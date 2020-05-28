import pysam
import tqdm
from statistics import (mode, StatisticsError)
import os
import sctools_imports

BAM_CMATCH = 0 # M
BAM_CREF_SKIP = 3 # N


def get_read_family(record):
    tag_keys = ['XC', 'XM', 'gn']
    tag_values = [sctools_imports.get_tag_or_default(record, key, "") for key in tag_keys]
    return tag_values


def umi_collapse_sorted_file(input_bam_filename, output_bam_filename, verbose=False, total_reads=None,
                             debug=False, synthetic_read_prefix="synthetic_read_",
                             family_temp_prefix="tmp/family_", save_collapsed_read=False,
                             save_collapsed_read_prefix="tmp/family_collapsed_"):
    current_family = None
    current_family_size = 0
    family_index = 0
    temp_bam = None
    input_record_count = 0

    if verbose and (total_reads is not None):
        show_progress_bar = True
    else:
        show_progress_bar = False

    with pysam.AlignmentFile(input_bam_filename, "rb") as input_bam:
        with pysam.AlignmentFile(output_bam_filename, "wb", header=input_bam.header) as output_bam:
            if show_progress_bar:
                pbar = tqdm.tqdm(total=total_reads)
            for input_record in input_bam:
                input_record_count = input_record_count + 1
                current_read_family = get_read_family(input_record)
                if current_family == current_read_family:
                    # same family
                    current_family_size = current_family_size + 1
                    temp_bam.write(input_record)
                else:
                    # next family
                    if temp_bam is not None:
                        temp_bam.close()
                        if current_family_size > 1:
                            new_read = call_consensus(temp_bam_filename,
                                                      new_read_name=f'{synthetic_read_prefix}{family_index}',
                                                      temp_sorted_file= f'{temp_bam_filename}.sorted.bam')
                            if save_collapsed_read:
                                with pysam.AlignmentFile(f'{save_collapsed_read_prefix}{family_index}', "wb",
                                                         header=input_bam.header) as family_collapsed_file:
                                    family_collapsed_file.write(new_read)
                            output_bam.write(new_read)
                        else:
                            # copy read
                            with pysam.AlignmentFile(temp_bam_filename, "rb") as family_file:
                                first_read = family_file.__next__()
                                output_bam.write(first_read)
                        if not debug:
                            os.remove(temp_bam_filename)
                    family_index = family_index + 1
                    # bam file for all reads in family
                    temp_bam_filename = f'{family_temp_prefix}{family_index}.bam'
                    temp_bam = pysam.AlignmentFile(temp_bam_filename, 'wb', header=input_bam.header)
                    current_family = current_read_family
                    temp_bam.write(input_record)
                    current_family_size = 1
                if show_progress_bar:
                    pbar.update(1)
            if temp_bam is not None:
                temp_bam.close()


def call_base(query_sequences, query_qualities):
    # TODO: Improve calling model
    try:
        base_call = mode([x.upper() for x in query_sequences])
    except StatisticsError:
        base_call = 'N'
    # TODO: Calculate quality (max observed)
    quality = 'F'
    return base_call, quality, 'M'


def call_consensus(temp_bam_filename, new_read_name=None, temp_sorted_file=None,
                   max_depth=10000, debug=False):
    assert temp_sorted_file is not None
    assert new_read_name is not None

    # sort and index the family file
    pysam.sort(temp_bam_filename, "-o", temp_sorted_file)
    pysam.index(temp_sorted_file)

    with pysam.AlignmentFile(temp_sorted_file, "rb") as family_file:
        first_read = family_file.__next__()
        reference_id = first_read.reference_id
        tags = first_read.get_tags(with_value_type=True)

    new_read_sequence_list = []
    new_read_quality_list = []
    new_read_cigar_tuple_list = []

    cigar_last = BAM_CMATCH
    cigar_last_count = 0

    last_pileup_position = None

    first_pileup_position = None
    with pysam.AlignmentFile(temp_sorted_file, "rb") as family_file:
        for pileupcolumn in family_file.pileup(stepper="nofilter", max_depth=max_depth):
            pos = pileupcolumn.pos

            if first_pileup_position is None:
                first_pileup_position = pos
                last_pileup_position = pos - 1

            position_delta = pos - last_pileup_position
            if position_delta > 1:
                # We have a gap
                if cigar_last == BAM_CREF_SKIP:
                    cigar_last_count += 1
                else:
                    new_read_cigar_tuple_list.append((cigar_last, cigar_last_count))
                    new_read_cigar_tuple_list.append((BAM_CREF_SKIP, position_delta - 1))
                    cigar_last = BAM_CREF_SKIP
                    cigar_last_count = position_delta

            qs = pileupcolumn.get_query_sequences()
            qq = pileupcolumn.get_query_qualities()

            called_base, called_quality, called_cigar = call_base(qs, qq)

            if called_base == "":
                if cigar_last == BAM_CREF_SKIP:
                    cigar_last_count += 1
                else:
                    new_read_cigar_tuple_list.append((cigar_last, cigar_last_count))
                    new_read_cigar_tuple_list.append((BAM_CREF_SKIP, position_delta - 1))
                    cigar_last = BAM_CREF_SKIP
                    cigar_last_count = position_delta
            else:
                new_read_sequence_list.append(called_base)
                new_read_quality_list.append(called_quality)

                if cigar_last == BAM_CMATCH:
                    cigar_last_count += 1
                else:
                    cigar_last = BAM_CMATCH
                    cigar_last_count = 1

            last_pileup_position = pos

        # append the final cigar tuple
        new_read_cigar_tuple_list.append((cigar_last, cigar_last_count))

    # Construct new AlignedSegment
    quality_string = ''.join(new_read_quality_list)
    quality_array = pysam.qualitystring_to_array(quality_string)

    nr = pysam.AlignedSegment()
    nr.query_name = new_read_name
    nr.query_sequence = ''.join(new_read_sequence_list)
    # TODO: Set the flags correctly
    nr.flag = 0
    nr.reference_id = reference_id
    nr.reference_start = first_pileup_position
    nr.mapping_quality = 255
    nr.cigartuples = new_read_cigar_tuple_list
    nr.query_qualities = quality_array
    nr.tags = tags

    if debug:
        print('---')
        print(new_read_name)
        print(nr.query_sequence)
        print(len(nr.query_sequence))
        print(new_read_cigar_tuple_list)
        print(quality_string)
        print(quality_array)

    return nr


def umi_collapse(input_file, output_file, verbose=False, input_is_sorted=False, debug=False):
    # TODO: Allow custom temp location and filename
    tag_sorted_tmp_filename = 'tag_sorted_tmp.bam'

    if not input_is_sorted:
        print('Sorting input bam file...')
        total_number_of_reads = tag_sort_bam(input_file, tag_sorted_tmp_filename, ['XC', 'XM', 'gn'])
        sorted_input_filename = tag_sorted_tmp_filename
    else:
        print('Counting input reads...')
        total_number_of_reads = int(pysam.view('-c', input_file).rstrip())
        sorted_input_filename = input_file

    # Do the collapse
    print('Collapsing reads...')
    umi_collapse_sorted_file(sorted_input_filename, output_file, verbose=verbose,
                             total_reads=total_number_of_reads, debug=debug,
                             save_collapsed_read=True)

    print('Cleaning up...')
    if not input_is_sorted and not debug:
        os.remove(tag_sorted_tmp_filename)
