import pysam
import tqdm
from statistics import (mode, StatisticsError)
import os
import sctools_imports


def get_read_family(record):
    tag_keys = ['XC', 'XM', 'gn']
    tag_values = [sctools_imports.get_tag_or_default(record, key, "") for key in tag_keys]
    return tag_values


def umi_collapse_sorted_file(input_bam_filename, output_bam_filename, verbose=False, total_reads=None,
                             debug=False):
    current_family = None
    current_family_size = 0
    family_index = 0
    temp_bam = None
    input_record_count = 0

    synthetic_read_prefix = "synthetic_read_"

    if verbose and (total_reads is not None):
        show_pbar = True
    else:
        show_pbar = False

    with pysam.AlignmentFile(input_bam_filename, "rb") as input_bam:
        with pysam.AlignmentFile(output_bam_filename, "wb", header=input_bam.header) as output_bam:
            if show_pbar:
                pbar = tqdm.tqdm(total=total_reads)
            for input_record in input_bam:
                input_record_count = input_record_count + 1
                current_read_family = get_read_family(input_record)
                if current_family == current_read_family:
                    # Continue outputting the same family
                    current_family_size = current_family_size + 1
                    temp_bam.write(input_record)
                else:
                    # Next family starting copy and reset
                    if temp_bam is not None:
                        temp_bam.close()
                        if current_family_size > 1:
                            call_consensus(output_bam, temp_bam_filename,
                                           new_read_name=f'{synthetic_read_prefix}{family_index}',
                                           temp_sorted_file= f'{temp_bam_filename}.sorted.bam')
                        else:
                            # copy read
                            with pysam.AlignmentFile(temp_bam_filename, "rb") as family_file:
                                first_read = family_file.__next__()
                                output_bam.write(first_read)
                        if not debug:
                            os.remove(temp_bam_filename)
                    family_index = family_index + 1
                    temp_bam_filename = f'tmp/family_{family_index}.bam'
                    temp_bam = pysam.AlignmentFile(temp_bam_filename, 'wb', header=input_bam.header)
                    current_family_size = 1
                    current_family = current_read_family
                    temp_bam.write(input_record)
                if show_pbar:
                    pbar.update(1)
            if temp_bam is not None:
                temp_bam.close()


def call_base(query_sequences, query_qualities):
    # TODO: Improve calling model
    try:
        base_call = mode(query_sequences)
    except StatisticsError:
        base_call = 'N'

    return base_call, 'F', 'M'


def call_consensus(output_bam, temp_bam_filename, new_read_name=None, temp_sorted_file=None):
    assert temp_sorted_file is not None
    assert new_read_name is not None

    # sort and index the family file
    pysam.sort(temp_bam_filename, "-o", temp_sorted_file)
    pysam.index(temp_sorted_file)

    with pysam.AlignmentFile(temp_sorted_file, "rb") as family_file:
        first_read = family_file.__next__()
        # reference_name = first_read.reference_name
        reference_id = first_read.reference_id
        tags = first_read.get_tags(with_value_type=True)

    new_read_sequence_list = []
    new_read_quality_list = []
    new_read_cigar_list = []

    first_pileup_position = None
    with pysam.AlignmentFile(temp_sorted_file, "rb") as family_file:
        for pileupcolumn in family_file.pileup():
            pos = pileupcolumn.pos
            qs = pileupcolumn.get_query_sequences()
            qq = pileupcolumn.get_query_qualities()

            if first_pileup_position is None:
                first_pileup_position = pos

            called_base, called_quality, called_cigar = call_base(qs, qq)

            new_read_sequence_list.append(called_base)
            new_read_quality_list.append(called_quality)
            new_read_cigar_list.append(called_cigar)


    # Construct new AlignedSegment
    nr = pysam.AlignedSegment()
    nr.query_name = new_read_name
    nr.query_sequence = ''.join(new_read_sequence_list)
    nr.flag = 0  # TODO
    nr.reference_id = reference_id
    nr.reference_start = first_pileup_position
    nr.mapping_quality = 255
    # TODO: Improve CIGAR String Handling
    #newcigarstring_debug = ''.join(new_read_cigar_list)
    #print(newcigarstring_debug)
    #nr.cigarstring = ''.join(new_read_cigar_list)
    n = len(new_read_cigar_list)
    nr.cigarstring = f'{n}M'
    nr.query_qualities = pysam.qualitystring_to_array(''.join(new_read_quality_list))
    nr.tags = tags

    # Write it to the output file
    output_bam.write(nr)


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
                             total_reads=total_number_of_reads, debug=debug)

    print('Cleaning up...')
    if not input_is_sorted and not debug:
        os.remove(tag_sorted_tmp_filename)


if __name__ == "__main__":
    #umi_collapse('final.EBOV.bam','output.bam',True, input_is_sorted = False)
    umi_collapse('tag_sorted_tmp.bam', 'output.bam', True, input_is_sorted=True)

