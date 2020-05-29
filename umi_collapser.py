import pysam
import tqdm
from statistics import (mode, StatisticsError)
import os
import sctools_imports
from typing import (List, Any)

BAM_CMATCH = 0  # M
BAM_CREF_SKIP = 3  # N


def get_read_family(record: pysam.AlignedSegment, tag_keys: List[str]) -> List[Any]:
    """Return a list of tag values that define which family this read belongs to

    Parameters
    ----------
    record : pysam.AlignedSegment
        the read to the get family of
    tag_keys : List[str]
        list of tag name strings

    """
    tag_values = [sctools_imports.get_tag_or_default(record, key, "") for key in tag_keys]
    return tag_values


def umi_collapse_sorted_file(input_bam_filename: str, output_bam_filename: str, family_tag_keys: List[str],
                             verbose: bool = False,
                             total_reads: bool = None, debug: bool = False,
                             synthetic_read_prefix: str = "synthetic_read_",
                             family_temp_prefix: str = "tmp/family_", save_collapsed_read: bool = False,
                             save_collapsed_read_prefix: str = "tmp/family_collapsed_",
                             ):
    # TODO: reparametrise, allow to specify directory for temp files and prefixes
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
                current_read_family = get_read_family(input_record, family_tag_keys)
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
                                                      temp_sorted_filename=f'{temp_bam_filename}.sorted.bam')
                            if save_collapsed_read:
                                with pysam.AlignmentFile(f'{save_collapsed_read_prefix}{family_index}.bam', "wb",
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


def call_base(query_sequences: List[str], query_qualities:List[int]) -> List[str, str, str]:
    """
    call a base based on the query sequences and query qualities
    :param query_sequences: base calls from pileup for given position
    :param query_qualities: base qualities from pileup for given position
    :return: list with base call, quality and cigar string
    """
    # TODO: Resolve ties by looking at quality scores
    try:
        base_call = mode([x.upper() for x in query_sequences])
    except StatisticsError:
        base_call = 'N'
    # TODO: Calculate quality (max observed)
    quality = 'F'
    # TODO: If base_call is '' return 'N' for cigar
    return base_call, quality, 'M'


def call_consensus(temp_bam_filename: str, new_read_name: str = None, temp_sorted_filename: str = None,
                   max_depth: int = 10000, debug: bool = False) -> pysam.AlignedSegment:
    """
    call a consensus read from a read family file
    :param temp_bam_filename: name of file containing the family reads
    :param new_read_name: name of new read
    :param temp_sorted_filename: name of temporary file in which to store family reads
    :param max_depth: max depth parameter for pileup
    :param debug: debug mode
    :return: consensus read
    """

    assert temp_sorted_filename is not None
    assert new_read_name is not None

    # sort and index the family file
    pysam.sort(temp_bam_filename, "-o", temp_sorted_filename)
    pysam.index(temp_sorted_filename)

    with pysam.AlignmentFile(temp_sorted_filename, "rb") as family_file:
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
    with pysam.AlignmentFile(temp_sorted_filename, "rb") as family_file:
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
    nr.flag = 0
    nr.reference_id = reference_id
    nr.reference_start = first_pileup_position
    nr.mapping_quality = 255
    nr.cigartuples = new_read_cigar_tuple_list
    nr.query_qualities = quality_array
    nr.tags = tags

    if debug:
        print('--- DEBUG INFORMATION ---')
        print(f'new_read_name: {new_read_name}')
        print(f'query_sequence: {nr.query_sequence}')
        print(f'length of query: {len(nr.query_sequence)}')
        print(f'cigar tuples: {new_read_cigar_tuple_list}')
        print(f'quality string: {quality_string}')
        print(f'quality array: {quality_array}')

    return nr


def umi_collapse(input_file: str, output_file: str, verbose: bool = False, input_is_sorted: bool = False,
                 debug: bool = False, tag_sorted_tmp_filename: str = 'tag_sorted_tmp.bam',
                 cell_barcode_tag: str = "XC", molecular_barcode_tag: str = "XM", gene_tag: str = "gn") -> None:
    """

    :param input_file: input bam file
    :param output_file: output bam file
    :param verbose: specify verbosity
    :param input_is_sorted: specify if the input bam is already sorted
    :param debug: debug flag
    :param tag_sorted_tmp_filename: name of temporary file for sorting into
    :param cell_barcode_tag: tag name for the cell barcode
    :param molecular_barcode_tag: tag name for the molecular barcode
    :param gene_tag: tag name for the gene tag
    :return: None
    """
    tag_ordering = [cell_barcode_tag, molecular_barcode_tag, gene_tag]

    if not input_is_sorted:
        print('Sorting input bam file...')
        total_number_of_reads = sctools_imports.tag_sort_bam(input_file,
                                                             tag_sorted_tmp_filename,
                                                             tag_ordering)
        sorted_input_filename = tag_sorted_tmp_filename
    else:
        print('Counting input reads...')
        total_number_of_reads = int(pysam.view('-c', input_file).rstrip())
        sorted_input_filename = input_file

    print('Collapsing reads...')
    umi_collapse_sorted_file(input_bam_filename=sorted_input_filename,
                             output_bam_filename=output_file,
                             family_tag_keys=tag_ordering,
                             verbose=verbose,
                             total_reads=total_number_of_reads,
                             debug=debug,
                             save_collapsed_read=debug)

    print('Cleaning up...')
    if not input_is_sorted and not debug:
        os.remove(tag_sorted_tmp_filename)

    return None