import pysam
import tqdm
import os
import sctools_imports
from typing import (List, Any)
import tempfile
from collections import Counter
from itertools import compress
from shutil import copyfile
import numpy as np
from scipy.special import logsumexp

# CIGAR String constants
BAM_CMATCH = 0  # M
BAM_CREF_SKIP = 3  # N


def get_read_family(record: pysam.AlignedSegment, tag_keys: List[str]) -> List[Any]:
    """
    Return a list of tag values that define which family this read belongs to
    :param record: the read to the get family of
    :param tag_keys: list of tag name strings
    :return: list of tag values
    """
    tag_values = [sctools_imports.get_tag_or_default(record, key, "") for key in tag_keys]
    return tag_values


def umi_collapse_sorted_file(input_bam_filename: str,
                             output_bam_filename: str,
                             family_tag_keys: List[str],
                             verbose: bool = False,
                             total_reads: bool = None,
                             synthetic_read_prefix: str = "synthetic_read_",
                             debug: bool = False,
                             debug_family_ids: List[int] = None,
                             debug_family_location: str = None,
                             ) -> None:
    """
    Perform family collapsing on a file that has been sorted
    :param input_bam_filename: input bam file name, must be sorted by tags
    :param output_bam_filename: output bam file name
    :param family_tag_keys: list of family tag keys
    :param verbose: verbosity
    :param total_reads: total number of reads in the file
    :param debug: debug flag
    :param synthetic_read_prefix: prefix for new reads
    :param debug_family_ids: list of ids of families to save for debug purposes
    :param debug_family_location: location where to save debug files
    :return: None
    """

    current_family = None
    current_family_size = 0
    current_family_forward_size = 0
    current_family_reverse_size = 0
    family_index = 0
    input_record_count = 0
    temp_bam_filename_forward = None
    temp_bam_forward = None
    family_file_prefix = "family_"

    with tempfile.TemporaryDirectory() as temp_directory_name:
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
                        if input_record.is_reverse:
                            temp_bam_reverse.write(input_record)
                            current_family_reverse_size += 1
                        else:
                            temp_bam_forward.write(input_record)
                            current_family_forward_size += 1
                        # same family
                        current_family_size += 1
                        temp_bam_forward.write(input_record)
                    else:
                        # next family
                        if temp_bam_forward is not None:
                            temp_bam_forward.close()
                            temp_bam_reverse.close()
                            call_family_consensus(current_family_forward_size, current_family_reverse_size,
                                                  current_family_size, debug, debug_family_ids, debug_family_location,
                                                  family_file_prefix, family_index, input_bam, output_bam,
                                                  synthetic_read_prefix, temp_bam_filename_forward,
                                                  temp_bam_filename_reverse)
                            # Cleanup
                            os.remove(temp_bam_filename_forward)
                            os.remove(temp_bam_filename_reverse)
                        family_index = family_index + 1
                        # create new bam for this family
                        temp_bam_filename_forward = f'{temp_directory_name}/{family_file_prefix}{family_index}_F.bam'
                        temp_bam_filename_reverse = f'{temp_directory_name}/{family_file_prefix}{family_index}_R.bam'
                        temp_bam_forward = pysam.AlignmentFile(temp_bam_filename_forward, 'wb', header=input_bam.header)
                        temp_bam_reverse = pysam.AlignmentFile(temp_bam_filename_reverse, 'wb', header=input_bam.header)
                        current_family = current_read_family
                        if input_record.is_reverse:
                            temp_bam_reverse.write(input_record)
                            current_family_forward_size = 0
                            current_family_reverse_size = 1
                        else:
                            temp_bam_forward.write(input_record)
                            current_family_forward_size = 1
                            current_family_reverse_size = 0
                        current_family_size = 1
                    if show_progress_bar:
                        pbar.update(1)
            if temp_bam_forward is not None:
                temp_bam_forward.close()
                temp_bam_reverse.close()
                call_family_consensus(current_family_forward_size, current_family_reverse_size,
                                      current_family_size, debug, debug_family_ids, debug_family_location,
                                      family_file_prefix, family_index, input_bam, output_bam,
                                      synthetic_read_prefix, temp_bam_filename_forward,
                                      temp_bam_filename_reverse)


# TODO: Add type hinting
def call_family_consensus(current_family_forward_size, current_family_reverse_size, current_family_size, debug,
                          debug_family_ids, debug_family_location, family_file_prefix, family_index, input_bam,
                          output_bam, synthetic_read_prefix, temp_bam_filename_forward, temp_bam_filename_reverse) -> None:
    """
    Call consensus read for family identifying if collapse is required and
    selecting if forward or reverse orientation should be used
    :param current_family_forward_size: number of reads in forward orientation for this family
    :param current_family_reverse_size: number of reads in reverse orientation for this family
    :param current_family_size: total family size
    :param debug: debug mode
    :param debug_family_ids: families for which to generate debug files
    :param debug_family_location: location where to save the debug files
    :param family_file_prefix: prefix of family files
    :param family_index: index for the family
    :param input_bam: input bam file (opened)
    :param output_bam: output bam files
    :param synthetic_read_prefix: prefix for synthetic reads
    :param temp_bam_filename_forward: temp filename of file with forward reads
    :param temp_bam_filename_reverse: temp filenem of file with reverse reads
    :return: None
    """
    if current_family_size > 1:
        if current_family_forward_size >= current_family_reverse_size:
            new_read = call_consensus(temp_bam_filename_forward,
                                      new_read_name=f'{synthetic_read_prefix}{family_index}',
                                      temp_sorted_filename=f'{temp_bam_filename_forward}'
                                      f'.sorted.bam')
            output_bam.write(new_read)
        else:
            new_read = call_consensus(temp_bam_filename_reverse,
                                      new_read_name=f'{synthetic_read_prefix}{family_index}',
                                      temp_sorted_filename=f'{temp_bam_filename_forward}'
                                      f'.sorted.bam')
            output_bam.write(new_read)
    else:
        if current_family_forward_size == 1:
            # copy read, could cache last read and avoid re-opening file
            with pysam.AlignmentFile(temp_bam_filename_forward, "rb") as family_file:
                first_read = family_file.__next__()
                output_bam.write(first_read)
        else:
            # copy read, could cache last read and avoid re-opening file
            with pysam.AlignmentFile(temp_bam_filename_reverse, "rb") as family_file:
                first_read = family_file.__next__()
                output_bam.write(first_read)
    if debug:
        # save information about specific families
        if family_index in debug_family_ids:
            save_family_debug(debug_family_location, family_file_prefix, family_index,
                              input_bam, new_read, temp_bam_filename_forward,
                              temp_bam_filename_reverse)

# TODO: Add type hinting
def save_family_debug(debug_family_location, family_file_prefix, family_index, input_bam, new_read,
                      temp_bam_filename_forward, temp_bam_filename_reverse) -> None:
    """
    Save debug information for family
    :param debug_family_location: location where to save the files
    :param family_file_prefix: prefix for files to save
    :param family_index: family index number
    :param input_bam: input bam file
    :param new_read: read to write in output for this family
    :param temp_bam_filename_forward: file name for forward reads in this family
    :param temp_bam_filename_reverse: file name for reverse reads in this family
    :return: None
    """
    # save information about this family in the target location
    copyfile(temp_bam_filename_forward, f"{debug_family_location}/"
    f"{family_file_prefix}{family_index}_F.bam")
    copyfile(temp_bam_filename_reverse, f"{debug_family_location}/"
    f"{family_file_prefix}{family_index}_R.bam")
    collapsed_read_bam_filename = f"{debug_family_location}" + \
                                  f"{family_file_prefix}{family_index}_collapsed.bam"
    with pysam.AlignmentFile(collapsed_read_bam_filename,
                             "wb", header=input_bam.header) as collapsed_read_bam:
        collapsed_read_bam.write(new_read)

def call_base(query_sequences: List[str], query_qualities: List[int], method="posterior") -> List[str]:
    if method == "majority":
        return call_base_majority_vote(query_sequences, query_qualities)
    elif method == "posterior":
        return call_base_posterior(query_qualities, query_qualities)
    else:
        raise Exception("Unknown method for base calling")


def get_log_prob_compl(ln_prob):
    """Get the natural base log complement"""
    ln1 = np.log(1)
    return ln1 + np.log(1 - np.exp(ln_prob - ln1))

LN_1_M_EXP_THRESHOLD = -np.log(2.)

def get_log_prob_compl_stable(ln_prob):
    """Get the natural base log complement"""
    return np.where(
        ln_prob > LN_1_M_EXP_THRESHOLD,
        np.log(-np.expm1(ln_prob)),
        np.log1p(-np.exp(ln_prob))
    )


def call_base_posterior(query_sequences: List[str], query_qualities: List[int], max_quality_score = 40) -> List[str]:

    assert len(query_sequences) == len(query_qualities)

    # k is the size of the pileup here
    k = len(query_sequences)

    # array of possible bases
    bases = ['A', 'T', 'G', 'C']

    # Calculate p of each base being correct or wrong in log_e
    q_i = np.array(query_qualities)
    ln_p_error = np.log(np.divide(np.power(10, np.divide(q_i, -10)), 3))
    ln_p_correct = get_log_prob_compl_stable(ln_p_error)

    # Construct arrays to manipulate all data together
    p_d_base = np.empty([k, 4])
    p_d_not_base = np.empty([k, 4])
    for j in range(len(bases)):
        for i in range(k):
            if query_sequences[i] == bases[j]:
                p_d_base[i, j] = ln_p_correct[i]
                p_d_not_base[i, j] = ln_p_error[i]
            else:
                p_d_base[i, j] = ln_p_error[i]
                p_d_not_base[i, j] = ln_p_correct[i]

    # Calculate p of data given each underlying base
    p_d_base_sum = np.sum(p_d_base, 0)
    p_d_not_base_sum = np.sum(p_d_not_base, 0)

    # Calculate Posterior
    nominator = p_d_base_sum + np.log(np.divide(1, 4))
    denominator = logsumexp(
        np.vstack(
            (
                nominator,
                np.add(
                    p_d_not_base_sum,
                    np.full(
                        [1, 4],
                        np.log(np.divide(3, 4))
                    )
                )
            )
        ),
        0
    )

    log_p = nominator - denominator
    log_p_norm = log_p - logsumexp(log_p)

    # Call the base with max p
    call_i = np.argmax(log_p_norm)

    if np.sum(log_p_norm == log_p_norm[call_i]) == 1:
        basecall = bases[call_i]

        # Probability of incorrect call is sum of p that another base is correct
        log_p_incorrect_call = logsumexp(log_p_norm[np.arange(len(log_p_norm)) != call_i])

        # Convert log_e probability of error to Phred scale
        quality_score = int(np.multiply(-10, np.log10(np.exp(log_p_incorrect_call))))

        # Cap quality score
        quality_score = (quality_score, 0)[quality_score < 0]
        quality_score = (quality_score, max_quality_score)[quality_score > max_quality_score]

        cigar_value = BAM_CMATCH
    else:
        # There is a tie, don't call base
        basecall = 'N'
        quality_score = 0
        cigar_value = BAM_CMATCH

    return [basecall, chr(quality_score + 33), BAM_CMATCH]


def call_base_majority_vote(query_sequences: List[str], query_qualities: List[int]) -> List[str]:
    """
    call a base based on the query sequences and query qualities
    :param query_sequences: base calls from pileup for given position
    :param query_qualities: base qualities from pileup for given position
    :return: list with base call, quality and cigar string
    """
    if len(query_sequences) == 0:
        return ['', chr(33), BAM_CREF_SKIP]
    else:
        query_sequences_upper = [x.upper() for x in query_sequences]
        base_frequencies = Counter(query_sequences_upper)
        most_common = base_frequencies.most_common(1)[0]
        most_common_base = most_common[0]
        most_common_freq = most_common[1]
        tie_exists = False

        # Find if tie exists
        for x in base_frequencies:
            if x != most_common_base and base_frequencies[x] == most_common_freq:
                tie_exists = True
                break

        if not tie_exists:
            base_call = most_common[0]
            select_vector = [qs == base_call for qs in query_sequences_upper]
            quality_score = max(list(compress(query_qualities, select_vector)))
        else:
            # Find which bases tie
            tied_bases = [x for x in base_frequencies if base_frequencies[x] == most_common_freq]
            tied_total_qualities = list()

            # Get total score for each tied base
            for tb in tied_bases:
                select_vector = [qs == tb for qs in query_sequences_upper]
                current_tied_base_qualities = list(compress(query_qualities, select_vector))
                tied_total_qualities.append((tb, sum(current_tied_base_qualities)))

            # Find if a base with maximal score exists
            max_quality = max([x[1] for x in tied_total_qualities])
            if sum([x[1] == max_quality for x in tied_total_qualities]) == 1:
                base_call = [x[0] for x in tied_total_qualities if x[1] == max_quality][0]
                select_vector = [qs == base_call for qs in query_sequences_upper]
                quality_score = max(list(compress(query_qualities, select_vector)))
            else:
                # qualities tie, can't call base
                base_call = 'N'
                quality_score = 0

        # Check if the most common base was ''
        cigar_value = BAM_CMATCH
        if base_call == "":
            cigar_value = BAM_CREF_SKIP

        return [base_call, chr(quality_score + 33), cigar_value]


def call_consensus(family_bam: str,
                   new_read_name: str = None,
                   temp_sorted_filename: str = None,
                   max_depth: int = 10000,
                   debug: bool = False,
                   debug_keep_families=False) -> pysam.AlignedSegment:
    """
    call a consensus read from a read family file
    :param family_bam: name of file containing the family reads
    :param new_read_name: name of new read
    :param temp_sorted_filename: name of temporary file in which to store family reads
    :param max_depth: max depth parameter for pileup
    :param debug: debug mode
    :return: consensus read
    """

    assert temp_sorted_filename is not None
    assert new_read_name is not None

    # sort and index the family file
    pysam.sort(family_bam, "-o", temp_sorted_filename)
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
        for pileup_column in family_file.pileup(stepper="nofilter", max_depth=max_depth, min_base_quality=0):
            pos = pileup_column.pos

            if first_pileup_position is None:
                first_pileup_position = pos
                last_pileup_position = pos - 1

            position_delta = pos - last_pileup_position
            if position_delta > 1:
                # We have a gap
                if cigar_last == BAM_CREF_SKIP:
                    # If we are already in a gap extend it
                    cigar_last_count += position_delta - 1
                else:
                    # If we are not in a gap close the previous segment and start a new
                    new_read_cigar_tuple_list.append((cigar_last, cigar_last_count))
                    cigar_last = BAM_CREF_SKIP
                    cigar_last_count = position_delta - 1

            query_sequences = pileup_column.get_query_sequences()
            query_qualities = pileup_column.get_query_qualities()

            called_base, called_quality, called_cigar = call_base(query_sequences, query_qualities)

            if called_cigar == BAM_CREF_SKIP:
                # No base could be called and we have a single skip
                if cigar_last == BAM_CREF_SKIP:
                    # If we are already in a skip extend it
                    cigar_last_count += 1
                else:
                    # otherwise close previous segment and start a skip
                    new_read_cigar_tuple_list.append((cigar_last, cigar_last_count))
                    cigar_last = BAM_CREF_SKIP
                    cigar_last_count = 1
            else:
                # we have a base call
                new_read_sequence_list.append(called_base)
                new_read_quality_list.append(called_quality)

                if cigar_last == BAM_CMATCH:
                    cigar_last_count += 1
                else:
                    new_read_cigar_tuple_list.append((cigar_last, cigar_last_count))
                    cigar_last = BAM_CMATCH
                    cigar_last_count = 1

            # update the last position for which we got a pileup
            last_pileup_position = pos

        # append the final cigar tuple
        new_read_cigar_tuple_list.append((cigar_last, cigar_last_count))

    # Construct new AlignedSegment
    quality_string = ''.join(new_read_quality_list)
    quality_array = pysam.qualitystring_to_array(quality_string)

    new_read = pysam.AlignedSegment()
    new_read.query_name = new_read_name
    new_read.query_sequence = ''.join(new_read_sequence_list)
    new_read.flag = 0
    new_read.reference_id = reference_id
    new_read.reference_start = first_pileup_position
    new_read.mapping_quality = 255
    new_read.cigartuples = new_read_cigar_tuple_list
    new_read.query_qualities = quality_array
    new_read.tags = tags

    return new_read


def umi_collapse(input_file: str,
                 output_file: str,
                 verbose: bool = False,
                 input_is_sorted: bool = False,
                 debug: bool = False,
                 tag_sorted_tmp_filename: str = 'tag_sorted_tmp.bam',
                 cell_barcode_tag: str = "XC",
                 molecular_barcode_tag: str = "XM",
                 gene_tag: str = "gn") -> None:
    """
    Collapse all umi families
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
                             )

    print('Cleaning up...')
    if not input_is_sorted and not debug:
        os.remove(tag_sorted_tmp_filename)

    return None
