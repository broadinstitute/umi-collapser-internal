import pysam
import functools
from typing import (Iterator, Iterable, Optional)
import tqdm
from statistics import (mode, StatisticsError)
import os


# Taken from sctools

def get_tag_or_default(
        alignment: pysam.AlignedSegment, tag_key: str, default: Optional[str] = None
) -> Optional[str]:
    """Extracts the value associated to `tag_key` from `alignment`, and returns a default value
    if the tag is not present."""
    try:
        return alignment.get_tag(tag_key)
    except KeyError:
        return default


@functools.total_ordering
class TagSortableRecord(object):
    """Wrapper for pysam.AlignedSegment that facilitates sorting by tags and query name."""

    def __init__(
            self,
            tag_keys: Iterable[str],
            tag_values: Iterable[str],
            query_name: str,
            record: pysam.AlignedSegment = None,
    ) -> None:
        self.tag_keys = tag_keys
        self.tag_values = tag_values
        self.query_name = query_name
        self.record = record

    @classmethod
    def from_aligned_segment(
            cls, record: pysam.AlignedSegment, tag_keys: Iterable[str]
    ) -> "TagSortableRecord":
        """Create a TagSortableRecord from a pysam.AlignedSegment and list of tag keys"""
        assert record is not None
        tag_values = [get_tag_or_default(record, key, "") for key in tag_keys]
        query_name = record.query_name
        return cls(tag_keys, tag_values, query_name, record)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, TagSortableRecord):
            return NotImplemented
        self.__verify_tag_keys_match(other)
        for (self_tag_value, other_tag_value) in zip(self.tag_values, other.tag_values):
            if self_tag_value < other_tag_value:
                return True
            elif self_tag_value > other_tag_value:
                return False
        return self.query_name < other.query_name

    def __eq__(self, other: object) -> bool:
        # TODO: Add more error checking
        if not isinstance(other, TagSortableRecord):
            return NotImplemented
        self.__verify_tag_keys_match(other)
        for (self_tag_value, other_tag_value) in zip(self.tag_values, other.tag_values):
            if self_tag_value != other_tag_value:
                return False
        return self.query_name == other.query_name

    def __verify_tag_keys_match(self, other) -> None:
        if self.tag_keys != other.tag_keys:
            format_str = "Cannot compare records using different tag lists: {0}, {1}"
            raise ValueError(format_str.format(self.tag_keys, other.tag_keys))

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        format_str = "TagSortableRecord(tags: {0}, tag_values: {1}, query_name: {2}"
        return format_str.format(self.tag_keys, self.tag_values, self.query_name)


def sort_reads_by_tags_and_query_name(
        records: Iterable[pysam.AlignedSegment], tag_keys: Iterable[str]
) -> Iterable[pysam.AlignedSegment]:
    tag_sortable_records = (
        TagSortableRecord.from_aligned_segment(r, tag_keys) for r in records
    )
    sorted_records = sorted(tag_sortable_records)
    aligned_segments = (r.record for r in sorted_records)
    return aligned_segments


def tag_sort_bam(input_bam_filename, output_bam_filename, tags):
    number_of_reads = 0
    with pysam.AlignmentFile(input_bam_filename, "rb") as f:
        header = f.header
        records = f.fetch(until_eof=True)
        sorted_records = sort_reads_by_tags_and_query_name(records, tags)
    with pysam.AlignmentFile(output_bam_filename, "wb", header=header) as f:
        for record in sorted_records:
            f.write(record)
            number_of_reads = number_of_reads + 1
    return number_of_reads


def get_read_family(record):
    tag_keys = ['XC', 'XM', 'gn']
    tag_values = [get_tag_or_default(record, key, "") for key in tag_keys]
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


def main():
    input_file = 'final.EBOV.bam'
    output_file = 'output.bam'
    verbose = True

    tag_sorted_tmp_filename = 'tag_sorted_tmp.bam'

    # Sort the input file by tags
    # TODO: Sort the file into a ramdisk location
    print('Sorting input bam file...')
    total_number_of_reads = tag_sort_bam(input_file, tag_sorted_tmp_filename, ['XC', 'XM', 'gn'])

    # Do the collapse
    print('Collapsing reads...')
    umi_collapse_sorted_file(tag_sorted_tmp_filename, output_file, verbose=verbose, total_reads=total_number_of_reads)

    # TODO: remove the sorted file


def main_nosort():
    # input_file = 'final.EBOV.bam'
    output_file = 'output.unsorted.bam'
    verbose = True

    tag_sorted_tmp_filename = 'tag_sorted_tmp.bam'

    # Sort the input file by tags
    # TODO: Sort the file into a ramdisk location
    print('Sorting input bam file...')
    # total_number_of_reads = tag_sort_bam(input_file,tag_sorted_tmp_filename, ['XC','XM','gn'])
    total_number_of_reads = 1153569

    # Do the collapse
    print('Collapsing reads...')
    umi_collapse_sorted_file(tag_sorted_tmp_filename, output_file, verbose=verbose, total_reads=total_number_of_reads,
                             debug=True)

    # TODO: remove the sorted file


main_nosort()
