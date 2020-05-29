import pysam
import functools
from typing import (Iterator, Iterable, Optional)



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


def tag_sort_bam(input_bam_filename: str, output_bam_filename: str, tags):
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
