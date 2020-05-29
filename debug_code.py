import pysam
import umi_collapser


def debug_single_family(family_bam,
                        output="consensus_read.bam",
                        temp_sorted_filename="temp.bam",
                        ):
    """
    Helper function to process a single family bam for debugging purposes
    :param family_bam: input family bam
    :param output: output file with single consensus read
    :param temp_sorted_filename: temporary file for sorting input family
    :return: None
    """
    nr = umi_collapser.call_consensus(family_bam=family_bam,
                                      new_read_name='new_read',
                                      temp_sorted_filename=temp_sorted_filename)
    with pysam.AlignmentFile(family_bam, "rb") as input_bam:
        with pysam.AlignmentFile(output, "wb", header=input_bam.header) as output_bam:
            output_bam.write(nr)

    return None
