import pysam
import umi_collapser


# Example call
# debug_code.debug_single_family(family_bam="/Users/barkasn/Desktop/debug_1/family_87617_F.bam",
#                               output="/Users/barkasn/Desktop/debug/family_87617_F_consensus.bam")

def debug_single_family(family_bam: pysam.AlignmentFile,
                        output: str = "consensus_read.bam",
                        temp_sorted_filename: str = "temp.bam",
                        calling_method: str = "posterior",
                        ):
    """
    Helper function to process a single family bam for debugging purposes
    :param family_bam: input family bam
    :param output: output file with single consensus read
    :param temp_sorted_filename: temporary file for sorting input family
    :param calling_method: method for calling individual bases
    :return: None
    """
    nr = umi_collapser.call_consensus(family_bam=family_bam,
                                      new_read_name='new_read',
                                      temp_sorted_filename=temp_sorted_filename,
                                      calling_method=calling_method)
    with pysam.AlignmentFile(family_bam, "rb") as input_bam:
        with pysam.AlignmentFile(output, "wb", header=input_bam.header) as output_bam:
            output_bam.write(nr)

    return None
