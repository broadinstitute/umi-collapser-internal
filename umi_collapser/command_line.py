import argparse
import umi_collapser

def main() -> int:
    description = "Collapse reads from the same Gene, UMI, Cell Barcode Triplet"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i',
                        "--input_bam",
                        required=True,
                        help="input bamfile")
    parser.add_argument('-o',
                        "--output_bam",
                        required=True,
                        help="output bamfile")
    parser.add_argument("--input_is_sorted",
                        dest="input_is_sorted",
                        action='store_true',
                        default=False,
                        help="flag indicating if the input bam is sorted by tags")
    parser.add_argument("--cell_barcode_tag",
                        default="XC",
                        help="tag name for the cell barcode")
    parser.add_argument("--molecular_barcode_tag",
                        default="XM",
                        help="tag name for the molecular barcode")
    parser.add_argument("--gene_tag",
                        default="gn",
                        help="tag name for the gene tag")
    parser.add_argument("--calling_method",
                        default="posterior",
                        choices=['posterior', 'majority'],
                        help="method to use to call individual bases")
    parser.add_argument("--verbose",
                        action="store_true",
                        help="verbosity level")
    parser.add_argument("--debug",
                        action="store_true",
                        help="flag for debug mode")

    args = parser.parse_args()

    umi_collapser.do_umi_collapse(input_file=args.input_bam,
                               output_file=args.output_bam,
                               verbose=args.verbose,
                               input_is_sorted=args.input_is_sorted,
                               cell_barcode_tag=args.cell_barcode_tag,
                               molecular_barcode_tag=args.molecular_barcode_tag,
                               gene_tag=args.gene_tag,
                               debug=args.debug,
                               calling_method=args.calling_method)
    return 0
