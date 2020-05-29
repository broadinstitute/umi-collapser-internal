import argparse
import umi_collapser


def main() -> int:
    description = "Collapse reads from the same Gene, UMI, Cell Barcode Triplet"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', "--input_bam", required=True, help="input bamfile")
    parser.add_argument('-o', "--output_bam", required=True, help="output bamfile")
    parser.add_argument("--input_is_sorted", dest="input_is_sorted", action='store_true', default=False)
    parser.add_argument("--cell_barcode_tag", default="XC")
    parser.add_argument("--molecular_barcode_tag", default="XM")
    parser.add_argument("--gene_tag", default="gn")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")

    args = parser.parse_args()

    umi_collapser.umi_collapse(input_file=args.input_bam,
                               output_file=args.output_bam,
                               verbose=args.verbose,
                               input_is_sorted=args.verbose,
                               cell_barcode_tag=args.cell_barcode_tag,
                               molecular_barcode_tag=args.molecular_barcode_tag,
                               gene_tag=args.gene_tag,
                               debug=args.debug)


if __name__ == "__main__":
    main()

