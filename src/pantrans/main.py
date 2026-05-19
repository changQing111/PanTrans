import sys
import argparse
import logging
from .pipeline import unit_construct, unit_append


def setup_logging(level=logging.INFO):
    """Configure logging format and level."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def read_parameters():
    parser = argparse.ArgumentParser(
        description="A Scalable Approach for Constructing Pan-Genome and Pan-Transcriptome in Polyploid Organisms"
    )
    # subcommands
    subparsers = parser.add_subparsers(
        title="subcommand", description="construct or append", dest="subcommand"
    )
    subparsers.required = True

    # construct subcommand
    construct_parser = subparsers.add_parser(
        "construct", help="Construct Pan-Genome and Pan-Transcriptome"
    )
    construct_parser.add_argument(
        "-n",
        "--name",
        required=True,
        help="variety names (file path, comma-separated list, or space-separated list)",
    )
    construct_parser.add_argument(
        "-c", "--cdna", required=True, default=None, help="all cdna sequences path"
    )
    construct_parser.add_argument(
        "-g", "--gdna", required=True, default=None, help="all gdna sequences path"
    )
    construct_parser.add_argument(
        "-b", "--bed", required=True, help="all bed path"
    )
    construct_parser.add_argument(
        "--bam",
        default=None,
        help="existing cDNA-to-gDNA BAM path; skip minimap2 alignment if provided",
    )
    construct_parser.add_argument(
        "--main-chroms",
        default=None,
        help="text file containing main chromosome names, one per line",
    )
    construct_parser.add_argument(
        "-r", "--reference", required=True, help="reference name"
    )
    construct_parser.add_argument(
        "-t", "--threads", type=int, default=8, help="number of threads for minimap2"
    )
    construct_parser.add_argument(
        "-p", "--prefix", default="Pan", required=True, help="output prefix"
    )
    construct_parser.add_argument(
        "-o", "--output", default="./", required=True, help="output dir"
    )

    # append subcommand
    append_parser = subparsers.add_parser(
        "append",
        help="Append New Variety to Existing Pan-Genome and Pan-Transcriptome",
    )
    append_parser.add_argument(
        "-n",
        "--name",
        required=True,
        help="new variety names (file path, comma-separated list, or space-separated list)",
    )
    append_parser.add_argument(
        "-c", "--cdna", required=True, help="new cdna sequences path"
    )
    append_parser.add_argument(
        "-g", "--gdna", required=True, help="new gdna sequences path"
    )
    append_parser.add_argument(
        "-b", "--bed", required=True, help="new bed file path"
    )
    append_parser.add_argument(
        "--refer_cdna", required=True, help="existing pan/reference cdna fasta path"
    )
    append_parser.add_argument(
        "--refer_gdna", required=True, help="existing pan/reference gdna fasta path"
    )
    append_parser.add_argument(
        "--refer_bed", required=True, help="existing pan/reference bed path"
    )
    append_parser.add_argument(
        "-p", "--prefix", default="Append", help="output prefix"
    )
    append_parser.add_argument(
        "-t", "--threads", type=int, default=8, help="number of threads for minimap2"
    )
    append_parser.add_argument(
        "-o", "--output", required=True, help="output dir"
    )

    args = parser.parse_args()
    
    return args


def main():
    setup_logging()
    args = read_parameters()

    if args.subcommand == "construct":
        unit_construct(
            all_cdna_path=args.cdna, 
            all_gdna_path=args.gdna, 
            all_bed_path=args.bed, 
            bam_path=args.bam,
            main_chrom_path=args.main_chroms,
            variety_li=args.name, 
            threads=args.threads, 
            out_dir=args.output, 
            prefix=args.prefix, 
            refer_name=args.reference
        )
        return 0

    elif args.subcommand == "append":
        unit_append(
            query_cdna_path=args.cdna,
            query_gdna_path=args.gdna,
            query_bed_path=args.bed,
            refer_cdna_path=args.refer_cdna,
            refer_gdna_path=args.refer_gdna,
            refer_bed_path=args.refer_bed,
            variety_name=args.name,
            threads=args.threads,
            out_dir=args.output,
            prefix=args.prefix,
        )
        return 0

    else:
        # No subcommand provided; log an error message.
        logger = logging.getLogger(__name__)
        logger.error("No subcommand specified. Use 'construct' or 'append'.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
