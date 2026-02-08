import sys
import argparse
import logging
from .pipeline import unit_construct


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

    # construct subcommand
    construct_parser = subparsers.add_parser(
        "construct", help="Construct Pan-Genome and Pan-Transcriptome"
    )
    construct_parser.add_argument("-n", "--name", required=True, help="all variety name")
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
    append_parser.add_argument("-n", "--name", required=True, help="new variety name")
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
            variety_li=args.name, 
            threads=args.threads, 
            out_dir=args.output, 
            prefix=args.prefix, 
            refer_name=args.reference
        )
        return 0

    elif args.subcommand == "append":
        # Future work: implement the logic to append a new variety
        # to an existing pan-genome and pan-transcriptome.
        raise NotImplementedError(
            "The 'append' subcommand is not implemented yet. "
            "Please use 'construct' to build a new pan-genome/transcriptome."
        )

    else:
        # No subcommand provided; log an error message.
        logger = logging.getLogger(__name__)
        logger.error("No subcommand specified. Use 'construct' or 'append'.")
        return 1


if __name__ == "__main__":
    sys.exit(main())


