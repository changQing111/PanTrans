import sys
import argparse
import logging
from .pipeline import unit_construct, unit_append
from .cnv import detect_cnv


def setup_logging(level=logging.INFO):
    """Configure logging format and level."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def read_parameters(argv=None):
    parser = argparse.ArgumentParser(
        description="A Scalable Approach for Constructing Pan-Genome and Pan-Transcriptome in Polyploid Organisms"
    )
    # subcommands
    subparsers = parser.add_subparsers(
        title="subcommand", description="construct, append, or cnv", dest="subcommand"
    )
    subparsers.required = True

    # construct subcommand
    construct_parser = subparsers.add_parser(
        "construct", help="Construct Pan-Genome and Pan-Transcriptome"
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
        "--rescue",
        action="store_true",
        help="rescue cluster representatives missing from the generated GTF",
    )
    construct_parser.add_argument(
        "-r", "--reference", required=True, help="reference name"
    )
    construct_parser.add_argument(
        "-t", "--threads", type=int, default=8, help="thread budget shared by minimap2 and samtools (minimum pipeline: 2; helper threads extra)"
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
        "-c",
        "--cdna",
        required=True,
        help="previous unrenamed final cDNA sequences path",
    )
    append_parser.add_argument(
        "--history-gtf",
        required=True,
        help="previous unrenamed final GTF path",
    )
    append_parser.add_argument(
        "--query-cdna",
        required=True,
        help="new-variety cDNA sequences path",
    )
    append_parser.add_argument(
        "-g", "--gdna", required=True, help="new gdna sequences path"
    )
    append_parser.add_argument(
        "-b",
        "--bed",
        required=True,
        help=(
            "combined BED containing historical pre-cluster representative gene IDs "
            "(the first gene of each pre cluster) plus new-variety genes"
        ),
    )
    append_parser.add_argument(
        "--history-graph",
        required=True,
        help="graph package emitted by the historical construct run",
    )
    append_parser.add_argument(
        "--rescue",
        action="store_true",
        help="rescue cluster representatives missing from the generated GTF",
    )
    append_parser.add_argument(
        "-p", "--prefix", default="Append", help="output prefix"
    )
    append_parser.add_argument(
        "-t", "--threads", type=int, default=8, help="thread budget shared by minimap2 and samtools (minimum pipeline: 2; helper threads extra)"
    )
    append_parser.add_argument(
        "--query-to-all-bam",
        default=None,
        help="reuse an existing new-cDNA to merged-gDNA BAM",
    )
    append_parser.add_argument(
        "--history-to-query-bam",
        default=None,
        help="reuse an existing historical-cDNA to new-gDNA BAM",
    )
    append_parser.add_argument(
        "-o", "--output", required=True, help="output dir"
    )

    cnv_parser = subparsers.add_parser(
        "cnv", help="Screen pre clusters for nearby same-variety gene copies"
    )
    cnv_parser.add_argument("--pre-cluster", required=True, help="pre.tmp.cluster file")
    cnv_parser.add_argument("-b", "--bed", required=True, help="complete gene BED with varieties in column 5")
    cnv_parser.add_argument("-r", "--reference", required=True, help="preferred representative variety")
    cnv_parser.add_argument("-o", "--output", required=True, help="output text file: one representative gene ID per hit cluster")
    cnv_parser.add_argument("--max-gap", type=int, default=1_000_000, help="maximum BED interval gap in bp (default: 1000000)")

    args = parser.parse_args(argv)
    if args.subcommand == "cnv" and args.max_gap < 0:
        cnv_parser.error("--max-gap must be nonnegative")
    
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
            rescue=args.rescue,
            threads=args.threads, 
            out_dir=args.output, 
            prefix=args.prefix, 
            refer_name=args.reference
        )
        return 0

    elif args.subcommand == "append":
        unit_append(
            history_cdna_path=args.cdna,
            history_gtf_path=args.history_gtf,
            query_cdna_path=args.query_cdna,
            query_gdna_path=args.gdna,
            all_bed_path=args.bed,
            history_graph_path=args.history_graph,
            rescue=args.rescue,
            threads=args.threads,
            out_dir=args.output,
            prefix=args.prefix,
            query_to_all_bam=args.query_to_all_bam,
            history_to_query_bam=args.history_to_query_bam,
        )
        return 0

    elif args.subcommand == "cnv":
        try:
            detect_cnv(args.pre_cluster, args.bed, args.reference, args.output, args.max_gap)
        except (ValueError, OSError) as exc:
            logging.getLogger(__name__).error("CNV screening failed: %s", exc)
            return 1
        return 0

    else:
        # No subcommand provided; log an error message.
        logger = logging.getLogger(__name__)
        logger.error("No subcommand specified. Use 'construct', 'append', or 'cnv'.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
