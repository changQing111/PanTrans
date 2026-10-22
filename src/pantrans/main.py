import sys
import argparse
import logging
from .pipeline import unit_construct, unit_append
from .cnv import detect_cnv
from .translocation import detect_translocations


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
        title="subcommand", description="construct, append, cnv, or translocation", dest="subcommand"
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
        help="rescue cluster representatives missing from the generated GFF3",
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
        "--history-gff3", "--history-gtf",
        required=True,
        help="previous unrenamed final GFF3 (legacy GTF also accepted; --history-gtf is a compatibility alias)",
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
        help="rescue cluster representatives missing from the generated GFF3",
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

    for coding_parser in (construct_parser, append_parser):
        coding_parser.add_argument("--transdecoder-parallel", action="store_true", help="parallelize CDS scoring using -t CPUs (default: serial); requires CDS annotation")
        coding_parser.add_argument("--annotate-cds", action="store_true", help="predict CDS and UTR features in final GFF3s using TransDecoder")
        coding_parser.add_argument("--transdecoder-dir", help="TransDecoder executable directory (default: search PATH); implies --annotate-cds")
        coding_parser.add_argument("--min-orf-length", type=int, default=50, help="minimum candidate ORF length in amino acids (default: 50)")

    cnv_parser = subparsers.add_parser(
        "cnv", help="Screen pre clusters for nearby same-variety gene copies"
    )
    cnv_parser.add_argument("--pre-cluster", required=True, help="pre.tmp.cluster file")
    cnv_parser.add_argument("-b", "--bed", required=True, help="complete gene BED with varieties in column 5")
    cnv_parser.add_argument("-r", "--reference", required=True, help="preferred representative variety")
    cnv_parser.add_argument("-o", "--output", required=True, help="output text file: one representative gene ID per hit cluster")
    cnv_parser.add_argument("--max-gap", type=int, default=1_000_000, help="maximum BED interval gap in bp (default: 1000000)")

    translocation_parser = subparsers.add_parser(
        "translocation", help="Screen reciprocal pre-cluster homologs for chromosome changes"
    )
    translocation_parser.add_argument("--pre-cluster", required=True, help="pre.tmp.cluster file")
    translocation_parser.add_argument("-b", "--bed", required=True, help="complete gene BED with varieties in column 5")
    translocation_parser.add_argument("--edges", required=True, help="graph edge TSV containing directed gene edges")
    translocation_parser.add_argument("--main-chroms", required=True, help="main chromosome names, one per line")
    translocation_parser.add_argument("-r", "--reference", required=True, help="reference variety")
    translocation_parser.add_argument("-o", "--output-prefix", required=True, help="prefix for .translocation.genes.tsv and .translocation.blocks.tsv")
    translocation_parser.add_argument("--min-anchor-genes", type=int, default=3, help="minimum independent pre-cluster anchors for a supported block (default: 3)")
    translocation_parser.add_argument("--max-gap", type=int, default=1_000_000, help="maximum interval gap on either chromosome in bp (default: 1000000)")
    translocation_parser.add_argument("--max-skipped-genes", type=int, default=5, help="maximum intervening BED genes allowed per side (default: 5)")
    translocation_parser.add_argument("--allow-reverse", action="store_true", help="also report collinear blocks in reverse orientation")
    translocation_parser.add_argument("--include-unplaced", action="store_true", help="include pairs with unplaced genes in reports and ambiguity counts (never as block anchors)")

    args = parser.parse_args(argv)
    if args.subcommand in ("construct", "append"):
        if args.threads < 1:
            parser.error("--threads must be positive")
        if args.min_orf_length < 1:
            parser.error("--min-orf-length must be positive")
        args.annotate_cds = args.annotate_cds or bool(args.transdecoder_dir)
    if args.subcommand == "cnv" and args.max_gap < 0:
        cnv_parser.error("--max-gap must be nonnegative")
    if args.subcommand == "translocation":
        if args.min_anchor_genes < 2:
            translocation_parser.error("--min-anchor-genes must be at least 2")
        if args.max_gap < 0 or args.max_skipped_genes < 0:
            translocation_parser.error("--max-gap and --max-skipped-genes must be nonnegative")
    
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
            annotate_cds=args.annotate_cds,
            transdecoder_dir=args.transdecoder_dir,
            min_orf_length=args.min_orf_length,
            transdecoder_parallel=args.transdecoder_parallel,
            threads=args.threads, 
            out_dir=args.output, 
            prefix=args.prefix, 
            refer_name=args.reference
        )
        return 0

    elif args.subcommand == "append":
        unit_append(
            history_cdna_path=args.cdna,
            history_gff3_path=args.history_gff3,
            query_cdna_path=args.query_cdna,
            query_gdna_path=args.gdna,
            all_bed_path=args.bed,
            history_graph_path=args.history_graph,
            rescue=args.rescue,
            annotate_cds=args.annotate_cds,
            transdecoder_dir=args.transdecoder_dir,
            min_orf_length=args.min_orf_length,
            transdecoder_parallel=args.transdecoder_parallel,
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

    elif args.subcommand == "translocation":
        try:
            detect_translocations(
                args.pre_cluster, args.bed, args.edges, args.main_chroms,
                args.reference, args.output_prefix, args.min_anchor_genes,
                args.max_gap, args.max_skipped_genes, args.allow_reverse,
                include_unplaced=args.include_unplaced,
            )
        except (ValueError, OSError) as exc:
            logging.getLogger(__name__).error("Translocation screening failed: %s", exc)
            return 1
        return 0

    else:
        # No subcommand provided; log an error message.
        logger = logging.getLogger(__name__)
        logger.error("No subcommand specified. Use 'construct', 'append', 'cnv', or 'translocation'.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
