import sys
import argparse
import logging
from .cnv import detect_cnv

from . import cds_construct, cds_append


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
    from .version import __version__
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    # subcommands
    subparsers = parser.add_subparsers(
        title="subcommand", description="preprocess, construct, append, cnv, synteny", dest="subcommand"
    )
    subparsers.required = True

    prepare_parser = subparsers.add_parser("preprocess", help="Prepare gene DNA, full cDNA, all CDS and longest valid CDS per gene from genomes/GFF3")
    for direct, file_option, label in (("genomes", "genome-list", "uncompressed genome FASTAs"),
                                       ("annotations", "annotation-list", "GFF3 annotations"),
                                       ("varieties", "variety-list", "variety names")):
        group = prepare_parser.add_mutually_exclusive_group(required=True)
        group.add_argument("--" + direct, nargs="+", help="Ordered list of " + label)
        group.add_argument("--" + file_option, help="Text file, one entry per line: " + label)
    prepare_parser.add_argument("-o", "--output", required=True, help="output directory")
    prepare_parser.add_argument("-p", "--prefix", default="Pan", help="output filename prefix")
    prepare_parser.add_argument("--flank", type=int, default=3000, help="flank on each side, clipped at sequence boundaries (default: 3000)")
    prepare_parser.add_argument("--namespace-ids", action="store_true", help="prefix every gene ID with VARIETY__; emit original/output ID map")

    cds_construct_parser = subparsers.add_parser("construct", help="Build a CDS homology graph and UTR-aware transcriptome")
    cds_construct.configure_parser(cds_construct_parser)
    cds_append_parser = subparsers.add_parser("append", help="Extend a CDS graph with new prepared varieties")
    cds_append.configure_parser(cds_append_parser)

    cnv_parser = subparsers.add_parser(
        "cnv", help="Screen pre clusters for nearby same-variety gene copies"
    )
    cnv_parser.add_argument("--pre-cluster", required=True, help="pre.tmp.cluster file")
    cnv_parser.add_argument("-b", "--bed", required=True, help="complete gene BED with varieties in column 5")
    cnv_parser.add_argument("-r", "--reference", required=True, help="preferred representative variety")
    cnv_parser.add_argument("-o", "--output", required=True, help="output text file: one representative gene ID per hit cluster")
    cnv_parser.add_argument("--max-gap", type=int, default=1_000_000, help="maximum BED interval gap in bp (default: 1000000)")

    from . import synteny
    synteny_parser = subparsers.add_parser('synteny', help='Build pairwise gene-order synteny from reciprocal pre-cluster homologs')
    synteny.configure_parser(synteny_parser)

    args = parser.parse_args(argv)
    if args.subcommand == "preprocess":
        from .preprocess import read_list
        try:
            for direct, file_attr, paths in (("genomes", "genome_list", True),
                                             ("annotations", "annotation_list", True),
                                             ("varieties", "variety_list", False)):
                if getattr(args, file_attr):
                    setattr(args, direct, read_list(getattr(args, file_attr), paths=paths))
        except (ValueError, OSError) as exc:
            prepare_parser.error(str(exc))
        if args.flank < 0:
            prepare_parser.error("--flank must be nonnegative")
        if not args.genomes or len(args.genomes) != len(args.annotations) or len(args.genomes) != len(args.varieties):
            prepare_parser.error("Genome, annotation and variety lists must have equal nonzero lengths")
    if args.subcommand == "construct":
        cds_construct.validate_args(args, cds_construct_parser)
    elif args.subcommand == "append":
        cds_append.validate_args(args, cds_append_parser)
    if args.subcommand == "cnv" and args.max_gap < 0:
        cnv_parser.error("--max-gap must be nonnegative")
    if args.subcommand == "synteny":
        synteny.validate_args(args, synteny_parser)

    return args


def main():
    setup_logging()
    args = read_parameters()

    if args.subcommand == "preprocess":
        from .preprocess import preprocess
        from pysam import SamtoolsError
        try:
            paths = preprocess(args.genomes, args.annotations, args.varieties, args.output,
                               args.prefix, args.flank, args.namespace_ids)
        except (ValueError, OSError, SamtoolsError) as exc:
            logging.getLogger(__name__).error("Preprocessing failed: %s", exc)
            return 1
        logging.getLogger(__name__).info("Preprocessing complete: %s", paths["manifest"])
        return 0

    elif args.subcommand in ("construct", "append"):
        try:
            (cds_construct.run if args.subcommand == "construct" else cds_append.run)(args)
        except (ValueError, OSError) as exc:
            logging.getLogger(__name__).error("%s failed: %s", args.subcommand, exc)
            return 1
        return 0

    elif args.subcommand == "cnv":
        try:
            detect_cnv(args.pre_cluster, args.bed, args.reference, args.output, args.max_gap)
        except (ValueError, OSError) as exc:
            logging.getLogger(__name__).error("CNV screening failed: %s", exc)
            return 1
        return 0

    elif args.subcommand == "synteny":
        from . import synteny
        try:
            synteny.run(args)
        except (ValueError, OSError, KeyError, IndexError) as exc:
            logging.getLogger(__name__).error("Synteny analysis failed: %s", exc)
            return 1
        return 0

    else:
        # No subcommand provided; log an error message.
        logger = logging.getLogger(__name__)
        logger.error("No subcommand specified. Use 'construct', 'append', 'cnv', or 'synteny'.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
