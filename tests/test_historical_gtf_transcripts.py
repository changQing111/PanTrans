import logging
import tempfile
import unittest
from pathlib import Path

import pysam

from pantrans import transcript_processor
from pantrans.pipeline import _write_cluster_outputs
from pantrans.transcript_processor import (
    load_gtf_transcript_models,
    transcript_dedup,
)


def _write_gtf(path, lines):
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _write_bam(path, alignments=()):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "Ref.g1", "LN": 100}]}
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for query_name, reference_start, cigar in alignments:
            read = pysam.AlignedSegment()
            read.query_name = query_name
            read.query_sequence = "A" * sum(
                length for op, length in cigar if op in (0, 1, 4, 7, 8)
            )
            read.flag = 0
            read.reference_id = 0
            read.reference_start = reference_start
            read.mapping_quality = 60
            read.cigartuples = cigar
            read.query_qualities = pysam.qualitystring_to_array(
                "I" * len(read.query_sequence)
            )
            read.set_tag("AS", len(read.query_sequence))
            read.set_tag("NM", 0)
            output.write(read)
    return str(path)


def _parse_output_gtf(path):
    transcripts = []
    exons = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        attrs = {
            item.strip().split(" ", 1)[0]: item.strip().split(" ", 1)[1].strip('"')
            for item in fields[8].split(";")
            if item.strip()
        }
        transcript_id = attrs["transcript_id"]
        if fields[2] == "transcript":
            transcripts.append(transcript_id)
        elif fields[2] == "exon":
            exons.setdefault(transcript_id, []).append(
                (int(fields[3]), int(fields[4]))
            )
    return transcripts, exons


class HistoricalGtfTranscriptTest(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.tmp_path = Path(self.temp_dir.name)

    def test_rename_gtf_ids_preserves_coordinates_and_transcript_suffix(self):
        source_gtf = _write_gtf(
            self.tmp_path / "source.gtf",
            [
                'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\ttranscript_id "Ref.g1.2"; gene_id "Ref.g1"; note "keep";',
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\ttranscript_id "Ref.g1.2"; gene_id "Ref.g1"; note "keep";',
                'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\ttranscript_id "Ref.g1.2"; gene_id "Ref.g1"; note "keep";',
            ],
        )
        output_gtf = self.tmp_path / "renamed.gtf"

        transcript_processor.rename_gtf_ids(
            source_gtf,
            output_gtf,
            {"Ref.g1": "Pan1A000001"},
        )

        with open(source_gtf, "rt") as handle:
            source_fields = [line.split("\t") for line in handle.read().splitlines()]
        output_fields = [
            line.split("\t") for line in output_gtf.read_text().splitlines()
        ]
        self.assertEqual(
            [fields[:8] for fields in output_fields],
            [fields[:8] for fields in source_fields],
        )
        self.assertTrue(
            all('gene_id "Pan1A000001"' in fields[8] for fields in output_fields)
        )
        self.assertTrue(
            all(
                'transcript_id "Pan1A000001.2"' in fields[8]
                for fields in output_fields
            )
        )
        self.assertTrue(all('note "keep"' in fields[8] for fields in output_fields))

    def test_rename_gtf_ids_rejects_non_numeric_transcript_suffix(self):
        source_gtf = _write_gtf(
            self.tmp_path / "invalid_suffix.gtf",
            [
                'Ref.g1\tPan\ttranscript\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.isoform";',
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.isoform";',
            ],
        )

        with self.assertRaisesRegex(ValueError, "numeric suffix"):
            transcript_processor.rename_gtf_ids(
                source_gtf,
                self.tmp_path / "renamed.gtf",
                {"Ref.g1": "Pan1A000001"},
            )

    def test_load_gtf_transcript_models_sorts_exons_and_derives_splice_sites(self):
        history_gtf = _write_gtf(
            self.tmp_path / "history.gtf",
            [
                'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            ],
        )

        models = load_gtf_transcript_models(history_gtf)

        self.assertEqual(models["transcript_gene"], {"Ref.g1.1": "Ref.g1"})
        self.assertEqual(models["exon_coords"], {"Ref.g1.1": [(1, 10), (21, 30)]})
        self.assertEqual(models["splice_sites"], {"Ref.g1.1": [(11, 20)]})

    def test_load_gtf_transcript_models_requires_transcript_feature(self):
        history_gtf = _write_gtf(
            self.tmp_path / "exon_only.gtf",
            [
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            ],
        )

        with self.assertRaisesRegex(ValueError, "no transcript feature"):
            load_gtf_transcript_models(history_gtf)

    def test_load_gtf_transcript_models_rejects_transcript_gene_id_mismatch(self):
        history_gtf = _write_gtf(
            self.tmp_path / "mismatched_transcript_gene.gtf",
            [
                'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Other.g1.1";',
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Other.g1.1";',
                'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Other.g1.1";',
            ],
        )

        with self.assertRaisesRegex(ValueError, "transcript_id.*gene_id"):
            load_gtf_transcript_models(history_gtf)

    def test_load_gtf_transcript_models_rejects_malformed_models(self):
        malformed_cases = [
            (
                [
                    'Ref.g1\tPan\ttranscript\t1\t10\t.\t+\t.\tgene_id "Other.g1"; transcript_id "Ref.g1.1";',
                    'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Other.g1"; transcript_id "Ref.g1.1";',
                ],
                "sequence name.*gene_id",
            ),
            (
                [
                    'Ref.g1\tPan\ttranscript\t1\t10\t.\t+\t.\tgene_id "Ref.g1";',
                    'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                ],
                "transcript_id",
            ),
            (
                [
                    'Ref.g1\tPan\ttranscript\t1\t20\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                    'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                    'Ref.g1\tPan\texon\t10\t20\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                ],
                "overlapping or adjacent exons",
            ),
            (
                [
                    'Ref.g1\tPan\ttranscript\t1\t20\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                    'Ref.g1\tPan\texon\t1\t5\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                    'Ref.g1\tPan\texon\t11\t20\t.\t-\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                ],
                "mixed sequence names or strands",
            ),
        ]

        for lines, message in malformed_cases:
            with self.subTest(message=message):
                history_gtf = _write_gtf(self.tmp_path / "malformed.gtf", lines)
                with self.assertRaisesRegex(ValueError, message):
                    load_gtf_transcript_models(history_gtf)

    def test_transcript_dedup_preserves_seed_transcript_absent_from_bam(self):
        history_gtf = _write_gtf(
            self.tmp_path / "history.gtf",
            [
                'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            ],
        )
        bam_path = _write_bam(self.tmp_path / "empty.bam")
        output_gtf = self.tmp_path / "output.gtf"

        transcript_dedup(
            bam_path,
            cluster_dic={"Ref.g1": ["Ref.g1", "JM22.g1"]},
            trans_len_dic={"Ref.g1.1": 20},
            gene_len_dic={"Ref.g1": 100},
            gene_strand_dic={"Ref.g1": "+"},
            rename_map=None,
            gtf_path=str(output_gtf),
            seed_gtf_path=history_gtf,
        )

        transcripts, exons = _parse_output_gtf(output_gtf)
        self.assertEqual(transcripts, ["Ref.g1.1"])
        self.assertEqual(exons, {"Ref.g1.1": [(1, 10), (21, 30)]})

    def test_transcript_dedup_combines_bam_and_seed_candidates(self):
        history_gtf = _write_gtf(
            self.tmp_path / "history.gtf",
            [
                'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Old.g2\tPan\ttranscript\t1\t10\t.\t+\t.\tgene_id "Old.g2"; transcript_id "Old.g2.1";',
                'Old.g2\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Old.g2"; transcript_id "Old.g2.1";',
            ],
        )
        bam_path = _write_bam(
            self.tmp_path / "new.bam",
            [("JM22.g1.1", 0, [(0, 8), (3, 5), (0, 12)])],
        )
        output_gtf = self.tmp_path / "output.gtf"

        with self.assertLogs(
            "pantrans.transcript_processor", level=logging.INFO
        ) as captured:
            transcript_dedup(
                bam_path,
                cluster_dic={"Ref.g1": ["Ref.g1", "JM22.g1"]},
                trans_len_dic={"Ref.g1.1": 20, "JM22.g1.1": 20},
                gene_len_dic={"Ref.g1": 100},
                gene_strand_dic={"Ref.g1": "+"},
                rename_map=None,
                gtf_path=str(output_gtf),
                seed_gtf_path=history_gtf,
            )

        transcripts, exons = _parse_output_gtf(output_gtf)
        self.assertEqual(len(transcripts), 2)
        self.assertEqual(
            sorted(exons.values()),
            [[(1, 8), (14, 25)], [(1, 10), (21, 30)]],
        )
        self.assertIn(
            "Skipped 1 historical GTF transcript",
            "\n".join(captured.output),
        )

    def test_transcript_dedup_prefers_seed_model_on_transcript_id_collision(self):
        history_gtf = _write_gtf(
            self.tmp_path / "history.gtf",
            [
                'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
                'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            ],
        )
        bam_path = _write_bam(
            self.tmp_path / "collision.bam",
            [("Ref.g1.1", 0, [(0, 8), (3, 5), (0, 12)])],
        )
        output_gtf = self.tmp_path / "output.gtf"

        transcript_dedup(
            bam_path,
            cluster_dic={"Ref.g1": ["Ref.g1"]},
            trans_len_dic={"Ref.g1.1": 20},
            gene_len_dic={"Ref.g1": 100},
            gene_strand_dic={"Ref.g1": "+"},
            rename_map=None,
            gtf_path=str(output_gtf),
            seed_gtf_path=history_gtf,
        )

        transcripts, exons = _parse_output_gtf(output_gtf)
        self.assertEqual(transcripts, ["Ref.g1.1"])
        self.assertEqual(exons, {"Ref.g1.1": [(1, 10), (21, 30)]})

    def test_construct_outputs_chain_into_append_historical_seed(self):
        gdna_path = self.tmp_path / "merged.gdna.fa"
        bed_path = self.tmp_path / "merged.bed"
        gdna_path.write_text(
            ">Ref.g1\n" + "A" * 50 + "\n>JM22.g1\n" + "C" * 50 + "\n"
        )
        bed_path.write_text(
            "chr1\t0\t50\tRef.g1\t0\t+\n"
            "chr1\t50\t100\tJM22.g1\t0\t+\n"
        )

        construct_bam = _write_bam(
            self.tmp_path / "construct.bam",
            [("Ref.g1.1", 0, [(0, 10), (3, 10), (0, 10)])],
        )
        construct_dir = self.tmp_path / "construct"
        construct_dir.mkdir()
        _write_cluster_outputs(
            cluster_dic={"Ref.g1": ["Ref.g1"]},
            all_gdna_path=str(gdna_path),
            all_bed_path=str(bed_path),
            filter_bam_path=construct_bam,
            trans_len_dic={"Ref.g1.1": 20},
            gene_len_dic={"Ref.g1": 50},
            gene_strand_dic={"Ref.g1": "+"},
            rename_map={"Ref.g1": "Pan1A000001"},
            out_dir=str(construct_dir),
            prefix="Pan",
            label="",
        )

        history_gtf = construct_dir / "Pan_unrenamed.gtf"
        history_cdna = construct_dir / "Pan_unrenamed_cdna.refer.fasta"
        self.assertTrue(history_gtf.is_file())
        self.assertTrue(history_cdna.is_file())

        append_bam = _write_bam(
            self.tmp_path / "append.bam",
            [("JM22.g1.1", 0, [(0, 8), (3, 5), (0, 12)])],
        )
        append_dir = self.tmp_path / "append"
        append_dir.mkdir()
        _write_cluster_outputs(
            cluster_dic={"Ref.g1": ["Ref.g1", "JM22.g1"]},
            all_gdna_path=str(gdna_path),
            all_bed_path=str(bed_path),
            filter_bam_path=append_bam,
            trans_len_dic={"Ref.g1.1": 20, "JM22.g1.1": 20},
            gene_len_dic={"Ref.g1": 50},
            gene_strand_dic={"Ref.g1": "+"},
            rename_map={"Ref.g1": "Pan1A000001"},
            out_dir=str(append_dir),
            prefix="Append",
            label="",
            seed_gtf_path=str(history_gtf),
        )

        append_unrenamed_gtf = append_dir / "Append_unrenamed.gtf"
        append_official_gtf = append_dir / "Append.gtf"
        unrenamed_transcripts, unrenamed_exons = _parse_output_gtf(
            append_unrenamed_gtf
        )
        official_transcripts, official_exons = _parse_output_gtf(
            append_official_gtf
        )
        unrenamed_fasta_ids = {
            line[1:].split()[0]
            for line in (append_dir / "Append_unrenamed_cdna.refer.fasta")
            .read_text()
            .splitlines()
            if line.startswith(">")
        }
        official_fasta_ids = {
            line[1:].split()[0]
            for line in (append_dir / "Append_cdna.refer.fasta")
            .read_text()
            .splitlines()
            if line.startswith(">")
        }

        self.assertEqual(set(unrenamed_transcripts), unrenamed_fasta_ids)
        self.assertEqual(set(official_transcripts), official_fasta_ids)
        self.assertEqual(len(unrenamed_transcripts), 2)
        self.assertEqual(len(official_transcripts), 2)
        self.assertTrue(
            all(transcript.startswith("Ref.g1.") for transcript in unrenamed_transcripts)
        )
        self.assertTrue(
            all(
                transcript.startswith("Pan1A000001.")
                for transcript in official_transcripts
            )
        )
        self.assertEqual(
            sorted(unrenamed_exons.values()),
            [[(1, 8), (14, 25)], [(1, 10), (21, 30)]],
        )
        self.assertEqual(
            sorted(official_exons.values()),
            sorted(unrenamed_exons.values()),
        )


if __name__ == "__main__":
    unittest.main()
