import logging

import pysam
import pytest

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


def test_load_gtf_transcript_models_sorts_exons_and_derives_splice_sites(tmp_path):
    history_gtf = _write_gtf(
        tmp_path / "history.gtf",
        [
            'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
        ],
    )

    models = load_gtf_transcript_models(history_gtf)

    assert models["transcript_gene"] == {"Ref.g1.1": "Ref.g1"}
    assert models["exon_coords"] == {"Ref.g1.1": [(1, 10), (21, 30)]}
    assert models["splice_sites"] == {"Ref.g1.1": [(11, 20)]}


def test_load_gtf_transcript_models_requires_transcript_feature(tmp_path):
    history_gtf = _write_gtf(
        tmp_path / "exon_only.gtf",
        [
            'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
        ],
    )

    with pytest.raises(ValueError, match="no transcript feature"):
        load_gtf_transcript_models(history_gtf)


def test_load_gtf_transcript_models_rejects_transcript_gene_id_mismatch(tmp_path):
    history_gtf = _write_gtf(
        tmp_path / "mismatched_transcript_gene.gtf",
        [
            'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Other.g1.1";',
            'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Other.g1.1";',
            'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Other.g1.1";',
        ],
    )

    with pytest.raises(ValueError, match="transcript_id.*gene_id"):
        load_gtf_transcript_models(history_gtf)


@pytest.mark.parametrize(
    "lines, message",
    [
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
    ],
)
def test_load_gtf_transcript_models_rejects_malformed_models(
    tmp_path, lines, message
):
    history_gtf = _write_gtf(tmp_path / "malformed.gtf", lines)

    with pytest.raises(ValueError, match=message):
        load_gtf_transcript_models(history_gtf)


def test_transcript_dedup_preserves_seed_transcript_absent_from_bam(tmp_path):
    history_gtf = _write_gtf(
        tmp_path / "history.gtf",
        [
            'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
        ],
    )
    bam_path = _write_bam(tmp_path / "empty.bam")
    output_gtf = tmp_path / "output.gtf"

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
    assert transcripts == ["Ref.g1.1"]
    assert exons == {"Ref.g1.1": [(1, 10), (21, 30)]}


def test_transcript_dedup_combines_bam_and_seed_candidates(tmp_path, caplog):
    history_gtf = _write_gtf(
        tmp_path / "history.gtf",
        [
            'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Old.g2\tPan\ttranscript\t1\t10\t.\t+\t.\tgene_id "Old.g2"; transcript_id "Old.g2.1";',
            'Old.g2\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Old.g2"; transcript_id "Old.g2.1";',
        ],
    )
    bam_path = _write_bam(
        tmp_path / "new.bam",
        [("JM22.g1.1", 0, [(0, 8), (3, 5), (0, 12)])],
    )
    output_gtf = tmp_path / "output.gtf"

    with caplog.at_level(logging.INFO, logger="pantrans.transcript_processor"):
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
    assert len(transcripts) == 2
    assert sorted(exons.values()) == [
        [(1, 8), (14, 25)],
        [(1, 10), (21, 30)],
    ]
    assert "Skipped 1 historical GTF transcript" in caplog.text


def test_transcript_dedup_prefers_seed_model_on_transcript_id_collision(tmp_path):
    history_gtf = _write_gtf(
        tmp_path / "history.gtf",
        [
            'Ref.g1\tPan\ttranscript\t1\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t1\t10\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
            'Ref.g1\tPan\texon\t21\t30\t.\t+\t.\tgene_id "Ref.g1"; transcript_id "Ref.g1.1";',
        ],
    )
    bam_path = _write_bam(
        tmp_path / "collision.bam",
        [("Ref.g1.1", 0, [(0, 8), (3, 5), (0, 12)])],
    )
    output_gtf = tmp_path / "output.gtf"

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
    assert transcripts == ["Ref.g1.1"]
    assert exons == {"Ref.g1.1": [(1, 10), (21, 30)]}
