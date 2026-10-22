import copy
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from pantrans.gff3 import (make_row, write_gff3, read_gff3, rename_gff3_ids,
                          load_gff3_transcript_models, load_history_transcript_models,
                          parse_attributes, format_attributes)
from pantrans.transcript_processor import generate_gff3, get_cdna_from_gff3
from pantrans.pipeline import (_append_rescue_records_from_gff3,
                               _extract_gff3_records_by_gene, _append_file_contents)
from pantrans.main import read_parameters


class Gff3Tests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root / "model.gff3"

    def rows(self, gene="g", strand="+"):
        tid = gene + ".12"
        return [make_row(gene, "gene", 1, 18, strand, gene),
                make_row(gene, "mRNA", 1, 18, strand, tid, gene),
                make_row(gene, "exon", 1, 6, strand, tid + ".exon1", tid),
                make_row(gene, "exon", 13, 18, strand, tid + ".exon2", tid)]

    def test_header_hierarchy_and_numeric_suffix(self):
        rows = self.rows()
        write_gff3(self.path, list(reversed(rows)))
        self.assertTrue(self.path.read_text().startswith("##gff-version 3\n"))
        models = load_gff3_transcript_models(self.path)
        self.assertEqual(models["transcript_gene"], {"g.12": "g"})
        self.assertEqual(models["splice_sites"], {"g.12": [(7, 12)]})
        self.assertEqual(read_gff3(self.path), rows)

    def test_transcripts_are_followed_by_all_their_children(self):
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                rows = self.rows(strand=strand)
                for suffix in ("2", "10"):
                    tid = "g." + suffix
                    rows.extend([
                        make_row("g", "mRNA", 1, 18, strand, tid, "g"),
                        make_row("g", "exon", 1, 6, strand, tid + ".exon1", tid),
                        make_row("g", "exon", 13, 18, strand, tid + ".exon2", tid),
                        make_row("g", "CDS", 3, 6, strand, tid + ".cds", tid, phase=0),
                        make_row("g", "CDS", 13, 17, strand, tid + ".cds", tid, phase=2),
                        make_row("g", "five_prime_UTR", 1, 2, strand, tid + ".utr", tid),
                    ])
                write_gff3(self.path, list(reversed(rows)))
                result = read_gff3(self.path)
                self.assertEqual(result[0][2], "gene")
                active = None
                order = []
                for row in result[1:]:
                    if row[2] == "mRNA":
                        active = row[8]["ID"]
                        order.append(active)
                    else:
                        self.assertEqual(row[8]["Parent"], [active])
                self.assertEqual(order, ["g.2", "g.10", "g.12"])
                self.assertCountEqual(result, rows)
                original = self.path.read_bytes()
                write_gff3(self.path, result)
                self.assertEqual(self.path.read_bytes(), original)

    def test_id_escaping_round_trip_and_rename(self):
        gid = "g;=,%20&"
        rows = self.rows(gid)
        rows[1][8]["note"] = "retain ;,=% text"
        write_gff3(self.path, rows)
        self.assertIn("%3B", self.path.read_text())
        self.assertEqual(read_gff3(self.path), rows)
        renamed = self.root / "official.gff3"
        rename_gff3_ids(self.path, renamed, {gid: "PanChr1000001"})
        result = read_gff3(renamed)
        self.assertEqual([r[:8] for r in result], [r[:8] for r in rows])
        self.assertEqual(result[1][8]["ID"], "PanChr1000001.12")
        self.assertEqual(result[1][8]["Parent"], ["PanChr1000001"])
        self.assertEqual(result[1][8]["note"], rows[1][8]["note"])
        self.assertTrue(all(r[8]["ID"].startswith("PanChr1000001") for r in result))

    def test_parent_lists_split_before_unescaping(self):
        attrs = {"ID": "x", "Parent": ["a,b", "c;d"]}
        self.assertEqual(parse_attributes(format_attributes(attrs)), attrs)

    def test_rename_rejects_invalid_suffix_and_missing_map(self):
        rows = self.rows()
        rows[1][8]["ID"] = "g.isoform"
        for row in rows[2:]:
            row[8]["Parent"] = ["g.isoform"]
        write_gff3(self.path, rows)
        with self.assertRaisesRegex(ValueError, "numeric suffix"):
            rename_gff3_ids(self.path, self.root / "out.gff3", {"g": "Pan1"})
        write_gff3(self.path, self.rows())
        with self.assertRaises(KeyError):
            rename_gff3_ids(self.path, self.root / "out.gff3", {})

    def test_invalid_models_do_not_replace_existing_output(self):
        good = self.rows()
        write_gff3(self.path, good)
        original = self.path.read_bytes()
        cases = []
        bad = copy.deepcopy(good); bad[2][8]["Parent"] = ["absent"]; cases.append(bad)
        bad = copy.deepcopy(good); bad[2][8]["ID"] = "g"; cases.append(bad)
        bad = copy.deepcopy(good); bad[2][6] = "-"; cases.append(bad)
        bad = copy.deepcopy(good); bad[3][3] = "5"; cases.append(bad)
        bad = copy.deepcopy(good); bad[3][4] = "25"; cases.append(bad)
        bad = copy.deepcopy(good); bad[1][8]["Parent"] = ["absent"]; cases.append(bad)
        bad = copy.deepcopy(good); bad[2][3] = "0"; cases.append(bad)
        cases.append(good + [make_row("g", "CDS", 3, 14, "+", "g.12.cds", "g.12", phase=0)])
        for bad in cases:
            with self.subTest(rows=bad), self.assertRaises(ValueError):
                write_gff3(self.path, bad)
            self.assertEqual(self.path.read_bytes(), original)

    def test_discontinuous_cds_and_renamed_parents(self):
        rows = self.rows() + [make_row("g", "CDS", 3, 6, "+", "g.12.cds", "g.12", phase=0),
                              make_row("g", "CDS", 13, 17, "+", "g.12.cds", "g.12", phase=2)]
        write_gff3(self.path, rows)
        renamed = self.root / "official.gff3"
        rename_gff3_ids(self.path, renamed, {"g": "Pan1"})
        cds = [r for r in read_gff3(renamed) if r[2] == "CDS"]
        self.assertEqual([r[7] for r in cds], ["0", "2"])
        self.assertEqual([r[8] for r in cds], [{"ID": "Pan1.12.cds", "Parent": ["Pan1.12"]}] * 2)

    def test_extract_cdna_both_strands_and_official_ids(self):
        fasta = self.root / "genes.fa"
        fasta.write_text(">g\nATGAAANNNNNNCCCTAA\n")
        for strand, expected in (("+", "ATGAAACCCTAA"), ("-", "TTAGGGTTTCAT")):
            write_gff3(self.path, self.rows(strand=strand))
            renamed = self.root / "official.gff3"
            rename_gff3_ids(self.path, renamed, {"g": "Pan1"})
            cdna = self.root / "cdna.fa"
            get_cdna_from_gff3(str(fasta), renamed, cdna)
            records = list(SeqIO.parse(cdna, "fasta"))
            self.assertEqual(records[0].id, "Pan1.12")
            self.assertEqual(str(records[0].seq), expected)

    def test_shared_exon_parent_is_extracted_for_both_transcripts(self):
        rows = self.rows()
        rows.insert(2, make_row("g", "mRNA", 1, 18, "+", "g.13", "g"))
        for row in rows:
            if row[2] == "exon":
                row[8]["Parent"] = ["g.12", "g.13"]
        write_gff3(self.path, rows)
        models = load_gff3_transcript_models(self.path)
        self.assertEqual(models["exon_coords"]["g.12"], models["exon_coords"]["g.13"])

    def test_rescue_merge_writes_single_header(self):
        pre = self.root / "pre.gff3"
        write_gff3(pre, self.rows("g") + self.rows("h"))
        write_gff3(self.path, self.rows("g"))
        records = _extract_gff3_records_by_gene(pre, ["h"])
        self.assertEqual(_append_rescue_records_from_gff3(self.path, records, None), {"h"})
        self.assertEqual(self.path.read_text().count("##gff-version"), 1)
        self.assertEqual(set(load_gff3_transcript_models(self.path)["gene_rows"]), {"g", "h"})
        rescue = self.root / "rescue.gff3"
        write_gff3(rescue, self.rows("i"))
        _append_file_contents(rescue, self.path)
        self.assertEqual(self.path.read_text().count("##gff-version"), 1)

    def test_legacy_history_and_cli_alias(self):
        path = self.root / "history.gtf"
        attrs = 'gene_id "g"; transcript_id "g.1";'
        path.write_text(f"g\tPan\ttranscript\t1\t6\t.\t+\t.\t{attrs}\n"
                        f"g\tPan\texon\t1\t6\t.\t+\t.\t{attrs}\n")
        self.assertEqual(load_history_transcript_models(path)["exon_coords"], {"g.1": [(1, 6)]})
        args = ["append", "-c", "c", "-g", "g", "-b", "b", "--query-cdna", "q",
                "--history-graph", "graph", "-o", "out"]
        for option in ("--history-gff3", "--history-gtf"):
            self.assertEqual(read_parameters(args + [option, str(path)]).history_gff3, str(path))

    def test_empty_file_has_valid_header(self):
        generate_gff3({}, {}, {}, {}, {}, None, self.path)
        self.assertEqual(self.path.read_text(), "##gff-version 3\n")
        self.assertEqual(load_gff3_transcript_models(self.path)["transcript_gene"], {})

    def test_rejects_gtf_disguised_as_gff3(self):
        self.path.write_text('g\tPan\texon\t1\t6\t.\t+\t.\tgene_id "g"; transcript_id "g.1";\n')
        with self.assertRaises(ValueError):
            read_gff3(self.path)


if __name__ == "__main__":
    unittest.main()
