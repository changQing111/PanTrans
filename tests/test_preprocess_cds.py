import csv
import json
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from pantrans.preprocess import preprocess
from pantrans.prepared_bundle import load_prepared_bundle, prepared_file, manifest
from pantrans.cds_preparation import coding_model
from pantrans.cds_merge import merge_prepared


class PreprocessCDSTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.sequence = "CC" + "ATGA" + "GGG" + "AACCCTAA" + "TT"
        self.genome = self.root / "genome.fa"
        self.genome.write_text(">p\n" + self.sequence + "\n>n\n" +
                               str(Seq(self.sequence).reverse_complement()) + "\n")
        lines = []
        for chrom, gene, strand, pieces in [("p", "positive", "+", [(3,6,0),(10,17,2)]),
                                            ("n", "negative", "-", [(14,17,0),(3,10,2)])]:
            lines += [f"{chrom}\tx\tgene\t3\t17\t.\t{strand}\t.\tID={gene}",
                      f"{chrom}\tx\tmRNA\t3\t17\t.\t{strand}\t.\tID={gene}.1;Parent={gene}"]
            for a,b,phase in pieces:
                for feature in ("exon", "CDS"):
                    lines.append(f"{chrom}\tx\t{feature}\t{a}\t{b}\t.\t{strand}\t{phase if feature=='CDS' else '.'}\tParent={gene}.1")
        self.gff = self.root / "models.gff3"
        self.gff.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")

    def prepare(self):
        return preprocess([self.genome], [self.gff], ["Var"], self.root / "prepared", prefix="custom", flank=3000)

    def sequences(self, path):
        with open(path) as handle:
            return {r.id: str(r.seq) for r in SeqIO.parse(handle, "fasta")}

    def test_positive_negative_split_phase_and_body(self):
        paths = self.prepare()
        expected = {"positive.1": "ATGAAACCCTAA", "negative.1": "ATGAAACCCTAA"}
        self.assertEqual(self.sequences(paths["cds"]), expected)
        self.assertEqual(self.sequences(paths["cds_valid"]), expected)
        self.assertEqual(self.sequences(paths["cds_longest"]), expected)
        self.assertTrue(Path(paths["cds_longest"] + ".fai").exists())
        self.assertEqual(self.sequences(paths["gene"]), {"positive": self.sequence[2:17], "negative": self.sequence[2:17]})
        self.assertEqual(self.sequences(paths["cdna"]), expected)
        self.assertTrue(Path(paths["cds_valid"] + ".fai").exists())
        import csv
        with open(paths["cds_regions"]) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual([r["phase"] for r in rows], ["0", "2", "0", "2"])
        for row in rows:
            seq = self.sequences(paths["gene"])[row["gene"]]
            self.assertEqual(seq[int(row["gene_start0"]):int(row["gene_end0"])],
                             self.sequence[int(row["gdna_start0"]):int(row["gdna_end0"])])

    def test_schema2_consumer_without_original_annotation_or_fixed_prefix(self):
        paths = self.prepare()
        self.gff.unlink()  # Consumption must not depend on re-reading original GFF3.
        out = self.root / "bundle"; out.mkdir()
        regions, models, tids = load_prepared_bundle(self.root / "prepared", out)
        self.assertEqual(len(regions), 2)
        self.assertTrue(all(m["valid"] for m in models.values()))
        self.assertEqual(models["negative.1"]["qcds"], [(0,4),(4,12)])
        self.assertEqual(self.sequences(out / "cds.fa"), self.sequences(paths["cds_valid"]))
        self.assertEqual(prepared_file(self.root / "prepared", "gdna"), Path(paths["gdna"]))

    def test_invalid_phase_is_reported_not_silently_dropped(self):
        self.gff.write_text(self.gff.read_text().replace("\t2\tParent=", "\t1\tParent="))
        paths = self.prepare()
        self.assertEqual(len(self.sequences(paths["cds"])), 2)
        self.assertEqual(Path(paths["cds_valid"]).read_text(), "")
        self.assertFalse(Path(paths["cds_valid"] + ".fai").exists())
        self.assertEqual(Path(paths["cds_longest"]).read_text(), "")
        self.assertFalse(Path(paths["cds_longest"] + ".fai").exists())
        with open(paths["cds_longest_selection"]) as f:
            self.assertEqual({r["status"] for r in csv.DictReader(f, delimiter="\t")}, {"no_valid_CDS"})
        data = json.loads(Path(paths["manifest"]).read_text())
        self.assertIn("cds_longest", data["empty_fastas"])
        self.assertEqual(data["longest_cds_selection"]["genes_without_valid_cds"], 2)
        self.assertEqual(len(self.sequences(paths["gene"])), 2)
        self.assertIn("invalid_phase", Path(paths["cds_status"]).read_text())
        out = self.root / "bundle"; out.mkdir()
        with self.assertRaisesRegex(ValueError, "No valid CDS"):
            load_prepared_bundle(self.root / "prepared", out)

    def test_no_CDS_keeps_gene_cdna_bed(self):
        self.gff.write_text("\n".join(l for l in self.gff.read_text().splitlines() if "\tCDS\t" not in l)+"\n")
        paths = self.prepare()
        self.assertEqual(len(self.sequences(paths["gene"])), 2)
        self.assertEqual(len(self.sequences(paths["cdna"])), 2)
        self.assertEqual(Path(paths["cds"]).read_text(), "")
        self.assertIn("no_CDS", Path(paths["cds_status"]).read_text())
        self.assertEqual(json.loads(Path(paths["manifest"]).read_text())["schema_version"], 2)

    def test_bad_CDS_structure_fails_without_partial_outputs(self):
        self.gff.write_text(self.gff.read_text().replace("CDS\t3\t6", "CDS\t3\t8"))
        with self.assertRaisesRegex(ValueError, "CDS outside"):
            self.prepare()
        self.assertFalse((self.root / "prepared/custom_gene.fa").exists())

    def test_ambiguous_stop_and_partial_ORF_classification(self):
        for seq, expected in [("ATGNNNTAA","ambiguous_CDS"), ("ATGTAATAA","internal_stop"),
                              ("CCCAAATAA","partial_or_invalid_ORF")]:
            _, status, model = coding_model(seq, [(0,len(seq))], [(0,len(seq),0)])
            self.assertEqual(status, expected)
            self.assertFalse(model["valid"])

    def test_selection_filters_both_graph_fastas(self):
        self.prepare();out=self.root/"bundle";out.mkdir()
        r,m,t=load_prepared_bundle(self.root/"prepared",out,{"negative"})
        self.assertEqual(set(r), {"negative"})
        self.assertEqual(set(self.sequences(out/"cds.fa")), {"negative.1"})
        self.assertEqual(set(self.sequences(out/"body.fa")), {"negative"})

    def test_longest_valid_selection_ties_invalid_longer_and_namespace(self):
        # The full 15bp CDS has an internal stop; skipping it yields valid 12bp.
        self.genome.write_text(">p\nATGAAATAACCCTAA\n")
        lines = ["p\tx\tgene\t1\t15\t.\t+\t.\tID=g"]
        for suffix, pieces in [("2", [(1, 6), (10, 15)]),
                               ("10", [(1, 6), (10, 15)]),
                               ("1", [(1, 6), (13, 15)]), ("3", [(1, 15)])]:
            tid = "g." + suffix
            lines.append(f"p\tx\tmRNA\t1\t15\t.\t+\t.\tID={tid};Parent=g")
            for a, b in pieces:
                lines.append(f"p\tx\texon\t{a}\t{b}\t.\t+\t.\tParent={tid}")
                lines.append(f"p\tx\tCDS\t{a}\t{b}\t.\t+\t0\tParent={tid}")
        self.gff.write_text("\n".join(lines) + "\n")
        paths = preprocess([self.genome, self.genome], [self.gff, self.gff],
                           ["V1", "V2"], self.root / "longest", namespace_ids=True)
        self.assertEqual(self.sequences(paths["cds_longest"]),
                         {"V1__g.10": "ATGAAACCCTAA", "V2__g.10": "ATGAAACCCTAA"})
        self.assertEqual(len(self.sequences(paths["cdna"])), 8)
        self.assertEqual(len(self.sequences(paths["cds"])), 8)
        self.assertEqual(len(self.sequences(paths["cds_valid"])), 6)
        with open(paths["cds_longest_selection"]) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        for row, variety in zip(rows, ["V1", "V2"]):
            self.assertEqual(row["gene"], variety + "__g")
            self.assertEqual(row["transcript"], variety + "__g.10")
            self.assertEqual((row["CDS_length"], row["valid_transcripts"], row["max_length_ties"]),
                             ("12", "3", "2"))
            self.assertEqual((row["variety"], row["chrom"], row["strand"],
                              row["gene_start0"], row["gene_end0"]), (variety, "p", "+", "0", "15"))
        self.gff.write_text("\n".join(reversed(lines)) + "\n")
        again = preprocess([self.genome, self.genome], [self.gff, self.gff],
                           ["V1", "V2"], self.root / "reordered", namespace_ids=True)
        for key in ("cds_longest", "cds_longest_selection"):
            self.assertEqual(Path(paths[key]).read_bytes(), Path(again[key]).read_bytes())
        data = json.loads(Path(paths["manifest"]).read_text())
        self.assertEqual(data["longest_cds_selection"]["genes_with_valid_cds"], 2)
        self.assertEqual(data["longest_cds_selection"]["genes_with_max_cds_length_ties"], 2)

    def test_gene_without_transcripts_retains_dna_coordinates_and_selection_status(self):
        with self.gff.open("a") as f:
            f.write("p\tx\tgene\t1\t2\t.\t+\t.\tID=unmodeled\n")
        paths = self.prepare()
        for key in ("gene", "gdna"):
            self.assertIn("unmodeled", self.sequences(paths[key]))
        self.assertNotIn("unmodeled", {t.rsplit(".", 1)[0] for t in self.sequences(paths["cdna"])})
        with open(paths["cds_longest_selection"]) as f:
            row = next(r for r in csv.DictReader(f, delimiter="\t") if r["gene"] == "unmodeled")
        self.assertEqual((row["status"], row["transcript"], row["CDS_length"]), ("no_transcripts", "", "0"))
        data = json.loads(Path(paths["manifest"]).read_text())
        self.assertEqual(data["varieties"][0]["genes"], 3)
        self.assertEqual(data["varieties"][0]["genes_without_transcripts"], 1)
        self.assertEqual(data["varieties"][0]["skipped_genes_without_transcripts"], 0)

    def test_gene_only_annotation_and_new_output_overwrite_guard(self):
        self.gff.write_text("p\tx\tgene\t1\t2\t.\t+\t.\tID=unmodeled\n")
        paths = self.prepare()
        self.assertEqual(self.sequences(paths["gene"]), {"unmodeled": "CC"})
        self.assertEqual(Path(paths["cdna"]).read_text(), "")
        self.assertEqual(Path(paths["cds_longest"]).read_text(), "")
        out = self.root / "blocked"
        out.mkdir()
        blocker = out / "custom_cds.longest.valid.fa.fai"
        blocker.write_text("existing index\n")
        with self.assertRaises(FileExistsError):
            preprocess([self.genome], [self.gff], ["Var"], out, prefix="custom")
        self.assertEqual(list(out.iterdir()), [blocker])

    def test_merge_preserves_longest_artifacts_and_accepts_older_schema2(self):
        first, second = self.root / "first", self.root / "second"
        for directory, variety in [(first, "A"), (second, "B")]:
            preprocess([self.genome], [self.gff], [variety], directory, namespace_ids=True)
        merged = self.root / "merged"
        merge_prepared(first, second, merged)
        mf = json.loads((merged / "merged_preprocess.json").read_text())
        self.assertEqual(set(self.sequences(merged / mf["files"]["cds_longest"])),
                         {"A__positive.1", "A__negative.1", "B__positive.1", "B__negative.1"})
        self.assertEqual(mf["longest_cds_selection"]["genes_with_valid_cds"], 4)
        with (merged / mf["files"]["cds_longest_selection"]).open() as f:
            self.assertEqual(len(list(csv.DictReader(f, delimiter="\t"))), 4)
        # Simulate a historical manifest predating the optional outputs.
        path = next(second.glob("*_preprocess.json"))
        old = json.loads(path.read_text())
        for key in ("cds_longest", "cds_longest_selection"):
            old["files"].pop(key)
        old.pop("longest_cds_selection")
        path.write_text(json.dumps(old))
        manifest.cache_clear()
        mixed = self.root / "mixed"
        merge_prepared(first, second, mixed)
        mixed_manifest = json.loads((mixed / "merged_preprocess.json").read_text())
        self.assertNotIn("cds_longest", mixed_manifest["files"])
        self.assertNotIn("longest_cds_selection", mixed_manifest)


if __name__ == "__main__":
    unittest.main()
