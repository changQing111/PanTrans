"""Check the publishable restart entry against the current pipeline interface."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import pysam

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("restart", REPO / "scripts/restart_cds_annotation.py")
restart = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(restart)


class RestartEntryTests(unittest.TestCase):
    def test_repository_entry_union_cache_and_stage_outputs(self):
        pipeline = restart.load_pipeline(REPO)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = root / "prepared"
            prepared.mkdir()
            reference = prepared / "five_variety_gene3000bp_gdna.fa"
            reference.write_text(">hub\n" + "A" * 100 + "\n")
            pysam.faidx(str(reference))
            output = root / "results"
            output.mkdir()
            calls = []

            def fake_cluster(cluster, prepared, models, tids, work):
                calls.append(cluster)
                candidates = [dict(source=g + ".1", exons=[(0, 100)], cds=[(20, 80)],
                                   complete=True, origin="test", pid=1.0) for g in cluster]
                return cluster[0], candidates, [], []

            pipeline.annotation_cluster = fake_cluster
            clusters = {"pre": [["hub", "q"]], "last": [["hub", "q", "r"]]}
            tids = {g: [g + ".1"] for g in ("hub", "q", "r")}
            models = {t: {"valid": True} for ts in tids.values() for t in ts}
            handles, db, counts = restart.install(pipeline, clusters, models, tids, output)
            try:
                for stage in ("pre", "last"):
                    pipeline.annotation_stage(stage, clusters[stage], prepared, models, tids, output, 1)
            finally:
                handles.close()
                db.close()
            self.assertEqual(calls, [["hub", "q", "r"]])
            self.assertEqual(counts["cache_hits"], 1)
            for stage, expected in (("pre", "hub.1,q.1"), ("last", "hub.1,q.1,r.1")):
                with (output / (stage + ".model_sources.tsv")).open() as handle:
                    row = next(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(row["sources"], expected)
                self.assertTrue(json.loads(row["source_relations_json"]))
                self.assertEqual(row["redundancy_status"], "observed")


if __name__ == "__main__":
    unittest.main()
