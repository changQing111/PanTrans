import copy
import itertools
import unittest

from pantrans.transcript_redundancy import deduplicate_models as dedup


def model(source, exons, cds=None, complete=True, **kwargs):
    return dict(source=source, exons=exons, cds=cds if cds is not None else [(200, 300)],
                complete=complete, pid=1., origin="test", **kwargs)


class TranscriptRedundancyTests(unittest.TestCase):
    def test_exact_chain_retains_longest_complete(self):
        a = model("a", [(0, 400)])
        b = model("b", [(20, 350)])
        out = dedup([a, b])
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["source"], "a")
        self.assertEqual(out[0]["sources"], ["a", "b"])

    def test_D_unique_with_extra_utr_exon(self):
        short = model("cds", [(200, 300)])
        full = model("full", [(0, 50), (100, 400)])
        out = dedup([short, full])
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["source"], "full")
        self.assertEqual(out[0]["sources"], ["cds", "full"])
        self.assertEqual(out[0]["source_relations"][0]["relation"], "unique_partial_compatible")

    def test_D_ambiguous_does_not_bridge_utr_isoforms(self):
        candidates = [model("cds", [(200, 300)]), model("f1", [(0, 50), (100, 400)]),
                      model("f2", [(0, 60), (150, 400)])]
        out = dedup(candidates)
        self.assertEqual(len(out), 3)
        short = next(m for m in out if m["source"] == "cds")
        self.assertEqual(short["redundancy_status"], "ambiguous")
        self.assertEqual(short["compatible_sources"], ["f1", "f2"])
        self.assertTrue(all("cds" not in m["sources"] for m in out if m["source"] != "cds"))

    def test_D_no_complete_kept(self):
        out = dedup([model("cds", [(200, 300)])])
        self.assertEqual(out[0]["utr_evidence"], "cds_only")
        self.assertEqual(out[0]["redundancy_status"], "unresolved")

    def test_E_missing_terminal_exon_unique_fold(self):
        full = model("f", [(0, 50), (100, 400)])
        short = model("p", [(150, 400)], complete=False, missing_5prime=True)
        self.assertEqual([m["source"] for m in dedup([full, short])], ["f"])

    def test_E_missing_three_prime_fold(self):
        full = model("f", [(0, 350), (400, 500)])
        short = model("p", [(0, 320)], complete=False, missing_3prime=True)
        self.assertEqual(len(dedup([full, short])), 1)

    def test_E_covered_utr_intron_retention_preserved(self):
        full = model("f", [(0, 400)])
        short = model("p", [(20, 70), (165, 400)], complete=False, missing_5prime=True)
        self.assertEqual(len(dedup([full, short])), 2)

    def test_E_covered_region_cannot_bridge_target_intron(self):
        full = model("f", [(0, 50), (100, 400)])
        short = model("p", [(20, 400)], complete=False, missing_5prime=True)
        self.assertEqual(len(dedup([full, short])), 2)

    def test_E_known_end_extension_rejected(self):
        full = model("f", [(0, 500)])
        short = model("p", [(20, 400)], complete=False, missing_5prime=True)
        self.assertEqual(len(dedup([full, short])), 2)

    def test_E_legacy_no_clipping_side_evidence_is_conservative(self):
        self.assertEqual(len(dedup([model("f", [(0, 400)]),
                                  model("p", [(20, 350)], complete=False)])), 2)

    def test_exact_coordinates_complete_and_partial_can_fold(self):
        self.assertEqual(len(dedup([model("f", [(0, 400)]),
                                  model("p", [(0, 400)], complete=False)])), 1)

    def test_different_CDS_remains_distinct(self):
        self.assertEqual(len(dedup([model("f", [(0, 400)]),
                                  model("cds", [(210, 300)], cds=[(210, 300)])])), 2)

    def test_no_CDS_mono_preserves_intervals(self):
        self.assertEqual(len(dedup([model("a", [(0, 100)], cds=[]),
                                  model("b", [(10, 90)], cds=[])])), 2)

    def test_complete_real_utr_splice_difference_preserved(self):
        self.assertEqual(len(dedup([model("a", [(0, 400)]),
                                  model("b", [(0, 50), (100, 400)])])), 2)

    def test_order_independent_and_no_mutation_with_json_lists(self):
        models = [model("cds", [[200, 300]]), model("a", [[0, 50], [100, 400]]),
                  model("b", [[0, 60], [150, 400]])]
        before = copy.deepcopy(models)
        reference = dedup(models)
        for order in itertools.permutations(models):
            self.assertEqual(dedup(order), reference)
        self.assertEqual(models, before)


if __name__ == "__main__":
    unittest.main()
