import unittest
import pysam
from pantrans.cds_filter_experiment import compensated_indels,compatible_junctions

class StructuralPolicyTests(unittest.TestCase):
    def read(self,cigar):
        r=pysam.AlignedSegment();r.reference_start=0;r.cigarstring=cigar;return r
    def test_local_compensating_indels(self):
        self.assertTrue(compensated_indels(self.read('30M1I3M2I40M'),[(0,200)],[(0,200)]))
        self.assertTrue(compensated_indels(self.read('30M1I3M1D40M'),[(0,200)],[(0,200)]))
        self.assertFalse(compensated_indels(self.read('30M1I30M2I40M'),[(0,200)],[(0,200)]))
        self.assertFalse(compensated_indels(self.read('30M1I3M40N2I40M'),[(0,200)],[(0,200)]))
        self.assertFalse(compensated_indels(self.read('30M1I40M'),[(0,200)],[(0,200)]))
        self.assertTrue(compensated_indels(self.read('30M3I40M'),[(0,200)],[(0,200)]))
    def test_rejects_block_boundary_and_long_disruptive_indel(self):
        self.assertFalse(compensated_indels(self.read('30M1I3M2I40M'),[(0,200)],[(0,31),(32,200)]))
        self.assertFalse(compensated_indels(self.read('30M14I1M1I40M'),[(0,200)],[(0,200)]))
    def test_shifted_junctions_must_match_one_to_one(self):
        self.assertTrue(compatible_junctions({(101,201)},{(100,200)},0,300))
        self.assertTrue(compatible_junctions({(103,203)},{(100,200)},0,300))
        self.assertFalse(compatible_junctions({(104,204)},{(100,200)},0,300))
        self.assertFalse(compatible_junctions({(101,202)},{(100,200)},0,300))
        self.assertFalse(compatible_junctions(set(),{(100,200)},0,300))
        self.assertFalse(compatible_junctions({(100,200)},set(),0,300))
        self.assertFalse(compatible_junctions({(101,201)},{(100,200),(102,202)},0,300))
        self.assertFalse(compatible_junctions({(99,199),(101,201)},{(100,200)},0,300))
