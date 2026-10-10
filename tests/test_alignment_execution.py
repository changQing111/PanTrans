import os
import shutil
import random
import subprocess
import tempfile
import unittest
from unittest import mock

import pysam

from pantrans import align_filter


class AlignmentExecutionTests(unittest.TestCase):
    def run_mock_mapping(self, directory, threads=8, map_status=0, view_status=0):
        calls = []

        def launch(command, **kwargs):
            calls.append(command)
            process = mock.Mock()
            status = map_status if command[0] == "minimap2" else view_status
            process.wait.return_value = status
            process.poll.return_value = status
            if command[0] == "minimap2":
                prefix = command[command.index("--split-prefix") + 1]
                self.assertTrue(os.path.isdir(os.path.dirname(prefix)))
            else:
                kwargs["stdout"].write(b"new BAM")
            return process

        output = os.path.join(directory, "output.bam")
        with mock.patch.object(align_filter.subprocess, "Popen", side_effect=launch):
            align_filter.minimap2_map("query with spaces.fa", "target.fa", threads, output)
        return calls, output

    def test_thread_allocation_and_unique_cleaned_directories(self):
        prefixes = set()
        with tempfile.TemporaryDirectory() as directory:
            for budget, mapping, extra in [(1, 1, 0), (2, 1, 0), (8, 6, 1), (16, 12, 3), (32, 28, 3)]:
                calls, output = self.run_mock_mapping(directory, budget)
                mapper, converter = calls
                self.assertEqual(mapper[mapper.index("-t") + 1], str(mapping))
                self.assertEqual(converter[converter.index("-@") + 1], str(extra))
                self.assertEqual(mapper[-1], "query with spaces.fa")
                prefix = mapper[mapper.index("--split-prefix") + 1]
                self.assertNotIn(prefix, prefixes)
                prefixes.add(prefix)
                self.assertFalse(os.path.exists(os.path.dirname(prefix)))
                with open(output, "rb") as handle:
                    self.assertEqual(handle.read(), b"new BAM")

    def test_failure_preserves_destination_and_cleans_temporary_files(self):
        for map_status, view_status in [(1, 0), (0, 1), (1, 1)]:
            with self.subTest(map_status=map_status, view_status=view_status):
                with tempfile.TemporaryDirectory() as directory:
                    output = os.path.join(directory, "output.bam")
                    with open(output, "wb") as handle:
                        handle.write(b"previous BAM")
                    with self.assertRaises(subprocess.CalledProcessError):
                        self.run_mock_mapping(directory, map_status=map_status, view_status=view_status)
                    self.assertEqual(os.listdir(directory), ["output.bam"])
                    with open(output, "rb") as handle:
                        self.assertEqual(handle.read(), b"previous BAM")

    def test_split_prefix_is_normalized_but_strategy_is_preserved(self):
        expected = align_filter.MINIMAP2_OPTIONS + align_filter.SPLIT_OPTIONS
        for prefix in ["--split-prefix '/tmp/a b/mm2'", "--split-prefix=/tmp/other/mm2"]:
            command = "minimap2 " + align_filter.MINIMAP2_OPTIONS + " " + prefix + " -t 6 target.fa query.fa"
            self.assertEqual(align_filter._minimap2_options_from_command(command), expected)
        plain = "minimap2 " + align_filter.MINIMAP2_OPTIONS + " target.fa query.fa"
        self.assertNotEqual(align_filter._minimap2_options_from_command(plain), expected)

    def test_invalid_threads_fail_before_launch(self):
        with mock.patch.object(align_filter.subprocess, "Popen") as launch:
            with self.assertRaises(ValueError):
                align_filter.minimap2_map("q.fa", "t.fa", 0, "out.bam")
            launch.assert_not_called()

    def test_converter_launch_failure_reaps_mapper(self):
        mapper = mock.Mock()
        mapper.poll.return_value = None
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(align_filter.subprocess, "Popen", side_effect=[mapper, FileNotFoundError()]):
                with self.assertRaises(FileNotFoundError):
                    align_filter.minimap2_map("q.fa", "t.fa", 4, os.path.join(directory, "out.bam"))
            mapper.kill.assert_called_once()
            mapper.wait.assert_called_once()
            self.assertEqual(os.listdir(directory), [])

    @unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "alignment tools unavailable")
    def test_real_split_index_has_complete_header_and_queries(self):
        rng = random.Random(42)
        sequences = {"gene%d" % i: "".join(rng.choices("ACGT", k=2000)) for i in range(4)}
        with tempfile.TemporaryDirectory() as directory:
            target = os.path.join(directory, "target.fa")
            query = os.path.join(directory, "query.fa")
            output = os.path.join(directory, "output.bam")
            with open(target, "w") as handle:
                for name, seq in sequences.items():
                    handle.write(">%s\n%s\n" % (name, seq))
            with open(query, "w") as handle:
                for name, seq in sequences.items():
                    handle.write(">%s.1\n%s\n" % (name, seq))
            options = align_filter.MINIMAP2_OPTIONS.replace("-I 2G", "-I 1k")
            with mock.patch.object(align_filter, "MINIMAP2_OPTIONS", options):
                align_filter.minimap2_map(query, target, 4, output)
                with pysam.AlignmentFile(output, "rb") as bam:
                    self.assertEqual(dict(zip(bam.references, bam.lengths)), {name: 2000 for name in sequences})
                    primaries = [read for read in bam if not read.is_secondary and not read.is_supplementary]
                    self.assertEqual(len(primaries), 4)
                    self.assertEqual({read.query_name for read in primaries}, {name + ".1" for name in sequences})
                    self.assertTrue(all(not read.is_unmapped for read in primaries))
                provenance = align_filter.bam_alignment_provenance(output)
                self.assertEqual(provenance["minimap2_options"], options + align_filter.SPLIT_OPTIONS)
            self.assertFalse(any(name.startswith(".pantrans-align-") for name in os.listdir(directory)))


if __name__ == "__main__":
    unittest.main()
