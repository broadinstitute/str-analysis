import os
import unittest

from str_analysis.utils.trf_runner import TRFRunner

class Tests(unittest.TestCase):
    def test_trf_runner(self):

        trf_path = os.path.expanduser("~/bin/trf409.macosx")
        if not os.path.exists(trf_path):
            print("WARNING: TRF binary not found at ~/bin/trf409.macosx")
            return

        trf_runner = TRFRunner(trf_path, mismatch_penalty=3, indel_penalty=5, minscore=8, html_mode=False)

        repeat_sequence = "ACACACACAC"
        dat_records = list(trf_runner.run_TRF_on_nucleotide_sequence(repeat_sequence))
        self.assertEqual(len(dat_records), 1)
        self.assertEqual(dat_records[0]["repeat_unit"], "AC")
        self.assertEqual(dat_records[0]["repeat_count"], 5)

        repeat_sequence = "ACGAACGAACGAACGAATGA"
        dat_records = list(trf_runner.run_TRF_on_nucleotide_sequence(repeat_sequence))
        self.assertEqual(len(dat_records), 1)
        self.assertEqual(dat_records[0]["repeat_unit"], "ACGA")
        self.assertEqual(dat_records[0]["repeat_count"], 5)
