import pandas as pd
import unittest

from generate_gnomad_json import compute_most_common_motif_lookup_dict


class Tests(unittest.TestCase):

    def setUp(self):
        self._df = pd.DataFrame([
            {"LocusId": "RFC1", "Motif: Allele 1": "AAAAG", "Motif: Allele 2": "GAAAG"},
            {"LocusId": "RFC1", "Motif: Allele 1": "AAAGA", "Motif: Allele 2": "AAAGG"},
            {"LocusId": "RFC1", "Motif: Allele 1": "AAAAG", "Motif: Allele 2": "GAAAG"},
        ])

    def test_compute_most_common_motif_lookup_dict(self):
        most_common_motif_lookup = compute_most_common_motif_lookup_dict(self._df)

        self.assertTrue(("RFC1", "AAAAG") in most_common_motif_lookup)
        self.assertTrue(("RFC1", "AAAGG") in most_common_motif_lookup)
        self.assertEqual(most_common_motif_lookup["RFC1", "AAAAG"], "AAAAG")
        self.assertEqual(most_common_motif_lookup["RFC1", "AAAGG"], "GAAAG")
