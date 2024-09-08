import unittest

import pandas as pd
from str_analysis.filter_to_largest_expansions import add_allele_histogram


class Tests(unittest.TestCase):

	def test_add_allele_histogram(self):
		df = pd.DataFrame({
			"LocusId": ["locus1", "locus1", "locus1", "locus2", "locus2", "locus3"],
			"Num Repeats: Allele 1": [1, 2, 3, 4, 5, 3],
			"Num Repeats: Allele 2": [2, 3, 4, 5, 6, None],
		})

		add_allele_histogram(df, "AlleleHist", "AlleleStdev")

		self.assertEqual("1x:1,2x:2,3x:2,4x:1", df[df["LocusId"] == "locus1"].iloc[0]["AlleleHist"])
		self.assertEqual("4x:1,5x:2,6x:1", df[df["LocusId"] == "locus2"].iloc[0]["AlleleHist"])
		self.assertEqual("3x:1", df[df["LocusId"] == "locus3"].iloc[0]["AlleleHist"])


		self.assertAlmostEqual(0.95742710775634, df[df["LocusId"] == "locus1"].iloc[0]["AlleleStdev"])
		self.assertAlmostEqual(0.70710678118655, df[df["LocusId"] == "locus2"].iloc[0]["AlleleStdev"])
		self.assertEqual(0, df[df["LocusId"] == "locus3"].iloc[0]["AlleleStdev"])