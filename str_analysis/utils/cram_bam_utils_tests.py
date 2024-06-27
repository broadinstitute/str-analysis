import os
import pkgutil
import tempfile
import unittest
import subprocess

from str_analysis.utils.cram_bam_utils import IntervalReader
from str_analysis.utils.file_utils import set_requester_pays_project


class TestCramBamUtils(unittest.TestCase):

	def setUp(self):
		self._FXN_intervals = [
			("chr9", 69037287, 69037304),
		]
		set_requester_pays_project("cmg-analysis")

	def test_interval_reader_non_overlapping_intervals(self):
		reader = IntervalReader("gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram")

		reader.add_interval("chr1", 1, 5)
		reader.add_interval("chr1", 6, 15)
		reader.add_interval("chr2", 16, 30)
		reader.add_interval("chr2", 1, 15)

		self.assertEqual(reader._get_merged_intervals(), [
			("chr1", 1, 5),
			("chr1", 6, 15),
			("chr2", 1, 15),
			("chr2", 16, 30),
		])

	def test_interval_reader_overlapping_intervals(self):
		reader = IntervalReader("gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram")

		reader.add_interval("chr1", 1, 10)
		reader.add_interval("chr1", 5, 15)
		reader.add_interval("chr1", 20, 30)
		reader.add_interval("chr2", 1, 10)
		reader.add_interval("chr2", 5, 15)
		reader.add_interval("chr2", 20, 30)

		self.assertEqual(reader._get_merged_intervals(), [
			("chr1", 1, 15),
			("chr1", 20, 30),
			("chr2", 1, 15),
			("chr2", 20, 30),
		])

	def test_interval_reader_adjacent_intervals(self):
		reader = IntervalReader("gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram")

		reader.add_interval("chr1", 1, 5)
		reader.add_interval("chr1", 5, 15)
		reader.add_interval("chr2", 15, 30)
		reader.add_interval("chr2", 1, 15)

		self.assertEqual(reader._get_merged_intervals(), [
			("chr1", 1, 15),
			("chr2", 1, 30),
		])


	def test_cram_reader_on_local_files(self):
		with tempfile.NamedTemporaryFile(suffix=".cram") as input_cram_file, \
			  tempfile.NamedTemporaryFile(suffix=".crai") as input_crai_file, \
			  tempfile.NamedTemporaryFile(suffix=".bam") as input_bam_file, \
			  tempfile.NamedTemporaryFile(suffix=".bai") as input_bai_file, \
			  tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file, \
			  tempfile.NamedTemporaryFile(suffix=".bam") as output_bam_file:

			input_cram_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.cram"))
			input_cram_file.flush()
			input_crai_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.cram.crai"))
			input_crai_file.flush()

			input_bam_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam"))
			input_bam_file.flush()
			input_bai_file.write(pkgutil.get_data("str_analysis", "data/tests/FXN.wgsim_HET_250xGAA.bam.bai"))
			input_bai_file.flush()

			cram_reader = IntervalReader(input_cram_file.name, input_crai_file.name)
			for interval in self._FXN_intervals:
				cram_reader.add_interval(*interval)
			cram_reads_counter = cram_reader.save_to_file(output_cram_file.name)

			bam_reader = IntervalReader(input_bam_file.name, input_bai_file.name)
			for interval in self._FXN_intervals:
				bam_reader.add_interval(*interval)
			bam_reads_counter = bam_reader.save_to_file(output_bam_file.name)

			#print(f"Retrieved {cram_reads_counter} reads from CRAM and {bam_reads_counter} reads from BAM")
			self.assertEqual(cram_reads_counter, bam_reads_counter)

	def test_cram_reader_on_google_storage_files(self):
		try:
			os.environ["GCS_OAUTH_TOKEN"] = subprocess.check_output("gcloud auth application-default print-access-token", shell=True).decode("utf-8").strip()
		except Exception as e:
			print(f"WARNING: Unable to set GCS_OAUTH_TOKEN: {e}. Skipping gs:// tests..")
			return

		with tempfile.NamedTemporaryFile(suffix=".cram") as output_cram_file, \
			tempfile.NamedTemporaryFile(suffix=".bam") as output_bam_file:
		
			# check cloud paths
			cram_reader = IntervalReader(
				"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram",
				"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram.crai",
				verbose=True)
			for interval in self._FXN_intervals:
				cram_reader.add_interval(*interval)

			cram_reads_counter = cram_reader.save_to_file(output_cram_file.name)
			#print(f"Retrieved {cram_reads_counter} reads from CRAM using default mode")

			cram_reader2 = IntervalReader(
				"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram",
				"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.cram.crai")
			for interval in self._FXN_intervals:
				cram_reader2.add_interval(*interval)
			cram_reads_counter2 = cram_reader2.save_to_file(output_bam_file.name)
			#print(f"Retrieved {cram_reads_counter2} reads from CRAM")
			self.assertEqual(cram_reads_counter, cram_reads_counter2)

			bam_reader = IntervalReader(
				"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.bam",
				"gs://str-analysis/tests/FXN.wgsim_HET_250xGAA.bam.bai")
			for interval in self._FXN_intervals:
				bam_reader.add_interval(*interval)
			bam_reads_counter = bam_reader.save_to_file(output_bam_file.name)

			#print(f"Retrieved {bam_reads_counter} reads from BAM")
			self.assertEqual(cram_reads_counter, bam_reads_counter)
