"""This script is a lighter-weight alternative to GATK PrintReads
(https://gatk.broadinstitute.org/hc/en-us/articles/360036883571-PrintReads).
It exracts data for genomic regions of a CRAM or BAM file.

Specifically for CRAM files stored in Google Cloud buckets, this script minimizes the number of bytes read from the input CRAM file compared to using 'samtools view' or other htslib-based approaches. 
For CRAM files stored in Nearline storage (as is currently the case with AllOfUs genomes) this substantially reduces Nearline access costs. 

Thanks to Ronen Mukamel for the original idea of how to determine and retrieve the subset of CRAM containers (ie. byte ranges) that corespond to specific genomic intervals. 

Example command-line:

python3 -u -m str_analysis.print_reads \
  -R gs://path/to/hg38.fa \
  -L chr4:12345-54321  -L chrM \
  -o sample1.minicram.cram \
  --output-data-transfer-stats \
  gs://bucket/path/to/sample1.cram

"""

import argparse
import os
import sys
import pysam

from str_analysis.utils.file_utils import set_requester_pays_project, file_exists, open_file, \
	get_local_file_stat, local_file_was_replaced
from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.cram_bam_utils import IntervalReader
pysam.set_verbosity(0)


def main():
	parser = argparse.ArgumentParser(description="Retrieve reads from a CRAM or BAM file that overlap one or more "
												 "genomic intervals, or are unmapped.")
	parser.add_argument("-u", "--gcloud-project", help="Google Cloud project to use for GCS requester pays buckets.")
	parser.add_argument("-o", "--output", help="Output file path. If not specified, it will be based on the input filename.")
	parser.add_argument("--read-index", help="Optional path of the input BAM or CRAM index file. This can be a local or a gs:// path")
	parser.add_argument("-R", "--reference-fasta", help="Optional reference genome FASTA file used when reading the CRAM file")
	parser.add_argument("-L", "--interval", action="append", default=[], help="The script will output aligned reads that "
						"overalp the given genomic interval(s). The value can be the path of a .bed file, .bed.gz file, "
						".interval_list file, an interval in the format \"chr:start-end\" with a 0-based start "
						"coordinate, or a chromosome name (ie. \"chrM\"). If a read within a given interval has an "
						"unmapped mate, the mate will be included in the output. If a mate is mapped somewhere beyond "
						"the specified interval(s), the mate will not be included. Read pairs where both mates are "
						"unmapped will not be included in the output unless --include-unmapped-read-pairs is specified.")
	parser.add_argument("--padding", type=int, default=0, help="Number of bases with which to pad each interval")
	parser.add_argument("--include-unmapped-read-pairs", action="store_true",
						help="Output read pairs where both mates are unmapped. This can be specified in addition to or "
							 "instead of -L intervals.")
	parser.add_argument("--enable-reference-compression-in-output-cram", action="store_true", help="Enable "
						"reference-based output CRAM compression. This produces a smaller output CRAM but can fail "
						"with reference MD5 errors when the supplied FASTA differs from the input CRAM header. By "
						"default the output CRAM is written self-contained (htslib no_ref=1) to avoid those failures.")
	parser.add_argument("--output-data-transfer-stats", action="store_true", help="Write out a TSV file with stats "
						"about the total number of bytes and containers downloaded from the CRAM")
	parser.add_argument("--verbose", action="store_true")
	parser.add_argument("--debug", action="store_true")
	parser.add_argument("input_bam_or_cram", help="Input BAM or CRAM file. This can be a local or a gs:// path")
	args = parser.parse_args()

	# validate args
	if args.debug:
		args.verbose = True

	if not args.interval and not args.include_unmapped_read_pairs:
		parser.error("At least one --interval or --include-unmapped-read-pairs arg must be specified")

	# checked here rather than left to IntervalReader, which rejects any other extension with a bare ValueError
	# traceback from its constructor instead of the usual argparse usage message
	if not args.input_bam_or_cram.endswith((".cram", ".bam")):
		parser.error(f"Input file must have a .cram or .bam extension: {args.input_bam_or_cram}")
	input_is_cram = args.input_bam_or_cram.endswith(".cram")

	set_requester_pays_project(args.gcloud_project)
	for path in args.input_bam_or_cram, args.read_index, args.reference_fasta:
		if path and not file_exists(path):
			parser.error(f"{path} not found")

	if not args.read_index:
		# IntervalReader falls back to "<input>.crai" / "<input>.bai" when --read-index is not specified, so check
		# that derived path too. Otherwise a missing index surfaced as a pysam or parse_crai_index error from deep
		# inside the reader rather than as a CLI error naming the file.
		read_index_path = f"{args.input_bam_or_cram}.crai" if input_is_cram else f"{args.input_bam_or_cram}.bai"
		if not file_exists(read_index_path):
			parser.error(f"Index file {read_index_path} not found")

	output_extension = ".cram" if input_is_cram else ".bam"
	if args.output is None:
		# sliced rather than str.removesuffix() so this keeps working on the Python 3.7 and 3.8 interpreters that
		# setup.py still declares support for. The extension is guaranteed present by the check above.
		args.output = os.path.basename(args.input_bam_or_cram)[:-len(output_extension)]
		args.output += f".print_reads{output_extension}"
	elif not args.output.endswith(output_extension):
		# save_to_file writes the format of the INPUT file regardless of this path's extension, so a mismatched
		# extension produces, for example, a CRAM named ".bam" that downstream tools then fail to open
		parser.error(f"Output file must have a {output_extension} extension because the input is a "
					 f"{output_extension[1:].upper()} file: {args.output}")

	# Refuse to write the output on top of the input. save_to_file replaces the output path with the extracted
	# subset and rewrites its index next to it, so "-o in.bam in.bam" would silently replace a full BAM with the
	# tiny subset and exit 0, destroying every read outside the requested intervals. samefile() catches symlinks
	# and hard links; the abspath comparison covers an output that does not exist yet, since samefile needs both
	# paths to exist. Skipped entirely for a gs:// input, which cannot name the same file as the always-local output.
	if not args.input_bam_or_cram.startswith("gs://"):
		if os.path.abspath(args.input_bam_or_cram) == os.path.abspath(args.output) or (
				os.path.isfile(args.output) and os.path.isfile(args.input_bam_or_cram)
				and os.path.samefile(args.input_bam_or_cram, args.output)):
			parser.error(f"Output path {args.output} is the same file as the input {args.input_bam_or_cram}. "
						 f"Choose a different -o path.")

	# the intervals are parsed and validated BEFORE the IntervalReader is constructed because, for a CRAM, the
	# constructor already does real I/O -- it downloads and parses the CRAI and then reads the CRAM header, which
	# for a gs:// input is network traffic (a real CRAI runs to megabytes). Validating first keeps a malformed -L
	# from costing that download before the CLI error is reported. make_minicram_for_expansion_hunter.py does the
	# same. The checks below cover both things IntervalReader.add_interval would reject -- an empty chromosome name
	# (reachable only from -L, since a .bed line's leading empty field is eaten by strip() and caught by the
	# column count) and an empty window -- so that either is reported here as an ordinary CLI error naming the
	# offending interval rather than as a raw ValueError traceback out of the reader.
	intervals = []
	for interval in args.interval:
		if interval.endswith(".bed") or interval.endswith(".bed.gz") or interval.endswith(".interval_list"):
			if not file_exists(interval):
				parser.error(f"{interval} file not found")

			with open_file(interval, is_text_file=True) as f:
				for line in f:
					if line.startswith("@"):
						continue
					fields = line.strip().split("\t")
					if len(fields) < 3:
						parser.error(f"Expected at least 3 columns in line {line}")
					chrom, start, end = fields[:3]
					start_offset = 1 if interval.endswith(".interval_list") else 0
					start = int(start) - start_offset  # convert to 0-based coordinates
					end = int(end)
					if start > end:
						parser.error(f"start coordinate {start} is greater than the end coordinate {end}")
					if end + args.padding <= max(0, start - args.padding):
						parser.error(f"Invalid interval in {interval}: {chrom}:{start}-{end} is empty after "
									 f"applying --padding {args.padding:,d}")
					intervals.append((chrom, start - args.padding, end + args.padding))
		else:
			try:
				if ":" in interval:
					chrom, start_0based, end = parse_interval(interval)
				else:
					chrom, start_0based, end = interval, 0, 10**9
			except ValueError as e:
				parser.error(f"Invalid interval {interval}  {e}")

			if not chrom:
				parser.error(f"Missing chromosome name in interval {interval}")
			if start_0based > end:
				parser.error(f"start coordinate {start_0based} is greater than end coordinate {end}")
			if ":" in interval:
				if end + args.padding <= max(0, start_0based - args.padding):
					parser.error(f"Invalid interval {interval}: it is empty after applying "
								 f"--padding {args.padding:,d}")
				intervals.append((chrom, start_0based - args.padding, end + args.padding))
			else:
				# a bare chromosome name is taken whole, so --padding does not apply to it
				intervals.append((chrom, start_0based, end))

	reader = IntervalReader(
		args.input_bam_or_cram,
		crai_or_bai_path=args.read_index,
		include_unmapped_read_pairs=args.include_unmapped_read_pairs,
		reference_fasta_path=args.reference_fasta,
		verbose=args.verbose,
		debug=args.debug,
	)

	for chrom, start, end in intervals:
		reader.add_interval(chrom, start, end)

	# stat'd before the run because save_to_file returns 0 both when it writes an output containing no reads and
	# when it returns without writing anything at all (which it does when no CRAM containers matched), and
	# os.path.isfile alone cannot tell those apart: a file left at this path by an earlier run passed that check,
	# so the run exited 0 while describing someone else's output as empty.
	output_stat_before = get_local_file_stat(args.output)
	try:
		read_count = reader.save_to_file(
			args.output,
			disable_reference_compression=not args.enable_reference_compression_in_output_cram)
	finally:
		# in a finally so that an exception anywhere above still releases the reader's GCS client, its open input
		# file handle, and its temp byte-range cache file
		reader.close()

	# The message says "matched the request" rather than naming intervals because --include-unmapped-read-pairs can
	# be given without any -L interval at all. The exit code differs between a CRAM and a BAM input holding the
	# same reads, and that is deliberate: for a CRAM with no matching containers nothing is written and this is a
	# failure, while the BAM path always writes a (valid, empty) output, which is reported as a warning below.
	if not local_file_was_replaced(args.output, output_stat_before):
		print(f"ERROR: Nothing was written to {args.output} because no data in {args.input_bam_or_cram} matched "
			  f"the request"
			  + (". The file already at that path is left over from a previous run"
				 if os.path.isfile(args.output) else ""))
		sys.exit(1)

	# written only once the run is known to have produced an output, matching make_minicram_for_expansion_hunter.py.
	# The counters it reports are read off the reader and stay valid after close().
	if args.output_data_transfer_stats:
		reader.save_data_transfer_stats()

	if read_count == 0:
		print(f"WARNING: No reads matched the request, so {args.output} is empty")
	else:
		print(f"Wrote {read_count:,d} reads to {args.output}")

if __name__ == "__main__":
	main()
