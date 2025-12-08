"""Convert TRFinder output .dat format to .bed

Example:

	trf409.macosx chr22.fa 2 10000000 10000000  80 10 2 2000 -h -l 6 -ngs > chr22.dat

	python3 -u -m str_analysis.convert_dat_to_bed chr22.dat
"""
import argparse
import os
import re

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.dat_utils import parse_dat_file


def main():
	p = argparse.ArgumentParser(description="Convert Tandem Repeat Finder output .dat format to .bed")
	p.add_argument("--min-motif-size", type=int, default=1, help="Minimum motif size")
	p.add_argument("--max-motif-size", type=int, default=1_000_000, help="Maximum motif size")
	p.add_argument("--min-repeat-count", type=int, default=3, help="Minimum number of repeats")
	p.add_argument("--min-base-pairs", type=int, default=9, help="Minimum base pairs spanned by repeat sequence")
	p.add_argument("-L", "--interval", action="append", help="Interval filter (eg. 'chr1:12345-54321'). "
															 "Can be specified multiple times.")
	p.add_argument("-o", "--output-path", help="Output BED file path")
	p.add_argument("-v", "--verbose", action="store_true")
	p.add_argument("dat_file_path", help="Input DAT file path")
	args = p.parse_args()

	if not os.path.isfile(args.dat_file_path):
		p.error(f"{args.dat_file_path} file not found")

	if not args.output_path:
		args.output_path = re.sub(".dat(.gz)?$", "", args.dat_file_path) + ".bed"

	with open(args.output_path, "wt") as f:
		for record in sorted(
			parse_dat_file(args.dat_file_path),
			key=lambda d: (d['chrom'], d['start_0based'], d['end_1based'])
		):
			if len(record['repeat_unit']) < args.min_motif_size or len(record['repeat_unit']) > args.max_motif_size:
				continue
			if record['repeat_count'] < args.min_repeat_count:
				continue
			if record['end_1based'] - record['start_0based'] < args.min_base_pairs:
				continue

			if args.interval:
				# Example: chr1:12345-54321
				for i in args.interval:
					chrom, start, end = parse_interval(i)
					if record['chrom'].replace("chr", "").upper() != chrom.replace("chr", "").upper():
						continue

					# check if current interval contains this record
					if record['start_0based'] >= start and record['end_1based'] <= end:
						break
				else:
					# none of the intervals contain this recrod
					continue

			f.write("\t".join(map(str, [
				record['chrom'],
				record['start_0based'],
				record['end_1based'],
				record['repeat_unit'],
				record['percent_matches'],
				".",
			])) + "\n")


if __name__ == "__main__":
	main()
