import argparse
import gzip
import os
import re


def main():
    p = argparse.ArgumentParser(description="""This script takes a .bed file and outputs a new .bed file where all loci
    have been trimmed so that their size is a multiple of their motif size. For example, if one of the input file rows 
    has coordinates chr1:12345-12358 and motif "CAG" (specified in the 4th column) then the output .bed file would 
    contain the same row but with the end coordinate 12358 changed to 12357 since this would make the interval 
    size = 12bp which is an exact multiple of the 3bp motif. The script assumes the .bed file uses 0-based coordinates 
    (as described in the .bed format spec), but also has a --one-based-coords option for non-standard .bed files 
    (such as the GangSTR repeat spec file). Also, it assumes the motif is specified in the 4th column by default, 
    but has a --motif-column arg which can be used to change this column #.
    """)
    p.add_argument("-1", "--one-based-coords", action="store_true", help="If specified, the .bed file will be parsed"
        "as having 1-based start coordinates (which is non-standard, but used in some cases such as the GangSTR input "
        "spec files)")
    p.add_argument("-o", "--output-path", help="Optional output file path")
    p.add_argument("--motif-column", help="Index of the motif column (counting from 1)", default=4, type=int)
    p.add_argument("input_bed", help=".bed or .bed.gz input file to trim")
    args = p.parse_args()

    if not os.path.isfile(args.input_bed):
        p.error(f"Input file not found: {args.input_bed}")

    if not args.output_path:
        args.output_path = re.sub(".bed(.gz)?$", "", args.input_bed) + ".trimmed.bed"

    print(f"Writing to {args.output_path}")
    open_func = gzip.open if args.input_bed.endswith(".gz") else open
    trim_counter = 0
    with open_func(args.input_bed, "rt") as f, open(args.output_path, "wt") as f2:
        for i, line in enumerate(f, start=1):
            fields = line.strip("\n").split("\t")
            if args.motif_column > len(fields):
                p.error(f"Unable to parse line #{i}. The motif column is column #{args.motif_column} but the "
                        f"line contains only {len(fields)} columns: {fields}")
            motif = fields[args.motif_column - 1]
            if not motif:
                p.error(f"Unable to parse line #{i}. The motif is an empty string:: {fields}")
            if not re.match("^[a-zA-Z]+$", motif):
                p.error(f"Unable to parse line #{i}. The motif contains characters other than letters: {motif}")

            try:
                start_0based = int(fields[1]) if not args.one_based_coords else (int(fields[1]) - 1)
                end_1based = int(fields[2])
            except ValueError as e:
                p.error(f"Unable to parse start or end coordinates in line #{i}: {fields}")

            trim_amount = (end_1based - start_0based) % len(motif)
            if trim_amount > 0:
                trim_counter += 1
                fields[2] = str(end_1based - trim_amount)
            f2.write("\t".join(fields) + "\n")

    print(f"Wrote {i} lines to {args.output_path}. Trimmed {trim_counter} lines ({100*trim_counter/i:0.1f}%)")

    if args.input_bed.endswith(".gz") and os.path.isfile(f"{args.input_bed}.tbi"):
        print(f"Running bgzip on {args.output_path}")
        os.system(f"bgzip -f {args.output_path}")
        print(f"Running tabix on {args.output_path}.gz")
        os.system(f"tabix {args.output_path}.gz")
    print("Done")


if __name__ == "__main__":
    main()
