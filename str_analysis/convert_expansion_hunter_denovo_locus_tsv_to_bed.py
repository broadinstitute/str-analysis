import argparse
import os
import pandas as pd

"""
$1              contig : chr6
$2               start : 38957549
$3                 end : 38957550
$4               motif : AAAAAAAAAAAAAAAAAAC
$5        num_anc_irrs : 1
$6   norm_num_anc_irrs : 0.59
$7        het_str_size : 8
"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_path")
    p.add_argument("output_path")
    args = p.parse_args()

    print(f"Reading from {args.input_path}")
    f = open(args.output_path, "wt")
    df = pd.read_table(args.input_path)
    df = df.sort_values(["contig", "start", "end"])
    for _, row in df.iterrows():
        fields = [
            row.contig,
            row.start,
            row.end,
            f"{len(row.motif)}bp:{row.motif}:{row.end-row.start}bp-ref-locus-size:{row.num_anc_irrs}-irrs:{row.het_str_size}x-het-str-size",
            ".",
        ]
        f.write("\t".join(map(str, fields)) + "\n")

    f.close()

    print(f"Wrote {len(df)} rows to {args.output_path}")


if __name__ == "__main__":
    main()