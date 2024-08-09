"""This script converts a TRGT repeat catalog BED file to an ExpansionHunter variant catalog JSON file.
Currently, this script will skip any loci that contain interruptions in the reference repeat sequence.
"""

import argparse
import collections
import gzip
import simplejson as json
import os
from pprint import pformat
import pyfaidx
import re
import tqdm

from str_analysis.utils.eh_catalog_utils import parse_motifs_from_locus_structure


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-r", "--reference-fasta", required=True, help="Reference FASTA file path")
    p.add_argument("-o", "--output-file", help="ExansionHunter variant catalog JSON output file path")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--set-locus-id", action="store_true", help="Ignore the ID field in the TRGT catalog and set the locus ID to {chrom}-{start}-{end}-{motif}")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("trgt_bed_path", help="path of the TRGT repeat catalog BED file")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".bed(.gz)?$", "", args.trgt_bed_path) + ".variant_catalog.json"

    if not os.path.isfile(args.trgt_bed_path):
        p.error(f"{args.trgt_bed_path} file not found")

    process_variant_catalog(args.trgt_bed_path, args.reference_fasta, args.output_file,
                            set_locus_id=args.set_locus_id,
                            verbose=args.verbose,
                            show_progress_bar=args.show_progress_bar)


def process_variant_catalog(trgt_bed_path_path, reference_fasta_path, output_file_path, set_locus_id=False,
                            verbose=False, show_progress_bar=False):

    fasta_obj = pyfaidx.Fasta(reference_fasta_path, one_based_attributes=False, as_raw=True)

    print(f"Parsing {trgt_bed_path_path}")
    output_json_records = []
    counter = collections.defaultdict(int)
    fopen = gzip.open if trgt_bed_path_path.endswith("gz") else open
    with fopen(trgt_bed_path_path, "rt") as f:
        iterator = f
        if show_progress_bar:
            iterator = tqdm.tqdm(iterator, unit=" records", unit_scale=True)

        for i, row in enumerate(iterator):
            fields = row.strip("\n").split("\t")
            if len(fields) < 4:
                print(f"WARNING: skipping invalid row which has fewer than 4 columns: {fields}")
                continue
            chrom = fields[0]
            start_0based = int(fields[1])
            end_1based = int(fields[2])
            info_fields = fields[3]  # Example: ID=chr1_153777445_153777638;MOTIFS=GAAGAGGAAGTGAGA,GGGGACA;STRUC=GCAGGGAAGGACTGTGTCAGTCCAG(GAAGAGGAAGTGAGA)n(GGGGACA)nGCACGGCTTCCTGACGGTGTGACTC
            info_fields_dict = {}
            for key_value in info_fields.split(";"):
                key_value = key_value.split("=")
                if len(key_value) != 2:
                    print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
                    continue
                key, value = key_value
                info_fields_dict[key] = value

            locus_structure = info_fields_dict.get("STRUC", "")
            motifs_from_locus_structure = parse_motifs_from_locus_structure(locus_structure)
            motifs = info_fields_dict.get("MOTIFS")
            if not motifs:
                raise ValueError(f"MOTIFS field not found in row #{i + 1}: {fields}")
            motifs = motifs.split(",")
            if set(motifs) != set(motifs_from_locus_structure):
                print(f"WARNING: MOTIFS {tuple(motifs)} != motifs in locus structure {tuple(motifs_from_locus_structure)} in row #{i + 1}: {fields}")
                continue

            motifs = motifs_from_locus_structure
            for motif in motifs:
                if not motif or (len(set(motif) - set("ACGTN")) > 0):
                    raise ValueError(f"Invalid repeat unit in row #{i + 1}: {repeat_unit}. Line: {row}")

            counter["total input loci"] += 1

            # get reference sequence
            if chrom not in fasta_obj:
                raise ValueError(f"Chromosome '{chrom}' not found in reference FASTA file which has chromosomes {list(fasta_obj.keys())}")

            if len(motifs) == 1:
                motif = motifs[0]
                if set_locus_id:
                    locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
                else:
                    locus_id = info_fields_dict.get("ID")
                    if not locus_id:
                        raise ValueError(f"ID field not found in row #{i + 1}: {fields}. One work-around is to rerun with --set-locus-id")

                if verbose:
                    print(f"Adding locus {locus_id} with motif {motif} and reference region {chrom}:{start_0based}-{end_1based}")
                output_json_records.append({
                    "LocusId": locus_id,
                    "ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
                    "LocusStructure": f"({motif})*",
                    "VariantType": "Repeat",
                })

                continue

            unprocessed_reference_sequence = fasta_obj[chrom][start_0based:end_1based]
            unprocessed_locus_structure = locus_structure
            current_start_coord_0based = start_0based

            added_locus = False
            while "(" in unprocessed_locus_structure:  # continue while there are more motifs
                current_end_coord = current_start_coord_0based
                left_flank_end_i = unprocessed_locus_structure.find("(")
                if left_flank_end_i != -1:
                    left_flank = unprocessed_locus_structure[:left_flank_end_i]
                    if not unprocessed_reference_sequence.startswith(left_flank):
                        print(f"WARNING: left flank {left_flank} doesn't match reference sequence "
                             f"{unprocessed_reference_sequence[:left_flank_end_i]} @ locus {locus_id}. Skipping...")
                        break

                    unprocessed_reference_sequence = unprocessed_reference_sequence[len(left_flank):]
                    unprocessed_locus_structure = unprocessed_locus_structure[left_flank_end_i:]
                    current_start_coord_0based += len(left_flank)

                motif_end_i = unprocessed_locus_structure.find(")n")
                if motif_end_i == -1:
                    break

                motif = unprocessed_locus_structure[1:motif_end_i]
                unprocessed_locus_structure = unprocessed_locus_structure[motif_end_i + 2:]

                # count number of repeats
                repeat_counter = 0
                for i in range(0, len(unprocessed_reference_sequence), len(motif)):
                    if not unprocessed_reference_sequence.startswith(motif):
                        break
                    repeat_counter += 1
                    unprocessed_reference_sequence = unprocessed_reference_sequence[len(motif):]
                    current_end_coord += len(motif)

                if current_start_coord_0based != current_end_coord and repeat_counter > 1:

                    if set_locus_id:
                        current_locus_id = f"{chrom}-{current_start_coord_0based}-{current_end_coord}"
                    else:
                        current_locus_id = info_fields_dict.get("ID")
                        if not current_locus_id:
                            raise ValueError(f"ID field not found in row #{i + 1}: {fields}. One work-around is to rerun with --set-locus-id")

                    current_locus_id = f"{current_locus_id}-{motif}"

                    added_locus=True
                    if verbose:
                        print(f"Adding locus {current_locus_id} with motif {motif} and reference region {chrom}:{current_start_coord_0based}-{current_end_coord}")

                    output_json_records.append({
                        "LocusId": current_locus_id,
                        "ReferenceRegion": f"{chrom}:{current_start_coord_0based}-{current_end_coord}",
                        "LocusStructure": f"({motif})*",
                        "VariantType": "Repeat",
                    })

                    current_start_coord_0based = current_end_coord

            if not added_locus:
                counter["skipped"] += 1
                if verbose:
                    print(f"WARNING: Skipped {locus_id} due to issues with parsing the reference sequence. This typically happens when reference repeat sequence contains interruptions. {motifs}, {fields}")

    print(f"Parsed {counter['total input loci']:,d} loci from {trgt_bed_path_path}")
    print(f"Writing {len(output_json_records):,d} records to {output_file_path}")
    fopen = gzip.open if output_file_path.endswith("gz") else open
    with fopen(output_file_path, "wt") as f:
        json.dump(output_json_records, f, indent=4, ignore_nan=True)

    print("Stats:")
    for key, value in sorted(counter.items(), key=lambda x: x[1], reverse=True):
        print(f"{value:10,d} loci {key}")

    print(f"Wrote {len(output_json_records):,d} records to {output_file_path}")


if __name__ == "__main__":
    main()
