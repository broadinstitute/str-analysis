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
    p.add_argument("trgt_bed_path", help="path of the TRGT repeat catalog BED file")
    args = p.parse_args()

    if not args.output_file:
        args.output_file = re.sub(".bed(.gz)?$", "", args.trgt_bed_path) + ".variant_catalog.json"

    if not os.path.isfile(args.trgt_bed_path):
        p.error(f"{args.trgt_bed_path} file not found")

    process_variant_catalog(args.trgt_bed_path, args.reference_fasta, args.output_file, verbose=args.verbose)


def process_variant_catalog(trgt_bed_path_path, reference_fasta_path, output_file_path, verbose=False):

    fasta_obj = pyfaidx.Fasta(reference_fasta_path, one_based_attributes=False, as_raw=True)

    print(f"Parsing {trgt_bed_path_path}")
    output_json_records = []
    counter = collections.defaultdict(int)
    fopen = gzip.open if trgt_bed_path_path.endswith("gz") else open
    with fopen(trgt_bed_path_path, "rt") as f:
        for i, row in tqdm.tqdm(enumerate(f), unit=" records", unit_scale=True):
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
                if len(key_value.split("=")) != 2:
                    print(f"WARNING: skipping invalid key-value pair '{key_value}' in line {fields}")
                    continue
                key, value = key_value.split("=")
                info_fields_dict[key] = value
            locus_id = info_fields_dict.get("ID", f"{chrom}-{start_0based}-{end_1based}")
            locus_structure = info_fields_dict.get("STRUC", "")
            motifs_from_locus_structure = parse_motifs_from_locus_structure(locus_structure)
            motifs = info_fields_dict.get("MOTIFS")
            motifs = motifs.split(",") if motifs else motifs_from_locus_structure
            if tuple(motifs) != tuple(motifs_from_locus_structure):
                print(f"WARNING: MOTIFS {tuple(motifs)} != motifs in locus structure {tuple(motifs_from_locus_structure)} in row #{i + 1}: {fields}")
                continue

            for motif in motifs:
                if not motif or (len(set(motif) - set("ACGTN")) > 0):
                    raise ValueError(f"Invalid repeat unit in row #{i + 1}: {repeat_unit}. Line: {row}")

            # TODO how do you convert to intervals
            counter["total input loci"] += 1

            # get reference sequence
            if chrom not in fasta_obj:
                raise ValueError(f"Chromosome '{chrom}' not found in reference FASTA file which has chromosomes {list(fasta_obj.keys())}")

            unprocessed_reference_sequence = fasta_obj[chrom][start_0based:end_1based]
            unprocessed_locus_structure = locus_structure
            current_start_coord_0based = start_0based

            #print("=======")
            #print("unprocessed_locus_structure", unprocessed_locus_structure)
            #print("unprocessed_reference_sequence:", unprocessed_reference_sequence)
            #print("start_coord_0based:", current_start_coord_0based)
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

                    #print("-------")
                    #print("processing left flank:", left_flank)
                    #print(f"current_start_0based/end_coord: {chrom}:{current_start_coord_0based:,d}-{current_end_coord:,d}")
                    unprocessed_reference_sequence = unprocessed_reference_sequence[len(left_flank):]
                    unprocessed_locus_structure = unprocessed_locus_structure[left_flank_end_i:]
                    current_start_coord_0based += len(left_flank)
                    #print("incrementing current_start_coord_0based to:", current_start_coord_0based)


                motif_end_i = unprocessed_locus_structure.find(")n")
                if motif_end_i == -1:
                    break

                motif = unprocessed_locus_structure[1:motif_end_i]
                #print("-------")
                #print("processing motif", motif)

                unprocessed_locus_structure = unprocessed_locus_structure[motif_end_i + 2:]
                #print("unprocessed_locus_structure:", unprocessed_locus_structure)

                # count number of repeats
                #print("-------")
                repeat_counter = 0
                for i in range(0, len(unprocessed_reference_sequence), len(motif)):
                    if not unprocessed_reference_sequence.startswith(motif):
                        break
                    repeat_counter += 1
                    unprocessed_reference_sequence = unprocessed_reference_sequence[len(motif):]
                    current_end_coord += len(motif)
                #print(f"found {repeat_counter} repeats of motif {motif}")
                #print("unprocessed_locus_structure is now:", unprocessed_locus_structure)

                if current_start_coord_0based != current_end_coord and repeat_counter > 1:
                    #print(f"incremented current_end_coord to {chrom}:{current_start_coord_0based:,d}-{current_end_coord:,d}")

                    current_locus_id = f"{locus_id}-{motif}" if len(motifs) > 1 else locus_id

                    #print("-------")
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
                    print(f"WARNING: No output loci written for {locus_id}. This is often because the reference repeat sequence contains interruptions")

    print(f"Writing {len(output_json_records):,d} records to {output_file_path}")
    fopen = gzip.open if output_file_path.endswith("gz") else open
    with fopen(output_file_path, "wt") as f:
        json.dump(output_json_records, f, indent=4, ignore_nan=True)

    print("Stats:", pformat(dict(counter)))
    print(f"Wrote {len(output_json_records):,d} records to {output_file_path}")


if __name__ == "__main__":
    main()
