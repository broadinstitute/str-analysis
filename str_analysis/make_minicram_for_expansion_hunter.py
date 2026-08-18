"""
This script takes an ExpansionHunter catalog and a CRAM file path - either local or in a Google Storage bucket (gs://) - and extracts the subset of reads that ExpansionHunter will need to genotype that catalog. 
It is designed so that ExpansionHunter's outputs when it runs on the extracted minicram will be identical to its outputs when it runs on the full CRAM.
One benefit of using this script is that it minimizes the number of bytes read from the input CRAM file, which, for CRAM files stored in Nearline storage 
(as is currently the case with AllOfUs genomes) also minimizes the Nearline access costs. 

Example command-lines:

python3 -u -m str_analysis.make_minicram_for_expansion_hunter \
  -R gs://path/to/hg38.fa \
  -c /path/to/variant_catalog.json \
  -i gs://bucket/path/to/sample1.cram.crai \
  -o sample1.minicram.cram \
  --output-data-transfer-stats \
  gs://bucket/path/to/sample1.cram

# The ExpansionHunter command below contains 3 flags that are only available in the https://github.com/bw2/ExpansionHunter fork, but the minicram works equally well with the original https://github.com/Illumina/ExpansionHunter
ExpansionHunter \
  --dont-output-consensus-sequences \
  --compress-output-files \
  --analysis-mode optimized-streaming \
  --reference gs://path/to/hg38.fa \
  --reads sample1.minicram.cram \
  --sex male \
  --variant-catalog /path/to/variant_catalog.json \
  --output-prefix sample1
"""

import argparse
import json
import os
import re
import sys
import pysam
import tempfile
import time

from str_analysis.make_bamlet import extract_region
from str_analysis.utils.cram_bam_utils import IntervalReader, normalize_chromosome_name
from str_analysis.utils.file_utils import set_requester_pays_project, file_exists, open_file, \
    get_local_file_stat, local_file_was_replaced
from str_analysis.utils.misc_utils import parse_interval

pysam.set_verbosity(0)


def main():
    parser = argparse.ArgumentParser(description="A script to generate BAMlets")
    # NOTE: make_bamlet.py defines this same arg with different help text, and that text is correct there:
    # make_bamlet passes a real bamlet, so it runs jump_for_mates and each merged region really does become a
    # single disk read. This tool passes bamlet=None, so that retrieval never happens here and only the merged
    # intervals themselves matter -- hence the different wording below.
    parser.add_argument("-d", "--merge-regions-distance", type=int, default=1000, help="Mate merge distance. "
        "Mates found within this distance of each other are merged into a single interval, and those intervals "
        "determine which CRAM containers get downloaded. Increasing this widens the merged intervals, which "
        "tends to increase the total number of bytes fetched; decreasing it keeps the intervals tighter.")
    parser.add_argument("-w", "--window-size", type=int, default=1000, help="Window size in bp to include around the "
                        "user-specified region(s). This is useful for including read pairs that may overlap the region(s)")
    parser.add_argument("-u", "--gcloud-project", help="Google Cloud project name to use when reading the input cram.")
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file used for reading the CRAM file")
    parser.add_argument("-o", "--output-cram", help="Output CRAM file path (must end in .cram)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-L", "--region", action="append", help="Region(s) for which to extract reads (chr:start-end). "
                       "For example, for the HTT repeat locus on hg38, specify chr4:3074877-3074933")
    group.add_argument("-c", "--variant-catalog", help="Variant catalog JSON path")

    parser.add_argument("-i", "--crai-index-path", help="Optional path of the input CRAM index file. This can be a "
                        "local or a gs:// path")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--enable-reference-compression-in-output-cram", action="store_true", help="Enable "
                        "reference-based output CRAM compression. This produces a smaller output CRAM but can fail "
                        "with reference MD5 errors when the supplied FASTA differs from the input CRAM header. By "
                        "default the output CRAM is written self-contained (htslib no_ref=1) to avoid those failures.")
    parser.add_argument("--output-data-transfer-stats", action="store_true", help="Write out a TSV file with stats "
                        "about the total number of bytes and containers downloaded from the CRAM")
    parser.add_argument("input_cram", help="Input CRAM file path. This can a local or a gs:// path")

    args = parser.parse_args()

    if not args.input_cram.endswith(".cram"):
        parser.error(f"Input CRAM file must have a .cram extension: {args.input_cram}")

    if not args.output_cram:
        args.output_cram = re.sub(".cram$", "", os.path.basename(args.input_cram))
        args.output_cram += ".subset.cram"
    elif not args.output_cram.endswith(".cram"):
        parser.error(f"Output CRAM file must have a .cram extension: {args.output_cram}")

    # Refuse to write the output on top of the input. save_to_file replaces the output path with the extracted
    # subset and rewrites its index next to it -- so "-o in.cram in.cram" silently replaced a full CRAM with the
    # tiny minicram subset and exited 0, destroying every read outside the requested loci.
    # samefile() catches symlinks and hard links; the abspath comparison covers a local output that does not
    # exist yet, since samefile needs both paths to exist. Skipped entirely for a gs:// input, which cannot
    # name the same file as the always-local output.
    if not args.input_cram.startswith("gs://"):
        if os.path.abspath(args.input_cram) == os.path.abspath(args.output_cram) or (
                os.path.isfile(args.output_cram) and os.path.isfile(args.input_cram)
                and os.path.samefile(args.input_cram, args.output_cram)):
            parser.error(f"Output CRAM path {args.output_cram} is the same file as the input CRAM "
                         f"{args.input_cram}. Choose a different -o path.")

    # set before the variant-catalog download below, not after it: download_local_copy reads this module global at
    # call time to build its "gsutil -u <project>" flag, so setting it later left the catalog fetch -- the first
    # gs:// read this tool does -- without the user's -u/--gcloud-project.
    set_requester_pays_project(args.gcloud_project)

    if args.variant_catalog:
        args.region = []
        with open_file(args.variant_catalog, download_local_copy_before_opening=True) as f:
            variant_catalog_records = json.load(f)
            for i, record in enumerate(variant_catalog_records):
                if "ReferenceRegion" not in record:
                    parser.error(f"Record #{i+1} does not have a ReferenceRegion field: {record}")

                reference_regions = record["ReferenceRegion"]
                if not isinstance(reference_regions, list):
                    reference_regions = [reference_regions]

                # OfftargetRegions are extracted too. ExpansionHunter recruits reads from a locus's off-target
                # regions when genotyping large expansions, so a minicram built from a catalog that defines them
                # but containing only the ReferenceRegion reads can yield different genotypes than running the
                # same catalog against the full CRAM -- which defeats the purpose of this tool. To keep the
                # minicram smaller, pass a catalog without them (the repo ships variant_catalog_without_offtargets).
                offtarget_regions = record.get("OfftargetRegions") or []
                if not isinstance(offtarget_regions, list):
                    offtarget_regions = [offtarget_regions]

                for region in reference_regions + offtarget_regions:
                    if args.debug:
                        print(f"DEBUG: Adding", record["LocusId"], "region", region)
                    args.region.append(region)

    if args.debug:
        args.verbose=True

    start_time = time.time()

    # validate args
    intervals = []
    for region in args.region:
        try:
            chrom, start, end = parse_interval(region)
        except ValueError as e:
            parser.error(f"Unable to parse region {region}: {e}")

        # parse_interval builds the chromosome from everything before the last colon, so a region written as
        # ":100-200" yields an empty name without raising. add_interval rejects it, but only from deep inside the
        # reader with a bare ValueError; report it here as an ordinary CLI error naming the region.
        if not chrom:
            parser.error(f"Missing chromosome name in region {region}")

        window_start = max(0, start - args.window_size)
        window_end = end + args.window_size
        # parse_interval does not check that start < end, so a zero-width or reversed region combined with a
        # small --window-size yields an empty window. add_interval would reject it, but only from deep inside
        # the reader with a bare ValueError; report it here as an ordinary CLI error naming the region.
        if window_end <= window_start:
            parser.error(f"Region {region} is empty after applying --window-size {args.window_size:,d} "
                         f"({chrom}:{window_start}-{window_end})")
        intervals.append((chrom, window_start, window_end))

    input_crai_path = args.crai_index_path if args.crai_index_path else f"{args.input_cram}.crai"
    if not file_exists(input_crai_path):
        parser.error(f"CRAM index path {input_crai_path} not found")

    if not file_exists(args.reference_fasta):
        parser.error(f"Reference file not found {args.reference_fasta}")

    # create a CramIntervalRreader and use it to generate a temp CRAM file containing the CRAM header and any reads
    # overlapping the user-specified region interval(s)
    # "regions", not "loci": a catalog record contributes one entry per ReferenceRegion AND one per
    # OfftargetRegions, so this count is much larger than the number of loci in the catalog
    print(f"Processing {len(args.region):,d} regions")
    print(f"Retrieving reads within {args.window_size:,d}bp of each region")
    cram_reader = IntervalReader(args.input_cram, input_crai_path, verbose=args.verbose, debug=args.debug,
                                 reference_fasta_path=args.reference_fasta, cache_byte_ranges=True)

    for chrom, start, end in intervals:
        cram_reader.add_interval(chrom, start, end)

    # stat'd before the run so the "was an output written" check below can tell a file this run wrote from one an
    # earlier run left at the same path -- os.path.isfile alone cannot, and would report success for a run that
    # wrote nothing. print_reads.py uses the same pair of calls around its own save_to_file.
    output_cram_stat_before = get_local_file_stat(args.output_cram)
    temporary_cram_file = tempfile.NamedTemporaryFile(suffix=".cram", delete=False)
    input_bam_file = None
    phase_name = "initial interval export"
    phase_start_time = time.time()
    try:
        # if no reads overlap the requested region(s), save_to_file returns 0 without writing the temp file,
        # so exit with a clear message instead of crashing on the empty file when pysam opens it below
        if cram_reader.save_to_file(temporary_cram_file.name) == 0:
            print(f"ERROR: No reads were found within {args.window_size:,d}bp of any of the requested "
                  f"regions in {args.input_cram}")
            sys.exit(1)
        print(f"Completed {phase_name} in {time.time() - phase_start_time:0.2f} seconds")
        temporary_cram_file.seek(0)

        # parse the temp CRAM file and get byte ranges for mates
        phase_name = "mate discovery"
        phase_start_time = time.time()
        input_bam_file = pysam.AlignmentFile(
            temporary_cram_file.name, reference_filename=args.reference_fasta)

        # map normalized chrom name -> the exact reference name in this file's header, so fetch() below always uses
        # a name valid for the header even when the user's -L/catalog used a different naming convention ("9" vs
        # "chr9"). IntervalReader.save_to_file does the same for its own fetch calls; without it, extract_region
        # passes the raw user-supplied name straight to pysam and raises "invalid contig".
        normalized_to_reference_name = {
            normalize_chromosome_name(name): name for name in input_bam_file.references
        }

        # deduplicated because mate discovery re-scans the reads in every entry it is given, and the raw list can
        # repeat the same coordinates -- a catalog record whose ReferenceRegion and OfftargetRegions overlap, or
        # two records sharing a region, made extract_region fetch the identical window more than once for no
        # additional mates. sorted() only to keep the traversal (and its verbose output) deterministic.
        deduplicated_intervals = sorted(set(intervals))
        for chrom, start, end in deduplicated_intervals:
            if args.verbose and len(deduplicated_intervals) > 1:
                print("-"*100)

            fetch_chrom = normalized_to_reference_name.get(normalize_chromosome_name(chrom))
            if fetch_chrom is None:
                # the contig matches no reference in the input CRAM's header at all -- a wrong naming convention
                # or a different genome build. A contig that is IN the header but has no CRAI entries does not
                # reach this branch: _load_cram_containers skips only its containers and still writes the full
                # original header (see the self._cram_header_bytes write), so it stays resolvable here.
                print(f"WARNING: {chrom} is not present in the header of {args.input_cram}; skipping mate "
                      f"discovery for {chrom}:{start}-{end}")
                continue

            genomic_regions = extract_region(
                fetch_chrom, start, end,
                input_bam=input_bam_file,
                bamlet=None,
                merge_regions_distance=args.merge_regions_distance,
                verbose=args.verbose)

            for genomic_region in genomic_regions:
                cram_reader.add_interval(*genomic_region)
        print(f"Completed {phase_name} in {time.time() - phase_start_time:0.2f} seconds")

        phase_name = "final primary-and-mate export"
        phase_start_time = time.time()
        print(f"Exporting data for {len(intervals)} intervals to {args.output_cram}")
        cram_reader.save_to_file(
            args.output_cram,
            disable_reference_compression=not args.enable_reference_compression_in_output_cram)
        print(f"Completed {phase_name} in {time.time() - phase_start_time:0.2f} seconds")
    except Exception:
        print(f"ERROR: Failed during {phase_name} after {time.time() - phase_start_time:0.2f} seconds "
              f"({time.time() - start_time:0.2f} seconds total). Loaded "
              f"{cram_reader.get_total_byte_ranges_loaded_from_cram():,d} CRAM byte ranges totaling "
              f"{cram_reader.get_total_bytes_loaded_from_cram()/10**6:0,.1f}Mb")
        raise
    finally:
        if input_bam_file is not None:
            input_bam_file.close()
        # temporary_cram_file was created with delete=False so it persists after its handle is closed;
        # remove it and the .crai index that save_to_file generated next to it
        temporary_cram_file.close()
        for temp_path in (temporary_cram_file.name, f"{temporary_cram_file.name}.crai"):
            if os.path.isfile(temp_path):
                os.remove(temp_path)
        cram_reader.close()

    if not local_file_was_replaced(args.output_cram, output_cram_stat_before):
        print(f"ERROR: No output CRAM was written to {args.output_cram} because none of the requested "
              f"regions had overlapping CRAM containers in {args.input_cram}"
              + (". The file already at that path is left over from a previous run"
                 if os.path.isfile(args.output_cram) else ""))
        sys.exit(1)

    total_bytes = cram_reader.get_total_bytes_loaded_from_cram()
    total_containers = cram_reader.get_total_containers_loaded_from_cram()
    total_requests = cram_reader.get_total_byte_ranges_loaded_from_cram()
    total_duration_seconds = time.time() - start_time
    print(f"Downloaded {total_containers:,d} containers ({total_requests:,d} read requests), "
          f"{total_bytes/10**6:0,.1f}Mb in {round(total_duration_seconds, 2)} seconds")
    if args.output_data_transfer_stats:
        cram_reader.save_data_transfer_stats()

if __name__ == "__main__":
    main()
