"""
This script takes an ExpansionHunter catalog and a CRAM or BAM file path - either local or in a Google Storage bucket (gs://) - and extracts the subset of reads that ExpansionHunter will need to genotype that catalog into a minicram (always written as CRAM, also from a BAM input).
A gs:// BAM is read through htslib, which needs the GCS_OAUTH_TOKEN environment variable for any bucket that isn't public (e.g. export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)); a gs:// CRAM is read with the google-cloud-storage client instead.
It is designed so that ExpansionHunter's outputs when it runs on the extracted minicram will be identical to its outputs when it runs on the full CRAM.
One benefit of using this script is that it minimizes the number of bytes read from the input CRAM file compared to using samtools or other htslib-based tools. For CRAM files stored in Nearline storage 
(as is currently the case with AllOfUs genomes) this substantially reduces Nearline access costs. 

Example command-lines:

python3 -u -m str_analysis.make_minicram_for_expansion_hunter \
  -R /path/to/hg38.fa \
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
  --reference /path/to/hg38.fa \
  --reads sample1.minicram.cram \
  --sex male \
  --variant-catalog /path/to/variant_catalog.json \
  --output-prefix sample1


Thanks to Ronen Mukamel for the original idea of how to determine and retrieve the subset of CRAM containers (ie. byte ranges) that correspond to specific genomic intervals. 
"""

import argparse
import collections
import intervaltree
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
    get_local_file_stat, local_file_was_replaced, download_local_copy
from str_analysis.utils.misc_utils import parse_interval

pysam.set_verbosity(0)


def merge_nearby_regions(regions, max_distance):
    """Returns the regions with those on the same chromosome at most max_distance apart merged into one.

    extract_region merges nearby mate regions only within one locus, but mates from many loci can map to the same
    place, e.g. the hg38 pileup at chr2:32,916,000-32,917,000 (24 separate mate regions within 400bp). Each separate
    fetch from such a region reads all of it again, since the BAM linear index only points to the start of a 16kb bin.

    Args:
        regions: list of (chrom, start_0based, end) tuples.
        max_distance: largest gap in bp between two regions that are merged.

    Returns:
        list of merged (chrom, start_0based, end) tuples, sorted by chromosome and start.
    """
    merged_regions = []
    for chrom, start, end in sorted(regions, key=lambda region: (normalize_chromosome_name(region[0]), region[1])):
        if merged_regions and normalize_chromosome_name(merged_regions[-1][0]) == normalize_chromosome_name(chrom) \
                and start - merged_regions[-1][2] <= max_distance:
            previous_chrom, previous_start, previous_end = merged_regions[-1]
            merged_regions[-1] = (previous_chrom, previous_start, max(previous_end, end))
        else:
            merged_regions.append((chrom, start, end))
    return merged_regions


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
    parser.add_argument("-u", "--gcloud-project", help="Google Cloud project name to use when reading the input "
                        "file (billed for requester-pays buckets).")
    parser.add_argument("-R", "--reference-fasta", required=True, help="Reference genome FASTA file used for reading the CRAM file")
    parser.add_argument("-o", "--output-cram", help="Output CRAM file path (must end in .cram)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-L", "--region", action="append", help="Region(s) for which to extract reads (chr:start-end). "
                       "For example, for the HTT repeat locus on hg38, specify chr4:3074877-3074933")
    group.add_argument("-c", "--variant-catalog", help="Variant catalog JSON path")

    # --crai-index-path is the option's original name, kept so existing command lines still work
    parser.add_argument("-i", "--bai-or-crai-index-path", "--crai-index-path", dest="bai_or_crai_index_path",
                        help="Optional path of the input BAM or CRAM index file (.bai or .crai). This can be a local or "
                        "a gs:// path")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--enable-reference-compression-in-output-cram", action="store_true", help="Enable "
                        "reference-based output CRAM compression. This produces a smaller output CRAM but can fail "
                        "with reference MD5 errors when the supplied FASTA differs from the input CRAM header. By "
                        "default the output CRAM is written self-contained (htslib no_ref=1) to avoid those failures.")
    parser.add_argument("--output-data-transfer-stats", action="store_true", help="Write out a TSV file with stats "
                        "about the total number of bytes and containers downloaded from the CRAM")
    parser.add_argument("input_bam_or_cram", help="Input BAM or CRAM file path. This can a local or a gs:// path")

    args = parser.parse_args()

    if not args.input_bam_or_cram.endswith((".cram", ".bam")):
        parser.error(f"Input file must have a .cram or .bam extension: {args.input_bam_or_cram}")
    input_is_bam = args.input_bam_or_cram.endswith(".bam")
    if input_is_bam and args.input_bam_or_cram.startswith("gs://"):
        # htslib reads a gs:// BAM, and takes its credentials and requester-pays project only from these environment
        # variables (the google-cloud-storage client used for CRAMs gets the project from -u instead).
        if not os.environ.get("GCS_OAUTH_TOKEN"):
            print("WARNING: GCS_OAUTH_TOKEN is not set, so only a BAM in a public bucket can be read. Set it with: "
                  "export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)")
        if args.gcloud_project:
            os.environ["GCS_REQUESTER_PAYS_PROJECT"] = args.gcloud_project

    if not args.output_cram:
        args.output_cram = re.sub(r"\.(cram|bam)$", "", os.path.basename(args.input_bam_or_cram))
        args.output_cram += ".subset.cram"
    elif not args.output_cram.endswith(".cram"):
        parser.error(f"Output CRAM file must have a .cram extension: {args.output_cram}")

    # Refuse to write the output on top of the input. save_to_file replaces the output path with the extracted
    # subset and rewrites its index next to it -- so "-o in.cram in.cram" silently replaced a full CRAM with the
    # tiny minicram subset and exited 0, destroying every read outside the requested loci.
    # samefile() catches symlinks and hard links; the abspath comparison covers a local output that does not
    # exist yet, since samefile needs both paths to exist. Skipped entirely for a gs:// input, which cannot
    # name the same file as the always-local output.
    if not args.input_bam_or_cram.startswith("gs://"):
        if os.path.abspath(args.input_bam_or_cram) == os.path.abspath(args.output_cram) or (
                os.path.isfile(args.output_cram) and os.path.isfile(args.input_bam_or_cram)
                and os.path.samefile(args.input_bam_or_cram, args.output_cram)):
            parser.error(f"Output CRAM path {args.output_cram} is the same file as the input CRAM "
                         f"{args.input_bam_or_cram}. Choose a different -o path.")

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

    input_bai_or_crai_path = args.bai_or_crai_index_path if args.bai_or_crai_index_path else \
        f"{args.input_bam_or_cram}.bai" if input_is_bam else f"{args.input_bam_or_cram}.crai"
    if not file_exists(input_bai_or_crai_path):
        parser.error(f"Index path {input_bai_or_crai_path} not found")
    if input_is_bam:
        # htslib saves a remote index into the working directory and later reuses any file there with the same
        # basename, even one that belongs to a different BAM, so use a local copy cached under a hash of the full URI
        input_bai_or_crai_path = download_local_copy(input_bai_or_crai_path)

    if not file_exists(args.reference_fasta):
        parser.error(f"Reference file not found {args.reference_fasta}")

    # create a CramIntervalRreader and use it to generate a temp CRAM file containing the CRAM header and any reads
    # overlapping the user-specified region interval(s)
    # "regions", not "loci": a catalog record contributes one entry per ReferenceRegion AND one per
    # OfftargetRegions, so this count is much larger than the number of loci in the catalog
    print(f"Processing {len(args.region):,d} regions")
    print(f"Retrieving reads within {args.window_size:,d}bp of each region")
    bam_or_cram_reader = IntervalReader(args.input_bam_or_cram, input_bai_or_crai_path, verbose=args.verbose, debug=args.debug,
                                 reference_fasta_path=args.reference_fasta, cache_byte_ranges=True)

    for chrom, start, end in intervals:
        bam_or_cram_reader.add_interval(chrom, start, end)

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
        # so exit with a clear message instead of crashing on the empty file when pysam opens it below.
        # write_cram=True so that a BAM input also produces a CRAM (with a .crai, which the cleanup below expects).
        if bam_or_cram_reader.save_to_file(temporary_cram_file.name, write_cram=True) == 0:
            print(f"ERROR: No reads were found within {args.window_size:,d}bp of any of the requested "
                  f"regions in {args.input_bam_or_cram}")
            sys.exit(1)
        print(f"Completed {phase_name} in {time.time() - phase_start_time:0.2f} seconds")
        temporary_cram_file.seek(0)

        # parse the temp CRAM file and get byte ranges for mates
        phase_name = "mate discovery"
        phase_start_time = time.time()
        input_bam_file = pysam.AlignmentFile(
            temporary_cram_file.name, reference_filename=args.reference_fasta)
        # measured here for the read filter below, before any fetch() moves the file position: fetch(until_eof=True)
        # reads from the current position rather than from the first read
        longest_read_length = max((read.infer_read_length() or read.query_length
                                   for read in input_bam_file.fetch(until_eof=True)), default=0)

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
        mate_regions = []
        for chrom, start, end in deduplicated_intervals:
            if args.verbose and len(deduplicated_intervals) > 1:
                print("-"*100)

            fetch_chrom = normalized_to_reference_name.get(normalize_chromosome_name(chrom))
            if fetch_chrom is None:
                # the contig matches no reference in the input CRAM's header at all -- a wrong naming convention
                # or a different genome build. A contig that is IN the header but has no CRAI entries does not
                # reach this branch: _load_cram_containers skips only its containers and still writes the full
                # original header (see the self._cram_header_bytes write), so it stays resolvable here.
                print(f"WARNING: {chrom} is not present in the header of {args.input_bam_or_cram}; skipping mate "
                      f"discovery for {chrom}:{start}-{end}")
                continue

            genomic_regions = extract_region(
                fetch_chrom, start, end,
                input_bam=input_bam_file,
                bamlet=None,
                merge_regions_distance=args.merge_regions_distance,
                verbose=args.verbose)

            # the first region is the window itself, which the reader already has
            mate_regions.extend(genomic_regions[1:])

        for mate_region in merge_nearby_regions(mate_regions, args.merge_regions_distance):
            bam_or_cram_reader.add_interval(*mate_region)
        print(f"Completed {phase_name} in {time.time() - phase_start_time:0.2f} seconds")

        # A mate region is fetched whole, so where the read depth is extreme (e.g. the hg38 pileup at
        # chr2:32,916,000-32,917,000 that mates from many loci map to) it brings in tens of thousands of reads that
        # ExpansionHunter never uses, up to 80% of a minicram. Keep only the reads in the windows around the regions
        # and the mates of those reads.
        window_widening_for_mates = longest_read_length * 3 // 2
        print(f"Keeping the mates of reads that start up to {window_widening_for_mates:,d}bp before a window (1.5 "
              f"times the longest read, {longest_read_length:,d}bp)")
        windows_by_chrom = collections.defaultdict(intervaltree.IntervalTree)
        windows_widened_on_the_left_by_chrom = collections.defaultdict(intervaltree.IntervalTree)
        for chrom, window_start, window_end in intervals:
            windows_by_chrom[normalize_chromosome_name(chrom)].addi(window_start, window_end)
            windows_widened_on_the_left_by_chrom[normalize_chromosome_name(chrom)].addi(
                max(0, window_start - window_widening_for_mates), window_end)

        # overlap() rather than overlaps() below, since overlaps() scans every window on the chromosome when the
        # position isn't in one, as for most of the reads this filter leaves out
        def is_read_in_a_window_or_a_mate_of_one(read):
            # an unmapped read placed at its mate's position has no reference_end and covers just its start
            read_end = read.reference_end if read.reference_end is not None else read.reference_start + 1
            if windows_by_chrom[normalize_chromosome_name(read.reference_name)].overlap(read.reference_start, read_end):
                return True
            return read.is_paired and read.next_reference_id >= 0 and bool(
                windows_widened_on_the_left_by_chrom[normalize_chromosome_name(read.next_reference_name)].overlap(
                    read.next_reference_start, read.next_reference_start + 1))

        phase_name = "final primary-and-mate export"
        phase_start_time = time.time()
        print(f"Exporting data for {len(intervals)} intervals to {args.output_cram}")
        bam_or_cram_reader.save_to_file(
            args.output_cram,
            disable_reference_compression=not args.enable_reference_compression_in_output_cram,
            write_cram=True,
            read_filter=is_read_in_a_window_or_a_mate_of_one)
        print(f"Completed {phase_name} in {time.time() - phase_start_time:0.2f} seconds")
    except Exception:
        # the byte counters are only kept for the CRAM byte-range path, so leave them out for a BAM input
        print(f"ERROR: Failed during {phase_name} after {time.time() - phase_start_time:0.2f} seconds "
              f"({time.time() - start_time:0.2f} seconds total)" + ("." if input_is_bam else
              f". Loaded {bam_or_cram_reader.get_total_byte_ranges_loaded_from_cram():,d} CRAM byte ranges totaling "
              f"{bam_or_cram_reader.get_total_bytes_loaded_from_cram()/10**6:0,.1f}Mb"))
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
        bam_or_cram_reader.close()

    if not local_file_was_replaced(args.output_cram, output_cram_stat_before):
        print(f"ERROR: No output CRAM was written to {args.output_cram} because none of the requested "
              f"regions had overlapping CRAM containers in {args.input_bam_or_cram}"
              + (". The file already at that path is left over from a previous run"
                 if os.path.isfile(args.output_cram) else ""))
        sys.exit(1)

    total_duration_seconds = time.time() - start_time
    if input_is_bam:
        # the byte counters below are only kept for the CRAM byte-range path, so they would all read 0 here
        print(f"Read the BAM through htslib in {round(total_duration_seconds, 2)} seconds")
    else:
        total_bytes = bam_or_cram_reader.get_total_bytes_loaded_from_cram()
        total_containers = bam_or_cram_reader.get_total_containers_loaded_from_cram()
        total_requests = bam_or_cram_reader.get_total_byte_ranges_loaded_from_cram()
        print(f"Downloaded {total_containers:,d} containers ({total_requests:,d} read requests), "
              f"{total_bytes/10**6:0,.1f}Mb in {round(total_duration_seconds, 2)} seconds")
    if args.output_data_transfer_stats:
        bam_or_cram_reader.save_data_transfer_stats()

if __name__ == "__main__":
    main()
