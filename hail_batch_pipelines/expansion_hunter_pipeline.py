import hailtop.fs as hfs
import logging
import os
import re

from str_analysis.utils.file_utils import file_exists
from step_pipeline import pipeline, Backend, Localize, Delocalize

# Built from bw2/ExpansionHunter fea497a, which embeds the 2026-10-08 genotype-quality model
DOCKER_IMAGE = "weisburd/str-analysis-with-expansion-hunter@sha256:10b3d26b5ae8b886688fb82bc65cf6e1bcaa3452b58c92ba3f188a0e1b7cdf0e"
#DOCKER_IMAGE = "gcr.io/bw2-rare-disease/ehd-bw2:latest"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

#READ_DATA_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.downsampled_to_10x.bam"
#READ_INDEX_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.downsampled_to_10x.bam.bai"
READ_DATA_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.downsampled_to_20x.bam"
READ_INDEX_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.downsampled_to_20x.bam.bai"
READ_DATA_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram"
READ_INDEX_PATH = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram.crai"

CATALOG = "https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json.gz"

OUTPUT_BASE_DIR = "gs://bw2-delete-after-10-days/ExpansionHunter"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--catalog", default=CATALOG, help="Path of variant catalog json file")
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--analysis-mode", default="streaming", choices=("seeking", "streaming"),
                        help="ExpansionHunter analysis mode")
    parser.add_argument("--cpu", type=int, help="Number of CPUs to use")
    parser.add_argument("--input-read-data", default=READ_DATA_PATH)
    parser.add_argument("--input-read-index", default=READ_INDEX_PATH)
    parser.add_argument("--sample-sex", default="male", choices=("male", "female"))
    parser.add_argument("--min-locus-coverage", type=int, help="ExpansionHunter --min-locus-coverage threshold")
    parser.add_argument("-r", "--run-reviewer-on-locus", action="append",
                        help="Generate REViewer read visualizations for a specific locus")
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    args = bp.parse_known_args()

    if not args.force:
        existing_json_paths = bp.precache_file_paths(os.path.join(args.output_dir, f"**/*.json"))
        logging.info(f"Precached {len(existing_json_paths)} json files")

    for path in args.input_read_data, args.reference_fasta:
        if not file_exists(path):
            parser.error(f"File not found: {path}")

    read_file_stats = hfs.ls(args.input_read_data)
    if len(read_file_stats) != 1:
        parser.error(f"Expected exactly one file at {args.input_read_data}, but found {len(read_file_stats)} files")

    read_file_stats = read_file_stats[0]

    if args.analysis_mode == "streaming":
        cpu = args.cpu if args.cpu else 16
        memory = "highmem"
    else:
        cpu = args.cpu if args.cpu else 2
        memory = "standard"

    bp.set_name(f"EH (cpu={cpu}, mem={memory}) / {os.path.basename(args.input_read_data)} / {os.path.basename(args.catalog)}")

    s1 = bp.new_step(
        f"Run EH",
        arg_suffix=f"step1",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=cpu,
        memory=memory,
        localize_by=Localize.COPY,
        storage=f"{int(read_file_stats.size/10**9) + 30}Gi",
        output_dir=args.output_dir)

    local_fasta = s1.input(args.reference_fasta)
    if args.reference_fasta_fai:
        s1.input(args.reference_fasta_fai)

    local_read_data = s1.input(args.input_read_data)
    if args.input_read_index:
        s1.input(args.input_read_index)

    s1.command("set -ex")

    local_catalog = s1.input(args.catalog)

    s1.command(f"echo Genotyping $(zcat {local_catalog} | grep LocusId | wc -l) loci")
    output_prefix = re.sub(".json(.gz)?$", "", local_catalog.filename)

    min_locus_coverage_arg = "" if args.min_locus_coverage is None else f"--min-locus-coverage {args.min_locus_coverage}"
    if args.analysis_mode == "streaming":
        s1.command(f"""/usr/bin/time --verbose ExpansionHunter {min_locus_coverage_arg} \
            --reference {local_fasta} \
            --reads {local_read_data} \
            --sex {args.sample_sex} \
            --variant-catalog {local_catalog} \
            --analysis-mode streaming \
            --threads 16 \
            --output-prefix {output_prefix}""")
    else:
        s1.command(f"""/usr/bin/time --verbose ExpansionHunter {min_locus_coverage_arg} \        
            --cache-mates \
            --reference {local_fasta} \
            --reads {local_read_data} \
            --sex {args.sample_sex} \
            --variant-catalog {local_catalog} \
            --output-prefix {output_prefix}""")

    s1.command("ls -lhrt")

    s1.command(f"gzip {output_prefix}.json")
    s1.output(f"{output_prefix}.json.gz", output_dir=os.path.join(args.output_dir, f"json"))

    bp.run()

    return


    step1_output_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

    if args.run_reviewer_on_locus:
        reviewer_remote_output_dir = os.path.join(output_dir, f"svg")
        reviewer_output_prefix = re.sub("(.bam|.cram)$", "", local_bam.filename)
        s1.command(f"samtools sort {output_prefix}_realigned.bam -o {output_prefix}_realigned.sorted.bam")
        s1.command(f"samtools index {output_prefix}_realigned.sorted.bam")
        s1.command(f"""/usr/bin/time --verbose REViewer \
            --reference {local_fasta}  \
            --catalog {local_variant_catalog_path} \
            --reads {output_prefix}_realigned.sorted.bam \
            --vcf {output_prefix}.vcf \
            --output-prefix {reviewer_output_prefix}
        """)

        done_file = f"done_generating_reviewer_images_for_{output_prefix}"
        s1.command(f"touch {done_file}")

        s1.output("*.svg", output_dir=reviewer_remote_output_dir, delocalize_by=Delocalize.GSUTIL_COPY)
        s1.output(done_file, output_dir=reviewer_remote_output_dir)

    # step2: combine json files
    s2 = bp.new_step(name=f"Combine EHv5 outputs",
                     step_number=2,
                     arg_suffix=f"combine-expansion-hunter-step",
                     image=DOCKER_IMAGE,
                     cpu=1,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=args.output_dir)

    s2.depends_on(s1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_paths:
        local_path = s2.input(json_path)
        s2.command(f"ln -s {local_path}")

    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-expansion-hunter-fields "
               f"--output-prefix {output_prefix}")
    s2.command(f"bgzip {output_prefix}.{len(step1_output_paths)}_json_files.bed")
    s2.command(f"tabix {output_prefix}.{len(step1_output_paths)}_json_files.bed.gz")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.variants.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.alleles.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.bed.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.bed.gz.tbi")

    return s2



if __name__ == "__main__":
    main()


