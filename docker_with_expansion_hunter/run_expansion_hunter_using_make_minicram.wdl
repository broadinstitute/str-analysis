version 1.0

# Runs make_minicram_for_expansion_hunter to localize the reads ExpansionHunter needs for genotyping 
# a specific TR catalog and input CRAM, and then runs ExpansionHunter (bw2 fork, optimized-streaming mode)
# on the minicram.
#
# While ExpansionHunter runs, a background helper process copies its incomplete output files to
# <output_bucket>/checkpoints/<cram name without .cram>.<hash of the inputs>/ every checkpoint_interval_seconds, so
# that, if the VM is preempted, the subsequent rerun resumes from the last checkpoint instead of starting over. The
# checkpoint folder is deleted once the json.gz is in output_bucket.
#
# The following output files are written to the output_bucket:
#   Step 1: <output_bucket>/<cram name without .cram>.minicram.cram (and .crai)
#   Step 2: <output_bucket>/<cram name without .cram>.json.gz
# 


workflow RunExpansionHunterUsingMakeMinicram {
    input {
        File cram
        String? crai_path
        File reference_fasta
        File reference_fasta_fai
        File tandem_repeat_catalog
        String sex
        String output_bucket
        String? sample_id
        String? requester_pays_project

        # The latest image built from docker_with_expansion_hunter/ (see sha256_dockerhub.txt)
        String docker = "weisburd/str-analysis-with-expansion-hunter@sha256:f6375d0e7ff3c35bbc0c1888c845f3c128703be23a6cb7ffef5d5f8b6f9dba78"
        Int checkpoint_interval_seconds = 600
        Int? make_minicram_disk_size_gb
        Int? expansion_hunter_disk_size_gb
        Int preemptible_tries = 3
    }

    parameter_meta {
        cram: "Input CRAM. Cromwell doesn't copy it to the VM: make_minicram_for_expansion_hunter reads only the byte ranges it needs straight from the bucket"
        crai_path: "gs:// path of the CRAM index. Defaults to <cram>.crai"
        reference_fasta: "Reference genome FASTA that the CRAM was aligned to"
        reference_fasta_fai: "Index of reference_fasta. It must be in the same bucket folder as reference_fasta"
        tandem_repeat_catalog: "ExpansionHunter variant catalog (.json or .json.gz)"
        sex: "'male' or 'female' (in any capitalization), passed to ExpansionHunter --sex"
        output_bucket: "gs:// folder where the minicram, the ExpansionHunter json.gz and the checkpoints are written"
        sample_id: "SampleId written to the ExpansionHunter output. Defaults to the CRAM file name without .cram"
        requester_pays_project: "Google Cloud project to bill when reading the CRAM from a requester-pays bucket"
        make_minicram_disk_size_gb: "Overrides the minicram step's disk size, which defaults to 3x the CRAM size plus the reference and catalog. That's enough for a catalog with loci spread across the genome, which touches most of the CRAM's containers"
        expansion_hunter_disk_size_gb: "Overrides the ExpansionHunter step's disk size, which defaults to the size of the inputs plus 10 GB. That's enough for catalogs of up to ~1M loci. The whole 5.6M-locus TRExplorer catalog needs ~25 GB more than the inputs"
    }

    String output_dir = sub(output_bucket, "/+$", "")
    String output_prefix = basename(cram, ".cram")

    call MakeMinicram {
        input:
            cram = cram,
            crai_path = crai_path,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            tandem_repeat_catalog = tandem_repeat_catalog,
            output_dir = output_dir,
            output_prefix = output_prefix,
            requester_pays_project = requester_pays_project,
            docker = docker,
            disk_size_gb_override = make_minicram_disk_size_gb,
            preemptible_tries = preemptible_tries,
    }

    call RunExpansionHunter {
        input:
            cram_path = "~{cram}",
            minicram = MakeMinicram.minicram,
            minicram_crai = MakeMinicram.minicram_crai,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            tandem_repeat_catalog = tandem_repeat_catalog,
            sex = sex,
            sample_id = select_first([sample_id, output_prefix]),
            output_dir = output_dir,
            output_prefix = output_prefix,
            docker = docker,
            checkpoint_interval_seconds = checkpoint_interval_seconds,
            disk_size_gb_override = expansion_hunter_disk_size_gb,
            preemptible_tries = preemptible_tries,
    }

    output {
        String expansion_hunter_json_path = RunExpansionHunter.json_path
    }
}

task MakeMinicram {
    input {
        File cram
        String? crai_path
        File reference_fasta
        File reference_fasta_fai
        File tandem_repeat_catalog
        String output_dir
        String output_prefix
        String? requester_pays_project
        String docker
        Int? disk_size_gb_override
        Int preemptible_tries
    }

    parameter_meta {
        # Cromwell passes the gs:// path instead of copying the CRAM to the VM
        cram: {
            localization_optional: true
        }
    }

    # make_minicram_for_expansion_hunter keeps every CRAM byte range it downloads in a temp file, writes them to a
    # second temp file, and keeps its first-pass CRAM until the minicram is written. With a catalog that touches most
    # of the CRAM's containers, each of these is about as large as the CRAM.
    Int disk_size_gb = select_first([disk_size_gb_override,
        ceil(3 * size(cram, "GB") + size(reference_fasta, "GB") + size(reference_fasta_fai, "GB") + size(tandem_repeat_catalog, "GB")) + 10])

    String minicram_filename = "~{output_prefix}.minicram.cram"

    command <<<
        set -euxo pipefail

        # make_minicram_for_expansion_hunter keeps the CRAM byte ranges it downloads in temp files. 
        # Put them on the task's disk rather than on the container's small boot disk.
        export TMPDIR="$(pwd)/tmp"
        mkdir -p "$TMPDIR"

        python3 -u -m str_analysis.make_minicram_for_expansion_hunter \
            -R ~{reference_fasta} \
            -c ~{tandem_repeat_catalog} \
            ~{"-i " + crai_path} \
            ~{"-u " + requester_pays_project} \
            -o ~{minicram_filename} \
            ~{cram}

        ls -l ~{minicram_filename} ~{minicram_filename}.crai

        # save the minicram for future reuse such as generating read visualizations for specific loci
        gcloud storage cp ~{minicram_filename} ~{minicram_filename}.crai ~{output_dir}/
    >>>

    output {
        File minicram = minicram_filename
        File minicram_crai = "~{minicram_filename}.crai"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "6 GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        preemptible: preemptible_tries
    }
}

task RunExpansionHunter {
    input {
        String cram_path
        File minicram
        File minicram_crai
        File reference_fasta
        File reference_fasta_fai
        File tandem_repeat_catalog
        String sex
        String sample_id
        String output_dir
        String output_prefix
        String docker
        Int checkpoint_interval_seconds
        Int? disk_size_gb_override
        Int preemptible_tries
    }

    # In addition to the CRAM, reference and catalog, the disk must store ExpansionHunter's per-contig temp output files,
    # the snapshot of them that the checkpoint uploader copies to the output bucket to enable --resume, and the
    # final output files. In previous tests, these took at most 5.4GB with a ~225k-locus catalog. They grow with the
    # number of loci, so much larger catalogs need disk_size_gb_override.
    Int disk_size_gb = select_first([disk_size_gb_override,
        ceil(size(minicram, "GB") + size(reference_fasta, "GB") + size(reference_fasta_fai, "GB") + size(tandem_repeat_catalog, "GB")) + 10])

    String minicram_filename = "~{output_prefix}.minicram.cram"
    String processed_loci_filename = "~{output_prefix}.processed_loci.txt"

    command <<<
        set -euxo pipefail

        # Cromwell localizes the minicram to a path that includes the workflow id, and ExpansionHunter's --resume run
        # signature includes the --reads path. Linking it to a fixed name keeps that path the same when the workflow is
        # resubmitted, so the new run can resume from the old run's checkpoint. It also puts the .crai where
        # ExpansionHunter looks for it.
        ln -s ~{minicram} ~{minicram_filename}
        ln -s ~{minicram_crai} ~{minicram_filename}.crai

        # - --output-motif-composition all-loci raises ExpansionHunter's default --max-depth from 150 to 500, so
        #   --max-depth 150 keeps the calls the same as a run without motif composition, which is what
        #   ExpansionHunter's genotype-quality model was trained on.
        # - Without --sample-id, ExpansionHunter would name the sample after the --reads file, which here is
        #   <cram name without .cram>.minicram
        # - --sex is lower-cased because ExpansionHunter accepts only 'male' or 'female', while sample tables often have
        #   'Male' or 'Female'.
        EXPANSION_HUNTER_ARGS=(
            --analysis-mode optimized-streaming
            --sample-id ~{sample_id}
            --sex "$(echo '~{sex}' | tr '[:upper:]' '[:lower:]')"
            --dont-output-consensus-sequences
            --output-motif-composition all-loci
            --max-depth 150
            --compress-output-files
        )

        # ExpansionHunter resumes from a checkpoint only if its run signature (reads, reference,
        # catalog path and size, sex and genotyping options) matches that of the previous run. Naming the checkpoint
        # folder after a hash of these inputs ensures that this condition is met. cram_path is in the hash so that two
        # CRAMs with the same file name in different bucket folders don't share a checkpoint, since ExpansionHunter's
        # run signature only sees the minicram's fixed file name.
        CHECKPOINT_KEY=$(printf '%s\n' "~{docker}" "~{cram_path}" "~{reference_fasta}" "~{tandem_repeat_catalog}" \
            "$(md5sum < ~{tandem_repeat_catalog})" "${EXPANSION_HUNTER_ARGS[*]}" | sha256sum | cut -c 1-12)
        CHECKPOINT_URL="~{output_dir}/checkpoints/~{output_prefix}.${CHECKPOINT_KEY}/"

        # Check whether an earlier attempt left a checkpoint, distinguishing "none exists" from "error while checking"
        CHECKPOINT_LS=""
        CHECKPOINT_FOUND=""
        CHECKPOINT_WAIT=10
        for try in 1 2 3 4 5; do
            if CHECKPOINT_LS=$(gcloud storage ls "${CHECKPOINT_URL}~{processed_loci_filename}" 2>&1); then
                CHECKPOINT_FOUND=yes
                break
            fi
            case "$CHECKPOINT_LS" in
                *"matched no objects"*)
                    CHECKPOINT_FOUND=no
                    break ;;
            esac
            echo "Could not check for a checkpoint (try $try of 5): $CHECKPOINT_LS"
            if [ "$try" -lt 5 ]; then
                sleep $CHECKPOINT_WAIT
                CHECKPOINT_WAIT=$(( CHECKPOINT_WAIT * 2 ))
            fi
        done
        if [ "$CHECKPOINT_FOUND" = "yes" ]; then
            echo "Found a checkpoint from an earlier attempt, resuming: $CHECKPOINT_URL"
            gcloud storage cp "${CHECKPOINT_URL}*" .
        elif [ "$CHECKPOINT_FOUND" = "no" ]; then
            echo "No checkpoint for this sample, ExpansionHunter will genotype the whole catalog"
        else
            echo "Gave up checking whether a checkpoint exists: $CHECKPOINT_LS"
            echo "Failing this attempt rather than starting over, which would overwrite it"
            exit 1
        fi

        # Keep copying ExpansionHunter's unfinished output to the checkpoint folder while it runs. Each round first copies
        # the processed-loci list into a local snapshot folder and then the per-contig temp files. cp -pu re-copies only the
        # files that changed since the last round, and keeps their modification times, so gcloud storage rsync skips
        # the contigs that are already finished. A round whose copies fail (e.g. ExpansionHunter deleting its temp
        # files at the end) uploads nothing.
        #
        # The upload keeps the same order: the temp files first, and the list only once all of them have uploaded, so
        # that the list in the bucket only ever describes bytes that every temp file in the bucket holds.
        #
        # The loop sleeps in 1-second steps so that it notices the stop file written once ExpansionHunter finishes
        # right away, while letting an upload that's already under way finish before the checkpoint folder is deleted.
        mkdir -p checkpoint_snapshot/contig_files
        (
            set +ex
            while true; do
                for second in $(seq ~{checkpoint_interval_seconds}); do
                    [ -e stop_checkpoint_uploader ] && exit 0
                    sleep 1
                done
                if [ -s ~{processed_loci_filename} ] \
                        && cp -p ~{processed_loci_filename} checkpoint_snapshot/ \
                        && cp -pu ~{output_prefix}.contig*.json ~{output_prefix}.contig*.vcf checkpoint_snapshot/contig_files/; then
                    gcloud storage rsync checkpoint_snapshot/contig_files "$CHECKPOINT_URL" \
                        && gcloud storage cp checkpoint_snapshot/~{processed_loci_filename} "${CHECKPOINT_URL}~{processed_loci_filename}" \
                        || echo "Checkpoint upload failed, will try again next round"
                fi
                echo "Checkpoint round disk usage: $(df -h . | tail -1)"
            done
        ) &
        CHECKPOINT_UPLOADER_PID=$!
        echo "Checkpoint uploader running as $CHECKPOINT_UPLOADER_PID, every ~{checkpoint_interval_seconds}s, to $CHECKPOINT_URL"

        # --resume is passed whether or not a checkpoint was restored: with nothing to resume from, ExpansionHunter
        # starts from the beginning and keeps its --resume files up to date as it goes.
        /usr/bin/time --verbose ExpansionHunter "${EXPANSION_HUNTER_ARGS[@]}" \
            --threads 4 \
            --resume \
            --reference ~{reference_fasta} \
            --reads ~{minicram_filename} \
            --variant-catalog ~{tandem_repeat_catalog} \
            --output-prefix ~{output_prefix}

        touch stop_checkpoint_uploader
        wait $CHECKPOINT_UPLOADER_PID

        ls -lhrt

        gcloud storage cp ~{output_prefix}.json.gz ~{output_dir}/

        # The json.gz is in the bucket, so the checkpoint is no longer needed
        gcloud storage rm "${CHECKPOINT_URL}**" || echo "Could not delete the checkpoint folder $CHECKPOINT_URL"
    >>>

    output {
        String json_path = "~{output_dir}/~{output_prefix}.json.gz"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "13 GB"
        disks: "local-disk ~{disk_size_gb} HDD"
        preemptible: preemptible_tries
    }
}
