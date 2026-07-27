#!/bin/bash

# pipefail so a cellranger failure is not masked by the exit code of `tee`
set -eo pipefail

# Parse the sample name from the CSV
CSV="demux.config.csv"

echo "Resolving relative links in \$CSV"
resolve_links.py "\$CSV"
cat "\$CSV"

echo
echo "Contents of the hashtags CSV:"
cat hashtags.csv
echo

echo "Running cellranger multi" | tee "demultiplexed_samples.log.txt"
cellranger --version 2>&1 | tee -a "demultiplexed_samples.log.txt"
cellranger multi \
            --id="demultiplexed_samples" \
            --csv="\${CSV}.resolved.csv" \
            --localcores=${task.cpus} \
            --localmem=${task.memory.toGiga()} \
    2>&1 | tee -a "demultiplexed_samples.log.txt"

echo "Finished running cellranger multi - " | tee -a "demultiplexed_samples.log.txt"

if [ -d "demultiplexed_samples" ]; then
    if [ -d "demultiplexed_samples/SC_MULTI_CS" ]; then
        echo "Cleaning up demultiplexed_samples/SC_MULTI_CS" | tee -a "demultiplexed_samples.log.txt"
        rm -r "demultiplexed_samples/SC_MULTI_CS"
    fi

    if [ -d "demultiplexed_samples/outs" ]; then
        echo "Cleaning up demultiplexed_samples/outs" | tee -a "demultiplexed_samples.log.txt"
        mv "demultiplexed_samples/outs/"* "demultiplexed_samples/"
        rmdir "demultiplexed_samples/outs"
    fi

    # Move the per-sample BAM files to the top level directory
    for sample in demultiplexed_samples/per_sample_outs/*; do
        BAM=\$sample/count/sample_alignments.bam
        if [ -s "\$BAM" ]; then

            # Samples with no cells called yield a header-only BAM, which
            # bamtofastq produces no output for. Leave those BAMs in
            # per_sample_outs/ so they are not picked up downstream.
            RC=0
            bam_has_alignments.py "\$BAM" || RC=\$?
            if [ "\$RC" -eq 1 ]; then
                echo Skipping \${sample##*/} - no aligned reads in \$BAM | tee -a "demultiplexed_samples.log.txt"
                continue
            elif [ "\$RC" -ne 0 ]; then
                exit "\$RC"
            fi

            DEST=demultiplexed_samples/\${sample##*/}.bam
            echo Moving BAM file from \$sample/count/sample_alignments.bam to \$DEST | tee -a "demultiplexed_samples.log.txt"
            mv "\$BAM" "\$DEST"
            mv "\${BAM}.bai" "\${DEST}.bai"
        fi
    done

fi
echo "Completed - demultiplexing samples" | tee -a "demultiplexed_samples.log.txt"
