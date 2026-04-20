/*
 * Demultiplex inDrop v3 Kotov libraries using cutadapt on the R3 (library barcode) read.
 * Outputs one fastq.gz per library named after the barcode entry in the fasta file.
 */
process demux_libraries {
    tag "demux on ${run_id}"
    label "cutadapt"

    publishDir "${params.output_dir}/cutadapt_stats", mode: 'copy', pattern: "*.cutadapt_stats.txt"

    input:
        tuple val(run_id), path(r3_fastq)
        path barcodes_fa

    output:
        tuple val(run_id), path("*.fastq.gz"),  emit: demuxed
        path "${run_id}.cutadapt_stats.txt",    emit: stats
        path "${run_id}.demux_counts.tsv",      emit: counts

    script:
    """
    set -euo pipefail

    cutadapt \
        -e 1 \
        -j ${task.cpus} \
        -g ^file:${barcodes_fa} \
        -o "{name}.fastq.gz" \
        --action=none \
        --discard-untrimmed \
        ${r3_fastq} \
        > ${run_id}.cutadapt_stats.txt

    # Count raw input and per-library demuxed reads
    raw_n=\$(zcat ${r3_fastq} | awk 'END{print NR/4}')
    echo -e "run\\tlibrary\\tstep\\treads" > ${run_id}.demux_counts.tsv
    echo -e "${run_id}\\tall\\traw_r3\\t\${raw_n}" >> ${run_id}.demux_counts.tsv
    input_base=\$(basename "${r3_fastq}" .fastq.gz)
    for f in *.fastq.gz; do
        lib=\$(basename "\$f" .fastq.gz)
        [[ "\$lib" == "\$input_base" ]] && continue
        n=\$(zcat "\$f" | awk 'END{print NR/4}')
        echo -e "${run_id}\\t\${lib}\\tdemuxed\\t\${n}" >> ${run_id}.demux_counts.tsv
    done
    """
}


/*
 * Extract reads from R1, R2, R4 matching the read IDs found in a demuxed R3 file.
 * Produces three gzipped fastq files (R1, R2, R4) for one library from one Kotov run.
 */
process extract_reads {
    tag "extract ${lib_name} from ${run_id}"
    label "seqtk"

    input:
        tuple val(lib_name), val(run_id), path(demuxed_r3), path(r1), path(r2), path(r4)

    output:
        tuple val(lib_name), val(run_id), \
              path("${lib_name}.${run_id}.R1.fastq.gz"), \
              path("${lib_name}.${run_id}.R2.fastq.gz"), \
              path("${lib_name}.${run_id}.R4.fastq.gz"),  emit: reads
        path "${lib_name}.${run_id}.extract_stats.tsv",   emit: stats

    script:
    """
    set -euo pipefail

    # Full read names from demuxed R3 — format-agnostic (works for SRA and raw-submitter headers)
    seqtk comp ${demuxed_r3} | cut -f1 > ids.txt

    seqtk subseq ${r1} ids.txt | gzip -c > ${lib_name}.${run_id}.R1.fastq.gz
    seqtk subseq ${r2} ids.txt | gzip -c > ${lib_name}.${run_id}.R2.fastq.gz
    seqtk subseq ${r4} ids.txt | gzip -c > ${lib_name}.${run_id}.R4.fastq.gz

    # Stats: count demuxed and extracted reads
    demuxed_n=\$(wc -l < ids.txt)
    extracted_n=\$(seqtk comp ${lib_name}.${run_id}.R1.fastq.gz | wc -l)
    echo -e "library\\trun\\tstep\\treads" > ${lib_name}.${run_id}.extract_stats.tsv
    echo -e "${lib_name}\\t${run_id}\\tdemux_ids\\t\${demuxed_n}" >> ${lib_name}.${run_id}.extract_stats.tsv
    echo -e "${lib_name}\\t${run_id}\\textracted\\t\${extracted_n}" >> ${lib_name}.${run_id}.extract_stats.tsv
    """
}


/*
 * Re-sync barcode reads (R2, R4) with a trimmed R1 that may have had reads
 * discarded by cutadapt -m. Extracts read IDs from trimmed R1 and filters
 * R2/R4 to only those IDs using seqtk subseq.
 */
process sync_reads {
    tag "sync ${sample_id}"
    label "seqtk"

    publishDir "${params.output_dir}/processed_fastqs", mode: 'copy', pattern: "*.synced_R*.fastq.gz"

    input:
        tuple val(sample_id), path(trimmed_r1), path(r2), path(r4)

    output:
        tuple val(sample_id), \
              path("${sample_id}.synced_R2.fastq.gz"), \
              path("${sample_id}.synced_R4.fastq.gz"),  emit: reads
        path "${sample_id}.sync_stats.tsv",             emit: stats

    script:
    """
    set -euo pipefail

    # Full read names from trimmed R1 — format-agnostic (works for SRA and raw-submitter headers)
    seqtk comp ${trimmed_r1} | cut -f1 > ids.txt

    seqtk subseq ${r2} ids.txt | gzip -c > ${sample_id}.synced_R2.fastq.gz
    seqtk subseq ${r4} ids.txt | gzip -c > ${sample_id}.synced_R4.fastq.gz

    # Stats: reads before (R2 input) and after sync (actual synced output)
    before_n=\$(seqtk comp ${r2} | wc -l)
    after_n=\$(seqtk comp ${sample_id}.synced_R2.fastq.gz | wc -l)
    echo -e "sample\\tstep\\treads" > ${sample_id}.sync_stats.tsv
    echo -e "${sample_id}\\tbefore_trim_sync\\t\${before_n}" >> ${sample_id}.sync_stats.tsv
    echo -e "${sample_id}\\tafter_trim_sync\\t\${after_n}" >> ${sample_id}.sync_stats.tsv
    """
}


/*
 * Re-sync a single target read file with a reference read that may have had
 * reads discarded (e.g. by cutadapt -m).  Used for v2 where only R1 (barcode)
 * needs syncing after R2 (gene read) is trimmed.
 */
process sync_single_read {
    tag "sync ${sample_id}"
    label "seqtk"

    publishDir "${params.output_dir}/processed_fastqs", mode: 'copy', pattern: "*.synced.fastq.gz"

    input:
        tuple val(sample_id), path(reference_read), path(target_read)

    output:
        tuple val(sample_id), path("${sample_id}.synced.fastq.gz"), emit: reads
        path "${sample_id}.sync_stats.tsv",                         emit: stats

    script:
    """
    set -euo pipefail

    # Full read names from the surviving (reference) read — format-agnostic (works for SRA and raw-submitter headers)
    seqtk comp ${reference_read} | cut -f1 > ids.txt

    seqtk subseq ${target_read} ids.txt | gzip -c > ${sample_id}.synced.fastq.gz

    # Stats: count actual synced output, not just input IDs
    before_n=\$(seqtk comp ${target_read} | wc -l)
    after_n=\$(seqtk comp ${sample_id}.synced.fastq.gz | wc -l)
    echo -e "sample\\tstep\\treads" > ${sample_id}.sync_stats.tsv
    echo -e "${sample_id}\\tbefore_trim_sync\\t\${before_n}" >> ${sample_id}.sync_stats.tsv
    echo -e "${sample_id}\\tafter_trim_sync\\t\${after_n}" >> ${sample_id}.sync_stats.tsv
    """
}
