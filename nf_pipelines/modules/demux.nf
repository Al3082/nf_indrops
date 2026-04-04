/*
 * Demultiplex inDrop v3 Kotov libraries using cutadapt on the R3 (library barcode) read.
 * Outputs one fastq.gz per library named after the barcode entry in the fasta file.
 */
process demux_libraries {
    tag "demux on ${run_id}"
    label "cutadapt"

    memory params.cutadapt_mem
    time params.cutadapt_time
    cpus params.cutadapt_threads

    publishDir "${params.output_dir}/cutadapt_stats", mode: 'copy', pattern: "*.cutadapt_stats.txt"

    input:
        tuple val(run_id), path(r3_fastq)
        path barcodes_fa

    output:
        tuple val(run_id), path("*.fastq.gz"), emit: demuxed
        path "${run_id}.cutadapt_stats.txt",   emit: stats

    script:
    """
    cutadapt \
        -e 1 \
        -j ${task.cpus} \
        -g ^file:${barcodes_fa} \
        -o "{name}.fastq.gz" \
        --action=none \
        --discard-untrimmed \
        ${r3_fastq} \
        > ${run_id}.cutadapt_stats.txt
    """
}


/*
 * Extract reads from R1, R2, R4 matching the read IDs found in a demuxed R3 file.
 * Produces three gzipped fastq files (R1, R2, R4) for one library from one Kotov run.
 */
process extract_reads {
    tag "extract ${lib_name} from ${run_id}"
    label "seqtk"

    memory params.seqtk_mem
    time params.seqtk_time
    cpus 1

    input:
        tuple val(lib_name), val(run_id), path(demuxed_r3), path(r1), path(r2), path(r4)

    output:
        tuple val(lib_name), val(run_id), \
              path("${lib_name}.${run_id}.R1.fastq.gz"), \
              path("${lib_name}.${run_id}.R2.fastq.gz"), \
              path("${lib_name}.${run_id}.R4.fastq.gz")

    script:
    """
    seqtk comp ${demuxed_r3} | cut -f1 | sort > ids.txt
    seqtk subseq ${r1} ids.txt | gzip -c > ${lib_name}.${run_id}.R1.fastq.gz
    seqtk subseq ${r2} ids.txt | gzip -c > ${lib_name}.${run_id}.R2.fastq.gz
    seqtk subseq ${r4} ids.txt | gzip -c > ${lib_name}.${run_id}.R4.fastq.gz
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

    memory params.seqtk_mem
    time params.seqtk_time
    cpus 1

    input:
        tuple val(sample_id), path(trimmed_r1), path(r2), path(r4)

    output:
        tuple val(sample_id), \
              path("${sample_id}.synced_R2.fastq.gz"), \
              path("${sample_id}.synced_R4.fastq.gz")

    script:
    """
    seqtk comp ${trimmed_r1} | cut -f1 | sort > ids.txt
    seqtk subseq ${r2} ids.txt | gzip -c > ${sample_id}.synced_R2.fastq.gz
    seqtk subseq ${r4} ids.txt | gzip -c > ${sample_id}.synced_R4.fastq.gz
    """
}
