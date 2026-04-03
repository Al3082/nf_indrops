/*
 * Trim adapter sequences, poly-A tails, and low-quality 3' bases from gene reads
 * using cutadapt, before STARsolo alignment.
 */
process trim_gene_read {
    tag "trim ${sample_id}"
    label "cutadapt"

    publishDir "${params.output_dir}/cutadapt_trim_stats", mode: 'copy', pattern: "*.cutadapt_trim_stats.txt"

    input:
        tuple val(sample_id), path(r1)
        path adapter_fasta

    output:
        tuple val(sample_id), path("${sample_id}.trimmed_R1.fastq.gz"), emit: trimmed
        path "${sample_id}.cutadapt_trim_stats.txt",                    emit: stats

    script:
    """
    cutadapt \
        -a file:${adapter_fasta} \
        --nextseq-trim=20 \
        --poly-a \
        -j ${task.cpus} \
        -o ${sample_id}.trimmed_R1.fastq.gz \
        ${r1} \
        > ${sample_id}.cutadapt_trim_stats.txt
    """
}
