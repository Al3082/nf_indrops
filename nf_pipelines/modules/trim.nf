/*
 * Trim adapter sequences, poly-A tails, and low-quality 3' bases from gene reads
 * using cutadapt, before STARsolo alignment.
 *
 * Reads shorter than 20 bp after trimming are discarded (-m 20).
 * Downstream processes must re-sync barcode reads (R2/R4) to match
 * the surviving R1 read IDs.
 */
process trim_gene_read {
    tag "trim ${sample_id}"
    label "cutadapt"

    publishDir "${params.output_dir}/cutadapt_trim_stats", mode: 'copy', pattern: "*.cutadapt_trim_stats.txt"
    publishDir "${params.output_dir}/processed_fastqs", mode: 'copy', pattern: "*.trimmed_R1.fastq.gz"

    input:
        tuple val(sample_id), path(r1)
        path adapter_fasta

    output:
        tuple val(sample_id), path("${sample_id}.trimmed_R1.fastq.gz"), emit: trimmed
        path "${sample_id}.cutadapt_trim_stats.txt",                    emit: stats
        path "${sample_id}.trim_counts.tsv",                            emit: counts

    script:
    """
    input_n=\$(zcat ${r1} | awk 'END{print NR/4}')

    cutadapt \
        -a file:${adapter_fasta} \
        --nextseq-trim=20 \
        --poly-a \
        -m 20 \
        -j ${task.cpus} \
        -o ${sample_id}.trimmed_R1.fastq.gz \
        ${r1} \
        > ${sample_id}.cutadapt_trim_stats.txt

    output_n=\$(zcat ${sample_id}.trimmed_R1.fastq.gz | awk 'END{print NR/4}')
    echo -e "sample\\tstep\\treads" > ${sample_id}.trim_counts.tsv
    echo -e "${sample_id}\\tbefore_trim\\t\${input_n}" >> ${sample_id}.trim_counts.tsv
    echo -e "${sample_id}\\tafter_trim\\t\${output_n}" >> ${sample_id}.trim_counts.tsv
    """
}
