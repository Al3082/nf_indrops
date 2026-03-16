/*
 * FastQC + MultiQC quality control on per-library reads (post-demultiplexing).
 *
 * fastqc  : runs FastQC on all three read types (R1, R2, R4) for one library
 * multiqc : aggregates all FastQC reports into a single MultiQC report
 */

process fastqc {
    tag "FastQC on ${lib_name}"
    label "fastqc"

    publishDir "${params.output_dir}/qc/fastqc", mode: 'copy'

    input:
        tuple val(lib_name), path(r1_files), path(r2_files), path(r4_files)

    output:
        path "*.{html,zip}", emit: reports

    script:
    all_files = ([r1_files, r2_files, r4_files].flatten()).join(' ')
    """
    fastqc --threads ${task.cpus} --outdir . ${all_files}
    """
}


process multiqc {
    label "multiqc"

    publishDir "${params.output_dir}/qc", mode: 'copy'

    input:
        path reports

    output:
        path "multiqc_report.html"
        path "multiqc_data"

    script:
    """
    multiqc .
    """
}
