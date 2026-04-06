/*
 * dropEst pipeline processes for inDrop v2/v3 scRNA-seq data.
 *
 * Replaces STARsolo with: dropTag (barcode extraction) → plain STAR → dropEst (quantification).
 * See https://github.com/kharchenkolab/dropEst
 */


/*
 * dropTag for inDrop v3 libraries.
 *
 * Input reads arrive as lists (one entry per run/source: Kotov runA, runB, Briggs).
 * dropTag accepts a single file per read type, so we concatenate first.
 *
 * Read order for droptag v3:  R2 (CB1, 8 bp)  R4 (CB2+UMI, 14 bp)  R1 (gene read)
 */
process droptag_v3 {
    tag "dropTag_v3 on ${lib_name}"
    label "dropest"
    stageInMode 'copy'

    memory params.dropest_mem
    time params.dropest_time
    cpus params.dropest_threads

    input:
        tuple val(lib_name), val(r1_files), val(r2_files), val(r4_files)
        path droptag_xml

    output:
        tuple val(lib_name), path("${lib_name}.*.fastq.gz"), path("*.params.gz"), emit: tagged

    script:
    r1_list = r1_files instanceof List ? r1_files : [r1_files]
    r2_list = r2_files instanceof List ? r2_files : [r2_files]
    r4_list = r4_files instanceof List ? r4_files : [r4_files]

    cat_r1 = r1_list.size() > 1 ? "cat ${r1_list.join(' ')} > merged_R1.fastq.gz" : "ln -s ${r1_list[0]} merged_R1.fastq.gz"
    cat_r2 = r2_list.size() > 1 ? "cat ${r2_list.join(' ')} > merged_R2.fastq.gz" : "ln -s ${r2_list[0]} merged_R2.fastq.gz"
    cat_r4 = r4_list.size() > 1 ? "cat ${r4_list.join(' ')} > merged_R4.fastq.gz" : "ln -s ${r4_list[0]} merged_R4.fastq.gz"

    """
    ${cat_r1}
    ${cat_r2}
    ${cat_r4}

    droptag \
        -c ${droptag_xml} \
        -S -s \
        -p ${task.cpus} \
        -n ${lib_name} \
        merged_R2.fastq.gz \
        merged_R4.fastq.gz \
        merged_R1.fastq.gz
    """
}


/*
 * dropTag for inDrop v3 from pre-merged FASTQs.
 *
 * Skips the merge step — takes already-concatenated R1/R2/R4 files directly.
 */
process droptag_v3_premerged {
    tag "dropTag_v3 on ${lib_name}"
    label "dropest"

    memory params.dropest_mem
    time params.dropest_time
    cpus params.dropest_threads

    input:
        tuple val(lib_name), path(r1), path(r2), path(r4)
        path droptag_xml

    output:
        tuple val(lib_name), path("${lib_name}.*.fastq.gz"), path("*.params.gz"), emit: tagged

    script:
    """
    echo "=== Staged XML (${droptag_xml}) ==="
    cat ${droptag_xml}
    echo "=== End XML ==="

    droptag \
        -c ${droptag_xml} \
        -S -s \
        -p ${task.cpus} \
        -n ${lib_name} \
        ${r2} \
        ${r4} \
        ${r1}
    """
}


/*
 * dropTag for inDrop v2 libraries.
 *
 * Read order for droptag v2:  R1 (barcode: CB1|adapter|CB2|UMI)  R2 (gene read)
 */
process droptag_v2 {
    tag "dropTag_v2 on ${lib_name}"
    label "dropest"
    stageInMode 'copy'

    memory params.dropest_mem
    time params.dropest_time
    cpus params.dropest_threads

    input:
        tuple val(lib_name), path(r1), path(r2)
        path droptag_xml

    output:
        tuple val(lib_name), path("${lib_name}.*.fastq.gz"), path("*.params.gz"), emit: tagged

    script:
    """
    droptag \
        -c ${droptag_xml} \
        -S -s \
        -p ${task.cpus} \
        -n ${lib_name} \
        ${r1} ${r2}
    """
}


/*
 * Plain STAR alignment on dropTag-tagged FASTQs.
 *
 * No --soloType: barcode information is already embedded in read names by dropTag.
 * The params file is passed through to dropEst.
 */
process star_plain {
    tag "STAR on ${lib_name}"
    label "star"

    memory params.star_mem
    time params.star_time
    cpus params.star_threads

    input:
        tuple val(lib_name), path(tagged_reads), path(params_file)
        val genome_dir

    output:
        tuple val(lib_name), path("STAR/Aligned.sortedByCoord.out.bam"), path(params_file), emit: aligned

    script:
    reads_str = tagged_reads instanceof List ? tagged_reads.join(',') : tagged_reads
    """
    mkdir -p STAR

    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${reads_str} \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 10 \
        --outFilterScoreMinOverLread 0.66 \
        --outFilterMatchNminOverLread 0.66 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.05 \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 100000 \
        --outSJfilterIntronMaxVsReadN 50000 100000 200000 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix STAR/
    """
}


/*
 * dropEst: UMI deduplication and gene expression quantification.
 *
 * Produces count matrices and cell-level statistics from the aligned BAM.
 * Gene parts "eiEIBA" counts reads from exons and introns.
 */
process dropest_quant {
    tag "dropEst on ${lib_name}"
    label "dropest"

    memory params.dropest_mem
    time params.dropest_time
    cpus 1

    publishDir "${params.output_dir}/${lib_name}/dropest", mode: 'copy'

    input:
        tuple val(lib_name), path(bam), path(params_file)
        path gtf
        path droptag_xml

    output:
        path "*"

    script:
    """
    dropest \
        -r ${params_file} \
        -w -V -m \
        -L eiEIBA \
        -g ${gtf} \
        -c ${droptag_xml} \
        ${bam}
    """
}
