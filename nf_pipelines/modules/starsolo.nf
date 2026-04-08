/*
 * STARsolo alignment for inDrop v3 libraries.
 *
 * Reads for each library are passed as lists (one entry per run/source):
 *   r1_files : gene reads  (R1 from Kotov runA, runB, Briggs)
 *   r2_files : CB part 1  (R2, 8 bp)
 *   r4_files : CB part 2 + UMI (R4, 14 bp: 8 bp CB2 + 6 bp UMI)
 *
 * R2 and R4 are merged into a single 22 bp barcode read before alignment:
 *   positions 0-7  : CB1 (from R2)
 *   positions 8-15 : CB2 (first 8 bp of R4)
 *   positions 16-21: UMI (last 6 bp of R4)
 *
 * STARsolo position notation (anchor_pos_anchor_pos, anchor 0 = read start):
 *   CB1 : 0_0_0_7
 *   CB2 : 0_8_0_15
 *   UMI : 0_16_0_21
 */
process starsolo_v3 {
    tag "STARsolo_v3 on ${lib_name}"
    label "star"

    publishDir "${params.output_dir}/${lib_name}", mode: 'copy'

    input:
        tuple val(lib_name), val(r1_files), val(r2_files), val(r4_files)
        val genome_dir
        val cb_whitelist2_sense
        val cb_whitelist2_rc

    output:
        path "STAR/*"

    script:
    r1_list = r1_files instanceof List ? r1_files : [r1_files]
    r2_list = r2_files instanceof List ? r2_files : [r2_files]
    r4_list = r4_files instanceof List ? r4_files : [r4_files]

    merge_cmds = [r2_list, r4_list].transpose().withIndex().collect { pair, i ->
        "paste <(zcat ${pair[0]}) <(zcat ${pair[1]}) | awk 'NR%4==1{print \$1} NR%4==2{print \$1\$2} NR%4==3{print \"+\"} NR%4==0{print \$1\$2}' | gzip -c > bc_${i}.fastq.gz"
    }.join('\n')

    r1_str = r1_list.join(',')
    bc_str = (0..<r2_list.size()).collect { "bc_${it}.fastq.gz" }.join(',')

    """
    set -euo pipefail

    mkdir -p STAR

    ${merge_cmds}

    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${r1_str} ${bc_str} \
        --readFilesCommand zcat \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype EditDist_2 \
        --soloCBwhitelist ${cb_whitelist2_sense} ${cb_whitelist2_rc} \
        --soloCBposition 0_0_0_7 0_8_0_15 \
        --soloUMIposition 0_16_0_21 \
        --soloFeatures Gene GeneFull Velocyto \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM CR CY UR UY CB UB sM \
        --soloMultiMappers EM \
        --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix STAR/ \
        --soloOutFormatFeaturesGeneField3 gene_name
    """
}


/*
 * STARsolo alignment for inDrop v3 libraries — 1MM barcode matching variant.
 *
 * Same as starsolo_v3 but uses --soloCBmatchWLtype 1MM instead of EditDist_2.
 */
process starsolo_v3_1mm {
    tag "STARsolo_v3_1mm on ${lib_name}"
    label "star"

    publishDir "${params.output_dir}/${lib_name}/starsolo_1mm", mode: 'copy'

    input:
        tuple val(lib_name), val(r1_files), val(r2_files), val(r4_files)
        val genome_dir
        val cb_whitelist2_sense
        val cb_whitelist2_rc

    output:
        path "STAR/*"

    script:
    r1_list = r1_files instanceof List ? r1_files : [r1_files]
    r2_list = r2_files instanceof List ? r2_files : [r2_files]
    r4_list = r4_files instanceof List ? r4_files : [r4_files]

    merge_cmds = [r2_list, r4_list].transpose().withIndex().collect { pair, i ->
        "paste <(zcat ${pair[0]}) <(zcat ${pair[1]}) | awk 'NR%4==1{print \$1} NR%4==2{print \$1\$2} NR%4==3{print \"+\"} NR%4==0{print \$1\$2}' | gzip -c > bc_${i}.fastq.gz"
    }.join('\n')

    r1_str = r1_list.join(',')
    bc_str = (0..<r2_list.size()).collect { "bc_${it}.fastq.gz" }.join(',')

    """
    set -euo pipefail

    mkdir -p STAR

    ${merge_cmds}

    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${r1_str} ${bc_str} \
        --readFilesCommand zcat \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype 1MM \
        --soloCBwhitelist ${cb_whitelist2_sense} ${cb_whitelist2_rc} \
        --soloCBposition 0_0_0_7 0_8_0_15 \
        --soloUMIposition 0_16_0_21 \
        --soloFeatures Gene GeneFull Velocyto \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM CR CY UR UY CB UB sM \
        --soloMultiMappers EM \
        --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix STAR/ \
        --soloOutFormatFeaturesGeneField3 gene_name
    """
}


/*
 * STARsolo alignment for inDrop v2 libraries (Briggs only, no Kotov counterpart).
 *
 * Read structure:
 *   reads[0] = R1: barcode read (CB1 variable | adapter | CB2 8 bp | UMI 6 bp)
 *   reads[1] = R2: gene read
 *
 * Adapter GAGTGATTGCTTGTGACGCCTT is used to locate CB1 end.
 */
process starsolo_v2 {
    tag "STARsolo_v2 on ${lib_name}"
    label "star"

    publishDir "${params.output_dir}/${lib_name}", mode: 'copy'

    input:
        tuple val(lib_name), path(r1), path(r2)
        val genome_dir
        val cb_whitelist1_rc
        val cb_whitelist2_sense

    output:
        path "STAR/*"

    script:
    """
    set -euo pipefail

    mkdir -p STAR

    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${r2} ${r1} \
        --readFilesCommand zcat \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype EditDist_2 \
        --soloCBwhitelist ${cb_whitelist1_rc} ${cb_whitelist2_sense} \
        --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT \
        --soloAdapterMismatchesNmax 2 \
        --soloCBposition 0_0_2_-1 3_1_3_8 \
        --soloUMIposition 3_9_3_14 \
        --soloFeatures Gene GeneFull Velocyto \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM CR CY UR UY CB UB sM \
        --soloMultiMappers EM \
        --limitBAMsortRAM 30000000000 \
        --outFileNamePrefix STAR/
    """
}
