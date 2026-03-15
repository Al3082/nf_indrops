/*
 * starsolo_indrops.nf
 *
 * Processes inDrop v2/v3 scRNA-seq data from Briggs (SRP142544) and Kotov (SRP363793)
 * using cutadapt demultiplexing (v3 only) + STARsolo alignment.
 *
 * Usage:
 *   nextflow run starsolo_indrops.nf --batch v3_s8_14  [options]
 *   nextflow run starsolo_indrops.nf --batch v3_s16_22 [options]
 *   nextflow run starsolo_indrops.nf --batch v2        [options]
 *
 * Batches:
 *   v3_s8_14  : Kotov GSM5949469 (stages 8-14)  + Briggs Pool 3 (21 v3 libraries)
 *   v3_s16_22 : Kotov GSM5949468 (stages 16-22) + Briggs Pool 4 (23 v3 libraries)
 *   v2        : Briggs Pool 1+2 (31 v2 libraries, no Kotov counterpart)
 *
 * Kotov read structure (4 reads × 2 sequencing runs per batch, stored as individual SRRs):
 *   R1 (95 bp) : gene read
 *   R2 ( 8 bp) : cell barcode part 1
 *   R3 ( 8 bp) : library barcode  ->  used for cutadapt demultiplexing
 *   R4 (14 bp) : cell barcode part 2 (8 bp) + UMI (6 bp)
 *
 * Briggs v3 read structure (4 files per SRR, already demultiplexed):
 *   _1 (R1) : gene read
 *   _2 (R2) : cell barcode part 1
 *   _3 (R3) : library barcode (demux already done)
 *   _4 (R4) : cell barcode part 2 + UMI
 *
 * Briggs v2 read structure (2 files per SRR, already demultiplexed):
 *   _1 (R1) : barcode read (CB1 | adapter | CB2 | UMI)
 *   _2 (R2) : gene read
 */

nextflow.enable.dsl = 2

// ── Parameters ────────────────────────────────────────────────────────────────

params.batch      = null    // required: v3_s8_14 | v3_s16_22 | v2
params.fastq_dir  = null    // required: directory containing all fastq.gz files
params.output_dir = null    // required: output root directory
params.genome_dir = null    // required: STAR genome index directory

// Cell barcode whitelists
params.cb_whitelist1       = null  // required: CB1 whitelist (gel_barcode1_list.txt) – v2 & v3
params.cb_whitelist2_sense = null  // required: CB2 forward whitelist (gel_barcode2_list.txt) – v3
params.cb_whitelist2_rc    = null  // required: CB2 reverse-complement whitelist (bc2_rc.txt) – v3

// Samplesheets  (defaults point to bundled files)
params.samplesheet = null   // overrides default per-batch samplesheet

// Batch-specific defaults (can be overridden)
params.barcodes_v3_s8_14  = "${projectDir}/../inputs/pool2.fa"
params.barcodes_v3_s16_22 = "${projectDir}/../inputs/pool1.fa"

// Kotov SRR accessions — v3_s8_14 (GSM5949469)
params.kotov_s8_14_runA_r1 = "SRR18313232.fastq.gz"
params.kotov_s8_14_runA_r2 = "SRR18313233.fastq.gz"
params.kotov_s8_14_runA_r3 = "SRR18313234.fastq.gz"
params.kotov_s8_14_runA_r4 = "SRR18313235.fastq.gz"
params.kotov_s8_14_runB_r1 = "SRR18313236.fastq.gz"
params.kotov_s8_14_runB_r2 = "SRR18313237.fastq.gz"
params.kotov_s8_14_runB_r3 = "SRR18313238.fastq.gz"
params.kotov_s8_14_runB_r4 = "SRR18313239.fastq.gz"

// Kotov SRR accessions — v3_s16_22 (GSM5949468)
params.kotov_s16_22_runA_r1 = "SRR18313240.fastq.gz"
params.kotov_s16_22_runA_r2 = "SRR18313241.fastq.gz"
params.kotov_s16_22_runA_r3 = "SRR18313242.fastq.gz"
params.kotov_s16_22_runA_r4 = "SRR18313243.fastq.gz"
params.kotov_s16_22_runB_r1 = "SRR18313244.fastq.gz"
params.kotov_s16_22_runB_r2 = "SRR18313245.fastq.gz"
params.kotov_s16_22_runB_r3 = "SRR18313246.fastq.gz"
params.kotov_s16_22_runB_r4 = "SRR18313247.fastq.gz"



// ── Imports ───────────────────────────────────────────────────────────────────

include { demux_libraries; extract_reads } from './modules/demux.nf'
include { starsolo_v3; starsolo_v2       } from './modules/starsolo.nf'
include { check_fastq                    } from './modules/check_fastq.nf'
include { fastqc; multiqc                } from './modules/qc.nf'


// ── Helper: resolve fastq path ────────────────────────────────────────────────

def fastq(filename) {
    return file("${params.fastq_dir}/${filename}", checkIfExists: true)
}


// ── Workflow: v3 (shared for both v3 batches) ─────────────────────────────────

workflow RUN_V3 {
    take:
        samplesheet    // path to CSV (library_name, briggs_srr)
        barcodes_fa    // path to pool fasta for cutadapt demux
        runA           // tuple: [r1, r2, r3, r4] paths for Kotov run A
        runB           // tuple: [r1, r2, r3, r4] paths for Kotov run B

    main:

        // ── 0. Pre-flight checks on raw Kotov and first Briggs file ───────────

        kotov_check_ch = Channel.of(
            tuple('kotov_runA', runA[0], runA[1], runA[3]),
            tuple('kotov_runB', runB[0], runB[1], runB[3])
        )

        // ── 4. Load Briggs samplesheet, build file paths ──────────────────────

        briggs_ch = Channel.fromPath(samplesheet)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row.library_name,
                    fastq("${row.briggs_srr}_1.fastq.gz"),
                    fastq("${row.briggs_srr}_2.fastq.gz"),
                    fastq("${row.briggs_srr}_4.fastq.gz")
                )
            }

        check_fastq(
            kotov_check_ch.mix(
                briggs_ch.first().map { lib, r1, r2, r4 -> tuple("briggs_${lib}", r1, r2, r4) }
            ),
            file(params.cb_whitelist1),
            file(params.cb_whitelist2_sense),
            file(params.cb_whitelist2_rc)
        )

        // ── 1. Cutadapt demultiplexing of Kotov R3 reads ──────────────────────
        //    Input:  (run_id, r3_fastq)
        //    Output: (run_id, [lib1.fastq.gz, lib2.fastq.gz, ...])

        kotov_runs_r3 = Channel.of(
            tuple('runA', runA[2]),
            tuple('runB', runB[2])
        )

        demux_out = demux_libraries(kotov_runs_r3, barcodes_fa)

        // Flatten: one tuple per (run_id, library) -> (run_id, lib_name, demuxed_file)
        demux_flat = demux_out.demuxed
            .transpose()
            .map { run_id, f -> tuple(run_id, f.simpleName, f) }

        // ── 2. Build channel of Kotov R1/R2/R4 per run ───────────────────────

        kotov_reads = Channel.of(
            tuple('runA', runA[0], runA[1], runA[3]),
            tuple('runB', runB[0], runB[1], runB[3])
        )

        // Join demuxed IDs with the corresponding run's R1/R2/R4
        // demux_flat : (run_id, lib_name, demuxed_r3)
        // kotov_reads: (run_id, r1, r2, r4)
        // -> (lib_name, run_id, demuxed_r3, r1, r2, r4)

        extract_in = demux_flat
            .map { run_id, lib_name, demuxed_r3 -> tuple(run_id, lib_name, demuxed_r3) }
            .combine(kotov_reads, by: 0)
            .map { run_id, lib_name, demuxed_r3, r1, r2, r4 ->
                tuple(lib_name, run_id, demuxed_r3, r1, r2, r4)
            }

        // ── 3. Seqtk: extract R1/R2/R4 reads matching demuxed IDs ────────────
        //    Output: (lib_name, run_id, r1_out, r2_out, r4_out)

        extracted = extract_reads(extract_in)

        // Group both runs per library
        // -> (lib_name, [run_ids], [r1s], [r2s], [r4s])
        kotov_grouped = extracted
            .groupTuple(by: 0, size: 2)

        // ── 4. Join Kotov (grouped) + Briggs per library ──────────────────────

        starsolo_in = kotov_grouped
            .join(briggs_ch, by: 0)
            .map { lib_name, run_ids, kotov_r1s, kotov_r2s, kotov_r4s,
                   briggs_r1, briggs_r2, briggs_r4 ->
                all_r1 = (kotov_r1s + [briggs_r1]).collect { it.toString() }
                all_r2 = (kotov_r2s + [briggs_r2]).collect { it.toString() }
                all_r4 = (kotov_r4s + [briggs_r4]).collect { it.toString() }
                tuple(lib_name, all_r1, all_r2, all_r4)
            }

        // ── 5. FastQC on per-library reads ────────────────────────────────────

        fastqc(starsolo_in)

        multiqc(fastqc.out.reports.collect())

        // ── 6. STARsolo alignment ─────────────────────────────────────────────

        starsolo_v3(
            starsolo_in,
            params.genome_dir,
            params.cb_whitelist2_sense,
            params.cb_whitelist2_rc
        )
}


// ── Workflow: v2 ──────────────────────────────────────────────────────────────

workflow RUN_V2 {
    take:
        samplesheet    // path to CSV (library_name, srr)

    main:

        briggs_v2 = Channel.fromPath(samplesheet)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row.library_name,
                    fastq("${row.srr}_1.fastq.gz"),
                    fastq("${row.srr}_2.fastq.gz")
                )
            }

        starsolo_v2(
            briggs_v2,
            params.genome_dir,
            params.cb_whitelist1_rc,
            params.cb_whitelist2_sense
        )
}


// ── Entry point ───────────────────────────────────────────────────────────────

workflow {

    if (!params.batch) error "Please specify --batch (v3_s8_14 | v3_s16_22 | v2)"
    if (!params.fastq_dir)  error "Please specify --fastq_dir"
    if (!params.output_dir) error "Please specify --output_dir"
    if (!params.genome_dir) error "Please specify --genome_dir"
    if (params.batch == 'v2' && (!params.cb_whitelist1_rc || !params.cb_whitelist2_sense)) {
        error "Please specify --cb_whitelist1_rc and --cb_whitelist2_sense (v2 whitelists)"
    }
    if (params.batch != 'v2' && (!params.cb_whitelist1 || !params.cb_whitelist2_sense || !params.cb_whitelist2_rc)) {
        error "Please specify --cb_whitelist1, --cb_whitelist2_sense and --cb_whitelist2_rc (v3 whitelists)"
    }

    if (params.batch == 'v3_s8_14') {

        samplesheet = params.samplesheet
            ?: "${projectDir}/../inputs/samplesheet_v3_s8_14.csv"

        RUN_V3(
            file(samplesheet),
            file(params.barcodes_v3_s8_14),
            [fastq(params.kotov_s8_14_runA_r1), fastq(params.kotov_s8_14_runA_r2),
             fastq(params.kotov_s8_14_runA_r3), fastq(params.kotov_s8_14_runA_r4)],
            [fastq(params.kotov_s8_14_runB_r1), fastq(params.kotov_s8_14_runB_r2),
             fastq(params.kotov_s8_14_runB_r3), fastq(params.kotov_s8_14_runB_r4)]
        )

    } else if (params.batch == 'v3_s16_22') {

        samplesheet = params.samplesheet
            ?: "${projectDir}/../inputs/samplesheet_v3_s16_22.csv"

        RUN_V3(
            file(samplesheet),
            file(params.barcodes_v3_s16_22),
            [fastq(params.kotov_s16_22_runA_r1), fastq(params.kotov_s16_22_runA_r2),
             fastq(params.kotov_s16_22_runA_r3), fastq(params.kotov_s16_22_runA_r4)],
            [fastq(params.kotov_s16_22_runB_r1), fastq(params.kotov_s16_22_runB_r2),
             fastq(params.kotov_s16_22_runB_r3), fastq(params.kotov_s16_22_runB_r4)]
        )

    } else if (params.batch == 'v2') {

        samplesheet = params.samplesheet
            ?: "${projectDir}/../inputs/samplesheet_v2.csv"

        RUN_V2(file(samplesheet))

    } else {
        error "Unknown --batch '${params.batch}'. Use: v3_s8_14 | v3_s16_22 | v2"
    }
}
