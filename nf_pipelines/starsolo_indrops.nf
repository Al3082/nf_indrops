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
 * Preprocessing options:
 *   --preprocess_only      Stop after trim+sync, publish FASTQs to <output_dir>/processed_fastqs/
 *   --kotov_runs both|runA|runB   Which Kotov sequencing runs to process (default: both)
 *   --skip_briggs          Skip Briggs preprocessing (use when already processed in a prior run)
 *   --processed_fastq_dir  Skip all preprocessing and load pre-processed FASTQs from this directory
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

// Adapter fasta for gene-read trimming (cutadapt)
params.adapter_fasta = null  // required: e.g. illumina_nextseq_p7.fasta

// Optional QC steps
params.skip_preflight = false  // --skip_preflight to skip pre-flight FASTQ checks
params.skip_qc        = false  // --skip_qc to skip FastQC/MultiQC

// Preprocessing options
params.preprocess_only = false  // --preprocess_only to stop after trim+sync (publish FASTQs, skip alignment)
params.kotov_runs      = 'both' // --kotov_runs both|runA|runB — which Kotov sequencing runs to process
params.skip_briggs     = false  // --skip_briggs to skip Briggs preprocessing (use when Briggs was already processed)

// Test mode: run starsolo EditDist_2 + 1MM in parallel on a single library (optionally also dropest)
params.test           = false  // --test to enable test mode
params.test_library   = 'S11_1_1'  // library to test on
params.test_dropest   = false  // --test_dropest to also run dropest alongside starsolo in test mode
params.merged_fastq_dir    = null  // --merged_fastq_dir to skip merge and use pre-merged FASTQs for dropest test
params.processed_fastq_dir = null  // --processed_fastq_dir to skip preprocessing and load pre-processed FASTQs directly

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

include { demux_libraries; extract_reads; sync_reads as sync_kotov; sync_reads as sync_briggs; sync_single_read as sync_v2_barcode } from './modules/demux.nf'
include { consolidate_stats                              } from './modules/stats.nf'
include { starsolo_v3; starsolo_v3_1mm; starsolo_v2      } from './modules/starsolo.nf'
include { droptag_v3; star_plain as star_plain_v3; dropest_quant as dropest_v3 } from './modules/dropest.nf'
include { droptag_v2; star_plain as star_plain_v2; dropest_quant as dropest_v2 } from './modules/dropest.nf'
include { droptag_v3 as droptag_v3_test; droptag_v3_premerged as droptag_v3_premerged_test; star_plain as star_plain_test; dropest_quant as dropest_test } from './modules/dropest.nf'
include { trim_gene_read as trim_kotov_r1                } from './modules/trim.nf'
include { trim_gene_read as trim_briggs_r1               } from './modules/trim.nf'
include { trim_gene_read as trim_v2_gene                 } from './modules/trim.nf'
include { check_fastq                                    } from './modules/check_fastq.nf'
include { fastqc; multiqc                                } from './modules/qc.nf'


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

        if (params.processed_fastq_dir) {

            // ── Skip preprocessing: load pre-processed FASTQs directly ───────
            //    Expected naming in processed_fastq_dir/:
            //      {lib}__{run}.trimmed_R1.fastq.gz   (Kotov gene reads)
            //      {lib}__{run}.synced_R2.fastq.gz    (Kotov CB1)
            //      {lib}__{run}.synced_R4.fastq.gz    (Kotov CB2+UMI)
            //      {lib}.trimmed_R1.fastq.gz          (Briggs gene reads)
            //      {lib}.synced_R2.fastq.gz           (Briggs CB1)
            //      {lib}.synced_R4.fastq.gz           (Briggs CB2+UMI)

            pdir = params.processed_fastq_dir
            if (!file(pdir).isDirectory()) error "--processed_fastq_dir '${pdir}' does not exist or is not a directory"

            r1_ch = Channel.fromPath("${pdir}/*.trimmed_R1.fastq.gz", checkIfExists: true)
                .map { f ->
                    def base = f.name.replace('.trimmed_R1.fastq.gz', '')
                    def lib = base.contains('__') ? base.split('__')[0] : base
                    tuple(lib, f)
                }
                .groupTuple()

            r2_ch = Channel.fromPath("${pdir}/*.synced_R2.fastq.gz")
                .map { f ->
                    def base = f.name.replace('.synced_R2.fastq.gz', '')
                    def lib = base.contains('__') ? base.split('__')[0] : base
                    tuple(lib, f)
                }
                .groupTuple()

            r4_ch = Channel.fromPath("${pdir}/*.synced_R4.fastq.gz")
                .map { f ->
                    def base = f.name.replace('.synced_R4.fastq.gz', '')
                    def lib = base.contains('__') ? base.split('__')[0] : base
                    tuple(lib, f)
                }
                .groupTuple()

            starsolo_in = r1_ch.join(r2_ch).join(r4_ch)
                .map { lib, r1s, r2s, r4s ->
                    tuple(lib,
                          r1s.sort { it.name }.collect { it.toString() },
                          r2s.sort { it.name }.collect { it.toString() },
                          r4s.sort { it.name }.collect { it.toString() })
                }

        } else {

        // ── 0. Pre-flight checks on raw Kotov and first Briggs file ───────────

        if (params.kotov_runs == 'runA') {
            kotov_check_ch = Channel.of(tuple('kotov_runA', runA[0], runA[1], runA[3]))
        } else if (params.kotov_runs == 'runB') {
            kotov_check_ch = Channel.of(tuple('kotov_runB', runB[0], runB[1], runB[3]))
        } else {
            kotov_check_ch = Channel.of(
                tuple('kotov_runA', runA[0], runA[1], runA[3]),
                tuple('kotov_runB', runB[0], runB[1], runB[3])
            )
        }

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

        if (!params.skip_preflight) {
            preflight_ch = params.skip_briggs
                ? kotov_check_ch
                : kotov_check_ch.mix(
                    briggs_ch.first().map { lib, r1, r2, r4 -> tuple("briggs_${lib}", r1, r2, r4) }
                )
            check_fastq(
                preflight_ch,
                file(params.cb_whitelist1),
                file(params.cb_whitelist2_sense),
                file(params.cb_whitelist2_rc)
            )
        }

        // ── 1. Cutadapt demultiplexing of Kotov R3 reads ──────────────────────
        //    Input:  (run_id, r3_fastq)
        //    Output: (run_id, [lib1.fastq.gz, lib2.fastq.gz, ...])

        if (params.kotov_runs == 'runA') {
            kotov_runs_r3 = Channel.of(tuple('runA', runA[2]))
        } else if (params.kotov_runs == 'runB') {
            kotov_runs_r3 = Channel.of(tuple('runB', runB[2]))
        } else {
            kotov_runs_r3 = Channel.of(
                tuple('runA', runA[2]),
                tuple('runB', runB[2])
            )
        }

        demux_out = demux_libraries(kotov_runs_r3, barcodes_fa)

        // Flatten: one tuple per (run_id, library) -> (run_id, lib_name, demuxed_file)
        demux_flat = demux_out.demuxed
            .transpose()
            .map { run_id, f -> tuple(run_id, f.simpleName, f) }

        // ── 2. Build channel of Kotov R1/R2/R4 per run ───────────────────────

        if (params.kotov_runs == 'runA') {
            kotov_reads = Channel.of(tuple('runA', runA[0], runA[1], runA[3]))
        } else if (params.kotov_runs == 'runB') {
            kotov_reads = Channel.of(tuple('runB', runB[0], runB[1], runB[3]))
        } else {
            kotov_reads = Channel.of(
                tuple('runA', runA[0], runA[1], runA[3]),
                tuple('runB', runB[0], runB[1], runB[3])
            )
        }

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

        extract_out = extract_reads(extract_in)
        extracted = extract_out.reads

        // ── 3b. Trim Kotov R1 (gene reads) ───────────────────────────────────
        //    cutadapt -m 20 discards short reads, breaking sync with R2/R4.
        //    sync_kotov re-filters R2/R4 to match surviving R1 read IDs.

        adapter_fa = file(params.adapter_fasta, checkIfExists: true)

        kotov_trim_in = extracted.map { lib, run_id, r1, r2, r4 ->
            tuple("${lib}__${run_id}", r1)
        }
        trimmed_kotov = trim_kotov_r1(kotov_trim_in, adapter_fa).trimmed

        // Join trimmed R1 back with original R2/R4, then re-sync
        kotov_sync_in = extracted
            .map { lib, run_id, r1, r2, r4 -> tuple("${lib}__${run_id}", lib, run_id, r2, r4) }
            .join(trimmed_kotov)
            .map { key, lib, run_id, r2, r4, tr1 -> tuple("${lib}__${run_id}", tr1, r2, r4) }

        kotov_sync_out = sync_kotov(kotov_sync_in)
        kotov_synced = trimmed_kotov
            .join(kotov_sync_out.reads)
            .map { key, tr1, r2, r4 ->
                def (lib, run_id) = key.split('__')
                tuple(lib, run_id, tr1, r2, r4)
            }

        // Group both runs per library
        // -> (lib_name, [run_ids], [r1s], [r2s], [r4s])
        kotov_grouped = kotov_synced
            .groupTuple(by: 0)

        // ── 3c. Trim Briggs R1 (gene reads) ──────────────────────────────────

        if (!params.skip_briggs) {
            briggs_trim_in = briggs_ch.map { lib, r1, r2, r4 -> tuple(lib, r1) }
            trimmed_briggs = trim_briggs_r1(briggs_trim_in, adapter_fa).trimmed

            // Re-sync Briggs R2/R4 with trimmed R1
            briggs_sync_in = briggs_ch
                .map { lib, r1, r2, r4 -> tuple(lib, r2, r4) }
                .join(trimmed_briggs)
                .map { lib, r2, r4, tr1 -> tuple(lib, tr1, r2, r4) }

            briggs_sync_out = sync_briggs(briggs_sync_in)
            briggs_trimmed = trimmed_briggs
                .join(briggs_sync_out.reads)
                .map { lib, tr1, r2, r4 -> tuple(lib, tr1, r2, r4) }
        }

        // ── 3d. Consolidate read statistics ──────────────────────────────────

        stats_id = params.kotov_runs == 'both'
            ? params.batch
            : "${params.batch}_${params.kotov_runs}" + (params.skip_briggs ? '_kotov_only' : '')

        if (params.skip_briggs) {
            consolidate_stats(
                stats_id,
                demux_out.counts
                    .mix(extract_out.stats)
                    .mix(trim_kotov_r1.out.counts)
                    .mix(kotov_sync_out.stats)
                    .collect()
            )
        } else {
            consolidate_stats(
                stats_id,
                demux_out.counts
                    .mix(extract_out.stats)
                    .mix(trim_kotov_r1.out.counts)
                    .mix(trim_briggs_r1.out.counts)
                    .mix(kotov_sync_out.stats)
                    .mix(briggs_sync_out.stats)
                    .collect()
            )
        }

        // ── 4. Join Kotov (grouped) + Briggs per library ──────────────────────

        if (params.preprocess_only) {
            // Stop here — processed FASTQs are published via publishDir
        } else {
            starsolo_in = kotov_grouped
                .join(briggs_trimmed, by: 0)
                .map { lib_name, run_ids, kotov_r1s, kotov_r2s, kotov_r4s,
                       briggs_r1, briggs_r2, briggs_r4 ->
                    all_r1 = (kotov_r1s + [briggs_r1]).collect { it.toString() }
                    all_r2 = (kotov_r2s + [briggs_r2]).collect { it.toString() }
                    all_r4 = (kotov_r4s + [briggs_r4]).collect { it.toString() }
                    tuple(lib_name, all_r1, all_r2, all_r4)
                }
        }

        } // end else (preprocessing path)

        // ── 5. Alignment & quantification (shared by both paths) ─────────────
        //    starsolo_in is built either from --processed_fastq_dir or from
        //    the preprocessing pipeline above.

        if (!params.preprocess_only) {

            if (!params.skip_qc && !params.processed_fastq_dir) {
                fastqc(starsolo_in)
                multiqc(fastqc.out.reports.collect())
            }

            if (params.test) {
            // Test mode: run STARsolo EditDist_2 and 1MM on a single library
            test_ch = starsolo_in.filter { lib_name, r1, r2, r4 ->
                lib_name == params.test_library
            }

            // 1) STARsolo EditDist_2
            starsolo_v3(
                test_ch,
                params.genome_dir,
                params.cb_whitelist2_sense,
                params.cb_whitelist2_rc
            )

            // 2) STARsolo 1MM
            starsolo_v3_1mm(
                test_ch,
                params.genome_dir,
                params.cb_whitelist2_sense,
                params.cb_whitelist2_rc
            )

            // 3) dropEst (optional, enabled with --test_dropest)
            if (params.test_dropest) {
                droptag_xml = params.droptag_v3_xml
                gtf_file    = params.gtf

                if (params.merged_fastq_dir) {
                    premerged_ch = Channel.of(
                        tuple(
                            params.test_library,
                            file("${params.merged_fastq_dir}/merged_R1.fastq.gz", checkIfExists: true),
                            file("${params.merged_fastq_dir}/merged_R2.fastq.gz", checkIfExists: true),
                            file("${params.merged_fastq_dir}/merged_R4.fastq.gz", checkIfExists: true)
                        )
                    )
                    tagged  = droptag_v3_premerged_test(premerged_ch, droptag_xml)
                } else {
                    tagged  = droptag_v3_test(test_ch, droptag_xml)
                }
                aligned = star_plain_test(tagged, params.genome_dir)
                dropest_test(aligned, gtf_file, droptag_xml)
            }

        } else if (params.aligner == 'starsolo') {
            starsolo_v3(
                starsolo_in,
                params.genome_dir,
                params.cb_whitelist2_sense,
                params.cb_whitelist2_rc
            )
        } else if (params.aligner == 'dropest') {
            droptag_xml = params.droptag_v3_xml
            gtf_file    = params.gtf

            tagged  = droptag_v3(starsolo_in, droptag_xml)
            aligned = star_plain_v3(tagged, params.genome_dir)
            dropest_v3(aligned, gtf_file, droptag_xml)
        }
        } // end !preprocess_only
}


// ── Workflow: v2 ──────────────────────────────────────────────────────────────

workflow RUN_V2 {
    take:
        samplesheet    // path to CSV (library_name, srr)

    main:

        briggs_v2_raw = Channel.fromPath(samplesheet)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row.library_name,
                    fastq("${row.srr}_1.fastq.gz"),
                    fastq("${row.srr}_2.fastq.gz")
                )
            }

        // ── Trim gene reads (R2) and re-sync barcode reads (R1) ─────────

        adapter_fa = file(params.adapter_fasta, checkIfExists: true)

        // Trim R2 (gene read): adapters, poly-A, low-quality bases
        v2_trim_in = briggs_v2_raw.map { lib, r1, r2 -> tuple(lib, r2) }
        trimmed_v2 = trim_v2_gene(v2_trim_in, adapter_fa).trimmed

        // Sync R1 (barcode read) with surviving trimmed R2 read IDs
        v2_sync_in = briggs_v2_raw
            .map { lib, r1, r2 -> tuple(lib, r1) }
            .join(trimmed_v2)
            .map { lib, r1, tr2 -> tuple(lib, tr2, r1) }

        synced_v2 = sync_v2_barcode(v2_sync_in)

        // Reassemble (lib, r1_barcode, r2_gene) with trimmed/synced reads
        briggs_v2 = synced_v2.reads
            .join(trimmed_v2)
            .map { lib, sr1, tr2 -> tuple(lib, sr1, tr2) }

        // ── Consolidate read statistics ──────────────────────────────────

        consolidate_stats(
            params.batch,
            trim_v2_gene.out.counts
                .mix(synced_v2.stats)
                .collect()
        )

        // ── Alignment ────────────────────────────────────────────────────

        if (params.preprocess_only) {
            // Stop here — processed FASTQs are published via publishDir
        } else if (params.aligner == 'starsolo') {
            starsolo_v2(
                briggs_v2,
                params.genome_dir,
                params.cb_whitelist1_rc,
                params.cb_whitelist2_sense
            )
        } else if (params.aligner == 'dropest') {
            droptag_xml = params.droptag_v2_xml
            gtf_file    = params.gtf

            tagged  = droptag_v2(briggs_v2, droptag_xml)
            aligned = star_plain_v2(tagged, params.genome_dir)
            dropest_v2(aligned, gtf_file, droptag_xml)
        }
}


// ── Entry point ───────────────────────────────────────────────────────────────

workflow {

    if (!params.batch) error "Please specify --batch (v3_s8_14 | v3_s16_22 | v2)"
    if (!params.processed_fastq_dir && !params.fastq_dir) error "Please specify --fastq_dir (or --processed_fastq_dir to skip preprocessing)"
    if (!params.output_dir) error "Please specify --output_dir"
    if (!params.preprocess_only && !params.genome_dir) error "Please specify --genome_dir"
    if (params.kotov_runs != 'both' && params.kotov_runs != 'runA' && params.kotov_runs != 'runB') {
        error "Unknown --kotov_runs '${params.kotov_runs}'. Use: both | runA | runB"
    }
    if (params.batch == 'v2' && params.kotov_runs != 'both') {
        log.warn "--kotov_runs is ignored for v2 batch (no Kotov data)"
    }
    if (params.batch == 'v2' && params.processed_fastq_dir) {
        error "--processed_fastq_dir is not supported for v2 batch"
    }
    if (params.processed_fastq_dir && params.preprocess_only) {
        error "--processed_fastq_dir and --preprocess_only are mutually exclusive"
    }
    if (params.skip_briggs && !params.preprocess_only) {
        error "--skip_briggs can only be used with --preprocess_only"
    }
    if (!params.preprocess_only && params.aligner != 'starsolo' && params.aligner != 'dropest') {
        error "Unknown --aligner '${params.aligner}'. Use: starsolo | dropest"
    }
    if (!params.preprocess_only && (params.aligner == 'dropest' || params.test_dropest)) {
        if (!params.dropest_container) error "Please specify --dropest_container when using --aligner dropest or --test_dropest"
        if (!params.gtf) error "Please specify --gtf (GTF annotation) when using --aligner dropest or --test_dropest"
    }
    if (!params.processed_fastq_dir && !params.adapter_fasta) {
        error "Please specify --adapter_fasta (e.g. illumina_nextseq_p7.fasta) for gene-read trimming"
    }
    if (params.batch == 'v2' && (!params.cb_whitelist1_rc || !params.cb_whitelist2_sense)) {
        error "Please specify --cb_whitelist1_rc and --cb_whitelist2_sense (v2 whitelists)"
    }
    if (params.batch != 'v2' && !params.processed_fastq_dir && (!params.cb_whitelist1 || !params.cb_whitelist2_sense || !params.cb_whitelist2_rc)) {
        error "Please specify --cb_whitelist1, --cb_whitelist2_sense and --cb_whitelist2_rc (v3 whitelists)"
    }
    if (params.batch != 'v2' && params.processed_fastq_dir && (!params.cb_whitelist2_sense || !params.cb_whitelist2_rc)) {
        error "Please specify --cb_whitelist2_sense and --cb_whitelist2_rc (needed for alignment)"
    }

    if (params.batch == 'v3_s8_14') {

        samplesheet = params.samplesheet
            ?: "${projectDir}/../inputs/samplesheet_v3_s8_14.csv"

        if (params.processed_fastq_dir) {
            // Skip raw FASTQ resolution — preprocessing will be bypassed
            RUN_V3(file(samplesheet), file(params.barcodes_v3_s8_14), [null,null,null,null], [null,null,null,null])
        } else {
            RUN_V3(
                file(samplesheet),
                file(params.barcodes_v3_s8_14),
                [fastq(params.kotov_s8_14_runA_r1), fastq(params.kotov_s8_14_runA_r2),
                 fastq(params.kotov_s8_14_runA_r3), fastq(params.kotov_s8_14_runA_r4)],
                [fastq(params.kotov_s8_14_runB_r1), fastq(params.kotov_s8_14_runB_r2),
                 fastq(params.kotov_s8_14_runB_r3), fastq(params.kotov_s8_14_runB_r4)]
            )
        }

    } else if (params.batch == 'v3_s16_22') {

        samplesheet = params.samplesheet
            ?: "${projectDir}/../inputs/samplesheet_v3_s16_22.csv"

        if (params.processed_fastq_dir) {
            RUN_V3(file(samplesheet), file(params.barcodes_v3_s16_22), [null,null,null,null], [null,null,null,null])
        } else {
            RUN_V3(
                file(samplesheet),
                file(params.barcodes_v3_s16_22),
                [fastq(params.kotov_s16_22_runA_r1), fastq(params.kotov_s16_22_runA_r2),
                 fastq(params.kotov_s16_22_runA_r3), fastq(params.kotov_s16_22_runA_r4)],
                [fastq(params.kotov_s16_22_runB_r1), fastq(params.kotov_s16_22_runB_r2),
                 fastq(params.kotov_s16_22_runB_r3), fastq(params.kotov_s16_22_runB_r4)]
            )
        }

    } else if (params.batch == 'v2') {

        samplesheet = params.samplesheet
            ?: "${projectDir}/../inputs/samplesheet_v2.csv"

        RUN_V2(file(samplesheet))

    } else {
        error "Unknown --batch '${params.batch}'. Use: v3_s8_14 | v3_s16_22 | v2"
    }
}
