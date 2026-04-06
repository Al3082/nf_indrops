# nf_cutadapt_starsolo

Nextflow DSL2 pipeline for processing **inDrop v2/v3 single-cell RNA-seq data** from *Xenopus tropicalis* developmental samples.

Combines two public datasets:
- **Kotov** (SRP363793) — Raw multiplexed reads from stages 8-22, requires demultiplexing
- **Briggs** (SRP142544) — Pre-demultiplexed reference libraries

## Pipeline overview

### v3 workflow (batches `v3_s8_14`, `v3_s16_22`)

```
Kotov R3 (library barcode)
  |-> cutadapt demux (pool1.fa / pool2.fa, 1 mismatch)
  |-> seqtk: extract matching R1/R2/R4 read IDs
  |-> cutadapt: trim adapters, poly-A, low-quality 3' bases from R1 (gene read)
  |-> seqtk: re-sync R2/R4 with surviving R1 reads
  |-> group by library (runA + runB)
  |-> join with Briggs R1/R2/R4 (also trimmed + synced)
  |-> STARsolo (CB_UMI_Complex) or dropEst
```

### v2 workflow (batch `v2`)

```
Briggs R1 (barcode) + R2 (gene read)
  |-> cutadapt: trim adapters, poly-A, low-quality 3' bases from R2 (gene read)
  |-> seqtk: re-sync R1 with surviving R2 reads
  |-> STARsolo (CB_UMI_Complex, adapter-delimited) or dropEst
```

### Read structures

| Version | Read | Length | Content |
|---------|------|--------|---------|
| v3 | R1 | 95 bp | cDNA (gene read) |
| v3 | R2 | 8 bp | Cell barcode part 1 (CB1) |
| v3 | R3 | 8 bp | Library barcode (demux only) |
| v3 | R4 | 14 bp | CB2 (8 bp) + UMI (6 bp) |
| v2 | R1 | variable | CB1 \| adapter \| CB2 (8 bp) \| UMI (6 bp) |
| v2 | R2 | variable | cDNA (gene read) |

## Batches

Processing is split into three batches due to HPC disk constraints:

| Batch | Kotov sample | Briggs libraries | Aligner options |
|-------|-------------|------------------|-----------------|
| `v3_s8_14` | GSM5949469 (stages 8-14) | 21 v3 libraries (Pool 2) | starsolo, dropest |
| `v3_s16_22` | GSM5949468 (stages 16-22) | 23 v3 libraries (Pool 1) | starsolo, dropest |
| `v2` | none | 31 v2 libraries (Pools 1+2) | starsolo, dropest |

## Requirements

- **Nextflow** >= 21.04 (DSL2)
- **Singularity** (for containerised tools)
- **Conda/Mamba** (for the `pystats` stats consolidation step; the `hpc` profile expects miniforge at `/mnt/beegfs/common/apps/miniforge/25.3.1/`)
- **STAR 2.7.11b** (loaded via HPC module or available in PATH)

### Containers (pulled automatically)

| Tool | Version | Container |
|------|---------|-----------|
| Cutadapt | 5.2 | `quay.io/biocontainers/cutadapt:5.2--py310h1fe012e_0` |
| Seqtk | 1.5 | `quay.io/biocontainers/seqtk:1.5--h577a1d6_1` |
| FastQC | 0.12.1 | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| MultiQC | 1.33 | `quay.io/biocontainers/multiqc:1.33--pyhdfd78af_0` |

dropEst requires a separate container specified via `--dropest_container`.

## Usage

All commands are run from `nf_pipelines/`.

### v3 batches

```bash
nextflow run starsolo_indrops.nf \
  --batch v3_s8_14 \
  --fastq_dir /path/to/fastq_files \
  --output_dir /path/to/output \
  --genome_dir /path/to/star_index \
  -profile hpc
```

### v2 batch

```bash
nextflow run starsolo_indrops.nf \
  --batch v2 \
  --fastq_dir /path/to/fastq_files \
  --output_dir /path/to/output \
  --genome_dir /path/to/star_index \
  -profile hpc
```

Cell barcode whitelists, adapter FASTA, and dropEst XML configs default to bundled files under `inputs/`. Only `--batch`, `--fastq_dir`, `--output_dir`, and `--genome_dir` are required.

### dropEst aligner

```bash
nextflow run starsolo_indrops.nf \
  --batch v3_s8_14 \
  --aligner dropest \
  --dropest_container /path/to/dropest.sif \
  --gtf /path/to/annotation.gtf \
  [other flags] \
  -profile hpc
```

### Optional flags

| Flag | Default | Description |
|------|---------|-------------|
| `--aligner` | `starsolo` | Aligner: `starsolo` or `dropest` |
| `--skip_preflight` | `false` | Skip pre-flight FASTQ sanity checks |
| `--skip_qc` | `false` | Skip FastQC/MultiQC |
| `--test` | `false` | Test mode: run all aligners on a single library |
| `--test_library` | `S11_1_1` | Library to use in test mode |
| `--samplesheet` | per-batch default | Override the bundled samplesheet CSV |
| `--cb_whitelist1` | bundled | CB1 forward whitelist |
| `--cb_whitelist1_rc` | bundled | CB1 reverse-complement whitelist (v2) |
| `--cb_whitelist2_sense` | bundled | CB2 forward whitelist |
| `--cb_whitelist2_rc` | bundled | CB2 reverse-complement whitelist (v3) |
| `--adapter_fasta` | bundled | Adapter sequences for gene-read trimming |
| `--dropest_container` | none | Singularity image for dropEst (required with `--aligner dropest`) |
| `--gtf` | none | GTF annotation (required with `--aligner dropest`) |
| `-resume` | — | Resume a previously failed run |

### Profiles

| Profile | Executor | Container runtime |
|---------|----------|-------------------|
| `hpc` | SLURM (`dev` queue) | Singularity (auto-mounts) + Conda/Mamba |
| `local` | Local | Singularity |

## Project structure

```
nf_pipelines/
  starsolo_indrops.nf           # Entry point — routes to RUN_V3 or RUN_V2 workflow
  nextflow.config               # Parameters, resource labels, profiles
  modules/
    demux.nf                    # Cutadapt demux + seqtk read extraction + sync
    trim.nf                     # Adapter/poly-A/quality trimming (cutadapt)
    starsolo.nf                 # STARsolo alignment (v2, v3, v3_1mm)
    dropest.nf                  # dropTag + plain STAR + dropEst quantification
    check_fastq.nf              # Pre-flight FASTQ integrity & barcode match checks
    qc.nf                       # FastQC + MultiQC
    stats.nf                    # Consolidate per-step read count statistics

inputs/
  samplesheet_v3_s8_14.csv      # library_name -> briggs_srr (21 libraries)
  samplesheet_v3_s16_22.csv     # library_name -> briggs_srr (23 libraries)
  samplesheet_v2.csv            # library_name -> srr (31 libraries)
  pool1.fa / pool2.fa           # 8 bp library barcodes for cutadapt demux
  gel_barcode1_list.txt         # CB1 whitelist (384 x 8 bp)
  gel_barcode1_rc_list.txt      # CB1 whitelist reverse complement
  gel_barcode2_list.txt         # CB2 whitelist forward (384 x 8 bp)
  bc2_rc.txt                    # CB2 whitelist reverse complement
  illumina_nextseq_p7.fasta     # Adapter sequences for gene-read trimming
  indrop_v3.xml / indrop_v1_2.xml  # dropEst/dropTag config files
  metadata_briggs.csv           # Briggs sample metadata
  metadata_kotov.csv            # Kotov sample metadata
```

## Resource allocations (per task)

| Tool | CPUs | RAM | Time | Retry behaviour |
|------|------|-----|------|-----------------|
| Cutadapt (demux + trim) | 20 | 12 GB | 2h | RAM & time scale with attempt |
| Seqtk (extract + sync) | 1 | 32 GB | 4h | RAM & time scale with attempt |
| STAR / STARsolo | 20 | 32 GB | 4h | RAM & time scale with attempt |
| dropEst | 20 | 80 GB | 8h | RAM & time scale with attempt |
| FastQC | 6 | 4 GB | 2h | RAM & time scale with attempt |
| MultiQC | 1 | 4 GB | 1h | RAM & time scale with attempt |
| Stats consolidation | 1 | 1 GB | 30m | — (conda Python env) |

Processes retry up to 2 times on SLURM OOM (exit 137), timeout (140), or segfault (139).

## Outputs

| Directory | Contents |
|-----------|----------|
| `<output_dir>/<library>/STAR/` | STARsolo count matrices (Gene, GeneFull, Velocyto), sorted BAM, logs |
| `<output_dir>/<library>/dropest/` | dropEst count matrices and cell statistics |
| `<output_dir>/cutadapt_stats/` | Per-run demultiplexing statistics |
| `<output_dir>/cutadapt_trim_stats/` | Per-library trimming statistics |
| `<output_dir>/pipeline_stats/` | Consolidated per-step read count summary TSV |
| `<output_dir>/preflight/` | Pre-flight FASTQ integrity reports |
| `<output_dir>/qc/` | MultiQC report and per-library FastQC reports |

## Notes

- The `hpc` profile has site-specific Singularity bind mounts, workDir, and miniforge paths that may need adjusting for other clusters.
- Conda environments are cached in `/mnt/beegfs/home/asalamit/persistent/nf_conda_envs/` (created automatically on first run).
- Spot-number extraction in `extract_reads` and `sync_reads` assumes SRA-style read names (`SRR*.N`). Non-SRA FASTQs would need adapter changes.
