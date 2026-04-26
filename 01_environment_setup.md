# Module 01 — Environment Setup & Reference Preparation

**HCS 7004 — Genome Analytics | The Ohio State University**

*← [Module 00: Conceptual Framework](00_conceptual_framework.md) | [Module 02: WGS Variant Calling](02_wgs_variant_calling.md) →*

---

> ### 🔑 Before you begin — username and environment variables
> Replace `<your_username>` with your OSC username in **every script** in this
> module. For example, if your username is `jsmith`, every `<your_username>`
> becomes `jsmith`.
>
> If you are starting here after a break or working through this module
> independently, set these variables in your terminal before running any
> interactive commands:
>
> ```bash
> export user_name=<your_username>
> export VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
> export SHARED=/fs/scratch/PAS3260/Variant_Calling/Data
> export CONTAINERS=/fs/scratch/PAS3260/Variant_Calling/Containers
> ```
>
> These four variables are referenced by every interactive command in this
> tutorial series. SLURM batch scripts re-declare all paths internally — you
> still need to replace `<your_username>` on the `user_name=` line inside each
> script. Setting the variables above only affects your interactive terminal.

---

## Learning Objectives

By the end of this module you will be able to:

1. Set up a reproducible personal working environment on OSC for the full
   tutorial series
2. Copy and inspect the *Arabidopsis thaliana* TAIR10 reference genome
3. Build the alignment indexes required by BWA-MEM2, STAR, and GATK4
4. Copy all raw sequencing data from the shared directory into your workspace
5. Perform raw read quality assessment with FastQC and MultiQC
6. Apply adapter trimming and quality filtering with fastp, correctly
   handling the different library types (paired-end WGS, single-end GBS,
   single-end RNA-seq) present in this tutorial

---

## Biological Background

Before any variant can be called, every sequencing read must be placed in its
genomic context by aligning it to a reference genome. The reference serves as
the coordinate system for the entire analysis — every SNP position you report
refers to a position in the reference, and every gene annotation you use to
interpret your variants is anchored to that same coordinate system.

For this tutorial we use the **TAIR10** assembly of *Arabidopsis thaliana*
(ecotype Columbia-0, Col-0). This assembly is chromosome-level (5 chromosomes
plus mitochondrion and chloroplast), nearly complete at ~135 Mb of nuclear
sequence, and supported by the most extensively curated plant gene annotation
in existence.

> **Why does reference quality matter for variant calling?** Errors, gaps, and
> misassemblies in the reference are directly inherited by your variant calls. A
> collapsed repeat causes reads from multiple loci to pile up at one position,
> creating false heterozygosity. Using the best available reference is the single
> most impactful quality decision in the entire pipeline.

---

## Overview of This Module

```
Step 1  OSC login and directory structure
Step 2  Acquire containers (copy from shared directory)
Step 3  Copy reference FASTA and GFF3
Step 4  Inspect reference chromosome structure
Step 5  Build alignment indexes
        ├── samtools fai index       (all modules)
        ├── GATK sequence dictionary (Modules 02 and 04)
        ├── BWA-MEM2 index           (Modules 02 and 03)
        └── STAR splice-aware index  (Module 04)
Step 6  Copy raw sequencing data
        ├── WGS  — 5 accessions, paired-end 101 bp
        ├── GBS  — 5 accessions, single-end 101 bp (ApeKI simulated)
        └── RNA-seq — 5 accessions, single-end 100 bp
Step 7  FastQC on all raw FASTQ files
Step 8  fastp adapter trimming
Step 9  MultiQC aggregate report
Step 10 Completion checklist
```

---

## Step 1 — Log In and Create Your Directory Structure

### 1.1 Connect to OSC

```bash
ssh <your_username>@owens.osc.edu
```

### 1.2 Set your username variable

> ⚠️ **Do this first, every session.** Replace `<your_username>` with your
> actual OSC username (e.g. `jsmith`). This variable is referenced by every
> command in this module.

```bash
export user_name=<your_username>
echo "Username is set to: ${user_name}"
```

### 1.3 Create the working directory tree

```bash
# Move to scratch — never run jobs from your home directory
cd /fs/scratch/PAS3260/

# Create the full working directory tree
mkdir -p /fs/scratch/PAS3260/${user_name}/variant_calling/{00_reference/{fasta,gff,indexes/{bwa,star}},01_raw_reads/{wgs,gbs,rnaseq},02_qc/{fastqc_raw,fastp},03_wgs/{aligned,gvcf,genotyped,filtered},04_gbs/{aligned,vcf},05_rnaseq/{aligned,vcf},06_annotation,07_applications,containers,logs,scripts}

# Verify
ls /fs/scratch/PAS3260/${user_name}/variant_calling/
```

### 1.4 Set persistent environment variables

The three variables below are used in every module. Add them to `.bashrc` so
they are available in every session and every SLURM job you submit.

```bash
cat >> ~/.bashrc << EOF

# ---- HCS 7004 Variant Calling Tutorial ----
export user_name=<your_username>
export VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
export SHARED=/fs/scratch/PAS3260/Variant_Calling/Data
export CONTAINERS=/fs/scratch/PAS3260/Variant_Calling/Containers
EOF

source ~/.bashrc

# Verify all three are set correctly
echo "user_name:  ${user_name}"
echo "VARCALL:    ${VARCALL}"
echo "SHARED:     ${SHARED}"
echo "CONTAINERS: ${CONTAINERS}"
```

> **Note the escaped `\${user_name}` inside the heredoc.** Without the
> backslash, bash expands `$user_name` immediately (to its current value)
> rather than storing the variable reference. The escaped form re-evaluates
> correctly every time you log in.

---

## Step 2 — Acquire Containers

All 16 containers are pre-pulled in the shared directory. Copy them to your
personal `containers/` directory so your scripts are self-contained and
independent of other students.

```bash
cp ${CONTAINERS}/*.sif ${VARCALL}/containers/

# Verify — expect 16 .sif files
ls -lh ${VARCALL}/containers/
echo "Container count: $(ls ${VARCALL}/containers/*.sif | wc -l)  (expected: 16)"
```

### Container inventory

| Container | Tool | Modules |
|---|---|---|
| `sratools_3.2.1.sif` | SRA Toolkit | 01 |
| `ncbi_datasets_18.24.sif` | NCBI Datasets CLI | 01 |
| `fastp_1.3.2.sif` | fastp (QC + trimming) | 01 |
| `fastqc_0.12.1.sif` | FastQC | 01 |
| `multiqc_1.34.sif` | MultiQC | 01 |
| `bwamem2_2.2.3.sif` | BWA-MEM2 | 01, 02, 03 |
| `samtools_1.23.1.sif` | SAMtools | 01, 02, 03, 04 |
| `star_2.7.11b.sif` | STAR | 01, 04 |
| `gatk4_4.6.2.0.sif` | GATK4 | 02, 03, 04 |
| `picard_3.4.0.sif` | Picard | 02, 04 |
| `stacks_2.68.sif` | Stacks | 03 |
| `bcftools_1.23.1.sif` | BCFtools | 02, 03, 04, 05 |
| `vcftools_0.1.17.sif` | VCFtools | 05, 06 |
| `snpeff_5.4.0c.sif` | SnpEff | 05 |
| `plink_1.90b7.sif` | PLINK | 05, 06 |
| `bedtools_2.31.1.sif` | BEDTools | 03 |

> **R is not containerised.** All R-based steps in Modules 05–06 use the
> native OSC installation. Load it when needed with:
> ```bash
> module load gcc/12.3.0 R/4.5.2
> ```

---

## Step 3 — Copy the Reference Genome and Annotation

The TAIR10 reference FASTA and GFF3 annotation are pre-staged in the shared
directory. Copy both to your personal reference directories.

```bash
# Copy reference FASTA
cp ${SHARED}/Reference/fasta/Athaliana_TAIR10.fasta \
   ${VARCALL}/00_reference/fasta/

# Copy GFF3 annotation
cp ${SHARED}/Reference/gff/Athaliana_TAIR10.gff3 \
   ${VARCALL}/00_reference/gff/

# Verify
ls -lh ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
ls -lh ${VARCALL}/00_reference/gff/Athaliana_TAIR10.gff3
```

> **Why copy rather than symlink?** Index files created by BWA-MEM2, STAR,
> and GATK4 are written alongside the FASTA. The shared directory is read-only,
> so indexing must happen in your personal copy. The FASTA is 116 MB — well
> within your scratch quota.

---

## Step 4 — Inspect the Reference

Always verify chromosome names and structure before building indexes. Mismatches
between FASTA chromosome names and GFF3 chromosome names silently break SnpEff
annotation in Module 05.

```bash
# Check chromosome names in FASTA
grep "^>" ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
```

Expected output (Ensembl Plants naming convention — bare numbers):

```
>1 dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF
>2 dna:chromosome chromosome:TAIR10:2:1:19698289:1 REF
>3 dna:chromosome chromosome:TAIR10:3:1:23459830:1 REF
>4 dna:chromosome chromosome:TAIR10:4:1:18585056:1 REF
>5 dna:chromosome chromosome:TAIR10:5:1:26975502:1 REF
>Mt dna:chromosome chromosome:TAIR10:Mt:1:366924:1 REF
>Pt dna:chromosome chromosome:TAIR10:Pt:1:154478:1 REF
```

Build the `.fai` index and confirm chromosome sizes:

```bash
apptainer exec \
  --bind ${VARCALL}/00_reference/fasta:${VARCALL}/00_reference/fasta \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools faidx ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta

cat ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta.fai
```

Expected `.fai` output:

```
1    30427671   54         60  61
2    19698289   30943144   60  61
3    23459830   50969786   60  61
4    18585056   74821736   60  61
5    26975502   93715842   60  61
Mt   366924     121182345  60  61
Pt   154478     121574413  60  61
```

Confirm GFF3 chromosome names match the FASTA:

```bash
grep -v "^#" ${VARCALL}/00_reference/gff/Athaliana_TAIR10.gff3 | \
  cut -f1 | sort -u
```

Both files must use the same names (`1`, `2`, `3`, `4`, `5`, `Mt`, `Pt`).

> ⚠️ **Chromosome naming is non-negotiable.** If your FASTA uses `1` and your
> GFF3 uses `Chr1`, SnpEff annotation will silently produce zero results.
> Verify now — it takes 30 seconds and saves hours later.

---

## Step 5 — Build Alignment Indexes

### 5.1 BWA-MEM2 index (for WGS and GBS alignment)

```bash
cat > ${VARCALL}/scripts/01a_bwa_index.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/bwa_index_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/bwa_index_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

echo "=== BWA-MEM2 index: $(date) ==="

IDX_DIR=${VARCALL}/00_reference/indexes/bwa
FASTA=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta

# Copy FASTA to index directory so all index files are co-located
cp ${FASTA}     ${IDX_DIR}/Athaliana_TAIR10.fasta
cp ${FASTA}.fai ${IDX_DIR}/Athaliana_TAIR10.fasta.fai

apptainer exec \
  --bind ${IDX_DIR}:${IDX_DIR} \
  ${VARCALL}/containers/bwamem2_2.2.3.sif \
  bwa-mem2 index ${IDX_DIR}/Athaliana_TAIR10.fasta

echo "=== BWA-MEM2 index complete: $(date) ==="
ls -lh ${IDX_DIR}/
EOF

sbatch ${VARCALL}/scripts/01a_bwa_index.sh
```

Expected index files after completion (7 total):

```
Athaliana_TAIR10.fasta
Athaliana_TAIR10.fasta.fai
Athaliana_TAIR10.fasta.0123
Athaliana_TAIR10.fasta.amb
Athaliana_TAIR10.fasta.ann
Athaliana_TAIR10.fasta.bwt.2bit.64
Athaliana_TAIR10.fasta.pac
```

### 5.2 STAR genome index (for RNA-seq alignment)

STAR requires the GFF3 at index time to encode known splice junctions, which
dramatically improves alignment accuracy across exon-exon boundaries.

```bash
cat > ${VARCALL}/scripts/01b_star_index.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --account=PAS3260
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/star_index_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/star_index_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

echo "=== STAR genome index: $(date) ==="

apptainer exec \
  --bind ${VARCALL}/00_reference:${VARCALL}/00_reference \
  ${VARCALL}/containers/star_2.7.11b.sif \
  STAR \
    --runMode genomeGenerate \
    --genomeDir ${VARCALL}/00_reference/indexes/star \
    --genomeFastaFiles ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta \
    --sjdbGTFfile ${VARCALL}/00_reference/gff/Athaliana_TAIR10.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbOverhang 99 \
    --genomeSAindexNbases 12 \
    --runThreadN 8

# --sjdbOverhang: read_length - 1  (RNA-seq reads are 100 bp → 99)
# --genomeSAindexNbases 12: appropriate for ~135 Mb genome
#   Formula: min(14, floor(log2(GenomeLength)/2 - 1))
# --sjdbGTFtagExonParentTranscript Parent: required for GFF3 format
#   (STAR defaults expect GTF 'gene_id'; GFF3 uses 'Parent' attribute)

echo "=== STAR index complete: $(date) ==="
ls -lh ${VARCALL}/00_reference/indexes/star/
EOF

sbatch ${VARCALL}/scripts/01b_star_index.sh
```

### 5.3 GATK sequence dictionary

GATK4 requires a `.dict` file alongside the `.fai` index. The `.dict` is a
SAM-format header that defines chromosome order for all GATK tools.

```bash
cat > ${VARCALL}/scripts/01c_seq_dict.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=seq_dict
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/seq_dict_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/seq_dict_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

echo "=== GATK sequence dictionary: $(date) ==="

REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta

apptainer exec \
  --bind ${VARCALL}/00_reference/fasta:${VARCALL}/00_reference/fasta \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk CreateSequenceDictionary \
    -R ${REF} \
    -O ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.dict

echo "=== Dictionary complete: $(date) ==="
ls -lh ${VARCALL}/00_reference/fasta/
EOF

sbatch ${VARCALL}/scripts/01c_seq_dict.sh
```

### 5.4 Verify all indexes

Run this after all three indexing jobs complete:

```bash
echo "--- samtools .fai ---"
ls -lh ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta.fai

echo "--- GATK .dict ---"
ls -lh ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.dict

echo "--- BWA-MEM2 index files (expect 7) ---"
ls ${VARCALL}/00_reference/indexes/bwa/ | wc -l

echo "--- STAR index (key file SA) ---"
ls -lh ${VARCALL}/00_reference/indexes/star/SA
```

---

## Step 6 — Copy Raw Sequencing Data

All raw data are pre-staged in the shared directory. Copy each track to your
personal working directories.

### 6.1 WGS reads — 5 accessions, paired-end 101 bp

Five geographically diverse accessions sequenced by the 1001 Genomes Consortium
(Illumina HiSeq 2000, PRJNA273563):

| Accession | Origin | WGS SRR | Coverage |
|---|---|---|---|
| Vas-0 | Västerås, Sweden | SRR1946455 | ~32× |
| Bez-9 | Besançon, France | SRR1946475 | ~25× |
| Gen-8 | Geneva, Switzerland | SRR1946461 | ~14× |
| Mah-6 | Mahón, Menorca, Spain | SRR1946459 | ~17× |
| Usa-0 | Utsunomiya, Japan | SRR1946457 | ~18× |

```bash
cp ${SHARED}/WGS/*.fastq.gz ${VARCALL}/01_raw_reads/wgs/

# Verify — expect 10 files (5 samples × R1 + R2)
ls -lh ${VARCALL}/01_raw_reads/wgs/
echo "WGS file count: $(ls ${VARCALL}/01_raw_reads/wgs/*.fastq.gz | wc -l)  (expected: 10)"
```

### 6.2 GBS reads — 5 accessions, single-end 101 bp (ApeKI simulated)

In silico ApeKI-digested reads derived from the WGS data. See Module 03 for
the full biological rationale and simulation methodology.

```bash
cp ${SHARED}/GBS/GBS_multiplexed.fastq.gz ${VARCALL}/01_raw_reads/gbs/
cp ${SHARED}/GBS/barcodes.tsv             ${VARCALL}/01_raw_reads/gbs/
```

### 6.3 RNA-seq reads — 5 accessions, single-end 100 bp

Transcriptome data from the 1001 Epigenomes Project (Illumina HiSeq 2500,
PRJNA319904). Multiple sequencing runs per accession were pre-merged into one
file per sample.

| Accession | SRX | File | Size |
|---|---|---|---|
| Vas-0 | SRX1735056 | `Vas-0_rnaseq.fastq.gz` | 763 MB |
| Bez-9 | SRX1735072 | `Bez-9_rnaseq.fastq.gz` | 2.6 GB |
| Gen-8 | SRX1735063 | `Gen-8_rnaseq.fastq.gz` | 3.2 GB |
| Mah-6 | SRX1735060 | `Mah-6_rnaseq.fastq.gz` | 2.9 GB |
| Usa-0 | SRX1735058 | `Usa-0_rnaseq.fastq.gz` | 3.1 GB |

```bash
cp ${SHARED}/RNA-seq/*.fastq.gz ${VARCALL}/01_raw_reads/rnaseq/

# Verify — expect 5 files
ls -lh ${VARCALL}/01_raw_reads/rnaseq/
echo "RNA-seq file count: $(ls ${VARCALL}/01_raw_reads/rnaseq/*.fastq.gz | wc -l)  (expected: 5)"
```

> ⚠️ **Three different library types in this tutorial.** WGS is paired-end;
> GBS and RNA-seq are single-end. Every downstream step — fastp trimming,
> alignment, GATK pre-processing — must be parameterised correctly for the
> library type. The scripts in each module are written for the correct type;
> pay attention to flags like `--in1`/`--in2` (paired-end) vs. `--in1` only
> (single-end) in fastp.

### 6.4 Confirm all raw data

```bash
echo "=== Raw data inventory ==="

echo ""
echo "WGS (expect 10 files):"
ls -lh ${VARCALL}/01_raw_reads/wgs/*.fastq.gz

echo ""
echo "GBS (expect 5 files):"
ls -lh ${VARCALL}/01_raw_reads/gbs/*.fastq.gz

echo ""
echo "RNA-seq (expect 5 files):"
ls -lh ${VARCALL}/01_raw_reads/rnaseq/*.fastq.gz

echo ""
TOTAL=$(find ${VARCALL}/01_raw_reads -name "*.fastq.gz" | wc -l)
echo "Total FASTQ files: ${TOTAL}  (expected: 20)"
```

---

## Step 7 — FastQC on All Raw Reads

Assess quality across all 20 FASTQ files before any trimming. We run FastQC
as a SLURM array — one task per file.

```bash
cat > ${VARCALL}/scripts/01d_fastqc_raw.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-20
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastqc_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastqc_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

echo "=== FastQC job ${SLURM_ARRAY_TASK_ID}: $(date) ==="

QC_DIR=${VARCALL}/02_qc/fastqc_raw

# Build sorted list of all raw FASTQ files across all three tracks
mapfile -t ALL_FILES < <(find \
  ${VARCALL}/01_raw_reads/wgs/ \
  ${VARCALL}/01_raw_reads/gbs/ \
  ${VARCALL}/01_raw_reads/rnaseq/ \
  -name "*.fastq.gz" | sort)

echo "Total FASTQ files found: ${#ALL_FILES[@]}  (expected: 20)"

FILE="${ALL_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "Processing: $(basename ${FILE})"

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/fastqc_0.12.1.sif \
  fastqc \
    --outdir ${QC_DIR} \
    --threads 8 \
    --extract \
    ${FILE}

echo "=== FastQC complete for $(basename ${FILE}): $(date) ==="
EOF

# Verify file count before submitting — adjust --array if needed
NFILES=$(find ${VARCALL}/01_raw_reads/ -name "*.fastq.gz" | wc -l)
echo "Found ${NFILES} FASTQ files. Array set to 1-20."
echo "If this differs, edit --array=1-${NFILES} before submitting."

sbatch ${VARCALL}/scripts/01d_fastqc_raw.sh
```

### Interpreting FastQC results by library type

Because this tutorial uses three different library types, some FastQC
"warnings" or "failures" are expected and do not indicate problems:

| FastQC module | WGS expectation | GBS expectation | RNA-seq expectation |
|---|---|---|---|
| Per base sequence quality | Green throughout | Green; slight 3' drop | Green; slight 3' drop |
| Sequence duplication | < 20% | **High (> 50%) — expected** | 20–50% — expected |
| Per sequence GC content | ~36%, smooth | May be non-normal | Non-normal (bimodal OK) |
| Overrepresented sequences | Rare | Cut-site sequences present | rRNA fragments possible |
| Adapter content | Near zero | Near zero | Near zero |

> **GBS duplication is not a quality problem.** Because every GBS library cuts
> the genome at the same ApeKI sites, reads from the same site across all
> individuals will be identical or nearly identical. FastQC reports this as
> high duplication — it is a biological property of the experiment, not a
> library preparation failure.

---

## Step 8 — Adapter Trimming with fastp

We trim all three tracks with fastp. The key difference is library type:
WGS uses paired-end mode (`--in1` + `--in2`); GBS and RNA-seq use
single-end mode (`--in1` only).

```bash
cat > ${VARCALL}/scripts/01e_fastp_trim.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-15
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastp_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastp_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
FASTP_DIR=${VARCALL}/02_qc/fastp

echo "=== fastp trim job ${SLURM_ARRAY_TASK_ID}: $(date) ==="

# Array layout:
#   Tasks  1– 5: WGS  (paired-end 101 bp)
#   Tasks  6–10: GBS  (single-end 101 bp)
#   Tasks 11–15: RNA-seq (single-end 100 bp)

WGS_SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
GBS_SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
RNA_SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)

TASK=${SLURM_ARRAY_TASK_ID}

if (( TASK <= 5 )); then
    # ---- WGS — paired-end ----
    SAMPLE=${WGS_SAMPLES[$((TASK - 1))]}
    R1=${VARCALL}/01_raw_reads/wgs/${SAMPLE}_R1.fastq.gz
    R2=${VARCALL}/01_raw_reads/wgs/${SAMPLE}_R2.fastq.gz
    echo "Track: WGS | Sample: ${SAMPLE}"

    apptainer exec \
      --bind ${VARCALL}:${VARCALL} \
      ${VARCALL}/containers/fastp_1.3.2.sif \
      fastp \
        --in1  ${R1} \
        --in2  ${R2} \
        --out1 ${FASTP_DIR}/${SAMPLE}_wgs_R1.fastq.gz \
        --out2 ${FASTP_DIR}/${SAMPLE}_wgs_R2.fastq.gz \
        --html ${FASTP_DIR}/${SAMPLE}_wgs_fastp.html \
        --json ${FASTP_DIR}/${SAMPLE}_wgs_fastp.json \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread 8 \
        --report_title "${SAMPLE} WGS"

elif (( TASK <= 10 )); then
    # ---- GBS — single-end ----
    SAMPLE=${GBS_SAMPLES[$((TASK - 6))]}
    SE=${VARCALL}/01_raw_reads/gbs/${SAMPLE}_GBS.fastq.gz
    echo "Track: GBS | Sample: ${SAMPLE}"

    apptainer exec \
      --bind ${VARCALL}:${VARCALL} \
      ${VARCALL}/containers/fastp_1.3.2.sif \
      fastp \
        --in1  ${SE} \
        --out1 ${FASTP_DIR}/${SAMPLE}_gbs.fastq.gz \
        --html ${FASTP_DIR}/${SAMPLE}_gbs_fastp.html \
        --json ${FASTP_DIR}/${SAMPLE}_gbs_fastp.json \
        --disable_adapter_trimming \
        --qualified_quality_phred 20 \
        --length_required 30 \
        --trim_front1 5 \
        --cut_right \
        --cut_right_mean_quality 20 \
        --thread 8 \
        --report_title "${SAMPLE} GBS"

else
    # ---- RNA-seq — single-end ----
    SAMPLE=${RNA_SAMPLES[$((TASK - 11))]}
    SE=${VARCALL}/01_raw_reads/rnaseq/${SAMPLE}_rnaseq.fastq.gz
    echo "Track: RNA-seq | Sample: ${SAMPLE}"

    apptainer exec \
      --bind ${VARCALL}:${VARCALL} \
      ${VARCALL}/containers/fastp_1.3.2.sif \
      fastp \
        --in1  ${SE} \
        --out1 ${FASTP_DIR}/${SAMPLE}_rnaseq.fastq.gz \
        --html ${FASTP_DIR}/${SAMPLE}_rnaseq_fastp.html \
        --json ${FASTP_DIR}/${SAMPLE}_rnaseq_fastp.json \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 30 \
        --thread 8 \
        --report_title "${SAMPLE} RNA-seq"
fi

echo "=== fastp complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/01e_fastp_trim.sh
```

### fastp parameter rationale by track

| Parameter | WGS | GBS | RNA-seq | Reason |
|---|---|---|---|---|
| `--detect_adapter_for_pe` | ✓ | — | ✓ (SE mode) | Auto-detect TruSeq/Nextera adapters |
| `--disable_adapter_trimming` | — | ✓ | — | GBS reads (simulated from WGS) have no adapter contamination |
| `--length_required` | 50 bp | 30 bp | 30 bp | GBS/RNA reads are informative even at shorter lengths |
| `--trim_front1 5` | — | ✓ | — | Removes cut-site artefact at 5' end of GBS reads |
| `--cut_right` | — | ✓ | — | Sliding-window 3' trimming; better for variable-quality GBS |
| `--thread` | 8 | 8 | 8 | Match to `--cpus-per-task` |

> **Expected read retention:** 93–98% for WGS; 90–97% for GBS; 88–95% for
> RNA-seq. Values below 85% should be investigated before proceeding.

---

## Step 9 — MultiQC Aggregate Report

Once all FastQC and fastp jobs complete, run MultiQC to produce a single
interactive HTML summary:

```bash
cat > ${VARCALL}/scripts/01f_multiqc.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/multiqc_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/multiqc_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

echo "=== MultiQC: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/multiqc_1.34.sif \
  multiqc \
    ${VARCALL}/02_qc/ \
    --outdir ${VARCALL}/02_qc/ \
    --filename multiqc_report \
    --title "HCS7004 Variant Calling — Raw Read QC (WGS + GBS + RNA-seq)"

echo "=== MultiQC complete: $(date) ==="
ls -lh ${VARCALL}/02_qc/multiqc_report.html
EOF

sbatch ${VARCALL}/scripts/01f_multiqc.sh
```

---

## Step 10 — Completion Checklist

Run this before proceeding to Module 02:

```bash
echo "=== Module 01 completion check ==="

echo ""
echo "--- Reference files ---"
ls -lh ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
ls -lh ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta.fai
ls -lh ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.dict

echo ""
echo "--- BWA-MEM2 index (expect 7 files) ---"
ls ${VARCALL}/00_reference/indexes/bwa/ | wc -l

echo ""
echo "--- STAR index (key file SA) ---"
ls -lh ${VARCALL}/00_reference/indexes/star/SA

echo ""
echo "--- Containers (expect 16) ---"
ls ${VARCALL}/containers/*.sif | wc -l

echo ""
echo "--- Raw reads ---"
echo "WGS    (expect 10): $(ls ${VARCALL}/01_raw_reads/wgs/*.fastq.gz | wc -l)"
echo "GBS    (expect  5): $(ls ${VARCALL}/01_raw_reads/gbs/*.fastq.gz | wc -l)"
echo "RNAseq (expect  5): $(ls ${VARCALL}/01_raw_reads/rnaseq/*.fastq.gz | wc -l)"

echo ""
echo "--- Trimmed reads (expect 20: 10 WGS + 5 GBS + 5 RNA-seq) ---"
ls ${VARCALL}/02_qc/fastp/*.fastq.gz | wc -l

echo ""
echo "--- MultiQC report ---"
ls -lh ${VARCALL}/02_qc/multiqc_report.html
```

---

## Discussion Questions

1. The samtools `.fai` index stores chromosome lengths. Using your `.fai`
   output, calculate: (a) the total nuclear genome size in Mb; (b) what fraction
   of the nuclear genome is on chromosome 1. How does this compare to the human
   genome, where chromosome 1 is ~8% of total genome size?

2. BWA-MEM2's `.bwt.2bit.64` index file is substantially larger than the
   original FASTA. What data structure does BWA-MEM2 use to enable rapid k-mer
   matching, and why does it require more storage space than the raw sequence?

3. STAR requires the GFF3 at index build time, but BWA-MEM2 indexes only the
   genome FASTA. Why does splice-aware RNA-seq alignment require knowledge of
   gene structures at index time, while DNA alignment does not?

4. The GBS fastp run uses `--disable_adapter_trimming` and `--trim_front1 5`,
   while the WGS run uses `--detect_adapter_for_pe` and no front-trimming. Why
   do these two tracks need different trimming strategies, even though the reads
   are the same length and come from the same original WGS data?

5. Vas-0 RNA-seq file is ~763 MB while the other four RNA-seq files are 2.6–3.2
   GB. Before proceeding to Module 04, how would you assess whether Vas-0 has
   sufficient depth for RNA-seq variant calling? What minimum read count would
   you expect to be necessary, and why?

---

*← [Module 00: Conceptual Framework](00_conceptual_framework.md) | [Module 02: WGS Variant Calling](02_wgs_variant_calling.md) →*