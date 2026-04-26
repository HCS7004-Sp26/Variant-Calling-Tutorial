# Module 03 — GBS Variant Calling

**HCS 7004 — Genome Analytics | The Ohio State University**

*← [Module 02: WGS Variant Calling](02_wgs_variant_calling.md) | [Module 04: RNA-seq Variant Calling](04_rnaseq_variant_calling.md) →*

---

> ### 🔑 Before you begin — username and environment variables
> Replace `<your_username>` with your OSC username in **every script** in this
> module. If you are starting here after a break or working through this module
> independently, also re-set these variables in your terminal before running any
> interactive commands:
>
> ```bash
> export user_name=<your_username>
> export VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
> export SHARED=/fs/scratch/PAS3260/Variant_Calling/Data
> export CONTAINERS=/fs/scratch/PAS3260/Variant_Calling/Containers
> ```
>
> SLURM batch scripts re-declare all paths internally — you still need to
> replace `<your_username>` on the `user_name=` line inside each script.
> Setting the variables above only affects your interactive terminal session.

---

## Learning Objectives

By the end of this module you will be able to:

1. Explain the biology of ApeKI-based GBS library preparation and how the
   restriction enzyme determines which portions of the genome are sampled
2. Demultiplex a pooled GBS library using Stacks `process_radtags`, correctly
   specifying enzyme, barcode file, and quality filters
3. Run the complete Stacks reference-based pipeline (`ref_map.pl` →
   `populations`) to produce a population-level VCF
4. Run the GATK4 variant calling pipeline on the same demultiplexed reads
   and compare the variant sets produced by both approaches
5. Explain the key differences between the Stacks genotype likelihood model
   and the GATK HaplotypeCaller model and when each is more appropriate
6. Evaluate GBS-specific quality metrics: locus coverage, missing data rate,
   and the relationship between restriction site proximity and read depth

---

## Biological Background

### Genotyping-by-Sequencing and the ApeKI protocol

Genotyping-by-Sequencing (GBS) is a reduced-representation sequencing strategy
that exploits restriction enzymes to reproducibly sample the same subset of the
genome across all individuals in a population. Rather than sequencing the entire
genome, GBS sequences only the regions immediately adjacent to restriction enzyme
cut sites — typically within 200–300 bp of each cut.

**ApeKI** recognises the degenerate sequence **GCWGC** (W = A or T), cutting
between the G and C to leave a 3-base 5' overhang:

```
5'...G  CWGC...3'       ApeKI cuts here ↓
3'...CGWC  G...5'

Resulting fragment ends:
  5'- CWG C -3'    (the "sticky end" that gets ligated to the sequencing adapter)
```

After ligation of barcoded adapters and size selection, sequencing begins from
the adapter and reads into the genomic sequence immediately following the cut
site. Every read therefore begins with the restriction site remnant (`CWG`) —
either `CAG` (W=A) or `CTG` (W=T) — immediately after the barcode is removed.
This is the signature that `process_radtags` uses to validate that a read
genuinely originated from an ApeKI cut site.

### Why ApeKI for *Arabidopsis*?

ApeKI's GCWGC recognition site occurs approximately every 1.5–2 kb in
the *Arabidopsis* genome, producing ~80,000–100,000 cut sites genome-wide.
After size selection for fragments in the 200–500 bp range, approximately
15–25% of the nuclear genome is represented. Crucially, ApeKI cut sites are
distributed throughout both gene-rich and intergenic regions in *Arabidopsis*,
unlike some enzymes that are strongly biased toward specific GC contexts.

### The multiplexed GBS library in this tutorial

The GBS data in this tutorial was **constructed in silico** from the WGS reads
by:
1. Identifying all ApeKI cut sites (GCAGC and GCTGC) in TAIR10
2. Extracting WGS reads overlapping those sites (±200 bp windows)
3. Prepending a unique 6 bp barcode and the CAGC overhang to each read

The result is a single multiplexed FASTQ file indistinguishable in structure
from a real GBS library output from an Illumina sequencer. Each read is 101 bp:

```
[6 bp barcode][4 bp overhang: CAGC][91 bp genomic sequence]
   AAGCTA        CAGC              ATGCTAGCTA...
```

| Sample | Barcode | Full prefix (removed by process_radtags) |
|---|---|---|
| Vas-0 | `AAGCTA` | `AAGCTACAGC` |
| Bez-9 | `CCGTAC` | `CCGTACCAGC` |
| Gen-8 | `TTGCAT` | `TTGCATCAGC` |
| Mah-6 | `GGACTT` | `GGACTTCAGC` |
| Usa-0 | `AACCGG` | `AACCGGCAGC` |

After `process_radtags` removes the 10 bp prefix, each read is 91 bp of
genuine genomic sequence originating from an ApeKI locus.

---

## Pipeline Overview

```
Shared data
  └── GBS/
        ├── GBS_multiplexed.fastq.gz   ─┐
        └── barcodes.tsv               ─┘→ Step 1: Copy to workspace

Step 2: FastQC on raw multiplexed library
         ↓
Step 3: process_radtags
         Demultiplex by barcode
         Validate ApeKI overhang (CAGC)
         Quality trim
         ↓
   Vas-0.fq.gz, Bez-9.fq.gz, Gen-8.fq.gz, Mah-6.fq.gz, Usa-0.fq.gz
         ↓
Step 4: FastQC on demultiplexed reads
         ↓
         ┌──────────────────────────────┐
         │   TRACK A — Stacks pipeline  │
         │   Step 5: BWA-MEM2 align     │
         │   Step 6: ref_map.pl +       │
         │           populations        │
         │   Output: stacks.vcf.gz      │
         └──────────────────────────────┘
                      ↕
         ┌──────────────────────────────┐
         │   TRACK B — GATK pipeline    │
         │   Step 7: BWA-MEM2 align     │
         │   Step 8: HaplotypeCaller    │
         │   Step 9: Joint genotyping   │
         │   Step 10: Hard filter       │
         │   Output: gatk_GBS.vcf.gz    │
         └──────────────────────────────┘
                      ↓
Step 11: Cross-pipeline comparison
         Which variants are shared? Which are unique?
```

---

## Step 1 — Copy GBS Data and Extend Directory Structure

The GBS data in the shared directory consists of just two files. Module 01
created a `04_gbs/` directory — we now add the subdirectories needed for this
module's outputs.

```bash
# Extend directory structure for Module 03
mkdir -p ${VARCALL}/04_gbs/{demux,aligned/{stacks,gatk},stacks,vcf}

# Copy multiplexed FASTQ and barcode file
cp ${SHARED}/GBS/GBS_multiplexed.fastq.gz ${VARCALL}/01_raw_reads/gbs/
cp ${SHARED}/GBS/barcodes.tsv             ${VARCALL}/01_raw_reads/gbs/

# Verify
ls -lh ${VARCALL}/01_raw_reads/gbs/
echo ""
echo "Barcode file contents:"
cat ${VARCALL}/01_raw_reads/gbs/barcodes.tsv
```

Expected barcode file format (tab-separated, no header):

```
AAGCTA  Vas-0
CCGTAC  Bez-9
TTGCAT  Gen-8
GGACTT  Mah-6
AACCGG  Usa-0
```

---

## Step 2 — FastQC on the Raw Multiplexed Library

Before any processing, assess the quality of the multiplexed library as a whole.
The multiplexed file is treated as a single sample at this stage.

```bash
cat > ${VARCALL}/scripts/03a_fastqc_gbs_raw.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=fastqc_gbs_raw
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastqc_gbs_raw_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastqc_gbs_raw_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

echo "=== FastQC: raw multiplexed GBS library: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/fastqc_0.12.1.sif \
  fastqc \
    --outdir ${VARCALL}/02_qc/fastqc_raw \
    --threads 4 \
    --extract \
    ${VARCALL}/01_raw_reads/gbs/GBS_multiplexed.fastq.gz

echo "=== FastQC complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/03a_fastqc_gbs_raw.sh
```

### What to expect in the raw multiplexed FastQC report

| Module | Expected pattern | Explanation |
|---|---|---|
| Per base sequence content | Strong non-random pattern in positions 1–10 | Positions 1–6 are one of 5 barcodes; positions 7–10 are always `CAGC` |
| Sequence duplication | Very high | Reads from the same ApeKI locus across 5 samples are near-identical |
| Per sequence GC content | Non-normal; may be multi-modal | A mix of reads from five different samples with slightly different GC backgrounds |
| Overrepresented sequences | Present | The `CAGC` overhang sequence is present in every read |

> **All of these patterns are expected** and are not indicators of library
> failure. They are inherent properties of a multiplexed GBS library.
> The per-sample FastQC reports after demultiplexing (Step 4) will look
> much closer to a normal single-end library.

---

## Step 3 — Demultiplexing with `process_radtags`

`process_radtags` is the Stacks tool that converts a pooled, barcoded GBS
library into one FASTQ file per sample. It simultaneously:

1. Identifies each read's barcode (allowing 1 bp mismatch by default)
2. Validates that the restriction site overhang is present immediately
   after the barcode (rejects reads without a genuine RAD site)
3. Performs quality trimming (sliding window, similar to `--cut_right` in fastp)
4. Removes the barcode and overhang, outputting clean genomic sequence

```bash
cat > ${VARCALL}/scripts/03b_process_radtags.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/process_radtags_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/process_radtags_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

MULTIPLEX=${VARCALL}/01_raw_reads/gbs/GBS_multiplexed.fastq.gz
BARCODES=${VARCALL}/01_raw_reads/gbs/barcodes.tsv
DEMUX_DIR=${VARCALL}/04_gbs/demux

echo "=== process_radtags demultiplexing: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/stacks_2.68.sif \
  process_radtags \
    -f ${MULTIPLEX} \
    -b ${BARCODES} \
    -o ${DEMUX_DIR} \
    -e apeKI \
    --quality \
    --rescue \
    --clean \
    --filter-illumina \
    --window-size 0.15 \
    --score-limit 10 \
    --len-limit 30 \
    --disable-rad-check

echo ""
echo "=== Demultiplexing complete: $(date) ==="
echo ""
echo "--- Output files ---"
ls -lh ${DEMUX_DIR}/

echo ""
echo "--- Read counts per sample ---"
for SAMPLE in Vas-0 Bez-9 Gen-8 Mah-6 Usa-0; do
    FQ=${DEMUX_DIR}/${SAMPLE}.fq.gz
    if [[ -f ${FQ} ]]; then
        N=$(apptainer exec \
              --bind ${VARCALL}:${VARCALL} \
              ${VARCALL}/containers/samtools_1.23.1.sif \
              bash -c "zcat ${FQ} | awk 'NR%4==1' | wc -l")
        echo "  ${SAMPLE}: ${N} reads"
    else
        echo "  ${SAMPLE}: FILE NOT FOUND"
    fi
done

echo ""
echo "--- Reads retained vs. discarded ---"
cat ${DEMUX_DIR}/process_radtags.*.log 2>/dev/null | \
  grep -E "Total|Retained|Discarded|Barcode" | head -20
EOF

sbatch ${VARCALL}/scripts/03b_process_radtags.sh
```

> **A note on `--disable-rad-check` for this tutorial.** In a real GBS
> experiment you should **not** use this flag — `process_radtags` validates
> that each read begins with the expected restriction site remnant, which is
> an important quality filter that catches contamination and chimeric reads.
> We disable it here because our simulated reads were constructed from
> internal WGS reads that do not all begin at exactly the cut site; a real
> ApeKI library would pass this check cleanly. The Module 03 discussion
> questions ask you to think about what removing this flag means for your
> false-positive rate.

### Expected output from `process_radtags`

After completion the `04_gbs/demux/` directory should contain:

```
Vas-0.fq.gz    ← clean single-end reads, ~91 bp, Vas-0 only
Bez-9.fq.gz
Gen-8.fq.gz
Mah-6.fq.gz
Usa-0.fq.gz
remainder.fq.gz                ← reads with no barcode match
process_radtags.XXXXX.log      ← detailed per-barcode statistics
```

Read retention per sample should be approximately 85–95% of the reads
attributed to that sample's barcode.

---

## Step 4 — FastQC on Demultiplexed Reads

```bash
cat > ${VARCALL}/scripts/03c_fastqc_demux.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=fastqc_demux
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastqc_demux_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/fastqc_demux_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=== FastQC demultiplexed: ${SAMPLE}: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/fastqc_0.12.1.sif \
  fastqc \
    --outdir ${VARCALL}/02_qc/fastqc_raw \
    --threads 4 \
    --extract \
    ${VARCALL}/04_gbs/demux/${SAMPLE}.fq.gz

echo "=== FastQC complete for ${SAMPLE}: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/03c_fastqc_demux.sh
```

### Comparing pre- and post-demultiplexing FastQC reports

After demultiplexing, per-sample FastQC reports should show:

- **Per base sequence content:** Positions 1–4 will show a strong fixed pattern
  (`CAG` or `CTG` from the ApeKI overhang remnant left after barcode removal)
  because `process_radtags` removes the barcode but retains the overhang. This
  is normal and expected for real GBS data.
- **Sequence duplication:** Remains elevated — reads from the same ApeKI locus
  within one sample are genuinely repetitive
- **Per base quality:** Should be uniformly high (the synthetic overhang was
  assigned Phred 40)

> **Comparison with Module 01 GBS FastQC.** If you ran FastQC on the
> per-sample GBS files in Module 01, compare those reports to the demultiplexed
> files here. The demultiplexed reads are 91 bp (10 bp shorter) and the
> barcode pattern in positions 1–6 has been replaced by the overhang pattern
> in positions 1–4.

---

## Track A — Stacks Reference-Based Pipeline

The Stacks reference-based pipeline consists of two main stages:

1. **`ref_map.pl`** — aligns reads to the reference and runs `gstacks` to
   call haplotypes at each RAD locus using a Stacks-specific population
   genetics model
2. **`populations`** — filters loci across samples and outputs a VCF

### Step 5 — Align Demultiplexed Reads for Stacks

`ref_map.pl` in Stacks 2.x expects pre-aligned, sorted BAM files as input.
We align each demultiplexed sample with BWA-MEM2 in single-end mode.

```bash
cat > ${VARCALL}/scripts/03d_bwa_align_gbs_stacks.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=bwa_gbs_stacks
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/bwa_gbs_stacks_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/bwa_gbs_stacks_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
THREADS=8

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

REF=${VARCALL}/00_reference/indexes/bwa/Athaliana_TAIR10.fasta
FQ=${VARCALL}/04_gbs/demux/${SAMPLE}.fq.gz
OUT=${VARCALL}/04_gbs/aligned/stacks/${SAMPLE}.bam

echo "=== BWA-MEM2 GBS alignment (Stacks track): ${SAMPLE}: $(date) ==="

# Single-end alignment — no R2 file, no -p flag
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bwamem2_2.2.3.sif \
  bwa-mem2 mem \
    -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}_GBS\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_GBS_lib1" \
    ${REF} \
    ${FQ} | \
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools sort \
    -@ ${THREADS} \
    -T ${VARCALL}/04_gbs/aligned/stacks/${SAMPLE}_tmp \
    -o ${OUT}

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools index ${OUT}

echo ""
echo "Alignment summary for ${SAMPLE}:"
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools flagstat ${OUT}

echo "=== ${SAMPLE} alignment complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/03d_bwa_align_gbs_stacks.sh
```

### Step 6 — ref_map.pl and populations

`ref_map.pl` takes the aligned BAMs and a **population map** (a tab-separated
file mapping each sample to a population) and runs `gstacks` internally to
call variant loci. We then run `populations` to filter and export the VCF.

**The population map** assigns samples to groups. For this tutorial all five
accessions are assigned to a single population because we are not testing for
population differentiation — we simply want a VCF comparable to the WGS output.

```bash
# ---- Create the population map ----
cat > ${VARCALL}/04_gbs/stacks/popmap.tsv << 'EOF'
Vas-0	1
Bez-9	1
Gen-8	1
Mah-6	1
Usa-0	1
EOF

echo "Population map:"
cat ${VARCALL}/04_gbs/stacks/popmap.tsv
```

```bash
cat > ${VARCALL}/scripts/03e_refmap_populations.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=refmap_populations
#SBATCH --account=PAS3260
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/refmap_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/refmap_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

BAM_DIR=${VARCALL}/04_gbs/aligned/stacks
STACKS_DIR=${VARCALL}/04_gbs/stacks
POPMAP=${STACKS_DIR}/popmap.tsv

echo "=== ref_map.pl: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/stacks_2.68.sif \
  ref_map.pl \
    --samples ${BAM_DIR} \
    --popmap ${POPMAP} \
    --out-path ${STACKS_DIR} \
    -T 8 \
    -X "gstacks: --max-clipped 0.5"

# Removed --rm-pcr-duplicates: this flag requires paired-end reads to
#   estimate insert sizes for duplicate detection. With single-end GBS
#   data, gstacks classifies every forward read as "unpaired" and
#   discards the entire dataset. PCR duplicates in single-end GBS are
#   instead handled at the locus level by Stacks's stacking model.
#
# Added --max-clipped 0.5: raises the soft-clip tolerance from the
#   default 20% (18 bp for 91 bp reads) to 50% (45 bp). Our simulated
#   reads come from throughout ApeKI windows rather than exactly from
#   cut sites, and divergent accessions accumulate soft-clipped bases
#   near indels relative to the Col-0 reference. The -- separator passes
#   this flag directly to gstacks rather than ref_map.pl.

echo ""
echo "=== ref_map.pl complete: $(date) ==="
echo "Stacks output files:"
ls -lh ${STACKS_DIR}/

# ---- Run populations to export VCF ----
echo ""
echo "=== populations: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/stacks_2.68.sif \
  populations \
    --in-path ${STACKS_DIR} \
    --out-path ${VARCALL}/04_gbs/vcf \
    --popmap ${POPMAP} \
    --vcf \
    --min-samples-overall 0.8 \
    --min-maf 0.05 \
    --threads 8

echo ""
echo "=== populations complete: $(date) ==="
echo "VCF output:"
ls -lh ${VARCALL}/04_gbs/vcf/

# ---- Index and sort the output VCF ----
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools sort \
    ${VARCALL}/04_gbs/vcf/populations.snps.vcf \
    -O z -o ${VARCALL}/04_gbs/vcf/stacks_GBS.vcf.gz

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${VARCALL}/04_gbs/vcf/stacks_GBS.vcf.gz

echo ""
echo "SNP count from Stacks:"
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${VARCALL}/04_gbs/vcf/stacks_GBS.vcf.gz | \
  grep "^SN" | grep "SNPs"
EOF

sbatch ${VARCALL}/scripts/03e_refmap_populations.sh
```

### What `ref_map.pl` / `gstacks` actually does

Unlike GATK's HaplotypeCaller, which uses local de novo assembly, Stacks
`gstacks` uses a **population genetics model** designed for RAD loci:

1. Groups reads by genomic position into **loci** (stacks of reads at each
   restriction site)
2. At each locus, models allele frequencies across the cohort using a
   multinomial likelihood
3. Calls genotypes under a diploid model, accounting for sequencing error
   and PCR duplicate rate
4. Applies a chi-square goodness-of-fit test for Hardy-Weinberg expectation
   as an additional quality filter

This model is optimised for GBS/RAD data where loci have distinctly different
depth profiles from WGS — very high depth immediately at the cut site, dropping
rapidly with distance.

---

## Track B — GATK Pipeline on Demultiplexed GBS Reads

Running GATK HaplotypeCaller on the same demultiplexed reads lets us directly
compare the two variant calling philosophies on identical input data.

### Step 7 — Align Demultiplexed Reads for GATK

We create a separate set of BAMs for the GATK track. Although the alignment
command is identical to Step 5, keeping them separate prevents the Stacks
`--rm-pcr-duplicates` and the GATK Picard duplicate-marking from interfering
with each other.

```bash
cat > ${VARCALL}/scripts/03f_bwa_align_gbs_gatk.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=bwa_gbs_gatk
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/bwa_gbs_gatk_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/bwa_gbs_gatk_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
THREADS=8

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

REF=${VARCALL}/00_reference/indexes/bwa/Athaliana_TAIR10.fasta
FQ=${VARCALL}/04_gbs/demux/${SAMPLE}.fq.gz
SORTED=${VARCALL}/04_gbs/aligned/gatk/${SAMPLE}_sorted.bam
MARKDUP=${VARCALL}/04_gbs/aligned/gatk/${SAMPLE}_markdup.bam
METRICS=${VARCALL}/04_gbs/aligned/gatk/${SAMPLE}_dup_metrics.txt

echo "=== BWA-MEM2 GBS alignment (GATK track): ${SAMPLE}: $(date) ==="

# Single-end alignment
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bwamem2_2.2.3.sif \
  bwa-mem2 mem \
    -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}_GBS\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}_GBS_lib1" \
    ${REF} \
    ${FQ} | \
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools sort \
    -@ ${THREADS} \
    -T ${VARCALL}/04_gbs/aligned/gatk/${SAMPLE}_tmp \
    -o ${SORTED}

# Mark duplicates
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/picard_3.4.0.sif \
  picard MarkDuplicates \
    -I ${SORTED} \
    -O ${MARKDUP} \
    -M ${METRICS} \
    --REMOVE_DUPLICATES false \
    --VALIDATION_STRINGENCY SILENT \
    --TMP_DIR ${VARCALL}/04_gbs/aligned/gatk/

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools index ${MARKDUP}

rm ${SORTED}

echo ""
echo "Duplicate metrics for ${SAMPLE}:"
grep -A 2 "ESTIMATED_LIBRARY_SIZE" ${METRICS} | tail -2

echo "=== ${SAMPLE} GATK GBS alignment complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/03f_bwa_align_gbs_gatk.sh
```

### Step 8 — HaplotypeCaller on GBS BAMs

```bash
cat > ${VARCALL}/scripts/03g_haplotypecaller_gbs.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hc_gbs
#SBATCH --account=PAS3260
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/hc_gbs_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/hc_gbs_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
BAM=${VARCALL}/04_gbs/aligned/gatk/${SAMPLE}_markdup.bam
GVCF=${VARCALL}/04_gbs/vcf/${SAMPLE}_GBS.g.vcf.gz

echo "=== HaplotypeCaller GBS: ${SAMPLE}: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk --java-options "-Xmx24g -Xms4g" HaplotypeCaller \
    --reference ${REF} \
    --input ${BAM} \
    --output ${GVCF} \
    --emit-ref-confidence GVCF \
    --sample-name ${SAMPLE} \
    --pcr-indel-model NONE \
    --native-pair-hmm-threads 4

echo "GVCF: $(ls -lh ${GVCF})"
echo "=== ${SAMPLE} HaplotypeCaller GBS complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/03g_haplotypecaller_gbs.sh
```

### Step 9 — Joint Genotyping of GBS GVCFs

```bash
cat > ${VARCALL}/scripts/03h_joint_genotype_gbs.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=joint_gbs
#SBATCH --account=PAS3260
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/joint_gbs_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/joint_gbs_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
DB=${VARCALL}/04_gbs/vcf/genomicsdb_gbs

echo "=== GenomicsDBImport + GenotypeGVCFs (GBS): $(date) ==="

SAMPLE_ARGS=""
for SAMPLE in Vas-0 Bez-9 Gen-8 Mah-6 Usa-0; do
    SAMPLE_ARGS="${SAMPLE_ARGS} -V ${VARCALL}/04_gbs/vcf/${SAMPLE}_GBS.g.vcf.gz"
done

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk --java-options "-Xmx24g -Xms4g" GenomicsDBImport \
    ${SAMPLE_ARGS} \
    --genomicsdb-workspace-path ${DB} \
    --intervals 1 --intervals 2 --intervals 3 \
    --intervals 4 --intervals 5 \
    --reader-threads 4 \
    --merge-input-intervals

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk --java-options "-Xmx24g -Xms4g" GenotypeGVCFs \
    --reference ${REF} \
    --variant gendb://${DB} \
    --output ${VARCALL}/04_gbs/vcf/gatk_GBS_raw.vcf.gz \
    --intervals 1 --intervals 2 --intervals 3 \
    --intervals 4 --intervals 5 \
    --merge-input-intervals

# Index the raw GBS cohort VCF — required by GATK SelectVariants in 03i
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk IndexFeatureFile \
    -I ${VARCALL}/04_gbs/vcf/gatk_GBS_raw.vcf.gz

echo ""
echo "Index created: ${VARCALL}/04_gbs/vcf/gatk_GBS_raw.vcf.gz.tbi"

echo ""
echo "Raw GATK GBS VCF:"
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${VARCALL}/04_gbs/vcf/gatk_GBS_raw.vcf.gz | grep "^SN"
EOF

sbatch ${VARCALL}/scripts/03h_joint_genotype_gbs.sh
```

### Step 10 — Hard Filter the GATK GBS VCF

```bash
cat > ${VARCALL}/scripts/03i_hard_filter_gbs.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=filter_gbs
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/filter_gbs_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/filter_gbs_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
VCF_DIR=${VARCALL}/04_gbs/vcf

echo "=== Hard filtering GATK GBS VCF: $(date) ==="

# Select SNPs
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk SelectVariants -R ${REF} \
    -V ${VCF_DIR}/gatk_GBS_raw.vcf.gz \
    --select-type-to-include SNP \
    -O ${VCF_DIR}/gatk_GBS_snps_raw.vcf.gz

# Apply SNP filters — identical thresholds to Module 02 for direct comparison
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk VariantFiltration -R ${REF} \
    -V ${VCF_DIR}/gatk_GBS_snps_raw.vcf.gz \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 60.0"             --filter-name "FS60" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${VCF_DIR}/gatk_GBS_snps_filtered.vcf.gz

# Extract PASS
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -f PASS \
    ${VCF_DIR}/gatk_GBS_snps_filtered.vcf.gz \
    -O z -o ${VCF_DIR}/gatk_GBS_PASS.vcf.gz

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${VCF_DIR}/gatk_GBS_PASS.vcf.gz

echo ""
echo "GATK GBS PASS SNP count:"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${VCF_DIR}/gatk_GBS_PASS.vcf.gz | grep "^SN" | grep "SNPs"

echo "=== Hard filtering complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/03i_hard_filter_gbs.sh
```

---

## Step 11 — Cross-Pipeline Comparison

This step is the pedagogical payoff of Module 03: comparing the variant sets
produced by Stacks and GATK on identical input reads, and overlaying both with
the WGS VCF from Module 02.

```bash
cat > ${VARCALL}/scripts/03j_compare_vcfs.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=compare_vcfs
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/compare_vcfs_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/compare_vcfs_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

STACKS_VCF=${VARCALL}/04_gbs/vcf/stacks_GBS.vcf.gz
GATK_GBS_VCF=${VARCALL}/04_gbs/vcf/gatk_GBS_PASS.vcf.gz
WGS_VCF=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
CMP_DIR=${VARCALL}/04_gbs/vcf/comparison
mkdir -p ${CMP_DIR}

echo "=== Cross-pipeline variant comparison: $(date) ==="

# ---- A. Stacks vs GATK GBS ----
echo ""
echo "--- Stacks vs. GATK GBS (same reads, different callers) ---"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools isec \
    ${STACKS_VCF} \
    ${GATK_GBS_VCF} \
    -p ${CMP_DIR}/stacks_vs_gatk_gbs \
    -O z

echo "Stacks-only SNPs (not in GATK GBS):"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/stacks_vs_gatk_gbs/0000.vcf.gz | grep "^SN" | grep "SNPs"

echo "GATK GBS-only SNPs (not in Stacks):"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/stacks_vs_gatk_gbs/0001.vcf.gz | grep "^SN" | grep "SNPs"

echo "SNPs in both Stacks and GATK GBS:"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/stacks_vs_gatk_gbs/0002.vcf.gz | grep "^SN" | grep "SNPs"

# ---- B. GBS (GATK) vs WGS (GATK) ----
echo ""
echo "--- GATK GBS vs. GATK WGS (same caller, different sequencing strategy) ---"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools isec \
    ${GATK_GBS_VCF} \
    ${WGS_VCF} \
    -p ${CMP_DIR}/gatk_gbs_vs_wgs \
    -O z

echo "GBS-only SNPs (not recovered by WGS):"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/gatk_gbs_vs_wgs/0000.vcf.gz | grep "^SN" | grep "SNPs"

echo "WGS-only SNPs (not recovered by GBS):"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/gatk_gbs_vs_wgs/0001.vcf.gz | grep "^SN" | grep "SNPs"

echo "SNPs in both GBS and WGS:"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/gatk_gbs_vs_wgs/0002.vcf.gz | grep "^SN" | grep "SNPs"

# ---- C. Summary table ----
echo ""
echo "=== Summary ==="
echo ""
python3 - << 'PYEOF'
import subprocess, gzip

def count_snps(vcf):
    try:
        result = subprocess.run(
            ['apptainer', 'exec',
             f'--bind', '/fs/scratch/PAS3260:/fs/scratch/PAS3260',
             '/fs/scratch/PAS3260/Variant_Calling/Containers/bcftools_1.23.1.sif',
             'bcftools', 'view', '-v', 'snps', '-H', vcf],
            capture_output=True, text=True)
        return sum(1 for line in result.stdout.split('\n') if line and not line.startswith('#'))
    except:
        return 'N/A'

import os
user = os.environ.get('user_name', '<user>')
base = f'/fs/scratch/PAS3260/{user}/variant_calling'
cmp  = f'{base}/04_gbs/vcf/comparison'

files = {
    'Stacks GBS':   f'{base}/04_gbs/vcf/stacks_GBS.vcf.gz',
    'GATK GBS':     f'{base}/04_gbs/vcf/gatk_GBS_PASS.vcf.gz',
    'GATK WGS':     f'{base}/03_wgs/filtered/cohort_PASS.vcf.gz',
    'Stacks∩GATK_GBS': f'{cmp}/stacks_vs_gatk_gbs/0002.vcf.gz',
    'GBS∩WGS':      f'{cmp}/gatk_gbs_vs_wgs/0002.vcf.gz',
}

print(f"{'Dataset':<22} {'SNP count':>12}")
print("-" * 36)
for label, vcf in files.items():
    if os.path.exists(vcf):
        n = count_snps(vcf)
        print(f"{label:<22} {str(n):>12}")
    else:
        print(f"{label:<22} {'FILE NOT FOUND':>12}")
PYEOF

echo ""
echo "=== Comparison complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/03j_compare_vcfs.sh
```

### Interpreting the comparison results

The comparison reveals three distinct perspectives on the same five genomes:

**Stacks vs. GATK GBS (same reads, different models):**
Variants unique to Stacks tend to be at loci with very high depth immediately
at the cut site — Stacks's RAD locus model is optimised for this depth profile.
Variants unique to GATK GBS tend to be at the edges of ApeKI windows where
depth is lower — GATK's local re-assembly handles low-depth sites better when
the full haplotype context is available.

**GATK GBS vs. GATK WGS (same caller, different strategies):**
This comparison directly answers the question: *what do you gain from WGS that
GBS cannot provide?* SNPs unique to WGS are predominantly in regions with no
ApeKI sites — centromeric regions, large intergenic stretches, and gene-dense
areas that happen to lack GCWGC motifs. SNPs shared between both are at ApeKI
loci and represent the highest-confidence common variants.

---

## Completion Checklist

```bash
echo "=== Module 03 completion check ==="

echo ""
echo "--- Demultiplexed FASTQ files (expect 5 + remainder) ---"
ls ${VARCALL}/04_gbs/demux/*.fq.gz | wc -l

echo ""
echo "--- Stacks output (key files) ---"
ls ${VARCALL}/04_gbs/stacks/gstacks.log 2>/dev/null && echo "gstacks: OK"
ls ${VARCALL}/04_gbs/vcf/stacks_GBS.vcf.gz 2>/dev/null && echo "Stacks VCF: OK"

echo ""
echo "--- GATK GBS output ---"
ls ${VARCALL}/04_gbs/vcf/gatk_GBS_PASS.vcf.gz 2>/dev/null && echo "GATK GBS VCF: OK"

echo ""
echo "--- Comparison output ---"
ls ${VARCALL}/04_gbs/vcf/comparison/stacks_vs_gatk_gbs/ 2>/dev/null
ls ${VARCALL}/04_gbs/vcf/comparison/gatk_gbs_vs_wgs/    2>/dev/null
```

---

## Discussion Questions

1. `process_radtags` was run with `--disable-rad-check` in this tutorial
   because our reads are simulated. In a real GBS experiment, removing this
   flag means that reads without a valid ApeKI overhang pass through. What
   types of reads would incorrectly enter your analysis, and how would this
   affect the resulting VCF?

2. The Stacks population map assigns all five accessions to a single population.
   How would the `populations` output change if you assigned European accessions
   (Vas-0, Bez-9, Gen-8, Mah-6) to population 1 and the Asian accession
   (Usa-0) to population 2? What new statistics would become available in the
   output, and how would the `--min-maf` filter behave differently?

3. Compare the duplicate rates reported by Picard MarkDuplicates for the GBS
   BAMs (Step 7) versus the WGS BAMs from Module 02. GBS duplicate rates are
   typically much higher. Why? Does this invalidate using MarkDuplicates
   for GBS data?

4. In the cross-pipeline comparison, some SNPs are found by Stacks but not
   by GATK on the same reads. If you examined the BAM alignments at one of
   these Stacks-unique sites using `samtools tview`, what alignment pattern
   might you expect to see that explains why GATK did not call the variant?

5. The GBS vs. WGS comparison shows that some variants are found only by GBS
   and not by WGS on the same accessions. Given that the GBS reads are a
   subset of the WGS reads, how is this possible? What property of the
   `populations` filtering in Stacks could explain the apparent presence of
   GBS-unique variants?

---

*← [Module 02: WGS Variant Calling](02_wgs_variant_calling.md) | [Module 04: RNA-seq Variant Calling](04_rnaseq_variant_calling.md) →*
