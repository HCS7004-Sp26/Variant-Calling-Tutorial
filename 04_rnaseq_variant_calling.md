# Module 04 — RNA-seq Variant Calling

**HCS 7004 — Genome Analytics | The Ohio State University**

*← [Module 03: GBS Variant Calling](03_gbs_variant_calling.md) | [Module 05: Variant Annotation & Filtering](05_variant_annotation.md) →*

---

> ### 🔑 Remember to replace `<your_username>`
> Every SLURM script contains `user_name=<your_username>`. Replace it with
> your actual OSC username (e.g. `jsmith`) before submitting any job.

---

## Learning Objectives

By the end of this module you will be able to:

1. Explain why RNA-seq variant calling requires a splice-aware aligner and
   why the N CIGAR operator must be resolved before variant calling
2. Run STAR in two-pass mode and interpret the key alignment statistics
   for single-end RNA-seq data
3. Apply the GATK RNA-seq best practices pipeline, including the
   `SplitNCigarReads` pre-processing step unique to spliced alignments
4. Explain which GATK hard filter thresholds differ between WGS and RNA-seq
   data and why
5. Compare RNA-seq variant calls to WGS calls at the same positions and
   interpret the sources of discordance
6. Use `ASEReadCounter` to quantify allelic expression at variant sites

---

## Biological Background

### Why RNA-seq variant calling is different from WGS

RNA-seq reads originate from processed mRNA, not genomic DNA. Three
properties make RNA-seq alignment and variant calling fundamentally
different from WGS:

**1. Splicing.** mRNA introns are removed during pre-mRNA processing.
A read that spans an exon-exon junction therefore aligns to two
discontinuous genomic regions, with a gap (encoded as `N` in the CIGAR
string) corresponding to the intron. Aligners like BWA-MEM2, designed for
DNA, cannot handle these gaps and will either misalign or discard such
reads entirely. STAR uses a splice-aware model that resolves junctions
correctly.

**2. Expression-level bias.** Highly expressed genes contribute many reads;
lowly expressed or unexpressed genes contribute none. Variant calling from
RNA-seq is therefore inherently biased toward variants in active genes.
A clean VCF from RNA-seq will never capture variants in intergenic regions
or silent chromatin — not because those variants don't exist, but because
no reads cover them.

**3. RNA editing and allele-specific expression.** RNA-to-DNA differences
can arise from adenosine-to-inosine (A-to-I) RNA editing, which is read
as A→G by the sequencer. These are genuine transcriptomic variants but not
genetic variants. Additionally, if a heterozygous locus is subject to
allele-specific expression, the read counts for the two alleles will be
unequal in RNA-seq data even though the underlying genotype is 1:1.
Both phenomena can produce variant calls that differ from the WGS truth.

### The GATK RNA-seq best practices pipeline

```
STAR (--twopassMode Basic)
        ↓
samtools sort + index
        ↓
Picard MarkDuplicates
        ↓
GATK SplitNCigarReads     ← key step unique to RNA-seq
        ↓
GATK HaplotypeCaller (GVCF mode, --dont-use-soft-clipped-bases)
        ↓
GenomicsDBImport → GenotypeGVCFs
        ↓
Hard filtering (FS threshold tightened for RNA-seq strand bias)
        ↓
Compare with WGS VCF → ASEReadCounter
```

### Why `SplitNCigarReads`?

After STAR alignment, reads spanning exon-exon junctions have CIGAR strings
like `50M1000N50M` — 50 bases aligned, 1000 base intronic gap (N), then 50
more aligned bases. GATK's variant caller cannot process the `N` operator
because it represents a genuine absence of sequence rather than a deletion,
and its pair-HMM model would interpret the gap as a large indel.

`SplitNCigarReads` resolves this by splitting each junction-spanning read
into two separate read segments at the `N` boundary:

```
Before:  ─────────[50M 1000N 50M]─────────
After:   ─────────[50M]          [50M]─────
                  ↑ split here
```

Each segment is then a contiguous alignment that HaplotypeCaller can
process normally.

---

## Pipeline Overview

```
Module 01 output
  └── 02_qc/fastp/
        ├── Vas-0_rnaseq.fastq.gz  ─┐
        ├── Bez-9_rnaseq.fastq.gz  ─┤
        ├── Gen-8_rnaseq.fastq.gz  ─┤→ Step 1: STAR 2-pass alignment
        ├── Mah-6_rnaseq.fastq.gz  ─┤
        └── Usa-0_rnaseq.fastq.gz  ─┘
                   ↓
          Step 2: Picard MarkDuplicates
                   ↓
          Step 3: GATK SplitNCigarReads
                   ↓
          Step 4: GATK HaplotypeCaller → per-sample GVCF
                   ↓
          Step 5: GenomicsDBImport → GenotypeGVCFs → raw VCF
                   ↓
          Step 6: Hard filtering → PASS VCF
                   ↓
          Step 7: Comparison with WGS VCF
                   ↓
          Step 8: ASEReadCounter → allelic read counts
```

---

## Step 0 — Extend Directory Structure

Module 01 created `05_rnaseq/{aligned,vcf}`. Add the subdirectories
needed for intermediate files in this module:

```bash
mkdir -p ${VARCALL}/05_rnaseq/aligned/{sorted,markdup,splitn}
mkdir -p ${VARCALL}/05_rnaseq/vcf/comparison
```

---

## Step 1 — STAR Two-Pass Alignment

STAR's `--twopassMode Basic` performs two alignment passes internally
per sample without requiring a separate index rebuild:

- **Pass 1:** Aligns reads using the pre-built splice junction database
  from Module 01 and discovers novel junctions not in the annotation
- **Pass 2:** Re-aligns reads using both the annotated and newly
  discovered junctions, improving alignment of reads at novel splice sites

```bash
cat > ${VARCALL}/scripts/04a_star_align.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --account=PAS3260
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/star_align_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/star_align_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
THREADS=12

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

STAR_IDX=${VARCALL}/00_reference/indexes/star
FQ=${VARCALL}/02_qc/fastp/${SAMPLE}_rnaseq.fastq.gz
OUT_PREFIX=${VARCALL}/05_rnaseq/aligned/sorted/${SAMPLE}_

echo "=== STAR alignment: ${SAMPLE}: $(date) ==="

[[ -f ${FQ} ]] || { echo "ERROR: FASTQ not found: ${FQ}"; exit 1; }

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/star_2.7.11b.sif \
  STAR \
    --runMode alignReads \
    --genomeDir ${STAR_IDX} \
    --readFilesIn ${FQ} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI NM MD AS XS \
    --outSAMattrRGline "ID:${SAMPLE}_RNA" "SM:${SAMPLE}" "PL:ILLUMINA" "LB:${SAMPLE}_RNA_lib1" \
    --twopassMode Basic \
    --outFileNamePrefix ${OUT_PREFIX} \
    --runThreadN ${THREADS} \
    --outBAMsortingThreadN ${THREADS} \
    --outSAMmapqUnique 60 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 10

# --outSAMtype BAM SortedByCoordinate: output sorted BAM directly
# --outSAMattributes: include standard tags needed by GATK
# --outSAMattrRGline: read group for GATK (must include SM: tag)
# --twopassMode Basic: two-pass alignment within a single run
# --outSAMmapqUnique 60: STAR assigns mapq=255 to uniquely mapped
#   reads; remap to 60 so Picard/GATK interpret quality correctly
# --outFilterIntronMotifs RemoveNoncanonical: remove reads with
#   non-canonical splice sites (GT-AG, GC-AG, AT-AC are canonical)
# --outFilterMismatchNmax 10: allow up to 10 mismatches; divergent
#   accessions need a relaxed threshold vs. the Col-0 reference

# Index the sorted BAM
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools index ${OUT_PREFIX}Aligned.sortedByCoord.out.bam

echo ""
echo "=== Alignment stats: ${SAMPLE} ==="
grep -E "Uniquely|splices|mismatch" \
  ${OUT_PREFIX}Log.final.out | head -15

echo ""
echo "Output BAM:"
ls -lh ${OUT_PREFIX}Aligned.sortedByCoord.out.bam
echo "=== ${SAMPLE} STAR complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/04a_star_align.sh
```

### Interpreting STAR alignment statistics

After completion, check `*Log.final.out` for each sample:

| Metric | Expected range | Flag if |
|---|---|---|
| Uniquely mapped reads | 75–90% | < 70% — check for contamination |
| Reads mapped to multiple loci | 5–15% | > 20% — reference quality issue |
| Unmapped reads (too short) | < 5% | > 10% — aggressive prior trimming |
| Splices: total | > 5M | Very low — index/annotation mismatch |
| Mismatch rate per base | < 1% | > 2% — library quality or wrong reference |

> **Why is the uniquely mapped rate lower for RNA-seq than WGS?**
> Multi-mapping reads are more common in RNA-seq because paralogous gene
> families often share highly similar exon sequences. A read from a
> member of a multigene family may align equally well to several loci.
> STAR reports these separately rather than forcing a mapping — which is
> the correct behaviour for variant calling where misassigned reads
> create false variants.

---

## Step 2 — Picard MarkDuplicates

Optical and PCR duplicates in RNA-seq are identified and flagged exactly
as in the WGS pipeline. Note that high duplicate rates are common and
expected for highly expressed transcripts — the same mRNA molecule is
genuinely sequenced many times. These are biological duplicates that
MarkDuplicates will also flag, slightly over-estimating the true PCR
duplicate rate.

```bash
cat > ${VARCALL}/scripts/04b_markduplicates_rna.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=markdup_rna
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/markdup_rna_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/markdup_rna_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

IN=${VARCALL}/05_rnaseq/aligned/sorted/${SAMPLE}_Aligned.sortedByCoord.out.bam
OUT=${VARCALL}/05_rnaseq/aligned/markdup/${SAMPLE}_markdup.bam
METRICS=${VARCALL}/05_rnaseq/aligned/markdup/${SAMPLE}_dup_metrics.txt

echo "=== Picard MarkDuplicates RNA: ${SAMPLE}: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/picard_3.4.0.sif \
  picard MarkDuplicates \
    -I ${IN} \
    -O ${OUT} \
    -M ${METRICS} \
    --REMOVE_DUPLICATES false \
    --VALIDATION_STRINGENCY SILENT \
    --TMP_DIR ${VARCALL}/05_rnaseq/aligned/markdup/

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools index ${OUT}

echo ""
echo "Duplicate metrics for ${SAMPLE}:"
grep -A 2 "ESTIMATED_LIBRARY_SIZE" ${METRICS} | tail -2
echo "=== ${SAMPLE} MarkDuplicates complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/04b_markduplicates_rna.sh
```

---

## Step 3 — SplitNCigarReads

This step is unique to RNA-seq variant calling. It splits reads that
span splice junctions into separate alignment segments, converting
`N`-containing CIGAR strings into clean continuous alignments that
HaplotypeCaller can process.

```bash
cat > ${VARCALL}/scripts/04c_splitncigar.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=splitncigar
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/splitncigar_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/splitncigar_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
IN=${VARCALL}/05_rnaseq/aligned/markdup/${SAMPLE}_markdup.bam
OUT=${VARCALL}/05_rnaseq/aligned/splitn/${SAMPLE}_splitn.bam

echo "=== SplitNCigarReads: ${SAMPLE}: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk --java-options "-Xmx24g -Xms4g" SplitNCigarReads \
    --reference ${REF} \
    --input ${IN} \
    --output ${OUT}

# No additional flags needed in GATK 4.x — SplitNCigarReads correctly
# handles the N→split conversion and re-assigns mapping quality scores.
# In GATK 3.x this required -rf ReassignOneMappingQuality; that flag
# is no longer needed because STAR's --outSAMmapqUnique 60 already
# handles the mapq reassignment at alignment time.

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools index ${OUT}

# Compare read counts before and after to verify the split worked
BEFORE=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools view -c -F 4 ${IN})

AFTER=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/samtools_1.23.1.sif \
  samtools view -c -F 4 ${OUT})

echo ""
echo "Read counts for ${SAMPLE}:"
echo "  Before SplitNCigarReads: ${BEFORE}"
echo "  After  SplitNCigarReads: ${AFTER}"
echo "  (After count is higher because junction reads become 2 segments)"
echo ""
echo "Output: $(ls -lh ${OUT})"
echo "=== ${SAMPLE} SplitNCigarReads complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/04c_splitncigar.sh
```

> **Why does read count increase after SplitNCigarReads?**
> Each junction-spanning read becomes two alignment segments, so the
> total number of records in the BAM increases. This is expected and
> correct. A read that spanned a 500 bp intron with 50 bp on each side
> becomes two 50 bp segments — each contributing independently to variant
> evidence at their respective positions.

---

## Step 4 — HaplotypeCaller on RNA-seq BAMs

HaplotypeCaller is run with two RNA-seq-specific flags that differ from
the WGS pipeline:

```bash
cat > ${VARCALL}/scripts/04d_haplotypecaller_rna.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=hc_rna
#SBATCH --account=PAS3260
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/hc_rna_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/hc_rna_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
BAM=${VARCALL}/05_rnaseq/aligned/splitn/${SAMPLE}_splitn.bam
GVCF=${VARCALL}/05_rnaseq/vcf/${SAMPLE}_rna.g.vcf.gz

echo "=== HaplotypeCaller RNA: ${SAMPLE}: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk --java-options "-Xmx24g -Xms4g" HaplotypeCaller \
    --reference ${REF} \
    --input ${BAM} \
    --output ${GVCF} \
    --emit-ref-confidence GVCF \
    --sample-name ${SAMPLE} \
    --dont-use-soft-clipped-bases \
    --pcr-indel-model NONE \
    --native-pair-hmm-threads 4

# --dont-use-soft-clipped-bases: RNA-seq reads accumulate soft-clipped
#   bases at splice junctions and read ends. These bases are often of
#   lower quality and can introduce false mismatches. This flag excludes
#   them from the variant calling model entirely.
#   Note: this flag is NOT used in the WGS pipeline — for WGS, soft-
#   clipped bases at indels contain real information.

echo ""
echo "GVCF: $(ls -lh ${GVCF})"
echo "=== ${SAMPLE} HaplotypeCaller RNA complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/04d_haplotypecaller_rna.sh
```

### WGS vs. RNA-seq HaplotypeCaller differences

| Flag | WGS | RNA-seq | Reason |
|---|---|---|---|
| `--dont-use-soft-clipped-bases` | Not used | **Required** | Splice junction soft-clips create false mismatches |
| `--pcr-indel-model` | `NONE` | `NONE` | Same — high-quality libraries |
| `--emit-ref-confidence` | `GVCF` | `GVCF` | Same — required for joint genotyping |
| `--native-pair-hmm-threads` | 4 | 4 | Same |

---

## Step 5 — Joint Genotyping

GenomicsDBImport and GenotypeGVCFs are run identically to the WGS
pipeline. The VCF index is created before the stats summary.

```bash
cat > ${VARCALL}/scripts/04e_joint_genotype_rna.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=joint_rna
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/joint_rna_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/joint_rna_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
DB=${VARCALL}/05_rnaseq/vcf/genomicsdb_rna
RAW_VCF=${VARCALL}/05_rnaseq/vcf/rna_cohort_raw.vcf.gz

echo "=== GenomicsDBImport RNA: $(date) ==="

SAMPLE_ARGS=""
for SAMPLE in Vas-0 Bez-9 Gen-8 Mah-6 Usa-0; do
    GVCF=${VARCALL}/05_rnaseq/vcf/${SAMPLE}_rna.g.vcf.gz
    [[ -f ${GVCF} ]] || { echo "ERROR: GVCF not found: ${GVCF}"; exit 1; }
    SAMPLE_ARGS="${SAMPLE_ARGS} -V ${GVCF}"
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

echo ""
echo "=== GenotypeGVCFs RNA: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk --java-options "-Xmx24g -Xms4g" GenotypeGVCFs \
    --reference ${REF} \
    --variant gendb://${DB} \
    --output ${RAW_VCF} \
    --intervals 1 --intervals 2 --intervals 3 \
    --intervals 4 --intervals 5 \
    --merge-input-intervals \
    --include-non-variant-sites false

# Index before any downstream tool reads this file
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk IndexFeatureFile \
    -I ${RAW_VCF}

echo ""
echo "Raw RNA VCF: $(ls -lh ${RAW_VCF})"
echo "Index: $(ls -lh ${RAW_VCF}.tbi)"

echo ""
echo "=== Basic site counts ==="
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${RAW_VCF} | grep "^SN"

echo "=== Joint genotyping RNA complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/04e_joint_genotype_rna.sh
```

---

## Step 6 — Hard Filtering

RNA-seq variant calling uses the same framework as the WGS pipeline
but with two threshold adjustments reflecting the different error profile
of spliced reads:

| Filter | WGS threshold | RNA-seq threshold | Reason for difference |
|---|---|---|---|
| `FS` SNPs | > 60.0 | > **30.0** | RNA-seq has stronger inherent strand bias because expression is often strand-specific; tighter FS removes more artefacts |
| `FS` indels | > 200.0 | > **30.0** | Same reasoning; indels at splice sites are especially prone to strand bias |
| `MQRankSum` | < -12.5 | **removed** | Mapping quality distributions differ between exonic and intronic reads; this filter is too aggressive for RNA-seq |
| All others | same | same | QD, SOR, MQ, ReadPosRankSum thresholds are appropriate for both |

```bash
cat > ${VARCALL}/scripts/04f_hard_filter_rna.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=filter_rna
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/filter_rna_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/filter_rna_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
VCF_DIR=${VARCALL}/05_rnaseq/vcf

echo "=== Hard filtering RNA VCF: $(date) ==="

# ---- Select SNPs ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk SelectVariants -R ${REF} \
    -V ${VCF_DIR}/rna_cohort_raw.vcf.gz \
    --select-type-to-include SNP \
    -O ${VCF_DIR}/rna_snps_raw.vcf.gz

# ---- Filter SNPs (RNA-seq thresholds) ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk VariantFiltration -R ${REF} \
    -V ${VCF_DIR}/rna_snps_raw.vcf.gz \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 30.0"             --filter-name "FS30" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${VCF_DIR}/rna_snps_filtered.vcf.gz

# ---- Select indels ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk SelectVariants -R ${REF} \
    -V ${VCF_DIR}/rna_cohort_raw.vcf.gz \
    --select-type-to-include INDEL \
    -O ${VCF_DIR}/rna_indels_raw.vcf.gz

# ---- Filter indels (RNA-seq thresholds) ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk VariantFiltration -R ${REF} \
    -V ${VCF_DIR}/rna_indels_raw.vcf.gz \
    --filter-expression "QD < 2.0"               --filter-name "QD2" \
    --filter-expression "FS > 30.0"              --filter-name "FS30" \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${VCF_DIR}/rna_indels_filtered.vcf.gz

# ---- Merge and extract PASS ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk MergeVcfs \
    -I ${VCF_DIR}/rna_snps_filtered.vcf.gz \
    -I ${VCF_DIR}/rna_indels_filtered.vcf.gz \
    -O ${VCF_DIR}/rna_cohort_filtered.vcf.gz

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -f PASS \
    ${VCF_DIR}/rna_cohort_filtered.vcf.gz \
    -O z -o ${VCF_DIR}/rna_cohort_PASS.vcf.gz

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${VCF_DIR}/rna_cohort_PASS.vcf.gz

echo ""
echo "=== Filtering summary ==="
echo "Raw variants:"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${VCF_DIR}/rna_cohort_raw.vcf.gz | \
  grep "^SN" | grep -E "SNPs|indels"

echo ""
echo "PASS variants:"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${VCF_DIR}/rna_cohort_PASS.vcf.gz | \
  grep "^SN" | grep -E "SNPs|indels"

echo "=== Hard filtering RNA complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/04f_hard_filter_rna.sh
```

### What to expect from the RNA VCF

For five *Arabidopsis* accessions with ~7–35M RNA-seq reads each:

| Metric | Expected range | Notes |
|---|---|---|
| Raw SNPs before filtering | 50,000–200,000 | Much fewer than WGS — only expressed regions covered |
| PASS SNPs after filtering | 20,000–80,000 | Ts/Tv should still be ~2.0–2.2 |
| Ts/Tv ratio | 2.0–2.3 | Lower than WGS would suggest a problem |
| Missing genotypes | Higher than WGS | Expected — unexpressed genes have no coverage |

> **Why so many fewer variants than WGS?** The RNA-seq reads cover only
> the ~5–10% of the genome that is actively transcribed in rosette leaf
> tissue under the conditions used. All intergenic SNPs, intronic SNPs
> outside splice sites, and variants in silent genes are invisible.

---

## Step 7 — Comparison with WGS Variants

The most informative analysis in this module is comparing the RNA-seq
variant calls to the WGS ground truth. Sites present in both represent
variants confidently called by two independent methods — the highest
confidence class. Sites unique to RNA-seq deserve scrutiny.

```bash
cat > ${VARCALL}/scripts/04g_rna_wgs_comparison.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=rna_wgs_compare
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/rna_wgs_compare_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/rna_wgs_compare_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

RNA_VCF=${VARCALL}/05_rnaseq/vcf/rna_cohort_PASS.vcf.gz
WGS_VCF=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
CMP_DIR=${VARCALL}/05_rnaseq/vcf/comparison

echo "=== RNA-seq vs WGS variant comparison: $(date) ==="

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools isec \
    ${RNA_VCF} \
    ${WGS_VCF} \
    -p ${CMP_DIR} \
    -O z

# isec output:
#   0000.vcf.gz = RNA-seq only (not in WGS)
#   0001.vcf.gz = WGS only (not in RNA-seq)
#   0002.vcf.gz = in both RNA-seq and WGS (concordant)

echo ""
echo "=== Comparison results ==="
echo ""
echo "RNA-seq only (not confirmed by WGS):"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/0000.vcf.gz | grep "^SN" | grep "SNPs"

echo ""
echo "WGS only (in genomic DNA but not detected in RNA):"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/0001.vcf.gz | grep "^SN" | grep "SNPs"

echo ""
echo "Concordant (confirmed by both methods):"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${CMP_DIR}/0002.vcf.gz | grep "^SN" | grep "SNPs"

# ---- Ts/Tv of each class ----
echo ""
echo "=== Ts/Tv by concordance class ==="
for IDX in 0000 0001 0002; do
    LABEL=""
    if [[ ${IDX} == "0000" ]]; then LABEL="RNA-only"; fi
    if [[ ${IDX} == "0001" ]]; then LABEL="WGS-only "; fi
    if [[ ${IDX} == "0002" ]]; then LABEL="Concordant"; fi

    TSTV=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
      ${VARCALL}/containers/bcftools_1.23.1.sif \
      bcftools stats ${CMP_DIR}/${IDX}.vcf.gz | \
      grep "^TSTV" | awk '{print $5}')

    echo "  ${LABEL}: Ts/Tv = ${TSTV}"
done

echo ""
echo "=== Comparison complete: $(date) ==="
echo ""
echo "Interpretation guide:"
echo "  Concordant SNPs:     highest confidence — validated by two methods"
echo "  WGS-only SNPs:       in unexpressed genes/regions — expected"
echo "  RNA-only SNPs:       investigate — potential RNA editing or"
echo "                       alignment artefacts at splice junctions"
EOF

sbatch ${VARCALL}/scripts/04g_rna_wgs_comparison.sh
```

> **The Ts/Tv test for RNA-only variants is powerful.** If RNA-only
> variants have a Ts/Tv ratio of ~2.0, they are likely real variants in
> lowly expressed genes that WGS missed by chance. If the ratio is well
> below 1.8, the RNA-only set is enriched for artefacts — most likely
> from misalignment at splice junctions or RNA editing events (A→I
> editing appears as A→G, a transition, but it would inflate the
> transition count beyond the expected ratio).

---

## Step 8 — Allele-Specific Expression with ASEReadCounter

`ASEReadCounter` counts reference and alternate allele reads at each
variant site in each RNA-seq BAM. This quantifies whether a variant
site shows balanced or biased allelic expression.

For *Arabidopsis*, which is nearly homozygous, the most informative
application is identifying sites that are **homozygous in the genome
but show unexpected read count imbalance in the transcriptome** —
a potential signature of RNA editing or post-transcriptional regulation.

```bash
cat > ${VARCALL}/scripts/04h_ase_readcounter.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=ase_counter
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/ase_%A_%a.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/ase_%A_%a.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

SAMPLES=(Vas-0 Bez-9 Gen-8 Mah-6 Usa-0)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
BAM=${VARCALL}/05_rnaseq/aligned/splitn/${SAMPLE}_splitn.bam

# Use concordant SNPs (confirmed by both WGS and RNA-seq) as the site list
SITES=${VARCALL}/05_rnaseq/vcf/comparison/0002.vcf.gz
OUT=${VARCALL}/05_rnaseq/vcf/${SAMPLE}_ASE_counts.table

echo "=== ASEReadCounter: ${SAMPLE}: $(date) ==="

[[ -f ${SITES} ]] || { echo "ERROR: Concordant VCF not found. Run 04g first."; exit 1; }

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/gatk4_4.6.2.0.sif \
  gatk --java-options "-Xmx24g -Xms4g" ASEReadCounter \
    --reference ${REF} \
    --input ${BAM} \
    --variant ${SITES} \
    --output ${OUT} \
    --min-depth-of-non-filtered-base 10 \
    --min-mapping-quality 10 \
    --min-base-quality 20

# --min-depth-of-non-filtered-base 10: require at least 10 reads
#   at a site to estimate allele ratios reliably
# --min-mapping-quality 10: exclude poorly mapped reads
# --min-base-quality 20: exclude low-quality base calls

echo ""
echo "ASE counts output for ${SAMPLE}:"
head -5 ${OUT}
echo "Total sites counted:"
wc -l < ${OUT}

echo "=== ${SAMPLE} ASEReadCounter complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/04h_ase_readcounter.sh
```

### Interpreting the ASEReadCounter output

The output is a tab-delimited table with one row per site per sample:

```
contig  position  variantID  refAllele  altAllele  refCount  altCount  totalCount  lowMAPQDepth  lowBaseQDepth  rawDepth  otherBases  improperPairs
1       142387    .          G          A          45        2         47          3             1              51        0           0
```

For a homozygous accession you expect either:
- `refCount` ≈ `totalCount` and `altCount` ≈ 0 (homozygous REF in this accession)
- `altCount` ≈ `totalCount` and `refCount` ≈ 0 (homozygous ALT in this accession)

Sites where both `refCount` and `altCount` are substantial in a sample
expected to be homozygous are candidates for RNA editing or mapping
artefacts and warrant further investigation.

---

## Completion Checklist

```bash
echo "=== Module 04 completion check ==="

echo ""
echo "--- STAR BAMs (expect 5) ---"
ls ${VARCALL}/05_rnaseq/aligned/sorted/*Aligned.sortedByCoord.out.bam | wc -l

echo ""
echo "--- SplitN BAMs (expect 5) ---"
ls ${VARCALL}/05_rnaseq/aligned/splitn/*_splitn.bam | wc -l

echo ""
echo "--- Per-sample GVCFs (expect 5) ---"
ls ${VARCALL}/05_rnaseq/vcf/*_rna.g.vcf.gz | wc -l

echo ""
echo "--- Raw cohort VCF + index ---"
ls -lh ${VARCALL}/05_rnaseq/vcf/rna_cohort_raw.vcf.gz
ls -lh ${VARCALL}/05_rnaseq/vcf/rna_cohort_raw.vcf.gz.tbi

echo ""
echo "--- PASS VCF ---"
ls -lh ${VARCALL}/05_rnaseq/vcf/rna_cohort_PASS.vcf.gz

echo ""
echo "--- SNP count in PASS VCF ---"
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${VARCALL}/05_rnaseq/vcf/rna_cohort_PASS.vcf.gz | \
  grep "^SN" | grep "SNPs"

echo ""
echo "--- Comparison output (expect 3 VCFs) ---"
ls ${VARCALL}/05_rnaseq/vcf/comparison/*.vcf.gz | wc -l

echo ""
echo "--- ASE count tables (expect 5) ---"
ls ${VARCALL}/05_rnaseq/vcf/*_ASE_counts.table | wc -l
```

---

## Discussion Questions

1. STAR's `--outSAMmapqUnique 60` flag remaps the default mapping quality
   of uniquely aligned reads from 255 to 60. Why does STAR use 255 as its
   default mapping quality, and why must this be changed before GATK
   processes the BAM? What would happen to HaplotypeCaller output if
   you forgot this flag?

2. After `SplitNCigarReads`, the total read count in the BAM is higher
   than before. A student claims this means the same bases are now being
   counted twice in the variant calling model, inflating confidence at
   junction-spanning variants. Is this concern valid? How does
   HaplotypeCaller handle split reads at the same locus?

3. The hard filter for Fisher Strand (`FS`) is set to 30 for RNA-seq
   versus 60 for WGS. Look at 5 random entries in your raw RNA VCF
   using `bcftools view` and examine the `FS` values. What is the
   distribution of `FS` in RNA-seq data compared to what you saw in
   Module 02? Does the tighter threshold seem justified?

4. In the WGS vs. RNA-seq comparison, a subset of SNPs is detected only
   by RNA-seq and not by WGS. Propose an experimental design using
   existing data from this tutorial to distinguish whether these
   RNA-only variants are: (a) genuine variants in lowly covered WGS
   regions, (b) RNA editing events, or (c) alignment artefacts at
   novel splice junctions.

5. The ASEReadCounter identifies sites where a homozygous accession shows
   mixed allele counts in RNA-seq. For Usa-0, which carries the most
   divergent haplotype from the Col-0 reference in this panel, would
   you expect more or fewer apparent ASE sites compared to Vas-0 (the
   most similar to Col-0)? Why?

---

*← [Module 03: GBS Variant Calling](03_gbs_variant_calling.md) | [Module 05: Variant Annotation & Filtering](05_variant_annotation.md) →*
