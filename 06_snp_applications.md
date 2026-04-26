# Module 06 — From SNPs to Genotyping Technologies

**HCS 7004 — Genome Analytics | The Ohio State University**

*← [Module 05: Variant Annotation & Filtering](05_variant_annotation.md)*

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
> R scripts in this tutorial use the same username as a plain string:
> ```r
> user_name <- "<your_username>"
> ```
>
> SLURM batch scripts re-declare all paths internally — you still need to
> replace `<your_username>` on the `user_name=` line inside each script.
> Setting the variables above only affects your interactive terminal session.

---

## Learning Objectives

By the end of this module you will be able to:

1. Explain the biological and technical criteria used to select SNPs for
   inclusion on a genotyping array or amplicon panel
2. Trace the complete pipeline from population-level VCF to a candidate
   marker set, applying MAF, LD, uniqueness, and annotation filters
3. Describe how the Illumina Infinium and Affymetrix Axiom array technologies
   work at the probe level, and why certain genomic contexts fail probe design
4. Design a minimal KASP assay panel from your own VCF data, applying the
   sequence context criteria that determine assay success rate
5. Explain the principles of Amplicon Sequencing (AmpSeq) and how it differs
   from GBS and SNP arrays in cost, flexibility, and throughput
6. Connect variant discovery data to the practical workflow of genomic
   selection and marker-assisted breeding in a diploid or polyploid crop

---

## Biological Background

### From VCF to genotyping tool — the fundamental question

Every SNP chip, KASP panel, and AmpSeq array in use in plant breeding today
began exactly where you are now: with a population-level VCF file from
whole-genome resequencing or GBS. The transition from VCF to genotyping
product is not a single step but a cascade of design decisions, each of which
shapes what the final tool can and cannot detect.

The central challenge is **selection**: a VCF for five *Arabidopsis*
accessions may contain 1–2 million SNPs, but a practical genotyping panel
might target 1,000–100,000. Every SNP excluded from the panel is a locus
where a QTL, introgressions, or causal variant might be invisible. Every
SNP included must be technically amenable to the chosen platform. The criteria
for selection differ between platforms but share a common logic:

```
Population VCF (1–10 million SNPs)
       ↓
Quality filters         Remove low-confidence calls (already done in Module 05)
       ↓
Population filters      MAF ≥ threshold, missingness ≤ threshold
       ↓
LD filters              Remove redundant markers (r² pruning)
       ↓
Genome context filters  Remove repetitive regions, low-complexity sequence
       ↓
Annotation filters      Prioritise by functional consequence if desired
       ↓
Platform-specific       Probe uniqueness (arrays), amplicon size (AmpSeq),
design filters          GC content (KASP), cut-site proximity (GBS)
       ↓
Candidate marker set    Ready for probe/primer design or chip submission
```

### Genotyping platforms — a comparative overview

| Platform | Principle | Throughput | Cost per sample | Flexibility | Best use case |
|---|---|---|---|---|---|
| SNP array (Illumina/Axiom) | Hybridisation to fixed probes | High (96–384/run) | Low at scale | None — fixed markers | Large breeding programs, GWAS |
| GBS / RADseq | Restriction enzyme + sequencing | High | Low | Moderate (enzyme choice) | Diversity panels, new crops |
| AmpSeq | PCR amplification + sequencing | High | Low–moderate | High (custom panel) | Targeted QTL regions, validation |
| KASP | Allele-specific PCR + FRET | Moderate | Moderate | High (any SNP) | MAS, genotype confirmation |
| Whole-genome sequencing | Direct sequencing | Low | High | Complete | Reference discovery, structural variants |

This module focuses on the design logic for **SNP arrays**, **KASP**, and
**AmpSeq** — the three platforms most commonly encountered in plant breeding
programs — using our *Arabidopsis* VCF data as the working example.

---

## Pipeline Overview

```
Module 05 outputs:
  ├── 03_wgs/filtered/cohort_PASS.vcf.gz           (full WGS SNP set)
  ├── 06_annotation/annotated/cohort_PASS_annotated.vcf.gz
  └── 07_population/filtered/wgs_filtered_final.vcf.gz

Step 1: Understand the 1001 Genomes → Arabidopsis array design
        (conceptual + data exploration)
        ↓
Step 2: Apply array design filters to our VCF
        (MAF, LD, annotation priority)
        ↓
Step 3: Genomic context filter — remove repetitive and low-complexity regions
        ↓
Step 4: Extract sequence context for probe/primer design
        ↓
Step 5: KASP assay design — sequence context criteria
        ↓
Step 6: AmpSeq panel design — amplicon size and primer placement
        ↓
Step 7: Visualise the final marker set
        ↓
Step 8: Connecting to genomic selection — what marker density is enough?
```

---

## Step 0 — Extend Directory Structure

```bash
mkdir -p ${VARCALL}/08_chip_design/{candidates,context,kasp,ampseq,plots}
```

---

## Step 1 — The 1001 Genomes → *Arabidopsis* Array: A Case Study

The **Arabidopsis 1001 Genomes SNP Array** (Affymetrix Axiom platform,
~1.1 million markers) was designed directly from the same VCF we have been
working with throughout this tutorial. Understanding its design criteria
makes the abstract design pipeline concrete.

### 1.1 Explore the population-level variant landscape

Before selecting markers, characterise the variant set you are working with.
This is also a sanity check that the upstream pipeline produced sensible results.

```bash
cat > ${VARCALL}/scripts/06a_characterise_variants.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=char_variants
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/char_variants_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/char_variants_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

VCF=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
ANN_VCF=${VARCALL}/06_annotation/annotated/cohort_PASS_annotated.vcf.gz
OUT=${VARCALL}/08_chip_design/candidates

echo "=== Variant characterisation: $(date) ==="

# ---- SNP counts per chromosome ----
echo ""
echo "--- SNP counts per chromosome ---"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -v snps ${VCF} 2>/dev/null | \
  grep -v "^#" | cut -f1 | sort | uniq -c | sort -k2,2V

# ---- MAF distribution ----
echo ""
echo "--- Computing allele frequencies ---"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/vcftools_0.1.17.sif \
  vcftools \
    --gzvcf ${VCF} \
    --remove-indels \
    --freq2 \
    --out ${OUT}/all_snps \
    2>/dev/null

echo "MAF distribution (frequency of minor allele):"
awk 'NR>1 {
    maf = ($5 < $6) ? $5 : $6
    if (maf < 0.05)        c1++
    else if (maf < 0.10)   c2++
    else if (maf < 0.20)   c3++
    else if (maf < 0.30)   c4++
    else                   c5++
    total++
}
END {
    print "  MAF < 0.05 (rare):     " c1 " (" int(c1/total*100) "%)"
    print "  MAF 0.05–0.10:         " c2 " (" int(c2/total*100) "%)"
    print "  MAF 0.10–0.20:         " c3 " (" int(c3/total*100) "%)"
    print "  MAF 0.20–0.30:         " c4 " (" int(c4/total*100) "%)"
    print "  MAF ≥ 0.30 (common):   " c5 " (" int(c5/total*100) "%)"
    print "  Total SNPs:            " total
}' ${OUT}/all_snps.frq

# ---- Count variants by impact category ----
echo ""
echo "--- Variants by SnpEff impact ---"
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools query -f '%INFO/ANN\n' ${ANN_VCF} 2>/dev/null | \
  awk -F'|' '{print $3}' | sort | uniq -c | sort -rn | head -10

echo "=== Characterisation complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06a_characterise_variants.sh
```

### 1.2 The 1001 Genomes array design criteria (conceptual)

The Affymetrix Axiom *Arabidopsis* array used the following criteria,
applied sequentially to the 10.7 million SNPs in the 1,135-accession VCF:

| Filter | Threshold | Rationale |
|---|---|---|
| Minor allele frequency | MAF ≥ 0.05 | Rare variants provide little information in GWAS panels |
| Missing data | ≤ 20% | Sites genotyped in fewer than 80% of accessions are unreliable |
| Probe uniqueness | BLAST e-value ≥ 10⁻⁵ at off-target sites | Probes that cross-hybridise produce unreliable calls |
| Sequence complexity | No low-complexity sequence within 35 bp | Poly-A/T runs prevent specific hybridisation |
| LD-based selection | Approximately uniform genome coverage | Ensure LD blocks are tagged without over-redundancy |
| Functional enrichment | Coding SNPs over-represented by ~2× | Missense and synonymous variants are more likely to be functional |

In our five-accession dataset the same logic applies, scaled down:

- With only five samples, MAF ≥ 0.20 is more meaningful (one minor allele
  in five diploid individuals at AF = 0.10 means only one copy of the
  alternate allele was observed)
- LD pruning criteria remain the same (r² < 0.2 in 50-SNP windows)
- Probe context criteria are platform-independent

---

## Step 2 — Apply Array Design Filters

Starting from the population-filtered VCF from Module 05, apply the
additional MAF and annotation-priority filters appropriate for array design.

```bash
cat > ${VARCALL}/scripts/06b_array_design_filter.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=array_filter
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/array_filter_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/array_filter_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

# Start from the population-filtered VCF (MAF≥0.05 already applied)
IN=${VARCALL}/07_population/filtered/wgs_filtered_final.vcf.gz
ANN_VCF=${VARCALL}/06_annotation/annotated/cohort_PASS_annotated.vcf.gz
OUT=${VARCALL}/08_chip_design/candidates

echo "=== Array design filtering: $(date) ==="

N_IN=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${IN} 2>/dev/null | wc -l)
echo "Input SNPs: ${N_IN}"

# ---- Step A: Raise MAF threshold to 0.20 for 5-sample panel ----
# With 5 diploid samples (10 alleles), MAF=0.10 means only 1 minor
# allele observed — too unreliable for array probe design.
# MAF=0.20 requires at least 2 minor allele copies.
echo ""
echo "[A] Applying MAF ≥ 0.20 threshold..."
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view \
    --min-af 0.20:minor \
    ${IN} \
    -O z -o ${OUT}/candidates_maf20.vcf.gz 2>/dev/null

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${OUT}/candidates_maf20.vcf.gz

N_MAF=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${OUT}/candidates_maf20.vcf.gz 2>/dev/null | wc -l)
echo "  After MAF ≥ 0.20: ${N_MAF} SNPs"

# ---- Step B: Apply LD pruning from Module 05 prune.in list ----
# Extract only SNPs in the LD-pruned set. The prune.in file contains
# variant IDs in CHR:POS format matching our --set-missing-var-ids @:#
echo ""
echo "[B] Applying LD pruning (r² < 0.2 from Module 05)..."
PRUNE_IN=${VARCALL}/07_population/plink/arabidopsis_wgs_pruned.prune.in

# Convert prune.in IDs to a regions BED for bcftools
awk 'BEGIN{OFS="\t"} {
    split($1, a, ":")
    chrom = a[1]
    pos   = a[2]
    print chrom, pos-1, pos
}' ${PRUNE_IN} | sort -k1,1V -k2,2n \
  > ${OUT}/ld_pruned_regions.bed

echo "  LD-pruned region count: $(wc -l < ${OUT}/ld_pruned_regions.bed)"

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view \
    --regions-file ${OUT}/ld_pruned_regions.bed \
    ${OUT}/candidates_maf20.vcf.gz \
    -O z -o ${OUT}/candidates_maf20_ldpruned.vcf.gz 2>/dev/null

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${OUT}/candidates_maf20_ldpruned.vcf.gz

N_LD=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${OUT}/candidates_maf20_ldpruned.vcf.gz 2>/dev/null | wc -l)
echo "  After LD pruning: ${N_LD} SNPs"

# ---- Step C: Annotated summary of candidate set ----
echo ""
echo "[C] Impact distribution in candidate set..."
# Use annotated VCF as reference, intersect with candidate positions
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools isec \
    ${OUT}/candidates_maf20_ldpruned.vcf.gz \
    ${ANN_VCF} \
    -p ${OUT}/isec_annotated \
    -O z 2>/dev/null

echo ""
echo "=== Array design filter summary ==="
echo "  Input SNPs:         ${N_IN}"
echo "  After MAF ≥ 0.20:   ${N_MAF}"
echo "  After LD pruning:   ${N_LD}"
python3 -c "
n0=${N_IN}; nf=${N_LD}
print(f'  Retained fraction:  {nf/n0*100:.1f}%')
print(f'  Candidate markers:  {nf:,}')
"
echo "=== Array design filtering complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06b_array_design_filter.sh
```

---

## Step 3 — Genomic Context Filter: Remove Repetitive Regions

Probes in repetitive or low-complexity regions cross-hybridise to multiple
genomic locations, producing unreliable genotype calls. We use the TAIR10
repeat annotation to remove candidates in or near annotated repeats.

```bash
cat > ${VARCALL}/scripts/06c_context_filter.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=context_filter
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/context_filter_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/context_filter_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

GFF=${VARCALL}/00_reference/gff/Athaliana_TAIR10.gff3
CAND=${VARCALL}/08_chip_design/candidates/candidates_maf20_ldpruned.vcf.gz
OUT=${VARCALL}/08_chip_design/context

echo "=== Genomic context filtering: $(date) ==="

# ---- Extract repeat regions from GFF3 ----
# Use awk -F'\t' instead of grep -E with $(printf '\t'), which does not
# reliably produce literal tab characters inside extended regex patterns.
echo "Extracting repeat annotations from GFF3..."

awk -F'\t' '
  $3 == "transposable_element"      ||
  $3 == "repeat_region"             ||
  $3 == "transposable_element_gene" ||
  $3 == "LTR_retrotransposon"       ||
  $3 == "DNA_transposon" {
      print $1 "\t" ($4-1) "\t" $5
  }
' ${GFF} | sort -k1,1V -k2,2n > ${OUT}/repeats.bed

NREP=$(wc -l < ${OUT}/repeats.bed)
echo "  Repeat regions found: ${NREP}"

# ---- Spot-check what feature types are in the GFF3 ----
echo "  Feature types present in GFF3 (top 20):"
awk -F'\t' 'NR>1 && !/^#/ {print $3}' ${GFF} | \
  sort | uniq -c | sort -rn | head -20

# ---- Extend and merge repeat windows ----
if [[ ${NREP} -gt 0 ]]; then
    apptainer exec --bind ${VARCALL}:${VARCALL} \
      ${VARCALL}/containers/bedtools_2.31.1.sif \
      bedtools slop \
        -i ${OUT}/repeats.bed \
        -g ${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta.fai \
        -b 50 | \
    apptainer exec --bind ${VARCALL}:${VARCALL} \
      ${VARCALL}/containers/bedtools_2.31.1.sif \
      bedtools merge -i - \
      > ${OUT}/repeats_extended.bed

    echo "  Extended + merged repeat windows: $(wc -l < ${OUT}/repeats_extended.bed)"
else
    echo "  WARNING: No repeat features found with expected types."
    echo "  Check the feature type list above and update the awk pattern if needed."
    echo "  Proceeding without repeat filtering (all candidates retained)."
    touch ${OUT}/repeats_extended.bed
fi

# ---- Convert candidate VCF to BED ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' \
  ${CAND} 2>/dev/null > ${OUT}/candidates.bed

N_IN=$(wc -l < ${OUT}/candidates.bed)
echo ""
echo "Candidate SNPs before context filter: ${N_IN}"

# ---- Remove candidates overlapping extended repeat windows ----
# Guard: if repeats_extended.bed is empty, skip subtract and copy directly
NREP_EXT=$(wc -l < ${OUT}/repeats_extended.bed)

if [[ ${NREP_EXT} -gt 0 ]]; then
    apptainer exec --bind ${VARCALL}:${VARCALL} \
      ${VARCALL}/containers/bedtools_2.31.1.sif \
      bedtools subtract \
        -a ${OUT}/candidates.bed \
        -b ${OUT}/repeats_extended.bed \
      > ${OUT}/candidates_nonrepeat.bed
else
    cp ${OUT}/candidates.bed ${OUT}/candidates_nonrepeat.bed
fi

N_OUT=$(wc -l < ${OUT}/candidates_nonrepeat.bed)
echo "Candidate SNPs after removing repeat-proximal sites: ${N_OUT}"
python3 -c "
n0=${N_IN}; n1=${N_OUT}
removed = n0 - n1
pct = removed/n0*100 if n0 > 0 else 0
print(f'  Removed: {removed:,} SNPs ({pct:.1f}%)')
"

# ---- Filter VCF to non-repeat candidates ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view \
    --regions-file ${OUT}/candidates_nonrepeat.bed \
    ${CAND} \
    -O z -o ${OUT}/candidates_final.vcf.gz 2>/dev/null

apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${OUT}/candidates_final.vcf.gz

echo ""
echo "Final candidate VCF: $(ls -lh ${OUT}/candidates_final.vcf.gz)"
echo "=== Context filtering complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06c_context_filter.sh
```

---

## Step 4 — Extract Sequence Context for Probe and Primer Design

Every genotyping platform needs the flanking sequence around each SNP.
For SNP arrays the probe is 35–70 bp; for KASP it is 20 bp on each side;
for AmpSeq it is 200–500 bp around the SNP.

We extract 200 bp on each side, which covers all three platforms.

```bash
cat > ${VARCALL}/scripts/06d_extract_context.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=extract_context
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/extract_context_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/extract_context_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
CAND=${VARCALL}/08_chip_design/context/candidates_final.vcf.gz
OUT=${VARCALL}/08_chip_design/context

echo "=== Extracting sequence context: $(date) ==="

# ---- Create ±200 bp windows around each candidate SNP ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools query -f '%CHROM\t%POS0\t%END\t%CHROM:%POS[%REF>%ALT]\n' \
  ${CAND} 2>/dev/null | \
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bedtools_2.31.1.sif \
  bedtools slop \
    -i stdin \
    -g ${REF}.fai \
    -b 200 \
  > ${OUT}/context_windows.bed

echo "Context windows created: $(wc -l < ${OUT}/context_windows.bed)"

# ---- Extract FASTA sequences ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bedtools_2.31.1.sif \
  bedtools getfasta \
    -fi ${REF} \
    -bed ${OUT}/context_windows.bed \
    -name \
    -fo ${OUT}/candidate_contexts.fa

echo "Context sequences extracted."
echo "Total sequences: $(grep -c '^>' ${OUT}/candidate_contexts.fa)"

# ---- Basic sequence quality checks ----
echo ""
echo "=== Sequence context quality summary ==="

# Define path at shell level so it resolves before Python reads it
CONTEXT_FA=${OUT}/candidate_contexts.fa

python3 - << PYEOF
gc_vals = []
n_vals  = []
at_runs = 0
gc_runs = 0
total   = 0

with open("${CONTEXT_FA}") as f:
    seq = ""
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if seq:
                seq_u = seq.upper()
                l     = len(seq_u)
                gc    = (seq_u.count('G') + seq_u.count('C')) / l * 100
                n_pct = seq_u.count('N') / l * 100
                gc_vals.append(gc)
                n_vals.append(n_pct)
                if seq_u.count('AAAA') > 0 or seq_u.count('TTTT') > 0:
                    at_runs += 1
                if seq_u.count('CCCC') > 0 or seq_u.count('GGGG') > 0:
                    gc_runs += 1
                total += 1
            seq = ""
        else:
            seq += line
    if seq:
        seq_u = seq.upper()
        l     = len(seq_u)
        gc    = (seq_u.count('G') + seq_u.count('C')) / l * 100
        n_pct = seq_u.count('N') / l * 100
        gc_vals.append(gc); n_vals.append(n_pct)
        if seq_u.count('AAAA') > 0 or seq_u.count('TTTT') > 0: at_runs += 1
        if seq_u.count('CCCC') > 0 or seq_u.count('GGGG') > 0: gc_runs += 1
        total += 1

mean_gc  = sum(gc_vals) / total
low_gc   = sum(1 for g in gc_vals if g < 25)
high_gc  = sum(1 for g in gc_vals if g > 65)
has_n    = sum(1 for n in n_vals if n > 0)

print(f"  Total candidate sequences:      {total:,}")
print(f"  Mean GC content:                {mean_gc:.1f}%")
print(f"  Low GC (<25%) — probe concern:  {low_gc:,} ({low_gc/total*100:.1f}%)")
print(f"  High GC (>65%) — probe concern: {high_gc:,} ({high_gc/total*100:.1f}%)")
print(f"  Contain Ns (gap sequence):      {has_n:,} ({has_n/total*100:.1f}%)")
print(f"  AT homopolymer runs (AAAA/TTTT):{at_runs:,} ({at_runs/total*100:.1f}%)")
print(f"  GC homopolymer runs (CCCC/GGGG):{gc_runs:,} ({gc_runs/total*100:.1f}%)")
PYEOF

echo "=== Context extraction complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06d_extract_context.sh
```

---

## Step 5 — KASP Assay Design Criteria

**KASP (Kompetitive Allele Specific PCR)** is the most widely used
assay format for targeted SNP genotyping in plant breeding. Each assay
uses three primers: two allele-specific forward primers (one per allele)
and one common reverse primer. The two forward primers are tagged with
different FRET (fluorescence resonance energy transfer) reporters so the
two alleles emit different fluorescent signals.

### What makes a good KASP assay?

The following criteria are applied to the ±20 bp flanking sequence around
each candidate SNP:

| Criterion | Threshold | Reason |
|---|---|---|
| GC content of flanking 20 bp | 40–60% | Ensures stable primer annealing; too low/high reduces efficiency |
| No SNPs within 5 bp of target | Required | Nearby SNPs disrupt primer-template complementarity |
| No homopolymer runs ≥ 4 bp | Required | Poly-A/T runs cause primer slippage and false calls |
| Target SNP is biallelic | Required | KASP assays one substitution at a time |
| 3' end of allele-specific primer ends at SNP | Design rule | The 3' mismatch between alleles is what drives allele specificity |

```bash
cat > ${VARCALL}/scripts/06e_kasp_candidates.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=kasp_design
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/kasp_design_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/kasp_design_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

CAND=${VARCALL}/08_chip_design/context/candidates_final.vcf.gz
WGS_VCF=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
OUT=${VARCALL}/08_chip_design/kasp

echo "=== KASP candidate evaluation: $(date) ==="

# ---- Step A: Convert candidates to BED ----
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools query -f '%CHROM\t%POS0\t%END\t%CHROM:%POS:%REF:%ALT\n' \
  ${CAND} 2>/dev/null > ${OUT}/kasp_candidates.bed

N_CAND=$(wc -l < ${OUT}/kasp_candidates.bed)
echo "Total candidates: ${N_CAND}"

# ---- Step B: Pre-generate WGS SNP BED as a proper file ----
# Avoids process substitution inside apptainer exec, which cannot
# access host file descriptors from inside the container namespace.
echo "Generating WGS SNP BED file..."
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools query -f '%CHROM\t%POS0\t%END\n' \
  ${WGS_VCF} 2>/dev/null > ${OUT}/wgs_all_snps.bed

N_WGS=$(wc -l < ${OUT}/wgs_all_snps.bed)
echo "WGS SNPs in BED: ${N_WGS}"

# ---- Step C: Extend each candidate by ±5 bp, count overlapping SNPs ----
echo "Checking for nearby SNPs (±5 bp) per candidate..."
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bedtools_2.31.1.sif \
  bedtools slop \
    -i ${OUT}/kasp_candidates.bed \
    -g ${REF}.fai \
    -b 5 | \
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bedtools_2.31.1.sif \
  bedtools intersect \
    -a stdin \
    -b ${OUT}/wgs_all_snps.bed \
    -c \
  > ${OUT}/kasp_nearby_snp_counts.bed

echo "Nearby SNP count distribution:"
awk '{print $NF}' ${OUT}/kasp_nearby_snp_counts.bed | \
  sort -n | uniq -c | \
  awk '{printf "  %s nearby SNPs: %s candidates\n", $2, $1}'

# ---- Step D: Extract clean candidates ----
# A candidate counts itself — the ±5 bp window of a candidate always
# overlaps its own position in the WGS VCF, so the minimum count is 1.
# "Clean" means only 1 overlap (the candidate itself), not 0.
MAX_NEARBY=1

awk -v max=${MAX_NEARBY} '$NF <= max {print $1"\t"$2"\t"$3"\t"$4}' \
  ${OUT}/kasp_nearby_snp_counts.bed \
  > ${OUT}/kasp_clean_candidates.bed

N_CLEAN=$(wc -l < ${OUT}/kasp_clean_candidates.bed)
N_FLAGGED=$((N_CAND - N_CLEAN))

echo ""
echo "=== KASP design summary ==="
echo "  Total candidates:                ${N_CAND}"
echo "  Flagged (extra SNP within ±5 bp): ${N_FLAGGED}"
echo "  Clean KASP candidates:           ${N_CLEAN}"

if [[ ${N_CLEAN} -eq 0 ]]; then
    echo ""
    echo "  NOTE: No candidates passed the strict ±5 bp filter."
    echo "  This is expected with a dense marker set — consider"
    echo "  relaxing to ±3 bp or using MAX_NEARBY=2 above."
    echo "  Saving all candidates as fallback for downstream steps."
    cp ${OUT}/kasp_candidates.bed ${OUT}/kasp_clean_candidates.bed
    N_CLEAN=$(wc -l < ${OUT}/kasp_clean_candidates.bed)
    echo "  Fallback candidate count: ${N_CLEAN}"
else
    echo ""
    echo "First 10 clean KASP candidates:"
    head -10 ${OUT}/kasp_clean_candidates.bed
fi

echo "=== KASP design complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06e_kasp_candidates.sh
```

---

## Step 6 — AmpSeq Panel Design

**Amplicon Sequencing (AmpSeq)** sequences short PCR amplicons (typically
100–500 bp) containing one or more target SNPs. Unlike KASP, a single
AmpSeq amplicon can capture multiple linked SNPs, providing haplotype
information rather than just genotypes at individual sites. Unlike GBS,
AmpSeq is entirely targeted — only the amplicons you design are sequenced.

### AmpSeq primer design criteria

| Criterion | Guideline |
|---|---|
| Amplicon size | 150–400 bp (for 2×150 bp paired-end sequencing) |
| Primer length | 20–25 bp |
| Primer Tm | 58–65°C, ±2°C between forward and reverse |
| GC content | 40–60% for primers; avoid 3'-end G/C runs longer than 3 bp |
| SNP placement | Target SNP in the central 60% of the amplicon |
| SNPs per amplicon | 1–5 informative SNPs per amplicon for haplotyping |

```bash
cat > ${VARCALL}/scripts/06f_ampseq_design.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=ampseq_design
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/ampseq_design_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/ampseq_design_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

CAND=${VARCALL}/08_chip_design/context/candidates_final.vcf.gz
WGS_VCF=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
OUT=${VARCALL}/08_chip_design/ampseq

echo "=== AmpSeq amplicon design: $(date) ==="

# ---- Create 300 bp amplicon windows around each candidate ----
# The target SNP will sit at position 150/300 = 50% (centre of amplicon)
echo "[A] Creating 300 bp amplicon windows..."
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools query -f '%CHROM\t%POS0\t%END\t%CHROM:%POS:%REF:%ALT\n' \
  ${CAND} 2>/dev/null | \
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bedtools_2.31.1.sif \
  bedtools slop \
    -i stdin \
    -g ${REF}.fai \
    -b 150 \
  > ${OUT}/amplicon_windows.bed

echo "  Amplicon windows: $(wc -l < ${OUT}/amplicon_windows.bed)"

# ---- Count additional SNPs within each amplicon window ----
# This tells us how many informative SNPs each amplicon captures beyond
# the primary target — more SNPs = richer haplotype information
echo "[B] Counting SNPs per amplicon..."
apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bedtools_2.31.1.sif \
  bedtools intersect \
    -a ${OUT}/amplicon_windows.bed \
    -b <(apptainer exec --bind ${VARCALL}:${VARCALL} \
           ${VARCALL}/containers/bcftools_1.23.1.sif \
           bcftools query -f '%CHROM\t%POS0\t%END\n' \
           ${WGS_VCF} 2>/dev/null) \
    -c \
  > ${OUT}/amplicons_with_snp_counts.bed

echo ""
echo "=== AmpSeq amplicon statistics ==="
awk '{print $NF}' ${OUT}/amplicons_with_snp_counts.bed | \
  python3 -c "
import sys
from collections import Counter
counts = [int(l.strip()) for l in sys.stdin if l.strip()]
c = Counter(counts)
total = len(counts)
print(f'  Total amplicons:              {total:,}')
print(f'  Single SNP amplicons:         {c[1]:,} ({c[1]/total*100:.1f}%)')
print(f'  2-3 SNP amplicons (haplotype):{c[2]+c[3]:,} ({(c[2]+c[3])/total*100:.1f}%)')
print(f'  4-5 SNP amplicons:            {c[4]+c[5]:,} ({(c[4]+c[5])/total*100:.1f}%)')
print(f'  >5 SNP amplicons (complex):   {sum(v for k,v in c.items() if k>5):,}')
print(f'  Mean SNPs per amplicon:       {sum(counts)/total:.1f}')
"

# ---- Extract high-value amplicons (2-5 SNPs = haplotype informative) ----
awk '$NF >= 2 && $NF <= 5' ${OUT}/amplicons_with_snp_counts.bed \
  > ${OUT}/haplotype_amplicons.bed

N_HAP=$(wc -l < ${OUT}/haplotype_amplicons.bed)
echo ""
echo "Haplotype-informative amplicons (2-5 SNPs): ${N_HAP}"

echo ""
echo "First 10 haplotype-informative amplicons:"
head -10 ${OUT}/haplotype_amplicons.bed

echo ""
echo "These amplicons are candidates for:"
echo "  1. Primer3 primer design within the window coordinates"
echo "  2. Multiplex PCR optimisation (pool up to 50-100 amplicons)"
echo "  3. Sequencing on MiSeq/NextSeq 2×150 bp paired-end"
echo "  4. Haplotype reconstruction from read pairs spanning the amplicon"
echo "=== AmpSeq design complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06f_ampseq_design.sh
```

---

## Step 7 — Visualise the Final Marker Set

Two visualisations help evaluate the candidate marker set: a chromosome-level
marker density plot and a comparison of marker counts by functional category.

```bash
# ---- Create the R script ----
cat > ${VARCALL}/scripts/06g_plot_marker_set.R << 'EOF'
# HCS 7004 — Module 06: Marker set visualisation
# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name <- "Jonathan"

library(ggplot2)

base_dir  <- paste0("/fs/scratch/PAS3260/", user_name, "/variant_calling")
cand_bed  <- file.path(base_dir, "08_chip_design/context/candidates_nonrepeat.bed")
kasp_bed  <- file.path(base_dir, "08_chip_design/kasp/kasp_clean_candidates.bed")
plot_dir  <- file.path(base_dir, "08_chip_design/plots")
dir.create(plot_dir, showWarnings = FALSE)

# ---- Chromosome sizes (from TAIR10 .fai) ----
chr_sizes <- data.frame(
  chrom = c("1","2","3","4","5"),
  size  = c(30427671, 19698289, 23459830, 18585056, 26975502),
  stringsAsFactors = FALSE
)

# ---- Load candidate SNP positions ----
cand <- tryCatch(
  read.table(cand_bed, header = FALSE,
             col.names = c("chrom","start","end","id")),
  error = function(e) {
    cat("WARNING: candidate BED not found, using empty dataset\n")
    data.frame(chrom=character(), start=integer(), end=integer(), id=character())
  }
)

cat("Total candidate SNPs:", nrow(cand), "\n")

# ---- Bin SNPs into 1 Mb windows per chromosome ----
cand_nuclear <- cand[cand$chrom %in% chr_sizes$chrom, ]
cand_nuclear$bin_mb <- floor(cand_nuclear$start / 1e6)

density <- aggregate(start ~ chrom + bin_mb, data = cand_nuclear, FUN = length)
colnames(density) <- c("chrom", "bin_mb", "count")
density$chrom <- factor(density$chrom, levels = as.character(1:5))

# ---- Plot 1: Marker density across chromosomes ----
p1 <- ggplot(density, aes(x = bin_mb, y = count, fill = chrom)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ chrom, ncol = 1, scales = "free_x",
             labeller = labeller(chrom = function(x) paste0("Chr", x))) +
  xlab("Position (Mb)") +
  ylab("SNPs per 1 Mb window") +
  ggtitle("Candidate marker density across Arabidopsis chromosomes\n(MAF≥0.20, LD-pruned, non-repetitive)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path(plot_dir, "marker_density.pdf"), plot = p1, width = 9, height = 10)
ggsave(file.path(plot_dir, "marker_density.png"), plot = p1, width = 9, height = 10, dpi = 150)
cat("Marker density plot saved.\n")

# ---- Plot 2: Marker counts per chromosome ----
chr_counts <- aggregate(start ~ chrom, data = cand_nuclear, FUN = length)
colnames(chr_counts) <- c("chrom", "n_markers")
chr_counts <- merge(chr_counts, chr_sizes, by = "chrom")
chr_counts$density_per_mb <- chr_counts$n_markers / (chr_counts$size / 1e6)
chr_counts$chrom <- factor(chr_counts$chrom, levels = as.character(1:5))

p2 <- ggplot(chr_counts, aes(x = chrom, y = density_per_mb, fill = chrom)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_markers), vjust = -0.3, size = 3.5) +
  xlab("Chromosome") +
  ylab("Markers per Mb") +
  ggtitle("Candidate marker density per Mb by chromosome") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_discrete(labels = paste0("Chr", 1:5))

ggsave(file.path(plot_dir, "markers_per_chr.pdf"), plot = p2, width = 7, height = 5)
ggsave(file.path(plot_dir, "markers_per_chr.png"), plot = p2, width = 7, height = 5, dpi = 150)
cat("Per-chromosome plot saved.\n")

cat("\nAll plots saved to:", plot_dir, "\n")
EOF

# ---- Create the SLURM wrapper ----
cat > ${VARCALL}/scripts/06g_plot_marker_set_slurm.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=marker_plots
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/marker_plots_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/marker_plots_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

module load gcc/12.3.0 R/4.5.2

echo "=== Marker set plots: $(date) ==="
Rscript ${VARCALL}/scripts/06g_plot_marker_set.R
echo "=== Plots complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06g_plot_marker_set_slurm.sh
```

---

## Step 8 — Connecting to Genomic Selection

### How many markers are enough?

Genomic selection (GS) uses genome-wide marker data to predict breeding
values without requiring prior knowledge of which markers are causal.
The accuracy of genomic selection depends on marker density relative to
the LD decay of the population.

The rule of thumb is that you need enough markers so that **every QTL is
in LD (r² > 0.2) with at least one marker**. Given the extended LD in
*Arabidopsis* (~100–200 kb), the number of markers needed is much lower
than for outcrossing species:

```
Effective genome size (nuclear)  = 119 Mb
Mean LD block size (Arabidopsis) ≈ 50 kb (accession-dependent)
Minimum markers for GS           ≈ 119,000 / 50 = ~2,400 markers

For outcrossing species (e.g. maize):
Mean LD block size               ≈ 2 kb
Minimum markers for GS           ≈ 119,000 / 2  = ~60,000 markers
```

This explains why *Arabidopsis* GWAS studies using 200K–500K markers are
more than adequate, while maize GS programs routinely use 50,000–650,000
markers.

```bash
cat > ${VARCALL}/scripts/06h_marker_density_summary.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=marker_summary
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/fs/scratch/PAS3260${user_name}variant_calling/logs/marker_summary_%j.out
#SBATCH --error=/fs/scratch/PAS3260${user_name}variant_calling/logs/marker_summary_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

CAND=${VARCALL}/08_chip_design/context/candidates_final.vcf.gz
KASP=${VARCALL}/08_chip_design/kasp/kasp_clean_candidates.bed
AMPSEQ=${VARCALL}/08_chip_design/ampseq/haplotype_amplicons.bed

echo "=== Tutorial marker set summary ==="
echo "Genome: Arabidopsis thaliana TAIR10 (5 nuclear chromosomes, 119 Mb)"
echo "Samples: 5 accessions from the 1001 Genomes Project"
echo ""

N_WGS=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz \
  2>/dev/null | wc -l)

N_FILT=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${VARCALL}/07_population/filtered/wgs_filtered_final.vcf.gz \
  2>/dev/null | wc -l)

N_CAND=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${CAND} 2>/dev/null | wc -l)

N_KASP=$(wc -l < ${KASP} 2>/dev/null || echo 0)
N_AMP=$(wc -l < ${AMPSEQ} 2>/dev/null || echo 0)

python3 - << PYEOF
n_wgs   = ${N_WGS}
n_filt  = ${N_FILT}
n_cand  = ${N_CAND}
n_kasp  = ${N_KASP}
n_amp   = ${N_AMP}
genome  = 119146348  # nuclear bp

print(f"{'Stage':<45} {'SNPs':>10} {'per Mb':>10}")
print("-" * 67)
print(f"{'WGS PASS (Modules 02 hard filters)':<45} {n_wgs:>10,} {n_wgs/(genome/1e6):>10.0f}")
print(f"{'Population filtered (MAF, depth, GQ)':<45} {n_filt:>10,} {n_filt/(genome/1e6):>10.0f}")
print(f"{'Array design candidates (MAF≥0.20, LD, non-repeat)':<45} {n_cand:>10,} {n_cand/(genome/1e6):>10.0f}")
print(f"{'KASP-ready candidates':<45} {n_kasp:>10,} {n_kasp/(genome/1e6):>10.0f}")
print(f"{'AmpSeq haplotype amplicons':<45} {n_amp:>10,}  (amplicons)")
print()
print("Comparison to published platforms:")
print(f"  Arabidopsis Axiom array (1001G design): ~1,100,000 markers")
print(f"  Typical GBS panel (Arabidopsis):          ~50,000–200,000 markers")
print(f"  This tutorial 5-accession panel:          {n_cand:,} markers")
print()
print("Note: the small sample size (5 accessions) limits marker discovery.")
print("A full 1,135-accession analysis recovers >10 million SNPs before")
print("filtering, yielding ~500K–1M array-designable markers.")
PYEOF

echo "=== Marker summary complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/06h_marker_density_summary.sh
```

---

## Completion Checklist

```bash
echo "=== Module 06 completion check ==="

echo ""
echo "--- Candidate VCFs ---"
ls -lh ${VARCALL}/08_chip_design/candidates/candidates_maf20.vcf.gz
ls -lh ${VARCALL}/08_chip_design/candidates/candidates_maf20_ldpruned.vcf.gz
ls -lh ${VARCALL}/08_chip_design/context/candidates_final.vcf.gz

echo ""
echo "--- Repeat filter ---"
ls -lh ${VARCALL}/08_chip_design/context/repeats_extended.bed

echo ""
echo "--- Sequence context ---"
ls -lh ${VARCALL}/08_chip_design/context/candidate_contexts.fa
echo "Sequences: $(grep -c '^>' ${VARCALL}/08_chip_design/context/candidate_contexts.fa 2>/dev/null)"

echo ""
echo "--- KASP candidates ---"
ls -lh ${VARCALL}/08_chip_design/kasp/kasp_clean_candidates.bed
echo "Clean KASP candidates: $(wc -l < ${VARCALL}/08_chip_design/kasp/kasp_clean_candidates.bed 2>/dev/null)"

echo ""
echo "--- AmpSeq amplicons ---"
ls -lh ${VARCALL}/08_chip_design/ampseq/haplotype_amplicons.bed
echo "Haplotype amplicons: $(wc -l < ${VARCALL}/08_chip_design/ampseq/haplotype_amplicons.bed 2>/dev/null)"

echo ""
echo "--- Plots ---"
ls -lh ${VARCALL}/08_chip_design/plots/marker_density.pdf
ls -lh ${VARCALL}/08_chip_design/plots/markers_per_chr.pdf
```

---

## Discussion Questions

1. The array design filter raised the MAF threshold from 0.05 (used in
   population filtering in Module 05) to 0.20 for this five-sample panel.
   If you were designing an array for a breeding program with 200 accessions
   instead of 5, would you lower or raise the MAF threshold? How does sample
   size interact with the reliability of MAF estimates, and what is the
   minimum MAF that ensures a variant is genuinely polymorphic rather than
   a genotype call error?

2. The LD decay analysis in Module 05 showed that *Arabidopsis* maintains
   high r² (>0.2) over distances of 100 kb or more. If you used the
   LD-pruned marker set from Module 05 as an AmpSeq panel for a genomic
   selection program, what fraction of the genome would be "tagged" (i.e.
   within 50 kb of a marker) assuming the panel contains 2,000 markers?
   Does this coverage meet the threshold calculated in Step 8?

3. The GBS data in Module 03 was generated from ApeKI digestion. The
   ApeKI cut sites preferentially occur in GC-rich regions (GCWGC), which
   means GBS misses AT-rich intergenic and centromeric regions. Looking at
   the marker density plot from Step 7, are there chromosomal regions with
   notably lower candidate marker density? How does this compare to the
   known locations of centromeres in *Arabidopsis* (approximately at the
   middle of each chromosome)?

4. A KASP assay fails in your laboratory for one of the candidates produced
   in Step 5 — both allele-specific reactions produce the same fluorescent
   signal, and the genotype calls cluster in only one region of the
   FAM/HEX plot. List three distinct molecular explanations for this
   failure, and for each, describe what you would check in the sequence
   context data to diagnose the cause.

5. The tutorial used a five-accession panel to simulate the SNP discovery
   → chip design pipeline. The 1001 Genomes Consortium used 1,135
   accessions and recovered 10.7 million SNPs before filtering. Extrapolate
   from your five-accession results: approximately what fraction of those
   10.7 million SNPs do you expect your five accessions to have recovered?
   What biological principle — related to the site frequency spectrum of
   natural variation — explains this fraction?

---

## Tutorial Completion — What You Have Built

Congratulations on completing the full tutorial series. Working from raw
sequencing reads through to genotyping technology design, you have:

| Module | What you produced |
|---|---|
| 01 | Indexed reference genome, trimmed reads from three sequencing strategies |
| 02 | WGS GVCF per sample → joint-genotyped, hard-filtered cohort VCF |
| 03 | Demultiplexed GBS library → Stacks VCF + GATK GBS VCF + cross-pipeline comparison |
| 04 | STAR 2-pass alignments → SplitNCigarReads → RNA VCF + ASE read counts |
| 05 | SnpEff-annotated VCF → population-filtered markers → PCA + LD decay |
| 06 | Candidate marker set for SNP array, KASP panel, and AmpSeq design |

The five *Arabidopsis* accessions you worked with — Vas-0, Bez-9, Gen-8,
Mah-6, and Usa-0 — represent a geographic transect from the Iberian
Peninsula to East Asia. The variants you discovered and the marker set
you designed reflect real genetic diversity documented by the 1001 Genomes
Consortium, and the pipeline you executed mirrors the one used to design
the Arabidopsis Axiom array that is in use in research laboratories today.

The same pipeline, with appropriate adjustments to reference genome,
sample panel, and platform parameters, applies directly to wheat, soybean,
maize, strawberry, and any other plant species with a reference genome and
a population of sequenced accessions.

---

*← [Module 05: Variant Annotation & Filtering](05_variant_annotation.md)*
