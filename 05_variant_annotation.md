# Module 05 — Variant Annotation & Filtering

**HCS 7004 — Genome Analytics | The Ohio State University**

*← [Module 04: RNA-seq Variant Calling](04_rnaseq_variant_calling.md) | [Module 06: From SNPs to Genotyping Technologies](06_snp_applications.md) →*

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

1. Build a SnpEff database from a reference genome and GFF3 annotation
   and annotate a VCF with predicted functional consequences
2. Interpret SnpEff consequence categories and their relevance to
   downstream applications in plant breeding and population genomics
3. Apply population-level filters to a multi-sample VCF using VCFtools,
   correctly specifying thresholds for missing data, minor allele
   frequency, and depth
4. Merge the three track-specific VCFs into a single unified callset
   and assess cross-track concordance quantitatively
5. Perform principal component analysis (PCA) and linkage disequilibrium
   (LD) pruning using PLINK, and interpret the results in the context of
   *Arabidopsis* population structure
6. Produce publication-quality PCA and LD decay plots in R submitted
   as SLURM batch jobs

---

## Biological Background

### Why annotate variants?

A VCF file after hard filtering contains positional information — chromosome,
position, reference allele, alternate allele — but no biological interpretation.
**Functional annotation** maps each variant onto the gene model and predicts
its consequence:

- Does it fall in a coding exon? If so, does it change the amino acid?
- Is it synonymous (same amino acid) or nonsynonymous (amino acid change)?
- Is it in a splice site, UTR, promoter, or intergenic region?
- Does it introduce a premature stop codon or disrupt a splice donor?

These predictions drive prioritisation: a breeder interested in selecting for
a trait associated with a candidate gene cares most about coding variants,
especially nonsynonymous and loss-of-function alleles. A population geneticist
estimating selection pressures needs to distinguish synonymous from
nonsynonymous rates (dN/dS). A chip designer wants to avoid splice site
variants that affect probe hybridisation.

### Why filter at the population level?

The GATK hard filters in Modules 02–04 removed technically poor variant calls.
Population-level filters remove variants that are **uninformative or
unreliable for a specific analytical purpose**, even if technically sound:

| Filter | Purpose | Typical threshold |
|---|---|---|
| Minor allele frequency (MAF) | Remove singletons and ultra-rare variants that are likely errors in small cohorts | MAF ≥ 0.05 |
| Missing data rate | Remove sites where too many samples lack a genotype | Max missingness ≤ 0.2 (≥ 80% genotyped) |
| Mean depth | Remove low-depth sites where genotype calls are unreliable | Min mean depth ≥ 5× |
| Hardy-Weinberg equilibrium | Remove sites with extreme genotype frequency deviations | p ≥ 1×10⁻⁶ |
| LD pruning | Reduce redundancy for PCA and genomic selection | Window-based r² < 0.2 |

---

## Pipeline Overview

```
Module 02 output: 03_wgs/filtered/cohort_PASS.vcf.gz      ─┐
Module 03 output: 04_gbs/vcf/gatk_GBS_PASS.vcf.gz         ─┤
Module 04 output: 05_rnaseq/vcf/rna_cohort_PASS.vcf.gz    ─┘
                              ↓
Step 1: Build SnpEff database for TAIR10
                              ↓
Step 2: Annotate WGS PASS VCF with SnpEff
                              ↓
Step 3: Parse annotation and summarise consequence categories
                              ↓
Step 4: Population-level filtering with VCFtools (WGS VCF)
                              ↓
Step 5: Merge all three track VCFs into unified callset
                              ↓
Step 6: LD pruning with PLINK
                              ↓
Step 7: PCA with PLINK
                              ↓
Step 8: Calculate pairwise LD with PLINK (SLURM)
                              ↓
Step 9: PCA plot in R (SLURM)
                              ↓
Step 10: LD decay plot in R (SLURM)
```

---

## Step 0 — Extend Directory Structure

```bash
mkdir -p ${VARCALL}/06_annotation/{snpeff_db,annotated,stats}
mkdir -p ${VARCALL}/07_population/{filtered,merged,plink,plots}
```

---

## Step 1 — Build the SnpEff Database for TAIR10

SnpEff ships with many pre-built databases, but *Arabidopsis thaliana* TAIR10
is not always current. Building from our own FASTA and GFF3 ensures the
chromosome names match the VCFs exactly — a mismatch would silently annotate
zero variants.

> **GFF3 prefix fix.** The Ensembl Plants GFF3 encodes IDs with `gene:` and
> `transcript:` prefixes (e.g. `ID=transcript:AT1G01010.1`). SnpEff strips
> these from `ID` fields but not from `Parent` references, which causes
> `WARNING_TRANSCRIPT_NOT_FOUND` for every feature and ~50% START codon errors.
> The `sed` command below strips all prefixes before the database build.

```bash
cat > ${VARCALL}/scripts/05a_snpeff_build.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=snpeff_build
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/snpeff_build_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/snpeff_build_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

DB_DIR=${VARCALL}/06_annotation/snpeff_db
REF=${VARCALL}/00_reference/fasta/Athaliana_TAIR10.fasta
GFF=${VARCALL}/00_reference/gff/Athaliana_TAIR10.gff3
DB_NAME=AtTAIR10

echo "=== Building SnpEff database: $(date) ==="

# ---- Create SnpEff config and data directory structure ----
mkdir -p ${DB_DIR}/data/${DB_NAME}

cat > ${DB_DIR}/snpEff.config << CONF
# SnpEff configuration for TAIR10
data.dir = ${DB_DIR}/data

# Arabidopsis thaliana TAIR10 custom database
${DB_NAME}.genome : Arabidopsis_thaliana
${DB_NAME}.chromosomes : 1, 2, 3, 4, 5, Mt, Pt

# Mitochondrial codon table (plant mt uses Vertebrate_Mitochondrial
# as the closest available approximation in SnpEff)
${DB_NAME}.Mt.codonTable : Vertebrate_Mitochondrial
CONF

echo "SnpEff config:"
cat ${DB_DIR}/snpEff.config

# ---- Fix Ensembl Plants GFF3 ID prefixes for SnpEff compatibility ----
echo "Fixing GFF3 ID prefixes for SnpEff..."
sed \
  -e 's/ID=gene:/ID=/g' \
  -e 's/ID=transcript:/ID=/g' \
  -e 's/Parent=gene:/Parent=/g' \
  -e 's/Parent=transcript:/Parent=/g' \
  ${GFF} > ${DB_DIR}/data/${DB_NAME}/genes.gff

# || true prevents set -e from treating the broken pipe from head as fatal
echo "  Done. Spot-checking fix:"
grep "^1" ${DB_DIR}/data/${DB_NAME}/genes.gff | head -3 | cut -f9 || true

cp ${REF} ${DB_DIR}/data/${DB_NAME}/sequences.fa

echo ""
echo "=== Building database ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/snpeff_5.4.0c.sif \
  snpEff build \
    -config ${DB_DIR}/snpEff.config \
    -dataDir ${DB_DIR}/data \
    -gff3 \
    -noCheckCds \
    -noCheckProtein \
    -v \
    ${DB_NAME}

# -noCheckCds / -noCheckProtein: skip the optional post-build validation
#   that requires cds.fa and protein.fa files. These files are not
#   required for annotation; the database is complete without them.

echo ""
echo "=== Database build complete: $(date) ==="
echo "Database files:"
ls -lh ${DB_DIR}/data/${DB_NAME}/
EOF

sbatch ${VARCALL}/scripts/05a_snpeff_build.sh
```

### Expected build summary

After a successful build the `.err` log should show:

```
# Has protein coding info    : true
# Genes                      : ~27,000
# START codon errors         : < 1%    ← confirms GFF3 fix worked
```

Two remaining warnings are harmless and expected:
- `WARNING_FILE_NOT_FOUND: protein.fa` — suppressed by `-noCheckProtein`
- Mitochondrion codon table — resolved by adding the `Mt.codonTable` line to the config above

---

## Step 2 — Annotate the WGS VCF with SnpEff

We annotate the WGS PASS VCF as the primary dataset because it has the
broadest genomic coverage and the most variants. The GBS and RNA-seq VCFs
can be annotated with the same command by substituting the input path.

```bash
cat > ${VARCALL}/scripts/05b_snpeff_annotate.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=snpeff_ann
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/snpeff_ann_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/snpeff_ann_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

DB_DIR=${VARCALL}/06_annotation/snpeff_db
DB_NAME=AtTAIR10
IN_VCF=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
ANN_VCF=${VARCALL}/06_annotation/annotated/cohort_PASS_annotated.vcf.gz
STATS_DIR=${VARCALL}/06_annotation/stats

echo "=== SnpEff annotation: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/snpeff_5.4.0c.sif \
  snpEff ann \
    -config ${DB_DIR}/snpEff.config \
    -dataDir ${DB_DIR}/data \
    -v \
    -stats ${STATS_DIR}/snpeff_summary.html \
    -csvStats ${STATS_DIR}/snpeff_summary.csv \
    ${DB_NAME} \
    ${IN_VCF} | \
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view \
    -O z -o ${ANN_VCF}

# -dataDir explicitly passed so SnpEff locates the database inside the
#   container environment without relying on config file path resolution.
#   Without this flag, SnpEff may attempt to download the database if
#   the config-relative path resolves differently inside the container.

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${ANN_VCF}

echo ""
echo "Annotated VCF: $(ls -lh ${ANN_VCF})"
echo "Summary report: ${STATS_DIR}/snpeff_summary.html"
echo "=== Annotation complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05b_snpeff_annotate.sh
```

### Understanding the ANN field

SnpEff adds an `ANN` INFO field to each variant record. A single variant
can have multiple annotations — one per transcript affected. Each annotation
is pipe-delimited:

```
ANN=A|missense_variant|MODERATE|AT1G01010|AT1G01010.1|transcript|...
     ↑ alt   ↑ effect      ↑ impact  ↑ gene    ↑ transcript
```

The four **impact categories** are:

| Impact | Examples | Relevance |
|---|---|---|
| `HIGH` | Stop gained, frameshift, splice site disruption | Loss-of-function candidates; highest priority in breeding |
| `MODERATE` | Missense, in-frame indel | Functional change likely; protein structure affected |
| `LOW` | Synonymous, splice region (non-essential) | Likely neutral; useful as neutral markers |
| `MODIFIER` | Intergenic, intronic, UTR, downstream | Regulatory potential; most variants fall here |

---

## Step 3 — Summarise Annotation Consequences

Parse the annotated VCF to count variants by consequence category.

```bash
cat > ${VARCALL}/scripts/05c_parse_annotations.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=parse_ann
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/parse_ann_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/parse_ann_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

ANN_VCF=${VARCALL}/06_annotation/annotated/cohort_PASS_annotated.vcf.gz
STATS_DIR=${VARCALL}/06_annotation/stats

echo "=== Parsing SnpEff annotations: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
    ${ANN_VCF} 2>/dev/null | \
python3 - << 'PYEOF'
import sys
from collections import defaultdict

consequence_counts = defaultdict(int)
impact_counts      = defaultdict(int)
total_variants     = 0

for line in sys.stdin:
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 5:
        continue
    total_variants += 1
    ann_field = parts[4]

    # Take the first annotation (most severe consequence)
    first_ann = ann_field.split(',')[0]
    fields = first_ann.split('|')

    if len(fields) >= 3:
        consequence = fields[1]
        impact      = fields[2]
        consequence_counts[consequence] += 1
        impact_counts[impact]           += 1

print(f"\nTotal annotated variants: {total_variants:,}")

print(f"\n{'Impact category':20s}  {'Count':>10s}  {'Percent':>8s}")
print("-" * 44)
for impact in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']:
    n   = impact_counts.get(impact, 0)
    pct = n / total_variants * 100 if total_variants > 0 else 0
    print(f"{impact:20s}  {n:>10,}  {pct:>7.1f}%")

print(f"\n{'Consequence':35s}  {'Count':>10s}")
print("-" * 48)
for cons, n in sorted(consequence_counts.items(),
                       key=lambda x: -x[1])[:20]:
    print(f"{cons:35s}  {n:>10,}")
PYEOF

echo ""
echo "=== Annotation summary complete: $(date) ==="
echo "Full HTML report: ${STATS_DIR}/snpeff_summary.html"
EOF

sbatch ${VARCALL}/scripts/05c_parse_annotations.sh
```

### Expected consequence distribution for *Arabidopsis* WGS

| Impact | Expected fraction | Notes |
|---|---|---|
| MODIFIER | ~75–80% | Intergenic and intronic — the majority |
| LOW | ~12–15% | Predominantly synonymous coding variants |
| MODERATE | ~8–10% | Predominantly missense variants |
| HIGH | ~0.5–1% | Stop gained, frameshift, splice disruption |

---

## Step 4 — Population-Level Filtering with VCFtools

Filters are applied sequentially so you can see how many sites each
threshold removes.

```bash
cat > ${VARCALL}/scripts/05d_vcftools_filter.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=vcftools_filter
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/vcftools_filter_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/vcftools_filter_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

IN=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
OUT=${VARCALL}/07_population/filtered
PREFIX=${OUT}/wgs_filtered

echo "=== VCFtools population filtering: $(date) ==="

# ---- Step A: SNPs only, MAF ≥ 0.05, max missingness 20% ----
echo "[A] Applying MAF and missingness filters..."
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/vcftools_0.1.17.sif \
  vcftools \
    --gzvcf ${IN} \
    --remove-indels \
    --maf 0.05 \
    --max-missing 0.8 \
    --recode \
    --recode-INFO-all \
    --stdout 2>/dev/null | \
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -O z -o ${PREFIX}_step1_maf_miss.vcf.gz

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${PREFIX}_step1_maf_miss.vcf.gz

N1=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${PREFIX}_step1_maf_miss.vcf.gz 2>/dev/null | wc -l)
echo "  After MAF≥0.05 + max-missing 20%: ${N1} SNPs"

# ---- Step B: Minimum mean depth ≥ 5 ----
echo "[B] Applying minimum mean depth filter..."
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/vcftools_0.1.17.sif \
  vcftools \
    --gzvcf ${PREFIX}_step1_maf_miss.vcf.gz \
    --min-meanDP 5 \
    --recode \
    --recode-INFO-all \
    --stdout 2>/dev/null | \
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -O z -o ${PREFIX}_step2_depth.vcf.gz

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${PREFIX}_step2_depth.vcf.gz

N2=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${PREFIX}_step2_depth.vcf.gz 2>/dev/null | wc -l)
echo "  After min mean depth ≥ 5×: ${N2} SNPs"

# ---- Step C: Per-sample depth and genotype quality ----
echo "[C] Applying per-sample depth and GQ filters..."
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/vcftools_0.1.17.sif \
  vcftools \
    --gzvcf ${PREFIX}_step2_depth.vcf.gz \
    --minDP 3 \
    --minGQ 20 \
    --recode \
    --recode-INFO-all \
    --stdout 2>/dev/null | \
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -O z -o ${PREFIX}_step3_pergeno.vcf.gz

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${PREFIX}_step3_pergeno.vcf.gz

N3=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H ${PREFIX}_step3_pergeno.vcf.gz 2>/dev/null | wc -l)
echo "  After per-sample minDP 3 + minGQ 20: ${N3} SNPs"

# ---- Rename final filtered VCF ----
cp ${PREFIX}_step3_pergeno.vcf.gz     ${PREFIX}_final.vcf.gz
cp ${PREFIX}_step3_pergeno.vcf.gz.tbi ${PREFIX}_final.vcf.gz.tbi

# ---- Filtering summary ----
N0=$(apptainer exec --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools view -H -v snps ${IN} 2>/dev/null | wc -l)

echo ""
echo "=== Filtering summary ==="
echo "  Input SNPs (PASS):                    ${N0}"
echo "  After MAF + missingness:              ${N1}"
echo "  After mean depth:                     ${N2}"
echo "  After per-sample depth + GQ (final):  ${N3}"
python3 -c "
n0=${N0}; n3=${N3}
print(f'  Retention rate: {n3/n0*100:.1f}%')
"
echo "=== VCFtools filtering complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05d_vcftools_filter.sh
```

---

## Step 5 — Merge All Three Track VCFs

```bash
cat > ${VARCALL}/scripts/05e_merge_vcfs.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=merge_vcfs
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/merge_vcfs_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/merge_vcfs_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

WGS=${VARCALL}/03_wgs/filtered/cohort_PASS.vcf.gz
GBS=${VARCALL}/04_gbs/vcf/gatk_GBS_PASS.vcf.gz
RNA=${VARCALL}/05_rnaseq/vcf/rna_cohort_PASS.vcf.gz
OUT=${VARCALL}/07_population/merged/all_tracks_merged.vcf.gz

echo "=== Merging three-track VCFs: $(date) ==="

for VCF in ${WGS} ${GBS} ${RNA}; do
    [[ -f ${VCF} ]]     || { echo "ERROR: VCF not found: ${VCF}"; exit 1; }
    [[ -f ${VCF}.tbi ]] || { echo "ERROR: Index not found: ${VCF}.tbi"; exit 1; }
    echo "  OK: $(basename ${VCF})"
done

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools merge \
    ${WGS} ${GBS} ${RNA} \
    --merge all \
    --force-samples \
    -O z -o ${OUT} 2>/dev/null

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools index --tbi ${OUT}

echo ""
echo "Merged VCF: $(ls -lh ${OUT})"
echo ""
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/bcftools_1.23.1.sif \
  bcftools stats ${OUT} 2>/dev/null | grep "^SN" | grep -E "SNPs|indels|samples"

echo "=== Merge complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05e_merge_vcfs.sh
```

---

## Step 6 — LD Pruning with PLINK

```bash
cat > ${VARCALL}/scripts/05f_plink_ld_prune.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=plink_ld
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/plink_ld_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/plink_ld_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
IN=${VARCALL}/07_population/filtered/wgs_filtered_final.vcf.gz
PLINK_DIR=${VARCALL}/07_population/plink
PREFIX=${PLINK_DIR}/arabidopsis_wgs

echo "=== PLINK LD pruning: $(date) ==="

# ---- Step A: Convert VCF to PLINK binary format ----
echo "[A] Converting VCF to PLINK binary..."
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/plink_1.90b7.sif \
  plink \
    --vcf ${IN} \
    --make-bed \
    --out ${PREFIX} \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --double-id \
    --vcf-half-call missing

echo ""
echo "PLINK conversion log (last 10 lines):"
tail -10 ${PREFIX}.log

echo ""
echo "PLINK binary files:"
ls -lh ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam

# ---- Step B: LD pruning ----
echo ""
echo "[B] LD pruning (window 50 SNPs, step 10, r² < 0.2)..."
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/plink_1.90b7.sif \
  plink \
    --bfile ${PREFIX} \
    --indep-pairwise 50 10 0.2 \
    --out ${PREFIX}_pruned \
    --allow-extra-chr

echo ""
echo "PLINK LD pruning log (last 5 lines):"
tail -5 ${PREFIX}_pruned.log

TOTAL=$(wc -l < ${PREFIX}.bim)
KEPT=$(wc -l < ${PREFIX}_pruned.prune.in)
REMOVED=$(wc -l < ${PREFIX}_pruned.prune.out)
echo ""
echo "  Total SNPs:          ${TOTAL}"
echo "  Kept (prune.in):     ${KEPT}"
echo "  Removed (prune.out): ${REMOVED}"

# ---- Step C: Extract pruned SNP set ----
echo ""
echo "[C] Extracting pruned SNP set..."
apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/plink_1.90b7.sif \
  plink \
    --bfile ${PREFIX} \
    --extract ${PREFIX}_pruned.prune.in \
    --make-bed \
    --out ${PREFIX}_LDpruned \
    --allow-extra-chr

echo ""
echo "LD-pruned dataset:"
ls -lh ${PREFIX}_LDpruned.bed ${PREFIX}_LDpruned.bim ${PREFIX}_LDpruned.fam
echo "=== LD pruning complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05f_plink_ld_prune.sh
```

---

## Step 7 — PCA with PLINK

```bash
cat > ${VARCALL}/scripts/05g_plink_pca.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=plink_pca
#SBATCH --account=PAS3260
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/plink_pca_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/plink_pca_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
PLINK_DIR=${VARCALL}/07_population/plink
PREFIX=${PLINK_DIR}/arabidopsis_wgs_LDpruned

echo "=== PLINK PCA: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/plink_1.90b7.sif \
  plink \
    --bfile ${PREFIX} \
    --pca 4 \
    --out ${PLINK_DIR}/arabidopsis_pca \
    --allow-extra-chr

echo ""
echo "Eigenvalues (variance explained per PC):"
cat ${PLINK_DIR}/arabidopsis_pca.eigenval

echo ""
echo "Eigenvectors (sample coordinates):"
cat ${PLINK_DIR}/arabidopsis_pca.eigenvec

echo "=== PCA complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05g_plink_pca.sh
```

---

## Step 8 — Calculate Pairwise LD with PLINK

This step calculates all pairwise r² values within a 500 kb window
on the five nuclear chromosomes. The output file (`arabidopsis_ld.ld`)
will be several GB in size — it is processed in Step 10 by the R script
and can be deleted afterwards to save space.

```bash
cat > ${VARCALL}/scripts/05h_plink_ld_calc.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=plink_ld_calc
#SBATCH --account=PAS3260
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/plink_ld_calc_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/plink_ld_calc_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling
PLINK_DIR=${VARCALL}/07_population/plink

echo "=== PLINK pairwise LD calculation: $(date) ==="

apptainer exec \
  --bind ${VARCALL}:${VARCALL} \
  ${VARCALL}/containers/plink_1.90b7.sif \
  plink \
    --bfile ${PLINK_DIR}/arabidopsis_wgs \
    --r2 \
    --ld-window 100 \
    --ld-window-kb 500 \
    --ld-window-r2 0 \
    --chr 1 2 3 4 5 \
    --out ${PLINK_DIR}/arabidopsis_ld \
    --allow-extra-chr

# --chr 1 2 3 4 5: nuclear chromosomes only; excludes Mt and Pt which
#   have distinct LD structure and very few SNPs
# --ld-window-r2 0: report all pairs including r² = 0 so the decay
#   curve extends to the zero-LD baseline
# --ld-window-kb 500: maximum physical distance of 500 kb

echo ""
echo "LD output file:"
ls -lh ${PLINK_DIR}/arabidopsis_ld.ld

echo ""
echo "First 5 lines of LD file:"
head -5 ${PLINK_DIR}/arabidopsis_ld.ld

echo "=== LD calculation complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05h_plink_ld_calc.sh
```

---

## Step 9 — PCA Plot in R

Two files are created: the R script itself and the SLURM wrapper that
submits it. Replace `<your_username>` in **both** files.

```bash
# ---- Create the R script ----
cat > ${VARCALL}/scripts/05i_plot_pca.R << 'EOF'
# HCS 7004 — Module 05: PCA plot
# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name <- "<your_username>"

library(ggplot2)

plink_dir <- paste0("/fs/scratch/PAS3260/", user_name,
                    "/variant_calling/07_population/plink")
plot_dir  <- paste0("/fs/scratch/PAS3260/", user_name,
                    "/variant_calling/07_population/plots")
dir.create(plot_dir, showWarnings = FALSE)

cat("plink_dir:", plink_dir, "\n")
cat("plot_dir: ", plot_dir,  "\n")

# ---- Load eigenvectors and eigenvalues ----
pca <- read.table(
  file.path(plink_dir, "arabidopsis_pca.eigenvec"),
  header = FALSE,
  col.names = c("FID", "IID", paste0("PC", 1:4))
)

eigenval <- read.table(
  file.path(plink_dir, "arabidopsis_pca.eigenval"),
  header = FALSE
)$V1

pct_var <- round(eigenval / sum(eigenval) * 100, 1)
cat("Variance explained: PC1 =", pct_var[1], "%, PC2 =", pct_var[2], "%\n")

# ---- Geographic metadata ----
geo <- data.frame(
  IID    = c("Vas-0",        "Bez-9",      "Gen-8",      "Mah-6",          "Usa-0"),
  Region = c("Scandinavia",  "W. Europe",  "C. Europe",  "Mediterranean",  "East Asia"),
  stringsAsFactors = FALSE
)

pca <- merge(pca, geo, by = "IID")
cat("Samples in PCA:\n")
print(pca[, c("IID", "Region", "PC1", "PC2")])

# ---- Plot PC1 vs PC2 ----
p <- ggplot(pca, aes(x = PC1, y = PC2, colour = Region, label = IID)) +
  geom_point(size = 5) +
  geom_text(vjust = -0.8, size = 3.5) +
  xlab(paste0("PC1 (", pct_var[1], "% variance)")) +
  ylab(paste0("PC2 (", pct_var[2], "% variance)")) +
  ggtitle("PCA of five Arabidopsis thaliana accessions\n(WGS SNPs, LD-pruned)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right") +
  scale_colour_brewer(palette = "Set1")

ggsave(file.path(plot_dir, "pca_wgs.pdf"),
       plot = p, width = 7, height = 6)
ggsave(file.path(plot_dir, "pca_wgs.png"),
       plot = p, width = 7, height = 6, dpi = 150)

cat("PCA plots saved to:", plot_dir, "\n")
EOF

# ---- Create the SLURM wrapper ----
cat > ${VARCALL}/scripts/05i_plot_pca_slurm.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=pca_plot
#SBATCH --account=PAS3260
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/pca_plot_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/pca_plot_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

module load gcc/12.3.0 R/4.5.2

echo "=== PCA plot: $(date) ==="
Rscript ${VARCALL}/scripts/05i_plot_pca.R
echo "=== PCA plot complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05i_plot_pca_slurm.sh
```

---

## Step 10 — LD Decay Plot in R

> **Note on file size.** The `arabidopsis_ld.ld` file produced by Step 8
> is several GB. Reading it in R requires ~32 GB of memory, which is why
> this step runs as a SLURM job rather than interactively. Once the plot
> is produced, the `.ld` file can be deleted to recover scratch space.

```bash
# ---- Create the R script ----
cat > ${VARCALL}/scripts/05j_plot_ld_decay.R << 'EOF'
# HCS 7004 — Module 05: LD decay plot
# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name <- "<your_username>"

library(ggplot2)

plink_dir <- paste0("/fs/scratch/PAS3260/", user_name,
                    "/variant_calling/07_population/plink")
plot_dir  <- paste0("/fs/scratch/PAS3260/", user_name,
                    "/variant_calling/07_population/plots")
dir.create(plot_dir, showWarnings = FALSE)

# ---- Load pairwise r² ----
cat("Reading LD file (this may take several minutes)...\n")
ld_file <- file.path(plink_dir, "arabidopsis_ld.ld")
cat("File size:", round(file.info(ld_file)$size / 1e9, 2), "GB\n")

ld <- read.table(ld_file, header = TRUE)
cat("Loaded", nrow(ld), "pairwise r² values\n")

# Physical distance in kb
ld$dist_kb <- abs(ld$BP_B - ld$BP_A) / 1000

# Remove same-position pairs (distance = 0)
ld <- ld[ld$dist_kb > 0, ]

# Bin by distance and compute mean r² with standard error
ld$bin_kb <- cut(ld$dist_kb,
                 breaks = c(0, 5, 10, 20, 50, 100, 200, 500),
                 labels = c("0-5", "5-10", "10-20", "20-50",
                            "50-100", "100-200", "200-500"))

ld_summary <- aggregate(R2 ~ bin_kb, data = ld,
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                            se   = sd(x, na.rm = TRUE) /
                                                   sqrt(sum(!is.na(x)))))
ld_summary <- do.call(data.frame, ld_summary)
colnames(ld_summary) <- c("bin_kb", "mean_r2", "se_r2")

cat("\nLD summary by distance bin:\n")
print(ld_summary)

# ---- Plot ----
p <- ggplot(ld_summary, aes(x = bin_kb, y = mean_r2, group = 1)) +
  geom_ribbon(aes(ymin = pmax(0, mean_r2 - se_r2),
                  ymax = pmin(1, mean_r2 + se_r2)),
              alpha = 0.2, fill = "steelblue") +
  geom_line(colour = "steelblue", linewidth = 1) +
  geom_point(colour = "steelblue", size = 2.5) +
  xlab("Physical distance (kb)") +
  ylab(expression(Mean ~ r^2)) +
  ggtitle("LD decay across five Arabidopsis accessions\n(WGS SNPs, nuclear chromosomes 1-5)") +
  theme_bw(base_size = 13) +
  ylim(0, 1)

ggsave(file.path(plot_dir, "ld_decay.pdf"),
       plot = p, width = 7, height = 5)
ggsave(file.path(plot_dir, "ld_decay.png"),
       plot = p, width = 7, height = 5, dpi = 150)

cat("LD decay plots saved to:", plot_dir, "\n")

# ---- Optional: remove the large .ld file to recover disk space ----
# Uncomment the line below once you have confirmed the plot looks correct
# file.remove(ld_file)
EOF

# ---- Create the SLURM wrapper ----
cat > ${VARCALL}/scripts/05j_plot_ld_decay_slurm.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=ld_decay_plot
#SBATCH --account=PAS3260
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/ld_decay_plot_%j.out
#SBATCH --error=/fs/scratch/PAS3260/${user_name}/variant_calling/logs/ld_decay_plot_%j.err

set -euo pipefail

# ↓ Replace <your_username> with your OSC username (e.g. jsmith)
user_name=<your_username>
VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling

module load gcc/12.3.0 R/4.5.2

echo "=== LD decay plot: $(date) ==="
Rscript ${VARCALL}/scripts/05j_plot_ld_decay.R
echo "=== LD decay plot complete: $(date) ==="
EOF

sbatch ${VARCALL}/scripts/05j_plot_ld_decay_slurm.sh
```

> **Interpreting LD decay in *Arabidopsis*.** Self-fertilising species like
> *A. thaliana* maintain high LD over long distances because selfing means
> effective recombination is rare. You should see r² remaining above
> 0.2–0.3 even at 100 kb distances — a stark contrast to outcrossing
> species like maize where r² drops below 0.1 within 2–5 kb. This extended
> LD is both a challenge (fewer independent markers per genomic region) and
> an advantage (fewer markers needed for GWAS) compared to outcrossers.

---

## Completion Checklist

```bash
echo "=== Module 05 completion check ==="

echo ""
echo "--- SnpEff database ---"
ls ${VARCALL}/06_annotation/snpeff_db/data/AtTAIR10/snpEffectPredictor.bin \
  2>/dev/null && echo "Database: OK"

echo ""
echo "--- Annotated WGS VCF ---"
ls -lh ${VARCALL}/06_annotation/annotated/cohort_PASS_annotated.vcf.gz

echo ""
echo "--- SnpEff summary report ---"
ls -lh ${VARCALL}/06_annotation/stats/snpeff_summary.html

echo ""
echo "--- Population-filtered VCF ---"
ls -lh ${VARCALL}/07_population/filtered/wgs_filtered_final.vcf.gz

echo ""
echo "--- Merged three-track VCF ---"
ls -lh ${VARCALL}/07_population/merged/all_tracks_merged.vcf.gz

echo ""
echo "--- PLINK LD-pruned files ---"
ls -lh ${VARCALL}/07_population/plink/arabidopsis_wgs_LDpruned.bed

echo ""
echo "--- PCA output ---"
ls -lh ${VARCALL}/07_population/plink/arabidopsis_pca.eigenvec
ls -lh ${VARCALL}/07_population/plink/arabidopsis_pca.eigenval

echo ""
echo "--- R plots ---"
ls -lh ${VARCALL}/07_population/plots/pca_wgs.pdf
ls -lh ${VARCALL}/07_population/plots/ld_decay.pdf
```

---

## Discussion Questions

1. SnpEff classified approximately 0.5–1% of your WGS variants as HIGH
   impact. If you were designing a follow-up experiment to validate
   candidate loss-of-function alleles in Usa-0, which of the HIGH impact
   consequence categories would you prioritise, and what additional
   information (beyond the SnpEff annotation) would you need to confirm
   functional impact?

2. The VCFtools `--max-missing 0.8` filter requires that a site be
   genotyped in at least 80% of samples — that is, at least 4 of 5
   accessions. Given that Gen-8 has the lowest WGS coverage (~14×),
   which chromosome regions do you expect to contribute disproportionately
   to sites failing this filter? How would you test this hypothesis using
   `vcftools --missing-site`?

3. PCA on five samples produces at most four non-zero principal components.
   PC1 and PC2 together likely explain >80% of the variance. What does it
   mean biologically if PC1 separates Usa-0 from the four European
   accessions, and PC2 separates the Mediterranean accession (Mah-6)
   from the northern European ones? How do these axes of variation relate
   to the geographic origins in the accession table?

4. LD pruning with r² < 0.2 in a 50-SNP window removes a large fraction
   of SNPs before PCA. If you repeated the PCA without LD pruning, which
   chromosomal regions would likely dominate the first few PCs, and why?
   What biological feature of the *Arabidopsis* genome makes this
   particularly important?

5. Compare the MODIFIER variant count between the WGS and RNA-seq
   annotated VCFs. The WGS set has far more MODIFIER variants. A student
   claims this proves that MODIFIER variants are non-functional. Explain
   why this comparison does not support that conclusion, and what it
   actually reflects about the two datasets.

---

*← [Module 04: RNA-seq Variant Calling](04_rnaseq_variant_calling.md) | [Module 06: From SNPs to Genotyping Technologies](06_snp_applications.md) →*
