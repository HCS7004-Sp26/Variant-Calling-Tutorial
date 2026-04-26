# Module 00 — Conceptual Framework

**HCS 7004 — Genome Analytics | The Ohio State University**

*[Module 01: Environment Setup →](01_environment_setup.md)*

---

> ### 📋 How to use this module
> Module 00 is a reading and discussion module — no commands are run here.
> It establishes the biological and computational framework that all
> subsequent modules build on. Work through it before starting Module 01.
>
> The four environment variables below are set up in Module 01 and used
> throughout the entire tutorial series. They are listed here so you know
> what to expect:
>
> ```bash
> export user_name=<your_username>                                    # your OSC login
> export VARCALL=/fs/scratch/PAS3260/${user_name}/variant_calling     # your working directory
> export SHARED=/fs/scratch/PAS3260/Variant_Calling/Data             # shared read-only data
> export CONTAINERS=/fs/scratch/PAS3260/Variant_Calling/Containers   # shared containers
> ```

---

## Learning Objectives

By the end of this tutorial series you will be able to:

1. Explain the biological origin of SNPs, indels, and structural variants and
   the statistical models used to distinguish true variants from sequencing error
2. Execute the GATK4 GVCF pipeline for WGS data, including per-sample
   HaplotypeCaller, joint genotyping, and hard-filter calibration
3. Contrast WGS, GBS, and RNA-seq variant calling in terms of genome coverage,
   cost, and appropriate use cases
4. Implement the GATK Best Practices workflow for RNA-seq variant calling,
   including the SplitNCigarReads correction step unique to spliced alignments
5. Annotate VCF files with SnpEff to predict the functional impact of each
   variant against the TAIR10 gene model set
6. Filter multi-sample VCFs for population genomic analysis using BCFtools
   and VCFtools
7. Explain how population-level variant data are used to design SNP genotyping
   arrays, AmpSeq panels, and KASP assays
8. Connect variant calling quality metrics to the performance of downstream
   genotyping technologies used in plant breeding programs

---

## Part 1 — The Biology of Genetic Variation

### What is a SNP — really?

A **single nucleotide polymorphism (SNP)** is a position in the genome where two
or more alleles exist in a population at a frequency above 1%. From a sequencing
perspective, a SNP is a position where aligned reads from one individual
consistently carry a different nucleotide than the reference genome.

The challenge is not finding mismatches — it is deciding which mismatches are
real. Every sequencing read carries errors, and every alignment carries
ambiguities. Variant callers must distinguish a **true variant** from:

- **Sequencing error** — Illumina base-call error rate is ~0.1–1% per base;
  at 15× depth, a random error will appear in ~1–15 reads
- **Mapping error** — reads from paralogous or repetitive regions may align
  to the wrong locus, creating apparent mismatches at the mapping site
- **PCR error and duplicate reads** — amplification during library preparation
  can fix errors and inflate allele counts
- **Systematic instrument bias** — position-specific quality degradation, often
  visible at the 3' end of reads

This is why variant callers do not simply count mismatches. They build a
**probabilistic model** of the evidence — incorporating base quality scores,
mapping quality, read strand, allele balance, and local haplotype context —
before emitting a genotype likelihood.

### SNPs, indels, and structural variants

Variant types form a size continuum:

| Type | Size | Detection strategy |
|---|---|---|
| SNP | 1 bp substitution | Read pile-up; local re-assembly |
| Small indel | 1–50 bp insertion or deletion | Local re-assembly (GATK); gap-aware alignment |
| Copy number variant | 1 kb – 1 Mb | Read-depth analysis; paired-end span |
| Structural variant | > 50 bp rearrangement | Split reads; paired-end discordance; long reads |

This tutorial focuses on **SNPs and small indels**, which are the variant
classes recovered by Illumina short-read sequencing and the primary input
to all genotyping array technologies.

### How *Arabidopsis thaliana* accumulates variants

*Arabidopsis thaliana* is a **predominantly self-fertilising** annual plant.
Because it reproduces almost entirely by selfing, each natural accession is
nearly completely homozygous — a natural inbred line. This has two important
consequences for variant calling:

**1. True heterozygous sites are rare.** When you observe a heterozygous
genotype call in your VCF, it deserves scrutiny. It may reflect a genuine
recent mutation, a residual heterozygous region, or — more commonly — a
mapping artefact at a repetitive or duplicated locus.

**2. Between-accession differences are large.** Although each accession is
internally homozygous, different accessions that have been geographically
isolated for thousands of generations have accumulated substantial divergence.
The 1001 Genomes Project identified ~10.7 million SNPs across 1,135 accessions,
with individual accessions carrying ~700,000–4,000,000 SNPs relative to the
Col-0 reference depending on their geographic and phylogenetic distance from
the reference ecotype.

> **Why does this matter for your analysis?** When you call variants by mapping
> reads from Bez-9 or Usa-0 against the TAIR10 reference (which represents
> Col-0), you are measuring the genetic distance between that accession and
> Col-0 — not within-accession diversity. Keep this in mind when interpreting
> heterozygosity statistics.

---

## Part 2 — The GATK4 HaplotypeCaller Model

GATK4's `HaplotypeCaller` is the variant caller used in all three sequencing
tracks in this tutorial (with modifications for the RNA-seq track). Understanding
its logic will help you interpret its output and make informed filtering decisions.

### Local re-assembly: why it matters

Earlier variant callers (e.g., SAMtools mpileup) operated in **pileup mode**:
at each genomic position, they counted the bases present in overlapping reads.
This approach fails near indels because indels shift the reading frame of
downstream bases, creating cascading mismatches that obscure the true variant.

HaplotypeCaller instead performs **local re-assembly**:

```
Step 1 — Identify "active regions"
         Regions where read evidence suggests variation (elevated mismatch rate,
         soft-clipped reads, reads with indels)

Step 2 — Build a De Bruijn graph
         Decompose reads in the active region into k-mers; the graph represents
         all possible haplotypes supported by the data

Step 3 — Determine candidate haplotypes
         Traverse the graph to enumerate distinct haplotype sequences

Step 4 — Re-align reads to each haplotype (pair-HMM)
         Assign each read a likelihood score for each candidate haplotype,
         accounting for base quality, indel penalties, and local context

Step 5 — Calculate genotype likelihoods
         For each possible genotype (e.g., AA, AB, BB), compute the probability
         of observing the read data given that genotype

Step 6 — Emit GVCF
         Report the most likely genotype and its likelihood at every position,
         including positions with no variant evidence (reference blocks)
```

This approach recovers indels that pileup callers miss entirely, reduces false
positives at repetitive loci by modeling the full local sequence context, and
produces calibrated genotype likelihoods that downstream tools (GenomicsDBImport,
GenotypeGVCFs) use for joint genotyping.

### The GVCF workflow and why we use it

For population-scale analysis, GATK4 uses a **two-step genotyping strategy**:

```
Per-sample step (parallelisable):
  Each sample → HaplotypeCaller → sample.g.vcf.gz (GVCF)

Joint genotyping step (run once):
  All GVCFs → GenomicsDBImport → GenotypeGVCFs → cohort.vcf.gz
```

The GVCF contains **reference confidence blocks** — compressed records stating
that positions in a region have no evidence of variation above a likelihood
threshold. This means that when samples are jointly genotyped, the caller can
distinguish a truly homozygous-reference site from a site that simply had no
reads (missing data). This distinction is critical for population genomic
analyses where missing data patterns can confound allele frequency estimates.

---

## Part 3 — Three Sequencing Strategies, One VCF

### Why three approaches?

```
WGS  ──────────────────────────────────────────────────── Full genome
     Every position covered (at sufficient depth).
     Dense, uniform SNP discovery. Most expensive per sample.
     Best for: small-N studies, structural variation, de novo mutation,
               candidate gene resequencing, reference-quality variant sets.

GBS  ────┬───────┬───────┬──────────────────────────── Sampled genome
         ↑cut    ↑cut    ↑cut site
     Sequences only near ApeKI (GCWGC) restriction sites.
     ~15–25% of the genome; consistent across individuals.
     Low per-sample cost; scalable to thousands of individuals.
     Best for: GWAS panels, genomic selection training sets,
               population structure, marker-assisted breeding.

RNA  ────────┬──┬──────────┬──┬───────────────────── Expressed regions
             │  │          │  │
             ↑exon         ↑exon
     Sequences only transcribed, polyadenylated RNA.
     Variants detected only where genes are expressed.
     No additional cost when RNA-seq already planned.
     Best for: eQTL discovery, allele-specific expression,
               connecting genetic variation to gene regulation.
```

All three strategies produce a **VCF (Variant Call Format)** file — the
universal representation of genetic variation. All downstream analyses (GWAS,
QTL mapping, population structure, genomic selection, chip design) consume VCF
files and are largely agnostic to how the VCF was produced.

### Key differences to keep in mind throughout the tutorial

| Property | WGS | GBS | RNA-seq |
|---|---|---|---|
| Genome coverage | ~100% (at depth) | ~15–25% | Expressed regions only |
| Read depth per locus | Uniform (function of sequencing depth) | High near cut sites; zero elsewhere | Proportional to expression level |
| Library type in this tutorial | Paired-end 101 bp | Single-end 101 bp | Single-end 100 bp |
| Aligner used | BWA-MEM2 | BWA-MEM2 | STAR (splice-aware) |
| Special pre-processing | Duplicate marking (Picard) | None | SplitNCigarReads (GATK) |
| Indel sensitivity | High | Moderate | Lower (splicing complicates indels) |
| Can detect structural variants | Yes (paired-end span) | No | No |

### The VCF → Genotyping Technology pipeline

The SNPs on any commercial or custom genotyping platform originate from exactly
the kind of analysis you will perform in this tutorial:

```
Population resequencing (WGS or GBS)
          ↓
  Joint variant calling (GATK4 GenotypeGVCFs)
          ↓
  Quality filtering (VQSR or hard filters)
          ↓
  Population-level filtering (MAF, missingness, depth)
          ↓
  LD pruning (remove redundant markers)
          ↓
  Functional annotation (prioritise coding / regulatory SNPs)
          ↓
  Probe design (uniqueness filter; Tm QC for KASP)
          ↓
  SNP chip  /  AmpSeq panel  /  KASP assay
```

In Module 06 you will walk through this logic using the 1001 Genomes variant
data to understand how the *Arabidopsis* Infinium array was designed — and apply
the same principles to a crop breeding scenario relevant to your own research.

---

## Part 4 — The Study Organism

### *Arabidopsis thaliana* genome statistics (TAIR10 reference)

| Feature | Value |
|---|---|
| Total assembly size | ~157 Mb (nuclear + Mt + Pt) |
| Nuclear genome size | ~135 Mb |
| Chromosomes | 5 nuclear + mitochondrion (Mt) + chloroplast (Pt) |
| Protein-coding genes | ~27,600 |
| Repeat content | ~14% |
| Mean gene density | ~1 gene per 4.5 kb |
| GC content | ~36% |
| Annotation | TAIR10 — community gold standard, continuously updated |

### Why *Arabidopsis* is ideal for this tutorial

**1. Uniform instrument and library type across all five accessions.** All WGS
data used here were sequenced on the Illumina HiSeq 2000 at 101 bp paired-end
by the 1001 Genomes Consortium. All RNA-seq data were generated on the Illumina
HiSeq 2500 at 100 bp single-end by the 1001 Epigenomes Project. Instrument
uniformity removes a major source of batch effect that complicates interpretation
of multi-sample variant calls.

**2. WGS and RNA-seq data exist for exactly the same five accessions from the
same biological material.** This means students can directly compare which SNPs
are recovered by WGS vs. GBS vs. RNA-seq at the same genomic positions in the
same genetic backgrounds — a comparison impossible with most organisms.

**3. The genome is small enough to compute in tutorial time.** At 135 Mb with
modest repeat content, alignment and GATK variant calling complete within
30–90 minutes of SLURM queue time for a single sample. The full five-sample
cohort is feasible within a single tutorial session.

**4. Real plant breeding relevance.** The population structure, environmental
adaptation signals, and variant density patterns visible in these five accessions
mirror what students will encounter in crop genomics work. The principles apply
directly to wheat, soybean, maize, and strawberry breeding programs.

---

## Part 5 — Datasets Used in This Tutorial

### Track A — WGS

Five natural *Arabidopsis thaliana* accessions resequenced by the **1001 Genomes
Consortium** (The 1001 Genomes Consortium, *Cell* 166:481–491, 2016).
All data are paired-end 101 bp from the Illumina HiSeq 2000.

| Accession | 1001G ID | SRR | Coverage | Geographic origin |
|---|---|---|---|---|
| Vas-0 | 9902 | SRR1946455 | ~32× | Västerås, Sweden |
| Bez-9 | 9928 | SRR1946475 | ~25× | Besançon, France |
| Gen-8 | 9909 | SRR1946461 | ~14× | Geneva, Switzerland |
| Mah-6 | 9906 | SRR1946459 | ~17× | Mahón, Menorca, Spain |
| Usa-0 | 9929 | SRR1946457 | ~18× | Utsunomiya, Japan |

The five accessions span from the Atlantic coast of the Iberian Peninsula
(Mah-6) across central Europe (Bez-9, Gen-8) to Scandinavia (Vas-0) and East
Asia (Usa-0), sampling the broad west-east geographic gradient of *Arabidopsis*
diversity. Coverage ranges from ~14× (Gen-8) to ~32× (Vas-0), reflecting
natural variation in library complexity rather than intentional design — this
is realistic for any real-world sequencing project and lets students observe
how sequencing depth influences genotype confidence.

> **WGS data files are pre-staged at:**
> `/fs/scratch/PAS3260/Variant_Calling/Data/WGS/`

### Track B — GBS (in silico ApeKI simulation)

Rather than downloading a separate GBS dataset, the GBS data in this tutorial
are **computationally simulated** from the WGS reads using an in silico
restriction enzyme digestion with ApeKI (recognition sequence: GCWGC, where
W = A or T).

**Why simulate rather than use real GBS data?**

*Arabidopsis* is a model organism primarily studied by WGS and RNA-seq; large-scale
GBS datasets comparable in quality to the 1001 Genomes WGS data do not exist
for the same accessions. More importantly, the simulation approach delivers
something no real GBS dataset can: GBS reads from the **identical five accessions**
used in the WGS track, enabling direct position-by-position comparison of which
variants are recovered by each strategy.

**How the simulation works:**

1. All ApeKI cut sites (GCAGC and GCTGC motifs) in TAIR10 are identified
   (~80,000–100,000 sites genome-wide)
2. A ±200 bp window is drawn around each cut site and overlapping windows merged
   (~15–25% of the nuclear genome is covered)
3. WGS reads falling within these windows are extracted and the R1 reads written
   as single-end FASTQ — mimicking the SE sequencing layout of real GBS libraries

The resulting files contain reads of realistic length (101 bp, inherited from the
WGS source data) that are uniformly distributed across ApeKI windows. The one
honest limitation is that simulated reads do not begin with the CWG cut-site
overhang that real ApeKI reads do after adapter trimming — students are encouraged
to think about what this means for the "Per base sequence content" FastQC module.

| Accession | File | Size |
|---|---|---|
| Vas-0 | `Vas-0_GBS.fastq.gz` | 275 MB |
| Bez-9 | `Bez-9_GBS.fastq.gz` | 232 MB |
| Gen-8 | `Gen-8_GBS.fastq.gz` | 152 MB |
| Mah-6 | `Mah-6_GBS.fastq.gz` | 222 MB |
| Usa-0 | `Usa-0_GBS.fastq.gz` | 496 MB |

> **GBS data files are pre-staged at:**
> `/fs/scratch/PAS3260/Variant_Calling/Data/GBS/`
> The ApeKI cut site BED file (`ApeKI_GBS_windows.bed`) is also present and
> will be used in Module 03 to visualise coverage patterns.

### Track C — RNA-seq

Transcriptome data for the same five accessions from the **1001 Epigenomes
Project** (Kawakatsu et al., *Cell* 166:492–505, 2016). Libraries were prepared
from rosette leaf tissue harvested just before bolting using the TruSeq Stranded
RNA kit (strand-specific). All data are single-end 100 bp from the Illumina
HiSeq 2500. Multiple sequencing runs per accession were concatenated into a
single merged file per sample.

| Accession | SRX | Merged file | Size | Approx. read count |
|---|---|---|---|---|
| Vas-0 | SRX1735056 | `Vas-0_rnaseq.fastq.gz` | 763 MB | ~7.7M |
| Bez-9 | SRX1735072 | `Bez-9_rnaseq.fastq.gz` | 2.6 GB | ~26M |
| Gen-8 | SRX1735063 | `Gen-8_rnaseq.fastq.gz` | 3.2 GB | ~35M |
| Mah-6 | SRX1735060 | `Mah-6_rnaseq.fastq.gz` | 2.9 GB | ~30M |
| Usa-0 | SRX1735058 | `Usa-0_rnaseq.fastq.gz` | 3.1 GB | ~32M |

> ⚠️ **Note on Vas-0 RNA-seq depth.** The Vas-0 RNA-seq file is substantially
> smaller than the other four (~763 MB vs. 2.6–3.2 GB). This reflects the
> original SRA submission, where several sequencing runs for this accession
> were either not deposited or failed the minimum read count filter applied
> during download. Vas-0 RNA-seq data will be analysed in Module 04 but may
> yield fewer variant calls than the other accessions — students should keep
> this in mind when interpreting cross-accession comparisons.

> **RNA-seq data files are pre-staged at:**
> `/fs/scratch/PAS3260/Variant_Calling/Data/RNA-seq/`

---

## Part 6 — Pre-Tutorial Discussion Questions

Consider the following before starting Module 01. These questions do not have
single correct answers — they are designed to surface assumptions and get you
thinking about the experimental design decisions that shape every variant
calling project.

1. All five accessions used in this tutorial were sequenced on the same
   instrument (HiSeq 2000 for WGS, HiSeq 2500 for RNA-seq) as part of
   coordinated consortium projects. Why does instrument uniformity matter
   for multi-sample variant calling? What problems arise when samples from
   different instruments or sequencing centres are combined in the same
   GATK cohort?

2. The GBS data in this tutorial is simulated from WGS reads by extracting
   reads that overlap ApeKI cut sites. A student argues this is not a fair
   comparison because "real GBS reads always start at the restriction site."
   Is this a valid criticism? What specific differences between simulated and
   real GBS data could affect variant calling outcomes?

3. *Arabidopsis thaliana* is nearly completely homozygous due to
   self-fertilisation. If you call variants on the five accessions used here
   and observe that 3% of sites are called as heterozygous in Mah-6, what
   are the most likely biological and technical explanations? How would you
   distinguish between them?

4. The 1001 Genomes Project found ~10.7 million SNPs across 1,135 accessions.
   If you were designing a 200,000-SNP genotyping array to capture maximum
   diversity across the species, what criteria would you use to select which
   200,000 SNPs to include? How would your criteria differ if the array were
   designed specifically for a breeding program in northern European germplasm?

5. RNA-seq variant calling is sometimes described as "opportunistic" — variants
   are a by-product of an experiment designed for expression analysis. What
   systematic biases would you expect in a VCF produced exclusively from
   RNA-seq data, and which classes of functional variants would be systematically
   under-represented?

---

*[Module 01: Environment Setup & Reference Preparation →](01_environment_setup.md)*
