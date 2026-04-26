# From Reads to Chips: Variant Calling Across Sequencing Technologies

**HCS 7004 — Genome Analytics | The Ohio State University**

---

## Overview

Genetic variants — single nucleotide polymorphisms (SNPs), insertions, deletions,
and structural rearrangements — are the raw material of plant breeding, population
genetics, and genomic selection. Every SNP chip, genotyping-by-sequencing (GBS)
assay, and amplicon panel begins with the same question: *how do you reliably
distinguish a true genetic difference from a sequencing error?*

This tutorial series walks through the complete variant discovery pipeline for
three complementary sequencing strategies using *Arabidopsis thaliana* as the
reference organism. All three strategies produce a **Variant Call Format (VCF)**
file, but they differ in genome coverage, read depth, and the biological questions
they are best suited to answer.

| Track | Data Type | Library type | Primary Use Case |
|---|---|---|---|
| **A** | Whole Genome Sequencing (WGS) | Paired-end 101 bp | High-confidence SNP discovery, structural variants |
| **B** | Genotyping-by-Sequencing (GBS) | Single-end 101 bp (ApeKI simulated) | Population structure, diversity panels, QTL mapping |
| **C** | RNA-seq | Single-end 100 bp | Allele-specific expression, eQTLs, expressed variants |

By the end of this tutorial you will understand not just *how* to call variants,
but *why* the pipeline works at each step — and how the SNPs you discover end up
in the genotyping tools used in real breeding programs.

---

## Why *Arabidopsis thaliana*?

- **Compact, well-annotated genome:** ~135 Mb across 5 chromosomes (TAIR10).
  Jobs finish in minutes to hours, not days.
- **All three data types exist for the same five accessions.** WGS and RNA-seq
  data for Vas-0, Bez-9, Gen-8, Mah-6, and Usa-0 come from two companion studies
  published simultaneously in *Cell* (July 2016), making direct cross-track
  comparison of variant calls meaningful.
- **SNP chip connection:** Affymetrix and Illumina have both produced *Arabidopsis*
  genotyping arrays designed from 1001 Genomes variants, making the bridge from
  discovery to application completely traceable.
- **Population biology is rich and documented.** Students can place their variants
  in the context of global geographic structure, admixture, and local adaptation
  using a reference panel of >1,100 accessions.

---

## Module Structure

| Module | File | Topics Covered |
|---|---|---|
| [00](modules/00_conceptual_framework.md) | Conceptual Framework | Biology of variants; VCF format; overview of all three sequencing strategies; study organism and datasets |
| [01](modules/01_environment_setup.md) | Environment Setup & Reference Preparation | OSC setup; reference indexing; container acquisition; raw data staging; QC and trimming |
| [02](modules/02_wgs_variant_calling.md) | WGS Variant Calling | BWA-MEM2 alignment; GATK4 GVCF pipeline; joint genotyping; hard filtering |
| [03](modules/03_gbs_variant_calling.md) | GBS Variant Calling | ApeKI simulation rationale; alignment-based calling; population-level filtering |
| [04](modules/04_rnaseq_variant_calling.md) | RNA-seq Variant Calling | STAR 2-pass alignment; GATK RNA-seq best practices; allele-specific expression |
| [05](modules/05_variant_annotation.md) | Variant Annotation & Filtering | SnpEff functional annotation; VCFtools population filtering; LD basics with PLINK |
| [06](modules/06_snp_applications.md) | From SNPs to Genotyping Technologies | SNP chip probe design; AmpSeq arrays; KASP assays; the 1001 Genomes → chip pipeline |

---

## Compute Environment

This tutorial runs on the **Ohio Supercomputer Center (OSC)** under allocation
**PAS3260**.

- All bioinformatics tools run inside **Apptainer containers** — no `module load`
  conflicts between tools.
- R-based analyses in Modules 05–06 use the **native OSC R installation**:
  `module load gcc/12.3.0 R/4.5.2`
- Each student works in their own scratch directory:
  `/fs/scratch/PAS3260/<your_username>/variant_calling/`
- **Shared read-only data** (reference genome, raw FASTQs, GBS simulation files)
  are pre-staged at: `/fs/scratch/PAS3260/Variant_Calling/Data/`
- **Containers** are pre-pulled at:
  `/fs/scratch/PAS3260/Variant_Calling/Containers/`

### Shared directory structure

```
/fs/scratch/PAS3260/Variant_Calling/          (32 GB total, 43 files)
├── Containers/                               (4.0 GB — 16 .sif files)
└── Data/                                     (28 GB)
    ├── Reference/                            (223 MB)
    │   ├── fasta/
    │   │   └── Athaliana_TAIR10.fasta        (116 MB)
    │   └── gff/
    │       └── Athaliana_TAIR10.gff3         (107 MB)
    ├── WGS/                                  (14 GB — PE 101 bp, HiSeq 2000)
    │   ├── Vas-0_R1.fastq.gz / _R2.fastq.gz (1.3 GB each)
    │   ├── Bez-9_R1.fastq.gz / _R2.fastq.gz (1.6 GB each)
    │   ├── Gen-8_R1.fastq.gz / _R2.fastq.gz (950 MB / 930 MB)
    │   ├── Mah-6_R1.fastq.gz / _R2.fastq.gz (1.0 GB / 1.1 GB)
    │   └── Usa-0_R1.fastq.gz / _R2.fastq.gz (2.2 GB each)
    ├── GBS/                                  (1.4 GB — SE 101 bp, ApeKI simulated)
    │   ├── Vas-0_GBS.fastq.gz               (275 MB)
    │   ├── Bez-9_GBS.fastq.gz               (232 MB)
    │   ├── Gen-8_GBS.fastq.gz               (152 MB)
    │   ├── Mah-6_GBS.fastq.gz               (222 MB)
    │   ├── Usa-0_GBS.fastq.gz               (496 MB)
    │   ├── ApeKI_sites.bed                  (3.1 MB — cut site positions)
    │   └── ApeKI_GBS_windows.bed            (1.2 MB — ±200 bp windows, merged)
    └── RNA-seq/                              (13 GB — SE 100 bp, HiSeq 2500)
        ├── Vas-0_rnaseq.fastq.gz            (763 MB)
        ├── Bez-9_rnaseq.fastq.gz            (2.6 GB)
        ├── Gen-8_rnaseq.fastq.gz            (3.2 GB)
        ├── Mah-6_rnaseq.fastq.gz            (2.9 GB)
        └── Usa-0_rnaseq.fastq.gz            (3.1 GB)
```

### Container inventory

| Container file | Tool | Version | Used in modules |
|---|---|---|---|
| `sratools_3.2.1.sif` | SRA Toolkit | 3.2.1 | 01 |
| `ncbi_datasets_18.24.sif` | NCBI Datasets CLI | 18.24.0 | 01 |
| `fastp_1.3.2.sif` | fastp | 1.3.2 | 01 |
| `fastqc_0.12.1.sif` | FastQC | 0.12.1 | 01 |
| `multiqc_1.34.sif` | MultiQC | 1.34 | 01 |
| `bwamem2_2.2.3.sif` | BWA-MEM2 | 2.2.3 | 01, 02, 03 |
| `samtools_1.23.1.sif` | SAMtools | 1.23.1 | 01, 02, 03, 04 |
| `star_2.7.11b.sif` | STAR | 2.7.11b | 01, 04 |
| `gatk4_4.6.2.0.sif` | GATK4 | 4.6.2.0 | 02, 03, 04 |
| `picard_3.4.0.sif` | Picard | 3.4.0 | 02, 04 |
| `stacks_2.68.sif` | Stacks | 2.68 | 03 |
| `bcftools_1.23.1.sif` | BCFtools | 1.23.1 | 02, 03, 04, 05 |
| `vcftools_0.1.17.sif` | VCFtools | 0.1.17 | 05, 06 |
| `snpeff_5.4.0c.sif` | SnpEff | 5.4.0c | 05 |
| `plink_1.90b7.sif` | PLINK | 1.90b7 | 05, 06 |
| `bedtools_2.31.1.sif` | BEDTools | 2.31.1 | 03, 06 |

> **R is not containerised.** Load the native OSC installation wherever R is
> needed: `module load gcc/12.3.0 R/4.5.2`

---

## Study Accessions

All five accessions are drawn from two companion studies published simultaneously
in *Cell* (14 July 2016). Because WGS and RNA-seq data originate from the same
physical plant material and experimental project, direct comparison of variants
across sequencing strategies is biologically meaningful.

| Accession | 1001G ID | WGS SRR | RNA-seq SRX | WGS instrument | WGS coverage | Origin |
|---|---|---|---|---|---|---|
| Vas-0 | 9902 | SRR1946455 | SRX1735056 | HiSeq 2000 | ~32× | Västerås, Sweden |
| Bez-9 | 9928 | SRR1946475 | SRX1735072 | HiSeq 2000 | ~25× | Besançon, France |
| Gen-8 | 9909 | SRR1946461 | SRX1735063 | HiSeq 2000 | ~14× | Geneva, Switzerland |
| Mah-6 | 9906 | SRR1946459 | SRX1735060 | HiSeq 2000 | ~17× | Mahón, Menorca, Spain |
| Usa-0 | 9929 | SRR1946457 | SRX1735058 | HiSeq 2000 | ~18× | Utsunomiya, Japan |

The five accessions span a broad geographic range from the Iberian Peninsula to
East Asia, representing both the central European diversity core of *Arabidopsis*
and peripheral relict populations. The coverage range (~14–32×) is natural
variation arising from library complexity differences across accessions and is
itself a teaching point: Module 02 demonstrates how depth affects genotype
confidence scores at heterozygous-looking sites in this naturally inbred species.

> **GBS data** for all five accessions is derived by in silico ApeKI digestion of
> the WGS reads. See Module 03 for the full biological and methodological
> rationale.

---

## Conventions Used in This Tutorial

- `code blocks` — commands to type exactly as written
- `<your_username>` — replace with your OSC username
- **Bold text** — key concept introduced for the first time
- > Blockquotes — biological or conceptual take-home messages
- ⚠️ — common pitfall worth pausing on before proceeding

---

## References

1. The 1001 Genomes Consortium (2016). 1,135 Genomes Reveal the Global Pattern
   of Polymorphism in *Arabidopsis thaliana*. *Cell* 166:481–491.
   https://doi.org/10.1016/j.cell.2016.05.063

2. Kawakatsu, T. et al. (2016). Epigenomic Diversity in a Global Collection of
   *Arabidopsis thaliana* Accessions. *Cell* 166:492–505.
   https://doi.org/10.1016/j.cell.2016.06.044

3. McKenna, A. et al. (2010). The Genome Analysis Toolkit: A MapReduce framework
   for analyzing next-generation DNA sequencing data. *Genome Research* 20:1297–1303.
   https://doi.org/10.1101/gr.107524.110

4. Elshire, R.J. et al. (2011). A Robust, Simple Genotyping-by-Sequencing (GBS)
   Approach for High Diversity Species. *PLOS ONE* 6:e19379.
   https://doi.org/10.1371/journal.pone.0019379

5. Van der Auwera, G.A. & O'Connor, B.D. (2020). *Genomics in the Cloud*.
   O'Reilly Media. (GATK Best Practices reference)

6. Dobin, A. et al. (2013). STAR: ultrafast universal RNA-seq aligner.
   *Bioinformatics* 29:15–21. https://doi.org/10.1093/bioinformatics/bts635
