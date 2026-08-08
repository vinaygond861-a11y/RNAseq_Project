# E. coli RNA-Seq Analysis: Ciprofloxacin Response

## Project Overview

This project presents an end-to-end RNA-Seq analysis of *Escherichia coli* to investigate transcriptional changes associated with Ciprofloxacin treatment compared with a Water control.

The analysis starts from raw paired-end RNA-Seq data and follows a reproducible bioinformatics workflow including quality control, read trimming, reference genome alignment, gene-level read quantification, and differential gene expression analysis using DESeq2.

The NGS processing pipeline was implemented using Nextflow, while downstream statistical analysis and visualization were performed in R.

---

## Aim

The main aim of this project is to identify genes that show significant changes in expression between Water-treated and Ciprofloxacin-treated *E. coli* samples.

### Research Question

Which genes are differentially expressed in *E. coli* in response to Ciprofloxacin treatment?

---

## Experimental Design

The dataset contains four RNA-Seq samples:

| Condition | Sample | SRA Accession |
|-----------|--------|---------------|
| Water | Water_1 | SRR22578513 |
| Water | Water_2 | SRR22578514 |
| Ciprofloxacin | Cipro_1 | SRR22578533 |
| Ciprofloxacin | Cipro_2 | SRR22578536 |

### Comparison

**Water vs Ciprofloxacin**

The experiment consists of:

- 2 Water control replicates
- 2 Ciprofloxacin treatment replicates
- Paired-end RNA-Seq data

---

# RNA-Seq Analysis Workflow

```text
Raw Paired-End FASTQ
        |
        v
     FastQC
        |
        v
      fastp
  Quality Control
   and Trimming
        |
        v
      HISAT2
Read Alignment to
Reference Genome
        |
        v
    SAMtools
BAM Processing
        |
        v
  featureCounts
Gene-Level Read
   Quantification
        |
        v
   Count Matrix
        |
        v
     DESeq2
Differential Gene
Expression Analysis
        |
        v
 Significant Genes
        |
        v
Visualization and
Biological Interpretation
