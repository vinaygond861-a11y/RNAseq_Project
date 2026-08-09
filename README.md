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

```
### Analysis Pipeline

The complete RNA-Seq analysis was performed through the following major steps:

1. **Raw RNA-Seq Data**
   - Paired-end sequencing data were obtained from the NCBI Sequence Read Archive (SRA).

2. **Quality Control**
   - FastQC was used to assess the quality of raw sequencing reads.

3. **Read Trimming**
   - fastp was used for adapter removal and quality filtering.

4. **Reference Genome Alignment**
   - Cleaned reads were aligned to the *E. coli* reference genome using HISAT2.

5. **BAM Processing**
   - SAMtools was used for processing and handling alignment files.

6. **Gene Quantification**
   - featureCounts was used to generate gene-level read count files.

7. **Count Matrix Construction**
   - Count files from all four samples were combined into a single count matrix.

8. **Differential Expression Analysis**
   - DESeq2 was used in R to identify genes with significant expression changes between Water and Ciprofloxacin conditions.

9. **Visualization**
   - MA plots and volcano plots were generated to visualize differential expression results.

---

# Results Summary

The final count matrix contained:

- **4,308 genes**
- **4 RNA-Seq samples**
- **2 Water control samples**
- **2 Ciprofloxacin-treated samples**

After DESeq2 analysis:

- **4,254 genes** had non-zero total read counts.
- **10 genes** were significantly differentially expressed.
- Significance threshold: **adjusted p-value < 0.05**

---

# Significant Differentially Expressed Genes

| Gene | Log2 Fold Change | Adjusted P-value |
|------|-----------------:|-----------------:|
| b0331 | 5.426 | 0.0111 |
| b0333 | 6.863 | 0.000739 |
| b0334 | 8.120 | 0.000136 |
| b0335 | 7.867 | 0.000136 |
| b0953 | 5.030 | 0.0215 |
| b1492 | 4.888 | 0.0215 |
| b1493 | 4.416 | 0.0463 |
| b1597 | 8.443 | 0.000064 |
| b4601 | 5.270 | 0.0158 |
| b2982 | 5.031 | 0.0215 |

The positive log2 fold-change values indicate higher expression in the **Water condition relative to Ciprofloxacin** for the current DESeq2 contrast.

---

# Differential Expression Visualization

## MA Plot

The MA plot was generated from the DESeq2 results to visualize the relationship between average normalized gene expression and log2 fold change.

It helps identify genes showing large expression changes between the two experimental conditions.
![MA Plot](plot/MA_plot.png
## Volcano Plot

The volcano plot combines:

- Log2 fold change
- Statistical significance

Genes with larger absolute fold changes and lower adjusted p-values represent stronger differential expression signals.

---

# Project Structure

```text
RNAseq_Project/
│
├── README.md
├── .gitignore
├── main.nf
├── main_old.nf
├── nextflow.config
├── SRR_Acc_List.txt
│
├── reference/
│
├── fastq/
│
├── results/
│   ├── counts/
│   ├── fastqc/
│   └── fastp/
│
└── workflow/
