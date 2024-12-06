# Bulk RNA-seq Analysis Workflow

This repository contains scripts and instructions for the bulk RNA-seq analysis workflow used in "Atezolizumab following definitive chemoradiotherapy in patients with unresectable locally advanced esophageal squamous cell carcinoma – a multicenter phase 2 trial (EPOC1802)" . 
This repository contains scripts for processing and analyzing bulk RNA sequencing data, including quality control, alignment, expression quantification, differential expression analysis, and enrichment analysis (GSEA and ssGSEA). 



## Overview
This pipeline was developed for analyzing paired tumor RNA-seq samples. It includes:
1. **Trimming and Quality Control**: Removal of adapter sequences and low-quality bases.
2. **rRNA Filtering**: Removal of rRNA-mapped reads using Bowtie2.
3. **Alignment**: Mapping reads to the human genome (GRCh38) using STAR.
4. **Gene Expression Quantification**: Calculating counts and TPM using StringTie.
5. **Differential Expression Analysis**: Conducting DESeq2 analysis for gene expression comparisons.
6. **GSEA and ssGSEA**: Enrichment analysis using ClusterProfiler and GSVA.

---

## Directory Structure

```plaintext
├── codes/                  # Analysis scripts
│   ├── trimming_fastq.sh     # Step 1: Adapter trimming
│   ├── filter_rRNA.sh        # Step 2: rRNA filtering
│   ├── STAR_Stringtie.sh     # Step 3: Alignment and expression quantification
│   ├── DESeq2.R              # Step 4: Differential expression analysis
│   ├── volcanoplot.R         # Step 5: Volcano plot generation
│   ├── GSEA.R                # Step 6: GSEA
│   └── ssGSEA.R              # Step 6: ssGSEA
└── README.md                 # This file
```
## Pipeline Steps

Step 1: Quality Control and Adapter Trimming

Run trimming_fastq.sh to remove adapters and low-quality bases:

```bash
bash codes/trimming_fastq.sh
```

Step 2: Filtering rRNA Reads

Run filter_rRNA.sh to filter out reads mapped to rRNA:

```bash
bash codes/filter_rRNA.sh
```

Step 3: Read Alignment

Run STAR_Stringtie.sh for genome alignment and expression quantification:

```bash
bash codes/STAR_Stringtie.sh
```

Step 4: Differential Expression Analysis

Run DESeq2.R to identify differentially expressed genes:

```
Rscript codes/DESeq2.R
```

Step 5: Visualization

Generate volcano plots for DESeq2 results:

```
Rscript codes/volcanoplot.R
```

Step 6: Gene Set Enrichment Analysis

Perform GSEA:

```
Rscript codes/GSEA.R
```

Perform ssGSEA:

```
Rscript codes/ssGSEA.R
```


