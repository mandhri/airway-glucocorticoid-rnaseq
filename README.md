# Snakemake RNA-seq Pipeline

A modular Snakemake workflow for RNA-seq differential expression analysis using DESeq2 and the Bioconductor airway dataset.

## Overview

This pipeline demonstrates a complete RNA-seq analysis workflow:
- Data import from Bioconductor packages
- DESeq2 differential expression analysis with proper design matrices
- Log fold change shrinkage using apeglm
- Gene set enrichment analysis (GSEA) with MSigDB Hallmark gene sets
- Modular rule organization for maintainability

## Dataset Information
The airway dataset contains RNA-seq data from:

Himes et al. (2014). "RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells.
" PLoS One. 9(6):e99625.

## Quick Start

```bash

# Create and activate environment
mamba create -n snakemake-7 -c bioconda -c conda-forge snakemake=7
mamba activate snakemake-7

# Run pipeline
snakemake -p -j1 --use-conda --conda-frontend mamba

```
## Quick Start

### System Requirements

Linux/macOS
Snakemake 7.x
Mamba/Conda package manager
4GB RAM minimum
~2GB disk space


## Environment Setup

The pipeline uses isolated conda environments specified in envs/rna.yml:

```yaml

channels:
  - conda-forge
  - bioconda
dependencies:
  - r-base=4.3.*
  - bioconductor-deseq2
  - bioconductor-apeglm
  - bioconductor-airway
  - msigdbr
  - fgsea

```

## R Packages (specified in envs/rna.yml)

- r-base 4.3.*
- bioconductor-deseq2
- bioconductor-apeglm
- bioconductor-airway
- bioconductor-summarizedexperiment
- msigdbr
- fgsea
- org.Hs.eg.db
- r-optparse, r-readr, r-tibble, r-dplyr, r-ggplot2

## Project Structure

```
.
├── Snakefile               # Main workflow orchestrator
├── config.yaml             # User-configurable parameters
├── envs/
│   └── rna.yml             # R/Bioconductor environment
├── rules/                  # Modular workflow components
│   ├── 00_export.smk       # Data export rules
│   ├── 01_make_dds.smk     # DESeq2 object creation
│   ├── 02_deseq.smk        # Differential expression
│   └── 03_gsea.smk         # Enrichment analysis
├── scripts/                # Analysis scripts
│   ├── export_airway.R    
│   ├── make_dds.R         
│   ├── run_deseq.R        
│   └── gsea.R             
└── results/                # Generated outputs
    ├── deseq2/
    └── gsea/
```

## Usage

### Basic Commands

``` bash

# Dry run (see what would be executed)
snakemake -n

# Run with single core
snakemake -p -j1 --use-conda --conda-frontend mamba

# Run with 4 cores
snakemake -p -j4 --use-conda --conda-frontend mamba

# Generate DAG visualization
snakemake --dag | dot -Tpng > workflow.png

# Run specific target
snakemake -p -j1 --use-conda results/deseq2/results.tsv

```


## Analysis Details
Parameters are currently hardcoded in the rule files:

DESeq2 Analysis (rules/01_make_dds.smk, rules/02_deseq.smk)

- Design formula: ~ cell + dex
- Reference level: untrt (untreated)
- Significance threshold (alpha): 0.05
- Minimum counts filter: ≥10 counts in ≥2 samples
- LFC shrinkage method: apeglm

GSEA Analysis (rules/03_gsea.smk)

- Species: Homo sapiens
- Gene sets: MSigDB Hallmark (category "H")
- Permutations: 10,000
- Gene ranking: Wald statistic from DESeq2

### Pipeline Steps

1. Data Export (rules/00_export.smk)

- Extracts count matrix and metadata from airway package
- Outputs: data/raw/airway_counts.tsv, data/raw/airway_coldata.tsv


2. Create DESeq2 Object (rules/01_make_dds.smk)

- Builds DESeqDataSet with design ~ cell + dex
- Output: results/deseq2/dds_raw.rds


3. Differential Expression (rules/02_deseq.smk)

- Filters low-count genes
- Runs DESeq2 normalization and testing
- Applies apeglm shrinkage
- Outputs: normalised counts, DE results, final dds object


Gene Set Enrichment (rules/03_gsea.smk)

- Maps Ensembl IDs to gene symbols
- Runs GSEA with Hallmark gene sets
- Outputs: enrichment table and visualisation

## Output Files


| File | Description |
|------|-------------|
| `data/raw/airway_counts.tsv` | Raw count matrix (64,102 genes × 8 samples) |
| `data/raw/airway_coldata.tsv` | Sample metadata |
| `results/deseq2/dds_raw.rds` | Initial DESeq2 object |
| `results/gsea/hallmark_gsea.csv` | GSEA results for 50 Hallmark gene sets |
| `figures/gsea_top20.pdf` | Top 20 enriched pathways visualization |













