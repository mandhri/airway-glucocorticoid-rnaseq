# Snakemake RNA-seq Pipeline

A modular Snakemake workflow for RNA-seq differential expression analysis using DESeq2 and the Bioconductor airway dataset.

## Overview

This pipeline demonstrates a complete RNA-seq analysis workflow:
- Data import from Bioconductor packages
- DESeq2 differential expression analysis with proper design matrices
- Log fold change shrinkage using apeglm
- Gene set enrichment analysis (GSEA) with MSigDB Hallmark gene sets
- Modular rule organization for maintainability

## Quick Start

```bash

# Create and activate environment
mamba create -n snakemake-7 -c bioconda -c conda-forge snakemake=7
mamba activate snakemake-7

# Run pipeline
snakemake -p -j1 --use-conda --conda-frontend mamba

```