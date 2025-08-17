#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(DESeq2)
})

opt <- OptionParser() |>
  add_option("--counts",     type="character") |>
  add_option("--coldata",    type="character") |>
  add_option("--design",     type="character", default="~ cell + dex") |>
  add_option("--condition",  type="character", default="dex") |>
  add_option("--reference",  type="character", default="untrt") |>
  add_option("--min_count",  type="integer",   default=10) |>
  add_option("--min_samples",type="integer",   default=2) |>
  add_option("--dds_out",    type="character") |>
  add_option("--summary_out",type="character") |>
  parse_args()

# Read counts table (first column is gene id)
cts <- readr::read_tsv(opt$counts, show_col_types = FALSE)
#Quick sanity check: there must be at least 2 columns — one for gene IDs + at least one sample column.
#If not, it stops with an error so you do not continue with a broken table.
stopifnot(ncol(cts) > 1)
rownames(cts) <- cts[[1]]
cts[[1]] <- NULL

# --- Read coldata, align to counts columns ---
coldata <- readr::read_tsv(opt$coldata, show_col_types = FALSE)
stopifnot("sample" %in% colnames(coldata))
rownames(coldata) <- coldata$sample
coldata <- coldata[colnames(cts), , drop = FALSE]

# --- Factors & design ---
stopifnot(opt$condition %in% colnames(coldata))
coldata[[opt$condition]] <- relevel(factor(coldata[[opt$condition]]), ref = opt$reference)
if ("cell" %in% names(coldata)) coldata$cell <- factor(coldata$cell)
form <- as.formula(opt$design)

# --- Construct DESeqDataSet ---
dds <- DESeqDataSetFromMatrix(countData = as.matrix(cts),
                              colData   = coldata,
                              design    = form)

saveRDS(dds, opt$dds_out)
