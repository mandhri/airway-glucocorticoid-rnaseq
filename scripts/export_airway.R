#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)              
  library(airway)                
  library(readr)
  library(SummarizedExperiment)
  library(tibble)
})

# 1) Define the command-line options we expect.

opt <- OptionParser() |>
  add_option("--counts",  type="character") |>
  add_option("--coldata", type="character") |>
  parse_args()

# 2) Load the dataset from the package into memory.
data("airway")
se <- airway

# 3) Extract the counts matrix: genes (rows) x samples (columns).

# assay() pulls the counts slot
cts <- as.data.frame(assay(se, "counts"))
# make gene IDs a real column
cts <- tibble::rownames_to_column(cts, "ENSEMBL")
# save to the path passed as --counts
readr::write_tsv(cts, opt$counts)

# 4) Extract the sample metadata (one row per sample).
cd  <- as.data.frame(colData(se))
cd  <- tibble::rownames_to_column(cd, "sample")
readr::write_tsv(cd, opt$coldata)
