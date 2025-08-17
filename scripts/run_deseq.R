#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(DESeq2)
  library(apeglm)
})

opt <- OptionParser() |>
  add_option("--dds_in",        type="character") |>
  add_option("--alpha",         type="double",   default=0.05) |>
  add_option("--lfc_threshold", type="double",   default=1.0) |>
  add_option("--norm_counts",   type="character") |>
  add_option("--res_out",       type="character") |>
  add_option("--dds_out",       type="character") |>
  parse_args()

# Load raw dds
dds <- readRDS(opt$dds_in)

# light in-object filtering
keep <- rowSums(counts(dds) >= 10) >=2
dds  <- dds[keep, ]

# Fit DESeq2
dds <- DESeq(dds)

# Normalised counts
# Normalised counts (library-size normalisation)
nc <- counts(dds, normalized = TRUE)

# Ensure gene IDs (row names) exist — use the dds row names as the source of truth
if (is.null(rownames(nc))) {
  rownames(nc) <- rownames(dds)
}

# Build a data frame with ENSEMBL as the first column
nc_df <- as.data.frame(nc)
nc_df$ENSEMBL <- rownames(nc_df)
# put ENSEMBL first
nc_df <- nc_df[, c("ENSEMBL", setdiff(colnames(nc_df), "ENSEMBL"))]

readr::write_tsv(nc_df, opt$norm_counts)


# 3) Results for dex: treated vs untreated (airway levels: untrt, trt)

res <- results(dds, name = "dex_trt_vs_untrt", alpha = opt$alpha)

# LFC shrinkage with apeglm (requires coef, not contrast)
res <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm", res = res)

# 5) Save table and object
res_df <- as.data.frame(res)
res_df$ENSEMBL <- rownames(res_df)
res_df <- res_df[, c("ENSEMBL", setdiff(colnames(res_df), "ENSEMBL"))]
readr::write_tsv(res_df, opt$res_out)

saveRDS(dds, opt$dds_out)
