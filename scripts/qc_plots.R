#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
  library(readr)
  library(dplyr)
  library(EnhancedVolcano)
})

# ---- CLI ----
opt <- OptionParser() |>
  add_option("--dds_in",  type="character") |>
  add_option("--res_in",  type="character") |>
  add_option("--topvar",  type="integer",   default=500) |>
  add_option("--figdir",  type="character") |>
  add_option("--pca",     type="character") |>
  add_option("--heat",    type="character") |>
  add_option("--ma",      type="character") |>
  add_option("--volc",    type="character") |>
  parse_args()

dir.create(opt$figdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load inputs
# fitted DESeqDataSet (~ cell + dex)
dds <- readRDS(opt$dds_in)          
res <- readr::read_tsv(opt$res_in, show_col_types = FALSE)

# ---- PCA on VST
# stabilise variance for distances/PCA
vsd <- vst(dds, blind = FALSE)  
p <- plotPCA(vsd, intgroup = c("dex","cell")) + ggtitle("PCA (VST)")
ggsave(opt$pca, p, width = 6, height = 5)

# ---- Sample-distance heatmap (top variable genes)
mat <- assay(vsd)
sel <- head(order(rowVars(mat), decreasing = TRUE), opt$topvar)
d   <- dist(t(mat[sel, ]))
pheatmap(as.matrix(d),
         annotation_col = as.data.frame(colData(vsd)[, c("dex","cell")]),
         filename = opt$heat)

# ---- MA plot (ggplot version using the results table)
# x-axis as log10(baseMean+1), y-axis as log2FC
p_ma <- ggplot(res, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "log10(baseMean + 1)", y = "log2FC", title = "MA (shrunken LFC)") +
  theme_minimal()
ggsave(opt$ma, p_ma, width = 6, height = 5)

# ---- Volcano (ENSEMBL labels; we’ll add SYMBOLs later) ----
res$LABEL <- res$ENSEMBL
pdf(opt$volc, width = 7, height = 6)
EnhancedVolcano(res,
  lab = res$LABEL,
  x = "log2FoldChange",
  y = "padj",
  title = "Dex (trt) vs Untreated",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  labSize = 3.5
)
dev.off()
