#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

# -------- CLI --------
opt <- OptionParser() |>
  add_option("--res_tsv", type="character", help="DESeq2 results.tsv (shrunken LFCs)") |>
  add_option("--species", type="character", default="Homo sapiens") |>
  add_option("--category", type="character", default="H") |>
  add_option("--out_csv", type="character", help="Output CSV for GSEA table") |>
  add_option("--out_fig", type="character", help="Output PDF for top-20 dotplot") |>
  parse_args()

dir.create(dirname(opt$out_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$out_fig), showWarnings = FALSE, recursive = TRUE)

# -------- 1) Load DE results --------
# Expected columns: ENSEMBL, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
res <- readr::read_tsv(opt$res_tsv, show_col_types = FALSE)

# Coerce numeric columns (some may be missing depending on shrinkage)
num_cols <- intersect(c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"), names(res))
res[num_cols] <- lapply(res[num_cols], function(x) suppressWarnings(as.numeric(x)))

# -------- 2) Map ENSEMBL -> SYMBOL (for MSigDB) --------
res$ENSEMBL_clean <- sub("\\..*$", "", res$ENSEMBL)

sym_map <- AnnotationDbi::mapIds(org.Hs.eg.db,
  keys = res$ENSEMBL_clean, keytype = "ENSEMBL",
  column = "SYMBOL", multiVals = "first"
)

res2 <- res %>%
  dplyr::mutate(SYMBOL = unname(sym_map[ENSEMBL_clean])) %>%
  dplyr::filter(!is.na(SYMBOL))

# Order by |stat| if available; otherwise by |log2FC|
order_val <- if ("stat" %in% names(res2) && any(is.finite(res2$stat))) {
  abs(res2$stat)
} else {
  abs(res2$log2FoldChange)
}

res2 <- res2 %>%
  dplyr::mutate(.ord = order_val) %>%
  dplyr::arrange(dplyr::desc(.ord)) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# -------- 3) Build the ranking --------
ranks <- if ("stat" %in% names(res2) && any(is.finite(res2$stat))) {
  res2$stat
} else {
  res2$log2FoldChange
}
names(ranks) <- res2$SYMBOL
ranks <- ranks[is.finite(ranks)]
stopifnot(length(ranks) > 0)

# -------- 4) Get Hallmark gene sets --------
msig <- msigdbr(species = opt$species, category = opt$category)
pathways <- split(msig$gene_symbol, msig$gs_name)

# -------- 5) Run GSEA --------
fg <- fgsea(pathways = pathways, stats = ranks, nperm = 10000)
fg <- arrange(fg, padj)

# Save table
readr::write_csv(fg, opt$out_csv)

# -------- 6) Plot top 20 pathways --------
topn <- head(fg, 20)
p <- ggplot(topn, aes(x = reorder(pathway, NES), y = NES, size = size, color = -log10(padj))) +
  geom_point() + coord_flip() +
  labs(x = "Pathway (Hallmarks)", y = "NES", size = "Genes", color = "-log10(padj)",
       title = "GSEA – Hallmarks (top 20)") +
  theme_minimal()

ggsave(opt$out_fig, p, width = 7, height = 6)
