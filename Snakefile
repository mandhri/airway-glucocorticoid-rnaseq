#Snakemake

#configfile: "config.yaml"

include: "rules/00_export.smk"
include: "rules/01_make_dds.smk"
include: "rules/02_deseq.smk"
include: "rules/03_qc.smk"

rule all:
    input:
        "results/deseq2/dds_raw.rds",
        "results/deseq2/results.tsv",
        "results/deseq2/norm_counts.tsv",
        "results/deseq2/dds.rds",
        "figures/pca.pdf",
        "figures/sample_distance_heatmap.pdf",
        "figures/MA.pdf",
        "figures/volcano.pdf"




