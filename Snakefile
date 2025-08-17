#Snakemake

#configfile: "config.yaml"

include: "rules/00_export.smk"
include: "rules/01_make_dds.smk"

rule all:
    input:
        "results/deseq2/dds_raw.rds"




