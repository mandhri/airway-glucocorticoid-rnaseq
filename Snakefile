#Snakemake

configfile: "config.yaml"

rule all:
    input:
        "data/raw/airway_counts.tsv",
        "data/raw/airway_coldata.tsv"

rule export_airway:
    output:
        counts  = "data/raw/airway_counts.tsv",
        coldata = "data/raw/airway_coldata.tsv"
    conda: "envs/rna.yml"
    shell:
        "mkdir -p $(dirname {output.counts}) && "
        "export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=; "
        "Rscript --vanilla scripts/export_airway.R "
        "--counts {output.counts} --coldata {output.coldata}"



