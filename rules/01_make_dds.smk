rule make_dds:
    input:
        counts  = "data/raw/airway_counts.tsv",
        coldata = "data/raw/airway_coldata.tsv"
    output:
        dds = "results/deseq2/dds_raw.rds"
    params:
        condition = "dex",
        reference = "untrt",
        design    = "~ cell + dex"
    conda: "../envs/rna.yml"
    shell:
        "mkdir -p $(dirname {output.dds}) && "
        "export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=; "
        "Rscript --vanilla scripts/make_dds.R "
        "--counts {input.counts} --coldata {input.coldata} "
        "--condition {params.condition} --reference {params.reference} "
        "--design '{params.design}' --dds_out {output.dds}"
