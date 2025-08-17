#Run DESeq2: fit model, normalised counts + shrunken results (simple).
rule deseq:
    input:
        dds = "results/deseq2/dds_raw.rds"
    output:
        dds_out = "results/deseq2/dds.rds",
        norm    = "results/deseq2/norm_counts.tsv",
        res     = "results/deseq2/results.tsv"
    params:
        alpha = 0.05
    conda: "../envs/rna.yml"
    shell:
        "mkdir -p $(dirname {output.res}) && "
        "export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=; "
        "Rscript --vanilla scripts/run_deseq.R "
        "--dds_in {input.dds} --alpha {params.alpha} "
        "--norm_counts {output.norm} --res_out {output.res} --dds_out {output.dds_out}"
