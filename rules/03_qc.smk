# QC plots: PCA, sample-distance heatmap, MA, Volcano
rule qc:
    input:
        dds = "results/deseq2/dds.rds",
        res = "results/deseq2/results.tsv"
    output:
        pca     = "figures/pca.pdf",
        heatmap = "figures/sample_distance_heatmap.pdf",
        ma      = "figures/MA.pdf",
        volc    = "figures/volcano.pdf"
    params:
        topvar = 500
    conda: "../envs/rna.yml"
    shell:
        "mkdir -p figures && "
        "export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=; "
        "Rscript --vanilla scripts/qc_plots.R "
        "--dds_in {input.dds} --res_in {input.res} "
        "--topvar {params.topvar} --figdir figures "
        "--pca {output.pca} --heat {output.heatmap} --ma {output.ma} --volc {output.volc}"
