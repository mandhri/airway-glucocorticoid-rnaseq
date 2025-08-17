# GSEA on Hallmark pathways using shrunken DE results
rule gsea:
    input:
        res = "results/deseq2/results.tsv"
    output:
        table  = "results/gsea/hallmark_gsea.csv",
        dotfig = "figures/gsea_top20.pdf"
    params:
        species  = "Homo sapiens",
        category = "H"                    # MSigDB Hallmarks
    conda: "../envs/rna.yml"
    shell:
        "mkdir -p results/gsea figures && "
        "export R_LIBS_USER= R_LIBS= R_PROFILE_USER= R_ENVIRON_USER=; "
        "Rscript --vanilla scripts/gsea.R "
        "--res_tsv {input.res} "
        "--species \"{params.species}\" "
        "--category {params.category} "
        "--out_csv {output.table} "
        "--out_fig {output.dotfig}"
