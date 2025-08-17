rna-seq-snakemake/
├─ Snakefile
├─ config.yaml
├─ envs/
│  └─ rna.yml
├─ scripts/
│  ├─ export_airway.R
│  ├─ deseq2.R
│  ├─ qc_plots.R
│  └─ gsea.R
├─ data/
│  ├─ raw/          # created by the pipeline
│  └─ processed/    # created by the pipeline
├─ results/
│  ├─ deseq2/       # created by the pipeline
│  └─ gsea/         # created by the pipeline
└─ figures/         # created by the pipeline
