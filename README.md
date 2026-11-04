# Cross-species foregut single-cell RNA-seq analysis

This repository contains an R script used to process and compare single-cell RNA-sequencing (scRNA-seq) data from **human, mouse, and chick foregut tissues**. It includes quality control of a human count matrix, one-to-one ortholog mapping, cross-species integration, cluster annotation, differential expression analysis, cell-cycle regression, and visualization. Contributed by Albert Z. Ma.

## Repository contents

| File | Description |
| --- | --- |
| `Cross_Species_Integration.R` | Analysis script. |
| `README.md` | This documentation. |


## What the script does

In approximate order, the code:

1. Imports human raw counts and cell/gene annotation tables from **GSE157329**, constructs a Seurat object, and filters cells using the thresholds specified in the script.
2. Retrieves mouse-to-human and chick-to-human **one-to-one orthologs** using `biomaRt`, converts gene symbols, and retains genes shared by all three species.
3. Subsets annotated human foregut-related cells and integrates human, mouse **GSE353356**, and chick **GSE307539** Seurat objects using **Seurat v5-style layers and CCA integration**.
4. Clusters cells, examines marker-gene expression, assigns cell-type identities, and performs species-level differential expression comparisons within epithelial and mesenchymal populations.
5. Runs the integration with **S.Score and G2M.Score regression**, then saves the annotated integrated object and generates UMAP figures.
6. Subsets selected mesenchymal populations for re-clustering and performs additional filtering of specified sample/stage metadata.


## Software requirements

- **R**, with a Seurat release supporting `IntegrateLayers()`, `JoinLayers()`, and assay `layer` access (i.e., a **Seurat v5-era workflow**).
- Packages loaded directly in the script:
  - CRAN/ecosystem: `Seurat`, `ggplot2`, `dplyr`, `RColorBrewer`, `Matrix`, `data.table`, `DoubletFinder`.
  - Bioconductor/ecosystem: `DropletUtils`, `SingleR`, `celldex`, `glmGamPoi`, `scDblFinder`, `biomaRt`.
- Internet access for the live **Ensembl BioMart** ortholog lookups.



Citation:
```bash
@article{Yan2026ForegutComparison,
  author  = {Yan, Rui and Ma, Albert Z. and Tabin, Clifford J.},
  title   = {Human foregut morphogenesis exhibits murine-like molecular patterning but avian-like epithelial architecture},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.09.17.752534},
  url     = {https://www.biorxiv.org/content/10.64898/2026.09.17.752534v1}
}

```

## License and citation

No license or definitive publication citation was supplied with the code. A repository owner should add the intended license and associated study citation/DOI before advertising the repository as a reusable published resource.
