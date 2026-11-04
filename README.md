# Cross-species foregut single-cell RNA-seq analysis

This repository contains an R script used to process and compare single-cell RNA-sequencing (scRNA-seq) data from **human, mouse, and chick foregut tissues**. It includes quality control of a human count matrix, one-to-one ortholog mapping, cross-species integration, cluster annotation, differential expression analysis, cell-cycle regression, and visualization.

## Repository contents

| File | Description |
| --- | --- |
| `HMS_Code_Repo_critical_fixes.R` | Analysis script. Critical code errors identified during a pre-publication review were corrected; other analytical choices and historical workflow sections were retained. |
| `README.md` | This documentation. |

**Scope and reproducibility:** This is an analysis record, **not a turnkey, independently tested pipeline**. The script contains earlier and later analysis variants, requires input files that are **not included** here, and uses computer-specific paths. The corrected script has not been executed end-to-end against the original datasets as part of this repository preparation.

## What the script does

In approximate order, the code:

1. Imports human raw counts and cell/gene annotation tables from **GSE157329**, constructs a Seurat object, and filters cells using the thresholds specified in the script.
2. Retrieves mouse-to-human and chick-to-human **one-to-one orthologs** using `biomaRt`, converts gene symbols, and retains genes shared by all three species.
3. Subsets annotated human foregut-related cells and integrates human, mouse, and chick Seurat objects using **Seurat v5-style layers and CCA integration**.
4. Clusters cells, examines marker-gene expression, assigns cell-type identities, and performs species-level differential expression comparisons within epithelial and mesenchymal populations.
5. Runs a later version of the integration with **S.Score and G2M.Score regression**, then saves the annotated integrated object and generates UMAP figures.
6. Subsets selected mesenchymal populations for re-clustering and performs additional filtering of specified sample/stage metadata.

The script retains an **earlier mouse input/ortholog mapping workflow**, followed by a later workflow using a different mouse object (`seurat_phase_Ms.rds`). The later cell-cycle-regressed analysis produces `seurat_HuChMs_regressed_Ruifinal.rds`. These branches should not be interpreted as a single fully automated workflow without reviewing their intended sequence.

## Software requirements

- **R**, with a Seurat release supporting `IntegrateLayers()`, `JoinLayers()`, and assay `layer` access (i.e., a **Seurat v5-era workflow**).
- Packages loaded directly in the script:
  - CRAN/ecosystem: `Seurat`, `ggplot2`, `dplyr`, `RColorBrewer`, `Matrix`, `data.table`, `DoubletFinder`.
  - Bioconductor/ecosystem: `DropletUtils`, `SingleR`, `celldex`, `glmGamPoi`, `scDblFinder`, `biomaRt`.
- Internet access for the live **Ensembl BioMart** ortholog lookups.

The exact R and package versions of the original analysis were **not provided**, and no lockfile or `sessionInfo()` snapshot is included. Some loaded libraries are not used explicitly in every section.

## Required input data

All file paths below are **as written in the script**. The leading `~/amdata/` refers to the original analyst's local directory and must be adjusted for another environment.

### Starting inputs

| Input file (relative to `~/amdata/`) | Role |
| --- | --- |
| `raw_data/GSE157329_raw_counts_human.mtx.gz` | Human gene-by-cell count matrix |
| `raw_data/GSE157329_cell_annotate_human.txt.gz` | Human cell annotation table, including `cell_id` |
| `raw_data/GSE157329_gene_annotate_human.txt.gz` | Human gene annotation table, including `gene_short_name` |
| `chick_seuratfiltered.rds` | Previously processed chick Seurat object |
| `seurat_phase_Ms.rds` | Mouse Seurat object used in the later ortholog-conversion section |
| `mousesub_seuratfiltered.rds` | Mouse object used in the **earlier, superseded** ortholog-conversion section |

The human input accession **GSE157329** is named in the script; accession numbers and provenance for the supplied mouse and chick Seurat inputs are not documented in the code and should be supplied separately if this repository accompanies a publication.

### Other file dependencies and manual handoffs

The script reads several `.rds` files that earlier sections write, including `hum_ortho.rds`, `chick_ortho.rds`, `*_orthocommon_order2.rds`, and `celltype_annotated.rds`.

Two additional handoffs need attention **before an end-to-end run**:

- The script reads `~/amdata/seurat_HuMsCh_finalRui.rds`. An earlier comment says `final_seurat2.rds` was **manually renamed** to this name; the script does not perform that rename.
- Near the end, it reads `~/amdata/For Rui/seurat_HuChMs_regressed_Ruifinal.rds`, although the earlier save command writes `~/amdata/seurat_HuChMs_regressed_Ruifinal.rds`. The file needs to be present at the read location, or the path needs to be reconciled.

Data files, intermediate Seurat objects, and final result files are **not distributed in this repository**.

## How to use the script

1. Obtain the input data and arrange the files in a working directory, preserving the metadata fields expected by the script.
2. Install the R packages listed above. Confirm you are using a compatible Seurat version.
3. Replace the local `~/amdata/` paths with paths on your machine. Also resolve the two manual file handoffs described above.
4. Open `HMS_Code_Repo_critical_fixes.R` in RStudio or another R editor. **Review and run the intended sections in sequence** rather than assuming the entire file can run unchanged from start to finish.
5. Check intermediate object dimensions, shared ortholog counts, cluster labels, and the generated figures against the original analysis outputs before using results for publication.

For an independently reproducible release, the analysis should additionally document the precise data accessions, original R/package versions, Ensembl release, and the connection between script sections and published figures.

## Outputs

The script saves multiple Seurat `.rds` objects and intermediate text tables under `~/amdata/`. Important named outputs include:

- `seurat_HuChMs_regressed_Ruifinal.rds` — later cell-cycle-regressed, integrated, annotated object.
- `cluster0_1subset.rds` — re-clustered selected mesenchymal populations.
- `celltypeDEGs.txt`, `top10markers.txt`, `top50markers.txt` — marker tables.
- `epi_*.txt` and `mes_*.txt` — species-level differential expression tables.
- UMAP figure PDFs, including `species.pdf`, `species_split.pdf`, `cluster.pdf`, `cluster_split.pdf`, `celltype.pdf`, and `celltype_split.pdf`.

**Important:** Some PDF filenames (for example `species.pdf`) occur in more than one analysis section and may be overwritten. PDF files are saved relative to the current R working directory, whereas most `.rds` and `.txt` files use `~/amdata/`.

## Analytical and interpretation notes

- Gene orthologs are fetched from the **live** Ensembl BioMart service; results may change when its underlying database changes. The original Ensembl release is not pinned in this script.
- Cell-type annotation uses ordered vectors associated with cluster levels; labels should be verified against the actual clusters when reproducing results.
- The `FindMarkers()` / `FindAllMarkers()` species comparisons operate on individual cells. Their significance tests should not automatically be interpreted as tests across independent biological replicates; the appropriate statistical unit depends on sample design.
- The script contains code from multiple dates and includes superseded processing steps. These were **not removed**, to preserve the provided analysis history.
- The critical fixes do not constitute independent verification of the biological results or a full computational reproducibility test.

## License and citation

No license or definitive publication citation was supplied with the code. A repository owner should add the intended license and associated study citation/DOI before advertising the repository as a reusable published resource.
