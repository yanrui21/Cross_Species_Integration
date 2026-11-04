# 15 cores, # of hours, 128 total memory in GB


library(DropletUtils)
library(Seurat)
library(SingleR)
library(celldex)
library(ggplot2)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(DoubletFinder)
library(glmGamPoi)
library(dplyr)
library(scDblFinder)

library(Matrix)
library(data.table)

mat <- readMM("~/amdata/raw_data/GSE157329_raw_counts_human.mtx.gz")
cells <- fread("~/amdata/raw_data/GSE157329_cell_annotate_human.txt.gz")
genes <- fread("~/amdata/raw_data/GSE157329_gene_annotate_human.txt.gz")

mat <- as(mat, "CsparseMatrix")

dim(mat)
nrow(genes) #nrow(mat) should = nrow(genes); ncol(mat) should = nrow(cells)
nrow(cells)

rownames(mat) <- make.unique(genes$gene_short_name)
colnames(mat) <- cells$cell_id

meta <- as.data.frame(cells)
rownames(meta) <- meta$cell_id

seu <- CreateSeuratObject(counts = mat, meta.data = meta, project   = "GSE157329", min.cells = 3, min.features = 200)

#27044 features left, i.e. 5307 genes that were present in less than 3 cells were removed. No cells were removed

VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

saveRDS(seu, file = "~/amdata/GSE157329_seurat.rds")




################################Basic QC######################################


hum <- seu
hum <- PercentageFeatureSet(hum, pattern = "^MT-", col.name = "percent.mito")
hum <- PercentageFeatureSet(hum, "^RP[SL]", col.name = "percent.ribo")
hum <- PercentageFeatureSet(hum, "^HB[^(P)]", col.name = "percent.HB")
hum <- PercentageFeatureSet(hum, "PECAM1|PF4", col.name = "percent.plate")

VlnPlot(hum, features = c("percent.mito","percent.ribo","nFeature_RNA","nCount_RNA"),
        pt.size = 0, ncol = 4)

hum <- subset(hum, subset = nFeature_RNA < 4500)
hum <- subset(hum, subset = percent.mito < 10)
hum <- subset(hum, subset = percent.ribo > 5)

saveRDS(hum, file = "~/amdata/hum_seuratfiltered.rds")




#########################################Converting gene orthologs with biomaRt###########################8/4/26

library(dplyr)
library(biomaRt)

mouse_obj <- readRDS("~/amdata/mousesub_seuratfiltered.rds") ###*****Replaced with Rui's own Mouse data later
chick_obj <- readRDS("~/amdata/chick_seuratfiltered.rds") ###Should be identical to Rui's original chick object (39250 features across 17765 samples within 3 assays)



human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
 
ortho <- getBM(
  attributes = c(
    "external_gene_name",                           
    "mmusculus_homolog_associated_gene_name",       
    "mmusculus_homolog_orthology_type",
    "ggallus_homolog_associated_gene_name",     
    "ggallus_homolog_orthology_type"
  ),
  mart = human
)


make_map <- function(df, sym_col, type_col) {
  df %>%
    dplyr::rename(source = all_of(sym_col), otype = all_of(type_col)) %>%
    dplyr::filter(otype == "ortholog_one2one",
                  external_gene_name != "", source != "") %>%
    dplyr::distinct(source, external_gene_name) %>%
    dplyr::group_by(source)             %>% dplyr::filter(n() == 1) %>% dplyr::ungroup() %>%
    dplyr::group_by(external_gene_name) %>% dplyr::filter(n() == 1) %>% dplyr::ungroup()
}
 
mouse_map <- make_map(ortho, "mmusculus_homolog_associated_gene_name",
                             "mmusculus_homolog_orthology_type")
chick_map <- make_map(ortho, "ggallus_homolog_associated_gene_name",
                             "ggallus_homolog_orthology_type")
 

convert_to_human <- function(obj, map, assay = "RNA") {
  counts <- GetAssayData(obj, assay = assay, layer = "counts")  
 
  keep    <- rownames(counts) %in% map$source
  counts  <- counts[keep, , drop = FALSE]
 

  new_names <- map$external_gene_name[match(rownames(counts), map$source)]
  rownames(counts) <- new_names
 

  if (any(duplicated(rownames(counts)))) {
    counts <- rowsum(as.matrix(counts), group = rownames(counts))
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
 
  CreateSeuratObject(counts = counts, meta.data = obj@meta.data,
                     project = Project(obj))
}
 
mouse_h <- convert_to_human(mouse_obj, mouse_map)
chick_h <- convert_to_human(chick_obj, chick_map)

human_obj <- readRDS("~/amdata/hum_seuratfiltered.rds")
human_h <- human_obj   # already in human symbols


 
#i.e. these each contain all the genes that biomaRt could find an ortholog for in humans
saveRDS(human_h, file = "~/amdata/hum_ortho.rds")
saveRDS(mouse_h, file = "~/amdata/mouse_ortho.rds")
saveRDS(chick_h, file = "~/amdata/chick_ortho.rds")


common <- Reduce(intersect, list(rownames(human_h),
                                 rownames(mouse_h),
                                 rownames(chick_h)))

message(length(common), " shared orthologous genes across all three species") ###10497 genes are shared between all 3

#Rebuilt seurat objects so that the genes are in the same order for all of them; unlike above(obj_h), these (obj_3) have ONLY genes that biomaRt could find human orthologs for in BOTH mouse and chick, i.e. are present in all three species

reorder_features <- function(obj, genes) {
  cm <- GetAssayData(obj, layer = "counts")[genes, ]  
  CreateSeuratObject(counts = cm, meta.data = obj@meta.data, project = Project(obj))
}
mouse_3 <- reorder_features(mouse_h, common)
chick_3 <- reorder_features(chick_h, common)
human_3 <- reorder_features(human_h, common)


stopifnot(
  all(rownames(human_3) == common),
  all(rownames(mouse_3) == common),
  all(rownames(chick_3) == common)
)
identical(rownames(human_3), rownames(mouse_3))   # should be TRUE
identical(rownames(human_3), rownames(chick_3))   # should be TRUE. Use this piece of code to check if all the genes are in the right order



saveRDS(human_3, file = "~/amdata/hum_orthocommon_order.rds")
saveRDS(mouse_3, file = "~/amdata/mouse_orthocommon_order.rds")
saveRDS(chick_3, file = "~/amdata/chick_orthocommon_order.rds")



human_h <- NormalizeData(human_h)
mouse_h <- NormalizeData(mouse_h)
chick_h <- NormalizeData(chick_h)



####(8/16/26) REDO WITH RUI's MOUSE SEURAT OBJECT; IGNORE ANY OF THE MOUSE CODE BEFORE

library(dplyr)
library(biomaRt)

mouse_obj <- readRDS("~/amdata/seurat_phase_Ms.rds")


human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
 

ortho <- getBM(
  attributes = c(
    "external_gene_name",                            # human symbol
    "mmusculus_homolog_associated_gene_name",        # mouse ortholog
    "mmusculus_homolog_orthology_type"
  ),
  mart = human
)
 

make_map <- function(df, sym_col, type_col) {
  df %>%
    dplyr::rename(source = all_of(sym_col), otype = all_of(type_col)) %>%
    dplyr::filter(otype == "ortholog_one2one",
                  external_gene_name != "", source != "") %>%
    dplyr::distinct(source, external_gene_name) %>%
    dplyr::group_by(source)             %>% dplyr::filter(n() == 1) %>% dplyr::ungroup() %>%
    dplyr::group_by(external_gene_name) %>% dplyr::filter(n() == 1) %>% dplyr::ungroup()
}
 
mouse_map <- make_map(ortho, "mmusculus_homolog_associated_gene_name",
                             "mmusculus_homolog_orthology_type")

 
convert_to_human <- function(obj, map, assay = "RNA") {
  counts <- GetAssayData(obj, assay = assay, layer = "counts")  
 
  keep    <- rownames(counts) %in% map$source
  counts  <- counts[keep, , drop = FALSE]

  new_names <- map$external_gene_name[match(rownames(counts), map$source)]
  rownames(counts) <- new_names
 

  if (any(duplicated(rownames(counts)))) {
    counts <- rowsum(as.matrix(counts), group = rownames(counts))
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
 
  CreateSeuratObject(counts = counts, meta.data = obj@meta.data,
                     project = Project(obj))
}
 
mouse_h <- convert_to_human(mouse_obj, mouse_map)

saveRDS(mouse_h, file = "~/amdata/mouse_ortho2.rds")


human_h <- readRDS("~/amdata/hum_ortho.rds")
chick_h <- readRDS("~/amdata/chick_ortho.rds")


 
common <- Reduce(intersect, list(rownames(human_h),
                                 rownames(mouse_h),
                                 rownames(chick_h)))

message(length(common), " shared orthologous genes across all three species") ###10497 genes are shared between all 3


#Rebuilt seurat objects so that the genes are in the same order for all of them; unlike above(obj_h), these (obj_3) have only genes that biomaRt could find human orthologs for in BOTH mouse and chick, i.e. are present in all three species

reorder_features <- function(obj, genes) {
  cm <- GetAssayData(obj, layer = "counts")[genes, ]  # use slot="counts" on older Seurat
  CreateSeuratObject(counts = cm, meta.data = obj@meta.data, project = Project(obj))
}
mouse_3 <- reorder_features(mouse_h, common)
chick_3 <- reorder_features(chick_h, common)
human_3 <- reorder_features(human_h, common)


stopifnot(
  all(rownames(human_3) == common),
  all(rownames(mouse_3) == common),
  all(rownames(chick_3) == common)
)
identical(rownames(human_3), rownames(mouse_3))   # should be TRUE
identical(rownames(human_3), rownames(chick_3))   # should be TRUE. 



saveRDS(human_3, file = "~/amdata/hum_orthocommon_order2.rds")
saveRDS(mouse_3, file = "~/amdata/mouse_orthocommon_order2.rds")
saveRDS(chick_3, file = "~/amdata/chick_orthocommon_order2.rds")



################################Subsetting to keep only the foregut epithelium and foregut mesenchyme cells###################################### (8/13/26)


##HUMAN

hum <- readRDS("~/amdata/hum_orthocommon_order2.rds")
Idents(hum) <- "developmental system"
foregutonly <- subset(hum, idents = c("epithelium", "splanchnic LPM", "endoderm"))

saveRDS(foregutonly, file = "~/amdata/humac_foregutONLY_ortho2.rds")




##################INTEGRATE OBJECTS TOGETHER (note that I overwrote this on 9/3/2026 w/ ScaleData that regressed out cell cycle markers)

chick <- readRDS("~/amdata/chick_orthocommon_order2.rds")
hum <- readRDS("~/amdata/humac_foregutONLY_ortho2.rds")
mouse <- readRDS("~/amdata/mouse_orthocommon_order2.rds")


hum$species   <- "human"
mouse$species <- "mouse"
chick$species <- "chick"

combined <- merge(hum, y = list(mouse, chick))

combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

## split the RNA counts into one layer per species
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$species)

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

combined <- IntegrateLayers(
  object = combined,
  method = CCAIntegration,           
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  verbose = TRUE
)

saveRDS(combined, file = "~/amdata/integrated_intermediate.rds")


combined <- FindNeighbors(combined, reduction = "integrated.cca", dims = 1:30)
combined <- FindClusters(combined, resolution = c(0.1, 0.2, 0.3, 0.4))
combined <- RunUMAP(combined, reduction = "integrated.cca", dims = 1:30)


combined <- JoinLayers(combined)


saveRDS(combined, file = "~/amdata/integrated.rds")


DimPlot(combined, group.by = "species",  reduction = "umap", raster = FALSE)
DimPlot(combined, group.by = "RNA_snn_res.0.4",  reduction = "umap", raster = FALSE, label = TRUE)
DimPlot(combined, group.by = "RNA_snn_res.0.3",  reduction = "umap", raster = FALSE, label = TRUE)
DimPlot(combined, group.by = "RNA_snn_res.0.2",  reduction = "umap", raster = FALSE, label = TRUE)
DimPlot(combined, group.by = "RNA_snn_res.0.1",  reduction = "umap", raster = FALSE, label = TRUE)




##################CLUSTER ANNOTATION

DefaultAssay(combined) <- "RNA"

Idents(combined) <- "RNA_snn_res.0.2"

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
     group_by(cluster) %>%
     slice_max(n = 50, order_by = avg_log2FC) -> top50

DoHeatmap(subset(combined, downsample = 100), features = top50$gene, size = 3.5)


markers %>%
     group_by(cluster) %>%
     slice_max(n = 10, order_by = avg_log2FC) -> top10

DoHeatmap(subset(combined, downsample = 100), features = top10$gene, size = 3.5)


write.table(markers, "~/amdata/celltypeDEGs.txt")
write.table(top10, "~/amdata/top10markers.txt")
write.table(top50, "~/amdata/top50markers.txt")

##For mesenchymal cells
FeaturePlot(combined, features = "PDGFRA")
FeaturePlot(combined, features = "VIM")
FeaturePlot(combined, features = "FOXF1")
FeaturePlot(combined, features = "PTCH1")

##For epithelial cells
FeaturePlot(combined, features = "EPCAM")
FeaturePlot(combined, features = "CDH1")
FeaturePlot(combined, features = "FOXA1")
FeaturePlot(combined, features = "FOXA2")
FeaturePlot(combined, features = "NKX2-1")


##For cluster 3 to check that its not pure mesenchymal cells (https://pmc.ncbi.nlm.nih.gov/articles/PMC7453027/)
FeaturePlot(combined, features = "WT1")
FeaturePlot(combined, features = "GATA4")



#####Annotating Clusters

Idents(combined) <- "RNA_snn_res.0.2"

celltypelist <- c("cycling", "Foregut Mesenchyme", "human_only", "Mesothelium", "Foregut Epithelium", "Neural Crest", "Vascular Endothelium", "Macrophages")
names(celltypelist) <- levels(combined)
combined$seurat_clusters <- Idents(combined)   # or confirm it already exists
combined <- RenameIdents(combined, celltypelist)
combined$celltype <- Idents(combined)


saveRDS(combined, file = "~/amdata/celltype_annotated.rds")


table(cluster = combined$celltype, species = combined$species)


#####DEG ANALYSIS (8/26/26)

##DEGs for each species against both other species
epi <- subset(combined, idents = "Foregut Epithelium")

Idents(epi) <- "species"
all_sp <- FindAllMarkers(epi, min.pct = 0.25)   # one block per species
write.table(all_sp, file = "~/amdata/epi_speciesvall.txt")


##DEGs for each species against each other species individually
Idents(epi) <- "species"
hm <- FindMarkers(epi, ident.1 = "human", ident.2 = "mouse", min.pct = 0.25)
hc <- FindMarkers(epi, ident.1 = "human", ident.2 = "chick", min.pct = 0.25)
mc <- FindMarkers(epi, ident.1 = "mouse", ident.2 = "chick", min.pct = 0.25)

write.table(hm, file = "~/amdata/epi_humanvmouse.txt")
write.table(hc, file = "~/amdata/epi_humanvchick.txt")
write.table(mc, file = "~/amdata/epi_mousevchick.txt")



##DEGs for each species against both other species
mes <- subset(combined, idents = "Foregut Mesenchyme")

Idents(mes) <- "species"
all_sp <- FindAllMarkers(mes, min.pct = 0.25)   # one block per species
write.table(all_sp, file = "~/amdata/mes_speciesvall.txt")


##DEGs for each species against each other species individually
Idents(mes) <- "species"
hmm <- FindMarkers(mes, ident.1 = "human", ident.2 = "mouse", min.pct = 0.25)
hcm <- FindMarkers(mes, ident.1 = "human", ident.2 = "chick", min.pct = 0.25)
mcm <- FindMarkers(mes, ident.1 = "mouse", ident.2 = "chick", min.pct = 0.25)

write.table(hmm, file = "~/amdata/mes_humanvmouse.txt")
write.table(hcm, file = "~/amdata/mes_humanvchick.txt")
write.table(mcm, file = "~/amdata/mes_mousevchick.txt")



####Figures for presentation (9/3/26):

srat <- readRDS("~/amdata/celltype_annotated.rds")
combined <- subset(srat, subset = celltype != "human_only")

saveRDS(combined, file = "~/amdata/finalseurat.rds")

p <- DimPlot(combined, group.by = "RNA_snn_res.0.2",  reduction = "umap", raster = FALSE, label = TRUE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("cluster_label.pdf", plot = p, width = 10, height = 8, units = "in")

a <- DimPlot(combined, group.by = "RNA_snn_res.0.2",  reduction = "umap", raster = FALSE, label = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("cluster_nolabel.pdf", plot = a, width = 10, height = 8, units = "in")

p <- DimPlot(combined, group.by = "celltype",  reduction = "umap", raster = FALSE, label = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL) + NoLegend()
ggsave("celltype_nolegend.pdf", plot = p, width = 10, height = 8, units = "in")

p <- DimPlot(combined, group.by = "celltype",  reduction = "umap", raster = FALSE, label = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("celltype_legend.pdf", plot = p, width = 10, height = 8, units = "in")


Idents(srat) <- "celltype"
srat <- RenameIdents(srat, "cycling" = "Other Mesenchyme")
srat$celltype <- Idents(srat)   # write the change back into the celltype column
table(srat$celltype)

saveRDS(srat, file = "~/amdata/final_seurat2.rds") ##renamed to seurat_HuMsCh_finalRui.rds


p <- DimPlot(combined, group.by = "celltype",  reduction = "umap", raster = FALSE, label = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL) + NoLegend()
ggsave("celltype_nolegend.pdf", plot = p, width = 10, height = 8, units = "in")

p <- DimPlot(srat, group.by = "celltype",  reduction = "umap", raster = FALSE, label = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("celltype_legend.pdf", plot = p, width = 10, height = 8, units = "in")

p <- DimPlot(combined, group.by = "species",  reduction = "umap", raster = FALSE, label = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("species.pdf", plot = p, width = 10, height = 8, units = "in")


##Make list of all 10k features that were kept after gene orthology
genes <- rownames(srat[["RNA"]])
head(genes)
length(genes)          # how many genes



#####Re-cluster w/ regressing out cell cycle markers (9/3/26)


srat <- readRDS("~/amdata/seurat_HuMsCh_finalRui.rds")

DefaultAssay(srat) <- "RNA"
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat)
srat <- CellCycleScoring(object = srat, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes, set.ident = FALSE)

srat[["RNA"]] <- split(srat[["RNA"]], f = srat$species)

srat <- ScaleData(srat, vars.to.regress = c("S.Score", "G2M.Score"))
srat <- RunPCA(srat)

srat <- IntegrateLayers(
  object = srat,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",  verbose = TRUE
)

srat <- FindNeighbors(srat, reduction = "integrated.cca", dims = 1:30)
srat <- FindClusters(srat, resolution = c(0.1, 0.2, 0.3, 0.4))
srat <- RunUMAP(srat, reduction = "integrated.cca", dims = 1:30)

srat <- JoinLayers(srat)

saveRDS(srat, file = "~/amdata/ccreg_seurat.rds")

##For mesenchymal cells
FeaturePlot(srat, features = "PDGFRA")
FeaturePlot(srat, features = "VIM")
FeaturePlot(srat, features = "FOXF1")
FeaturePlot(srat, features = "PTCH1")

##For epithelial cells
FeaturePlot(srat, features = "EPCAM")
FeaturePlot(srat, features = "CDH1")
FeaturePlot(srat, features = "FOXA1")
FeaturePlot(srat, features = "FOXA2")
FeaturePlot(srat, features = "NKX2-1")

FeaturePlot(srat, features = "WT1")
FeaturePlot(srat, features = "GATA4")



Idents(srat) <- "RNA_snn_res.0.2"

celltypelist <- c("Ventral Subepithelial Mesenchyme", "Dorsal Subepithelial Mesenchyme", "Dorsal Peripheral Mesenchyme", "Foregut Epithelium", "Mesothelium", "Neural Crest", "Vascular Endothelium", "Immune Cells")
names(celltypelist) <- levels(srat)
srat <- RenameIdents(srat, celltypelist)
srat$cell_type <- Idents(srat)


srat$celltype <- NULL


saveRDS(srat, file = "~/amdata/seurat_HuChMs_regressed_Ruifinal.rds")



####Figures

p <- DimPlot(srat, group.by = "species",  reduction = "umap", raster = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("species.pdf", plot = p, width = 10, height = 8, units = "in")

p <- DimPlot(srat, group.by = "species", split.by = "species", reduction = "umap", raster = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("species_split.pdf", plot = p, width = 20, height = 8, units = "in")

p <- DimPlot(srat, group.by = "RNA_snn_res.0.2",  reduction = "umap", raster = FALSE, label = TRUE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("cluster.pdf", plot = p, width = 10, height = 8, units = "in")

p <- DimPlot(srat, group.by = "RNA_snn_res.0.2",  reduction = "umap", split.by = "species", raster = FALSE, label = TRUE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("cluster_split.pdf", plot = p, width = 20, height = 8, units = "in")

p <- DimPlot(srat, group.by = "cell_type",  reduction = "umap", raster = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("celltype.pdf", plot = p, width = 10, height = 8, units = "in")

p <- DimPlot(srat, group.by = "cell_type",  reduction = "umap", split.by = "species", raster = FALSE) + coord_fixed() + theme(aspect.ratio = 1) + ggtitle(NULL)
ggsave("celltype_split.pdf", plot = p, width = 20, height = 8, units = "in")


###Subsetting and Re-clustering Cluster 0 and Cluster 1

srat <- readRDS("~/amdata/For Rui/seurat_HuChMs_regressed_Ruifinal.rds")

Idents(srat) <- "cell_type"
subs <- subset(srat, idents = c("Ventral Subepithelial Mesenchyme", "Dorsal Subepithelial Mesenchyme"))

DefaultAssay(subs) <- "RNA"
subs <- NormalizeData(subs)
subs <- FindVariableFeatures(subs)
subs <- CellCycleScoring(object = subs, g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes, set.ident = FALSE)

subs[["RNA"]] <- split(subs[["RNA"]], f = subs$species)

subs <- ScaleData(subs, vars.to.regress = c("S.Score", "G2M.Score"))
subs <- RunPCA(subs)

subs <- IntegrateLayers(
  object = subs,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",  verbose = TRUE
)

subs <- FindNeighbors(subs, reduction = "integrated.cca", dims = 1:30)
subs <- FindClusters(subs, resolution = c(0.1, 0.2, 0.3, 0.4))
subs <- RunUMAP(subs, reduction = "integrated.cca", dims = 1:30)

subs <- JoinLayers(subs)

saveRDS(subs, file = "~/amdata/cluster0_1subset.rds")

#########################Post processing
obj_filtered <- subset(
    srat,
    subset = orig.ident != "Ch_E4_Cyclop"
)
meta <- obj_filtered[[]]

keep_cells <- rownames(meta)[
    is.na(meta$stage) |
        meta$stage != "CS15-16"
]

obj_filtered <- subset(obj_filtered, cells = keep_cells)
