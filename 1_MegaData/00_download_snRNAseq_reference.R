################################################################################
#
# Xenium human hippocampus snRNA-seq reference (Seurat v5)
# From Thompson et al, 2025, NatNeuro
# downloaded from https://christinehou11.github.io/humanHippocampus2024/articles/humanHippocampus2024.html#package-tutorial

# Download and inspect ExperimentHub EH9606; optional independent clustering QC.
# Original analysis: Larissa Kraus
# Code cleanup and standardization: Larissa Kraus (October 2026)
######

### 0) Setup -------------------------------------------------------------------
suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(ggplot2)
        library(ExperimentHub)
})
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)
options(future.globals.maxSize = 100000 * 1024^2)

### IF NOT DONE, INSTALL PACKAGES
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

stopifnot(BiocManager::version() == "3.20")

BiocManager::install("SpatialExperiment")

BiocManager::install("humanHippocampus2024") # If it does not work, do developer version instead as below

# dowload developer version
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("christinehou11/humanHippocampus2024")

### 1) User-editable paths and load inputs -------------------------------------
output_rds = NA_character_

if (is.na(output_rds)) stop("Set output_rds")
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
# Install humanHippocampus2024 and ExperimentHub separately before running.


### 2) Download and reference QC ------------------------------------------------
###CONNECT TO DATAHUB, LOAD DATA
# Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()

## Load the datasets of the package
myfiles <- ExperimentHub::query(ehub, "humanHippocampus2024")

#Download single cell data
dat.sce <- myfiles[["EH9606"]]

#convert file to seurat object
sc.ref <- as.Seurat(dat.sce, counts = "counts", data = "logcounts")

saveRDS(sc.ref, output_rds)


###CHECK QC, CLUSTERING AND ANNOTATIONS
sc.ref[["percent.mt"]] <- PercentageFeatureSet(object = sc.ref, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(sc.ref, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3, pt.size = 0)
FeatureScatter(sc.ref, feature1 = "nCount_originalexp", feature2 = "nFeature_originalexp")

DimPlot(sc.ref, reduction = 'umap', group.by = 'broad.cell.class')
DimPlot(sc.ref, reduction = 'umap', group.by = 'mid.cell.class')
DimPlot(sc.ref, reduction = 'umap', group.by = 'fine.cell.class', label = T)
DimPlot(sc.ref, reduction = 'umap', group.by = 'Sample', shuffle=T)
FeaturePlot(sc.ref, features = 'percent.mt')
FeaturePlot(sc.ref, features = 'nCount_originalexp')
FeaturePlot(sc.ref, features = 'nFeature_originalexp')
FeatureScatter(sc.ref, feature1 = "nCount_originalexp", feature2 = "nFeature_originalexp")

# run additional clustering to sanity check
sc.ref <- NormalizeData(sc.ref)
sc.ref <- FindVariableFeatures(sc.ref)
sc.ref <- ScaleData(sc.ref)
sc.ref <- RunPCA(sc.ref)
ElbowPlot(sc.ref,ndims = 50)+geom_hline(yintercept = 2.5)

sc.ref <- FindNeighbors(sc.ref, dims = 1:40)
sc.ref <- FindClusters(sc.ref, resolution = 0.6)
sc.ref <- RunUMAP(sc.ref, dims = 1:40)

DimPlot(sc.ref, reduction = 'umap', shuffle=T)
DimPlot(sc.ref, reduction = 'umap', group.by = 'broad.cell.class', label = T)
DimPlot(sc.ref, reduction = 'umap', group.by = 'mid.cell.class', label = T)
DimPlot(sc.ref, reduction = 'umap', group.by = 'fine.cell.class', label = T)
DimPlot(sc.ref, reduction = 'umap', group.by = 'Sample', shuffle=T)

