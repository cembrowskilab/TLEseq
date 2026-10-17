################################################################################
#
# 01_microglia_subsetting_clustering_annotation.R
#
# Xenium analysis of microglia in human hippocampus samples
# Using MegaData_v8
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# Microglia subsetting, clustering, and cell-state annotation


### 0) Setup -------------------------------------------------------------------

# Core packages
suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(purrr)
        library(ggplot2)
        library(patchwork)
        library(harmony)
        library(clustree)
})

### PACKAGE CONFLICTS
conflicted::conflicts_prefer(dplyr::filter, 
                             dplyr::select, 
                             dplyr::rename, 
                             dplyr::lag, 
                             dplyr::count, 
                             dplyr::first)


# Memory for large objects (adjust if needed)
options(future.globals.maxSize = 100000 * 1024^2)


### 1) Helper functions --------------------------------------------------------
# Reuse an existing object in the global environment; otherwise load its RDS file.
# object_name: name to assign; file: user-configured RDS path; class: optional class check.
load_if_missing <- function(object_name, file, class = NULL) {
        if (exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
                obj <- get(object_name, envir = .GlobalEnv)
                if (is.null(class) || inherits(obj, class)) {
                        message("Reusing ", object_name, " from the environment")
                        return(obj)
                }
                warning("Existing ", object_name, " is not a ", class, "; reloading")
        }
        if (is.null(file) || length(file) != 1 || is.na(file) || !nzchar(file) || !file.exists(file))
                stop("Set a valid input path for '", object_name, "' in the user-editable configuration")
        obj <- readRDS(file)
        if (!is.null(class) && !inherits(obj, class)) stop("Input for '", object_name, "' must be a ", class)
        assign(object_name, obj, envir = .GlobalEnv)
        message("Loaded ", object_name, " from ", file)
        obj
}

# Prepare SCT model metadata for FindAllMarkers (original workflow).
prep_FindMarkers <- function(obj, num_slices = length(obj@assays$SCT@SCTModel.list)) {
        for (i in seq_len(num_slices)) slot(object = obj@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
        PrepSCTFindMarkers(obj)
}

### 2) User-editable paths and load inputs -------------------------------------
# Set input/output paths before running. Objects already in the session are reused.
cfg <- list(
        megadata_rds = NA_character_,
        region_metadata_rds = NA_character_,
        microglia_rds = NA_character_,
        markers_csv = NA_character_,
        humica_xlsx = NA_character_,
        output_dir_data = NA_character_,
        output_dir_plots = NA_character_)
if (is.na(cfg$output_dir_data)) stop("Set cfg$output_dir_data")
dir.create(cfg$output_dir_data, recursive = TRUE, showWarnings = FALSE)
if (is.na(cfg$output_dir_plots)) stop("Set cfg$output_dir_plots")
dir.create(cfg$output_dir_plots, recursive = TRUE, showWarnings = FALSE)

### 3) Color palettes -------------------------------------------------------

ilae_levels <- c("P", "0", "1", "2", "3")
pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00", "2" = "#00C0C5", "3" = "#E69F00")
pal_cluster_named <- c("Homeo" = "#ddcc77", "DIM" = "#117733", "BAM" = "#88ccee", 
                       "DAM" = "#cc6677", "Inflam.DAM" = "#aa4499", "T-cells" = "#332288", 
                       "PeriNeurM" = "#44aa99")
cluster_key <- data.frame(id = c("0", "1", "3", "6", "5", "4", "2"), label = names(pal_cluster_named), colour = unname(pal_cluster_named))
cluster_name_levels_noAll <- cluster_key$label
cluster_name_levels <- c(cluster_name_levels_noAll, "All")
box_colours <- c(pal_cluster_named, "All" = "white")

dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")

### 4) Subset and annotate source microglia ---------------------------------
# Load in MegaData_v8.rds and set the original cluster identities as the active Idents

Idents(dat) <- dat$seurat_clusters

# View data: Plot UMAPs coloured by cluster number and by annotated cell type, side by side
dat_dimplot1 <- DimPlot(dat, shuffle = TRUE, label = TRUE) +
        coord_fixed() + theme(legend.position = 'none') + ggtitle("MegaData Clustering - Ident")
dat_dimplot2 <- DimPlot(dat, shuffle = TRUE, group.by = 'celltype', label = TRUE) +
        coord_fixed() + theme(legend.position = 'none') + ggtitle("MegaData Clustering - Cell Type")
dat_dimplot1 + dat_dimplot2

# Sanity check that Cluster 2 is the microglia cluster: 
c(Cluster2 = sum(dat$seurat_clusters == 2), Micro = sum(dat$celltype == "Micro")) #Cluster 2 contains the same number of cells as the "Micro" annotated cluster = 297,653 cells 

# Subset the full dataset down to just the microglia-annotated cluster (cluster 2)
datM <- subset(dat, idents = 2)
nrow(datM@meta.data) # confirm datM contains 297,653 obs, matching the counts above

# Remove off-target cell types by further restricting the dataset to Micro/Cluster2 cells labeled "Micro/Macro/T" by SingleR mid.cell.class
sum(datM$SingleR_label_midfine == "Micro/Macro/T") #262,144 cells are "Micro/Macro/T"
datM_clean <- subset(datM, subset = SingleR_label_midfine == "Micro/Macro/T")
sum(datM_clean$SingleR_label_midfine == "Micro/Macro/T") #262,144 cells remain, matches above
class(datM_clean) # confirm this is still a Seurat object
rm(datM) # Remove intermediate file to save space

### 5) Add regional metadata ------------------------------------------------
# Load megadata_metadata_md_v8.rds, the file that contains anatomical domain/region information
dat_md <- load_if_missing("dat_md", cfg$region_metadata_rds)

# Filter that metadata down to just the cells present in the microglia subset
dat_md_subset <- dat_md[colnames(datM_clean), , drop = FALSE] #262,144 matches datM_clean above

# Attach the region metadata to the microglia subset
datM_clean$region_md <- dat_md_subset$region_md

# Sanity Check that merging was clean:
nrow(dat_md_subset) == ncol(datM_clean) # same number of cells
all(rownames(dat_md_subset) == colnames(datM_clean))  # same cells, same order
sum(is.na(dat_md_subset$region_md)) # cells with no region assigned

### 6) Normalize and integrate ----------------------------------------------
# Set the default assay to the raw Xenium counts before normalizing
DefaultAssay(datM_clean) <- 'Xenium'

# Normalize with SCTransform
datM_clean <- SCTransform(datM_clean, assay = "Xenium", min_cells = 0)

# Run PCA (100 PCs) on the SCT-normalized data
datM_clean <- RunPCA(datM_clean, npcs = 100)

# Batch-correct across samples with Harmony
datM_clean <- IntegrateLayers(
        object = datM_clean,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        verbose = TRUE,
        normalization.method = "SCT")

### 7) Parameter screening --------------------------------------------------
# Screen dimensions and resolutions using clustree plots.

# Create a new folder in the working directory for the output of the following tests:
rationale_dir <- file.path(cfg$output_dir_plots, "dim_res_clust")
dir.create(rationale_dir, showWarnings = FALSE)

# Candidate dimensions and clustering resolutions for parameter screening.
dims_list   <- list(1:10, 1:15, 1:20, 1:30, 1:35, 1:40)
resolutions <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)

# Make SCT the active assay to set the cluster column prefix, "SCT_snn_res."
DefaultAssay(datM_clean) <- "SCT"
res_prefix <- paste0(DefaultAssay(datM_clean), "_snn_res.")

# For each set of dimensions: cluster at every resolution, then save a clustree
cluster_counts <- list()

for (dims in dims_list) {
        max_dim <- max(dims)
        message("Clustering dims 1:", max_dim)

        dat_mod <- datM_clean %>%
                FindNeighbors(dims = dims, reduction = "harmony") %>%
                FindClusters(resolution = resolutions)   # one metadata column per resolution

        cluster_counts[[paste0("Dim_", max_dim)]] <-
                sapply(dat_mod@meta.data[paste0(res_prefix, resolutions)], n_distinct)

        ggsave(file.path(rationale_dir, paste0("Clustree_Dim_", max_dim, ".pdf")),
               clustree(dat_mod, prefix = res_prefix) + ggtitle(paste("Clustree - Dimensions:", max_dim)),
               width = 10, height = 10)

        rm(dat_mod)
        gc()  # Free memory before testing the next dimensionality.
}

# After viewing the clustrees above, the following dimensions and resolutions were chosen for further evaluation
chosen_combos <- data.frame(max_dim = c(10, 15, 20, 30, 35, 40),
                            res     = c(0.1, 0.1, 0.1, 0.05, 0.05, 0.05))

# FindAllMarkers on SCT data needs the SCT models prepared first

datM_clean <- prep_FindMarkers(datM_clean)

# For each combination: cluster, make the UMAP, find markers, and draw top-10 marker heatmaps
for (i in seq_len(nrow(chosen_combos))) {
        max_dim     <- chosen_combos$max_dim[i]
        res         <- chosen_combos$res[i]
        dims        <- 1:max_dim
        res_col     <- paste0(res_prefix, res)
        combo_label <- paste0("Dim", max_dim, "_Res", res)
        message("Running: ", combo_label)

        # Neighbours, clusters, and UMAP for this combination
        dat_mod <- datM_clean %>%
                FindNeighbors(dims = dims, reduction = "harmony") %>%
                FindClusters(resolution = res) %>%
                RunUMAP(dims = dims, reduction = "harmony")
        Idents(dat_mod) <- res_col

        # UMAP coloured by cluster
        ggsave(file.path(rationale_dir, paste0("UMAP_Dim_", max_dim, "_Res_", res, ".pdf")),
               DimPlot(dat_mod, shuffle = TRUE, label = TRUE) + coord_fixed() +
                       ggtitle(paste("Dims:", max_dim, "| Res:", res)),
               width = 7, height = 6)

        # Cluster markers (reused from the CSV if a previous run already made them; delete it to recompute)
        marker_file <- file.path(rationale_dir, paste0("FindAllMarkers_", combo_label, ".csv"))
        if (file.exists(marker_file)) {
                markers_combo <- read.csv(marker_file)
        } else {
                markers_combo <- FindAllMarkers(dat_mod, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
                write.csv(markers_combo, marker_file, row.names = FALSE)
        }

        # Top-10 markers per cluster, one heatmap page per cluster
        top10 <- markers_combo %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()

        pdf(file.path(rationale_dir, paste0(combo_label, "_ClusterHeatmaps.pdf")), width = 12, height = 8)
        for (cl in sort(unique(top10$cluster))) {
                print(DoHeatmap(dat_mod, features = top10$gene[top10$cluster == cl], group.by = res_col) +
                              NoLegend() + ggtitle(paste0(combo_label, " - Cluster ", cl, " - Top Markers")))
        }
        dev.off()

        rm(dat_mod)
        gc()  # Free memory before testing the next parameter combination.
}


### 8) Final microglia clustering -------------------------------------------
# Perform UMAP clustering using the dims and resolution chosen from outputs above
datM_clean <- datM_clean %>%
        FindNeighbors(dims = 1:10, reduction = "harmony") %>%
        FindClusters(resolution = 0.1) %>%
        RunUMAP(dims = 1:10, reduction = "harmony")

# Find cluster marker genes (SCT models were already prepared with prep_FindMarkers() above)
markers <- FindAllMarkers(datM_clean,
                          only.pos = TRUE,
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

# Checkpoint: save the processed object so this pipeline doesn't need to be re-run from scratch if the R session crashes
saveRDS(datM_clean, file.path(cfg$output_dir_data, "datM_clean_v8_unannotated.rds"))

### 9) Assign microglia states ----------------------------------------------
# Microglia subcluster-specific colour palette and labels (shared across all microglia figures)
# Cluster labels were assigned from marker heatmaps and refined using
# HuMicA reference comparisons and perineuronal microglia analyses.


# Add the cluster name column to the metadata
datM_clean$cluster_names <- factor(
        cluster_key$label[match(as.character(datM_clean$seurat_clusters), cluster_key$id)],
        levels = cluster_name_levels_noAll)
write.csv(markers, cfg$markers_csv, row.names = FALSE)
saveRDS(datM_clean, file.path(cfg$output_dir_data, "datM_clean_v8_annotated.rds"))
