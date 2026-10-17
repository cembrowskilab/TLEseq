# 01_inn_clustering_annotation.R
#
# Xenium analysis of inhibitory neurons in human hippocampal samples
# Input: annotated MegaData_v8 Seurat object (see cfg$megadata_rds)
# LK, October 2026
#
# This script:
# 1. Extracts the initial inhibitory-neuron cluster from the full Xenium dataset
# 2. Iteratively removes off-target clusters, re-normalizing and re-integrating after each subset
# 3. Compares clustering dimensions and resolutions
# 4. Performs final clustering and assigns fine and broad inhibitory-neuron subtypes
# 5. Saves the annotated inhibitory-neuron Seurat object
#
# Reproducibility notes:
# - Cluster exclusions and subtype assignments below are specific to MegaData_v8.
# - Each round's cluster IDs refer to the clustering from the preceding round.
# - Inspect the diagnostic UMAPs before reusing these choices
# - Run from a Seurat v5 environment

### Setup ----------------------------------------------------------------------

suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
        library(ggplot2)
        library(patchwork)
        library(beepr)
})

options(future.globals.maxSize = 100000 * 1024^2)

### User-editable configuration -------------------------------------------------

# Replace the example paths with local paths before running the script.
# Diagnostic PDFs are saved separately from the final annotated Seurat object.
cfg <- list(
        megadata_rds = "path/to/directory/MegaData_v8_annotated_final.rds",
        output_dir_data = "path/to/directory/data",
        output_dir_plots = "path/to/directory/plots/InN_v8/dimTest"
)

if (!file.exists(cfg$megadata_rds))
        stop("Input Seurat object not found: ", cfg$megadata_rds)

dir.create(cfg$output_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$output_dir_plots, recursive = TRUE, showWarnings = FALSE)

### Functions ------------------------------------------------------------------

# Run graph-based clustering and UMAP using an existing reduction (Harmony here).
# This function does not repeat normalization or integration.
Xenium_SeuratLK <- function(x, DIMs = 1:15, cluster.resolution = 0.5, red = "harmony", playMario = FALSE) {
        x <- FindNeighbors(x, dims = DIMs, reduction = red)
        x <- FindClusters(x, resolution = cluster.resolution)
        x <- RunUMAP(x, dims = DIMs, reduction = red)
        if (playMario) beep(sound = "mario")
        return(x)
}

# Compare candidate PC ranges and clustering resolutions using the same input object.
# Each PDF shows clusters followed by selected metadata overlays for interpretation.
# UMAP is recomputed for each dimension set, but reused across its resolutions.
dim_reso_test <- function(dataset, dims_list, resolutions, output_dir, red,
                          group_vars = c("orig.ident", "ILAE_score", "predicted.celltype", "FFPE_FF", "type")) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        missing_vars <- setdiff(group_vars, colnames(dataset@meta.data))
        if (length(missing_vars)) stop("Missing metadata columns for diagnostic plots: ", paste(missing_vars, collapse = ", "))
        
        for (dims in dims_list) {
                dat_mod <- FindNeighbors(dataset, dims = dims, reduction = red)
                dat_mod <- RunUMAP(dat_mod, dims = dims, reduction = red)
                
                for (res in resolutions) {
                        dat_res <- FindClusters(dat_mod, resolution = res)
                        plots <- list(DimPlot(dat_res, shuffle = TRUE, label = TRUE) + coord_fixed())
                        for (group_var in group_vars) {
                                plots <- append(plots, list(DimPlot(dat_res, group.by = group_var, shuffle = TRUE) + coord_fixed()))
                        }
                        
                        combined_plot <- wrap_plots(plots, ncol = 1) +
                                plot_annotation(title = paste("Dimensions:", paste(range(dims), collapse = "-"), "| Resolution:", res))
                        file_name <- paste0("Dim_", paste(range(dims), collapse = "-"), "_Res_", res, ".pdf")
                        pdf_path <- file.path(output_dir, file_name)
                        pdf(pdf_path, width = 8, height = 5 * length(plots))
                        tryCatch(print(combined_plot), finally = dev.off())
                }
        }
}

# Subset the current cluster identities, then repeat SCT normalization, PCA and
# Harmony integration. Optionally cluster/plot the subset before the next round.
# When run.SeuratLK = FALSE, the returned object is integrated but not reclustered;
# this is used immediately before the final dimension/resolution comparison.
subcluster_flow <- function(seurat_obj, cluster_id, DIMs = 1:40, cluster.resolution = 0.1,
                            red = "harmony", playMario = TRUE,
                            group.by.p2 = "SingleR_pruned_superfine", run.SeuratLK = TRUE) {
        missing_clusters <- setdiff(as.character(cluster_id), levels(Idents(seurat_obj)))
        if (length(missing_clusters)) stop("Requested clusters absent from current identities: ", paste(missing_clusters, collapse = ", "))
        message("Subsetting clusters: ", paste(cluster_id, collapse = ", "))
        sub_obj <- subset(seurat_obj, idents = cluster_id)
        message("Subset: ", ncol(sub_obj), " cells; ", nrow(sub_obj), " features")
        
        DefaultAssay(sub_obj) <- "Xenium"
        sub_obj <- SCTransform(sub_obj, assay = "Xenium", min_cells = 0)
        sub_obj <- RunPCA(sub_obj, npcs = 100)
        sub_obj <- IntegrateLayers(object = sub_obj, method = HarmonyIntegration,
                                   orig.reduction = "pca", new.reduction = red,
                                   verbose = TRUE, normalization.method = "SCT")
        
        if (!run.SeuratLK) {
                message("Returning normalized and integrated object without new clusters/UMAP")
                return(sub_obj)
        }
        
        sub_obj <- FindNeighbors(sub_obj, dims = DIMs, reduction = red)
        sub_obj <- FindClusters(sub_obj, resolution = cluster.resolution)
        sub_obj <- RunUMAP(sub_obj, dims = DIMs, reduction = red)
        if (playMario) beep(sound = "mario")
        
        if (!group.by.p2 %in% colnames(sub_obj@meta.data)) stop("Metadata column not found: ", group.by.p2)
        p1 <- DimPlot(sub_obj, shuffle = TRUE, label = TRUE) + coord_fixed() + theme(legend.position = "none")
        p2 <- DimPlot(sub_obj, shuffle = TRUE, group.by = group.by.p2, label = TRUE) +
                coord_fixed() + theme(legend.position = "none")
        print(p1 + p2)
        return(sub_obj)
}

### Load full Xenium dataset ---------------------------------------------------

# The initial inhibitory-neuron population corresponds to cluster 6 in MegaData_v8.
# Set identities explicitly because subsequent subsetting uses Idents(), not celltype.
dat <- readRDS(cfg$megadata_rds)

if (!"seurat_clusters" %in% colnames(dat@meta.data)) 
        stop("Input object requires seurat_clusters metadata")

if (!"SingleR_pruned_superfine" %in% colnames(dat@meta.data)) 
        stop("Input object requires SingleR_pruned_superfine metadata")

Idents(dat) <- "seurat_clusters"

### Inhibitory-neuron subsetting ---------------------------------------

# Round 1: extract cluster 6, then normalize/integrate and recluster with 40
# Harmony dimensions at resolution 0.5 to reveal off-target populations.
datIN <- subcluster_flow(seurat_obj = dat, cluster_id = 6, DIMs = 1:40,
                         cluster.resolution = 0.4, red = "harmony", playMario = TRUE,
                         group.by.p2 = "SingleR_pruned_superfine", run.SeuratLK = TRUE)

# Round 2: remove off target clusters (adjust the number below as needed)
# Retained cells are reprocessed and reclustered (over-clustering to better subset) 
toKeep<-setdiff(0:20, c(0, 8, 9, 12,13, 14, 19,20))

datIN <- subcluster_flow(seurat_obj = datIN, 
                         cluster_id = toKeep, 
                         DIMs = 1:40,
                         cluster.resolution = 0.4, 
                         red = "harmony", 
                         playMario = TRUE,
                         group.by.p2 = "SingleR_pruned_superfine",
                         run.SeuratLK = TRUE)

# Round 3: remove off target clusters (adjust the number below as needed)
# Retained cells are reprocessed; defer reclustering until parameter testing.
toKeep<-setdiff(0:18, c(14, 18))

datIN <- subcluster_flow(seurat_obj = datIN, cluster_id = toKeep,
                         playMario = TRUE, group.by.p2 = "SingleR_pruned_superfine",
                         run.SeuratLK = FALSE)

### Clustering parameter optimization ------------------------------------------

# Compare multiple dimensions and resolutions
# PDFs show cluster identities, SingleR predictions, ILAE score and patient ID.
dim_reso_test(dataset = datIN, dims_list = list(1:10, 1:20, 1:30),
              resolutions = c(0.2, 0.3, 0.4, 0.5), 
              output_dir = here('plots', 'InN_v8_github', 'dim_reso_test'),
              red = "harmony", 
              group_vars = c("SingleR_pruned_superfine"))

### Final clustering -----------------------------------------------------------

# Selected parameters: 20 dims and resolution 0.2
# IMPORTANT: subtype annotations below depend on these final cluster IDs.
datIN <- Xenium_SeuratLK(x = datIN, DIMs = 1:15, cluster.resolution = 0.4,
                         red = "harmony", playMario = FALSE)

p1 <- DimPlot(datIN, shuffle = TRUE, label = TRUE) + coord_fixed() + theme(legend.position = "none")
p2 <- DimPlot(datIN, shuffle = TRUE, group.by = "SingleR_pruned_superfine", label = TRUE) +
        coord_fixed() + theme(legend.position = "none")

print(p1 + p2)

### Optional plotting approaches to identify cluster identity -------------------
# Collection of different ways to look at the data and decide which subtypes underly clusters
# overall clusters were identified by SingleR label, marker gene expression and spatial location

##plot SingleR label clusters individually to find what clusters are
# Make a copy of the vector and assign "Unlabeled" to NAs
datIN2<-datIN
labels <- datIN2$SingleR_pruned_superfine
labels[is.na(labels)] <- "Unlabeled"
datIN2$SingleR_labeled <- labels  # Add new metadata column

# Set identities to new column with "Unlabeled"
Idents(datIN2) <- datIN2$SingleR_labeled

# Get unique labels (now including "Unlabeled")
unique_labels <- unique(datIN2$SingleR_labeled)

# Initialize list for plots
plot_list <- list()

for (label in unique_labels) {
        cells <- WhichCells(datIN2, idents = label)
        
        if (length(cells) > 0) {
                p <- DimPlot(datIN2,
                             cells.highlight = cells,
                             cols.highlight = "red",
                             sizes.highlight = 0.2, pt.size=0.1) +
                        ggtitle(label) + coord_fixed()+
                        theme(legend.position='none')+
                        theme_minimal()
                plot_list[[label]] <- p
        }
}

###check proportions of SingleR label transfer for each cluster
# Ensure metadata is a proper data frame
meta <- as.data.frame(datIN@meta.data)

# Convert columns to character and replace NAs with "Unassigned"
meta$seurat_clusters <- as.character(meta$seurat_clusters)
meta$SingleR_pruned_superfine <- as.character(meta$SingleR_pruned_superfine)
meta$SingleR_pruned_superfine[is.na(meta$SingleR_pruned_superfine)] <- "Unassigned"

# Count number of cells per SingleR type per Seurat cluster
count_table <- dplyr::count(meta, seurat_clusters, SingleR_pruned_superfine, name = "n")

# Calculate proportions
prop_table <- count_table %>%
        group_by(seurat_clusters) %>%
        mutate(proportion = n / sum(n)) %>%
        ungroup()

# Get top 5 SingleR labels per cluster
top5_table <- prop_table %>%
        arrange(seurat_clusters, desc(proportion)) %>%
        group_by(seurat_clusters) %>%
        slice_head(n = 5) %>%
        ungroup() %>%
        mutate(proportion = round(proportion, 5))

###Check canonical marker gene expression 
genes_to_plot<-c("SLC17A7","GAD1","SST",'NPY',"PVALB",'LHX6',"LAMP5",'KIT','CXCL14',"CCK","VIP", "ADARB2")
genes_to_plot<-c('TAC1','SST','PVALB', 'KIT', 'ADARB2', 'VIP', 'LAMP5', 'LHX6', 'TAC1', 'NPY')


Idents(datIN)<-'broad_subtype_InhN'
p2<-DotPlot(datIN, features = genes_to_plot) +
        scale_color_gradientn(colors = c("cyan", "grey", "red")) +
        theme(legend.position = "top",
              legend.justification = "center",
              axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 13))

### Annotate subtypes ----------------------------------------------------------

# Manual mapping of the 15 final clusters to fine inhibitory-neuron subtypes.
recode_vector <- c(
        "0" = "VIP_1",
        "1" = "LAMP5.MGE",
        "2" = "VIP_2",
        "3" = "PV_1",
        "4" = "SST_1",
        "5" = "CXCL14",
        "6" = "LAMP5.CGE",
        "7" = "SST_2",
        "8" = "HTR3A_1",
        "9" = "HTR3A_2",
        "10" = "PV_2",
        "11" = "HTR3A_3",
        "12" = "C1QL1",
        "13" = "VIP_3",
        "14" = "Cajal")

datIN$fine_subtype_InhN <- recode(datIN$seurat_clusters, !!!recode_vector)

# optional check of new label
# p1<-DimPlot(datIN, shuffle = TRUE, group.by = 'fine_subtype_InhN',label=T)+
#         coord_fixed()+theme(legend.position='none')
# p2<-DimPlot(datIN, shuffle = TRUE, group.by = 'SingleR_pruned_superfine', label=T)+
#         coord_fixed()+theme(legend.position='none')
# p1+p2

# Manual mapping of the 15 final clusters to broad inhibitory-neuron subtypes.
recode_vector <- c(
        "0" = "VIP",
        "1" = "LAMP5.MGE",
        "2" = "VIP",
        "3" = "PV",
        "4" = "SST",
        "5" = "CXCL14",
        "6" = "LAMP5.CGE",
        "7" = "SST",
        "8" = "HTR3A",
        "9" = "HTR3A",
        "10" = "PV",
        "11" = "HTR3A",
        "12" = "C1QL1",
        "13" = "VIP",
        "14" = "Cajal")

datIN$broad_subtype_InhN <- recode(datIN$seurat_clusters, !!!recode_vector)

#optional check of new label
# p1<-DimPlot(datIN, shuffle = TRUE, group.by = 'broad_subtype_InhN',label=T)+
#         coord_fixed()+theme(legend.position='none')
# p2<-DimPlot(datIN, shuffle = TRUE, group.by = 'fine_subtype_InhN', label=T)+
#         coord_fixed()+theme(legend.position='none')
# p1+p2


### Save annotated object ------------------------------------------------------

# Downstream InN scripts should load this object, not repeat clustering.
output_rds <- file.path(cfg$output_dir_data, "datIN_v8.rds")
saveRDS(datIN, output_rds)
