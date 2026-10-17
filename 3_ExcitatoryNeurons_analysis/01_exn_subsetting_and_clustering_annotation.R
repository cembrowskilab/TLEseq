
# 01_exn_subsetting_and_clustering_annotation.R
#
# Xenium analysis of excitatory neurons in human hippocampus samples
# Using MegaData_v8
# Author: Larissa Kraus
# (October 2026)
#
# This script:
# 1. Subsets excitatory neurons from the full Xenium dataset (Megadata)
# 2. Iteratively removes off-target cell populations
# 3. Tests clustering dimensions and resolutions
# 4. Performs final clustering and assigns ExN subtype annotations
# 5. Saves the annotated excitatory-neuron Seurat object


### LOAD PACKAGES ##############################################################
suppressPackageStartupMessages({
        library(dplyr)
        library(patchwork)
        library(ggplot2)
        library(Seurat)
        library(beepr)
        library(conflicted)
})

conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")
conflicted::conflict_prefer("lag", "dplyr")
conflicted::conflict_prefer("first", "dplyr")

# Increase the maximum size of global objects exported to future workers
# when processing large Seurat objects.
options(future.globals.maxSize = 100000 * 1024^2)


### HELPER FUNCTIONS ###########################################################
# Reuse an existing object; otherwise load it from a user-configured RDS path.
# Returns the object without overwriting an existing object in the calling environment.
load_if_missing <- function(object_name, path, envir = .GlobalEnv) {
        if (exists(object_name, envir = envir, inherits = FALSE)) {
                message("Reusing ", object_name, " from the environment.")
                return(get(object_name, envir = envir, inherits = FALSE))
        }
        if (is.null(path) || !nzchar(path) || !file.exists(path)) {
                stop("Missing ", object_name, ". Set its input RDS path in cfg: ", path)
        }
        message("Loading ", object_name, " from: ", path)
        readRDS(path)
}



# Subset a Seurat object using a metadata column, then repeat normalization,
# dimensionality reduction, and Harmony integration.
# Clustering and UMAP can be skipped for subsequent parameter testing.
subcluster_flow <- function(seurat_obj, subsetIdents = "seurat_clusters",
                            cluster_id, DIMs = 1:40, cluster.resolution = 0.1,
                            red = "harmony", playMario = TRUE,
                            group.by.p2 = "SingleR_pruned_superfine",
                            run.SeuratLK = TRUE) {
        
        # Select identities and subset the requested cells
        Idents(seurat_obj) <- seurat_obj@meta.data[[subsetIdents]]
        print("Using Idents")
        print(unique(Idents(seurat_obj)))
        
        print("Subsetting takes a sec...")
        sub_obj <- subset(seurat_obj, idents = cluster_id)
        cat("Subset dimensions:\n")
        print(dim(sub_obj))
        
        # Recalculate SCT normalization and PCA for the subset
        DefaultAssay(sub_obj) <- "Xenium"
        sub_obj <- SCTransform(sub_obj, assay = "Xenium", min_cells = 0)
        sub_obj <- RunPCA(sub_obj, npcs = 100)
        
        # Integrate using Harmony
        sub_obj <- IntegrateLayers(object = sub_obj,
                                   method = HarmonyIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = red,
                                   verbose = TRUE,
                                   normalization.method = "SCT")
        
        # Return the integrated object without clustering when requested
        if (!run.SeuratLK) {
                message("Skipping clustering and UMAP. Returning integrated object only.")
                return(sub_obj)
        }
        
        # Cluster cells and calculate UMAP
        sub_obj <- FindNeighbors(sub_obj, dims = DIMs, reduction = red)
        sub_obj <- FindClusters(sub_obj, resolution = cluster.resolution)
        sub_obj <- RunUMAP(sub_obj, dims = DIMs, reduction = red)
        
        if (playMario) beep(sound = "mario")
        
        # Check that the requested annotation exists
        if (!(group.by.p2 %in% colnames(sub_obj@meta.data))) {
                stop(paste0("'", group.by.p2, "' not found in metadata columns."))
        }
        
        # Compare unsupervised clusters with existing cell-type annotations
        p1 <- DimPlot(sub_obj, shuffle = TRUE, label = TRUE) +
                coord_fixed() + theme(legend.position = "none")
        
        p2 <- DimPlot(sub_obj, shuffle = TRUE, group.by = group.by.p2,
                      label = TRUE) +
                coord_fixed() + theme(legend.position = "none")
        
        print(p1 + p2)
        
        return(sub_obj)
}


# Test combinations of principal components and clustering resolutions.
# For each combination, save a PDF showing clusters and selected metadata from group_vars (can be any metadata).
dim_reso_test <- function(dataset, dims_list, resolutions, output_dir, red,
                          playMario=TRUE,
                          group_vars = c("orig.ident", "ILAE_score",
                                         "SingleR_pruned_superfine")) {
        
        if (!dir.exists(output_dir)) {
                dir.create(output_dir, recursive = TRUE)
        }
        
        for (dims in dims_list) {
                max_dim <- max(dims)
                
                dat_mod <- dataset
                dat_mod <- FindNeighbors(dat_mod, dims = dims, reduction = red)
                dat_mod <- RunUMAP(dat_mod, dims = dims, reduction = red)
                
                for (res in resolutions) {
                        dat_res <- FindClusters(dat_mod, resolution = res)
                        
                        plots <- list(
                                DimPlot(dat_res, shuffle = TRUE, label = TRUE) +
                                        coord_fixed()
                        )
                        
                        for (group_var in group_vars) {
                                plots <- append(
                                        plots,
                                        list(DimPlot(dat_res, group.by = group_var,
                                                     shuffle = TRUE, raster = T) + coord_fixed())
                                )
                        }
                        
                        combined_plot <- wrap_plots(plots, ncol = 1)
                        
                        plot_title <- paste("Dimensions:", max_dim,
                                            "| Resolution:", res)
                        combined_plot <- combined_plot +
                                plot_annotation(title = plot_title)
                        
                        file_name <- paste0("Dim_", max_dim, "_Res_", res, ".pdf")
                        pdf_path <- file.path(output_dir, file_name)
                        
                        pdf_width <- 8
                        pdf_height <- 5 * length(plots)
                        
                        pdf(pdf_path, width = pdf_width, height = pdf_height)
                        print(combined_plot)
                        dev.off()
                }
        }
       
         if (playMario) beep(sound = "mario")
}


# Run clustering and UMAP on an already normalized and integrated object.
Xenium_SeuratLK <- function(x, DIMs = 1:15, cluster.resolution = 0.5,
                            red = "harmony", playMario = FALSE) {
        
        x <- FindNeighbors(x, dims = DIMs, reduction = red)
        x <- FindClusters(x, resolution = cluster.resolution)
        x <- RunUMAP(x, dims = DIMs, reduction = red)
        
        if (playMario) beep(sound = "mario")
        
        return(x)
}


### USER-EDITABLE PATHS ########################################################
cfg <- list(
        # Annotated MegaData object available from GEO
        megadata_rds = "path/to/directory/MegaData_v8_annotated.rds",
        
        # outputs
        output_dir_plots = "path/to/directory/plots",
        output_dir_data = "path/to/directory/data"
        )

dir.create(cfg$output_dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$output_dir_data, recursive = TRUE, showWarnings = FALSE)

### SUBSET EXCITATORY NEURONS ##################################################

# Load the full Xenium dataset and inspect unsupervised clusters
# alongside previously assigned broad cell-type annotations.
dat <- load_if_missing("dat", cfg$megadata_rds)

# p1 <- DimPlot(dat, shuffle = TRUE, label = TRUE) +
#         coord_fixed() + theme(legend.position = "none")
# 
# p2 <- DimPlot(dat, shuffle = TRUE, group.by = "celltype", label = TRUE) +
#         coord_fixed() + theme(legend.position = "none")
# 
# p1 + p2
# 

# Initial ExN subset: cluster 3 of the full dataset
Idents(dat) <- dat$seurat_clusters
datEx <- subset(dat, idents = 3)

# Retain cells annotated as ExcN by SingleR to remove remaining
# non-excitatory cells. Recalculate normalization, PCA, Harmony,
# clustering, and UMAP. Choose high DIMs and resolution to better identify off target cells
datEx <- subcluster_flow(seurat_obj = datEx,
                         subsetIdents = "SingleR_pruned_midfine",
                         cluster_id = "ExcN",
                         red = "harmony",
                         run.SeuratLK = TRUE,
                         DIMs = 1:25,
                         cluster.resolution = 0.3,
                         playMario = TRUE,
                         group.by.p2 = "SingleR_pruned_superfine")

# Remove cluster 5, identified as off-target during intermediate inspection.
# Recalculate SCT, PCA, and Harmony without clustering or UMAP.
toKeep <- setdiff(0:13, 5)

datEx <- subcluster_flow(seurat_obj = datEx,
                         subsetIdents = "seurat_clusters",
                         cluster_id = toKeep,
                         playMario = FALSE,
                         run.SeuratLK = FALSE)

### CLUSTERING PARAMETER TESTING ##############################################

# Compare 20, 25, 30, and 35 dimensions across resolutions 0.1-0.3.
# PDFs are saved for visual evaluation of cluster separation and annotations.
dim_reso_test(dataset = datEx,
              dims_list = list(1:20, 1:25, 1:30, 1:35),
              resolutions = c(0.1, 0.2, 0.3),
              output_dir = file.path(cfg$output_dir_plots, "ExN_newTest", "dimresTest"),
              red = "harmony", playMario=TRUE,
              group_vars = c("SingleR_pruned_superfine"))


### FINAL EXN CLUSTERING #######################################################

# Selected clustering parameters: 25 dimensions, resolution 0.2
datEx <- Xenium_SeuratLK(x = datEx, DIMs = 1:25,
                         cluster.resolution = 0.2,
                         red = "harmony", playMario = TRUE)

# Inspect final clusters alongside SingleR annotations
p1 <- DimPlot(datEx, shuffle = TRUE, label = TRUE) +
        coord_fixed() + theme(legend.position = "none")

p2 <- DimPlot(datEx, shuffle = TRUE, group.by = "SingleR_pruned_superfine",
              label = TRUE) +
        coord_fixed() + theme(legend.position = "none")

p1 + p2


### ASSIGN EXN SUBTYPE ANNOTATIONS #############################################

# Assign fine ExN subtypes to the final unsupervised clusters.
# Granule cells (GC) are retained as three separate populations.
recode_vector <- c(
        "0" = "GC_1",
        "1" = "dCA1",
        "2" = "CA2-3",
        "3" = "sCA1",
        "4" = "GC_2",
        "5" = "dSub",
        "6" = "sSub",
        "7" = "GC_3",
        "8" = "mSub",
        "9" = "L6b",
        "10" = "GC_4",
        "11" = "MC",
        "12" = "dCA3"
)

datEx$fine_subtype_ExN <- dplyr::recode(
        as.character(datEx$seurat_clusters), !!!recode_vector
)


# Assign broad ExN subtypes by combining the GC sub-populations.
# All other subtype annotations remain unchanged.
recode_vector <- c(
        "0" = "GC",
        "1" = "dCA1",
        "2" = "CA2-3",
        "3" = "sCA1",
        "4" = "GC",
        "5" = "dSub",
        "6" = "sSub",
        "7" = "GC",
        "8" = "mSub",
        "9" = "L6b",
        "10" = "GC",
        "11" = "MC",
        "12" = "dCA3"
)

datEx$broad_subtype_ExN <- dplyr::recode(
        as.character(datEx$seurat_clusters), !!!recode_vector
)


# Visualize fine and broad ExN subtype annotations
p1 <- DimPlot(datEx, shuffle = TRUE, label = TRUE,
              group.by = "fine_subtype_ExN") +
        coord_fixed() + theme(legend.position = "none")

p2 <- DimPlot(datEx, shuffle = TRUE, label = TRUE,
              group.by = "broad_subtype_ExN") +
        coord_fixed() + theme(legend.position = "none")

p1 + p2


### SAVE ANNOTATED EXN OBJECT ##################################################

# Save for use in subsequent ExN analysis scripts.
saveRDS(datEx, file.path(cfg$output_dir_data, "datEX_v8.rds"))


### CLEAN UP ###################################################################

# Remove temporary plotting and annotation objects.
# Retain datEx in the environment for optional interactive inspection or next steps
rm(p1, p2, recode_vector, toKeep)
