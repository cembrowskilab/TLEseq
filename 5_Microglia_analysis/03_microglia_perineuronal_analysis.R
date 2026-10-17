################################################################################
#
# 03_microglia_perineuronal_analysis.R
#
# Xenium analysis of microglia in human hippocampus samples
# Using MegaData_v8
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# Perineuronal microglia analysis


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
        library(FNN)
})

### PACKAGE CONFLICTS
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag, dplyr::count, dplyr::first)


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

# Round axis limits to the next 0.5.
ceiling_to_half <- function(x) ceiling(x * 2) / 2

### 2) User-editable paths and load inputs -------------------------------------
# Set input/output paths before running. Objects already in the session are reused.
cfg <- list(
        megadata_rds = NA_character_,
        region_metadata_rds = NA_character_,
        microglia_rds = NA_character_,
        markers_csv = NA_character_,
        humica_xlsx = NA_character_,
        output_dir_data = NA_character_)
if (is.na(cfg$output_dir_data)) stop("Set cfg$output_dir_data")
dir.create(cfg$output_dir_data, recursive = TRUE, showWarnings = FALSE)
datM_clean <- load_if_missing("datM_clean", cfg$microglia_rds, class = "Seurat")
if (!"cluster_names" %in% colnames(datM_clean@meta.data)) stop("Annotated microglia object required (run Script 01).")
dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")
if (!"cellsubtype" %in% colnames(dat@meta.data)) stop("Annotated cellsubtypes required (run Script 01_make_megasata_v8).")

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



### 4) nCount and nfeature QC --------------------------------------------------

# QC violin plot (transcript / feature counts per cell) across clusters
PeriNeurM_featplot <- VlnPlot(datM_clean, features = c("nCount_Xenium", "nFeature_Xenium"),
                              group.by = "cluster_names", pt.size = 0, cols = pal_cluster_named)
PeriNeurM_featplot
### 5) Distance to Nearest Neuron --------------------------------------------

# Pull cell coordinates + cell type for one representative FOV (sample D9), for a
#    nearest-neuron distance analysis of the microglia in that FOV
coords_full <- GetTissueCoordinates(dat, image = "D9")
coords_full$celltype <- dat$celltype[match(coords_full$cell, colnames(dat))]

neuron_coords <- coords_full %>% filter(celltype %in% c("ExN", "InN"))
micro_coords  <- coords_full %>% filter(celltype == "Micro")

# Nearest-neuron distance for every microglia cell in this FOV
nn <- get.knnx(data = neuron_coords[, c("x", "y")], query = micro_coords[, c("x", "y")], k = 1)
micro_coords$dist_to_nearest_neuron <- nn$nn.dist[, 1]

# Merge in each microglia cell's cluster identity
micro_coords$cluster_names <- datM_clean$cluster_names[match(micro_coords$cell, colnames(datM_clean))]

# Compare distance-to-nearest-neuron across clusters.
#    Note: cluster_names = NA corresponds to cells that were labeled "Micro" in dat$celltype
#    but were filtered out of datM_clean earlier (i.e. did not pass subsetting
#    in "01_microglia_subsetting_clustering_annotation.R").
DNN_violin_D9 <- ggplot(micro_coords, aes(x = cluster_names, y = dist_to_nearest_neuron, fill = cluster_names)) +
        geom_violin() +
        scale_fill_manual(values = pal_cluster_named) +
        theme_classic() +
        labs(title = "Distance to nearest neuron by microglia cluster (D9)",
            y = "Distance to Nearest Neuron (coordinate units)") +
        theme(legend.position = "none")
DNN_violin_D9


# Build the D9-only microglia object and attach each cell's distance to its nearest neuron.
#    intersect() guards against D9 "Micro" cells that were filtered out of datM_clean earlier
#    (see note above) so subset() doesn't fail on cell names that aren't in the object.
D9_cells <- intersect(micro_coords$cell, colnames(datM_clean))
datM_clean_D9 <- subset(datM_clean, cells = D9_cells)
datM_clean_D9$dist_to_nearest_neuron <- micro_coords$dist_to_nearest_neuron[
        match(colnames(datM_clean_D9), micro_coords$cell)]

# Round a value up to the nearest 0.5, used for setting y-axis maxima below
# SNAP25 vs. distance-to-neuron scatter
snap25_vals <- FetchData(datM_clean_D9, "SNAP25")[, 1]
dist_vals   <- datM_clean_D9$dist_to_nearest_neuron
snap25_max  <- ceiling_to_half(max(snap25_vals))
snap25_r    <- cor(dist_vals, snap25_vals)

scatter_snap25 <- FeatureScatter(datM_clean_D9, feature1 = "dist_to_nearest_neuron", feature2 = "SNAP25",
                                 group.by = "cluster_names", cols = pal_cluster_named) +
        scale_y_continuous(limits = c(0, snap25_max)) +
        ggtitle(bquote(italic(r) == .(round(snap25_r, 2)))) +
        ylab(bquote(italic("SNAP25"))) +
        xlab("Distance to Nearest Neuron (coordinate units)")

# GRIN1 vs. distance-to-neuron scatter
grin1_vals <- FetchData(datM_clean_D9, "GRIN1")[, 1]
grin1_max  <- ceiling_to_half(max(grin1_vals))
grin1_r    <- cor(dist_vals, grin1_vals)

scatter_grin1 <- FeatureScatter(datM_clean_D9, feature1 = "dist_to_nearest_neuron", feature2 = "GRIN1",
                                group.by = "cluster_names", cols = pal_cluster_named) +
        scale_y_continuous(limits = c(0, grin1_max)) +
        ggtitle(bquote(italic(r) == .(round(grin1_r, 2)))) +
        ylab(bquote(italic("GRIN1"))) +
        xlab("Distance to Nearest Neuron (coordinate units)")

# Combine both scatter plots onto one page and save
combined_scatter <- scatter_snap25 | scatter_grin1
combined_scatter
### 6) Cropped Spatialplots with selected features -----------------------------

# Make a metadata column that makes it easier to find PeriNeurM without changing 
# original MegaData
dat1<-dat
dat1$cellsubtype_clean <- as.character(dat1$cellsubtype)

# Collapse all InN subtypes to "InN"
dat1$cellsubtype_clean[grepl("^InN", dat1$cellsubtype_clean)] <- "InN"

# Collapse all ExN subtypes to "ExN"
dat1$cellsubtype_clean[grepl("^ExN", dat1$cellsubtype_clean)] <- "ExN"

# Collapse all Micro_ subtypes to "Micro_other", EXCEPT PeriNeurM
dat1$cellsubtype_clean[grepl("^Micro", dat1$cellsubtype_clean) &
                               dat1$cellsubtype_clean != "Micro_PeriNeurM"] <- "Micro/T"

# Check result
unique(dat1$cellsubtype_clean)
unique(dat1$cellsubtype)

#make a new palette for plotting
pal_cellsubtype_clean <- c(
        # Excitatory neurons
        "ExN"         = "#E69F00",
        "InN"         = "#56B4E9",
        "Micro/T"       = "#332288",
        "Micro_PeriNeurM"   = "#44AA99",
        "Vascular"    = "#F0E442",
        "OPC"         = "#0072B2",
        "Astro"       = "#D55E00",
        "Oligo"       = "#CC79A7",
        "CSF"         = "#666666"
)


###Crop to CA2 for representative sections
#plot the sample you want to look at, Scale helps you find the coordinates for crop function
ImageDimPlot(dat1, fov = 'I1', group.by = 'cellsubtype_clean',size = 1,  border.size = NA,
             axes = TRUE, dark.background = FALSE, cols=pal_cellsubtype_clean)+
        scale_x_continuous(breaks = seq(0, 12000, by = 1000)) +
        scale_y_continuous(breaks = seq(0, 12000, by = 1000))

#make a new fov/image only for the cropped view: crop CA2 small piece
#important to not have '_' in new fov name, otherwise molecules don't plot
dat1[["D9.CA2"]] <- Crop(dat1[["D9"]],
                         x = c(1300, 1800),
                         y = c(4400, 4800),
                         coords = "plot")

dat1[["L15.CA2"]] <- Crop(dat1[["L15"]],
                          x = c(6000, 6400), 
                          y = c(2800, 3200),
                          coords = "plot")

dat1[["L1.CA2"]] <- Crop(dat1[["L1"]],
                         x = c(5000, 5400),
                         y = c(8200, 8600),
                         coords = "plot")

dat1[["M70.CA2"]] <- Crop(dat1[["M70"]],
                          x = c(6200, 6600),
                          y = c(1600, 2000),
                          coords = "plot")

dat1[["I1.CA2"]] <- Crop(dat1[["I1"]],
                         x = c(8000, 8400),
                         y = c(1500, 1900),
                         coords = "plot")

# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(dat1[["I1.CA2"]]) <- "segmentation"
ImageDimPlot(dat1, fov = 'I1.CA2', group.by = 'cellsubtype_clean', border.color = "#A9A9A9", 
             border.size = 0.2,axes = TRUE, dark.background = FALSE, cols = pal_cellsubtype_clean,
             molecules = c("SLC17A7","SNAP25", "GRIN1", 'P2RY12'),nmols = 1000, mols.size = 0.5) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

DefaultBoundary(dat1[["D9.CA2"]]) <- "segmentation"
ImageDimPlot(dat1, fov = 'D9.CA2', group.by = 'cellsubtype_clean', border.color = "#A9A9A9", 
             border.size = 0.2,axes = TRUE, dark.background = FALSE, cols = pal_cellsubtype_clean,
             molecules = c("SLC17A7","SNAP25", "GRIN1", 'P2RY12'),nmols = 1000, mols.size = 0.5) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

DefaultBoundary(dat1[["L15.CA2"]]) <- "segmentation"
ImageDimPlot(dat1, fov = 'L15.CA2', group.by = 'cellsubtype_clean', border.color = "#A9A9A9", 
             border.size = 0.2,axes = TRUE, dark.background = FALSE, cols = pal_cellsubtype_clean,
             molecules = c("SLC17A7","SNAP25", "GRIN1", 'P2RY12'),nmols = 1000, mols.size = 0.5) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

DefaultBoundary(dat1[["L1.CA2"]]) <- "segmentation"
ImageDimPlot(dat1, fov = 'L1.CA2', group.by = 'cellsubtype_clean', border.color = "#A9A9A9", 
             border.size = 0.2,axes = TRUE, dark.background = FALSE, cols = pal_cellsubtype_clean,
             molecules = c("SLC17A7","SNAP25", "GRIN1", 'P2RY12'),nmols = 1000, mols.size = 0.5) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

DefaultBoundary(dat1[["M70.CA2"]]) <- "segmentation"
ImageDimPlot(dat1, fov = 'M70.CA2', group.by = 'cellsubtype_clean', border.color = "#A9A9A9", 
             border.size = 0.2,axes = TRUE, dark.background = FALSE, cols = pal_cellsubtype_clean,
             molecules = c("SLC17A7","SNAP25", "GRIN1", 'P2RY12'),nmols = 1000, mols.size = 0.5) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




