################################################################################
#
# 02_plot_megadata_v8.R
#
# Plotting and overview analyses for the annotated MegaData v8 Seurat object.
#
# This script:
#   1) Loads the annotated MegaData v8 object available from GEO
#   2) Defines plotting order and color palettes
#   3) Plots QC metrics and UMAP overviews
#   4) Summarizes broad cell-type composition
#   5) Plots broad cell-type marker expression
#   6) Generates spatial cell-type maps for all FOVs
#   7) Generates cropped representative spatial maps used for Fig. 1
#
# Input:
#   MegaData_v8_annotated.rds
#
# Author: Larissa Kraus (October 2026)
#
################################################################################

### 0) Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(patchwork)
        library(scCustomize)
        
})

### PACKAGE CONFLICTS
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)

# Memory for large objects (adjust if needed)
options(future.globals.maxSize = 100000 * 1024^2)

### 1) User-editable config -----------------------------------------------------

cfg <- list(
        # Annotated MegaData object available from GEO
        megadata_rds = NA_character_,
        
        # Plot outputs
        out_imagedim_dir = NA_character_,
        
        # Sample/FOV plotting order
        sample_order = c("L10", "L5", "M70", "P28", "P60", "P71", "P51", "M56",
                         "B1", "F1", "H1", "I1", "J1", "L20", "D11", "D8",
                         "C1", "C2", "E008", "L14", "L15", "L18",
                         "L1", "L21", "L3", "D12", "D13", "E015", "L2", "D9", "D14"),
        
        # Marker genes for broad cell types
        marker_features = c("SLC17A7", "GAD1", "AQP4", "P2RY12",
                            "MOBP", "PDGFRA", "FLT1", "FBLN1")
)

dir.create(cfg$out_imagedim_dir, recursive = TRUE, showWarnings = FALSE)

### 2) Helper functions ---------------------------------------------------------

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


### 3) Load and prepare MegaData ------------------------------------------------

dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")

# Factor ordering used throughout plotting
dat$type <- factor(dat$type, levels = c("epilepsy", "postmortem"))
dat$Sample <- factor(dat$Sample, levels = cfg$sample_order)
dat$ILAE_score <- factor(dat$ILAE_score, levels = c("P", "0", "1", "2", "3"))
dat$celltype <- factor(dat$celltype,
                       levels = c("ExN", "InN", "Astro", "Micro",
                                  "Oligo", "OPC", "Vascular", "CSF"))

# make palette for plotting per sample
dat$Sample <- factor(dat$Sample, levels = cfg$sample_order)
pal_Sample <- c("L10" = "#4477AA", "L5" = "#EE6677", "M70" = "#228833", "P28" = "#CCBB44", 
                "P60" = "#66CCEE", "P71" = "#AA3377", "P51" = "#BBBBBB", "M56" = "#332288", 
                "B1" = "#88CCEE", "F1" = "#44AA99", "H1" = "#999933", "I1" = "#882255", 
                "J1" = "#661100", "L20" = "#6699CC", "D11" = "#AA4466", "D8" = "#4477AA", 
                "C1" = "#117733", "C2" = "#DDCC77", "E008" = "#CC6677", "L14" = "#88CCAA", 
                "L15" = "#BBCC33", "L18" = "#AA4499", "L1" = "#DDDDDD", "L21" = "#77AADD", 
                "L3" = "#EE8866", "D12" = "#EEDD88", "D13" = "#FFAABB", "E015" = "#99DDFF", 
                "L2" = "#44BB99", "D9" = "#BBCCEE", "D14" = "#CC99BB")

# Cell type palette used throughout the manuscript
pal_celltype <- c(
        "ExN"      = "#E69F00",
        "InN"      = "#56B4E9",
        "Astro"    = "#009E73",
        "Micro"    = "#F0E442",
        "Oligo"    = "#0072B2",
        "OPC"      = "#D55E00",
        "Vascular" = "#CC79A7",
        "CSF"      = "#666666")

# palette for sample groups (ILAE_score)
pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00",
                    "2" = "#00C0C5", "3" = "#E69F00")

# Sample colors corresponding to each sample's ILAE score
sample_ilae <- dat@meta.data %>%
        distinct(Sample, ILAE_score)

pal_Sample_ILAE <- setNames(
        pal_ILAE_score[as.character(sample_ilae$ILAE_score[
                match(cfg$sample_order, sample_ilae$Sample)])],
        cfg$sample_order
)
### 4) QC plots -----------------------------------------------------------------

p_qc_features <- VlnPlot_scCustom(
        seurat_object = dat,
        features = "nFeature_Xenium",
        group.by = "Sample",
        pt.size = 0,
        colors_use = pal_Sample_ILAE,
        log = TRUE) +
        labs(title = paste0("nFeature_Xenium (", ncol(dat), " total cells)"),
             y = "nFeature (log)")

p_qc_counts <- VlnPlot_scCustom(
        seurat_object = dat,
        features = "nCount_Xenium",
        group.by = "Sample",
        pt.size = 0,
        colors_use = pal_Sample_ILAE,
        log = TRUE) +
        labs(title = paste0("nCount_Xenium (", ncol(dat), " total cells)"),
             y = "nCount (log)")

p_qc_features
p_qc_counts

### 5) UMAP overview ------------------------------------------------------------

Idents(dat) <- dat$seurat_clusters

p_cluster <- DimPlot(dat, shuffle = TRUE, label = TRUE) + coord_fixed()
p_sample <- DimPlot(dat, group.by = "orig.ident", shuffle = TRUE, cols = pal_Sample) + coord_fixed()
p_ilae <- DimPlot(dat, group.by = "ILAE_score", shuffle = TRUE, cols = pal_ILAE_score) + coord_fixed()
p_celltype <- DimPlot(dat, group.by = "celltype", shuffle = TRUE, cols = pal_celltype) + coord_fixed()
p_singler <- DimPlot(dat, group.by = "SingleR_pruned_midfine", shuffle = TRUE) + coord_fixed()
p_fixation <- DimPlot(dat, group.by = "FFPE_FF", shuffle = TRUE) + coord_fixed()
p_xoa <- DimPlot(dat, group.by = "XOA_version", shuffle = TRUE) + coord_fixed()
p_segmentation <- DimPlot(dat, group.by = "cell_segmentation", shuffle = TRUE) + coord_fixed()
p_type <- DimPlot(dat, group.by = "type", shuffle = TRUE) + coord_fixed()

wrap_plots(
        p_cluster, p_sample, p_ilae, p_celltype, p_singler,
        p_fixation, p_xoa, p_segmentation, p_type,
        ncol = 1
) +
        plot_annotation(title = "Dimensions: 45 | Resolution: 0.05")

### 6) Broad cell-type composition ---------------------------------------------

# Cell-type proportions within each sample
celltype_by_sample <- dat@meta.data %>%
        count(Sample, celltype, name = "cell_count") %>%
        group_by(Sample) %>%
        mutate(proportion = cell_count / sum(cell_count)) %>%
        ungroup()

celltype_by_sample$Sample <- factor(celltype_by_sample$Sample, levels = cfg$sample_order)
celltype_by_sample$celltype <- factor(
        celltype_by_sample$celltype,
        levels = levels(dat$celltype)
)

ggplot(celltype_by_sample, aes(x = Sample, y = proportion, fill = celltype)) +
        geom_col() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
        labs(title = "% Cells per Sample for each celltype",
             x = "Sample", y = "Proportion") +
        scale_fill_manual(values = pal_celltype)

# Cell-type proportions pooled within each ILAE score
celltype_by_ilae <- dat@meta.data %>%
        count(ILAE_score, celltype, name = "cell_count") %>%
        group_by(ILAE_score) %>%
        mutate(proportion = cell_count / sum(cell_count)) %>%
        ungroup()

celltype_by_ilae$ILAE_score <- factor(celltype_by_ilae$ILAE_score,
                                      levels = c("P", "0", "1", "2", "3"))
celltype_by_ilae$celltype <- factor(
        celltype_by_ilae$celltype,
        levels = levels(dat$celltype)
)

ggplot(celltype_by_ilae, aes(x = ILAE_score, y = proportion, fill = celltype)) +
        geom_col() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
        labs(title = "% Cells per ILAE score for each celltype",
             x = "ILAE score", y = "Proportion") +
        scale_fill_manual(values = pal_celltype)

### 7) Broad cell-type marker expression ---------------------------------------

DotPlot(dat, features = cfg$marker_features, group.by = "celltype") +
        scale_color_gradientn(colors = c("cyan", "white", "red")) +
        theme(
                legend.position = "top",
                legend.justification = "center",
                axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                axis.text.y = element_text(size = 12),
                axis.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 13)
        )

### 8) Spatial cell-type maps ---------------------------------------------------

# Generate one spatial cell-type map for each FOV
for (fov in cfg$sample_order) {
        pdf_file <- file.path(cfg$out_imagedim_dir, paste0(fov, "_ImageDimPlot.pdf"))
        
        pdf(pdf_file, width = 10, height = 10)
        print(
                ImageDimPlot(
                        dat,
                        fov = fov,
                        group.by = "celltype",
                        size = 1.2,
                        border.size = NA,
                        axes = FALSE,
                        dark.background = FALSE,
                        cols = pal_celltype
                ) +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank())
        )
        dev.off()
}

### 9) Representative spatial crops --------------------------------------------

# Create a copy so cropped FOVs can be added without modifying dat
dat_sub <- dat

# Coordinates correspond to representative ~1000 x 1000 regions used for Fig. 1
crop_coordinates <- list(
        I1_CA2  = list(fov = "I1",  x = c(7700, 8700), y = c(3000, 4000)),
        L15_CA2 = list(fov = "L15", x = c(4500, 5500), y = c(6700, 7700)),
        L1_CA2  = list(fov = "L1",  x = c(5000, 6000), y = c(7200, 8200)),
        M70_CA2 = list(fov = "M70", x = c(3500, 4500), y = c(4000, 5000)),
        D9_CA2  = list(fov = "D9",  x = c(2800, 3800), y = c(4000, 5000))
)

for (crop_name in names(crop_coordinates)) {
        crop_info <- crop_coordinates[[crop_name]]
        dat_sub[[crop_name]] <- Crop(
                dat_sub[[crop_info$fov]],
                x = crop_info$x,
                y = crop_info$y,
                coords = "plot"
        )
}

# Use cell segmentations for the cropped spatial maps
for (fov in names(crop_coordinates)) {
        DefaultBoundary(dat_sub[[fov]]) <- "segmentation"
}

# Plot representative cropped regions
for (fov in names(crop_coordinates)) {
        print(ImageDimPlot(dat_sub,
                           fov = fov,
                           group.by = "celltype",
                           border.color = "#A9A9A9",
                           border.size = 0.2,
                           axes = TRUE,
                           dark.background = FALSE,
                           cols = pal_celltype) +
                      theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank()))
}

                