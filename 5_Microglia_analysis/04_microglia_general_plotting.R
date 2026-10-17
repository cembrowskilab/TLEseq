################################################################################
#
# 04_microglia_UMAP_markers_proportions.R
#
# Xenium analysis of microglia in human hippocampus samples
# Using MegaData_v8
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# Microglia UMAPs, markers, and cell-state proportions


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
        library(gridExtra)
        library(dittoSeq)
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

### 3) Color palettes -------------------------------------------------------
pal_sample <- c("D9"   = "#E69F00",
                "L15"  = "#56B4E9",
                "L14"  = "#009E73",
                "D11"  = "#F0E442",
                "E015" = "#0072B2",
                "L1"   = "#D55E00",
                "L3"   = "#CC79A7",
                "L2"   = "#666666",
                "L21"  = "#AD7700",
                "L20"  = "#1C91D4",
                "E008" = "#007756",
                "M70"  = "#D5C711",
                "D13"  = "#005685",
                "P28"  = "#A04700",
                "B1"   = "#B14380",
                "L18"  = "#4D4D4D",
                "C1"   = "#FFBE2D",
                "F1"   = "#80C7EF",
                "P71"  = "#00F6B3",
                "P51"  = "#F4EB71",
                "C2"   = "#06A5FF",
                "J1"   = "#FF8320",
                "I1"   = "#D99BBD",
                "D12"  = "#8C8C8C",
                "L10"  = "#FFCB57",
                "P60"  = "#9AD2F2",
                "H1"   = "#2CFFC6",
                "L5"   = "#F6EF8E",
                "D8"   = "#38B7FF",
                "D14"  = "#FF9B4D",
                "M56"  = "#E0AFCA")


ilae_levels <- c("P", "0", "1", "2", "3")
pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00", "2" = "#00C0C5", "3" = "#E69F00")
pal_cluster_named <- c("Homeo" = "#ddcc77", "DIM" = "#117733", "BAM" = "#88ccee", 
                       "DAM" = "#cc6677", "Inflam.DAM" = "#aa4499", "T-cells" = "#332288", 
                       "PeriNeurM" = "#44aa99")
cluster_key <- data.frame(id = c("0", "1", "3", "6", "5", "4", "2"), 
                          label = names(pal_cluster_named), 
                          colour = unname(pal_cluster_named))
cluster_name_levels_noAll <- cluster_key$label
cluster_name_levels <- c(cluster_name_levels_noAll, "All")
box_colours <- c(pal_cluster_named, "All" = "white")


datM_clean <- load_if_missing("datM_clean", cfg$microglia_rds, class = "Seurat")
if (!"cluster_names" %in% colnames(datM_clean@meta.data)) stop("Annotated microglia object required (run Script 01).")


### 4) Plot Umaps --------------------------------------------------------------
umap10_0.1 <- DimPlot(datM_clean, group.by = 'cluster_names', shuffle = TRUE, cols = pal_cluster_named,
        label = TRUE, raster = TRUE, raster.dpi = c(1024, 1024), pt.size = 1.2) +
        coord_fixed() +
        ggtitle("Microglia") +
        xlab(expression(UMAP[1])) +
        ylab(expression(UMAP[2])) +
        theme(legend.position = 'right',
              axis.line = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3))

umapILAE <- DimPlot(datM_clean, group.by = 'ILAE_score', shuffle = TRUE, cols = pal_ILAE_score,
        label = FALSE, raster = TRUE, pt.size = 1) +
        coord_fixed() +
        ggtitle("ILAE Score") +
        xlab(expression(UMAP[1])) +
        ylab(expression(UMAP[2])) +
        theme(legend.position = 'right',
              axis.line = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3))

umapSAMPLE <- DimPlot(datM_clean, group.by = 'Sample', cols = pal_sample, shuffle = TRUE,
        label = FALSE, raster = TRUE, pt.size = 1) +
        coord_fixed() +
        ggtitle("Sample") +
        xlab(expression(UMAP[1])) +
        ylab(expression(UMAP[2])) +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
        theme(legend.position = 'right',
              legend.text = element_text(size = 7),
              legend.title = element_text(size = 9),
              legend.key.size = unit(0.3, "cm"),
              legend.spacing.y = unit(0.05, "cm"))

print(umap10_0.1)
print(umapILAE)
print(umapSAMPLE)


### 5) Marker dotplot -------------------------------------------------------
# Marker gene panel to display 
DotPlot_Genes <- c('P2RY12', 'CX3CR1', 'CCL4', 'CD83', 'CD163', 'LYVE1', 'GPNMB', 'APOE', 'CD68',
                   'IFITM3', 'FCGR3A', 'LRRK2', 'TRAC', 'CD2', 'SNAP25', 'SLC17A7', 'GRIN1')

# Plot, grouped by cluster name
DotPlot_micro <- DotPlot(datM_clean, features = DotPlot_Genes, group.by = 'cluster_names') +
        scale_color_gradientn(colors = c("cyan", "white", "red")) +
        theme(legend.position = "right",
              legend.justification = "center",
              axis.text.x = element_text(size = 9, angle = 90, hjust = 1, face = "italic"),
              axis.text.y = element_text(size = 9),
              axis.title = element_text(size = 13),
              legend.text = element_text(size = 10),
              axis.line = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3),
              legend.title = element_text(size = 11)) +
        labs(x = NULL, y = "Cluster")

DotPlot_micro


### 6) Cell-state proportions -----------------------------------------------
# Calculate proportions within each cluster (plot 1, gets an extra "All" group) and within each ILAE score (plot 2)
# ILAE score composition within each cluster
plot_data1 <- bind_rows(datM_clean@meta.data, mutate(datM_clean@meta.data, cluster_names = "All")) %>%
        dplyr::count(cluster_names, ILAE_score) %>%
        group_by(cluster_names) %>%
        mutate(proportion = n / sum(n),
               label_text = ifelse(proportion > 0.03, paste0(round(proportion * 100, 1), "%"), "")) %>%
        ungroup() %>%
        mutate(cluster_names = factor(cluster_names, levels = cluster_name_levels))
# cluster composition within each ILAE score
plot_data2 <- datM_clean@meta.data %>%
        dplyr::count(ILAE_score, cluster_names) %>%
        group_by(ILAE_score) %>%
        mutate(proportion = n / sum(n)) %>%
        ungroup() %>%
        mutate(ILAE_score = factor(ILAE_score, levels = ilae_levels))

# Plot
p1_bars <- ggplot(plot_data1, aes(cluster_names, proportion, fill = ILAE_score)) +
        geom_col(width = 0.95) +
        geom_text(aes(label = label_text), position = position_stack(vjust = 0.5), size = 3) +
        scale_fill_manual(values = pal_ILAE_score) +
        scale_y_continuous(breaks = 0:1, labels = c("0%", "100%"), expand = c(0, 0)) + 
        theme_minimal() +
        theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 13), plot.margin = margin(b = 2)) +
        labs(title = "ILAE Score Composition by Cluster", x = NULL, y = "Proportion of cells", fill = "ILAE Score")

p1_strip <- ggplot(data.frame(x = factor(cluster_name_levels, cluster_name_levels)), aes(x, 1, fill = x)) +
        geom_tile(aes(colour = x == "All"), linewidth = 0.8, width = 0.95, height = 0.9) +
        geom_text(aes(label = x), fontface = "bold", size = 3.5) +
        scale_fill_manual(values = box_colours, guide = "none") +
        scale_colour_manual(values = c("TRUE" = "black", "FALSE" = NA), guide = "none") +
        labs(x = "Cluster") +
        theme_void() +
        theme(axis.title.x = element_text(size = 12, margin = margin(t = 4)), plot.margin = margin(t = 2))

p1 <- p1_bars / p1_strip + plot_layout(heights = c(10, 1.5))


p2_bars <- ggplot(plot_data2, aes(ILAE_score, proportion, fill = cluster_names)) +
        geom_col(width = 0.95) +
        scale_fill_manual(values = pal_cluster_named) +
        scale_y_continuous(breaks = 0:1, labels = c("0%", "100%"), expand = c(0, 0)) + 
        theme_minimal() +
        theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 13), plot.margin = margin(b = 2)) +
        labs(title = "Cluster Composition by ILAE Score", x = NULL, y = "Proportion of cells", fill = "Cluster")

p2_strip <- ggplot(data.frame(x = factor(ilae_levels, ilae_levels)), aes(x, 1, fill = x)) +
        geom_tile(width = 0.95, height = 0.9) +
        geom_text(aes(label = x), fontface = "bold", size = 4) +
        scale_fill_manual(values = pal_ILAE_score, guide = "none") +
        labs(x = "ILAE Score") + 
        theme_void() +
        theme(axis.title.x = element_text(size = 12, margin = margin(t = 4)), plot.margin = margin(t = 2))

p2 <- p2_bars / p2_strip + plot_layout(heights = c(10, 1))

print(p1)
print(p2)
