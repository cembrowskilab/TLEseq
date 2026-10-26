
################################################################################
#
# 04_microglia_UMAP_markers_proportions.R
#
# Xenium analysis of microglia in human hippocampus samples
# Uses datM_clean_v8_annotated.rds generated in Script 01
#
# Microglia UMAPs, marker expression, and cell-state proportions
#
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
### 0) Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(ggplot2)
        library(patchwork)
})

conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename,
                             dplyr::lag, dplyr::count, dplyr::first)

### 1) Helper functions --------------------------------------------------------

# Reuse an existing Seurat object or load it from an RDS file
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

### 2) User-editable paths -----------------------------------------------------

cfg <- list(microglia_rds = NA_character_)

### 3) Color palettes ----------------------------------------------------------

# Sample colors
pal_sample <- c(
        "D9" = "#E69F00", "L15" = "#56B4E9", "L14" = "#009E73",
        "D11" = "#F0E442", "E015" = "#0072B2", "L1" = "#D55E00",
        "L3" = "#CC79A7", "L2" = "#666666", "L21" = "#AD7700",
        "L20" = "#1C91D4", "E008" = "#007756", "M70" = "#D5C711",
        "D13" = "#005685", "P28" = "#A04700", "B1" = "#B14380",
        "L18" = "#4D4D4D", "C1" = "#FFBE2D", "F1" = "#80C7EF",
        "P71" = "#00F6B3", "P51" = "#F4EB71", "C2" = "#06A5FF",
        "J1" = "#FF8320", "I1" = "#D99BBD", "D12" = "#8C8C8C",
        "L10" = "#FFCB57", "P60" = "#9AD2F2", "H1" = "#2CFFC6",
        "L5" = "#F6EF8E", "D8" = "#38B7FF", "D14" = "#FF9B4D",
        "M56" = "#E0AFCA")

# ILAE hippocampal sclerosis grades
ilae_levels <- c("P", "0", "1", "2", "3")
pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C",
                    "1" = "#7BAD00", "2" = "#00C0C5", "3" = "#E69F00")

# Microglia cell states defined in Script 01
pal_cluster_named <- c("Homeo" = "#ddcc77", "DIM" = "#117733",
                       "BAM" = "#88ccee", "DAM" = "#cc6677",
                       "Inflam.DAM" = "#aa4499", "T-cells" = "#332288",
                       "PeriNeurM" = "#44aa99")

# Marker expression color gradient
pal_marker <- c("cyan", "white", "red")

### 4) Load and prepare data ---------------------------------------------------

datM_clean <- load_if_missing("datM_clean", cfg$microglia_rds, class = "Seurat")

# Check required metadata and UMAP
required_meta <- c("cellstates", "ILAE_score", "Sample")
if (!all(required_meta %in% colnames(datM_clean@meta.data)))
        stop("Missing required metadata: ", paste(setdiff(required_meta, colnames(datM_clean@meta.data)), collapse = ", "))

if (!"umap" %in% Reductions(datM_clean))
        stop("UMAP reduction not found in datM_clean")

# Validate annotation and palette coverage
if (anyNA(datM_clean$cellstates))
        stop("Missing cell-state annotations")

if (!all(unique(as.character(datM_clean$cellstates)) %in% names(pal_cluster_named)))
        stop("Cell-state palette does not cover all annotations")

if (!all(unique(as.character(datM_clean$ILAE_score)) %in% names(pal_ILAE_score)))
        stop("ILAE palette does not cover all groups")

if (!all(unique(as.character(datM_clean$Sample)) %in% names(pal_sample)))
        stop("Sample palette does not cover all samples")

# Apply consistent ordering to metadata
datM_clean$cellstates <- factor(as.character(datM_clean$cellstates),
                                levels = names(pal_cluster_named))
datM_clean$ILAE_score <- factor(as.character(datM_clean$ILAE_score),
                                levels = ilae_levels)

### 5) UMAP plots --------------------------------------------------------------

# UMAP 1: Microglia cell-state annotations
# Shows the spatial organization of annotated states in UMAP space
umap10_0.1 <- DimPlot(datM_clean, group.by = "cellstates", shuffle = TRUE,
                      cols = pal_cluster_named, label = TRUE,
                      raster = TRUE, raster.dpi = c(1024, 1024), pt.size = 1.2) +
        coord_fixed() +
        ggtitle("Microglia") +
        xlab(expression(UMAP[1])) +
        ylab(expression(UMAP[2])) +
        theme(legend.position = "right",
              axis.line = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3))

# UMAP 2: ILAE score
# Shows which disease groups contribute cells to the UMAP
umapILAE <- DimPlot(datM_clean, group.by = "ILAE_score", shuffle = TRUE,
                    cols = pal_ILAE_score, label = FALSE,
                    raster = TRUE, pt.size = 1) +
        coord_fixed() +
        ggtitle("ILAE Score") +
        xlab(expression(UMAP[1])) +
        ylab(expression(UMAP[2])) +
        theme(legend.position = "right",
              axis.line = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3))

# UMAP 3: Sample identity
# Shows the contribution of individual samples to the UMAP
umapSAMPLE <- DimPlot(datM_clean, group.by = "Sample", shuffle = TRUE,
                      cols = pal_sample, label = FALSE,
                      raster = TRUE, pt.size = 1) +
        coord_fixed() +
        ggtitle("Sample") +
        xlab(expression(UMAP[1])) +
        ylab(expression(UMAP[2])) +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
        theme(legend.position = "right",
              legend.text = element_text(size = 7),
              legend.title = element_text(size = 9),
              legend.key.size = grid::unit(0.3, "cm"),
              legend.spacing.y = grid::unit(0.05, "cm"))

print(umap10_0.1)
print(umapILAE)
print(umapSAMPLE)

### 6) Marker gene DotPlot -----------------------------------------------------

# Selected genes for identifying and characterizing microglia cell states
DotPlot_Genes <- c("P2RY12", "CX3CR1", "CCL4", "CD83", "CD163", "LYVE1",
                   "GPNMB", "APOE", "CD68", "IFITM3", "FCGR3A", "LRRK2",
                   "TRAC", "CD2", "SNAP25", "SLC17A7", "GRIN1")

# Dot size represents the percentage of cells expressing a gene
# Dot color represents scaled average expression within each cell state
DotPlot_micro <- DotPlot(datM_clean, features = DotPlot_Genes,
                         group.by = "cellstates") +
        scale_color_gradientn(colors = pal_marker) +
        theme(legend.position = "right",
              legend.justification = "center",
              axis.text.x = element_text(size = 9, angle = 90, hjust = 1, face = "italic"),
              axis.text.y = element_text(size = 9),
              axis.title = element_text(size = 13),
              legend.text = element_text(size = 10),
              axis.line = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3),
              legend.title = element_text(size = 11)) +
        labs(x = NULL, y = "Cell state")

print(DotPlot_micro)

### 7) Cell-state proportions --------------------------------------------------

# Both analyses use pooled cell counts across all samples.
# Each cell contributes equally, so samples with more cells have more weight.
# These are descriptive cell-level proportions, not patient-level averages.

# Prepare metadata used for both plots
meta_prop <- datM_clean@meta.data %>%
        select(cellstates, ILAE_score) %>%
        filter(!is.na(cellstates), !is.na(ILAE_score))

### 7.1) Plot 1: ILAE score composition within each cell state ------------------

# Question: Of all cells assigned to a particular cell state,
# what percentage comes from each ILAE group?
#
# Numerator: Number of cells in a cell state from a specific ILAE group
# Denominator: Total number of cells in that cell state, across all groups
#
# "All" pools all microglia cell states and shows the overall ILAE
# composition of the complete microglia dataset.
#
# Each bar sums to 100%.
# This does NOT measure cell-state abundance within an ILAE group.

cellstate_levels <- names(pal_cluster_named)
cellstate_levels_all <- c(cellstate_levels, "All")

plot_data1 <- bind_rows(meta_prop,
                        mutate(meta_prop, cellstates = factor("All"))) %>%
        count(cellstates, ILAE_score, name = "n") %>%
        group_by(cellstates) %>%
        mutate(proportion = n / sum(n),
               label_text = ifelse(proportion > 0.03,
                                   paste0(round(proportion * 100, 1), "%"), "")) %>%
        ungroup() %>%
        mutate(cellstates = factor(as.character(cellstates), levels = cellstate_levels_all))

# Stacked bar graph: each bar is one cell state
# Fill color indicates the ILAE group of origin
p1_bars <- ggplot(plot_data1, aes(cellstates, proportion, fill = ILAE_score)) +
        geom_col(width = 0.95) +
        geom_text(aes(label = label_text),
                  position = position_stack(vjust = 0.5), size = 3) +
        scale_fill_manual(values = pal_ILAE_score, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        scale_y_continuous(breaks = c(0, 1), labels = c("0%", "100%"),
                           expand = c(0, 0)) +
        theme_minimal() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 13),
              plot.margin = margin(b = 2)) +
        labs(title = "ILAE score composition by cell state",
             x = NULL, y = "Proportion of cells", fill = "ILAE score")

p1_bars

### 7.2) Plot 2: Cell-state composition within each ILAE group ------------------

# Question: Within an ILAE group, what percentage of microglia
# belongs to each cell state?
#
# Numerator: Number of cells belonging to a specific cell state
#            within one ILAE group
# Denominator: Total number of microglia cells in that ILAE group
#
# Each bar sums to 100%.
# This plot describes relative cell-state abundance by ILAE group.
# It does not account for unequal patient contributions.

plot_data2 <- meta_prop %>%
        count(ILAE_score, cellstates, name = "n") %>%
        group_by(ILAE_score) %>%
        mutate(proportion = n / sum(n)) %>%
        ungroup()

# Stacked bar graph: each bar is one ILAE group
# Fill color indicates the microglia cell state
p2_bars <- ggplot(plot_data2, aes(ILAE_score, proportion, fill = cellstates)) +
        geom_col(width = 0.95) +
        scale_fill_manual(values = pal_cluster_named, drop = FALSE) +
        scale_x_discrete(drop = FALSE) +
        scale_y_continuous(breaks = c(0, 1), labels = c("0%", "100%"),
                           expand = c(0, 0)) +
        theme_minimal() +
        theme(panel.grid = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 13),
              plot.margin = margin(b = 2)) +
        labs(title = "Cell-state composition by ILAE score",
             x = NULL, y = "Proportion of cells", fill = "Cell state")
