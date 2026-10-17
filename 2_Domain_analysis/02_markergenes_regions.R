################################################################################
# 02_markergenes_regions.R
# Marker gene analysis for PRECAST-defined hippocampal domains
# Input: 07_seumerged.rds created by 01_domainIdentification.R (Section 8)
#
# Original author: Mathias Delhaye, October 2026
# Code cleanup and standardization: Larissa Kraus (October 2026)

### 0) Setup -------------------------------------------------------------------
suppressPackageStartupMessages({
        library(conflicted)
        library(tidyverse)
        library(Seurat)
})

### PACKAGE CONFLICTS ---------------------------------------------------------
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)


### 1) User-editable config ----------------------------------------------------
cfg <- list(
        # Output of 01_domainIdentification.R, Section 8 (cfg$rds_dir/07_seumerged.rds).
        # Only used if 'seumerged' is not already loaded in the R environment.
        seumerged_file = NA_character_,
        
        # Folder for marker-gene plots and session information.
        fig_dir = NA_character_,
        
        # Reuse an existing Seurat object named 'seumerged' when available.
        # Set FALSE to force loading the RDS file instead.
        
        # Marker-detection thresholds
        min_log2fc = 1,
        top_n_markers = 10,
        
        # Curated marker panel for region dot plot
        marker_features_dotplot = c("GRM2", "NTNG2", "KIT", "WIF1", "GRM4",
                                    "GAD2", "RELN", "MAG", "TGFBI", "TP73"),
        
        # Region order and palette (matching script 01)
        region_levels = c("GCL", "ML", "CA2_CA4", "CA1", "SUB", "GABA",
                          "SL.SR.SLM", "WM", "vascular", "pia", "void"),
        col_regions = setNames(
                c("#3D0F99", "#C7B8E6", "#0F8299", "#B8DEE6", "#54990F", "#CFE6B8",
                  "#333333", "#990F26", "#CCAA7A", "#99600F", "#FFD300"),
                c("GCL", "ML", "CA2_CA4", "CA1", "SL.SR.SLM", "GABA", "SUB", "WM", "pia", "vascular", "void")
        ),
        f_fig_dotplot = "marker_genes_dotplot.png"
)

# All output directories are created automatically.
marker_fig_dir <- file.path(cfg$fig_dir, "04_marker_genes_regions")
dir.create(marker_fig_dir, recursive = TRUE, showWarnings = FALSE)

### 2) Helper functions --------------------------------------------------------

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

plot_marker_violin <- function(seu, gene, group_by, palette, out_file) {
        png(out_file)
        tryCatch(print(VlnPlot(seu, features = gene, group.by = group_by,
                               pt.size = 0, cols = palette)),
                 finally = dev.off())
}

### 3) Load input from script 01, unless already in memory --------------------
# Script 01 creates 'seumerged' after PRECAST clustering and manual region edits.
# It saves this object as 07_seumerged.rds; no MegaData load is required here.
seumerged <- load_if_missing("seumerged", cfg$seumerged_file, class = "Seurat")
stopifnot(inherits(seumerged, "Seurat"))
if (!all(c("region", "region_md") %in% colnames(seumerged@meta.data))) {
        stop("'seumerged' must contain 'region' and 'region_md' metadata from script 01.")
}
seumerged$region <- factor(seumerged$region, levels = cfg$region_levels)
seumerged$region_md <- factor(seumerged$region_md, levels = cfg$region_levels)

### 4) Prep seumerged ----------------------------------------------------------
seumerged_joined <- JoinLayers(seumerged)
seumerged_joined <- subset(seumerged_joined, subset = !is.na(region) & region != "void")

### 5) Marker genes per region -------------------------------------------------
message("Finding region marker genes (FindAllMarkers, group.by = 'region')")
markers_regions <- FindAllMarkers(seumerged_joined, group.by = "region", only.pos = TRUE)

# Rank by fold change before taking the top N markers per region.
top10_regions <- markers_regions %>%
        group_by(cluster) %>%
        filter(avg_log2FC > cfg$min_log2fc) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        slice_head(n = cfg$top_n_markers) %>%
        ungroup()

top10_regions_markers_only <- top10_regions %>%
        select(cluster, gene) %>%
        rename(regions = cluster)

### 6) Per-gene violin plots ---------------------------------------------------
for (region in unique(top10_regions_markers_only$regions)) {
        out_dir <- file.path(marker_fig_dir, as.character(region))
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        markers <- top10_regions %>% filter(cluster == region) %>% pull(gene)
        for (g in markers) {
                plot_marker_violin(seu = seumerged_joined, gene = g, group_by = "region",
                                   palette = cfg$col_regions,
                                   out_file = file.path(out_dir, paste0(g, "_vln_regions.png")))
        }
}
message("Saved per-gene violin plots to: ", marker_fig_dir)

### 7) Curated marker panel dot plot -------------------------------------------
missing_features <- setdiff(cfg$marker_features_dotplot, rownames(seumerged_joined))
if (length(missing_features)) warning("Missing dot-plot genes: ", paste(missing_features, collapse = ", "))
features_to_plot <- intersect(cfg$marker_features_dotplot, rownames(seumerged_joined))
if (!length(features_to_plot)) stop("None of the requested dot-plot genes are in seumerged.")

p <- DotPlot(seumerged_joined, features = features_to_plot, group.by = "region") +
        scale_color_gradientn(colors = c("cyan", "white", "red")) +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              legend.position = "top", legend.justification = "right",
              legend.text = element_text(size = 9), legend.title = element_text(size = 9))

ggsave(filename = file.path(marker_fig_dir, cfg$f_fig_dotplot), plot = p,
       dpi = 300, width = 8, height = 6)
message("Saved marker dot plot to: ", file.path(marker_fig_dir, cfg$f_fig_dotplot))

