################################################################################
#
# 02_microglia_HuMicA_validation.R
#
# Xenium analysis of microglia in human hippocampus samples
# Using MegaData_v8
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# Microglia HuMicA reference validation


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
        library(readxl)
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


datM_clean <- load_if_missing("datM_clean", cfg$microglia_rds, class = "Seurat")
if (!"cluster_names" %in% colnames(datM_clean@meta.data)) stop("Annotated microglia object required (run Script 01).")
markers <- read.csv(cfg$markers_csv)

### 4) Humica reference comparison ------------------------------------------

# Load and clean the HuMicA reference marker table (Supplementary Data 2; Martins-Ferreira et al, 2025, NatComm)
ref <- readxl::read_excel(cfg$humica_xlsx, skip = 1)
# Some numeric columns were exported with commas instead of decimal points — fix those
ref <- ref %>% mutate(avg_log2FC = as.numeric(gsub(",", ".", avg_log2FC)),
                      p_val_adj = as.numeric(gsub(",", ".", p_val_adj)))

# Add cluster names to our own marker table
markers <- markers %>%
  mutate(cluster_name = cluster_key$label[match(as.character(cluster), cluster_key$id)])

# Rank marker genes within each cluster by a combined score that rewards genes that are
#    both strongly upregulated AND statistically robust. A gene's importance-weighted overlap
#    with the HuMicA reference is high only if it ranks highly in BOTH marker lists.

#(p_val_adj of exactly 0 is floored to avoid -log10(0) = Inf.)
p_floor <- 1e-300

# Rank our marker genes within each cluster
rank_weighted <- markers %>%
        mutate(combined_score = -log10(pmax(p_val_adj, p_floor)) * avg_log2FC) %>%
        group_by(cluster_name) %>%
        mutate(ours_rank = rank(-combined_score)) %>%
        ungroup() %>%
        select(Gene = gene, ours_cluster = cluster_name, ours_rank)

# Rank the HuMicA reference marker genes within each of their clusters
ref_rank <- ref %>%
        filter(!is.na(cluster)) %>%
        mutate(combined_score = -log10(pmax(p_val_adj, p_floor)) * avg_log2FC) %>%
        group_by(cluster) %>%
        mutate(humica_rank = rank(-combined_score)) %>%
        ungroup() %>%
        select(Gene, humica_cluster = cluster, humica_rank)

# Join on shared genes and compute a per-cluster-pair overlap score. Top-ranked shared
#    genes (low rank in both lists) contribute the most to the score.
rank_overlap <- inner_join(rank_weighted, ref_rank, by = "Gene") %>%
        mutate(rank_score = 1 / (ours_rank * humica_rank)) %>%
        group_by(ours_cluster, humica_cluster) %>%
        summarise(weighted_score = sum(rank_score), n_shared = n(), .groups = "drop")

# Heatmap of cluster-to-cluster overlap, with the number of shared genes overlaid as text
RankPlot <- ggplot(rank_overlap,
                   aes(x = factor(ours_cluster), y = factor(humica_cluster), fill = weighted_score)) +
        geom_tile() +
        geom_text(aes(label = n_shared,
                     color = weighted_score > 0.5 * max(rank_overlap$weighted_score, na.rm = TRUE)),
                 size = 3, show.legend = FALSE) +
        scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
        scale_fill_viridis_c(option = "viridis", direction = -1, trans = "sqrt",
                             limits = c(0, max(rank_overlap$weighted_score, na.rm = TRUE))) +
        theme_test() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Present Study Cluster", y = "HuMicA Cluster")

ggsave(file.path(cfg$output_dir_plots, "HuMicAMarkerPlot.pdf"), plot = RankPlot, width = 7, height = 5)
