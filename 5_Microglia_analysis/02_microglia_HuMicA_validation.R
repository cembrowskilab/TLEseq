
################################################################################
#
# 02_microglia_HuMicA_validation.R
#
# Xenium analysis of microglia in human hippocampus samples
# Uses datM_clean_v8_annotated.rds generated in Script 01
#
# Reference: Supplementary Data 2, Martins-Ferreira et al., 2025,
# Nature Communications (HuMicA dataset)
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
        library(readxl)
})

conflicted::conflicts_prefer(dplyr::filter, 
                             dplyr::select, 
                             dplyr::rename,
                             dplyr::lag, 
                             dplyr::count, 
                             dplyr::first)

### 1) Helper functions --------------------------------------------------------

# Reuse an existing object in the global environment; otherwise load its RDS file
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
NA_character_
cfg <- list(
        microglia_rds = NA_character_,
        markers_csv =  'Z:/CembrowskiLab/Larissa/1. Projects/Human brain tissue/3_Xenium/1. Human_samples_analysis/Analysis/github_micro_test/micro_marker.csv',
        humica_xlsx =  'Z:/CembrowskiLab/Larissa/1. Projects/Human brain tissue/3_Xenium/1. Human_samples_analysis/Analysis/github_micro_test/41467_2025_56124_MOESM4_ESM.xlsx')
        
### 3) Load data ---------------------------------------------------------------

datM_clean <- load_if_missing("datM_clean", cfg$microglia_rds, class = "Seurat")
markers <- read.csv(cfg$markers_csv)
ref <- readxl::read_excel(cfg$humica_xlsx, skip = 1)

if (!all(c("seurat_clusters", "cellstates") %in% colnames(datM_clean@meta.data)))
        stop("Microglia object requires seurat_clusters and cellstates metadata from Script 01_microglia")

if (!all(c("cluster", "gene", "avg_log2FC", "p_val_adj") %in% colnames(markers)))
        stop("Marker table is missing required columns")

if (!all(c("cluster", "Gene", "avg_log2FC", "p_val_adj") %in% colnames(ref)))
        stop("HuMicA reference is missing required columns")

### 4) Prepare marker tables ---------------------------------------------------

# Convert decimal commas in HuMicA reference
ref <- ref %>%
        mutate(avg_log2FC = as.numeric(gsub(",", ".", avg_log2FC)),
               p_val_adj = as.numeric(gsub(",", ".", p_val_adj)))

# Extract cluster-to-cell-state annotations directly from Seurat metadata
cluster_annotations <- datM_clean@meta.data %>%
        distinct(seurat_clusters, cellstates) %>%
        mutate(seurat_clusters = as.character(seurat_clusters),
               cellstates = as.character(cellstates))

# Each Seurat cluster must correspond to exactly one cell state
if (anyDuplicated(cluster_annotations$seurat_clusters))
        stop("Multiple cellstates annotations found within a Seurat cluster")

if (anyNA(cluster_annotations$cellstates))
        stop("Missing cellstates annotations")

# Add cell-state names to marker table
markers <- markers %>%
        mutate(cluster = as.character(cluster)) %>%
        left_join(cluster_annotations, by = c("cluster" = "seurat_clusters")) %>%
        rename(cluster_name = cellstates)

if (anyNA(markers$cluster_name))
        stop("Some marker clusters could not be matched to Seurat metadata")

### 5) HuMicA reference comparison --------------------------------------------

# Combined score rewards strong positive fold changes and statistical significance
# Adjusted p-values of zero are floored to avoid infinite scores
p_floor <- 1e-300

# Rank marker genes within each microglia cell state
rank_weighted <- markers %>%
        mutate(combined_score = -log10(pmax(p_val_adj, p_floor)) * avg_log2FC) %>%
        group_by(cluster_name) %>%
        mutate(ours_rank = rank(-combined_score)) %>%
        ungroup() %>%
        select(Gene = gene, ours_cluster = cluster_name, ours_rank)

# Rank HuMicA reference marker genes within each cluster
ref_rank <- ref %>%
        filter(!is.na(cluster)) %>%
        mutate(combined_score = -log10(pmax(p_val_adj, p_floor)) * avg_log2FC) %>%
        group_by(cluster) %>%
        mutate(humica_rank = rank(-combined_score)) %>%
        ungroup() %>%
        select(Gene, humica_cluster = cluster, humica_rank)

# Calculate weighted marker overlap between cluster pairs
rank_overlap <- inner_join(rank_weighted, ref_rank, by = "Gene") %>%
        mutate(rank_score = 1 / (ours_rank * humica_rank)) %>%
        group_by(ours_cluster, humica_cluster) %>%
        summarise(weighted_score = sum(rank_score),
                  n_shared = n(), .groups = "drop")

if (nrow(rank_overlap) == 0)
        stop("No shared marker genes found between Xenium and HuMicA")

### 6) Plot weighted marker overlap -------------------------------------------

RankPlot <- ggplot(rank_overlap,
                   aes(x = factor(ours_cluster), y = factor(humica_cluster),
                       fill = weighted_score)) +
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

print(RankPlot)
