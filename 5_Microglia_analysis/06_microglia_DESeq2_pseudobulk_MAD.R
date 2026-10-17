################################################################################
#
# 06_microglia_DESeq2_pseudobulk_MAD.R
#
# Xenium analysis of microglia in human hippocampus samples
# Using MegaData_v8
# Original analysis: Larissa Kraus (2026)
#
# Patient-level microglia pseudobulk DESeq2 and MAD filtering

# DESeq2 approach adapted from the ExN pseudobulk workflow.
# MAD filtering approach adapted from Adi Swaro (July 2025).


### 0) Setup -------------------------------------------------------------------

# Core packages
suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(purrr)
        library(readr)
        library(ggplot2)
        library(EnhancedVolcano)
        library(DESeq2)
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

# Build all available pairwise contrasts in the same order as the ExN DEG workflow.
generate_ilae_pairs <- function(groups, order = c("P", "0", "1", "2", "3")) {
        groups <- order[order %in% as.character(groups)]
        if (length(groups) < 2) stop("Fewer than two ILAE groups present.")
        combn(groups, 2, simplify = FALSE)
}

# MAD = stats::mad() default (normal-consistency factor 1.4826).
# Per-gene threshold = median(mean raw counts per cell across full-data clusters) + 2*MAD.
# The full-data clusters are used ONLY as a background reference for this filter.
calculate_mad_thresholds <- function(full_obj, micro_obj, z = 2, assay = "Xenium") {
        if (!"seurat_clusters" %in% colnames(full_obj@meta.data)) stop("Full MegaData needs seurat_clusters.")
        full_obj <- JoinLayers(full_obj, assay = assay)
        micro_obj <- JoinLayers(micro_obj, assay = assay)
        full_counts <- GetAssayData(full_obj, assay = assay, layer = "counts")
        micro_counts <- GetAssayData(micro_obj, assay = assay, layer = "counts")
        genes <- intersect(rownames(full_counts), rownames(micro_counts))
        if (!length(genes)) stop("No common genes in full and microglia Xenium counts.")
        clusters <- as.character(full_obj$seurat_clusters)
        if (anyNA(clusters)) stop("Missing full-data cluster labels.")
        cluster_means <- vapply(unique(clusters), function(cl) {
                Matrix::rowMeans(full_counts[genes, which(clusters == cl), drop = FALSE])
        }, numeric(length(genes)))
        if (is.null(dim(cluster_means))) cluster_means <- matrix(cluster_means, ncol = 1)
        med <- apply(cluster_means, 1, median, na.rm = TRUE)
        mad_value <- apply(cluster_means, 1, stats::mad, na.rm = TRUE)
        micro_mean <- Matrix::rowMeans(micro_counts[genes, , drop = FALSE])
        tibble(gene = genes, cluster_median = unname(med), cluster_MAD = unname(mad_value),
               MAD_threshold = unname(med + z * mad_value), microglia_mean = unname(micro_mean),
               passes_MAD = microglia_mean > med + z * mad_value)
}

### 2) User-editable paths and load inputs -------------------------------------

# Set input/output paths before running. Objects already in the session are reused.
cfg <- list(
        microglia_rds = NA_character_,
        megadata_rds = NA_character_,
        output_dir_de = NA_character_)

if (is.na(cfg$output_dir_de)) stop("Set cfg$output_dir_de")

dir.create(cfg$output_dir_de, recursive = TRUE, showWarnings = FALSE)

raw_dir <- file.path(cfg$output_dir_de, "DESeq2_unfiltered")

filtered_dir <- file.path(cfg$output_dir_de, "DESeq2_MAD_filtered")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

### 3) Color palette --------------------------------------------------------

pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00", "2" = "#00C0C5", "3" = "#E69F00")

### 4) Load microglia and aggregate by patient ------------------------------

# Same microglia subset as the clustering workflow: full-data cluster 2, restricted
# to SingleR_label_midfine == "Micro/Macro/T". No state or subcluster restriction.
datM_clean <- load_if_missing("datM_clean", cfg$microglia_rds, class = "Seurat")

required <- c("orig.ident", "ILAE_score")

if (!all(required %in% colnames(datM_clean@meta.data))) stop("Missing metadata: ", paste(setdiff(required, colnames(datM_clean@meta.data)), collapse = ", "))

if (anyNA(datM_clean$orig.ident) || anyNA(datM_clean$ILAE_score)) stop("Missing patient or ILAE labels.")

patient_groups <- datM_clean@meta.data %>% distinct(orig.ident, ILAE_score)

if (anyDuplicated(patient_groups$orig.ident)) stop("Some patients have multiple ILAE scores.")

# Aggregate raw counts once per patient. ILAE is patient-level metadata, not an
# additional biological replicate; one pseudobulk column represents one patient.
pseudo_datM <- AggregateExpression(datM_clean, assays = "Xenium", return.seurat = TRUE,
                                   group.by = c("ILAE_score", "orig.ident"))

DefaultAssay(pseudo_datM) <- "Xenium"

pseudo_datM$ILAE_score <- factor(sub("^g(?=[0-9]+$)", "", as.character(pseudo_datM$ILAE_score), perl = TRUE),
                                  levels = c("P", "0", "1", "2", "3"))

if (anyNA(pseudo_datM$ILAE_score)) stop("Could not recover ILAE group metadata in pseudobulk object.")

Idents(pseudo_datM) <- "ILAE_score"

print(table(pseudo_datM$ILAE_score, useNA = "ifany"))

### 5) Pairwise deseq2 comparisons ------------------------------------------

# Match the ExN pipeline: Seurat::FindMarkers(test.use = "DESeq2") on patient
# pseudobulk profiles; positive avg_log2FC indicates higher expression in group 1.
# Adjusted p-values are from the DESeq2/Seurat output for each comparison.
all_results <- list()

failed_comparisons <- list()

for (pair in generate_ilae_pairs(pseudo_datM$ILAE_score)) {
        group1 <- pair[1]; group2 <- pair[2]
        contrast <- paste0(group1, "_vs_", group2)
        message("DESeq2: ", group1, " - ", group2)
        markers <- tryCatch(FindMarkers(pseudo_datM, ident.1 = group1, ident.2 = group2,
                                        test.use = "DESeq2", assay = "Xenium"),
                            error = function(e) {message("DESeq2 failed: ", conditionMessage(e)); NULL})
        if (is.null(markers)) {
                failed_comparisons[[contrast]] <- tibble(comparison = contrast, status = "DESeq2 error")
                next
        }
        markers$gene <- rownames(markers)
        markers$group_1 <- group1
        markers$group_2 <- group2
        markers$comparison <- contrast
        readr::write_csv(as_tibble(markers), file.path(raw_dir, paste0("Markergenes_", contrast, ".csv")))
        all_results[[contrast]] <- as_tibble(markers)
        if (!nrow(markers)) next
        plot_data <- markers %>% mutate(p_val_adj_plot = pmax(p_val_adj, .Machine$double.xmin))
        if (any(is.finite(plot_data$avg_log2FC) & is.finite(plot_data$p_val_adj_plot))) {
                p <- EnhancedVolcano::EnhancedVolcano(plot_data, lab = plot_data$gene,
                        x = "avg_log2FC", y = "p_val_adj_plot", pCutoff = 0.05,
                        title = paste("Microglia:", group1, "-", group2))
                ggsave(file.path(raw_dir, paste0("Volcano_", contrast, ".pdf")), p, width = 10, height = 8)
        }
}

if (!length(all_results)) stop("No DESeq2 comparisons completed successfully.")

all_deg_results <- bind_rows(all_results)

readr::write_csv(all_deg_results, file.path(raw_dir, "Microglia_DESeq2_all_comparisons.csv"))

if (length(failed_comparisons)) readr::write_csv(bind_rows(failed_comparisons), file.path(raw_dir, "failed_comparisons.csv"))

### 6) Mad-based crossover filter -------------------------------------------

# Preserve the original filter: full MegaData unsupervised clusters provide
# background distributions of mean raw Xenium counts. Retain a gene only when
# its mean raw count per microglial cell exceeds median + 2*MAD across clusters.
# Note: this selects microglia-enriched genes relative to the full dataset;
# it is not an outlier-removal filter for individual pseudobulk samples.
dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")

mad_table <- calculate_mad_thresholds(dat, datM_clean, z = 2)

readr::write_csv(mad_table, file.path(filtered_dir, "Microglia_MAD_gene_thresholds.csv"))

filtered_genes <- mad_table$gene[mad_table$passes_MAD]

all_deg_results_MAD <- all_deg_results %>% filter(gene %in% filtered_genes)

readr::write_csv(all_deg_results_MAD, file.path(filtered_dir, "Microglia_DESeq2_all_comparisons_MAD_filtered.csv"))

for (contrast in names(all_results)) {
        filtered <- all_results[[contrast]] %>% filter(gene %in% filtered_genes)
        readr::write_csv(filtered, file.path(filtered_dir, paste0("Markergenes_", contrast, "_MAD_filtered.csv")))
}

### 7) Significant deg summaries --------------------------------------------

# Same cutoffs as the ExN DEG overview (adjusted P < 0.05 and |log2FC| >= 1).
padj_cutoff <- 0.05

log2fc_cutoff <- 1

sig_deg <- all_deg_results %>% filter(!is.na(p_val_adj), p_val_adj < padj_cutoff, abs(avg_log2FC) >= log2fc_cutoff)

sig_deg_MAD <- all_deg_results_MAD %>% filter(!is.na(p_val_adj), p_val_adj < padj_cutoff, abs(avg_log2FC) >= log2fc_cutoff)

readr::write_csv(sig_deg, file.path(raw_dir, "Microglia_DESeq2_significant_DEGs.csv"))
readr::write_csv(sig_deg_MAD, file.path(filtered_dir, "Microglia_DESeq2_significant_DEGs_MAD_filtered.csv"))

summary_table <- tibble(comparison = names(all_results)) %>%
        mutate(n_tested = vapply(all_results, nrow, integer(1)),
               n_significant = vapply(names(all_results), function(x) sum(sig_deg$comparison == x), integer(1)),
               n_MAD_retained = vapply(names(all_results), function(x) sum(all_deg_results_MAD$comparison == x), integer(1)),
               n_significant_MAD = vapply(names(all_results), function(x) sum(sig_deg_MAD$comparison == x), integer(1)))

print(summary_table, n = Inf)

readr::write_csv(summary_table, file.path(cfg$output_dir_de, "Microglia_DESeq2_MAD_summary.csv"))
