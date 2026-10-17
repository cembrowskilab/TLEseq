################################################################################
#
# 07_microglia_exploratory_pseudobulk_boxplots.R
#
# Xenium analysis of microglia in human hippocampus samples
# Using MegaData_v8
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# Exploratory microglia pseudobulk expression plots


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
        library(FSA)
        library(ggpubr)
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


plot_gene <- function(gene) {

        df <- FetchData(pseudo_datM, vars = c(gene, "ILAE_score")) %>%
                setNames(c("expression", "ILAE_score")) %>%
                mutate(ILAE_score = factor(ILAE_score, levels = ilae_levels)) %>%
                drop_na() %>%
                droplevels()

        posthoc <- NULL
        if (isTRUE(kruskal.test(expression ~ ILAE_score, data = df)$p.value < 0.05)) {
                posthoc <- FSA::dunnTest(expression ~ ILAE_score, data = df, method = "bh")$res %>%
                        filter(P.adj < 0.05) %>%
                        separate(Comparison, c("group1", "group2"), sep = " - ") %>%
                        mutate(significance = case_when(P.adj < 0.001 ~ "***", P.adj < 0.01 ~ "**", TRUE ~ "*"))
        }

        p <- ggplot(df, aes(ILAE_score, expression, fill = ILAE_score)) +
                geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
                geom_jitter(width = 0.15, alpha = 0.6, color = "black", size = 1.5) +
                scale_fill_manual(values = pal_ILAE_score, drop = FALSE, guide = "none") +
                scale_x_discrete(expand = expansion(add = 0.6)) +
                labs(title = gene, x = NULL, y = "Pseudobulk Expression") +
                theme_classic() +
                theme(plot.title = element_text(face = "italic", size = 9),
                      axis.text.y = element_text(size = 20),
                      plot.margin = margin(2, 2, 2, 2))

        if (isTRUE(nrow(posthoc) > 0)) {
                p <- p + ggpubr::stat_pvalue_manual(
                        mutate(posthoc, y.position = max(df$expression) * seq(1.1, 1.3, length.out = n())),
                        label = "significance", xmin = "group1", xmax = "group2", label.size = 12)
        }
        p
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

pseudo_datM <- AggregateExpression(datM_clean, assays = "Xenium", group.by = c("orig.ident", "ILAE_score"), return.seurat = TRUE)

Idents(pseudo_datM) <- "ILAE_score"

### 4) Plot genes of interest and test kw/dunn pseudobulk ----------------------
# Plot one gene across ILAE scores, with brackets on significant pairs (Kruskal-Wallis, then Dunn's test)

# Plot the genes of interest side by side and save
genes <- c("CCL4", "CD83", "P2RY12", "CX3CR1", "GPNMB", "APOE", "TREM2")

ggsave(file.path(cfg$output_dir_plots, "Microglia_pseudobulk_boxplots.pdf"),
       wrap_plots(lapply(genes, plot_gene), nrow = 1),
       width = 1.8 * length(genes), height = 7)
