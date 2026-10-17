################################################################################
#
# 05_microglia_spatial_expression.R
#
# Xenium analysis of microglia in human hippocampus samples
# Using MegaData_v8
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# Microglia spatial expression and regional analysis


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
        library(ggrastr)
        library(cowplot)
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


fov_figure <- function(groups, n_col, show_label, body_h) {

        fill <- scale_fill_manual(values = pal_cluster_named[keep_states],
                                  limits = keep_states, na.value = NA, name = NULL)

        # one FOV; cell layer rasterized for file size, everything else stays vector
        panel <- function(fov) {
                p <- ImageDimPlot(FOV_plots, fov = fov, group.by = "cluster_names",
                                  size = 0.5, border.size = NA, axes = FALSE, dark.background = FALSE) +
                        fill +
                        theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))
                p$layers <- lapply(p$layers, rasterise, dpi = 300)
                if (show_label) p <- p + annotate("text", x = -Inf, y = Inf, label = fov,
                                                  hjust = -0.15, vjust = 1.5, size = 3.2, fontface = "bold")
                p
        }

        legend <- cowplot::get_legend(
                ImageDimPlot(FOV_plots, fov = "M70", group.by = "cluster_names", size = 1.2,
                             border.size = NA, axes = FALSE, dark.background = FALSE) +
                        fill + theme(legend.position = "right"))

        blocks <- imap(groups, function(fovs, score) {
                title <- ggplot(data.frame(x = 1, y = 1,
                                           label = ifelse(score == "P", "Postmortem", paste("ILAE", score))),
                                aes(x, y)) +
                        geom_tile(fill = pal_ILAE_score[score], width = 1, height = 0.5) +
                        geom_text(aes(label = label), fontface = "bold", size = 4.5) +
                        theme_void() +
                        theme(plot.margin = margin(2, 2, 2, 2))
                title / wrap_plots(lapply(fovs, panel), ncol = n_col[score]) +
                        plot_layout(heights = c(1, body_h))
        })

        (wrap_plots(blocks, nrow = 1) | wrap_elements(legend)) +
                plot_layout(widths = c(unname(n_col[names(blocks)]), 0.3))
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


### 4) Spatial cell-state maps ----------------------------------------------
# Restrict to the clusters of interest for this figure 
keep_states <- c("Homeo", "DIM", "BAM", "Inflam.DAM", "DAM")
FOV_plots <- subset(datM_clean, subset = cluster_names %in% keep_states)

# Format titles for each ILAE score group
group_titles <- c("P" = "Postmortem", "0" = "ILAE 0", "1" = "ILAE 1", "2" = "ILAE 2", "3" = "ILAE 3")

# Group which samples (FOVs) belong to each ILAE score group
sample_groups <- split(as.character(FOV_plots$Sample), as.character(FOV_plots$ILAE_score)) %>%
        lapply(function(x) sort(unique(x)))

# Sanity checks: 
anyDuplicated(unlist(sample_groups)) == 0  # each sample is in exactly one ILAE group (TRUE)
all(unlist(sample_groups) %in% Images(FOV_plots)) # each sample has a matching FOV image (TRUE)


# Build a full FOV figure for each ILAE group with a title strip and a shared legend.
representative_fovs <- list("P" = "M70", "0" = "I1", "1" = "L15", "2" = "L1", "3" = "D9")
fov_figure(representative_fovs, n_col = c("P" = 1, "0" = 1, "1" = 1, "2" = 1, "3" = 1),
           show_label = FALSE, body_h = 10)

fov_figure(sample_groups, n_col = c("P" = 2, "0" = 2, "1" = 2, "2" = 2, "3" = 1),
           show_label = TRUE, body_h = 20)

### 5) Spatial expression maps ----------------------------------------------
## ImageFeaturePlots – Marker Genes

# Choose genes to plot spatially, one per PDF page
genes_of_interest <- c("CCL4", "CD83")

# Ensure SCT assay so expression values reflect normalized data
DefaultAssay(datM_clean) <- "SCT"

# Plot
pdf(file.path(cfg$output_dir_plots, "Microglia_ImageFeaturePlots_allgenes.pdf"), width = 16, height = 4.5)

for (gene in genes_of_interest) {

        # Per-gene colour scale: min fixed at 0, max rounded up to the nearest whole number
        gene_min <- 0
        gene_max <- ceiling(max(datM_clean[["SCT"]]@data[gene, ], na.rm = TRUE))

        # One block per ILAE score: title strip over a single-FOV expression panel,
        #    all on this gene's scale
        blocks <- imap(representative_fovs, function(fov, score) {
                title <- ggplot(data.frame(x = 1, y = 1,
                                           label = ifelse(score == "P", "Postmortem", paste("ILAE", score))),
                                aes(x, y)) +
                        geom_tile(fill = pal_ILAE_score[score], width = 1, height = 0.5) +
                        geom_text(aes(label = label), fontface = "bold", size = 4.5) +
                        theme_void() +
                        theme(plot.margin = margin(2, 2, 2, 2))

                panel <- ImageFeaturePlot(
                        datM_clean, features = gene, fov = fov,
                        size = 0.9, border.size = NA, axes = FALSE, dark.background = FALSE) +
                        scale_fill_viridis_c(option = "viridis", limits = c(gene_min, gene_max)) +
                        ggtitle(NULL) +
                        theme(legend.position = "none",
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              plot.margin = margin(2, 2, 2, 2))

                title / panel + plot_layout(heights = c(1, 10))
        })

        # Legend for this gene's scale, labeled with its min and rounded max
        legend <- cowplot::get_legend(
                ImageFeaturePlot(datM_clean, features = gene, fov = representative_fovs[["P"]],
                                 size = 1.2, border.size = NA, axes = FALSE, dark.background = FALSE) +
                        scale_fill_viridis_c(option = "viridis",
                                             limits = c(gene_min, gene_max),
                                             breaks = c(gene_min, gene_max),
                                             name = parse(text = paste0("italic('", gene, "')"))) +
                        theme(legend.position = "right"))

        # Combine blocks + legend and add this gene's page to the PDF
        print((wrap_plots(blocks, nrow = 1) | wrap_elements(legend)) +
                      plot_layout(widths = c(rep(1, length(blocks)), 0.3)))
}

dev.off()


### 6) Regional expression --------------------------------------------------
# Pull normalized expression values for the genes of interest chosen above and attach to metadata
expr_df <- FetchData(datM_clean, vars = genes_of_interest, layer = "data")
meta <- datM_clean@meta.data
meta[genes_of_interest] <- expr_df[genes_of_interest]

# Average expression per brain region within each ILAE score
tile_df_ilae <- meta %>%
        filter(!is.na(region_md)) %>%
        group_by(ILAE_score, region_md) %>%
        summarise(across(all_of(genes_of_interest), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
                 .groups = "drop")

# Reshape to long format for faceted tile plotting
tile_long_ilae <- tile_df_ilae %>%
        pivot_longer(cols = starts_with("mean_"), names_to = "gene", values_to = "mean_expr") %>%
        mutate(gene = sub("mean_", "", gene))

tile_long_ilae$region_md <- factor(tile_long_ilae$region_md, levels = levels(datM_clean$region_md))
tile_long_ilae$ILAE_score <- factor(tile_long_ilae$ILAE_score, levels = ilae_levels)
tile_long_ilae$gene <- factor(tile_long_ilae$gene, levels = genes_of_interest)

# Shared colour scale limits, computed across all genes so panels are comparable
max_val <- max(tile_long_ilae$mean_expr, na.rm = TRUE)
min_val <- min(tile_long_ilae$mean_expr, na.rm = TRUE)

# Build and save the tile plot (one row per gene, region x ILAE score grid)
tileplot <- ggplot(tile_long_ilae, aes(x = region_md, y = ILAE_score, fill = mean_expr)) +
        geom_tile(color = "grey90") +
        scale_fill_gradient(low = "white", high = "blue",
                            limits = c(min_val, max_val),
                            breaks = c(min_val, max_val),
                            labels = scales::number(c(min_val, max_val), accuracy = 0.1)) +
        facet_wrap(~ gene, ncol = 1, strip.position = "right") +
        labs(x = "Region", y = "ILAE Score", fill = "Average expression") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              strip.text = element_text(face = "italic", size = 11))

tileplot