################################################################################
#
# 03_megadata_deg_v8.R
#
# Broad-cell-type differential expression analysis using the annotated
# MegaData v8 Seurat object.
#
# This script:
#   1) Loads the annotated MegaData v8 object available from GEO
#   2) Aggregates raw Xenium counts by patient and broad cell type
#   3) Runs pseudobulk DESeq2 across ILAE groups within each cell type
#   4) Saves DE tables and volcano plots
#   5) Summarizes the number and direction of significant DE genes
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
        library(ggplot2)
        library(stringr)
        library(EnhancedVolcano)
        library(ggtext)
})

### PACKAGE CONFLICTS
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)

options(future.globals.maxSize = 25000 * 1024^2)

### 1) User-editable config -----------------------------------------------------

cfg <- list(
        megadata_rds = NA_character_,
        out_de_dir = NA_character_,
        out_plot_dir = NA_character_,
        celltype_order = c("ExN", "InN", "Astro", "Micro", "Oligo", "OPC", "Vascular", "CSF"),
        ilae_order = c("P", "0", "1", "2", "3")
)

dir.create(cfg$out_de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$out_plot_dir, recursive = TRUE, showWarnings = FALSE)

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


# Shared patient-level pseudobulk DEG helpers | October 2026
# Identity format: ILAE_group_cluster (cluster may contain underscores).
# Positive avg_log2FC always means higher expression in the first group.

ilae_order <- c("P", "0", "1", "2", "3")
comparison_order <- unlist(lapply(seq_len(length(ilae_order) - 1), function(i) paste(ilae_order[i], "vs", ilae_order[(i + 1):length(ilae_order)])))

parse_deg_filename <- function(file) {
  stem <- tools::file_path_sans_ext(basename(file))
  match <- stringr::str_match(stem, "^Markergenes_(.+?)_vs_(.+)$")
  id1 <- match[, 2]; id2 <- match[, 3]
  group1 <- sub("_.*$", "", id1); group2 <- sub("_.*$", "", id2)
  group1 <- sub("^g(?=[0-9]+$)", "", group1, perl = TRUE)
  group2 <- sub("^g(?=[0-9]+$)", "", group2, perl = TRUE)
  cluster1 <- sub("^[^_]+_", "", id1); cluster2 <- sub("^[^_]+_", "", id2)
  valid <- !is.na(id1) & !is.na(id2) & grepl("_", id1) & grepl("_", id2) &
    group1 %in% ilae_order & group2 %in% ilae_order & cluster1 == cluster2 &
    match(group1, ilae_order) < match(group2, ilae_order)
  tibble::tibble(full_path = file, source_file = basename(file), group_1 = group1,
                 group_2 = group2, cluster = cluster1, comparison = paste(group1, "vs", group2),
                 filename_parsed = valid)
}

read_deg_file <- function(file) {
  x <- readr::read_csv(file, show_col_types = FALSE, name_repair = "unique")
  names(x)[1] <- "gene"
  required <- c("gene", "p_val", "avg_log2FC", "p_val_adj")
  if (!all(required %in% names(x))) stop("Missing columns in ", basename(file))
  dplyr::mutate(x, gene = as.character(gene), source_file = basename(file), full_path = file)
}

generate_combinations <- function(idents, clusters = NULL) {
  idents <- unique(as.character(idents))
  groups <- sub("_.*$", "", idents)
  groups <- sub("^g(?=[0-9]+$)", "", groups, perl = TRUE)
  if (any(!grepl("_", idents)) || any(!groups %in% ilae_order)) stop("Unexpected pseudobulk identity format or ILAE group")
  cluster_names <- sub("^[^_]+_", "", idents)
  by_cluster <- split(idents, cluster_names)
  if (!is.null(clusters)) by_cluster <- by_cluster[names(by_cluster) %in% clusters]
  pairs <- lapply(by_cluster, function(ids) {
    scores <- sub("_.*$", "", ids)
    scores <- sub("^g(?=[0-9]+$)", "", scores, perl = TRUE)
    if (anyDuplicated(scores)) stop("Duplicate ILAE group identities within a cluster")
    ids <- ids[order(match(scores, ilae_order))]
    if (length(ids) < 2) return(NULL)
    combn(ids, 2, simplify = FALSE)
  })
  pairs[lengths(pairs) > 0]
}

run_BulkDE_and_volcano_plots <- function(pseudo.data, output_dir, test_method = "DESeq2", clusters = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  score_combinations <- generate_combinations(Seurat::Idents(pseudo.data), clusters)
  if (!length(score_combinations)) stop("No valid pairwise comparisons found")
  for (cluster in names(score_combinations)) {
    for (scores in score_combinations[[cluster]]) {
      ident1 <- scores[1]; ident2 <- scores[2]
      message("Running ", ident1, " - ", ident2)
      de_markers <- tryCatch(Seurat::FindMarkers(pseudo.data, ident.1 = ident1, ident.2 = ident2,
                                                 test.use = test_method),
                             error = function(e) {message("DE failed: ", conditionMessage(e)); NULL})
      if (is.null(de_markers) || !nrow(de_markers)) next
      file_id <- paste0(ident1, "_vs_", ident2)
      write.csv(de_markers, file.path(output_dir, paste0("Markergenes_", file_id, ".csv")), row.names = TRUE)
      plot_data <- de_markers
      plot_data$p_val_adj_plot <- pmax(plot_data$p_val_adj, .Machine$double.xmin)
      valid <- is.finite(plot_data$avg_log2FC) & is.finite(plot_data$p_val_adj_plot)
      if (!any(valid)) next
      x_min <- min(plot_data$avg_log2FC[valid]); x_max <- max(plot_data$avg_log2FC[valid])
      y_max <- max(-log10(plot_data$p_val_adj_plot[valid]))
      arrow_length <- max((x_max - x_min) * 0.08, 0.25)
      p <- EnhancedVolcano::EnhancedVolcano(plot_data, lab = rownames(plot_data),
                                            x = "avg_log2FC", y = "p_val_adj_plot", drawConnectors = TRUE, labSize = 4, pCutoff = 0.05) +
        ggplot2::ggtitle(paste("DE between", ident1, "and", ident2, "in", cluster)) +
        ggplot2::annotate("text", x = x_max, y = y_max + 4, label = paste("Upregulated in", ident1), hjust = 1, size = 5) +
        ggplot2::annotate("text", x = x_min, y = y_max + 4, label = paste("Upregulated in", ident2), hjust = 0, size = 5) +
        ggplot2::annotate("segment", x = x_min + arrow_length, xend = x_min, y = y_max + 3.2, yend = y_max + 3.2,
                          arrow = grid::arrow(type = "closed", length = grid::unit(0.1, "inches"))) +
        ggplot2::annotate("segment", x = x_max - arrow_length, xend = x_max, y = y_max + 3.2, yend = y_max + 3.2,
                          arrow = grid::arrow(type = "closed", length = grid::unit(0.1, "inches"))) +
        ggplot2::coord_cartesian(clip = "off")
      ggplot2::ggsave(file.path(output_dir, paste0("Volcano_", file_id, ".pdf")), p, width = 10, height = 8)
    }
  }
  invisible(score_combinations)
}

### 3) Load and prepare MegaData ------------------------------------------------

dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")

dat$ILAE_score <- factor(dat$ILAE_score, levels = cfg$ilae_order)
dat$celltype <- factor(dat$celltype, levels = cfg$celltype_order)

### 4) Pseudobulk DE across ILAE groups ----------------------------------------

# Aggregate raw Xenium counts within patient x cell type
pseudo_dat <- AggregateExpression(dat, assays = "Xenium", return.seurat = TRUE,
                                  group.by = c("ILAE_score", "orig.ident", "celltype"))

# Compare ILAE groups within each broad cell type
pseudo_dat$ILAE_celltype <- paste(pseudo_dat$ILAE_score, pseudo_dat$celltype, sep = "_")
Idents(pseudo_dat) <- "ILAE_celltype"

#Function will run on Idents() defined before
run_BulkDE_and_volcano_plots(pseudo.data = pseudo_dat, output_dir = cfg$out_de_dir,
                             test_method = "DESeq2")

### 5) Summarize number and direction of DE genes -------------------------------

# Read DE tables generated above
csv_files <- list.files(cfg$out_de_dir, pattern = "^Markergenes_.*\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) stop("No DE result files found in cfg$out_de_dir.")

file_metadata <- purrr::map_dfr(csv_files, parse_deg_filename)
if (any(!file_metadata$filename_parsed)) stop("Invalid or old DEG filenames: ", paste(file_metadata$source_file[!file_metadata$filename_parsed], collapse = ", "))
all_deg_results <- purrr::map_dfr(csv_files, read_deg_file) %>%
        left_join(file_metadata, by = c("full_path", "source_file")) %>%
        mutate(significant = !is.na(p_val_adj) & p_val_adj < 0.05 & abs(avg_log2FC) > 1)
deg_summary <- all_deg_results %>%
        group_by(cluster, comparison) %>%
        summarise(Up = sum(significant & avg_log2FC > 0, na.rm = TRUE),
                  Down = sum(significant & avg_log2FC < 0, na.rm = TRUE), .groups = "drop") %>%
        mutate(Total = Up + Down, Direction_balance = if_else(Total > 0, (Up - Down) / Total, 0),
               Comparison = factor(comparison, levels = comparison_order),
               Cluster = factor(cluster, levels = rev(cfg$celltype_order)))
comparison_labels <- setNames(vapply(strsplit(comparison_order, " vs ", fixed = TRUE), function(pair) {
        paste0("<span style='color:#1B7837;'>", pair[1], "</span><br>-<br>",
               "<span style='color:#762A83;'>", pair[2], "</span>")
}, character(1)), comparison_order)

# Plot total number and direction of significant DE genes
p_deg_bubble <- ggplot(deg_summary,
                       aes(x = Comparison, y = Cluster, size = Total, color = Direction_balance)) +
        geom_point(alpha = 0.9) +
        scale_size_area(max_size = 16, breaks = scales::pretty_breaks(n = 4), name = "Total DEGs") +
        scale_color_gradient2(low = "#762A83", mid = "#D8C99B", high = "#1B7837",
                              midpoint = 0, limits = c(-1, 1),
                              breaks = c(-1, -0.5, 0, 0.5, 1),
                              name = "Direction of\nregulation") +
        scale_x_discrete(labels = comparison_labels) +
        labs(x = NULL, y = NULL,
             title = "Differential expression across ILAE groups",
             subtitle = "Bubble size indicates total DEGs; color indicates direction of regulation") +
        theme_classic(base_size = 12) +
        theme(axis.text.x = ggtext::element_markdown(hjust = 0.5, vjust = 1),
              axis.text.y = element_text(color = "black"),
              axis.line = element_line(color = "black", linewidth = 0.5),
              axis.ticks = element_line(color = "black", linewidth = 0.5),
              panel.grid = element_blank(),
              legend.position = "right",
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(color = "grey30"))

p_deg_bubble
