# 04_exn_deg.R
#
# Pseudobulk differential expression across ILAE groups within ExN subtypes
# Using MegaData_v8
# Author: Larissa Kraus
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# 1. Aggregate raw Xenium counts per patient and broad ExN subtype
# 2. Run pairwise DESeq2 comparisons and save volcano plots
# 3. Summarize significant DEGs and plot direction-of-change bubble chart
# 4. Explore selected-gene expression across ILAE groups

### LOAD PACKAGES ##############################################################
suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(purrr)
        library(readr)
        library(stringr)
        library(ggplot2)
        library(ggtext)
        library(EnhancedVolcano)
        library(DESeq2)
        library(patchwork)
        library(FSA)
        library(conflicted)
})

### PACKAGE CONFLICTS
conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")
conflicted::conflict_prefer("lag", "dplyr")
conflicted::conflict_prefer("first", "dplyr")

options(future.globals.maxSize = 100000 * 1024^2)

### COLOR PALETTES #############################################################
pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00",
                    "2" = "#00C0C5", "3" = "#E69F00")

### HELPER FUNCTIONS ###########################################################

# Reuse an existing object; otherwise load it from a user-configured RDS path.
# Returns the object without overwriting an existing object in the calling environment.
load_if_missing <- function(object_name, path, envir = .GlobalEnv) {
        if (exists(object_name, envir = envir, inherits = FALSE)) {
                message("Reusing ", object_name, " from the environment.")
                return(get(object_name, envir = envir, inherits = FALSE))
        }
        if (is.null(path) || !nzchar(path) || !file.exists(path)) {
                stop("Missing ", object_name, ". Set its input RDS path in cfg: ", path)
        }
        message("Loading ", object_name, " from: ", path)
        readRDS(path)
}


# Shared patient-level pseudobulk DEG helpers | October 2026
# Identity format: ILAE_group_cluster (cluster may contain underscores).
# Positive avg_log2FC always means higher expression in the first group.

# DEG helpers use the fixed ILAE ordering specified below.
ilae_order <- c("P", "0", "1", "2", "3")

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

# For each subtype, plot patient-level pseudobulk expression across ILAE groups.
# The omnibus Kruskal-Wallis p-value is unadjusted across subtypes/genes.
# When significant, post hoc comparisons use BH-adjusted Dunn p-values
# (or a two-group Wilcoxon test when only two ILAE groups are represented).
Test_Plot_DEGs <- function(seurat_obj, gene, cluster_col = "broad_subtype_ExN",
                           ilae_col = "ILAE_score", assay = "Xenium",
                           posthoc_p_cutoff = 0.05, cluster_list = NULL, ncol = 3) {
  
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) stop("Missing cluster column: ", cluster_col)
  if (!ilae_col %in% colnames(seurat_obj@meta.data)) stop("Missing ILAE column: ", ilae_col)
  if (!gene %in% rownames(seurat_obj[[assay]])) {
    warning("Gene not found in assay ", assay, ": ", gene)
    return(NULL)
  }
  if (is.null(cluster_list)) cluster_list <- unique(as.character(seurat_obj@meta.data[[cluster_col]]))
  cluster_list <- cluster_list[!is.na(cluster_list)]
  
  # Fetch once rather than for every subtype; explicitly use Xenium assay.
  expr <- FetchData(seurat_obj[[assay]], vars = gene, layer = "data")
  meta <- seurat_obj@meta.data
  meta$expression <- expr[rownames(meta), gene]
  meta$ILAE_group <- factor(sub("^g(?=[0-9]+$)", "", as.character(meta[[ilae_col]]), perl = TRUE),
                            levels = c("P", "0", "1", "2", "3"))
  
  plot_list <- list()
  test_results <- list()
  posthoc_results <- list()
  
  for (cl in cluster_list) {
    message("Processing ", gene, " in ", cl)
    df <- meta %>% filter(as.character(.data[[cluster_col]]) == cl,
                          !is.na(ILAE_group), is.finite(expression)) %>%
      droplevels()
    n_groups <- n_distinct(df$ILAE_group)
    if (n_groups < 2 || n_distinct(df$expression) < 2) {
      warning("Insufficient groups or expression variation for ", gene, " in ", cl)
      next
    }
    
    main_test <- kruskal.test(expression ~ ILAE_group, data = df)
    posthoc_df <- NULL
    if (is.finite(main_test$p.value) && main_test$p.value < posthoc_p_cutoff) {
      if (n_groups == 2) {
        groups <- levels(droplevels(df$ILAE_group))
        w <- wilcox.test(expression ~ ILAE_group, data = df,
                         exact = FALSE, correct = FALSE)
        posthoc_df <- data.frame(group1 = groups[1], group2 = groups[2], P.adj = w$p.value)
      } else {
        dunn <- FSA::dunnTest(expression ~ ILAE_group, data = df, method = "bh")$res
        split_groups <- strsplit(dunn$Comparison, " - ", fixed = TRUE)
        posthoc_df <- data.frame(group1 = vapply(split_groups, `[`, "", 1),
                                 group2 = vapply(split_groups, `[`, "", 2),
                                 P.adj = dunn$P.adj)
      }
      posthoc_df$significance <- ifelse(posthoc_df$P.adj < 0.001, "***",
                                        ifelse(posthoc_df$P.adj < 0.01, "**",
                                               ifelse(posthoc_df$P.adj < 0.05, "*", "ns")))
    }
    
    subtitle_text <- paste0("Kruskal-Wallis P = ", format(main_test$p.value, digits = 2, scientific = TRUE))
    if (!is.null(posthoc_df)) {
      sig <- posthoc_df[!is.na(posthoc_df$P.adj) & posthoc_df$P.adj < posthoc_p_cutoff, ]
      if (nrow(sig)) subtitle_text <- paste0(subtitle_text, "\nPost hoc: ",
                                             paste0(sig$group1, " - ", sig$group2, " (", sig$significance, ")", collapse = "; "))
    }
    
    p <- ggplot(df, aes(x = ILAE_group, y = expression)) +
      geom_boxplot(aes(fill = ILAE_group), outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
      scale_fill_manual(values = pal_ILAE_score, drop = FALSE) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      labs(title = paste(gene, "in", cl), subtitle = subtitle_text,
           x = "ILAE score", y = "Pseudobulk expression (Xenium)") +
      theme_classic() +
      theme(legend.position = "none", plot.title = element_text(size = 8),
            plot.subtitle = element_text(size = 7))
    
    plot_list[[cl]] <- p
    test_results[[cl]] <- main_test
    posthoc_results[[cl]] <- posthoc_df
  }
  
  combined_plot <- if (length(plot_list)) wrap_plots(plot_list, ncol = ncol) else NULL
  list(plot_list = plot_list, combined_plot = combined_plot,
       test_results = test_results, posthoc_results = posthoc_results)
}


### USER-EDITABLE PATHS ########################################################
# Use the same saved ExN object as Scripts 01-03.
cfg <- list(
        datEx_rds = "path/to/directory/data/datEX_v8.rds",
        output_dir_plots = "path/to/directory/plots",
        output_dir_de = "path/to/directory/de_marker"
)

output_dir <- file.path(cfg$output_dir_de, "ExN_v8", "broad")
overview_dir <- file.path(output_dir, "DEG_overview")
plot_dir <- file.path(cfg$output_dir_plots, "ExN_v8")
dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

### LOAD EXN OBJECT ############################################################
datEx <- load_if_missing("datEx", cfg$datEx_rds) else {
        message("Using datEx already in the environment.")
}

### PSEUDOBULK AGGREGATION #####################################################
# Sum raw Xenium counts per patient, ILAE group and broad ExN subtype.
# Patients (not individual cells) are the pseudobulk replicates.
pseudo_datEx <- AggregateExpression(datEx, assays = "Xenium", return.seurat = TRUE,
                                    group.by = c("ILAE_score", "orig.ident", "broad_subtype_ExN"))
DefaultAssay(pseudo_datEx) <- "Xenium"

# Seurat may prefix numeric group labels with 'g' when constructing identities.
# Preserve actual metadata group values when building identities.
pseudo_datEx$ILAE_clust <- paste(pseudo_datEx$ILAE_score, pseudo_datEx$broad_subtype_ExN, sep = "_")
Idents(pseudo_datEx) <- "ILAE_clust"
print(table(pseudo_datEx$ILAE_score, useNA = "ifany"))
print(head(Cells(pseudo_datEx)))

### PAIRWISE DESEQ2 ############################################################
# DESeq2 is applied to pseudobulk samples, not individual Xenium cells.
run_BulkDE_and_volcano_plots(pseudo.data = pseudo_datEx, output_dir = output_dir,
                             test_method = "DESeq2")

### DEG OVERVIEW ###############################################################
padj_cutoff <- 0.05
log2fc_cutoff <- 1
subtype_order <- c("GC", "MC", "CA2-3", "dCA3", "sCA1", "dCA1", "sSub", "mSub", "dSub", "L6b")
ilae_order <- c("P", "0", "1", "2", "3")
comparison_order <- c("P vs 0", "P vs 1", "P vs 2", "P vs 3", "0 vs 1", "0 vs 2",
                      "0 vs 3", "1 vs 2", "1 vs 3", "2 vs 3")

deg_files <- list.files(output_dir, pattern = "^Markergenes_.*\\.csv$", full.names = TRUE)
message("Found ", length(deg_files), " DEG files.")
if (!length(deg_files)) stop("No Markergenes CSV files found in ", output_dir)

file_metadata <- purrr::map_dfr(deg_files, parse_deg_filename) %>%
        mutate(subtype = cluster)
failed_files <- file_metadata %>% filter(!filename_parsed | !subtype %in% subtype_order)
if (nrow(failed_files)) {
        print(failed_files)
        stop("Unexpected DEG filenames, group labels, or subtype names")
}
all_deg_results <- map_dfr(deg_files, read_deg_file) %>%
        left_join(file_metadata %>% select(source_file, subtype, comparison, group_1, group_2), by = "source_file") %>%
        mutate(significant = !is.na(p_val_adj) & p_val_adj < padj_cutoff & abs(avg_log2FC) >= log2fc_cutoff,
               direction = case_when(significant & avg_log2FC > 0 ~ "Up",
                                     significant & avg_log2FC < 0 ~ "Down", TRUE ~ "Not significant"),
               higher_in = case_when(significant & avg_log2FC > 0 ~ group_1,
                                     significant & avg_log2FC < 0 ~ group_2, TRUE ~ NA_character_),
               subtype = factor(subtype, levels = subtype_order),
               comparison = factor(comparison, levels = comparison_order))
write_csv(all_deg_results, file.path(overview_dir, "all_DE_results_combined.csv"))

deg_summary <- all_deg_results %>%
        group_by(subtype, comparison, group_1, group_2) %>%
        summarise(Up = sum(significant & direction == "Up", na.rm = TRUE),
                  Down = sum(significant & direction == "Down", na.rm = TRUE), .groups = "drop") %>%
        mutate(total_DEGs = Up + Down,
               direction_balance = if_else(total_DEGs > 0, (Up - Down) / total_DEGs, 0))

# Bubble size = significant DEGs; green = more genes higher in first group;
# purple = more genes higher in second group (standardized comparison order).
comparison_labels <- setNames(vapply(strsplit(comparison_order, " vs "), function(pair) {
        paste0("<span style='color:#1B7837;'>", pair[1], "</span><br>vs<br>",
               "<span style='color:#762A83;'>", pair[2], "</span>")
}, character(1)), comparison_order)

p_bubble <- ggplot(deg_summary, aes(x = comparison, y = subtype, size = total_DEGs, color = direction_balance)) +
        geom_point(alpha = 0.9) +
        scale_size_area(max_size = 16, breaks = scales::pretty_breaks(n = 4), name = "Total DEGs") +
        scale_color_gradient2(low = "#762A83", mid = "#D8C99B", high = "#1B7837",
                              midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                              name = "Direction of\nregulation") +
        scale_x_discrete(labels = comparison_labels) +
        scale_y_discrete(limits = rev(subtype_order)) +
        labs(x = NULL, y = NULL, title = "Differential expression across ILAE groups",
             subtitle = "Bubble size indicates total DEGs; color indicates direction of regulation") +
        theme_classic(base_size = 12) +
        theme(axis.text.x = ggtext::element_markdown(hjust = 0.5, vjust = 1),
              axis.text.y = element_text(color = "black"),
              axis.line = element_line(color = "black", linewidth = 0.5),
              axis.ticks = element_line(color = "black", linewidth = 0.5),
              panel.grid = element_blank(), legend.position = "right",
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(color = "grey30"))
print(p_bubble)

### SELECTED GENES: BROAD EXN SUBTYPES #########################################
# Exploratory patient-level pseudobulk expression analysis, distinct from DESeq2.
# SNCG: original broad-subtype comparison retained from Script 05.
pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00",
                    "2" = "#00C0C5", "3" = "#E69F00")

SNCG_ExN <- Test_Plot_DEGs(pseudo_datEx, gene = "SNCG", ncol = 10)

SNCG_ExN$combined_plot

### CLEAN UP ###################################################################
# Retain datEx, pseudo_datEx, DEG summaries and SNCG_ExN for inspection.
rm(file_metadata, failed_files, deg_files, comparison_labels)
