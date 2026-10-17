# 03_inn_region_deg.R
#
# Regional pseudobulk differential expression in inhibitory neurons | MegaData_v8
# Input: datIN_v8.rds from script 01, and per-cell anatomical region annotations
# LK, October 2026
#
# This script:
# 1. Aggregates SST, PV and VIP expression by patient and anatomical region
# 2. Performs pairwise ILAE-group comparisons using DESeq2
# 3. Saves DEG tables and volcano plots
# 4. Summarizes regional DEG counts and recurrent genes

### Setup ----------------------------------------------------------------------

suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(purrr)
        library(stringr)
        library(readr)
        library(EnhancedVolcano)
})

options(future.globals.maxSize = 100000 * 1024^2)

### Helper Functions -----------------------------------------------------------

# Keep the same DEG functions used in the ExN analysis.
generate_combinations <- function(idents, get_cluster) {
        idents <- unique(as.character(idents))
        clusters <- get_cluster(idents)
        idents_by_cluster <- split(idents, clusters)
        
        combinations <- lapply(idents_by_cluster, function(ids) {
                ids <- unique(ids)
                if (length(ids) < 2) return(list())
                
                # Ensure comparisons follow the intended ILAE order.
                scores <- sub("_.*$", "", ids)
                ids <- ids[order(match(scores, c("P", "0", "1", "2", "3")))]
                
                combn(ids, 2, simplify = FALSE)
        })
        
        combinations[lengths(combinations) > 0]
}

run_BulkDE_and_volcano_plots_v2 <- function(pseudo.data, output_dir, test_method = "DESeq2", get_cluster = NULL) {
        if (is.null(get_cluster)) stop("You must supply a function to extract cluster names from identities.")
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        
        idents <- unique(as.character(Seurat::Idents(pseudo.data)))
        # The pair order determines the sign of avg_log2FC; retain the
        # original combination logic for consistency with the ExN scripts.
        score_combinations <- generate_combinations(idents, get_cluster = get_cluster)
        
        if (length(score_combinations) == 0) {
                message("No valid pairwise comparisons found.")
                return(invisible(NULL))
        }
        
        for (cluster in names(score_combinations)) {
                for (scores in score_combinations[[cluster]]) {
                        ident1 <- scores[1]
                        ident2 <- scores[2]
                        message("Running ", ident1, " vs ", ident2)
                        
                        de_markers <- tryCatch(FindMarkers(pseudo.data, ident.1 = ident1, ident.2 = ident2, test.use = test_method), error = function(e) {
                                message("Error in DE analysis for ", ident1, " vs ", ident2, ": ", e$message)
                                NULL
                        })
                        
                        if (is.null(de_markers) || nrow(de_markers) == 0) {
                                message("No DE results for ", ident1, " vs ", ident2, ", skipping.")
                                next
                        }
                        
                        safe_ident1 <- gsub("[^A-Za-z0-9_-]", "_", ident1)
                        safe_ident2 <- gsub("[^A-Za-z0-9_-]", "_", ident2)
                        safe_cluster <- gsub("[^A-Za-z0-9_-]", "_", cluster)
                        
                        csv_file <- file.path(output_dir, paste0("Markergenes_", safe_ident1, "_vs_", safe_ident2, "_Group_", safe_cluster, ".csv"))
                        write.csv(de_markers, csv_file, row.names = TRUE)
                        
                        de_markers$p_val_adj_plot <- pmax(de_markers$p_val_adj, .Machine$double.xmin)
                        finite_fc <- de_markers$avg_log2FC[is.finite(de_markers$avg_log2FC)]
                        finite_p <- -log10(de_markers$p_val_adj_plot[is.finite(de_markers$p_val_adj_plot)])
                        
                        if (length(finite_fc) == 0 || length(finite_p) == 0) {
                                message("No finite values available for volcano plot for ", ident1, " vs ", ident2, ", skipping plot.")
                                next
                        }
                        
                        x_min <- min(finite_fc, na.rm = TRUE)
                        x_max <- max(finite_fc, na.rm = TRUE)
                        y_max <- max(finite_p, na.rm = TRUE)
                        arrow_length <- max(diff(range(finite_fc, na.rm = TRUE)) * 0.08, 0.25)
                        
                        p <- EnhancedVolcano(de_markers, lab = rownames(de_markers), x = "avg_log2FC", y = "p_val_adj_plot", drawConnectors = TRUE, labSize = 4, pCutoff = 0.05) +
                                ggtitle(paste("DE between", ident1, "and", ident2, "in", cluster)) +
                                annotate("text", x = x_max, y = y_max + 4, label = paste("Upregulated in", ident1), hjust = 1, vjust = 0, size = 5) +
                                annotate("text", x = x_min, y = y_max + 4, label = paste("Upregulated in", ident2), hjust = 0, vjust = 0, size = 5) +
                                annotate("segment", x = x_min + arrow_length, xend = x_min, y = y_max + 3.2, yend = y_max + 3.2, arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
                                annotate("segment", x = x_max - arrow_length, xend = x_max, y = y_max + 3.2, yend = y_max + 3.2, arrow = arrow(type = "closed", length = unit(0.1, "inches")))
                        
                        pdf_file <- file.path(output_dir, paste0("Volcano_", safe_ident1, "_vs_", safe_ident2, "_Group_", safe_cluster, ".pdf"))
                        grDevices::pdf(pdf_file, width = 10, height = 8)
                        print(p)
                        grDevices::dev.off()
                }
        }
        
        invisible(score_combinations)
}

# Parse DEG filenames
parse_deg_filename <- function(file) {
        filename <- tools::file_path_sans_ext(basename(file))
        subtype <- basename(dirname(file))
        
        # Remove filename prefix and trailing _Group_...
        identity_part <- str_remove(filename, "^Markergenes_")
        identity_part <- str_remove(identity_part, "_Group_.*$")
        
        # Split comparison
        identities <- str_split_fixed(identity_part, "_vs_", 2)
        id1 <- identities[, 1]
        id2 <- identities[, 2]
        
        # Extract ILAE groups
        group_1_raw <- str_extract(id1, "^[^_]+")
        group_2_raw <- str_extract(id2, "^[^_]+")
        group_1 <- str_remove(group_1_raw, "^g")
        group_2 <- str_remove(group_2_raw, "^g")
        
        # Extract region by removing group + subtype
        region_1 <- str_remove(id1, paste0("^", group_1_raw, "_", subtype, "_"))
        region_2 <- str_remove(id2, paste0("^", group_2_raw, "_", subtype, "_"))
        
        tibble(full_path = file, source_file = basename(file), subtype = subtype, region = region_1,
               group_1 = group_1, group_2 = group_2,
               comparison = paste(group_1, "vs", group_2),
               filename_parsed = !is.na(group_1) & !is.na(group_2) & !is.na(region_1) & region_1 != "" & region_1 == region_2 & subtype %in% subtype_order)
}

# Read DEG tables
# First unnamed column contains gene names
read_deg_file <- function(file) {
        x <- read_csv(file, show_col_types = FALSE, name_repair = "unique")
        names(x)[1] <- "gene"
        required_columns <- c("gene", "p_val", "avg_log2FC", "p_val_adj")
        if (!all(required_columns %in% names(x))) stop("Missing required columns in ", 
                                                       basename(file), ": ",
                                                       paste(setdiff(required_columns, names(x)),
                                                             collapse = ", "))
        x %>% mutate(gene = as.character(gene), source_file = basename(file), full_path = file) %>% 
                select(gene, p_val, avg_log2FC, p_val_adj, source_file, full_path)
}

### User-editable configuration -------------------------------------------------

# Supply the annotated InN object from script 01 and the full-dataset regional
# metadata generated in script XXXX (to be added after domain analysis is done). 
# Row names of the regional metadata must be Xenium cell identifiers.

cfg <- list(
        datIN_rds = "path/data/datIN_v8.rds",
        region_metadata_rds = "path/data/megadata_metadata_md_v8.rds",
        output_dir_deg = "path/de_marker/InN_v8/TableS3_broad_region_md"
)

for (path in c(cfg$datIN_rds, cfg$region_metadata_rds)) {
        if (!file.exists(path)) stop("Required input not found: ", path)
}

dir.create(cfg$output_dir_deg, recursive = TRUE, showWarnings = FALSE)

### Load annotated data ------------------------------------------

datIN <- readRDS(cfg$datIN_rds)

required_meta <- c("seurat_clusters", "fine_subtype_InhN", "broad_subtype_InhN",
                   "ILAE_score", "orig.ident", "region_md")

missing_meta <- setdiff(required_meta, colnames(datIN@meta.data))

if (length(missing_meta)) stop("Missing metadata in datIN: ", paste(missing_meta, collapse = ", "))

### Region-specific pseudobulk differential expression ------------------------

# Aggregate Xenium counts by patient, ILAE group, inhibitory subtype and region.
# DESeq2 comparisons are made between patient-level pseudobulk profiles.
# Run for main IN subtype and their respective regions

subtypes_keep <- c("SST", "PV", "VIP")
regions_keep <- c("SL.SR.SLM", "CA1", "GCL", "SUB", "CA2_CA4", "GABA")

for (subtype_i in subtypes_keep) {
        message("Processing ", subtype_i)
        
        md <- datIN@meta.data
        cells_keep <- rownames(md)[!is.na(md$broad_subtype_InhN) & !is.na(md$region_md) & as.character(md$broad_subtype_InhN) == subtype_i & as.character(md$region_md) %in% regions_keep]
        
        if (length(cells_keep) == 0) {
                message("No cells for subtype ", subtype_i, " in selected regions, skipping.")
                next
        }
        
        datIN_subset <- subset(datIN, cells = cells_keep)
        message("Found ", ncol(datIN_subset), " cells for ", subtype_i)
        
        # Each aggregated profile represents one patient in a subtype-region group.
        # Check the printed patient counts before interpreting DESeq2 results;
        # the function below catches DESeq2 errors but does not enforce a
        # minimum patient count per group.
        pseudo_subset <- AggregateExpression(datIN_subset, assays = "Xenium", return.seurat = TRUE, group.by = c("ILAE_score", "orig.ident", "broad_subtype_InhN", "region_md"))
        
        # Compare ILAE groups only within the same subtype and region.
        pseudo_subset$ILAE_type_region <- paste(as.character(pseudo_subset$ILAE_score), as.character(pseudo_subset$broad_subtype_InhN), as.character(pseudo_subset$region_md), sep = "_")
        Idents(pseudo_subset) <- "ILAE_type_region"
        
        print(table(Idents(pseudo_subset)))
        
        run_BulkDE_and_volcano_plots_v2(pseudo.data = pseudo_subset, output_dir = file.path(cfg$output_dir_deg, subtype_i), test_method = "DESeq2", get_cluster = function(ids) sub("^[^_]+_", "", ids))
}


### Combine and summarize DEG results -----------------------------------------

# Combine saved DESeq2 tables, annotate each comparison, and summarize the
# number and direction of significant genes by subtype and region.
output_dir <- cfg$output_dir_deg

padj_cutoff <- 0.05
log2fc_cutoff <- 1

# Desired plotting order
subtype_order <- c("SST", "PV", "VIP")
region_order <- c("GCL","CA2_CA4","CA1", "SUB", "SL.SR.SLM", "GABA")
comparison_order <- c("P vs 0", "P vs 1", "P vs 2", "P vs 3", "0 vs 1", "0 vs 2",
                      "0 vs 3", "1 vs 2", "1 vs 3", "2 vs 3")

# Locate subtype-specific DEG tables
deg_files <- list.files(output_dir, pattern = "^Markergenes_.*\\.csv$", 
                        full.names = TRUE, recursive = TRUE)

message("Found ", length(deg_files), " DEG files.")

if (length(deg_files) == 0) stop("No Markergenes CSV files were found in: ", output_dir)

### Load and validate DEG files ------------------------------------------------

file_metadata <- map_dfr(deg_files, parse_deg_filename) %>%
        mutate(region = recode(region,
                               "CA2-CA4" = "CA2_CA4",
                               "SL_SR_SLM" = "SL.SR.SLM"))

# Check parsing
failed_files <- file_metadata %>% filter(!filename_parsed)

if (nrow(failed_files) > 0) {
        print(failed_files)
        stop("Some filenames could not be parsed. Check the filenames printed above.")
}

# Combine DEG results
# Positive avg_log2FC = higher in group_1
# Negative avg_log2FC = higher in group_2
all_deg_results <- map_dfr(deg_files, read_deg_file) %>%
        left_join(file_metadata, by = c("full_path", "source_file")) %>%
        mutate(significant = !is.na(p_val_adj) & p_val_adj < padj_cutoff & abs(avg_log2FC) >= log2fc_cutoff,
               direction = case_when(significant & avg_log2FC > 0 ~ "Up", significant & avg_log2FC < 0 ~ "Down", TRUE ~ "Not significant"),
               subtype = factor(subtype, levels = subtype_order),
               region = factor(region, levels = c(region_order, setdiff(unique(region), region_order))),
               comparison = factor(comparison, levels = c(comparison_order, setdiff(unique(comparison), comparison_order))))

# Count up- and downregulated genes
# Positive log2FC indicates higher expression in the first ILAE group.
deg_summary <- all_deg_results %>%
        filter(significant) %>%
        count(subtype, region, comparison, group_1, group_2, direction, name = "n_genes") %>%
        pivot_wider(names_from = direction, values_from = n_genes, values_fill = 0)

if (!"Up" %in% names(deg_summary)) deg_summary$Up <- integer(nrow(deg_summary))
if (!"Down" %in% names(deg_summary)) deg_summary$Down <- integer(nrow(deg_summary))

deg_template <- all_deg_results %>%
        distinct(subtype, region, comparison, group_1, group_2)

deg_summary <- deg_template %>%
        left_join(deg_summary, by = c("subtype", "region", "comparison", "group_1", "group_2")) %>%
        mutate(Up = replace_na(Up, 0),
               Down = replace_na(Down, 0),
               total_DEGs = Up + Down,
               direction_balance = if_else(total_DEGs > 0, (Up - Down) / total_DEGs, 0),
               comparison = factor(comparison, levels = c(comparison_order, setdiff(unique(comparison), comparison_order))))


### Plot DEG counts and direction ----------------------------------------------

p_bubble <- ggplot(deg_summary, aes(x = comparison, y = region, 
                                    size = total_DEGs, color = direction_balance)) +
        geom_point(alpha = 0.9) +
        facet_wrap(~subtype, ncol=1)+
        scale_size_area(max_size = 16, breaks = scales::pretty_breaks(n = 4), name = "Total DEGs") +
        scale_color_gradient2(low = "#762A83", mid = "#D8C99B", high = "#1B7837",
                              midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                              name = "Direction of\nregulation") +
        labs(x = NULL, y = NULL, size = "Total DEGs",
             title = "Differential expression across inhibitory-neuron subtypes and regions", 
             subtitle = paste0("Adjusted p < ", padj_cutoff, "; |avg log2FC| ≥ ", log2fc_cutoff)) +
        theme_classic(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(color = "black"),
              axis.line = element_line(color = "black", linewidth = 0.5),
              axis.ticks = element_line(color = "black", linewidth = 0.5),
              panel.grid = element_blank(),
              legend.position = "right",
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(color = "grey30"))
p_bubble
