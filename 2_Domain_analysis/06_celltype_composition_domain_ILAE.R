################################################################################
#
# Xenium regional cell-type composition by ILAE grade (Seurat v5)
#
# Mean within-sample cell-type proportions by spatial domain and ILAE grade; KW BH tests.
# Original analysis: Hope Luo
# Code cleanup and standardization: Larissa Kraus (October 2026)
######

### 0) Setup -------------------------------------------------------------------
suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(broom)
        library(tibble)
})
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)
options(future.globals.maxSize = 100000 * 1024^2)

### 1) Helper functions --------------------------------------------------------
# Reuse an existing object in the global environment; otherwise load its RDS file.
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
dat_rds = NA_character_

dat <- load_if_missing("dat", dat_rds, "Seurat")

# Cell type palette used throughout the manuscript
pal_celltype <- c(
  "ExN"      = "#E69F00",
  "InN"      = "#56B4E9",
  "Astro"    = "#009E73",
  "Micro"    = "#F0E442",
  "Oligo"    = "#0072B2",
  "OPC"      = "#D55E00",
  "Vascular" = "#CC79A7",
  "CSF"      = "#666666")


### 3) Regional composition and Kruskal-Wallis tests --------------------------
### STACKED BAR PLOT FOR GENERAL CELL TYPES ACROSS SPATIAL DOMAINS AND ILAE GRADES

# Defines celltype_prop_df: cell type proportion differences
celltype_prop_df <- dat@meta.data %>%
        # Filter out missing annotations if present
        filter(!is.na(celltype), !is.na(region_md), !is.na(ILAE_score), !is.na(orig.ident)) %>%
        # Group by sample, spatial domain, disease grade, and cell type
        group_by(orig.ident, region_md, ILAE_score, celltype) %>%
        summarise(n_cells = n(), .groups = "drop_last") %>%
        # Calculate relative proportion per sample within each domain
        mutate(proportion = n_cells / sum(n_cells)) %>%
        ungroup()

# Aggregate mean proportions for stacked bar plotting
mean_celltype_df <- celltype_prop_df %>%
        filter(tolower(as.character(region_md)) != "void") %>%
        group_by(region_md, ILAE_score, celltype) %>%
        summarise(mean_proportion = mean(proportion), .groups = "drop") %>%
        mutate(ILAE_score = factor(ILAE_score, levels = c("P", "0", "1", "2", "3")))

# Render Stacked Bar Plot, all stacked bars adjusted to 100%
# x axis is ILAE grades across spatial domains
# y axis is proportion of each cell type
p_celltype_bar <- ggplot(mean_celltype_df, aes(x = ILAE_score, y = mean_proportion, fill = celltype)) +
        geom_bar(stat = "identity", position = "fill", width = 0.72) +
        facet_wrap(~ region_md, nrow = 1) +
        scale_y_continuous(
                labels = percent_format(accuracy = 1),
                expand = expansion(mult = c(0, 0.08))
        ) +
        # Apply custom palette
        scale_fill_manual(
                values = pal_celltype
        ) +
        theme_classic(base_size = 11) +
        theme(
                axis.text.x = element_text(size = 8, face = "bold", color = "black"),
                axis.text.y = element_text(size = 8, color = "black"),
                axis.title = element_text(size = 10, face = "bold", color = "black"),
                strip.background = element_rect(fill = "#ECEFF1", color = NA),
                strip.text = element_text(face = "bold", size = 10, color = "black"),
                panel.spacing = unit(0.8, "lines"),
                legend.title = element_text(face = "bold", size = 10),
                legend.text = element_text(size = 8),
                legend.position = "right",
                plot.title = element_text(face = "bold", size = 12, color = "black", margin = margin(b = 10))
        ) +
        labs(
                x = "ILAE Grade",
                y = "Relative Cell Type Composition",
                fill = "Cell Type",
                title = "General Cell Type Composition Across ILAE Grades by Spatial Domain"
        )

print(p_celltype_bar)



# Run Kruskal-Wallis test per spatial domain and per cell type across ILAE grades
# and applies the Benjamini-Hochberg False Discovery Rate adjustment across all tests
kw_celltype_results <- celltype_prop_df %>%
        filter(tolower(as.character(region_md)) != "void") %>%
        mutate(
                region_md = droplevels(factor(region_md)), # drops “void” spatial category
                ILAE_score = factor(ILAE_score, levels = c("P", "0", "1", "2", "3")) # enforces order of ILAE grades
        ) %>%
        group_by(region_md, celltype) %>%
        filter(n_distinct(ILAE_score) > 1) %>%
        do(tidy(kruskal.test(proportion ~ ILAE_score, data = .))) %>%
        ungroup() %>%
        # Correct p-values using Benjamini-Hochberg (FDR)
        mutate(
                p.adj = p.adjust(p.value, method = "BH"),
                neg_log_padj = -log10(pmax(p.adj, 1e-10)),
                significance = case_when(
                        p.adj < 0.001 ~ "***",
                        p.adj < 0.01  ~ "**",
                        p.adj < 0.05  ~ "*",
                        TRUE          ~ ""
                )
        )

# View the full statistical results table
print(kw_celltype_results)

# Build complete grid matrix
heatmap_data <- kw_celltype_results %>%
        complete(region_md, celltype, fill = list(p.adj = 1, neg_log_padj = 0, significance = ""))

# Define linear gradient limits, with maximum limit being at least p = 0.05
max_val <- max(heatmap_data$neg_log_padj, na.rm = TRUE)
upper_limit <- max(max_val, 1.3)

# Render Heatmap
# Color intensity of tiles represents -log10(FDR-p),
# significance asterisks is overlaid and centered inside significant tiles
p_kw_heatmap <- ggplot(heatmap_data, aes(x = region_md, y = celltype, fill = neg_log_padj)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = significance), color = "white", size = 8, vjust = 0.75, fontface = "bold") +
        scale_fill_gradient(
                low = "#F5F5F5",
                high = "#004D40",
                limits = c(0, upper_limit),
                name = expression(-log[10](p[adj]))
        ) +
        theme_minimal(base_size = 11) +
        theme(
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold", color = "black"),
                axis.text.y = element_text(size = 10, face = "bold", color = "black"),
                axis.title = element_text(size = 10, face = "bold", color = "black"),
                panel.grid = element_blank(),
                legend.title = element_text(size = 10, face = "bold"),
                legend.position = "right",
                plot.title = element_text(face = "bold", size = 12, color = "black", margin = margin(b = 10))
        ) +
        labs(
                x = "Spatial Domain (Region)",
                y = "Cell Type",
                title = "Kruskal-Wallis Significance Heatmap Across Spatial Domains",
                subtitle = "Continuous color intensity represents -log10(FDR-adjusted p-value)",
                caption = "* p.adj < 0.05 | ** p.adj < 0.01 | *** p.adj < 0.001"
        )

print(p_kw_heatmap)
