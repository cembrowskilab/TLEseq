################################################################################
#
# Xenium Xenium-run batch correlations (Seurat v5)
#
# Across-sample cell-type expression correlations, split by run and ILAE composition.
# Original analysis: Larissa Kraus
# Code cleanup and standardization: Larissa Kraus (October 2026)
######

### 0) Setup -------------------------------------------------------------------
suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(Matrix)
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
cfg <- list(dat_rds = NA_character_, 
            meta_csv = NA_character_)

dat <- load_if_missing("dat", cfg$dat_rds, "Seurat")

if (is.na(cfg$meta_csv) || !file.exists(cfg$meta_csv)) stop("Set cfg$meta_csv")

meta <- read.csv(cfg$meta_csv)

dat$Xenium_run <- meta$Xenium_run[match(as.character(dat$orig.ident), as.character(meta$Sample))]

if (anyNA(dat$Xenium_run)) warning("Missing Xenium_run assignments; check sample IDs")


### 3) Sample and batch correlations ------------------------------------------
assay_use <- "SCT"
layer_use <- "data"
sample_col <- "orig.ident"
celltype_col <- "celltype"
batch_col <- "Xenium_run"

# Extract expression and metadata
expr_batch <- GetAssayData(dat, assay = assay_use, layer = layer_use)
meta_batch <- dat[[]]

# Keep metadata in the same order as expression matrix columns
meta_batch <- meta_batch[colnames(expr_batch), , drop = FALSE]

### SAMPLE X CELL-TYPE COUNTS
sample_celltype_counts <- meta_batch %>%
        count(.data[[sample_col]], .data[[celltype_col]], name = "n_cells")

sample_celltype_counts

# Calculate mean expression across selected cells
get_mean_expression_batch <- function(cells) {
        Matrix::rowMeans(expr_batch[, cells, drop = FALSE])
}

### SAMPLE X CELL-TYPE MEAN EXPRESSION
sample_expression <- lapply(seq_len(nrow(sample_celltype_counts)), function(i) {
  
  sample_i <- as.character(sample_celltype_counts[[sample_col]][i])
  celltype_i <- as.character(sample_celltype_counts[[celltype_col]][i])
  
  # Cells belonging to this sample x cell-type combination
  cells <- rownames(meta_batch)[
    as.character(meta_batch[[sample_col]]) == sample_i &
      as.character(meta_batch[[celltype_col]]) == celltype_i
  ]
  
  data.frame(
    sample = sample_i,
    celltype = celltype_i,
    gene = rownames(expr_batch),
    expression = get_mean_expression_batch(cells)
  )
}) %>%
  bind_rows()

### SAMPLE BATCH METADATA
sample_batch <- meta_batch %>%
        transmute(sample = as.character(.data[[sample_col]]),
                  batch = as.character(.data[[batch_col]])) %>%
        distinct()

sample_batch

### SAMPLE PAIRS
sample_pairs <- combn(unique(sample_batch$sample), 2, simplify = FALSE) %>%
        lapply(function(x) {
                data.frame(sample_A = x[1],sample_B = x[2])
        }) %>%
        bind_rows()

dim(sample_pairs)
head(sample_pairs)

#Add Xenium run for A and B
sample_pairs <- sample_pairs %>%
        left_join(sample_batch %>%
                        rename(sample_A = sample,batch_A = batch),by = "sample_A") %>%
        left_join(sample_batch %>%rename(sample_B = sample,batch_B = batch),by = "sample_B") %>%
        mutate(batch_comparison = if_else(batch_A == batch_B,"Within batch","Across batch"))



### ACROSS-SAMPLE CORRELATIONS
across_sample_cor <- lapply(seq_len(nrow(sample_pairs)), function(i) {
        
        sample_A_i <- sample_pairs$sample_A[i]
        sample_B_i <- sample_pairs$sample_B[i]
        
        # Gene x cell-type matrix for sample A
        mat_A <- sample_expression %>%
                filter(sample == sample_A_i) %>%
                select(gene, celltype, expression) %>%
                pivot_wider(names_from = celltype, values_from = expression) %>%
                tibble::column_to_rownames("gene") %>%
                as.matrix()
        
        # Gene x cell-type matrix for sample B
        mat_B <- sample_expression %>%
                filter(sample == sample_B_i) %>%
                select(gene, celltype, expression) %>%
                pivot_wider(names_from = celltype, values_from = expression) %>%
                tibble::column_to_rownames("gene") %>%
                as.matrix()
        
        # Correlate cell types between samples
        cor_mat <- cor(mat_A, mat_B, method = "pearson")
        
        # Convert matrix to long format
        as.data.frame(as.table(cor_mat)) %>%
                rename(celltype_A = Var1,
                       celltype_B = Var2,
                       correlation = Freq) %>%
                mutate(sample_A = sample_A_i,
                       sample_B = sample_B_i)
}) %>%
        bind_rows()


### ADD BATCH INFORMATION
across_sample_cor <- across_sample_cor %>%
        left_join(sample_pairs %>%
                          select(sample_A, sample_B, batch_A, batch_B, batch_comparison),
                  by = c("sample_A", "sample_B"))


### MEAN CORRELATIONS BY BATCH
batch_heatmap_mean <- across_sample_cor %>%
        group_by(batch_comparison, celltype_A, celltype_B) %>%
        summarise(mean_correlation = mean(correlation, na.rm = TRUE),
                  n_comparisons = sum(!is.na(correlation)),
                  .groups = "drop")

batch_heatmap_mean


### PLOT WITHIN- AND ACROSS-BATCH CORRELATIONS
# Set cell-type order
celltype_order <- sort(unique(c(batch_heatmap_mean$celltype_A,
                                batch_heatmap_mean$celltype_B)))

batch_heatmap_mean <- batch_heatmap_mean %>%
        mutate(batch_comparison = factor(batch_comparison,
                                         levels = c("Within batch", "Across batch")),
               celltype_A = factor(celltype_A, levels = rev(celltype_order)),
               celltype_B = factor(celltype_B, levels = celltype_order))

# Plot
p_batch <- ggplot(batch_heatmap_mean,
                  aes(x = celltype_B, y = celltype_A,
                      fill = mean_correlation)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.3f", mean_correlation)), size = 2) +
        scale_fill_viridis_c(name = "Mean Pearson r",
                             limits = c(0, 1)) +
        facet_wrap(~batch_comparison, nrow = 1) +
        coord_equal() +
        theme_classic() +
        labs(x = "Cell type B",
             y = "Cell type A") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 10),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"))

p_batch


### WITHIN-BATCH COR BY XENIUM RUN WITH LABELS FOR MIXED VS NON_MIXED RUN
# Get one ILAE score per sample
sample_ILAE <- meta_batch %>%
        transmute(sample = as.character(.data[[sample_col]]),
                  Xenium_run = as.character(.data[[batch_col]]),
                  ILAE_score = as.character(ILAE_score)) %>%
        distinct()

# Determine whether each run contains the same or mixed ILAE scores
run_ILAE <- sample_ILAE %>%
        group_by(Xenium_run) %>%
        summarise(n_ILAE = n_distinct(ILAE_score),
                  ILAE_composition = if_else(n_ILAE == 1, "same ILAE", "mixed"),
                  .groups = "drop") %>%
        mutate(run_label = paste0("Run ", Xenium_run, " (", ILAE_composition, ")"))

# Add labels to heatmap data
within_run_mean <- across_sample_cor %>%
        filter(batch_comparison == "Within batch") %>%
        group_by(batch_A, celltype_A, celltype_B) %>%
        summarise(mean_correlation = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
        rename(Xenium_run = batch_A)
within_run_mean <- within_run_mean %>%
        left_join(run_ILAE %>%
                          select(Xenium_run, run_label),
                  by = "Xenium_run")

# Plot one correlation matrix per Xenium run
p_within_run <- ggplot(within_run_mean,
                       aes(x = celltype_B, y = celltype_A,
                           fill = mean_correlation)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.3f", mean_correlation)), size = 2) +
        scale_fill_viridis_c(name = "Mean Pearson r",
                             limits = c(0, 1)) +
        facet_wrap(~run_label, ncol = 5) +
        coord_equal() +
        theme_classic() +
        labs(x = "Cell type B",
             y = "Cell type A",
             title = "Within-batch correlations by Xenium run") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 10),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"))

p_within_run

### BATCH CORRELATION HEATMAPS
# Add ILAE composition of the Xenium run to within-batch comparisons
batch_cor_plot <- across_sample_cor %>%
        left_join(run_ILAE %>%
                          select(Xenium_run, ILAE_composition),
                  by = c("batch_A" = "Xenium_run"))

# Across-batch mean correlations
cor_across <- batch_cor_plot %>%
        filter(batch_comparison == "Across batch") %>%
        group_by(celltype_A, celltype_B) %>%
        summarise(mean_correlation = mean(correlation, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(comparison = "Across batch")

# All within-batch mean correlations
cor_within <- batch_cor_plot %>%
        filter(batch_comparison == "Within batch") %>%
        group_by(celltype_A, celltype_B) %>%
        summarise(mean_correlation = mean(correlation, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(comparison = "Within batch")

# Combine the four matrices
batch_cor_heatmap <- bind_rows(cor_across,
                               cor_within) %>%
        mutate(comparison = factor(comparison,
                                   levels = c("Across batch",
                                              "Within batch")),
               celltype_A = factor(celltype_A, levels = rev(celltype_order)),
               celltype_B = factor(celltype_B, levels = celltype_order))

# Get minimum mean correlation across all four matrices
min_cor <- min(batch_cor_heatmap$mean_correlation, na.rm = TRUE)
min_cor

# Plot all four matrices using the same color scale
p_batch_comparison <- ggplot(batch_cor_heatmap,
                             aes(x = celltype_B, y = celltype_A,
                                 fill = mean_correlation)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.3f", mean_correlation)), size = 2) +
        scale_fill_viridis_c(name = "Mean Pearson r",
                             limits = c(min_cor, 1))+
        facet_wrap(~comparison, ncol = 2) +
        coord_equal() +
        theme_classic() +
        labs(x = "Cell type (sample 2)",
             y = "Cell type (sample 1)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 10),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"))

p_batch_comparison




