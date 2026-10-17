################################################################################
#
# Xenium split-half correlations by neuronal subtype and microglial state (Seurat v5)
#
# Split-half expression reproducibility by ILAE grade.
# this is run on individual seurat files from subset analysis, so those need to 
# be run first before running this script
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

split_half_heatmap <- function(dat,
                               group_col,
                               celltype_col = "celltype",
                               assay_use = "SCT",
                               layer_use = "data",
                               min_cells = 20,
                               group_order = NULL,
                               seed = 123,
                               plot_title = NULL,
                               facet_ncol = 6) {
        
        ### SETUP #####
        expr <- GetAssayData(dat, assay = assay_use, layer = layer_use)
        meta <- dat[[]]
        meta <- meta[colnames(expr), , drop = FALSE]
        
        # Remove cells missing grouping information
        keep_cells <- !is.na(meta[[group_col]]) & !is.na(meta[[celltype_col]])
        meta <- meta[keep_cells, , drop = FALSE]
        expr <- expr[, rownames(meta), drop = FALSE]
        
        # Mean expression helper
        get_mean_expression <- function(cells) {
                Matrix::rowMeans(expr[, cells, drop = FALSE])
        }
        
        ### RANDOM SPLIT #####
        meta$split_group <- paste(meta[[group_col]], meta[[celltype_col]], sep = "___")
        
        set.seed(seed)
        meta$half <- ave(seq_len(nrow(meta)), meta$split_group, FUN = function(x) {
                sample(rep(c("A", "B"), length.out = length(x)))
        })
        
        # Keep groups with sufficient cells
        group_counts <- meta %>%
                count(.data[[group_col]], .data[[celltype_col]], name = "n_cells")
        
        valid_groups <- group_counts %>%
                filter(n_cells >= min_cells) %>%
                transmute(split_group = paste(.data[[group_col]],
                                              .data[[celltype_col]], sep = "___")) %>%
                pull(split_group)
        
        meta$keep <- meta$split_group %in% valid_groups
        
        ### MEAN EXPRESSION #####
        split_ids <- unique(paste(meta$split_group[meta$keep],
                                  meta$half[meta$keep], sep = "___"))
        
        split_expression <- lapply(split_ids, function(id) {
                
                parts <- strsplit(id, "___", fixed = TRUE)[[1]]
                group_i <- parts[1]
                celltype_i <- parts[2]
                half_i <- parts[3]
                
                cells <- rownames(meta)[
                        as.character(meta[[group_col]]) == group_i &
                                as.character(meta[[celltype_col]]) == celltype_i &
                                meta$half == half_i
                ]
                
                data.frame(
                        group = group_i,
                        celltype = celltype_i,
                        half = half_i,
                        gene = rownames(expr),
                        expression = get_mean_expression(cells)
                )
        }) %>%
                bind_rows()
        
        ### WITHIN-CELL-TYPE CORRELATIONS #####
        within_cor <- split_expression %>%
                select(group, celltype, half, gene, expression) %>%
                pivot_wider(names_from = half, values_from = expression) %>%
                group_by(group, celltype) %>%
                summarise(correlation = cor(A, B, method = "pearson"),
                          .groups = "drop")
        
        ### BETWEEN-CELL-TYPE CORRELATIONS #####
        between_cor <- split_expression %>%
                group_split(group) %>%
                lapply(function(x) {
                        
                        group_i <- unique(x$group)
                        
                        mat_A <- x %>%
                                filter(half == "A") %>%
                                select(gene, celltype, expression) %>%
                                pivot_wider(names_from = celltype, values_from = expression) %>%
                                tibble::column_to_rownames("gene") %>%
                                as.matrix()
                        
                        mat_B <- x %>%
                                filter(half == "B") %>%
                                select(gene, celltype, expression) %>%
                                pivot_wider(names_from = celltype, values_from = expression) %>%
                                tibble::column_to_rownames("gene") %>%
                                as.matrix()
                        
                        cor_mat <- cor(mat_A, mat_B, method = "pearson")
                        
                        as.data.frame(as.table(cor_mat)) %>%
                                rename(celltype1 = Var1,
                                       celltype2 = Var2,
                                       correlation = Freq) %>%
                                mutate(group = group_i) %>%
                                filter(celltype1 != celltype2)
                }) %>%
                bind_rows()
        
        ### HEATMAP DATA #####
        heatmap_df <- bind_rows(
                between_cor %>%
                        transmute(group,
                                  celltype_A = as.character(celltype1),
                                  celltype_B = as.character(celltype2),
                                  correlation),
                within_cor %>%
                        transmute(group,
                                  celltype_A = as.character(celltype),
                                  celltype_B = as.character(celltype),
                                  correlation)
        )
        
        # Cell-type order
        celltype_order <- sort(unique(c(heatmap_df$celltype_A,
                                        heatmap_df$celltype_B)))
        
        # Optional group order
        if (is.null(group_order)) {
                group_order <- unique(heatmap_df$group)
        }
        
        heatmap_df <- heatmap_df %>%
                mutate(group = factor(group, levels = group_order),
                       celltype_A = factor(celltype_A, levels = rev(celltype_order)),
                       celltype_B = factor(celltype_B, levels = celltype_order))
        
        ### GROUP HEATMAPS #####
        p_group <- ggplot(heatmap_df,
                          aes(x = celltype_B, y = celltype_A,
                              fill = correlation)) +
                geom_tile() +
                scale_fill_viridis_c(name = "Pearson r", limits = c(0, 1)) +
                facet_wrap(~group, ncol = facet_ncol) +
                coord_equal() +
                theme_classic() +
                labs(x = "Half B",
                     y = "Half A",
                     title = plot_title) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.grid = element_blank(),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "bold"))
        
        ### MEAN HEATMAP #####
        heatmap_mean <- heatmap_df %>%
                group_by(celltype_A, celltype_B) %>%
                summarise(mean_correlation = mean(correlation, na.rm = TRUE),
                          .groups = "drop")
        
        p_mean <- ggplot(heatmap_mean,
                         aes(x = celltype_B, y = celltype_A,
                             fill = mean_correlation)) +
                geom_tile() +
                geom_text(aes(label = sprintf("%.3f", mean_correlation)), size = 2) +
                scale_fill_viridis_c(name = "Mean Pearson r", limits = c(0, 1)) +
                coord_equal() +
                theme_classic() +
                labs(x = "Half B",
                     y = "Half A",
                     title = paste("Mean:", plot_title)) +
                theme(axis.text = element_text(size = 6),
                      axis.title = element_text(size = 10))
        
        ### RETURN EVERYTHING #####
        return(list(
                group_counts = group_counts,
                split_expression = split_expression,
                within_cor = within_cor,
                between_cor = between_cor,
                heatmap_df = heatmap_df,
                heatmap_mean = heatmap_mean,
                p_group = p_group,
                p_mean = p_mean
        ))
}

### 2) User-editable paths and load inputs -------------------------------------
cfg <- list(datEx_rds = NA_character_, 
            datIN_rds = NA_character_, 
            datM_clean_rds = NA_character_,
            out_dir = NA_character_, 
            min_cells = 20, 
            seed = 123)
datEx <- load_if_missing("datEx", cfg$datEx_rds, "Seurat")
datIN <- load_if_missing("datIN", cfg$datIN_rds, "Seurat")
datM_clean <- load_if_missing("datM_clean", cfg$datM_clean_rds, "Seurat")
if (is.na(cfg$out_dir)) stop("Set cfg$out_dir")


### 3) Split-half analyses ----------------------------------------------------
#for ExN subtypes by ILAE
ILAE_cor_datEx <- split_half_heatmap(
        dat = datEx,
        group_col = "ILAE_score",
        celltype_col = "broad_subtype_ExN_LK",
        assay_use = "SCT",
        layer_use = "data",
        group_order = c("P", "0", "1", "2", "3"),
        min_cells = cfg$min_cells, seed = cfg$seed,
        plot_title = "Grouped by ILAE grade",
        facet_ncol = 5)

ILAE_cor_datEx$p_group
ILAE_cor_datEx$p_mean
newExN_plot<- ILAE_cor_datEx$p_mean +
        scale_fill_viridis_c(name = "Pearson r", 
                             limits = c(min(ILAE_cor_datEx$heatmap_mean$mean_correlation), 1))



#for InN subtypes by ILAE
ILAE_cor_datIN <- split_half_heatmap(
        dat = datIN,
        group_col = "ILAE_score",
        celltype_col = "broad_subtype_InhN",
        assay_use = "SCT",
        layer_use = "data",
        group_order = c("P", "0", "1", "2", "3"),
        min_cells = cfg$min_cells, seed = cfg$seed,
        plot_title = "Grouped by ILAE grade",
        facet_ncol = 5)

ILAE_cor_datIN$p_group
ILAE_cor_datIN$p_mean
newInN_plot<- ILAE_cor_datIN$p_mean +
        scale_fill_viridis_c(name = "Pearson r", 
        limits = c(min(ILAE_cor_datIN$heatmap_mean$mean_correlation), 1))


#for datM subtypes by ILAE
ILAE_cor_datM <- split_half_heatmap(
        dat = datM_clean,
        group_col = "ILAE_score",
        celltype_col = "cellstates",
        assay_use = "SCT",
        layer_use = "data",
        group_order = c("P", "0", "1", "2", "3"),
        min_cells = cfg$min_cells, seed = cfg$seed,
        plot_title = "Grouped by ILAE grade",
        facet_ncol = 5)

ILAE_cor_datM$p_group
ILAE_cor_datM$p_mean
newM_plot<- ILAE_cor_datM$p_mean +
        scale_fill_viridis_c(name = "Pearson r", 
                             limits = c(min(ILAE_cor_datM$heatmap_mean$mean_correlation), 1))



# Save separate subtype and cell-state heatmaps.
for (nm in c("ExN", "InN", "Micro")) dir.create(file.path(cfg$out_dir, nm), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(cfg$out_dir, "ExN", "split_half_mean.pdf"), newExN_plot, width = 5, height = 5)
ggsave(file.path(cfg$out_dir, "ExN", "split_half_by_ILAE.pdf"), ILAE_cor_datEx$p_group, width = 9, height = 4)
ggsave(file.path(cfg$out_dir, "InN", "split_half_mean.pdf"), newInN_plot, width = 5, height = 5)
ggsave(file.path(cfg$out_dir, "InN", "split_half_by_ILAE.pdf"), ILAE_cor_datIN$p_group, width = 9, height = 4)
ggsave(file.path(cfg$out_dir, "Micro", "split_half_mean.pdf"), newM_plot, width = 5, height = 5)
ggsave(file.path(cfg$out_dir, "Micro", "split_half_by_ILAE.pdf"), ILAE_cor_datM$p_group, width = 9, height = 4)
