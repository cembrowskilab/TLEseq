################################################################################
#
# 04_patient_PCA_by_celltype.R
#
# Patient-level PCA and ILAE regression separately for each broad cell type.
#
# Original analysis: Larissa Kraus (September–October 2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
################################################################################

### 0) Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(Matrix)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(emmeans)
        library(broom)
})

### PACKAGE CONFLICTS
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag, dplyr::count)

options(future.globals.maxSize = 25000 * 1024^2)

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

### 2) User-editable paths and load inputs -------------------------------------
cfg <- list(megadata_rds = NA_character_,
            out_results_dir = NA_character_,
            patient_column = "orig.ident", celltype_column = "celltype", assay = "Xenium",
            ilae_order = c("P", "0", "1", "2", "3"))

if (is.na(cfg$out_results_dir)) stop("Set cfg$out_results_dir")
dir.create(cfg$out_results_dir, recursive = TRUE, showWarnings = FALSE)

### 3) Palettes -----------------------------------------------------------------

pal_ILAE_score <- c(
        "P" = "#C67BFF",
        "0" = "#F8766C",
        "1" = "#7BAD00",
        "2" = "#00C0C5",
        "3" = "#E69F00"
)

### 4) Load and prepare MegaData -----------------------------------------------
dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")

DefaultAssay(dat) <- cfg$assay
dat_joinLayers <- JoinLayers(dat, assay = cfg$assay)
counts <- GetAssayData(dat_joinLayers, assay = cfg$assay, layer = "counts")
meta <- dat_joinLayers@meta.data
meta$ILAE_score <- factor(meta$ILAE_score, levels = cfg$ilae_order)

### 5) PCA and adjusted PC1 analysis for each broad cell type -------------------
# Patient profiles are mean raw counts across cells of the same type, followed
# by library-size normalization, log1p transformation and scaled PCA.
patient_column <- cfg$patient_column
celltype_column <- cfg$celltype_column
celltypes <- sort(unique(as.character(meta[[celltype_column]])))
celltypes <- celltypes[!is.na(celltypes) & nzchar(celltypes)]
output_dir <- cfg$out_results_dir

for (current_celltype in celltypes) {
        
        cat("\n\n############################################################\n")
        cat("RUNNING:", current_celltype, "\n")
        cat("############################################################\n")
        
        current_output_dir <- file.path(output_dir, gsub("[^A-Za-z0-9_-]", "_", current_celltype))
        dir.create(current_output_dir, recursive = TRUE, showWarnings = FALSE)
        
        
        # Subset current cell type
        cells_current <- which(meta[[celltype_column]] == current_celltype)
        
        counts_ct <- counts[, cells_current, drop = FALSE]
        meta_ct <- meta[cells_current, , drop = FALSE]
        
        cat("Cells:", ncol(counts_ct), "\n")
        cat("Patients:", length(unique(meta_ct[[patient_column]])), "\n")
        
        
        # Patient metadata
        # Median nCount is calculated from the original cells of this cell type
        patient_metadata <- meta_ct %>%
                mutate(patient = as.character(.data[[patient_column]])) %>%
                group_by(patient) %>%
                summarise(ILAE_score = first(ILAE_score),
                          sex = first(sex),
                          patient_age = first(patient_age),
                          median_nCount_Xenium = median(nCount_Xenium, na.rm = TRUE),
                          n_cells = n(),
                          .groups = "drop") %>%
                mutate(patient_age = as.numeric(as.character(patient_age)))
        
        
        # Mean raw counts per patient
        # Calculate mean raw expression across cells for each patient
        patient_labels <- as.character(meta_ct[[patient_column]])
        patient_ids <- unique(patient_labels)
        
        patient_profiles <- sapply(patient_ids, function(patient) {
                Matrix::rowMeans(counts_ct[, patient_labels == patient, drop = FALSE])
        })
        
        patient_profiles <- as.matrix(patient_profiles)
        colnames(patient_profiles) <- patient_ids
        
        
        # Log normalization
        # Library-size normalize each patient mean profile and log-transform
        if (any(colSums(patient_profiles) <= 0)) stop("Zero-sum patient profile in ", current_celltype)
        normalized_matrix <- sweep(patient_profiles, 2,
                                   Matrix::colSums(patient_profiles), "/") * 10000
        
        log_normalized_matrix <- log1p(normalized_matrix)
        
        
        # PCA
        # Remove genes with no variance across patients
        gene_variance <- apply(log_normalized_matrix, 1, var)
        genes_for_pca <- names(gene_variance[is.finite(gene_variance) & gene_variance > 0])
        
        if (length(genes_for_pca) < 2 || ncol(patient_profiles) < 3) stop("Insufficient genes/patients for PCA in ", current_celltype)
        # Run PCA across patients
        pca <- prcomp(t(log_normalized_matrix[genes_for_pca, , drop = FALSE]),
                      center = TRUE, scale. = TRUE)
        
        variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
        
        # Combine PCA scores with patient metadata
        pca_df <- data.frame(patient = rownames(pca$x),
                             PC1 = pca$x[, 1],
                             PC2 = pca$x[, 2]) %>%
                left_join(patient_metadata, by = "patient")
        
        
        # PCA plots
        
        # PCA colored by ILAE score
        pca_ilae_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = ILAE_score)) +
                geom_point(size = 3) +
                coord_fixed() +
                scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
                theme_classic() +
                theme(legend.position = "bottom") +
                labs(title = current_celltype,
                     x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
                     y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
                     color = "ILAE score")
        
        # PC1 gene loadings
        # Extract PC1 gene loadings
        pc1_loadings <- data.frame(gene = rownames(pca$rotation),
                                   loading = pca$rotation[, "PC1"])
        
        # Select 20 strongest negative and positive PC1 loadings
        top_pc1_genes <- bind_rows(
                pc1_loadings %>% arrange(loading) %>% slice_head(n = 20),
                pc1_loadings %>% arrange(desc(loading)) %>% slice_head(n = 20)
        ) %>%
                arrange(loading) %>%
                mutate(gene = factor(gene, levels = gene),
                       sign = ifelse(loading > 0, "positive", "negative"))
        
        # Regression adjusting for age and transcript count
        # Center continuous covariates
        pc1_regression_df <- pca_df %>%
                filter(!is.na(PC1),
                       !is.na(ILAE_score),
                       !is.na(patient_age),
                       !is.na(median_nCount_Xenium)) %>%
                mutate(ILAE_score = factor(ILAE_score,
                                           levels = c("P", "0", "1", "2", "3")),
                       patient_age_centered = patient_age - mean(patient_age),
                       patient_count_centered = median_nCount_Xenium -
                               mean(median_nCount_Xenium))
        
        # Model PC1 while adjusting for age and median nCount
        pc1_regression_model <- lm(
                PC1 ~ ILAE_score + patient_age_centered + patient_count_centered,
                data = pc1_regression_df
        )
        
        # Adjusted ILAE group means
        pc1_adjusted_means <- emmeans(pc1_regression_model, ~ ILAE_score)
        
        # Comparisons against postmortem controls
        pc1_vs_P <- contrast(pc1_adjusted_means,
                             method = "trt.vs.ctrl",
                             ref = "P",
                             adjust = "mvt")
        
        pc1_adjusted_means_df <- as.data.frame(pc1_adjusted_means)
        
        
        # Adjusted PC1 scores
        # Adjust PC1 to mean age and mean nCount while retaining residual variation
        adjustment_data <- pc1_regression_df %>%
                mutate(patient_age_centered = 0,
                       patient_count_centered = 0)
        
        pc1_regression_df$PC1_adjusted <- predict(
                pc1_regression_model,
                newdata = adjustment_data
        ) + residuals(pc1_regression_model)
        
        
        # Adjusted PC1 plot
        pc1_adjusted_plot <- ggplot(pc1_regression_df,
                                    aes(x = ILAE_score, y = PC1_adjusted)) +
                geom_jitter(aes(color = ILAE_score),
                            width = 0.2, size = 2.5) +
                geom_point(data = pc1_adjusted_means_df,
                           aes(x = ILAE_score, y = emmean),
                           size = 3, shape = 21,
                           fill = "white", color = "black",
                           inherit.aes = FALSE) +
                geom_errorbar(data = pc1_adjusted_means_df,
                              aes(x = ILAE_score,
                                  ymin = lower.CL,
                                  ymax = upper.CL),
                              width = 0.12,
                              color = "black",
                              inherit.aes = FALSE) +
                scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
                scale_x_discrete(limits = c("P", "0", "3", "2", "1")) +
                theme_classic() +
                labs(title = current_celltype,
                     x = "ILAE score",
                     y = "Adjusted PC1 score") +
                theme(legend.position = "none")
        
        print(pc1_adjusted_plot)
        
        
        # Save plots
        ggsave(file.path(current_output_dir, "patient_PCA_ILAE.pdf"),
               pca_ilae_plot, width = 4, height = 4)
        
        ggsave(file.path(current_output_dir, "patient_PC1_adjusted.pdf"),
               pc1_adjusted_plot, width = 3, height = 5)
        
        
        # Save tables
        # PCA scores and patient metadata
        write.csv(pca_df,
                  file.path(current_output_dir, "patient_PCA_scores.csv"),
                  row.names = FALSE)
        
        # All PC1 gene loadings
        write.csv(pc1_loadings,
                  file.path(current_output_dir, "PC1_gene_loadings.csv"),
                  row.names = FALSE)
        
        # Adjusted patient PC1 scores
        write.csv(pc1_regression_df,
                  file.path(current_output_dir, "patient_PC1_adjusted_scores.csv"),
                  row.names = FALSE)
        
        # Regression coefficients
        write.csv(broom::tidy(pc1_regression_model),
                  file.path(current_output_dir, "PC1_regression_coefficients.csv"),
                  row.names = FALSE)
        
        # Partial F tests
        write.csv(broom::tidy(drop1(pc1_regression_model, test = "F")),
                  file.path(current_output_dir, "PC1_regression_partial_F_tests.csv"),
                  row.names = FALSE)
        
        # Adjusted ILAE means
        write.csv(pc1_adjusted_means_df,
                  file.path(current_output_dir, "PC1_adjusted_ILAE_means.csv"),
                  row.names = FALSE)
        
        # Adjusted comparisons against postmortem controls
        write.csv(as.data.frame(pc1_vs_P),
                  file.path(current_output_dir, "PC1_ILAE_vs_P.csv"),
                  row.names = FALSE)
}

