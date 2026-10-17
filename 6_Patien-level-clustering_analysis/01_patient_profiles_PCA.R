################################################################################
#
# 01_patient_profiles_PCA.R
#
# Generate equally weighted cell-type patient profiles, PCA, PC1 loadings and covariate correlations.
#
# Original analysis: Larissa Kraus (October 2026)
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
        library(patchwork)
        library(broom)
})

### PACKAGE CONFLICTS
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)
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

### 2) User-editable paths and load inputs ------------------------------------
# Set paths before running. A loaded 'dat' object will be reused if available.
cfg <- list(megadata_rds = NA_character_,
            out_results_dir = NA_character_,
            celltype_order = c("ExN", "InN", "Astro", "Micro", "Oligo", "OPC", "Vascular", "CSF"),
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
dat$ILAE_score <- factor(dat$ILAE_score, levels = cfg$ilae_order)
dat$celltype <- factor(dat$celltype, levels = cfg$celltype_order)

DefaultAssay(dat) <- "Xenium"
dat_joined <- JoinLayers(dat, assay = "Xenium")

counts <- GetAssayData(dat_joined, assay = "Xenium", layer = "counts")
meta <- dat_joined@meta.data

### 5) Generate patient-level expression profiles ------------------------------

# Extract patient-level metadata
patient_metadata <- meta %>%
        mutate(patient = as.character(orig.ident)) %>%
        group_by(patient) %>%
        summarise(ILAE_score = first(ILAE_score),
                  sex = first(sex),
                  patient_age = first(patient_age),
                  median_nCount_Xenium = median(nCount_Xenium, na.rm = TRUE),
                  .groups = "drop")

# Calculate mean raw expression within each patient x cell type
patient_celltype <- interaction(meta$orig.ident, meta$celltype, drop = TRUE)

celltype_mean_profiles <- sapply(levels(patient_celltype), function(group) {
        cells <- which(patient_celltype == group)
        Matrix::rowMeans(counts[, cells, drop = FALSE])
})

# Identify the patient corresponding to each cell-type mean profile
group_patients <- as.character(meta$orig.ident[match(levels(patient_celltype), patient_celltype)])
patient_ids <- unique(group_patients)

# Aggregate cell-type mean profiles within each patient
patient_profiles <- sapply(patient_ids, function(patient) {
        Matrix::rowSums(celltype_mean_profiles[, group_patients == patient, drop = FALSE])
})
colnames(patient_profiles) <- patient_ids

# Library-size normalize and log-transform patient profiles
normalized_matrix <- sweep(patient_profiles, 2, colSums(patient_profiles), "/") * 10000
log_normalized_matrix <- log1p(normalized_matrix)

### 6) Patient-level PCA -------------------------------------------------------

# Remove genes with no variance across patients
gene_variance <- apply(log_normalized_matrix, 1, var)
genes_for_pca <- names(gene_variance[is.finite(gene_variance) & gene_variance > 0])

# Run PCA across patients
pca <- prcomp(t(log_normalized_matrix[genes_for_pca, , drop = FALSE]),
              center = TRUE, scale. = TRUE)

variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

pca_df <- data.frame(patient = rownames(pca$x),
                     PC1 = pca$x[, 1],
                     PC2 = pca$x[, 2]) %>%
        left_join(patient_metadata, by = "patient") %>%
        mutate(ILAE_score = factor(ILAE_score, levels = cfg$ilae_order))

# PCA colored by ILAE score
pca_ilae_plot <- ggplot(pca_df, aes(PC1, PC2, color = ILAE_score)) +
        geom_point(size = 3) +
        coord_fixed() +
        scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
             y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
             color = "ILAE score")

# PCA colored by patient age
pca_age_plot <- ggplot(pca_df, aes(PC1, PC2, color = patient_age)) +
        geom_point(size = 3) +
        coord_fixed() +
        scale_color_viridis_c() +
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
             y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
             color = "Patient age")

# PCA colored by median transcript count
pca_nCount_plot <- ggplot(pca_df, aes(PC1, PC2, color = median_nCount_Xenium)) +
        geom_point(size = 3) +
        coord_fixed() +
        scale_color_viridis_c() +
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
             y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
             color = "Median nCount_Xenium")

pca_ilae_plot
pca_age_plot
pca_nCount_plot

### 7) PC1 gene loadings ------------------------------------------------------

pc1_loadings <- data.frame(gene = rownames(pca$rotation),
                           loading = pca$rotation[, "PC1"])

top_pc1_genes <- bind_rows(
        pc1_loadings %>% arrange(loading) %>% slice_head(n = 20),
        pc1_loadings %>% arrange(desc(loading)) %>% slice_head(n = 20)
) %>%
        arrange(loading) %>%
        mutate(gene = factor(gene, levels = gene),
               sign = ifelse(loading > 0, "positive", "negative"))

pc1_loading_plot <- ggplot(top_pc1_genes, aes(gene, loading, fill = sign)) +
        geom_col(width = 0.6) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = c("negative" = "#377EB8", "positive" = "#E41A1C")) +
        coord_flip() +
        theme_classic() +
        theme(axis.text = element_text(size = 10, color = "black"),
              legend.position = "none") +
        labs(x = "Gene", y = "PC1 loading")

pc1_loading_plot

### 8) PC1 associations with age and transcript count -------------------------

pc1_age_cor <- cor.test(pca_df$PC1, pca_df$patient_age,
                        method = "pearson", use = "complete.obs")

pc1_age_label <- paste0("r = ", sprintf("%.2f", pc1_age_cor$estimate),
                        "\np = ", format.pval(pc1_age_cor$p.value, digits = 2, eps = 0.001))

pc1_age_cor_plot <- ggplot(pca_df, aes(patient_age, PC1, color = ILAE_score)) +
        geom_point(size = 1.5) +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
        annotate("text", x = Inf, y = Inf, label = pc1_age_label,
                 hjust = 1.1, vjust = 1.1, size = 3.2) +
        scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
        theme_classic() +
        labs(x = "Patient age", y = "PC1 score", color = "ILAE score") +
        theme(legend.position = "bottom")

pc1_nCount_cor <- cor.test(pca_df$PC1, pca_df$median_nCount_Xenium,
                           method = "pearson", use = "complete.obs")

pc1_nCount_label <- paste0("r = ", sprintf("%.2f", pc1_nCount_cor$estimate),
                           "\np = ", format.pval(pc1_nCount_cor$p.value, digits = 2, eps = 0.001))

pc1_nCount_cor_plot <- ggplot(pca_df, aes(median_nCount_Xenium, PC1, color = ILAE_score)) +
        geom_point(size = 1.5) +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
        annotate("text", x = Inf, y = Inf, label = pc1_nCount_label,
                 hjust = 1.1, vjust = 1.1, size = 3.2) +
        scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
        theme_classic() +
        labs(x = "Median nCount_Xenium", y = "PC1 score", color = "ILAE score") +
        theme(legend.position = "bottom")

pc1_age_cor_plot
pc1_nCount_cor_plot
pc1_age_cor
pc1_nCount_cor


### 9) Save patient-level PCA checkpoint ----------------------------------------
# The checkpoint allows regression and LOPO scripts to run without reloading cells.
pca_checkpoint <- list(patient_profiles = patient_profiles, patient_metadata = patient_metadata,
                       pca = pca, pca_df = pca_df, variance_explained = variance_explained,
                       pal_ILAE_score = pal_ILAE_score)
saveRDS(pca_checkpoint, file.path(cfg$out_results_dir, "patient_level_pca_checkpoint.rds"))
write.csv(pca_df, file.path(cfg$out_results_dir, "patient_level_pca_scores.csv"), row.names = FALSE)
write.csv(pc1_loadings, file.path(cfg$out_results_dir, "patient_level_pc1_loadings.csv"), row.names = FALSE)
