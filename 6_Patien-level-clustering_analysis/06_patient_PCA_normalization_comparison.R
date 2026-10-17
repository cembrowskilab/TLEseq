################################################################################
#
# 06_patient_PCA_normalization_comparison.R
#
# Compare patient LogN/SCT PCA and batch-corrected Harmony embeddings.
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
        library(harmony)
        library(patchwork)
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

# Fraction of variance across Harmony dimensions (not PCA variance explained).
variance_harmony <- function(x) {
        dimension_variance <- apply(x, 2, var)
        dimension_variance / sum(dimension_variance) * 100
}

# Align PC1 signs to the patient-matched LogN reference.
align_pc1 <- function(embedding, reference_pc1, common_patients) {
        r <- cor(reference_pc1, embedding[common_patients, 1], method = "pearson")
        if (is.finite(r) && r < 0) embedding[, 1] <- -embedding[, 1]
        embedding
}

# Regress patient PC1 on ILAE score, age and median transcript count.
run_pc1_regression <- function(embedding, method_name, patient_metadata) {
        df <- data.frame(patient = rownames(embedding), PC1 = embedding[, 1]) %>%
                left_join(patient_metadata, by = "patient") %>%
                mutate(patient_age = as.numeric(as.character(patient_age)),
                       median_nCount_Xenium = as.numeric(as.character(median_nCount_Xenium)),
                       ILAE_score = factor(ILAE_score, levels = c("P", "0", "1", "2", "3"))) %>%
                filter(!is.na(PC1), !is.na(ILAE_score), !is.na(patient_age), !is.na(median_nCount_Xenium)) %>%
                mutate(patient_age_centered = patient_age - mean(patient_age),
                       nCount_centered = median_nCount_Xenium - mean(median_nCount_Xenium))
        model <- lm(PC1 ~ ILAE_score + patient_age_centered + nCount_centered, data = df)
        df$PC1_adjusted <- df$PC1 - coef(model)["patient_age_centered"] * df$patient_age_centered -
                coef(model)["nCount_centered"] * df$nCount_centered
        df$method <- method_name
        list(data = df, model = model)
}

# Return emmeans and multiple-comparison-adjusted ILAE vs postmortem contrasts.
get_emmeans <- function(result, method_name) {
        means <- emmeans::emmeans(result$model, ~ ILAE_score)
        pairwise <- as.data.frame(emmeans::contrast(means, method = "trt.vs.ctrl", ref = "P", adjust = "mvt"))
        pairwise$method <- method_name
        pairwise$significance <- ifelse(pairwise$p.value < 0.001, "***",
                                ifelse(pairwise$p.value < 0.01, "**",
                                ifelse(pairwise$p.value < 0.05, "*", "ns")))
        list(emmeans = as.data.frame(means), pairwise = pairwise)
}

### 2) User-editable paths and load inputs -------------------------------------
cfg <- list(megadata_rds = NA_character_, out_results_dir = NA_character_,
            batch_column = "Xenium_run", ilae_order = c("P", "0", "3", "2", "1"),
            distance_color_min = 0.2)
if (is.na(cfg$out_results_dir)) stop("Set cfg$out_results_dir")
dir.create(cfg$out_results_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- cfg$out_results_dir
method_order <- c("LogN", "LogN + Harmony", "SCT", "SCT + Harmony")
ilae_order <- cfg$ilae_order
batch_column <- cfg$batch_column

### 3) Palettes -----------------------------------------------------------------

pal_ILAE_score <- c(
        "P" = "#C67BFF",
        "0" = "#F8766C",
        "1" = "#7BAD00",
        "2" = "#00C0C5",
        "3" = "#E69F00"
)

### 4) Load data and construct patient profiles -------------------------------
dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")
DefaultAssay(dat) <- "Xenium"
dat2 <- JoinLayers(dat, assay = "Xenium")

# Extract raw Xenium counts and metadata
counts <- GetAssayData(dat2, assay = "Xenium", layer = "counts")
meta <- dat2@meta.data
if (!cfg$batch_column %in% names(meta)) stop("Missing batch column: ", cfg$batch_column)

# Patient-level metadata, including median nCount from original cells
patient_metadata <- meta %>%
        mutate(patient = as.character(orig.ident)) %>%
        group_by(patient) %>%
        summarise(ILAE_score = first(ILAE_score),
                  sex = first(sex),
                  patient_age = first(patient_age),
                  median_nCount_Xenium = median(nCount_Xenium, na.rm = TRUE),
                  Xenium_run = first(.data[[cfg$batch_column]]),
                  .groups = "drop")

# Mean raw counts within patient x cell type
# Define patient x cell-type groups
patient_celltype <- interaction(meta$orig.ident, meta$celltype, drop = TRUE)

# Mean raw expression across cells within each patient x cell-type group
group_means <- sapply(levels(patient_celltype), function(group) {
        cells <- which(patient_celltype == group)
        Matrix::rowMeans(counts[, cells, drop = FALSE])
})

# Keep patient identity for each resulting column
group_patients <- as.character(meta$orig.ident[
        match(levels(patient_celltype), patient_celltype)
])

# Sum cell-type means for each patient to generate one expression profile per patient
patient_ids <- unique(group_patients)

patient_profiles <- sapply(patient_ids, function(patient) {
        Matrix::rowSums(group_means[, group_patients == patient, drop = FALSE])
})

patient_profiles <- as.matrix(patient_profiles)
colnames(patient_profiles) <- patient_ids


### 5) LogN and SCT PCA (w/o Harmony)-------------------------------------------
# Log normalization
if (any(colSums(patient_profiles) <= 0)) stop("Zero-sum patient profile")
normalized_matrix <- sweep(patient_profiles, 2, colSums(patient_profiles), "/") * 10000
log_normalized_matrix <- log1p(normalized_matrix)

gene_variance_log <- apply(log_normalized_matrix, 1, var)
genes_log <- names(gene_variance_log[is.finite(gene_variance_log) & gene_variance_log > 0])

cat("LogN genes available:", nrow(log_normalized_matrix), "\n")
cat("LogN genes used:", length(genes_log), "\n")

if (length(genes_log) < 2) stop("Too few variable genes for LogN PCA")
pca_log <- prcomp(t(log_normalized_matrix[genes_log,,drop = FALSE]),
                  center = TRUE, scale. = TRUE)

variance_log <- (pca_log$sdev^2 / sum(pca_log$sdev^2)) * 100

# Sct Normalization
DefaultAssay(dat2) <- "Xenium"
dat2 <- SCTransform(dat2, assay = "Xenium", new.assay.name = "SCT",
                    vst.flavor = "v2", return.only.var.genes = FALSE, verbose = TRUE)

# Sct Group Identifiers
# Reuse the patient x cell-type grouping for SCT
# Create group IDs if not already present
if (!"patient_celltype_group" %in% colnames(dat2@meta.data)) {
        dat2$patient_celltype_group <- as.character(patient_celltype)
}

# Recalculate using the defined groups
group_means_sct <- AverageExpression(dat2, assays = "SCT",
                                     group.by = "patient_celltype_group",
                                     layer = "data", return.seurat = FALSE,
                                     verbose = FALSE)[["SCT"]]

# Match group columns back to patients
sct_group_lookup <- data.frame(
        group = levels(patient_celltype),
        patient = group_patients
)

# Match Seurat's sanitized group names
sct_group_lookup$group_clean <- gsub("[^[:alnum:]_]", ".", sct_group_lookup$group)

group_patients_sct <- sct_group_lookup$patient[
        match(colnames(group_means_sct), sct_group_lookup$group_clean)
]

if (any(is.na(group_patients_sct))) stop("Could not match SCT groups to patients.")

# Sum Sct Cell-Type Means Within Each Patient
patient_profiles_sct <- sapply(patient_ids, function(patient) {
        Matrix::rowSums(group_means_sct[, group_patients_sct == patient, drop = FALSE])
})

patient_profiles_sct <- as.matrix(patient_profiles_sct)
colnames(patient_profiles_sct) <- patient_ids

# Sct Pca
gene_variance_sct <- apply(patient_profiles_sct, 1, var)
genes_sct <- names(gene_variance_sct[is.finite(gene_variance_sct) & gene_variance_sct > 0])

cat("SCT genes available:", nrow(patient_profiles_sct), "\n")
cat("SCT genes used:", length(genes_sct), "\n")

if (length(genes_sct) < 2) stop("Too few variable genes for SCT PCA")

# Run PCA but without scaling, sicne SCTransform() already performs scaling
pca_sct <- prcomp(t(patient_profiles_sct[genes_sct,,drop = FALSE]),
                  center = TRUE, scale. = FALSE)

variance_sct <- (pca_sct$sdev^2 / sum(pca_sct$sdev^2)) * 100

# Harmony Integration
harmony_metadata_log <- patient_metadata[match(rownames(pca_log$x), patient_metadata$patient), ]
harmony_metadata_sct <- patient_metadata[match(rownames(pca_sct$x), patient_metadata$patient), ]

harmony_log <- harmony::HarmonyMatrix(pca_log$x, harmony_metadata_log,
                             vars_use = batch_column, do_pca = FALSE, verbose = TRUE)

harmony_sct <- harmony::HarmonyMatrix(pca_sct$x, harmony_metadata_sct,
                             vars_use = batch_column, do_pca = FALSE, verbose = TRUE)

rownames(harmony_log) <- rownames(pca_log$x)
rownames(harmony_sct) <- rownames(pca_sct$x)

# Variance Explained
# Original PCA variance explained
variance_log <- (pca_log$sdev^2 / sum(pca_log$sdev^2)) * 100
variance_sct <- (pca_sct$sdev^2 / sum(pca_sct$sdev^2)) * 100

# Relative variance across Harmony dimensions


variance_log_harmony <- variance_harmony(harmony_log)
variance_sct_harmony <- variance_harmony(harmony_sct)

# Combine variance information
variance_list <- list("LogN" = variance_log,
                      "LogN + Harmony" = variance_log_harmony,
                      "SCT" = variance_sct,
                      "SCT + Harmony" = variance_sct_harmony)

# Print first two dimensions
lapply(variance_list, function(x) round(x[1:2], 2))

# Align Pc1 Directions
# Use LogN PC1 as the reference
common_patients <- Reduce(intersect, list(rownames(pca_log$x), rownames(pca_sct$x),
                                          rownames(harmony_log), rownames(harmony_sct)))

reference_pc1 <- pca_log$x[common_patients, 1]



pca_sct$x <- align_pc1(pca_sct$x, reference_pc1, common_patients)
harmony_log <- align_pc1(harmony_log, reference_pc1, common_patients)
harmony_sct <- align_pc1(harmony_sct, reference_pc1, common_patients)

# Prepare Embeddings
embeddings <- list("LogN" = pca_log$x,
                   "LogN + Harmony" = harmony_log,
                   "SCT" = pca_sct$x,
                   "SCT + Harmony" = harmony_sct)

### 6) Compare PCA and Harmony embeddings ------------------------------------
# Construct long and wide patient-level score tables (absent from the source script).
embedding_df <- dplyr::bind_rows(lapply(names(embeddings), function(method_name) {
        emb <- embeddings[[method_name]]
        data.frame(patient = rownames(emb), PC1 = emb[, 1], PC2 = emb[, 2], method = method_name) %>%
                left_join(patient_metadata %>% select(patient, ILAE_score), by = "patient")
}))
pc1_wide <- embedding_df %>% select(patient, ILAE_score, method, PC1) %>%
        tidyr::pivot_wider(names_from = method, values_from = PC1)

p_embeddings_list <- list()

for (method in method_order) {
        plot_df <- embedding_df[embedding_df$method == method, , drop = FALSE]
        var_exp <- variance_list[[method]]
        is_harmony <- grepl("Harmony", method)
        
        x_label <- sprintf("%s (%.1f%%)", ifelse(is_harmony, "Harmony 1", "PC1"), var_exp[1])
        y_label <- sprintf("%s (%.1f%%)", ifelse(is_harmony, "Harmony 2", "PC2"), var_exp[2])
        
        p_embeddings_list[[method]] <- ggplot(plot_df, aes(PC1, PC2, color = ILAE_score)) +
                geom_point(size = 1) +
                coord_fixed() +
                scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
                theme_classic() +
                theme(legend.position = 'none')+
                labs(title = method, x = x_label, y = y_label, color = "ILAE score")
}

p_embeddings <- wrap_plots(p_embeddings_list, ncol = 4)

p_embeddings


# Plot pairwise pc1 correlations
pc1_comparisons <- combn(method_order, 2, simplify = FALSE)
pc1_scatterplots <- list()

for (comparison in pc1_comparisons) {
        method_A <- comparison[1]
        method_B <- comparison[2]
        
        plot_df <- pc1_wide %>%
                transmute(patient, ILAE_score,
                          PC1_A = .data[[method_A]], PC1_B = .data[[method_B]])
        
        test <- cor.test(plot_df$PC1_A, plot_df$PC1_B, method = "spearman", exact = FALSE)
        
        p <- ggplot(plot_df, aes(PC1_A, PC1_B, color = ILAE_score)) +
                geom_point(size = 1) +
                annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
                         label = sprintf("rho = %.3f\nP = %.2e",
                                         unname(test$estimate), test$p.value)) +
                scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
                theme_classic() +
                coord_fixed()+
                theme(legend.position = 'none')+
                labs(x = paste0(method_A, " PC1"), y = paste0(method_B, " PC1"))
        
        pc1_scatterplots[[paste(method_A, method_B, sep = "_")]] <- p
}

p_pc1_comparisons <- wrap_plots(pc1_scatterplots, ncol = 6)
p_pc1_comparisons

# Adjust spearman p values
pc1_correlation_results <- bind_rows(lapply(pc1_comparisons, function(comparison) {
        test <- cor.test(pc1_wide[[comparison[1]]], pc1_wide[[comparison[2]]],
                         method = "spearman", exact = FALSE)
        data.frame(Method_A = comparison[1], Method_B = comparison[2],
                   rho = unname(test$estimate), p.value = test$p.value)
})) %>%
        mutate(p.adj = p.adjust(p.value, method = "BH"))

pc1_correlation_results

# Distance correlation: pc1
embeddings_PC1 <- lapply(embeddings, function(x) x[common_patients,1,drop = FALSE])
distance_vectors_PC1 <- lapply(embeddings_PC1, function(x) as.vector(dist(x)))

distance_cor_matrix_PC1 <- cor(do.call(cbind, distance_vectors_PC1), method = "pearson")
distance_cor_matrix_PC1

# Plot pc1 distance correlation heatmap
distance_cor_df_PC1 <- as.data.frame(as.table(distance_cor_matrix_PC1))
colnames(distance_cor_df_PC1) <- c("Method_A", "Method_B", "r")

distance_cor_df_PC1$Method_A <- factor(distance_cor_df_PC1$Method_A, levels = method_order)
distance_cor_df_PC1$Method_B <- factor(distance_cor_df_PC1$Method_B, levels = method_order)

p_distance_PC1 <- ggplot(distance_cor_df_PC1, aes(Method_B, Method_A, fill = r)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.3f", r)), size = 3.5) +
        scale_fill_viridis_c(limits = c(cfg$distance_color_min, 1), option = "viridis") +
        coord_fixed() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "PC1 patient-distance correlation", x = NULL, y = NULL, fill = "Pearson r")

p_distance_PC1

### 7) Run regression model for logn and sct (no Harmony) ----------------------

# Regression: adjust pc1 for age and ncount
# Run regression for logn and sct
reg_methods <- c("LogN", "SCT")

reg_results <- lapply(reg_methods, function(method) {
        run_pc1_regression(embeddings[[method]], method, patient_metadata)
})
names(reg_results) <- reg_methods

for (method in names(reg_results)) {
        cat("\n", method, "\n")
        print(summary(reg_results[[method]]$model))
        print(drop1(reg_results[[method]]$model, test = "F"))
}

# Emmeans: compare each ilae group against p
emm_results <- lapply(names(reg_results), function(method) {
        get_emmeans(reg_results[[method]], method)
})
names(emm_results) <- names(reg_results)

pairwise_results <- bind_rows(lapply(emm_results, function(x) x$pairwise))
pairwise_results %>% select(method, contrast, estimate, SE, df, t.ratio, p.value, significance)

# Prepare adjusted patient scores and group means
adjusted_pc1 <- bind_rows(lapply(reg_results, function(x) x$data))

pc1_adjusted_means <- bind_rows(lapply(names(emm_results), function(method) {
        emm_results[[method]]$emmeans %>% mutate(method = method)
}))

adjusted_pc1 <- adjusted_pc1 %>%
        mutate(method = factor(method, levels = reg_methods),
               ILAE_score = factor(ILAE_score, levels = ilae_order))

pc1_adjusted_means <- pc1_adjusted_means %>%
        mutate(method = factor(method, levels = reg_methods),
               ILAE_score = factor(ILAE_score, levels = ilae_order))

# Plot adjusted pc1 scores and emmeans
p_adjusted_PC1 <- ggplot(adjusted_pc1, aes(x = ILAE_score, y = PC1_adjusted)) +
        geom_jitter(aes(color = ILAE_score), width = 0.2, size = 2) +
        geom_errorbar(data = pc1_adjusted_means, aes(x = ILAE_score, ymin = lower.CL, ymax = upper.CL),
                      width = 0.12, color = "black", inherit.aes = FALSE) +
        geom_point(data = pc1_adjusted_means, aes(x = ILAE_score, y = emmean),
                   size = 3, shape = 21, fill = "white", color = "black", inherit.aes = FALSE) +
        facet_wrap(~method, scales = "free_y", nrow = 1) +
        scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
        scale_x_discrete(limits = ilae_order) +
        theme_classic() +
        theme(legend.position = "none", strip.background = element_blank()) +
        labs(x = "ILAE score", y = "Adjusted PC1 score")

p_adjusted_PC1

# Save results
write.csv(pairwise_results, file.path(output_dir, "PC1_pairwise_MVT.csv"), row.names = FALSE)
write.csv(pc1_adjusted_means, file.path(output_dir, "PC1_adjusted_means.csv"), row.names = FALSE)
write.csv(adjusted_pc1, file.path(output_dir, "PC1_adjusted_patients.csv"), row.names = FALSE)

# Save reusable embeddings and model outputs.
saveRDS(list(embeddings = embeddings, patient_metadata = patient_metadata,
             variance_list = variance_list, regression = reg_results),
        file.path(output_dir, "normalization_comparison_checkpoint.rds"))
