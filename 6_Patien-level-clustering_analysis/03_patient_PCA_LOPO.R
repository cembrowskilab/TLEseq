################################################################################
#
# 03_patient_PCA_LOPO.R
#
# Leave-one-patient-out PCA and regression sensitivity analysis, including KCNQ3 loading stability.
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
cfg <- list(pca_checkpoint_rds = NA_character_,
            regression_checkpoint_rds = NA_character_,
            out_results_dir = NA_character_,
            ilae_order = c("P", "0", "1", "2", "3"))

if (is.na(cfg$out_results_dir)) stop("Set cfg$out_results_dir")
dir.create(cfg$out_results_dir, recursive = TRUE, showWarnings = FALSE)
pca_checkpoint <- load_if_missing("pca_checkpoint", cfg$pca_checkpoint_rds, class = "list")

regression_checkpoint <- load_if_missing("regression_checkpoint", cfg$regression_checkpoint_rds, class = "list")

if (!all(c("pca", "patient_profiles", "variance_explained") %in% names(pca_checkpoint))) stop("Invalid PCA checkpoint")
if (!all(c("pc1_regression_model", "pc1_regression_df") %in% names(regression_checkpoint))) stop("Invalid regression checkpoint")
pca <- pca_checkpoint$pca
patient_profiles <- pca_checkpoint$patient_profiles
variance_explained <- pca_checkpoint$variance_explained
pc1_regression_model <- regression_checkpoint$pc1_regression_model
pc1_regression_df <- regression_checkpoint$pc1_regression_df

### 3) Leave-one-patient-out analysis ----------------------------------------
# Test stability of the patient-level PCA and ILAE effects by removing one
# patient at a time and rerunning normalization, PCA and regression

gene_of_interest <- "KCNQ3"

# Use results from the full-data patient-level analysis as reference
full_pc1_loadings <- pca$rotation[, "PC1"]
full_pc1_variance <- variance_explained[1]

full_ilae_coefficients <- broom::tidy(pc1_regression_model, conf.int = TRUE) %>%
        filter(grepl("^ILAE_score", term)) %>%
        select(term, full_estimate = estimate,
               full_conf_low = conf.low,
               full_conf_high = conf.high,
               full_p_value = p.value)

full_partial_F <- broom::tidy(drop1(pc1_regression_model, test = "F")) %>%
        filter(term == "ILAE_score")

# KCNQ3 loading in the full-data PCA
if (!gene_of_interest %in% names(full_pc1_loadings)) {
        stop(paste0(gene_of_interest, " is not present in the PCA."))
}

full_kcnq3_loading <- full_pc1_loadings[gene_of_interest]
full_kcnq3_rank <- rank(-abs(full_pc1_loadings), ties.method = "min")[gene_of_interest]

cat("Full-data PC1 variance:", round(full_pc1_variance, 2), "%\n")
cat("Full-data KCNQ3 PC1 loading:", round(full_kcnq3_loading, 4), "\n")
cat("Full-data KCNQ3 absolute loading rank:", full_kcnq3_rank,
    "of", length(full_pc1_loadings), "\n")

# Prepare result tables
patients_lopo <- colnames(patient_profiles)

lopo_pca_summary <- data.frame()
lopo_ilae_results <- data.frame()

# Remove each patient individually and rerun the patient-level analysis
for (patient_out in patients_lopo) {
        
        cat("\nLOPO:", which(patients_lopo == patient_out), "/", length(patients_lopo),
            "- removing", patient_out, "\n")
        
        # Remove patient from the already constructed patient profiles
        patient_profiles_lopo <- patient_profiles[, colnames(patient_profiles) != patient_out, drop = FALSE]
        
        # Normalize remaining patient profiles
        normalized_lopo <- sweep(patient_profiles_lopo, MARGIN = 2,
                                 STATS = Matrix::colSums(patient_profiles_lopo), FUN = "/") * 10000
        log_normalized_lopo <- log1p(normalized_lopo)
        
        # Remove genes without variance and rerun PCA
        gene_variance_lopo <- apply(log_normalized_lopo, 1, var)
        genes_for_pca_lopo <- names(gene_variance_lopo[
                is.finite(gene_variance_lopo) & gene_variance_lopo > 0])
        
        pca_lopo <- prcomp(t(log_normalized_lopo[genes_for_pca_lopo, , drop = FALSE]),
                           center = TRUE, scale. = TRUE)
        
        variance_lopo <- (pca_lopo$sdev^2 / sum(pca_lopo$sdev^2)) * 100
        lopo_loadings <- pca_lopo$rotation[, "PC1"]
        
        # Align LOPO PC1 direction to the full-data PC1
        common_genes <- intersect(names(full_pc1_loadings), names(lopo_loadings))
        loading_cor_before <- cor(full_pc1_loadings[common_genes],
                                  lopo_loadings[common_genes])
        sign_flip <- loading_cor_before < 0
        
        if (sign_flip) {
                pca_lopo$x[, "PC1"] <- -pca_lopo$x[, "PC1"]
                pca_lopo$rotation[, "PC1"] <- -pca_lopo$rotation[, "PC1"]
                lopo_loadings <- pca_lopo$rotation[, "PC1"]
        }
        
        loading_cor <- cor(full_pc1_loadings[common_genes],
                           lopo_loadings[common_genes])
        
        # Record KCNQ3 loading and absolute loading rank
        kcnq3_loading <- if (gene_of_interest %in% names(lopo_loadings)) {
                lopo_loadings[gene_of_interest]
        } else NA_real_
        
        kcnq3_rank <- if (gene_of_interest %in% names(lopo_loadings)) {
                rank(-abs(lopo_loadings), ties.method = "min")[gene_of_interest]
        } else NA_real_
        
        # Add metadata from the original regression data
        pca_df_lopo <- data.frame(patient = rownames(pca_lopo$x),
                                  PC1 = pca_lopo$x[, "PC1"]) %>%
                left_join(pc1_regression_df %>%
                                  select(patient, ILAE_score, patient_age, median_nCount_Xenium) %>%
                                  distinct(),
                          by = "patient")
        
        # Recenter covariates after patient removal and rerun regression
        regression_df_lopo <- pca_df_lopo %>%
                filter(!is.na(PC1), !is.na(ILAE_score),
                       !is.na(patient_age), !is.na(median_nCount_Xenium)) %>%
                mutate(ILAE_score = factor(ILAE_score, levels = c("P", "0", "1", "2", "3")),
                       patient_age_centered = patient_age - mean(patient_age),
                       median_nCount_centered = median_nCount_Xenium - mean(median_nCount_Xenium))
        
        model_lopo <- lm(PC1 ~ ILAE_score + patient_age_centered + median_nCount_centered,
                         data = regression_df_lopo)
        
        model_tidy_lopo <- broom::tidy(model_lopo, conf.int = TRUE) %>%
                filter(grepl("^ILAE_score", term)) %>%
                left_join(full_ilae_coefficients, by = "term") %>%
                mutate(patient_removed = patient_out,
                       influence = estimate - full_estimate)
        
        lopo_ilae_results <- bind_rows(lopo_ilae_results, model_tidy_lopo)
        
        partial_F_lopo <- broom::tidy(drop1(model_lopo, test = "F")) %>%
                filter(term == "ILAE_score")
        
        # Store PCA and regression stability metrics
        lopo_pca_summary <- bind_rows(
                lopo_pca_summary,
                data.frame(patient_removed = patient_out,
                           ILAE_group_removed = as.character(pc1_regression_df$ILAE_score[
                                   match(patient_out, pc1_regression_df$patient)]),
                           n_patients = ncol(patient_profiles_lopo),
                           n_genes = length(genes_for_pca_lopo),
                           PC1_variance = variance_lopo[1],
                           PC1_variance_change = variance_lopo[1] - full_pc1_variance,
                           PC1_loading_correlation = loading_cor,
                           PC1_sign_flipped = sign_flip,
                           KCNQ3_loading = kcnq3_loading,
                           KCNQ3_loading_change = kcnq3_loading - full_kcnq3_loading,
                           KCNQ3_loading_rank = kcnq3_rank,
                           ILAE_partial_F = partial_F_lopo$statistic,
                           ILAE_partial_F_p = partial_F_lopo$p.value)
        )
}

# Add readable ILAE comparison labels
lopo_ilae_results <- lopo_ilae_results %>%
        mutate(comparison = recode(term,
                                   "ILAE_score0" = "0 - P",
                                   "ILAE_score1" = "1 - P",
                                   "ILAE_score2" = "2 - P",
                                   "ILAE_score3" = "3 - P"))

# Plot PC1 structural stability
lopo_loading_plot <- ggplot(lopo_pca_summary,
                            aes(x = reorder(patient_removed, PC1_loading_correlation),
                                y = PC1_loading_correlation)) +
        geom_point(size = 2.5) +
        geom_text(aes(label = sprintf("%.3f", PC1_loading_correlation)),
                  hjust = 1.5, size = 3) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
        scale_y_continuous(limits = c(0, 1.03),
                           breaks = seq(0, 1, 0.2),
                           expand = expansion(mult = c(0, 0.01))) +
        coord_flip() +
        theme_classic() +
        labs(x = "Patient removed",
             y = "Correlation with full-data PC1 loadings",
             title = "Leave-one-patient-out PC1 stability")

lopo_loading_plot

# Prepare full-data ILAE estimates for comparison
full_ilae_plot <- full_ilae_coefficients %>%
        mutate(comparison = recode(term,
                                   "ILAE_score0" = "0 - P",
                                   "ILAE_score1" = "1 - P",
                                   "ILAE_score2" = "2 - P",
                                   "ILAE_score3" = "3 - P"),
               patient_removed = "Full dataset") %>%
        transmute(patient_removed, comparison,
                  estimate = full_estimate,
                  conf.low = full_conf_low,
                  conf.high = full_conf_high)

patient_order <- c("L10", "L5", "M70", "P28", "P60", "P71", "P51", "M56",
                   "B1", "F1", "H1", "I1", "J1", "L20", "D11", "D8",
                   "C1", "C2", "E008", "L14", "L15", "L18", "L1", "L21",
                   "L3", "D12", "D13", "E015", "L2", "D9", "D14")

# Plot KCNQ3 PC1 loading stability
kcnq3_loading_mean <- mean(lopo_pca_summary$KCNQ3_loading, na.rm = TRUE)
kcnq3_loading_sd <- sd(lopo_pca_summary$KCNQ3_loading, na.rm = TRUE)

loading_order <- lopo_pca_summary %>%
        arrange(KCNQ3_loading) %>%
        pull(patient_removed)

kcnq3_loading_plot_df <- lopo_pca_summary %>%
        select(patient_removed, KCNQ3_loading) %>%
        bind_rows(data.frame(patient_removed = "Full dataset",
                             KCNQ3_loading = full_kcnq3_loading)) %>%
        mutate(patient_removed = factor(patient_removed,
                                        levels = c("Full dataset", loading_order)))

lopo_kcnq3_loading_plot <- ggplot(kcnq3_loading_plot_df,
                                  aes(x = KCNQ3_loading, y = patient_removed)) +
        annotate("rect",
                 xmin = kcnq3_loading_mean - kcnq3_loading_sd,
                 xmax = kcnq3_loading_mean + kcnq3_loading_sd,
                 ymin = -Inf, ymax = Inf,
                 fill = "red", alpha = 0.12) +
        geom_vline(xintercept = kcnq3_loading_mean,
                   color = "red", linetype = "dotted") +
        geom_point(aes(color = patient_removed == "Full dataset"), size = 2.5) +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                           guide = "none") +
        theme_classic() +
        labs(x = "KCNQ3 PC1 loading",
             y = "Patient removed",
             title = "Leave-one-patient-out KCNQ3 loading",
             subtitle = "Red = full dataset; dotted line and shading = LOPO mean ± SD")

lopo_kcnq3_loading_plot

# Plot KCNQ3 absolute PC1 loading-rank stability
kcnq3_rank_mean <- mean(lopo_pca_summary$KCNQ3_loading_rank, na.rm = TRUE)
kcnq3_rank_sd <- sd(lopo_pca_summary$KCNQ3_loading_rank, na.rm = TRUE)

rank_order <- lopo_pca_summary %>%
        arrange(KCNQ3_loading_rank) %>%
        pull(patient_removed)

kcnq3_rank_plot_df <- lopo_pca_summary %>%
        select(patient_removed, KCNQ3_loading_rank) %>%
        bind_rows(data.frame(patient_removed = "Full dataset",
                             KCNQ3_loading_rank = full_kcnq3_rank)) %>%
        mutate(patient_removed = factor(patient_removed,
                                        levels = c("Full dataset", rank_order)))

lopo_kcnq3_rank_plot <- ggplot(kcnq3_rank_plot_df,
                               aes(x = KCNQ3_loading_rank, y = patient_removed)) +
        annotate("rect",
                 xmin = kcnq3_rank_mean - kcnq3_rank_sd,
                 xmax = kcnq3_rank_mean + kcnq3_rank_sd,
                 ymin = -Inf, ymax = Inf,
                 fill = "red", alpha = 0.12) +
        geom_vline(xintercept = kcnq3_rank_mean,
                   color = "red", linetype = "dotted") +
        geom_point(aes(color = patient_removed == "Full dataset"), size = 2.5) +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                           guide = "none") +
        scale_x_continuous(limits = c(1, length(full_pc1_loadings))) +
        theme_classic() +
        labs(x = "KCNQ3 absolute PC1 loading rank",
             y = "Patient removed",
             title = "Leave-one-patient-out KCNQ3 loading rank",
             subtitle = "Rank 1 = strongest absolute loading")

lopo_kcnq3_rank_plot

# Print summary
cat("\n===== LOPO PC1 SUMMARY =====\n")
cat("Full-data PC1 variance:", round(full_pc1_variance, 2), "%\n")
cat("LOPO PC1 variance range:",
    round(min(lopo_pca_summary$PC1_variance), 2), "-",
    round(max(lopo_pca_summary$PC1_variance), 2), "%\n")
cat("PC1 loading correlation range:",
    round(min(lopo_pca_summary$PC1_loading_correlation), 3), "-",
    round(max(lopo_pca_summary$PC1_loading_correlation), 3), "\n")
cat("Full-data KCNQ3 loading:", round(full_kcnq3_loading, 4), "\n")
cat("LOPO KCNQ3 loading range:",
    round(min(lopo_pca_summary$KCNQ3_loading), 4), "-",
    round(max(lopo_pca_summary$KCNQ3_loading), 4), "\n")
cat("Full-data KCNQ3 absolute loading rank:", full_kcnq3_rank, "\n")
cat("LOPO KCNQ3 absolute loading rank range:",
    min(lopo_pca_summary$KCNQ3_loading_rank), "-",
    max(lopo_pca_summary$KCNQ3_loading_rank), "\n")
cat("Full-data overall ILAE p:", signif(full_partial_F$p.value, 3), "\n")
cat("LOPO overall ILAE p range:",
    signif(min(lopo_pca_summary$ILAE_partial_F_p), 3), "-",
    signif(max(lopo_pca_summary$ILAE_partial_F_p), 3), "\n")

### 4) Save LOPO outputs -------------------------------------------------------
write.csv(lopo_pca_summary, file.path(cfg$out_results_dir, "LOPO_PCA_summary.csv"), row.names = FALSE)
write.csv(lopo_ilae_results, file.path(cfg$out_results_dir, "LOPO_ILAE_coefficients.csv"), row.names = FALSE)
