################################################################################
#
# 02_patient_PC1_regression.R
#
# Model patient PC1 against ILAE score, adjusting for age and median transcript counts.
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
            out_results_dir = NA_character_,
            ilae_order = c("P", "0", "1", "2", "3"))

if (is.na(cfg$out_results_dir)) stop("Set cfg$out_results_dir")
dir.create(cfg$out_results_dir, recursive = TRUE, showWarnings = FALSE)

pca_checkpoint <- load_if_missing("pca_checkpoint", cfg$pca_checkpoint_rds, class = "list")
required <- c("pca_df", "pca", "patient_profiles", "variance_explained")
if (!all(required %in% names(pca_checkpoint))) stop("PCA checkpoint missing required fields")
pca_df <- pca_checkpoint$pca_df
pca <- pca_checkpoint$pca
patient_profiles <- pca_checkpoint$patient_profiles
variance_explained <- pca_checkpoint$variance_explained

### 3) Palettes -----------------------------------------------------------------

pal_ILAE_score <- c(
        "P" = "#C67BFF",
        "0" = "#F8766C",
        "1" = "#7BAD00",
        "2" = "#00C0C5",
        "3" = "#E69F00"
)

### 4) Regression: PC1 ~ ILAE + age + transcript count ------------------------

pc1_regression_df <- pca_df %>%
        filter(!is.na(PC1), !is.na(ILAE_score), !is.na(patient_age), !is.na(median_nCount_Xenium)) %>%
        mutate(ILAE_score = factor(ILAE_score, levels = cfg$ilae_order),
               patient_age_centered = patient_age - mean(patient_age),
               patient_count_centered = median_nCount_Xenium - mean(median_nCount_Xenium))

pc1_regression_model <- lm(PC1 ~ ILAE_score + patient_age_centered + patient_count_centered,
                           data = pc1_regression_df)

summary(pc1_regression_model)
drop1(pc1_regression_model, test = "F")

# Estimate adjusted ILAE group means and compare each epilepsy group with postmortem controls
pc1_adjusted_means <- emmeans(pc1_regression_model, ~ ILAE_score)
pc1_vs_P <- contrast(pc1_adjusted_means, method = "trt.vs.ctrl", ref = "P", adjust = "mvt")

pc1_adjusted_means
pc1_vs_P

pc1_adjusted_means_df <- as.data.frame(pc1_adjusted_means)

# Calculate patient PC1 scores adjusted to mean age and mean transcript count
adjustment_data <- pc1_regression_df %>%
        mutate(patient_age_centered = 0,
               patient_count_centered = 0)

pc1_regression_df$PC1_adjusted <- predict(pc1_regression_model, newdata = adjustment_data) +
        residuals(pc1_regression_model)

pc1_adjusted_plot <- ggplot(pc1_regression_df, aes(ILAE_score, PC1_adjusted)) +
        geom_jitter(aes(color = ILAE_score), width = 0.2, size = 2.5) +
        geom_point(data = pc1_adjusted_means_df, aes(ILAE_score, emmean),
                   size = 3, shape = 21, fill = "white", color = "black", inherit.aes = FALSE) +
        geom_errorbar(data = pc1_adjusted_means_df, aes(ILAE_score, ymin = lower.CL, ymax = upper.CL),
                      width = 0.12, color = "black", inherit.aes = FALSE) +
        scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
        scale_x_discrete(limits = c("P", "0", "2", "3", "1")) +
        theme_classic() +
        labs(x = "ILAE score", y = "Adjusted PC1 score") +
        theme(legend.position = "none")

pc1_adjusted_plot


### 5) Save regression checkpoint and tables -----------------------------------
regression_checkpoint <- list(pc1_regression_model = pc1_regression_model,
                              pc1_regression_df = pc1_regression_df,
                              pc1_adjusted_means_df = pc1_adjusted_means_df)
saveRDS(regression_checkpoint, file.path(cfg$out_results_dir, "patient_level_regression_checkpoint.rds"))
write.csv(pc1_regression_df, file.path(cfg$out_results_dir, "patient_level_adjusted_PC1.csv"), row.names = FALSE)
write.csv(as.data.frame(pc1_vs_P), file.path(cfg$out_results_dir, "patient_level_ILAE_vs_P.csv"), row.names = FALSE)
