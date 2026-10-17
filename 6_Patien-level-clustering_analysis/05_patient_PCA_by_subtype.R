################################################################################
#
# 05_patient_PCA_by_subtype.R
#
# ExN, InN and Micro patient PCA using equally weighted subtype mean profiles.
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

### 2) User-editable paths and load inputs -------------------------------------
cfg <- list(megadata_rds = NA_character_, 
            out_results_dir = NA_character_,
            subtype_column = "cellsubtype", 
            assay = "Xenium",
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
dat_joined <- JoinLayers(dat, assay = cfg$assay)
dat_joined$ILAE_score <- factor(dat_joined$ILAE_score, levels = cfg$ilae_order)

### 5) Build subtype-weighted patient profiles and run PCA ---------------------
subtype_column <- cfg$subtype_column
assay_to_use <- cfg$assay

# Only use cells assigned to an ExN, InN or Micro subtype
dat_joined$major_celltype <- case_when(
        grepl("^ExN_", dat_joined[[subtype_column, drop = TRUE]]) ~ "ExN",
        grepl("^InN_", dat_joined[[subtype_column, drop = TRUE]]) ~ "InN",
        grepl("^Micro_", dat_joined[[subtype_column, drop = TRUE]]) ~ "Micro",
        TRUE ~ NA_character_
)

cells_use <- colnames(dat_joined)[!is.na(dat_joined$major_celltype)]
dat_joined <- subset(dat_joined, cells = cells_use)

# Run PCA for each major cell type using equally weighted subtype means
celltypes_pca <- c("ExN", "InN", "Micro")
pca_scores_subtypes <- data.frame()
pca_subtype_results <- list()
patient_column <- "orig.ident"

for (current_celltype in celltypes_pca) {
        
        cat("\n\n############################################################\n")
        cat("RUNNING:", current_celltype, "\n")
        cat("############################################################\n")
        # SUBSET CELL TYPE
        cells_current <- colnames(dat_joined)[dat_joined$major_celltype == current_celltype]
        dat_ct <- subset(dat_joined, cells = cells_current)
        
        cat("Cells:", ncol(dat_ct), "\n")
        cat("Patients:", length(unique(dat_ct[[patient_column, drop = TRUE]])), "\n")
        cat("Subtypes:", paste(sort(unique(dat_ct[[subtype_column, drop = TRUE]])), collapse = ", "), "\n")
        # PATIENT METADATA
        patient_metadata <- dat_ct@meta.data %>%
                mutate(patient = trimws(as.character(.data[[patient_column]]))) %>%
                select(patient, ILAE_score, sex, patient_age) %>%
                distinct(patient, .keep_all = TRUE)
        # RAW COUNTS AND PATIENT × SUBTYPE GROUPS
        counts <- LayerData(dat_ct, assay = assay_to_use, layer = "counts")
        
        meta <- dat_ct@meta.data %>%
                mutate(patient = trimws(as.character(.data[[patient_column]])),
                       subtype = trimws(as.character(.data[[subtype_column]])),
                       patient_subtype = paste(patient, subtype, sep = "___"))
        # MEAN RAW COUNTS PER PATIENT × SUBTYPE
        patient_subtype_groups <- split(seq_len(ncol(counts)), meta$patient_subtype)
        
        subtype_mean_matrix <- sapply(patient_subtype_groups, function(idx) {
                Matrix::rowMeans(counts[, idx, drop = FALSE])
        })
        # IDENTIFY PATIENT FOR EACH SUBTYPE PROFILE
        subtype_info <- data.frame(
                patient_subtype = colnames(subtype_mean_matrix),
                stringsAsFactors = FALSE
        ) %>%
                separate(
                        patient_subtype,
                        into = c("patient", "subtype"),
                        sep = "___",
                        extra = "merge",
                        remove = FALSE
                )
        # AGGREGATE SUBTYPE MEANS PER PATIENT
        patients <- unique(subtype_info$patient)
        
        patient_profiles_raw <- sapply(patients, function(current_patient) {
                cols <- which(subtype_info$patient == current_patient)
                Matrix::rowSums(subtype_mean_matrix[, cols, drop = FALSE])
        })
        
        patient_profiles_raw <- as.matrix(patient_profiles_raw)
        colnames(patient_profiles_raw) <- patients
        # CHECK PATIENT PROFILES
        cat("Patient-subtype profiles:", ncol(subtype_mean_matrix), "\n")
        cat("Final patient profiles:", ncol(patient_profiles_raw), "\n")
        cat("Expected patients:", length(unique(dat_ct[[patient_column, drop = TRUE]])), "\n")
        # LOG NORMALIZE PATIENT PROFILES
        normalized_matrix <- sweep(
                patient_profiles_raw,
                MARGIN = 2,
                STATS = Matrix::colSums(patient_profiles_raw),
                FUN = "/"
        ) * 10000
        
        patient_profiles_matrix <- log1p(normalized_matrix)
        # PCA
        gene_variance <- apply(patient_profiles_matrix, 1, var)
        genes_for_pca <- names(gene_variance[is.finite(gene_variance) & gene_variance > 0])
        
        if (length(genes_for_pca) < 2 || ncol(patient_profiles_raw) < 3) stop("Insufficient genes/patients for ", current_celltype)
        pca <- prcomp(t(patient_profiles_matrix[genes_for_pca, , drop = FALSE]),
                      center = TRUE, scale. = TRUE)
        
        variance_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
        
        pca_df <- data.frame(patient = rownames(pca$x),
                             PC1 = pca$x[, 1],
                             PC2 = pca$x[, 2]) %>%
                left_join(patient_metadata, by = "patient") %>%
                mutate(celltype = current_celltype)
        
        cat("PC1 variance:", round(variance_explained[1], 1), "%\n")
        cat("PC2 variance:", round(variance_explained[2], 1), "%\n")
        # PCA PLOT
        pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = ILAE_score)) +
                geom_point(size = 3) +
                scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
                coord_fixed() +
                theme_classic() +
                labs(title = paste0(current_celltype, " equal-subtype-weighted PCA"),
                     x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
                     y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
                     color = "ILAE score") +
                theme(legend.position = "bottom")
        
        print(pca_plot)
        # Save results
        pca_scores_subtypes <- bind_rows(pca_scores_subtypes, pca_df)
        
        pca_subtype_results[[current_celltype]] <- list(
                pca = pca,
                scores = pca_df,
                patient_profiles = patient_profiles_matrix,
                patient_profiles_raw = patient_profiles_raw,
                subtype_profiles = subtype_mean_matrix,
                subtype_info = subtype_info,
                variance_explained = variance_explained
        )
}


# Clean temporary objects
rm(dat_ct, cells_current, patient_metadata, counts, meta,
   patient_subtype_groups, subtype_mean_matrix, subtype_info,
   patients, patient_profiles_raw, normalized_matrix,
   patient_profiles_matrix, gene_variance, genes_for_pca,
   pca, variance_explained, pca_df, pca_plot)


### 6) Align PC1 directions across major cell types ----------------------------
reference_celltype <- "ExN"
reference_loadings <- pca_subtype_results[[reference_celltype]]$pca$rotation[, "PC1"]

for (current_celltype in setdiff(celltypes_pca, reference_celltype)) {
        
        current_loadings <- pca_subtype_results[[current_celltype]]$pca$rotation[, "PC1"]
        common_genes <- intersect(names(reference_loadings), names(current_loadings))
        
        loading_r <- cor(reference_loadings[common_genes],
                         current_loadings[common_genes],
                         method = "pearson")
        
        cat(current_celltype, "vs", reference_celltype,
            "PC1 loading correlation:", round(loading_r, 3), "\n")
        
        if (loading_r < 0) {
                pca_subtype_results[[current_celltype]]$pca$x[, "PC1"] <-
                        -pca_subtype_results[[current_celltype]]$pca$x[, "PC1"]
                
                pca_subtype_results[[current_celltype]]$pca$rotation[, "PC1"] <-
                        -pca_subtype_results[[current_celltype]]$pca$rotation[, "PC1"]
                
                pca_subtype_results[[current_celltype]]$scores$PC1 <-
                        -pca_subtype_results[[current_celltype]]$scores$PC1
                
                cat("Flipped", current_celltype, "PC1\n")
        }
}

pca_scores_subtypes <- bind_rows(
        lapply(pca_subtype_results, function(x) x$scores)
)

### 7) Compare PC1 scores across ExN, InN and Micro ----------------------------
pc1_scatterplots <- list()

for (i in 1:(length(celltypes_pca) - 1)) {
        for (j in (i + 1):length(celltypes_pca)) {
                
                celltype_A <- celltypes_pca[i]
                celltype_B <- celltypes_pca[j]
                
                plot_df <- pca_scores_subtypes %>%
                        filter(celltype %in% c(celltype_A, celltype_B)) %>%
                        select(patient, celltype, PC1, ILAE_score) %>%
                        pivot_wider(names_from = celltype, values_from = PC1) %>%
                        distinct(patient, .keep_all = TRUE) %>%
                        filter(!is.na(.data[[celltype_A]]),
                               !is.na(.data[[celltype_B]]))
                
                r <- cor(plot_df[[celltype_A]], plot_df[[celltype_B]],
                         method = "pearson")
                
                p <- cor.test(plot_df[[celltype_A]], plot_df[[celltype_B]],
                              method = "pearson")$p.value
                
                current_plot <- ggplot(plot_df,
                                       aes(x = .data[[celltype_A]],
                                           y = .data[[celltype_B]],
                                           color = ILAE_score)) +
                        geom_point(size = 1) +
                        geom_smooth(aes(group = 1), method = "lm",
                                    se = TRUE, color = "black") +
                        scale_color_manual(values = pal_ILAE_score, drop = FALSE) +
                        theme_classic() +
                        labs(x = paste0(celltype_A, " PC1"),
                             y = paste0(celltype_B, " PC1"),
                             title = paste0(celltype_A, " - ", celltype_B,
                                            " (r = ", sprintf("%.2f", r),
                                            ", p = ", format.pval(p, digits = 2), ")"),
                             color = "ILAE score") +
                        theme(
                                axis.text = element_text(size = 6),
                                axis.title = element_text(size = 10),
                                plot.title = element_text(size = 8),
                                legend.title = element_text(size = 8),
                                legend.text = element_text(size = 6),
                                legend.position = "bottom"
                        )
                
                pc1_scatterplots[[paste0(celltype_A, "_", celltype_B)]] <- current_plot
        }
}

subtype_combine_plot <- wrap_plots(pc1_scatterplots, ncol = 3)
subtype_combine_plot

# Save subtype-weighted PCA results for reuse in downstream analyses.
saveRDS(pca_subtype_results, file.path(cfg$out_results_dir, "subtype_PCA_results.rds"))
write.csv(pca_scores_subtypes, file.path(cfg$out_results_dir, "subtype_PCA_scores.csv"), row.names = FALSE)
