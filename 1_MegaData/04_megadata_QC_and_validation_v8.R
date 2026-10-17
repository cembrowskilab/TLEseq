################################################################################
# Author: Larissa Kraus (October 2026)
# 05_megadata_QC_and_validation_v8.R
#
# Purpose:
#   Quality-control, covariate and technical-validation analyses for MegaData.
#
# Analyses:
#   1. PMI effects on quality and cell-type proportions
#   2. Age effects on quality and cell-type proportions
#   3. Age differences across ILAE groups
#   4. Power analysis
#   5. Split-half expression correlations
#   6. Across-sample Xenium-run correlations
#   7. Combined correlation figure
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
        library(pwr)
        library(FSA)
})

### PACKAGE CONFLICTS ---------------------------------------------------------
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)


### 1) User-editable config -----------------------------------------------------

cfg <- list(
        megadata_rds = NA_character_,
        celltype_order = c("ExN", "InN", "Astro", "Micro", "Oligo", "OPC", "Vascular", "CSF"),
        ilae_order = c("P", "0", "1", "2", "3"),
        split_half_seed = 123
)

pal_ILAE_score <- c(
        "P" = "#C67BFF",
        "0" = "#F8766C",
        "1" = "#7BAD00",
        "2" = "#00C0C5",
        "3" = "#E69F00"
)

pal_celltype <- c(
        "ExN"      = "#E69F00",
        "InN"      = "#56B4E9",
        "Astro"    = "#009E73",
        "Micro"    = "#F0E442",
        "Oligo"    = "#0072B2",
        "OPC"      = "#D55E00",
        "Vascular" = "#CC79A7",
        "CSF"      = "#666666"
)


### 2) Helper functions ---------------------------------------------------------

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


# Split cells within group x cell-type combinations and correlate mean expression
# Plot mean pairwise correlations using the shared cell-type order and scale.
make_final_cor_plot <- function(df, title, x_lab, y_lab, show_legend = FALSE) {
        df %>%
                mutate(celltype_A = factor(as.character(celltype_A),
                                           levels = rev(cfg$celltype_order)),
                       celltype_B = factor(as.character(celltype_B),
                                           levels = cfg$celltype_order)) %>%
                ggplot(aes(x = celltype_B, y = celltype_A,
                           fill = mean_correlation)) +
                geom_tile() +
                geom_text(aes(label = sprintf("%.3f", mean_correlation)), size = 2) +
                scale_fill_viridis_c(name = "Mean Pearson r",
                                     limits = c(min_cor_final, 1)) +
                coord_equal() +
                theme_classic() +
                labs(x = x_lab, y = y_lab, title = title) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      axis.text.y = element_text(angle = 45, hjust = 1),
                      axis.text = element_text(size = 6),
                      axis.title = element_text(size = 10),
                      legend.position = if (show_legend) "right" else "none")
}

split_half_heatmap <- function(dat, group_col, celltype_col = "celltype",
                               assay_use = "SCT", layer_use = "data",
                               group_order = NULL, celltype_order = NULL,
                               seed = 123, plot_title = NULL, facet_ncol = 6) {
        
        expr <- GetAssayData(dat, assay = assay_use, layer = layer_use)
        meta <- dat[[]][colnames(expr), , drop = FALSE]
        
        keep_cells <- !is.na(meta[[group_col]]) & !is.na(meta[[celltype_col]])
        meta <- meta[keep_cells, , drop = FALSE]
        expr <- expr[, rownames(meta), drop = FALSE]
        
        meta$split_group <- paste(meta[[group_col]], meta[[celltype_col]], sep = "___")
        
        set.seed(seed)
        meta$half <- ave(seq_len(nrow(meta)), meta$split_group, FUN = function(x) {
                sample(rep(c("A", "B"), length.out = length(x)))
        })
        
        split_ids <- unique(paste(meta$split_group, meta$half, sep = "___"))
        
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
                
                data.frame(group = group_i,
                           celltype = celltype_i,
                           half = half_i,
                           gene = rownames(expr),
                           expression = Matrix::rowMeans(expr[, cells, drop = FALSE]))
        }) %>%
                bind_rows()
        
        within_cor <- split_expression %>%
                select(group, celltype, half, gene, expression) %>%
                pivot_wider(names_from = half, values_from = expression) %>%
                group_by(group, celltype) %>%
                summarise(correlation = cor(A, B, method = "pearson"), .groups = "drop")
        
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
                        
                        common_celltypes <- intersect(colnames(mat_A), colnames(mat_B))
                        mat_A <- mat_A[, common_celltypes, drop = FALSE]
                        mat_B <- mat_B[, common_celltypes, drop = FALSE]
                        
                        cor_mat <- cor(mat_A, mat_B, method = "pearson")
                        
                        as.data.frame(as.table(cor_mat)) %>%
                                rename(celltype1 = Var1,
                                       celltype2 = Var2,
                                       correlation = Freq) %>%
                                mutate(group = group_i) %>%
                                filter(celltype1 != celltype2)
                }) %>%
                bind_rows()
        
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
        
        if (is.null(celltype_order)) {
                celltype_order <- sort(unique(c(heatmap_df$celltype_A,
                                                heatmap_df$celltype_B)))
        }
        
        if (is.null(group_order)) group_order <- unique(heatmap_df$group)
        
        heatmap_df <- heatmap_df %>%
                mutate(group = factor(group, levels = group_order),
                       celltype_A = factor(celltype_A, levels = rev(celltype_order)),
                       celltype_B = factor(celltype_B, levels = celltype_order))
        
        p_group <- ggplot(heatmap_df,
                          aes(x = celltype_B, y = celltype_A, fill = correlation)) +
                geom_tile() +
                scale_fill_viridis_c(name = "Pearson r", limits = c(0, 1)) +
                facet_wrap(~group, ncol = facet_ncol) +
                coord_equal() +
                theme_classic() +
                labs(x = "Half B", y = "Half A", title = plot_title) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.grid = element_blank(),
                      strip.background = element_blank(),
                      strip.text = element_text(face = "bold"))
        
        heatmap_mean <- heatmap_df %>%
                group_by(celltype_A, celltype_B) %>%
                summarise(mean_correlation = mean(correlation, na.rm = TRUE), .groups = "drop")
        
        p_mean <- ggplot(heatmap_mean,
                         aes(x = celltype_B, y = celltype_A, fill = mean_correlation)) +
                geom_tile() +
                geom_text(aes(label = sprintf("%.3f", mean_correlation)), size = 2) +
                scale_fill_viridis_c(name = "Mean Pearson r", limits = c(0, 1)) +
                coord_equal() +
                theme_classic() +
                labs(x = "Half B", y = "Half A", title = paste("Mean:", plot_title)) +
                theme(axis.text = element_text(size = 6),
                      axis.title = element_text(size = 10))
        
        list(split_expression = split_expression,
             within_cor = within_cor,
             between_cor = between_cor,
             heatmap_df = heatmap_df,
             heatmap_mean = heatmap_mean,
             p_group = p_group,
             p_mean = p_mean)
}


### 3) Load MegaData ------------------------------------------------------------

dat <- load_if_missing("dat", cfg$megadata_rds, class = "Seurat")

dat$ILAE_score <- factor(dat$ILAE_score, levels = cfg$ilae_order)
dat$celltype <- factor(dat$celltype, levels = cfg$celltype_order)


### 4) PMI effects on quality and cell-type proportions -------------------------

# Extract metadata needed for PMI analyses
meta_pmi <- dat@meta.data %>%
        select(orig.ident, celltype, type, PMI, nCount_Xenium)

# Create one row of patient-level metadata per sample
patient_meta_pmi <- meta_pmi %>%
        select(orig.ident, type, PMI) %>%
        distinct()

# Calculate cell-type proportions and retain zero-count combinations
cell_prop_pmi <- meta_pmi %>%
        count(orig.ident, celltype, name = "n_cells") %>%
        complete(orig.ident, celltype, fill = list(n_cells = 0)) %>%
        group_by(orig.ident) %>%
        mutate(n_total = sum(n_cells),
               proportion = n_cells / n_total) %>%
        ungroup() %>%
        left_join(patient_meta_pmi, by = "orig.ident")

# Test association between cell-type proportion and PMI in postmortem samples
pmi_results <- lapply(cfg$celltype_order, function(ct) {
        tmp <- cell_prop_pmi %>%
                filter(celltype == ct, type == "postmortem", !is.na(PMI))
        
        if (nrow(tmp) < 3) return(NULL)
        
        test <- cor.test(tmp$proportion, tmp$PMI, method = "pearson")
        
        data.frame(celltype = ct,
                   n = nrow(tmp),
                   r = unname(test$estimate),
                   p_value = test$p.value)
}) %>%
        bind_rows() %>%
        mutate(p_adjusted = p.adjust(p_value, method = "BH"),
               label = paste0("r = ", round(r, 2),
                              "\nPadj = ", signif(p_adjusted, 2)),
               celltype = factor(celltype, levels = cfg$celltype_order))

cell_prop_pmi$celltype <- factor(cell_prop_pmi$celltype, levels = cfg$celltype_order)

# Plot cell-type proportion against PMI
pmi_prop_plot <- ggplot(cell_prop_pmi %>%
                                filter(type == "postmortem", !is.na(PMI)),
                        aes(x = PMI, y = proportion)) +
        geom_point(size = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
        facet_wrap(~celltype, scales = "free_y", nrow = 1) +
        geom_text(data = pmi_results,
                  aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.2, vjust = 1, size = 1,
                  inherit.aes = FALSE) +
        scale_y_continuous(labels = scales::label_percent(accuracy = 1, suffix = "", big.mark = ""),
                           limits = c(-0.01, NA),
                           expand = expansion(mult = c(0, 0.05))) +
        labs(x = "PMI (days)", y = "Cell proportion (%)") +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 10),
              panel.spacing.x = unit(0.1, "lines"))

pmi_prop_plot

# Plot cell-type composition across postmortem samples
postmortem_prop_plot <- ggplot(cell_prop_pmi %>%
                                       filter(type == "postmortem"),
                               aes(x = orig.ident, y = proportion, fill = celltype)) +
        geom_col(width = 0.8) +
        scale_fill_manual(values = pal_celltype) +
        scale_y_continuous(labels = scales::percent) +
        labs(x = "Patient", y = "Cell proportion", fill = "Cell type") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

postmortem_prop_plot

# Calculate patient-level median transcripts/cell and number of cells
pmi_quality <- meta_pmi %>%
        filter(type == "postmortem", !is.na(PMI)) %>%
        group_by(orig.ident, PMI) %>%
        summarise(median_nCount_Xenium = median(nCount_Xenium, na.rm = TRUE),
                  n_cells = n(),
                  .groups = "drop")

# Test associations with PMI
test_pmi_nCount <- cor.test(pmi_quality$median_nCount_Xenium, pmi_quality$PMI,
                            method = "pearson")

test_pmi_ncells <- cor.test(pmi_quality$n_cells, pmi_quality$PMI,
                            method = "pearson")

label_pmi_nCount <- paste0("r = ", round(unname(test_pmi_nCount$estimate), 2),
                           "\nP = ", signif(test_pmi_nCount$p.value, 2))

label_pmi_ncells <- paste0("r = ", round(unname(test_pmi_ncells$estimate), 2),
                           "\nP = ", signif(test_pmi_ncells$p.value, 2))

# Plot median transcripts/cell against PMI
pmi_nCount_plot <- ggplot(pmi_quality,
                          aes(x = PMI, y = median_nCount_Xenium)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        annotate("text", x = -Inf, y = Inf, label = label_pmi_nCount,
                 hjust = -0.2, vjust = 1.5) +
        labs(x = "PMI (days)", y = "Median transcripts/cell") +
        theme_classic()

pmi_nCount_plot

# Plot number of cells against PMI
pmi_ncells_plot <- ggplot(pmi_quality,
                          aes(x = PMI, y = n_cells)) +
        geom_point(size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        annotate("text", x = -Inf, y = Inf, label = label_pmi_ncells,
                 hjust = -0.2, vjust = 1.5) +
        scale_y_continuous(labels = function(x) x / 1000) +
        labs(x = "PMI (days)",
             y = expression("Number of cells ("*10^3*")")) +
        theme_classic()

pmi_ncells_plot


### 5) Age effects on quality and cell-type proportions -------------------------

# Extract metadata needed for age analyses
meta_age <- dat@meta.data %>%
        select(orig.ident, celltype, patient_age, type, ILAE_score, nCount_Xenium)

# Create one row of patient-level metadata per sample
patient_meta_age <- meta_age %>%
        select(orig.ident, patient_age, type, ILAE_score) %>%
        distinct()

# Calculate cell-type proportions and retain zero-count combinations
cell_prop_age <- meta_age %>%
        count(orig.ident, celltype, name = "n_cells") %>%
        complete(orig.ident, celltype, fill = list(n_cells = 0)) %>%
        group_by(orig.ident) %>%
        mutate(n_total = sum(n_cells),
               proportion = n_cells / n_total) %>%
        ungroup() %>%
        left_join(patient_meta_age, by = "orig.ident")

# Test association between cell-type proportion and age
age_results <- lapply(cfg$celltype_order, function(ct) {
        tmp <- cell_prop_age %>%
                filter(celltype == ct, !is.na(patient_age))
        
        test <- cor.test(tmp$proportion, tmp$patient_age, method = "pearson")
        
        data.frame(celltype = ct,
                   n = nrow(tmp),
                   r = unname(test$estimate),
                   p_value = test$p.value)
}) %>%
        bind_rows() %>%
        mutate(p_adjusted = p.adjust(p_value, method = "BH"),
               label = paste0("r = ", round(r, 2),
                              "\nPadj = ", signif(p_adjusted, 2)),
               celltype = factor(celltype, levels = cfg$celltype_order))

cell_prop_age$celltype <- factor(cell_prop_age$celltype, levels = cfg$celltype_order)

# Plot cell-type proportion against age
age_prop_plot <- ggplot(cell_prop_age %>%
                                filter(!is.na(patient_age)),
                        aes(x = patient_age, y = proportion)) +
        geom_point(aes(color = ILAE_score), size = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
        scale_color_manual(values = pal_ILAE_score, name = "ILAE score") +
        facet_wrap(~celltype, scales = "free_y", nrow = 1) +
        geom_text(data = age_results,
                  aes(x = -Inf, y = Inf, label = label),
                  hjust = -0.2, vjust = 1, size = 1,
                  inherit.aes = FALSE) +
        scale_y_continuous(labels = scales::label_percent(accuracy = 1, suffix = "", big.mark = ""),
                           limits = c(-0.01, NA),
                           expand = expansion(mult = c(0, 0.05))) +
        labs(x = "Age (years)", y = "Cell proportion (%)") +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 10),
              panel.spacing.x = unit(0.1, "lines"))

age_prop_plot

# Calculate patient-level median transcripts/cell and number of cells
age_quality <- meta_age %>%
        filter(!is.na(patient_age)) %>%
        group_by(orig.ident, patient_age, ILAE_score) %>%
        summarise(median_nCount_Xenium = median(nCount_Xenium, na.rm = TRUE),
                  n_cells = n(),
                  .groups = "drop")

age_quality$ILAE_score <- factor(age_quality$ILAE_score, levels = cfg$ilae_order)

# Test associations with age
test_age_nCount <- cor.test(age_quality$median_nCount_Xenium, age_quality$patient_age,
                            method = "spearman", exact = FALSE)

test_age_ncells <- cor.test(age_quality$n_cells, age_quality$patient_age,
                            method = "spearman", exact = FALSE)

label_age_nCount <- paste0("r = ", round(unname(test_age_nCount$estimate), 2),
                           "\nP = ", signif(test_age_nCount$p.value, 2))

label_age_ncells <- paste0("r = ", round(unname(test_age_ncells$estimate), 2),
                           "\nP = ", signif(test_age_ncells$p.value, 2))

# Plot median transcripts/cell against age
age_nCount_plot <- ggplot(age_quality,
                          aes(x = patient_age, y = median_nCount_Xenium)) +
        geom_point(aes(color = ILAE_score), size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        scale_color_manual(values = pal_ILAE_score, name = "ILAE score") +
        annotate("text", x = -Inf, y = Inf, label = label_age_nCount,
                 hjust = -0.2, vjust = 1.5) +
        labs(x = "Age (years)", y = "Median transcripts/cell") +
        theme_classic()

age_nCount_plot

# Plot number of cells against age
age_ncells_plot <- ggplot(age_quality,
                          aes(x = patient_age, y = n_cells)) +
        geom_point(aes(color = ILAE_score), size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        scale_color_manual(values = pal_ILAE_score, name = "ILAE score") +
        annotate("text", x = -Inf, y = Inf, label = label_age_ncells,
                 hjust = -0.2, vjust = 1.5) +
        scale_y_continuous(labels = function(x) x / 1000) +
        labs(x = "Age (years)",
             y = expression("Number of cells ("*10^3*")")) +
        theme_classic()

age_ncells_plot


### 6) Age across ILAE groups ---------------------------------------------------

# Create one row per patient
age_groups <- meta_age %>%
        select(orig.ident, patient_age, ILAE_score) %>%
        distinct() %>%
        filter(!is.na(patient_age), !is.na(ILAE_score)) %>%
        mutate(ILAE_score = factor(ILAE_score, levels = cfg$ilae_order))

# Test overall difference in age between ILAE groups
age_kw <- kruskal.test(patient_age ~ ILAE_score, data = age_groups)
age_kw

# Dunn post-hoc comparisons with BH correction
age_dunn <- FSA::dunnTest(patient_age ~ ILAE_score,
                          data = age_groups,
                          method = "bh")$res
age_dunn

# Plot age by ILAE group
age_group_plot <- ggplot(age_groups,
                         aes(x = ILAE_score, y = patient_age, fill = ILAE_score)) +
        geom_boxplot(width = 0.6, outlier.shape = NA) +
        geom_jitter(width = 0.12, size = 0.7) +
        scale_fill_manual(values = pal_ILAE_score) +
        labs(x = "ILAE group", y = "Age (years)") +
        theme_classic() +
        theme(legend.position = "none",
              axis.text = element_text(size = 6),
              axis.title = element_text(size = 10))

age_group_plot


### 7) Power analysis -----------------------------------------------------------

# Determine number of patients in each ILAE group
patient_groups <- dat@meta.data %>%
        transmute(patient = as.character(orig.ident),
                  ILAE_score = as.character(ILAE_score)) %>%
        distinct(patient, .keep_all = TRUE)

group_sizes <- patient_groups %>%
        count(ILAE_score, name = "n")

group_sizes

n_group <- setNames(group_sizes$n, group_sizes$ILAE_score)

stopifnot(all(cfg$ilae_order %in% names(n_group)))
stopifnot(n_group["P"] == n_group["0"])

# P and ILAE 0 are represented together because they have the same sample size
comparisons <- data.frame(
        group1 = c("P/0", "P/0", "P/0", "1", "1", "2"),
        group2 = c("1", "2", "3", "2", "3", "3"),
        n1 = c(n_group["P"], n_group["P"], n_group["P"],
               n_group["1"], n_group["1"], n_group["2"]),
        n2 = c(n_group["1"], n_group["2"], n_group["3"],
               n_group["2"], n_group["3"], n_group["3"])
) %>%
        mutate(comparison_label = paste0(group1, " (n=", n1, ") - ",
                                         group2, " (n=", n2, ")"))

effect_sizes <- seq(0.01, 3, by = 0.01)

power_results <- comparisons %>%
        crossing(d = effect_sizes) %>%
        rowwise() %>%
        mutate(power = pwr.t2n.test(n1 = n1, n2 = n2, d = d,
                                    sig.level = 0.05,
                                    alternative = "two.sided")$power) %>%
        ungroup() %>%
        mutate(comparison_label = factor(comparison_label,
                                         levels = comparisons$comparison_label))

power_80 <- comparisons %>%
        rowwise() %>%
        mutate(d_80 = pwr.t2n.test(n1 = n1, n2 = n2,
                                   power = 0.80,
                                   sig.level = 0.05,
                                   alternative = "two.sided")$d) %>%
        ungroup()

power_80

power_colors <- c("#0072B2", "#E69F00", "#009E73",
                  "#CC79A7", "#56B4E9", "#D55E00")
names(power_colors) <- comparisons$comparison_label

p_power <- ggplot(power_results,
                  aes(x = d, y = power,
                      group = comparison_label,
                      color = comparison_label)) +
        geom_line(linewidth = 0.8) +
        geom_hline(yintercept = 0.80, linetype = "dashed") +
        scale_color_manual(values = power_colors) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        theme_classic() +
        labs(x = "Cohen's d",
             y = "Statistical power",
             color = "Comparison")

p_power


### 8) Split-half correlation analysis -----------------------------------------

# Split cells within each sample and correlate mean expression profiles
sample_cor <- split_half_heatmap(dat = dat,
                                 group_col = "orig.ident",
                                 celltype_col = "celltype",
                                 assay_use = "SCT",
                                 layer_use = "data",
                                 celltype_order = cfg$celltype_order,
                                 seed = cfg$split_half_seed,
                                 plot_title = "Grouped by sample")

sample_cor$p_group
sample_cor$p_mean

# Split cells within each ILAE group and correlate mean expression profiles
ILAE_cor <- split_half_heatmap(dat = dat,
                               group_col = "ILAE_score",
                               celltype_col = "celltype",
                               assay_use = "SCT",
                               layer_use = "data",
                               group_order = cfg$ilae_order,
                               celltype_order = cfg$celltype_order,
                               seed = cfg$split_half_seed,
                               plot_title = "Grouped by ILAE score",
                               facet_ncol = 5)

ILAE_cor$p_group
ILAE_cor$p_mean


### 9) Across-sample Xenium-run correlation analysis ----------------------------

# Add Xenium run from metadata if not already present
if (!"Xenium_run" %in% colnames(dat@meta.data)) {
        xenium_meta <- read.csv(cfg$metadata_csv)
        
        if (!all(c("Sample", "Xenium_run") %in% colnames(xenium_meta))) {
                stop("metadata_csv must contain Sample and Xenium_run columns.")
        }
        
        dat$Xenium_run <- xenium_meta$Xenium_run[
                match(dat$orig.ident, xenium_meta$Sample)]
}

if (anyNA(dat$Xenium_run)) warning("Some cells have no Xenium_run assignment.")

expr_batch <- GetAssayData(dat, assay = "SCT", layer = "data")
meta_batch <- dat[[]][colnames(expr_batch), , drop = FALSE]

# Identify observed sample x cell-type combinations
sample_celltype_groups <- meta_batch %>%
        filter(!is.na(orig.ident), !is.na(celltype)) %>%
        distinct(orig.ident, celltype)

# Calculate mean expression for each sample x cell-type combination
sample_expression <- lapply(seq_len(nrow(sample_celltype_groups)), function(i) {
        sample_i <- as.character(sample_celltype_groups$orig.ident[i])
        celltype_i <- as.character(sample_celltype_groups$celltype[i])
        
        cells <- rownames(meta_batch)[
                as.character(meta_batch$orig.ident) == sample_i &
                        as.character(meta_batch$celltype) == celltype_i
        ]
        
        data.frame(sample = sample_i,
                   celltype = celltype_i,
                   gene = rownames(expr_batch),
                   expression = Matrix::rowMeans(expr_batch[, cells, drop = FALSE]))
}) %>%
        bind_rows()

# Assign each sample to its Xenium run
sample_batch <- meta_batch %>%
        transmute(sample = as.character(orig.ident),
                  batch = as.character(Xenium_run)) %>%
        filter(!is.na(sample), !is.na(batch)) %>%
        distinct()

# Define all sample pairs as within or across Xenium runs
sample_pairs <- combn(sample_batch$sample, 2, simplify = FALSE) %>%
        lapply(function(x) data.frame(sample_A = x[1], sample_B = x[2])) %>%
        bind_rows() %>%
        left_join(sample_batch %>%
                          rename(sample_A = sample, batch_A = batch),
                  by = "sample_A") %>%
        left_join(sample_batch %>%
                          rename(sample_B = sample, batch_B = batch),
                  by = "sample_B") %>%
        mutate(batch_comparison = if_else(batch_A == batch_B,
                                          "Within batch", "Across batch"))

table(sample_pairs$batch_comparison)

# Correlate cell-type expression profiles between each pair of samples
across_sample_cor <- lapply(seq_len(nrow(sample_pairs)), function(i) {
        sample_A_i <- sample_pairs$sample_A[i]
        sample_B_i <- sample_pairs$sample_B[i]
        
        mat_A <- sample_expression %>%
                filter(sample == sample_A_i) %>%
                select(gene, celltype, expression) %>%
                pivot_wider(names_from = celltype, values_from = expression) %>%
                tibble::column_to_rownames("gene") %>%
                as.matrix()
        
        mat_B <- sample_expression %>%
                filter(sample == sample_B_i) %>%
                select(gene, celltype, expression) %>%
                pivot_wider(names_from = celltype, values_from = expression) %>%
                tibble::column_to_rownames("gene") %>%
                as.matrix()
        
        common_celltypes <- intersect(colnames(mat_A), colnames(mat_B))
        
        cor_mat <- cor(mat_A[, common_celltypes, drop = FALSE],
                       mat_B[, common_celltypes, drop = FALSE],
                       method = "pearson")
        
        as.data.frame(as.table(cor_mat)) %>%
                rename(celltype_A = Var1,
                       celltype_B = Var2,
                       correlation = Freq) %>%
                mutate(sample_A = sample_A_i,
                       sample_B = sample_B_i)
}) %>%
        bind_rows() %>%
        left_join(sample_pairs %>%
                          select(sample_A, sample_B, batch_A, batch_B,
                                 batch_comparison),
                  by = c("sample_A", "sample_B"))

# Calculate mean correlations across Xenium runs
cor_across <- across_sample_cor %>%
        filter(batch_comparison == "Across batch") %>%
        group_by(celltype_A, celltype_B) %>%
        summarise(mean_correlation = mean(correlation, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(comparison = "Across batch")

# Calculate mean correlations within Xenium runs
cor_within <- across_sample_cor %>%
        filter(batch_comparison == "Within batch") %>%
        group_by(celltype_A, celltype_B) %>%
        summarise(mean_correlation = mean(correlation, na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(comparison = "Within batch")

# Plot within- and across-run correlations
batch_cor_heatmap <- bind_rows(cor_across, cor_within) %>%
        mutate(comparison = factor(comparison,
                                   levels = c("Across batch", "Within batch")),
               celltype_A = factor(celltype_A,
                                   levels = rev(cfg$celltype_order)),
               celltype_B = factor(celltype_B,
                                   levels = cfg$celltype_order))

min_cor_batch <- min(batch_cor_heatmap$mean_correlation, na.rm = TRUE)

p_batch_comparison <- ggplot(batch_cor_heatmap,
                             aes(x = celltype_B, y = celltype_A,
                                 fill = mean_correlation)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.3f", mean_correlation)), size = 2) +
        scale_fill_viridis_c(name = "Mean Pearson r",
                             limits = c(min_cor_batch, 1)) +
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


### 10) Combined correlation figure --------------------------------------------

# Use one color scale and cell-type order across correlation heatmaps
min_cor_final <- min(cor_across$mean_correlation,
                     cor_within$mean_correlation,
                     sample_cor$heatmap_mean$mean_correlation,
                     ILAE_cor$heatmap_mean$mean_correlation,
                     na.rm = TRUE)


p_sample_final <- make_final_cor_plot(sample_cor$heatmap_mean,
                                      "Within sample",
                                      "Half B", "Half A")

p_ILAE_final <- make_final_cor_plot(ILAE_cor$heatmap_mean,
                                    "Within ILAE group",
                                    "Half B", "Half A")

p_within_final <- make_final_cor_plot(cor_within,
                                      "Within batch",
                                      "Cell type (sample 2)",
                                      "Cell type (sample 1)")

p_across_final <- make_final_cor_plot(cor_across,
                                      "Across batch",
                                      "Cell type (sample 2)",
                                      "Cell type (sample 1)",
                                      show_legend = TRUE)

p_extfig2 <- wrap_plots(p_sample_final,
                        p_within_final,
                        p_across_final,
                        nrow = 1)

p_extfig2