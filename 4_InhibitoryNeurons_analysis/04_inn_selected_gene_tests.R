# 04_inn_selected_gene_tests.R
#
# Selected-gene comparisons in inhibitory-neuron pseudobulk profiles | MegaData_v8
# Input: datIN_v8.rds from script 01 and regional metadata
# genes were selected from DEG analysis in 03_inn_region_deg.R
# LK, October 2026
#
# This script:
# 1. Aggregates SST, PV and VIP counts per patient and anatomical region
# 2. Tests GABRA1, GABRA2 and SNCG across ILAE groups (Kruskal-Wallis + Dunn)
# 3. Displays boxplots and inspects significant post hoc comparisons
#
### Setup ----------------------------------------------------------------------

suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
        library(ggplot2)
        library(patchwork)
})

options(future.globals.maxSize = 100000 * 1024^2)

### Function -------------------------------------------------------------------

# Test one gene across ILAE groups within each anatomical region.
# Each point represents a patient-level pseudobulk profile, not a single cell.
# The Kruskal-Wallis test is followed by BH-adjusted Dunn tests when significant.
Test_Plot_DEGs <- function(seurat_obj, gene, cluster_col = "broad_subtype_ExN",
                           ilae_col = "ILAE_score", assay = "Xenium",
                           posthoc_p_cutoff = 0.05, cluster_list = NULL, ncol = 3) {
        
        if (!cluster_col %in% colnames(seurat_obj@meta.data)) stop("Missing cluster column: ", cluster_col)
        if (!ilae_col %in% colnames(seurat_obj@meta.data)) stop("Missing ILAE column: ", ilae_col)
        if (!gene %in% rownames(seurat_obj[[assay]])) {
                warning("Gene not found in assay ", assay, ": ", gene)
                return(NULL)
        }
        if (is.null(cluster_list)) cluster_list <- unique(as.character(seurat_obj@meta.data[[cluster_col]]))
        cluster_list <- cluster_list[!is.na(cluster_list)]
        
        # Fetch once rather than for every subtype; explicitly use Xenium assay.
        expr <- FetchData(seurat_obj[[assay]], vars = gene, layer = "data")
        meta <- seurat_obj@meta.data
        meta$expression <- expr[rownames(meta), gene]
        meta$ILAE_group <- factor(sub("^g(?=[0-9]+$)", "", as.character(meta[[ilae_col]]), perl = TRUE),
                                  levels = c("P", "0", "1", "2", "3"))
        
        plot_list <- list()
        test_results <- list()
        posthoc_results <- list()
        
        for (cl in cluster_list) {
                message("Processing ", gene, " in ", cl)
                df <- meta %>% filter(as.character(.data[[cluster_col]]) == cl,
                                      !is.na(ILAE_group), is.finite(expression)) %>%
                        droplevels()
                n_groups <- n_distinct(df$ILAE_group)
                if (n_groups < 2 || n_distinct(df$expression) < 2) {
                        warning("Insufficient groups or expression variation for ", gene, " in ", cl)
                        next
                }
                
                main_test <- kruskal.test(expression ~ ILAE_group, data = df)
                posthoc_df <- NULL
                if (is.finite(main_test$p.value) && main_test$p.value < posthoc_p_cutoff) {
                        if (n_groups == 2) {
                                groups <- levels(droplevels(df$ILAE_group))
                                w <- wilcox.test(expression ~ ILAE_group, data = df,
                                                 exact = FALSE, correct = FALSE)
                                posthoc_df <- data.frame(group1 = groups[1], group2 = groups[2], P.adj = w$p.value)
                        } else {
                                dunn <- FSA::dunnTest(expression ~ ILAE_group, data = df, method = "bh")$res
                                split_groups <- strsplit(dunn$Comparison, " - ", fixed = TRUE)
                                posthoc_df <- data.frame(group1 = vapply(split_groups, `[`, "", 1),
                                                         group2 = vapply(split_groups, `[`, "", 2),
                                                         P.adj = dunn$P.adj)
                        }
                        posthoc_df$significance <- ifelse(posthoc_df$P.adj < 0.001, "***",
                                                          ifelse(posthoc_df$P.adj < 0.01, "**",
                                                                 ifelse(posthoc_df$P.adj < 0.05, "*", "ns")))
                }
                
                subtitle_text <- paste0("Kruskal-Wallis P = ", format(main_test$p.value, digits = 2, scientific = TRUE))
                if (!is.null(posthoc_df)) {
                        sig <- posthoc_df[!is.na(posthoc_df$P.adj) & posthoc_df$P.adj < posthoc_p_cutoff, ]
                        if (nrow(sig)) subtitle_text <- paste0(subtitle_text, "\nPost hoc: ",
                                                               paste0(sig$group1, " - ", sig$group2, " (", sig$significance, ")", collapse = "; "))
                }
                
                p <- ggplot(df, aes(x = ILAE_group, y = expression)) +
                        geom_boxplot(aes(fill = ILAE_group), outlier.shape = NA, alpha = 0.6) +
                        geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
                        scale_fill_manual(values = pal_ILAE_score, drop = FALSE) +
                        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
                        labs(title = paste(gene, "in", cl), subtitle = subtitle_text,
                             x = "ILAE score", y = "Pseudobulk expression (Xenium)") +
                        theme_classic() +
                        theme(legend.position = "none", plot.title = element_text(size = 8),
                              plot.subtitle = element_text(size = 7))
                
                plot_list[[cl]] <- p
                test_results[[cl]] <- main_test
                posthoc_results[[cl]] <- posthoc_df
        }
        
        combined_plot <- if (length(plot_list)) wrap_plots(plot_list, ncol = ncol) else NULL
        list(plot_list = plot_list, combined_plot = combined_plot,
             test_results = test_results, posthoc_results = posthoc_results)
}


### User-editable configuration -------------------------------------------------

# datIN_v8.rds is generated by script 01. The regional metadata file must have
# Xenium cell identifiers as row names and contain a region_md column.

cfg <- list(datIN_rds = "path/to/directory/data/datIN_v8.rds",
            region_metadata_rds = "path/to/directory/data/megadata_metadata_md_v8.rds")
        

for (path in c(cfg$datIN_rds, cfg$region_metadata_rds)) {
        if (!file.exists(path)) stop("Required input not found: ", path)
}

pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00",
                    "2" = "#00C0C5", "3" = "#E69F00")

### Load annotated data ------------------------------------------

datIN <- readRDS(cfg$datIN_rds)

required_meta <- c("seurat_clusters", "fine_subtype_InhN", "broad_subtype_InhN",
                   "ILAE_score", "orig.ident", "region_md")

missing_meta <- setdiff(required_meta, colnames(datIN@meta.data))

if (length(missing_meta)) stop("Missing metadata in datIN: ", paste(missing_meta, collapse = ", "))

### SST: CA1 and CA2-CA4 ------------------------------------------------------

# Subset cells before aggregating to obtain one pseudobulk profile per patient and region.
datIN_SST <- subset(datIN,subset = broad_subtype_InhN == "SST" & 
                            region_md %in% c("CA1", "CA2_CA4"))


# Pseudobulk per patient, ILAE score and region.
pseudo_datIN_SST <- AggregateExpression( datIN_SST, assays = "Xenium",
                                         return.seurat = TRUE,
                                         group.by = c("ILAE_score", "orig.ident", "region_md"))

# Inspect pseudobulk group sizes before interpreting statistical tests.
print(table(pseudo_datIN_SST$ILAE_score, pseudo_datIN_SST$region_md))


# Compare GABRA1, GABRA2 and SNCG across ILAE scores.
GABRA1_SST <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_SST,
                             gene = "GABRA1", cluster_col = "region_md",
                             ncol=4)

GABRA2_SST <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_SST,
                             gene = "GABRA2", cluster_col = "region_md",
                             ncol=4, type='SST')

SNCG_SST <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_SST,
                           gene = "SNCG", cluster_col = "region_md",
                           ncol=4, type='SST')


### PV: CA1 and CA2-CA4 -------------------------------------------------------

# Subset cells before aggregating to obtain one pseudobulk profile per patient and region.
datIN_PV <- subset(datIN,
                   subset = broad_subtype_InhN == "PV" & 
                           region_md %in% c("CA1", "CA2_CA4"))


# Pseudobulk per patient, ILAE score and region.
pseudo_datIN_PV <- AggregateExpression(
        datIN_PV,
        assays = "Xenium",
        return.seurat = TRUE,
        group.by = c("ILAE_score", "orig.ident", "region_md"))


print(table(pseudo_datIN_PV$ILAE_score, pseudo_datIN_PV$region_md))


GABRA1_PV <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_PV,
                            gene = "GABRA1", cluster_col = "region_md",
                            ncol=4, type='PV')

GABRA2_PV <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_PV,
                            gene = "GABRA2", cluster_col = "region_md",
                            ncol=4, type='PV')

SNCG_PV <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_PV,
                          gene = "SNCG", cluster_col = "region_md",
                          ncol=4, type='PV')




### VIP: SL/SR/SLM ------------------------------------------------------------

# Subset cells before aggregating to obtain one pseudobulk profile per patient and region.
datIN_VIP <- subset(datIN,
                    subset = broad_subtype_InhN == "VIP" & 
                            region_md == "SL.SR.SLM")


# Pseudobulk per patient, ILAE score and region.
pseudo_datIN_VIP <- AggregateExpression(
        datIN_VIP,
        assays = "Xenium",
        return.seurat = TRUE,
        group.by = c("ILAE_score", "orig.ident", "region_md"))


print(table(pseudo_datIN_VIP$ILAE_score, pseudo_datIN_VIP$region_md))


GABRA1_VIP <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_VIP,
                             gene = "GABRA1", cluster_col = "region_md",
                             ncol=4, type='VIP')

GABRA2_VIP <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_VIP,
                             gene = "GABRA2", cluster_col = "region_md",
                             ncol=4, type='VIP')

SNCG_VIP <- Test_Plot_DEGs(seurat_obj = pseudo_datIN_VIP,
                           gene = "SNCG", cluster_col = "region_md",
                           ncol=4, type='VIP')


### Combine selected-gene plots ------------------------------------------------

combined_plot <- wrap_plots(GABRA1_SST$combined_plot,
                            GABRA2_SST$combined_plot,
                            SNCG_SST$combined_plot,
                            GABRA1_PV$combined_plot,
                            GABRA2_PV$combined_plot,
                            SNCG_PV$combined_plot, 
                            GABRA1_VIP$combined_plot,
                            GABRA2_VIP$combined_plot,
                            SNCG_VIP$combined_plot, 
                            ncol = 3)
print(combined_plot)

### Inspect significant Dunn comparisons --------------------------------------

# The function stores BH-adjusted Dunn comparisons below the cutoff.
GABRA1_SST$posthoc_results
GABRA2_SST$posthoc_results
SNCG_SST$posthoc_results

GABRA1_PV$posthoc_results
GABRA2_PV$posthoc_results
SNCG_PV$posthoc_results

GABRA1_VIP$posthoc_results
GABRA2_VIP$posthoc_results
SNCG_VIP$posthoc_results

