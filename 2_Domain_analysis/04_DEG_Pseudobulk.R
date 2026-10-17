################################################################################
# 04_deg_pseudobulk.R
#
# Differential expression analysis on pseudobulked data, between ILAE scores
#  - Re-extracts "in-boundary" (somatic) vs "out-boundary" (extrasomatic)
#    transcripts per spatial bin from the raw Xenium molecule coordinates
#  - Normalizes the extrasomatic-only counts, merges into one Seurat object
#  - Pseudobulks by (ILAE_score, sample, region_md), runs DESeq2-based
#    FindMarkers for every pairwise ILAE comparison within each region
#  - Volcano plots + CSV export of significant DEGs, plus a summary bubble
#    plot of DEG counts/direction across neuropil regions
#  - Per-gene expression boxplots (with Kruskal-Wallis / Dunn's test) for a
#    curated gene list per region
#
# Original author: Mathias Delhaye, October 2026
# Code cleanup and standardization: Larissa Kraus (October 2026)

### 0) Setup -------------------------------------------------------------------

# Core packages (keep this list lean; add others only if used)
suppressPackageStartupMessages({
        library(conflicted)
        library(EnhancedVolcano)
        library(tidyverse)
        library(DESeq2)
        
        library(FSA)
        library(ggpubr)
        library(grid)
        library(colorspace)
        library(Matrix)
        library(future)
        library(future.apply)
        library(Seurat)
        library(ggh4x)
        library(ggtext)
        library(scales)
})

### PACKAGE CONFLICTS ---------------------------------------------------------
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)


### 1) User-editable config ----------------------------------------------------

cfg <- list(
        # INPUT: original annotated MegaData Seurat object (also used in script 01).
        # Required here for raw Xenium molecules, cell centroids and counts.
        megadata_file = NA_character_,
        
        # INPUT: binned Seurat object with region_md annotations, saved by
        # 01_domainIdentification.R (Section 8) as 07_seumerged.rds.
        seumerged_file = NA_character_,
        
        # OUTPUT: RDS results, figures and DEG CSV tables (directories created below).
        rds_dir = NA_character_,
        fig_dir = NA_character_,
        csv_dir = NA_character_,
        
        # Reuse valid Seurat objects already in the R environment if available.
        # Set FALSE to force loading from the input paths above.
        
        # Binning parameters (must match 01_domainIdentification.R)
        bin_width  = 150,
        bin_height = 150,
        
        # QC filter for the extrasomatic-only bins
        min_nCount_extrasoma = 5,
        
        # DEG thresholds
        fc_cutoff   = 1,
        padj_cutoff = 0.05,
        
        # Region + ILAE palettes (same as the other scripts)
        col_regions = setNames(
                c("#3D0F99", "#C7B8E6", "#0F8299", "#B8DEE6", "#54990F", "#CFE6B8",
                  "#333333", "#990F26", "#CCAA7A", "#99600F", "#FFD300"),
                c("GCL", "ML", "CA2_CA4", "CA1", "SL.SR.SLM", "GABA", "SUB", "WM", "pia", "vascular", "void")
        ),
        color_ILAE = setNames(
                c("#C67BFF", "#F8766C", "#7BAD00", "#00C0C5", "#E69F00"),
                c("P", "ILAE0", "ILAE1", "ILAE2", "ILAE3")
        ),
        ILAE_levels = c("P", "ILAE0", "ILAE1", "ILAE2", "ILAE3"),
        
        # Regions + comparison order for the DEG summary bubble plot
        neuropil_region  = c("ML", "CA2_CA4", "SL.SR.SLM", "GABA"),
        comparison_order = c(
                "P_vs_ILAE0", "P_vs_ILAE1", "P_vs_ILAE2", "P_vs_ILAE3",
                "ILAE0_vs_ILAE1", "ILAE0_vs_ILAE2", "ILAE0_vs_ILAE3",
                "ILAE1_vs_ILAE2", "ILAE1_vs_ILAE3",
                "ILAE2_vs_ILAE3"
        ),
        
        # Curated gene lists for the per-gene ILAE comparison plots (Section 8),
        # picked from the volcano plots / DEG table above
        list_genes_to_compare = list(
                ML        = c("ANO3", "MCTP2", "GABRG1"),
                SL.SR.SLM = c("ANO3", "MCTP2", "GABRG1"),
                CA2_CA4   = c("GABRG1"),
                GABA      = c("GABRG1")
        ),
        
        list_genes_to_compare_extended = list(
                ML = c("ANO3", "GABBR1", "GABRG1", "GPNMB", "GRIN2C", "NTS", "TRIL", "TSHZ2", "CALB1", "SLC17A6"),
                SL.SR.SLM = c("ANO3", "GABBR1", "GABRA2", "GABRG1", "GRIN2C", "GRM5", "MCTP2", "MEPE",
                              "MYO5B", "SLC17A6", "TOP2A", "TRHDE", "TRIL", "TRPC6", "TTYH1"),
                CA2_CA4 = c("ABCC9", "ANO3", "CALCRL", "GABBR1", "GABRA2", "GABRG1", "GRM5", "GRIN1",
                            "GRIN2B", "GRIN2C", "MCTP2", "PVALB", "RIT2", "SLC24A3", "TRHDE", "TRIL",
                            "TRPC6", "TSHZ2", "TTYH1"),
                GABA = c("ANO3", "GRIN2C", "GABRA2", "GABRG1", "MGST1", "MCTP2", "RIMS1", "RIMS2", "TRIL", "CALB2")
        ),
        
        # Output filenames
        f_insoma_persoma        = "10_inboundaries_transcripts_persoma.rds",
        f_insoma_perbin         = "11_inboundaries_transcripts_perbin.rds",
        f_seulist_extrasoma     = "12_seulist_extraboundaries.rds",
        f_seulist_extrasoma_log = "13_seulist_extraboundaries_log5.rds",
        f_seumerged_extrasoma   = "14_seumerged_extraboundaries.rds",
        f_degtable              = "15_degtable.rds",
        f_fig_deg_bubble        = "07_deg_bubble_summary.png",
        f_stats_table           = "ILAE_gene_comparison_stats.csv"
)

# Create output directories
for (dir in c(cfg$rds_dir, cfg$fig_dir, cfg$csv_dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        if (!dir.exists(dir)) stop("Cannot create output directory: ", dir)
}

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


# 2.1 Build a Seurat object of extrasomatic-only transcript counts for one
#     sample: (total transcripts in bin) - (transcripts already assigned to a
#     soma in that bin), floored at 0.
build_extrasoma_seurat <- function(im, gene_list, megadata,
                                   insoma_bin_counts,
                                   bin_width = 150, bin_height = 150) {
        
        fov <- megadata@images[[im]]@molecules$molecules
        
        # Pull every gene's coords once, combine with a single rbind (not one per gene)
        coords_list <- lapply(gene_list, function(g) fov[[g]]@coords)
        n_per_gene  <- vapply(coords_list, nrow, integer(1))
        all_coords  <- do.call(rbind, coords_list)
        gene_vec    <- rep(gene_list, n_per_gene)
        
        bin_x <- ceiling(all_coords[, 1] / bin_width)
        bin_y <- ceiling(all_coords[, 2] / bin_height)
        
        number_bin_x <- ceiling(max(all_coords[, 1]) / bin_width)
        number_bin_y <- ceiling(max(all_coords[, 2]) / bin_height)
        n_bins <- number_bin_x * number_bin_y
        
        bin_number <- (bin_y - 1) * number_bin_x + bin_x
        gene_idx   <- match(gene_vec, gene_list)
        bin_names  <- paste0(im, "_", seq_len(n_bins))
        
        # Total transcripts per gene x bin (sparse; duplicate i,j pairs are summed automatically)
        mat_total <- sparseMatrix(
                i = gene_idx, j = bin_number, x = 1,
                dims = c(length(gene_list), n_bins),
                dimnames = list(gene_list, bin_names)
        )
        
        # Matching in-soma counts, aligned to the same gene/bin dimensions
        insoma_sub <- insoma_bin_counts[insoma_bin_counts$bin_number %in% bin_names, ]
        mat_insoma <- sparseMatrix(
                i = match(insoma_sub$gene, gene_list),
                j = match(insoma_sub$bin_number, bin_names),
                x = insoma_sub$n_count,   # <-- confirm this is the right column name
                dims = c(length(gene_list), n_bins),
                dimnames = list(gene_list, bin_names)
        )
        
        mat_extrasoma <- mat_total - mat_insoma
        mat_extrasoma[mat_extrasoma < 0] <- 0  # safety guard, see note below
        
        # col fastest, row slowest -- matches your original bin_number formula
        bin_grid <- expand.grid(col = 1:number_bin_x, row = 1:number_bin_y)
        rownames(bin_grid) <- bin_names
        
        CreateSeuratObject(counts = mat_extrasoma, assay = "RNA", meta.data = bin_grid)
}

# 2.2 For every region, run FindMarkers (DESeq2) on every pairwise ILAE
#     comparison, and bind the results into one long table
compute_deg_table <- function(pseudo,
                              cat = "ILAE_score",
                              cluster = "region_md",
                              regions = NULL,
                              group_sep = "_",
                              fc_cutoff = 1,
                              padj_cutoff = 0.05,
                              verbose = TRUE) {
        
        ILAE_score  <- levels(pseudo[[cat]][, 1])
        if (is.null(regions)) regions <- unique(pseudo[[cluster]][, 1])
        comb_idents <- combn(ILAE_score, 2)
        
        results_list <- list()
        
        for (r in regions) {
                temp_mat <- apply(comb_idents, 2, function(col) paste(col, r, sep = group_sep))
                
                for (i in seq_len(ncol(comb_idents))) {
                        ident1 <- temp_mat[1, i]
                        ident2 <- temp_mat[2, i]
                        comparison <- paste0(comb_idents[1, i], "_vs_", comb_idents[2, i])
                        
                        if (verbose) message("DE: ", ident1, " vs ", ident2)
                        
                        de <- tryCatch(
                                FindMarkers(object = pseudo, ident.1 = ident1, ident.2 = ident2,
                                            test.use = "DESeq2"),
                                error = function(e) {
                                        warning("Skipped ", ident1, " vs ", ident2, ": ", conditionMessage(e))
                                        NULL
                                }
                        )
                        if (is.null(de)) next
                        
                        de <- de %>%
                                tibble::rownames_to_column(var = "gene") %>%
                                mutate(
                                        region     = r,
                                        group1     = ident1,
                                        group2     = ident2,
                                        comparison = comparison,
                                        direction  = case_when(
                                                avg_log2FC >  fc_cutoff & p_val_adj < padj_cutoff ~ "up",
                                                avg_log2FC < -fc_cutoff & p_val_adj < padj_cutoff ~ "down",
                                                TRUE ~ "ns"
                                        ),
                                        significant = direction != "ns"
                                )
                        
                        results_list[[paste(r, ident1, ident2, sep = "|")]] <- de
                }
        }
        
        bind_rows(results_list)
}

# 2.3 Volcano plot for every region x comparison in `deg_table`, saved to disk
plot_deg_volcanoes <- function(deg_table, output.dir,
                               pCutoff = 0.05, labSize = 4, w = 6, h = 6) {
        
        for (r in unique(deg_table$region)) {
                output_save <- file.path(output.dir, r)
                dir.create(output_save, recursive = TRUE, showWarnings = FALSE)
                
                combos <- deg_table %>% filter(region == r) %>% distinct(group1, group2)
                
                for (j in seq_len(nrow(combos))) {
                        g1 <- combos$group1[j]; g2 <- combos$group2[j]
                        
                        sub_de <- deg_table %>%
                                filter(region == r, group1 == g1, group2 == g2) %>%
                                column_to_rownames("gene")
                        
                        p <- EnhancedVolcano(sub_de,
                                             lab = rownames(sub_de),
                                             x = "avg_log2FC", y = "p_val_adj",
                                             drawConnectors = TRUE,
                                             subtitle = NULL,
                                             captionLabSize = 12,
                                             legendLabSize = 10,
                                             titleLabSize = 14,
                                             pCutoff = pCutoff) +
                                ggtitle(paste(g1, "vs", g2, " - ", r))
                        
                        filename <- file.path(output_save, paste0("Volcano_", g1, "_vs_", g2, "_region_", r, ".png"))
                        
                        ggsave(filename = filename, plot = p, dpi = 300, width = w, height = h)
                }
        }
}

# 2.4 Export significant DEGs (one CSV per region) from `deg_table`
export_deg_csv <- function(deg_table, output.dir) {
        dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
        
        for (r in unique(deg_table$region)) {
                df <- deg_table %>%
                        filter(region == r, significant) %>%
                        select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, group1, group2)
                
                write.csv(df, file.path(output.dir, paste0("table_deg_", r, ".csv")), row.names = FALSE)
        }
}

# 2.5 Expression of one gene, in one region, across ILAE scores: boxplot +
#     Kruskal-Wallis / Dunn's test (BH-adjusted), same pattern as
#     compare_groups_boxplot() in 03_region_proportions.R but pulling
#     pseudobulk expression via FetchData() instead of spot proportions.
ILAE_comparison_pseudobulk <- function(pseudo, region_col = "region_md", loc, gene, color, l.width = 0.2) {
        
        # get the data for the wanted gene only
        df <- FetchData(pseudo, vars = c(gene, region_col, "ILAE_score"), assay = Seurat::Assays(pseudo))
        names(df)[1] <- "expression"
        
        # isolate for the specific region
        df_sub <- filter(df, .data[[region_col]] == loc)
        
        kw <- kruskal.test(df_sub$expression, g = df_sub$ILAE_score)
        dunn <- dunnTest(expression ~ ILAE_score, data = df_sub, method = "bh")
        dunn <- dunn$res %>%
                mutate(significance = case_when(
                        P.adj < 0.001 ~ "***",
                        P.adj < 0.01  ~ "**",
                        P.adj < 0.05  ~ "*",
                        TRUE ~ "ns"
                )) %>%
                separate(Comparison, into = c("group1", "group2"), sep = " - ")
        dunn_plot <- dunn %>% filter(P.adj < 0.05)
        
        max_df <- max(df_sub$expression)
        y_limit <- max_df + 0.3 * 4  # should accommodate 4 significant comparisons
        kw_pval <- paste("p.val", as.character(round(kw$p.value, 4)))
        
        p <- ggplot(df_sub, aes(x = ILAE_score, y = expression, fill = ILAE_score)) +
                geom_boxplot(lwd = l.width, colour = "black", outlier.shape = NA) + # can be adjusted
                theme_minimal() +
                theme(
                        axis.line = element_line(colour = "black"),
                        axis.ticks = element_line(colour = "black"),
                        axis.text.x = element_text(colour = "black",
                                                   angle = 45,
                                                   vjust = 1,
                                                   hjust = 1),
                        axis.text.y = element_text(colour = "black", size = 9),
                        panel.grid = element_blank(),
                        legend.position = "none",
                        axis.title = element_blank(),
                        title = element_text(colour = "black", size = 9)
                ) +
                scale_fill_manual(values = color) +
                geom_jitter(color = "black", size = 0, alpha = 0.9, width = 0.15) +
                labs(title = paste0(loc, " - ", kw_pval)) +
                coord_cartesian(ylim = c(0, y_limit), clip = "off")
        
        if (nrow(dunn_plot) > 0) {
                position_comparison <- seq(from = max_df + 0.3, to = max_df + nrow(dunn_plot) * 0.3, by = 0.3)
                p <- p + stat_pvalue_manual(dunn_plot, label = "significance", tip.length = 0.01, y.position = position_comparison)
        }
        
        list(p = p, kw = kw, dunn = dunn)
}

# 2.6 Run ILAE_comparison_pseudobulk() for every gene in every region of
#     `list_genes`, saving one plot per gene
export_gene_comparison <- function(list_genes, pseudo, w = 2, h = 5, output_dir) {
        for (region in names(list_genes)) {
                save_path <- file.path(output_dir, region)
                dir.create(path = save_path, recursive = TRUE, showWarnings = FALSE)
                genes <- list_genes[[region]]
                
                for (g in genes) {
                        gene_comp <- ILAE_comparison_pseudobulk(pseudo = pseudo, loc = region, gene = g, color = cfg$color_ILAE)
                        filename <- file.path(save_path, paste0("expression_", g, "_", region, ".png"))
                        ggsave(filename = filename, plot = gene_comp$p, width = w, height = h, dpi = 300)
                }
        }
}

# 2.7 Same gene x region loop as export_gene_comparison(), but collecting the
#     Kruskal-Wallis / Dunn's-test stats into one table instead of saving plots
get_ILAE_stats_table <- function(list_genes, pseudo, region_col = "region_md") {
        
        results <- list()
        
        for (region in names(list_genes)) {
                genes <- list_genes[[region]]
                for (g in genes) {
                        
                        # reuse the same function as export_gene_comparison() (plot is discarded here)
                        res <- ILAE_comparison_pseudobulk(
                                pseudo = pseudo,
                                loc = region,
                                gene = g,
                                color = cfg$color_ILAE  # not used for stats, but required by function signature
                        )
                        
                        dunn_df <- res$dunn %>%
                                mutate(
                                        region = region,
                                        gene = g,
                                        kruskal_pval = res$kw$p.value
                                ) %>%
                                select(region, gene, group1, group2, P.unadj, P.adj, significance, kruskal_pval)
                        
                        results[[paste(region, g, sep = "_")]] <- dunn_df
                }
        }
        
        bind_rows(results)
}

### 3) Load inputs ---------------------------------------------------------------

# MegaData: original input also used by 01_domainIdentification.R.
# Script 04 requires raw molecule coordinates, cell centroids and Xenium counts.
megadata <- load_if_missing("megadata", cfg$megadata_file, class = "Seurat")

# seumerged: output 07_seumerged.rds from 01_domainIdentification.R.
# Contains 150 x 150 um spatial bins and region / manually corrected region_md.
seumerged <- load_if_missing("seumerged", cfg$seumerged_file, class = "Seurat")

required_cols <- c("orig.ident", "col", "row", "smpl_ILAE", "ILAE_score", "cluster_precast", "region", "region_md")
missing_cols <- setdiff(required_cols, colnames(seumerged@meta.data))
if (length(missing_cols)) stop("seumerged missing metadata columns: ", paste(missing_cols, collapse = ", "))
if (!"Xenium" %in% Assays(megadata)) stop("megadata is missing the Xenium assay")
if (!length(megadata@images)) stop("megadata has no images; raw Xenium coordinates are required")

samples_names <- names(megadata@images)

# Genes analyzed with Xenium (266 panel + 100 epilepsy-specific)
genes_list <- Features(megadata)

# Merge all per-sample layers into a single counts layer
megadata[["Xenium"]] <- JoinLayers(megadata[["Xenium"]])

### 4) Extract in-boundary (in-soma) transcripts per bin --------------------------

# Pull the full counts matrix once
xen_count <- LayerData(megadata, assay = "Xenium", layer = "counts")

# Build the transcript coordinate table in one shot
temp_df <- as.data.frame(summary(xen_count))
cell_names <- colnames(xen_count)  # use colnames/rownames directly, not gene_list
gene_names <- rownames(xen_count)

insoma_transcript_coord <- tibble(
        gene    = gene_names[temp_df$i],
        cell    = cell_names[temp_df$j],
        n_count = temp_df$x
)

metadata_seumerged <- seumerged@meta.data

df_bin_number <- metadata_seumerged %>%
        rownames_to_column(var = "bin_number") %>%
        select(orig.ident, col, row, bin_number)

cell_coord <- map_dfr(samples_names, function(smpl) {
        coord <- megadata@images[[smpl]]@boundaries[["centroids"]]@coords
        cell_names <- megadata@images[[smpl]]@boundaries[["centroids"]]@cells
        tibble(cell = cell_names, x = coord[, 1], y = coord[, 2])
})

# Add column for the corresponding col/row coordinates, using the same bin
# size the bins were originally built with (ceiling since bins are
# bin_width x bin_height and we want each soma's bin index, not its raw position)
cell_coord <- cell_coord %>%
        mutate(col = ceiling(x / cfg$bin_width),
               row = ceiling(y / cfg$bin_height))

# Extract the sample name for each soma using the metadata of megadata
metadata_megadata <- megadata@meta.data
cell_to_orig <- metadata_megadata %>%
        rownames_to_column(var = "cell") %>%
        select(cell, orig.ident)  # keep only the cell name and the sample name
cell_coord <- left_join(cell_coord, cell_to_orig, by = "cell")  # add the sample name to the cell_coord table
rm(cell_to_orig, metadata_megadata)

# Join cell_coord and metadata_seumerged by orig.ident, col and row to get the
# bin (and therefore region/region_md) where each cell belongs to
cell_coord <- cell_coord %>%
        left_join(df_bin_number, by = c("orig.ident", "col", "row"))

insoma_transcript_coord <- insoma_transcript_coord %>%
        left_join(cell_coord %>% select(cell, bin_number), by = "cell")

insoma_transcript_coord_bin <- insoma_transcript_coord %>%
        group_by(bin_number, gene) %>%
        summarise(n_count = sum(n_count), .groups = "drop")

saveRDS(insoma_transcript_coord, file.path(cfg$rds_dir, cfg$f_insoma_persoma))
saveRDS(insoma_transcript_coord_bin, file.path(cfg$rds_dir, cfg$f_insoma_perbin))
message("Saved in-boundary transcript tables to: ", cfg$rds_dir)

### 5) Isolate extrasomatic-only transcripts, normalize, merge -------------------

seulist_extrasoma <- setNames(
        lapply(samples_names, build_extrasoma_seurat,
               gene_list = genes_list,
               megadata = megadata,
               insoma_bin_counts = insoma_transcript_coord_bin,
               bin_width = cfg$bin_width, bin_height = cfg$bin_height),
        samples_names
)

saveRDS(seulist_extrasoma, file.path(cfg$rds_dir, cfg$f_seulist_extrasoma))

message("Normalizing extrasomatic-only counts (nCount_RNA > ", cfg$min_nCount_extrasoma, ") in parallel")

plan(multisession, workers = parallel::detectCores() - 1)  # adjust worker count as needed

seulist_extrasoma_log5 <- future_lapply(seulist_extrasoma[samples_names], function(seu) {
        seu %>%
                subset(subset = nCount_RNA > cfg$min_nCount_extrasoma) %>%
                NormalizeData(verbose = FALSE) %>%
                FindVariableFeatures(verbose = FALSE) %>%
                ScaleData(verbose = FALSE)
}, future.seed = TRUE)

names(seulist_extrasoma_log5) <- samples_names
plan(sequential)  # reset back to sequential when done
saveRDS(seulist_extrasoma_log5, file.path(cfg$rds_dir, cfg$f_seulist_extrasoma_log))

seumerged_extrasoma <- merge(seulist_extrasoma_log5[[1]], seulist_extrasoma_log5[-1])

metadata_seumerged_extrasoma <- seumerged_extrasoma@meta.data

df_ILAE_reg <- metadata_seumerged %>%
        select(smpl_ILAE, ILAE_score, cluster_precast, region, region_md) %>%
        rownames_to_column(var = "bin_name")

metadata_seumerged_extrasoma <- metadata_seumerged_extrasoma %>%
        rownames_to_column(var = "bin_name") %>%
        left_join(df_ILAE_reg, by = "bin_name") %>%
        column_to_rownames(var = "bin_name")

seumerged_extrasoma@meta.data <- metadata_seumerged_extrasoma

saveRDS(seumerged_extrasoma, file.path(cfg$rds_dir, cfg$f_seumerged_extrasoma))
message("Saved extrasomatic seumerged to: ", file.path(cfg$rds_dir, cfg$f_seumerged_extrasoma))

### 6) Pseudobulk aggregation ------------------------------------------------------

# seumerged_extrasoma <- readRDS(file.path(cfg$rds_dir, cfg$f_seumerged_extrasoma))  # optional: reload if resuming from here

# Join layers for Seurat functions
seumerged_extrasoma_joined <- JoinLayers(seumerged_extrasoma)

# Removing void (and NAs)
seumerged_extrasoma_joined <- subset(x = seumerged_extrasoma_joined, subset = region != "void")

# Script 01 stores ILAE_score as P/0/1/2/3; standardize to the labels
# used by this script BEFORE AggregateExpression creates group identities.
seumerged_extrasoma_joined$ILAE_score <- dplyr::recode(as.character(seumerged_extrasoma_joined$ILAE_score),
                                                       "0" = "ILAE0", "1" = "ILAE1", "2" = "ILAE2", "3" = "ILAE3",
                                                       "g0" = "ILAE0", "g1" = "ILAE1", "g2" = "ILAE2", "g3" = "ILAE3")
seumerged_extrasoma_joined$ILAE_score <- factor(seumerged_extrasoma_joined$ILAE_score, levels = cfg$ILAE_levels)

# Use ILAE_score (or ILAE_LK, which has Mathias' ILAE scores for FF samples)
pseudo_seu_extrasoma <- AggregateExpression(
        seumerged_extrasoma_joined, assays = "RNA", return.seurat = TRUE,
        group.by = c("ILAE_score", "orig.ident", "region_md")
)

pseudo_seu_extrasoma@meta.data <- pseudo_seu_extrasoma@meta.data %>%
        mutate(
                region_md = case_when(
                        region_md == "CA2-CA4" ~ "CA2_CA4",
                        TRUE ~ region_md
                ),
                ILAE_score = case_when(
                        ILAE_score %in% c("0", "g0") ~ "ILAE0",
                        ILAE_score %in% c("1", "g1") ~ "ILAE1",
                        ILAE_score %in% c("2", "g2") ~ "ILAE2",
                        ILAE_score %in% c("3", "g3") ~ "ILAE3",
                        TRUE ~ ILAE_score
                ),
                ILAE_score = factor(ILAE_score, levels = cfg$ILAE_levels)
        )

# Check data to see if grouping makes sense
message("Pseudobulk groups (tail): ", paste(tail(Cells(pseudo_seu_extrasoma)), collapse = ", "))

# Add metadata column which will be a combination of region and ILAE grade
pseudo_seu_extrasoma$ILAE_reg <- paste(pseudo_seu_extrasoma$ILAE_score, pseudo_seu_extrasoma$region_md, sep = "_")

# Set Idents() to the new ILAE_reg combo to test DE genes within one region between ILAE grades
Idents(pseudo_seu_extrasoma) <- "ILAE_reg"

### 7) DEG analysis on pseudobulk data ---------------------------------------------

df_deg <- compute_deg_table(
        pseudo = pseudo_seu_extrasoma,
        cat = "ILAE_score",
        cluster = "region_md",
        fc_cutoff = cfg$fc_cutoff,
        padj_cutoff = cfg$padj_cutoff
)

saveRDS(df_deg, file = file.path(cfg$rds_dir, cfg$f_degtable))
message("Saved DEG table to: ", file.path(cfg$rds_dir, cfg$f_degtable))

# Volcano plots, one per region x ILAE comparison
volcano_dir <- file.path(cfg$fig_dir,"01_volcano_plots")
plot_deg_volcanoes(df_deg, output.dir = volcano_dir)

# CSV of significant DEGs, one per region
export_deg_csv(df_deg, cfg$csv_dir)

# Summarize number of significant up/down DEGs per region x comparison, for
# the neuropil-rich regions only
deg_summary <- df_deg %>%
        filter(region %in% cfg$neuropil_region, significant) %>%
        group_by(region, comparison, direction) %>%
        summarise(DEG_count = n(), .groups = "drop") %>%
        mutate(direction = case_when(
                direction == "up"   ~ "Up",
                direction == "down" ~ "Down",
                TRUE ~ direction
        )) %>%
        pivot_wider(names_from = direction, values_from = DEG_count, values_fill = 0) %>%
        mutate(
                Total = Up + Down,
                Direction_balance = if_else(Total > 0, (Up - Down) / Total, 0)
        ) %>%
        mutate(
                region     = factor(region, levels = rev(cfg$neuropil_region)),
                comparison = factor(comparison, levels = cfg$comparison_order)
        )

# Rich-text axis labels: "P" in green, the comparison ILAE grade in purple
# (kept explicit, as in the original - a regex/sprintf version would be more
# "DRY" but is genuinely harder to verify at a glance than just listing them)
comparison_labels <- c(
        "P_vs_ILAE0" = "<b><span style='color:#1B7837;'>P</span><br>vs<br><span style='color:#762A83;'>0</span></b>",
        "P_vs_ILAE1" = "<b><span style='color:#1B7837;'>P</span><br>vs<br><span style='color:#762A83;'>1</span></b>",
        "P_vs_ILAE2" = "<b><span style='color:#1B7837;'>P</span><br>vs<br><span style='color:#762A83;'>2</span></b>",
        "P_vs_ILAE3" = "<b><span style='color:#1B7837;'>P</span><br>vs<br><span style='color:#762A83;'>3</span></b>",
        "ILAE0_vs_ILAE1" = "<b><span style='color:#1B7837;'>0</span><br>vs<br><span style='color:#762A83;'>1</span></b>",
        "ILAE0_vs_ILAE2" = "<b><span style='color:#1B7837;'>0</span><br>vs<br><span style='color:#762A83;'>2</span></b>",
        "ILAE0_vs_ILAE3" = "<b><span style='color:#1B7837;'>0</span><br>vs<br><span style='color:#762A83;'>3</span></b>",
        "ILAE1_vs_ILAE2" = "<b><span style='color:#1B7837;'>1</span><br>vs<br><span style='color:#762A83;'>2</span></b>",
        "ILAE1_vs_ILAE3" = "<b><span style='color:#1B7837;'>1</span><br>vs<br><span style='color:#762A83;'>3</span></b>",
        "ILAE2_vs_ILAE3" = "<b><span style='color:#1B7837;'>2</span><br>vs<br><span style='color:#762A83;'>3</span></b>"
)

p_bubble <- ggplot(
        deg_summary,
        aes(x = comparison, y = region, size = Total, color = Direction_balance)
) +
        geom_point(alpha = 0.9) +
        scale_size_area(max_size = 16, breaks = scales::pretty_breaks(n = 4), name = "Total DEGs") +
        scale_color_gradient2(
                low = "#762A83", mid = "#D8C99B", high = "#1B7837",
                midpoint = 0, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                name = "Direction of\nregulation"
        ) +
        scale_x_discrete(labels = comparison_labels) +
        guides(
                size = guide_legend(order = 1, override.aes = list(shape = 1, colour = "black")),
                color = guide_colorbar(order = 2)
        ) +
        theme_classic(base_size = 12) +
        theme(
                axis.title = element_blank(),
                axis.text.x = ggtext::element_markdown(hjust = 0.5, vjust = 1),
                axis.text.y = element_text(color = "black", angle = 90),
                axis.line = element_line(color = "black", linewidth = 0.5),
                axis.ticks = element_line(color = "black", linewidth = 0.5),
                panel.grid = element_blank(),
                legend.position = "right",
                plot.title = element_text(face = "bold"),
                plot.subtitle = element_text(color = "grey30")
        )

out_dir <- file.path(cfg$fig_dir,"02_deg_number")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(filename = file.path(out_dir, cfg$f_fig_deg_bubble), plot = p_bubble, width = 7, height = 6, dpi = 300)
message("Saved DEG summary bubble plot to: ", file.path(out_dir, cfg$f_fig_deg_bubble))

### 8) Per-gene ILAE comparisons for curated genes ---------------------------------

gene_comparison_dir <- file.path(cfg$fig_dir, "03_boxplot_gene_comparison")

export_gene_comparison(
        list_genes = cfg$list_genes_to_compare,
        pseudo = pseudo_seu_extrasoma,
        w = 2, h = 5,
        output_dir = gene_comparison_dir
)

stats_table <- get_ILAE_stats_table(cfg$list_genes_to_compare, pseudo = pseudo_seu_extrasoma)

write.csv(stats_table, file.path(cfg$csv_dir, cfg$f_stats_table), row.names = FALSE)
message("Saved gene-comparison stats table to: ", file.path(cfg$csv_dir, cfg$f_stats_table))
