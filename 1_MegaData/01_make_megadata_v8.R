################################################################################
#
# Xenium MegaData_v8 (Seurat v5)
#
# This script shows how we made the integrated full dataset with steps:
#   1) Load per-sample Xenium Seurat objects
#   2) Filter low-count cells
#   3) Add sample metadata (from a CSV)
#   4) Merge into one Seurat object with orig.ident-prefixed cell IDs
#   5) Run SingleR using Thompson snRNA-seq reference at:
#        a) mid.cell.class
#        b) superfine.cell.class
#   6) SCTransform + PCA + Harmony integration
#   7) test Clusters + UMAP
#   8) run final clustering and adds metadata columns 'celltypes'
#   9) Transfer subtype and region annotations
#
# this script is intended to help users understand how we made 'Megadata', not to be run as is
# Individual datasets available on Geo do not include all files needed for LoadXenium() 
# instead we encourage users to use 'Megadata_v8_annotated.rds' for analysis
# Author: Larissa Kraus (October 2026) 
######

### 0) Setup -------------------------------------------------------------------

# Core packages 
suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(patchwork)
        library(SingleR)
        library(future)
        library(rlang)
        library(scCustomize) 
        library(beepr)       # optional: mario sound
        
})

### PACKAGE CONFLICTS
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)

# Memory for large objects (adjust if needed)
options(future.globals.maxSize = 100000 * 1024^2)


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


# Cluster + UMAP convenience wrapper
Xenium_SeuratLK <- function(x,
                            DIMs = 1:15,
                            cluster.resolution = 0.5,
                            red = "harmony",
                            playMario = FALSE) {
        x <- FindNeighbors(x, dims = DIMs, reduction = red)
        x <- FindClusters(x, resolution = cluster.resolution)
        x <- RunUMAP(x, dims = DIMs, reduction = red)
        
        if (isTRUE(playMario)) {
                beepr::beep(sound = "mario")
        }
        x
}

# Load + annotate one sample Seurat object
load_one_xenium <- function(rds_path, sample_id, metadata_df, min_nCount_Xenium = 10) {
        obj <- readRDS(rds_path)
        
        # Ensure sample identity is set
        obj$orig.ident <- sample_id
        
        # Basic QC filter
        if (!"nCount_Xenium" %in% colnames(obj@meta.data)) {
                stop("Expected nCount_Xenium in meta.data but didn't find it in: ", rds_path)
        }
        obj <- subset(obj, subset = nCount_Xenium > min_nCount_Xenium)
        
        # Join metadata by Sample column (must match sample_id)
        if (!"Sample" %in% colnames(metadata_df)) {
                stop("metadata_df must contain a column named 'Sample'.")
        }
        sample_meta <- metadata_df %>% filter(Sample == sample_id)
        if (nrow(sample_meta) != 1) {
                warning("Metadata row for sample '", sample_id, "' is not exactly 1 row (n=", nrow(sample_meta), ").")
        }
        
        # Add each metadata column into Seurat meta.data (recycled across cells)
        for (col in colnames(sample_meta)) {
                obj[[col]] <- sample_meta[[col]][1]
        }
        
        obj
}

# Run SingleR per-sample, save per-sample results, and return full-label vectors
run_singler_by_sample <- function(dat,
                                  ref_seurat,
                                  ref_assay = "originalexp",
                                  ref_label_col,
                                  out_dir,
                                  query_assay = "Xenium") {
        
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Normalize reference
        ref_seurat <- NormalizeData(ref_seurat, assay = ref_assay, verbose = FALSE)
        
        # Pseudobulk reference by label column -> genes x labels matrix
        pseudobulk <- AggregateExpression(
                ref_seurat,
                group.by = ref_label_col,
                assays = ref_assay,
                slot = "data",
                return.seurat = FALSE
        )
        ref_expr <- pseudobulk[[ref_assay]]
        
        sample_names <- unique(dat$orig.ident)
        
        all_labels <- rep(NA_character_, ncol(dat))
        all_pruned <- rep(NA_character_, ncol(dat))
        names(all_labels) <- colnames(dat)
        names(all_pruned) <- colnames(dat)
        
        for (sample_name in sample_names) {
                message("[SingleR] Running sample: ", sample_name)
                
                sample_obj <- subset(dat, subset = orig.ident == sample_name)
                sample_obj <- NormalizeData(sample_obj, assay = query_assay, verbose = FALSE)
                
                query_expr <- GetAssayData(sample_obj, slot = "data", assay = query_assay)
                
                common_genes <- intersect(rownames(ref_expr), rownames(query_expr))
                ref_sub <- ref_expr[common_genes, , drop = FALSE]
                query_sub <- query_expr[common_genes, , drop = FALSE]
                
                singleR_results <- SingleR(
                        test = as.matrix(query_sub),
                        ref  = as.matrix(ref_sub),
                        labels = colnames(ref_expr)
                )
                
                saveRDS(singleR_results, file = file.path(out_dir, paste0("SingleR_", sample_name, "_", ref_label_col, ".rds")))
                
                sample_cells <- colnames(sample_obj)
                all_labels[sample_cells] <- singleR_results$labels
                all_pruned[sample_cells] <- singleR_results$pruned.labels
        }
        
        list(labels = all_labels, pruned = all_pruned)
}

# test clustering
dim_reso_test <- function(dataset, dims_list, 
                          resolutions, 
                          output_dir, 
                          red,
                          group_vars = c('orig.ident', 'ILAE_score', 'predicted.celltype', 'FFPE_FF', 'type')) {
        
        if (!dir.exists(output_dir)) {
                dir.create(output_dir, recursive = TRUE)
        }
        
        for (dims in dims_list) {
                max_dim <- max(dims)
                
                dat_mod <- dataset
                dat_mod <- FindNeighbors(dat_mod, dims = dims, reduction = red)
                dat_mod <- RunUMAP(dat_mod, dims = dims, reduction = red)
                
                for (res in resolutions) {
                        dat_res <- FindClusters(dat_mod, resolution = res)
                        
                        plots <- list(DimPlot(dat_res, shuffle = TRUE, label = TRUE) + coord_fixed())
                        
                        for (group_var in group_vars) {
                                plots <- append(
                                        plots,
                                        list(DimPlot(dat_res, group.by = group_var, shuffle = TRUE) + coord_fixed())
                                )
                        }
                        
                        combined_plot <- wrap_plots(plots, ncol = 1)
                        
                        plot_title <- paste("Dimensions:", max_dim, "| Resolution:", res)
                        combined_plot <- combined_plot + plot_annotation(title = plot_title)
                        
                        file_name <- paste0("Dim_", max_dim, "_Res_", res, ".pdf")
                        pdf_path <- file.path(output_dir, file_name)
                        
                        pdf_width <- 8
                        pdf_height <- 5 * length(plots)
                        
                        pdf(pdf_path, width = pdf_width, height = pdf_height)
                        print(combined_plot)
                        dev.off()
                }
        }
}


### 2) User-editable paths and load inputs -------------------------------------
# Fill in paths before running; existing objects are reused where applicable.
cfg <- list(xenium_rds_dir = NA_character_, 
            metadata_csv = NA_character_,
            thompson_rds = NA_character_, 
            out_singler_dir = NA_character_,
            exn_rds = NA_character_, 
            inn_rds = NA_character_, 
            micro_rds = NA_character_,
            domain_annotations_rds = NA_character_)
if (is.na(cfg$xenium_rds_dir) || !dir.exists(cfg$xenium_rds_dir)) stop("Set cfg$xenium_rds_dir")
xenium_files <- list.files(cfg$xenium_rds_dir, pattern = "\\.rds$", full.names = FALSE)
sample_ids <- sub("\\.rds$", "", xenium_files)
if (!length(sample_ids)) stop("No sample RDS files found in cfg$xenium_rds_dir")
metadata <- if (exists("metadata", envir = .GlobalEnv, inherits = FALSE)) metadata else {
        if (is.na(cfg$metadata_csv) || !file.exists(cfg$metadata_csv)) stop("Set cfg$metadata_csv")
        read.csv(cfg$metadata_csv)
}
if (!"Sample" %in% names(metadata)) stop("metadata must contain Sample")
missing_meta <- setdiff(sample_ids, metadata$Sample)
if (length(missing_meta)) warning("Missing sample metadata: ", paste(missing_meta, collapse = ", "))
sc.ref <- load_if_missing("sc.ref", cfg$thompson_rds, class = "Seurat")
if (is.na(cfg$out_singler_dir)) stop("Set cfg$out_singler_dir")
dir.create(cfg$out_singler_dir, recursive = TRUE, showWarnings = FALSE)
xenium_rds_dir <- cfg$xenium_rds_dir

### 3) Load + annotate all Xenium objects -------------------------------------
# Load individual rds files
message("Loading Xenium samples from: ", xenium_rds_dir)

obj_list <- purrr::map2(
        xenium_files,
        sample_ids,
        ~ load_one_xenium(
                rds_path = file.path(xenium_rds_dir, .x),
                sample_id = .y,
                metadata_df = metadata,
                min_nCount_Xenium = 10
        )
)
names(obj_list) <- sample_ids

# Optional: sort by total UMIs so the biggest object is first (merge speed/memory)
umi_counts <- purrr::map_dbl(obj_list, ~ sum(.x@meta.data$nCount_Xenium, na.rm = TRUE))
obj_list <- obj_list[names(sort(umi_counts, decreasing = TRUE))]

# Merge with informative cell ID prefixes (orig.ident)
add_ids <- names(obj_list)
dat <- merge(
        x = obj_list[[1]],
        y = obj_list[-1],
        add.cell.ids = add_ids
)

# Keep workspace clean (but don’t rm random sample objects by name)
rm(obj_list, umi_counts)

### 4) SingleR annotations -----------------------------------------------------

# mid.cell.class
out_mid <- file.path(cfg$out_singler_dir, "mid.cell.class")
res_mid <- run_singler_by_sample(
        dat = dat,
        ref_seurat = sc.ref,
        ref_assay = "originalexp",
        ref_label_col = "mid.cell.class",
        out_dir = out_mid,
        query_assay = "Xenium"
)
dat$SingleR_label_midfine  <- res_mid$labels
dat$SingleR_pruned_midfine <- res_mid$pruned

# superfine.cell.class
out_super <- file.path(cfg$out_singler_dir, "superfine")
res_super <- run_singler_by_sample(
        dat = dat,
        ref_seurat = sc.ref,
        ref_assay = "originalexp",
        ref_label_col = "superfine.cell.class",
        out_dir = out_super,
        query_assay = "Xenium"
)
dat$SingleR_label_superfine  <- res_super$labels
dat$SingleR_pruned_superfine <- res_super$pruned

# Clean up big reference if you want
rm(sc.ref, res_mid, res_super)

### 5) Normalize + integrate + cluster ----------------------------------------

# SCTransform + PCA
dat <- SCTransform(dat, assay = "Xenium", verbose = FALSE)
dat <- RunPCA(dat, npcs = 100, verbose = FALSE)

# Harmony integration (Seurat v5)
dat <- IntegrateLayers(
        object = dat,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        verbose = TRUE,
        normalization.method = "SCT"
)

# Quick elbow plot if needed, this will plot PCs not harmony, but can still be helpful
ElbowPlot(dat, ndims = 70) + geom_hline(yintercept = 2.5)

# run multiple reso and dims and look at output to decide final clustering
dim_reso_test(dataset=dat, dims_list=list(1:30, 1:35, 1:45, 1:50), 
              resolutions=c(0.05, 0.1, 0.2), 
              output_dir=NA_character_, red='harmony',
              group_vars = c('SingleR_pruned_midfine','ILAE_score','orig.ident'))

# cluster will run on 45PCs and 0.05reso
dat <- Xenium_SeuratLK(
        dat,
        DIMs = 1:45,
        cluster.resolution = 0.05,
        red = "harmony",
        playMario = FALSE
)

### 6) Coarse cell type labels --------------------------------------------------
p1<- DimPlot(dat, group.by= 'seurat_clusters', shuffle = TRUE, label=T) + coord_fixed() 
p2<- DimPlot(dat, group.by = 'SingleR_pruned_midfine', shuffle = TRUE, label=T) + coord_fixed() 
wrap_plots(p1,p2, nrow=1)

# IMPORTANT:
# This mapping is dataset/version-specific. If you change clustering parameters,
# you should re-validate this mapping.
recode_vector <- c(
        "0"  = "Oligo",
        "1"  = "Astro",
        "2"  = "Micro",
        "3"  = "ExN",
        "4"  = "Vascular",
        "5"  = "OPC",
        "6"  = "InN",
        "7"  = "CSF"
)

dat$celltype <- dplyr::recode(dat$seurat_clusters, !!!recode_vector)

