################################################################################
#
# 05_visium_domain_SingleR_transfer.R
#
# 1. Transfer reference Visium domain labels to Xenium domain labels
# 2. Compare transferred Visium domain labels with annotated Xenium regions.
#
# Original analysis: Larissa Kraus
# Code cleanup and standardization: Larissa Kraus (October 2026)
######

### 0) Setup -------------------------------------------------------------------
suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(ggplot2)
        library(SingleR)
        library(ExperimentHub)
})
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)
options(future.globals.maxSize = 100000 * 1024^2)

### 1) Helper functions --------------------------------------------------------
# Reuse an existing object in the global environment; otherwise load its RDS file.
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
cfg <- list(domain_rds = NA_character_, #should be "07_seumerged.rds" saved from '01_domainIdentification.R'
            spatial_reference_rds = NA_character_,
            results_dir = NA_character_, 
            output_rds = NA_character_)
if (anyNA(unlist(cfg))) stop("Set all paths in cfg")
dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cfg$output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(cfg$spatial_reference_rds), recursive = TRUE, showWarnings = FALSE)

### 3) Prepare reference and transfer labels ----------------------------------
###LOAD THOMPSON DATA (for more info on install, see '00_download_snRNAseq_reference.R')
##CONNECT TO DATAHUB, LOAD THOMPSON DATA
ehub <- ExperimentHub::ExperimentHub()

## Load the datasets of the package
myfiles <- ExperimentHub::query(ehub, "humanHippocampus2024")

#Download spatial data
spatial_hpc_spe <- myfiles[["EH9605"]]

#convert file to seurat object
spatial.ref <- as.Seurat(spatial_hpc_spe, counts = "counts", data = "logcounts")

saveRDS(spatial.ref, cfg$spatial_reference_rds)


#possible meta data for grouping:
unique(spatial.ref$domain)
# Levels: GCL CA2.4 CA1 SUB SUB.RHP RHP GABA SL.SR ML SR.SLM SLM.SGZ WM.1 WM.2 WM.3 Vascular Choroid

### 4) SINGLER - RUN USING domain ----------------------------------------------
# Set the output directory
output_dir <- cfg$results_dir
if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
}

# Prepare the reference expression matrix
#Load Thompson et al data
spatial.ref <- load_if_missing("spatial.ref", cfg$spatial_reference_rds, "Seurat")

# Log-Normalize spatial.ref
spatial.ref <- NormalizeData(spatial.ref, assay = "originalexp", verbose = FALSE)

#Check out some marker genes for domains to see if makes sense
domain_marker<- FindAllMarkers(spatial.ref, group.by = 'domain')
RNA.marker <- subset(domain_marker, domain_marker$pct.2<0.3) #!is very important! basically means that a marker gene is not existing in other clusters
RNA.marker <- RNA.marker %>%
        group_by(cluster) %>%
        slice_max(n = 10, order_by = pct.1) #also to order so that the most significant marker for specific cluster are on top

# Pseudobulk the spatial reference
pseudobulk <- AggregateExpression(
        spatial.ref,
        group.by = "domain",
        assays = "originalexp",
        slot = "data", 
        return.seurat = F)

# Matrix: genes x cell types
ref_expr <- pseudobulk$originalexp

###Load atlas data provided by mathias
domain.dat <- load_if_missing("domain.dat", cfg$domain_rds, "Seurat")

# Get the unique sample names of xenium data
sample_names <- unique(domain.dat$orig.ident)

# Initialize empty vectors to hold all labels
all_labels <- rep(NA, ncol(domain.dat))
all_pruned <- rep(NA, ncol(domain.dat))
names(all_labels) <- colnames(domain.dat)
names(all_pruned) <- colnames(domain.dat)

#run loop to run SingleR on each sample in domain.dat
for (sample_name in sample_names) {
        message(paste0("Running SingleR for sample: ", sample_name))
        
        # Subset Seurat object for current sample
        sample_obj <- subset(domain.dat, subset = orig.ident == sample_name)
        
        # Extract log-normalized expression
        query_expr <- GetAssayData(sample_obj, layer = "data", assay = "RNA")
        
        # Subset to common genes
        common_genes <- intersect(rownames(ref_expr), rownames(query_expr))
        if (length(common_genes) < 10) stop("Too few common genes for ", sample_name)
        ref_sub <- ref_expr[common_genes, , drop = FALSE]
        query_sub <- query_expr[common_genes, , drop = FALSE]
        
        # Run SingleR
        singleR_results <- SingleR(
                test = as.matrix(query_sub),
                ref = as.matrix(ref_sub),
                labels = colnames(ref_expr)
        )
        
        # Save results to the specified output directory
        saveRDS(singleR_results, file = file.path(output_dir, paste0("SingleR_", sample_name, "_domain_results.rds")))
        
        # Add results to the full metadata using cell barcodes
        sample_cells <- colnames(sample_obj)
        all_labels[sample_cells] <- singleR_results$labels
        all_pruned[sample_cells] <- singleR_results$pruned.labels
}

# Add to metadata of the full domain.dat object
domain.dat$SingleR_label_domain <- all_labels
domain.dat$SingleR_pruned_domain <- all_pruned

#make new metadata with simplified labels more similar to regions_md
recode_vector <- c(
        "Vascular" = "Vascular",
        "Choroid"  = "Choroid",
        "SR.SLM"   = "SL.SR.SLM",
        "SL.SR"    = "SL.SR.SLM",
        "SLM.SGZ"  = "SL.SR.SLM",
        "WM.1"     = "WM",
        "WM.2"     = "WM",
        "WM.3"     = "WM",
        "CA2.4"    = "CA2.4",
        "SUB.RHP"  = "SUB",
        "RHP"      = "SUB",
        "SUB"      = "SUB",
        "CA1"      = "CA1",
        "GCL"      = "GCL",
        "ML"       = "ML",
        "GABA"     = "GABA")

domain.dat$SingleR_pruned_domain <- recode(domain.dat$SingleR_pruned_domain, !!!recode_vector)


#make new metadata with more simplified labels for plotting
recode_vector <- c(
        "Vascular" = "Vasc.Chor",
        "Choroid"  = "Vasc.Chor",
        "SR.SLM"   = "SL.SR.SLM",
        "SL.SR"    = "SL.SR.SLM",
        "SLM.SGZ"  = "SL.SR.SLM",
        "WM.1"     = "WM",
        "WM.2"     = "WM",
        "WM.3"     = "WM",
        "CA2.4"    = "CA2.4",
        "SUB.RHP"  = "SUB",
        "RHP"      = "SUB",
        "SUB"      = "SUB",
        "CA1"      = "CA1",
        "GCL"      = "GCL",
        "ML"       = "ML",
        "GABA"     = "GABA")

domain.dat$SingleR_pruned_domain_simple <- recode(domain.dat$SingleR_pruned_domain, !!!recode_vector)


#save
saveRDS(domain.dat, cfg$output_rds)



### 3) Domain label comparisons and plotting-------------------------------------
# get data for plotting
domain.dat <- load_if_missing("domain.dat", cfg$domain_rds, "Seurat")
md <- domain.dat@meta.data

# make sure labels are character (or factor) and change NA labels for both region_md and singleR labels
md$SingleR_pruned_domain <- as.character(md$SingleR_pruned_domain)
md$SingleR_pruned_domain[is.na(md$SingleR_pruned_domain)] <- "Unassigned"
md$region_md <- as.character(md$region_md)
md$region_md[is.na(md$region_md)] <- "Unassigned"
md$SingleR_pruned_domain_simple <- as.character(md$SingleR_pruned_domain_simple)
md$SingleR_pruned_domain_simple[is.na(md$SingleR_pruned_domain_simple)] <- "Unassigned"

md$region_md <- factor(md$region_md,
                       levels = c("GCL","ML","CA2_CA4","CA1","SUB","GABA",
                                  "SL.SR.SLM","WM","pia","vascular","void","Unassigned"))

md$SingleR_pruned_domain <- factor(md$SingleR_pruned_domain,
                                      levels = c("GCL","ML","CA2.4","CA1","SUB","GABA",
                                                 "SL.SR.SLM","WM",
                                                 "Vascular","Choroid","Unassigned"),
                                      exclude = NULL)

md$SingleR_pruned_domain_simple <- factor(md$SingleR_pruned_domain_simple,
                                             levels = c("GCL","ML","CA2.4","CA1","SUB","GABA",
                                                        "SL.SR.SLM","WM",
                                                        "Vasc.Chor","Unassigned"),
                                             exclude = NULL)

# get label levels in the order ggplot will see them
labs <- unique(md$SingleR_pruned_domain)
labs_simple <- unique(md$SingleR_pruned_domain_simple)

# make a distinct palette with the right length
cols <- Polychrome::createPalette(length(labs),
                                  seedcolors = c("#000000", "#E69F00", "#56B4E9"))
cols_simple <- Polychrome::createPalette(length(labs_simple),
                                         seedcolors = c("#000000", "#E69F00", "#56B4E9"))

# IMPORTANT: name the colors by label
cols <- setNames(cols, labs)
cols_simple <- setNames(cols_simple, labs_simple)

#plot proportions for simple groups
print(ggplot(md, aes(x = region_md, fill = SingleR_pruned_domain_simple)) +
        geom_bar(position = "fill") +
        facet_wrap(~ ILAE_score, nrow=1) +
        scale_fill_manual(values = cols_simple, drop = FALSE) +
        theme_classic() +
        ylab("Proportion") +
        xlab("Regions_md") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))

##order with highest prop at bottom
# calculate overall proportions
prop_order <- md %>%
  count(SingleR_pruned_domain_simple) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(prop) %>%
  pull(SingleR_pruned_domain_simple)

# reorder factor so largest proportion is first (bottom of stack)
md_ordered<-md
md_ordered$SingleR_pruned_domain_simple <- factor(md_ordered$SingleR_pruned_domain_simple,
                                                     levels = prop_order)

# plot
print(ggplot(md_ordered, aes(x = region_md, fill = SingleR_pruned_domain_simple)) +
        geom_bar(position = "fill") +
        facet_wrap(~ ILAE_score, nrow=1) +
        scale_fill_manual(values = cols_simple, drop = FALSE) +
        theme_classic() +
        ylab("Proportion") +
        xlab("Region") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))


## make a table to see proportions and counts
tab_long <- md %>%
  dplyr::count(region_md, SingleR_pruned_domain, name = "n") %>%
  group_by(region_md) %>%
  mutate(prop = (n / sum(n))*100) %>%
  ungroup() %>%
  arrange(region_md, desc(prop))




