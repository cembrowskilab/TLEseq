################################################################################
#
# 05_megadata_final_annotations.R
#
# This script brings annotations of subtypes back into megadata to have one final 
# object which can be uploaded to GEO as 'processed data' and used more easily
#
# Author: Larissa Kraus (October 2026) 
######

### 0) Setup -------------------------------------------------------------------
uppressPackageStartupMessages({
  library(conflicted)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(future)
  library(rlang)
  library(beepr)       # optional: mario sound
  
})


### PACKAGE CONFLICTS ---------------------------------------------------------
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

### 2) User-editable paths and load inputs -------------------------------------
# Fill in paths before running; existing objects are reused where applicable.

cfg <- list(xenium_rds_dir = NA_character_, 
            metadata_csv = NA_character_,
            exn_rds = NA_character_, 
            inn_rds = NA_character_, 
            micro_rds = NA_character_,
            domain_annotations_rds = NA_character_)


### 3) Include subtype labels from subset analysis------------------------------
# after running all scripts in this Github, we have cell type annotations from subset analysis
# these are now added to the main data for easy plotting and access to our annotations for others

#For this script load metadata, domain annotations, final excitatory neuron dataset, inhibitory neuron dataset and microglia dataset
meta <- metadata
datEx <- load_if_missing("datEx", cfg$exn_rds)
datIN <- load_if_missing("datIN", cfg$inn_rds)
datM <- load_if_missing("datM", cfg$micro_rds)
dat_md <- load_if_missing("dat_md", cfg$domain_annotations_rds)

#add domain annotations
dat$region_md<-dat_md$region_md

#transfer following annotations to dat
# datEx$broad_subtype_ExN_LK
# datEx$fine_subtype_ExN_LK
# datIN$broad_subtype_InhN
# datIN$fine_subtype_InhN
# datM$cellstates

# start by copying the general celltype for everything
dat$cellsubtype <- as.character(dat$celltype)
dat$cellsubtype_fine <- as.character(dat$celltype)

# ExN annotations
cells_ex <- intersect(colnames(dat), colnames(datEx))

dat$cellsubtype[cells_ex] <- paste(
  dat$celltype[cells_ex],
  datEx$broad_subtype_ExN_LK[cells_ex],
  sep = "_"
)

dat$cellsubtype_fine[cells_ex] <- paste(
  dat$celltype[cells_ex],
  datEx$fine_subtype_ExN_LK[cells_ex],
  sep = "_"
)

# InN annotations
cells_in <- intersect(colnames(dat), colnames(datIN))

dat$cellsubtype[cells_in] <- paste(
  dat$celltype[cells_in],
  datIN$broad_subtype_InhN[cells_in],
  sep = "_"
)

dat$cellsubtype_fine[cells_in] <- paste(
  dat$celltype[cells_in],
  datIN$fine_subtype_InhN[cells_in],
  sep = "_"
)

# Microglia annotations
cells_m <- intersect(colnames(dat), colnames(datM))

dat$cellsubtype[cells_m] <- paste(dat$celltype[cells_m],
                                  datM$cellstates[cells_m],
                                  sep = "_")

# for microglia there is no separate fine annotation given,
# so use the same annotation for fine as well
dat$cellsubtype_fine[cells_m] <- paste(dat$celltype[cells_m],
                                       datM$cellstates[cells_m],
                                       sep = "_")