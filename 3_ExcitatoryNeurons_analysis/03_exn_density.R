# 03_exn_density.R
#
# Xenium analysis of excitatory neurons in human hippocampus samples
# Using MegaData_v8
# Author: Larissa Kraus
# (October 2026)
#
# This script:
# 1. Calculates sample surface areas from full MegaData FOV convex hulls
# 2. Calculates broad ExN subtype densities (cells/mm^2) per sample
# 3. Compares ILAE groups using Kruskal-Wallis and BH-adjusted Dunn tests
# 4. Prints statistical results and generates density plots without p-value annotations

### LOAD PACKAGES ##############################################################
suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(readr)
        library(purrr)
        library(ggplot2)
        library(patchwork)
        library(Seurat)
        library(FSA)
        library(scales)
        library(conflicted)
})

### PACKAGE CONFLICTS
conflicted::conflict_prefer("count", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")
conflicted::conflict_prefer("lag", "dplyr")
conflicted::conflict_prefer("first", "dplyr")

options(future.globals.maxSize = 100000 * 1024^2)

### HELPER FUNCTIONS ###########################################################
# Reuse an existing object; otherwise load it from a user-configured RDS path.
# Returns the object without overwriting an existing object in the calling environment.
load_if_missing <- function(object_name, path, envir = .GlobalEnv) {
        if (exists(object_name, envir = envir, inherits = FALSE)) {
                message("Reusing ", object_name, " from the environment.")
                return(get(object_name, envir = envir, inherits = FALSE))
        }
        if (is.null(path) || !nzchar(path) || !file.exists(path)) {
                stop("Missing ", object_name, ". Set its input RDS path in cfg: ", path)
        }
        message("Loading ", object_name, " from: ", path)
        readRDS(path)
}


# Compute polygon area using the shoelace formula (coordinates in um).
polygon_area <- function(x, y) {
        x <- c(x, x[1])
        y <- c(y, y[1])
        abs(sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1])) / 2
}

# Identify spatial x/y coordinate columns returned by Seurat.
find_xy_columns <- function(coords) {
        possible_pairs <- list(c("x", "y"), c("imagecol", "imagerow"),
                               c("x.centroid", "y.centroid"))
        pair_found <- purrr::keep(possible_pairs, ~ all(.x %in% colnames(coords)))
        if (length(pair_found) == 0) {
                stop("Could not identify coordinate columns. Available: ",
                     paste(colnames(coords), collapse = ", "))
        }
        pair_found[[1]]
}

# Convex hull of all tissue-cell coordinates in one full-dataset FOV.
# Coordinate units are um; hull_area is converted to mm^2.
calculate_hull_area <- function(fov_name, seurat_obj) {
        coords <- GetTissueCoordinates(seurat_obj[[fov_name]]) %>% as.data.frame()
        xy_cols <- find_xy_columns(coords)
        coords <- coords %>%
                transmute(x = .data[[xy_cols[1]]], y = .data[[xy_cols[2]]]) %>%
                filter(is.finite(x), is.finite(y)) %>%
                distinct(x, y)
        
        if (nrow(coords) < 3) {
                warning("Fewer than three coordinates for ", fov_name)
                return(tibble(orig.ident = fov_name, n_coordinates = nrow(coords),
                              hull_area_um2 = NA_real_, hull_area = NA_real_))
        }
        
        hull_coords <- coords[chull(coords$x, coords$y), ]
        area_um2 <- polygon_area(hull_coords$x, hull_coords$y)
        tibble(orig.ident = fov_name, n_coordinates = nrow(coords),
               hull_area_um2 = area_um2, hull_area = area_um2 / 1e6)
}

calculate_subtype_stats <- function(dat_sub, subtype_name, ilae_levels) {
        dat_sub <- dat_sub %>%
                filter(is.finite(density), !is.na(ILAE_score)) %>%
                droplevels()
        
        group_summary <- dat_sub %>%
                group_by(ILAE_score) %>%
                summarise(n = n(),median = median(density), mean = mean(density),
                          .groups = "drop") %>%
                complete(ILAE_score = factor(ilae_levels, levels = ilae_levels)) %>%
                pivot_wider(names_from = ILAE_score, values_from = c(n, median, mean),names_sep = ".")
        
        n_groups <- n_distinct(dat_sub$ILAE_score)
        
        kw_result <- if (n_groups >= 2) {
                tryCatch(kruskal.test(density ~ ILAE_score,data = dat_sub),
                         error = function(e) NULL)
        } else {
                NULL
        }
        
        if (!is.null(kw_result)) {
                kw_table <- tibble(
                        subtype = subtype_name,
                        n_total = nrow(dat_sub),
                        df = unname(kw_result$parameter),
                        statistic = unname(kw_result$statistic),
                        p.value = kw_result$p.value,
                        method = kw_result$method
                )
        } else {
                kw_table <- tibble(
                        subtype = subtype_name,
                        n_total = nrow(dat_sub),
                        df = NA_real_,
                        statistic = NA_real_,
                        p.value = NA_real_,
                        method = "Insufficient groups"
                )
        }
        
        kw_table <- bind_cols(kw_table, group_summary)
        
        dunn_result <- if (n_groups >= 2) {
                tryCatch(FSA::dunnTest(density ~ ILAE_score,
                                       data = dat_sub,
                                       method = "bh")$res,
                         error = function(e) NULL)
        } else {
                NULL
        }
        
        if (is.null(dunn_result) || nrow(dunn_result) == 0) {
                dunn_table <- tibble(
                        subtype = character(),
                        group1 = character(),
                        group2 = character(),
                        z = numeric(),
                        p.unadj = numeric(),
                        p.adj = numeric()
                )
        } else {
                dunn_table <- dunn_result %>%
                        rename(z = Z, p.unadj = P.unadj,p.adj = P.adj) %>%
                        separate(Comparison,
                                 into = c("group1", "group2"),
                                 sep = "\\s*-\\s*") %>%
                        mutate(subtype = subtype_name,
                               group1 = trimws(group1),
                               group2 = trimws(group2)) %>%
                        select(subtype, group1, group2, z,
                               p.unadj, p.adj)
        }
        
        list(kw = kw_table, dunn = dunn_table)
}

# Plot cell density by ILAE group without statistical annotations.
build_density_plot <- function(subtype_name, cell_props, pal_ILAE_score) {
        dat_sub <- cell_props %>%
                filter(subtype == subtype_name, is.finite(density), !is.na(ILAE_score))
        
        ggplot(dat_sub, aes(x = ILAE_score, y = density, fill = ILAE_score)) +
                geom_boxplot(outlier.shape = NA, alpha = 0.4) +
                geom_jitter(width = 0.2, alpha = 0.8, size = 1) +
                scale_fill_manual(values = pal_ILAE_score) +
                scale_y_continuous(limits = c(0, NA),
                                   labels = scales::label_number(accuracy = 1),
                                   expand = expansion(mult = c(0, 0.05))) +
                labs(title = subtype_name, x = NULL, y = expression("Cells per mm"^2)) +
                theme_minimal() +
                theme(legend.position = "none", panel.grid = element_blank(),
                      axis.line = element_line(colour = "black", linewidth = 0.4),
                      axis.ticks = element_line(colour = "black"))
}

### USER-EDITABLE PATHS ########################################################
# Match datEx_rds and output_dir_plots to Scripts 01 and 02.
cfg <- list(
        datEx_rds = "path/to/directory/data/datEX_v8.rds",
        megadata_rds = "path/to/directory/MegaData_v8_annotated.rds",
        output_dir_data = "path/to/directory/data",
        output_dir_plots = "path/to/directory/plots"
)

### LOAD EXN OBJECT ############################################################
# Reuse the object from Script 01/02 if already loaded.
datEx <- load_if_missing("datEx", cfg$datEx_rds) else {
        message("Using datEx already in the environment.")
}

### COLOR PALETTE ##############################################################
pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C", "1" = "#7BAD00",
                    "2" = "#00C0C5", "3" = "#E69F00")

### CALCULATE OR LOAD SAMPLE SURFACE AREAS ####################################
# The denominator must represent sampled tissue area, not the spatial extent
# of ExN alone. Therefore calculate convex hulls from the full dat object.
# If you have done this before, best reused to avoid loading the large full dataset again.
area_file <- file.path(cfg$output_dir_data, "sample_surface_areas_v8.csv")
dir.create(cfg$output_dir_data, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(area_file)) {
        dat <- load_if_missing("dat", cfg$megadata_rds)
        
        sample_surface_areas_v8 <- purrr::map_dfr(names(dat@images),
                                                  calculate_hull_area,
                                                  seurat_obj = dat) %>%
                arrange(orig.ident)
        
        print(sample_surface_areas_v8)
        print(summary(sample_surface_areas_v8$hull_area))
        readr::write_csv(sample_surface_areas_v8, area_file)
        message("Saved sample areas to: ", area_file)
} else {
        sample_surface_areas_v8 <- readr::read_csv(area_file, show_col_types = FALSE)
        message("Loaded existing sample areas from: ", area_file)
}

### EXCITATORY-NEURON DENSITY BY ILAE GROUP ###################################
subtype_col <- "broad_subtype_ExN"
ilae_levels <- c("P", "0", "1", "2", "3")
stats_dir <- file.path(cfg$output_dir_plots, "ExN_v8", "stats_density")
plot_dir <- file.path(cfg$output_dir_plots, "ExN_v8")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# hull_area is in mm^2 (converted from um^2 above).
area_table <- sample_surface_areas_v8 %>% select(orig.ident, hull_area)
if (anyDuplicated(area_table$orig.ident)) stop("Duplicate sample names in surface area table.")

df <- datEx@meta.data %>%
        as.data.frame() %>%
        filter(!is.na(orig.ident),
               !is.na(ILAE_score),
               !is.na(.data[[subtype_col]]))

sample_metadata <- df %>% distinct(orig.ident, ILAE_score)

if (anyDuplicated(sample_metadata$orig.ident)) {
        stop("Some samples have more than one ILAE_score.")
}

subtype_order <- c("GC", "MC", "CA2-3", "dCA3", "sCA1", "dCA1",
                   "sSub", "mSub", "dSub", "L6b")
subtype_order <- subtype_order[subtype_order %in% unique(as.character(df[[subtype_col]]))]

cell_counts <- df %>%
        count(orig.ident, subtype = .data[[subtype_col]],
              name = "celltype_count") %>%
        complete(orig.ident = sample_metadata$orig.ident,
                 subtype = subtype_order,
                 fill = list(celltype_count = 0))

cell_props <- cell_counts %>%
        left_join(sample_metadata, by = "orig.ident") %>%
        left_join(area_table, by = "orig.ident") %>%
        mutate(subtype = factor(subtype, levels = subtype_order),
               ILAE_score = factor(ILAE_score, levels = ilae_levels),
               density = celltype_count / hull_area)

missing_areas <- cell_props %>%
        filter(is.na(hull_area)) %>%
        distinct(orig.ident) %>%
        pull(orig.ident)

if (length(missing_areas) > 0) {
        stop("Missing sample areas for: ",
             paste(missing_areas, collapse = ", "))
}

if (any(cell_props$hull_area <= 0, na.rm = TRUE)) {
        stop("All hull areas must be greater than zero.")
}

### CALCULATE AND SAVE STATISTICS #############################################
stats_by_subtype <- map(subtype_order, function(subtype_name) {
        dat_sub <- cell_props %>%
                filter(subtype == subtype_name)
        
        calculate_subtype_stats(dat_sub, subtype_name, ilae_levels)
})

names(stats_by_subtype) <- subtype_order

kw_table <- map_dfr(stats_by_subtype, "kw") %>%
        relocate(subtype, n_total, df, statistic, p.value, method)

dunn_table <- map_dfr(stats_by_subtype, "dunn")
# Print all test results to the console (no p-values or stars on plots).
print(kw_table, n = Inf, width = Inf)
print(as_tibble(dunn_table), n = Inf, width = Inf)

### GENERATE PLOTS ###################################################
plots <- map(subtype_order, ~ build_density_plot(.x, cell_props, pal_ILAE_score)) %>%
        set_names(subtype_order)

combined_plot <- wrap_plots(plots[subtype_order], nrow = 2)

print(combined_plot)

### CLEAN UP ###################################################################
# Keep datEx and the statistical output tables available for inspection.
rm(area_table, df, sample_metadata, cell_counts, missing_areas, stats_by_subtype)



