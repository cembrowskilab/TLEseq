# 06_exn_CA1_spatial_analysis.R
#
# Spatial positions of deep and superficial CA1 excitatory neurons
# Using MegaData_v8
# Author: Larissa Kraus
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
# 1. Create manually defined CA1 image crops
# 2. Align deep/superficial CA1 positions across FOVs
# 3. Plot normalized positions and bootstrap median differences
# 4. Save combined figure and per-FOV summary

### LOAD PACKAGES ##############################################################
suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(purrr)
        library(readr)
        library(ggplot2)
        library(patchwork)
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


# Bootstrap 95% confidence interval for the median deep-minus-superficial difference.
boot_diff <- function(v_deep, v_super, R = 1000) {
        v_deep <- v_deep[is.finite(v_deep)]
        v_super <- v_super[is.finite(v_super)]
        if (length(v_deep) < 2 || length(v_super) < 2) return(c(ci_low = NA_real_, ci_high = NA_real_))
        b <- replicate(R, median(sample(v_deep, length(v_deep), replace = TRUE)) -
                               median(sample(v_super, length(v_super), replace = TRUE)))
        setNames(as.numeric(quantile(b, c(0.025, 0.975), na.rm = TRUE)), c("ci_low", "ci_high"))
}

# Extract cell centroids from a cropped Xenium FOV and match their metadata.
get_df2_for_fov <- function(fov, seurat_obj, depth_col, groups_of_interest) {
        img <- seurat_obj@images[[fov]]
        if (is.null(img)) stop("Missing cropped image: ", fov)
        if (!"centroids" %in% names(img@boundaries)) stop("Centroid boundaries missing in: ", fov)
        
        centroids <- img@boundaries[["centroids"]]
        coords <- as.data.frame(centroids@coords)
        if (ncol(coords) < 2) stop("Expected two centroid coordinates in: ", fov)
        names(coords)[1:2] <- c("x", "y")
        coords$cell <- centroids@cells
        if (nrow(coords) != length(coords$cell)) stop("Coordinate/cell mismatch in: ", fov)
        
        meta <- seurat_obj@meta.data
        if (!depth_col %in% colnames(meta)) stop("Missing metadata column: ", depth_col)
        matched <- match(coords$cell, rownames(meta))
        if (anyNA(matched)) warning(sum(is.na(matched)), " cropped centroids have no metadata match in ", fov)
        df <- cbind(coords, meta[matched, , drop = FALSE])
        df$fov <- fov
        df <- df[!is.na(df[[depth_col]]) & df[[depth_col]] %in% groups_of_interest &
                         is.finite(df$x) & is.finite(df$y), , drop = FALSE]
        if (!nrow(df)) return(NULL)
        df[[depth_col]] <- factor(as.character(df[[depth_col]]), levels = groups_of_interest)
        df
}

### USER-EDITABLE PATHS ########################################################
cfg <- list(
        datEx_rds = "path/to/directory/data/datEX_v8.rds",
        output_dir_data = "path/to/directory/data",
        output_dir_plots = "path/to/directory/plots"
)

plot_dir <- file.path(cfg$output_dir_plots, "ExN_v8", "deep_superCA1")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

### LOAD EXN OBJECT ############################################################
datEx <- load_if_missing("datEx", cfg$datEx_rds) else {
        message("Using datEx already in the environment.")
}

pal_exN_subtype_broad <- c("L6b" = "#E69F00", "dSub" = "#56B4E9",
                           "mSub" = "#009E73", "sCA1" = "#F0E442", "sSub" = "#0072B2",
                           "dCA1" = "#D55E00", "CA2-3" = "#CC79A7", "GC" = "#666666",
                           "dCA3" = "#007756", "MC" = "#1C91D4")

depth_col <- "broad_subtype_ExN"
groups_of_interest <- c("dCA1", "sCA1")
axis_used <- "y"  # Alignment axis; verify orientation for each manually selected crop.

if (!axis_used %in% c("x", "y")) stop("axis_used must be 'x' or 'y'.")
if (!depth_col %in% colnames(datEx@meta.data)) stop("Missing metadata column: ", depth_col)

### CROP CA1 REGIONS ###########################################################
datEx_sub <- datEx

# Original crop coordinates
datEx_sub[["B1_CA1"]] <- Crop(datEx_sub[["B1"]],
                              x = c(9900, 11200), y = c(5600, 8100), coords = "plot")

datEx_sub[["F1_CA1"]] <- Crop(datEx_sub[["F1"]],
                              x = c(7800, 9400), y = c(2400, 5100), coords = "plot")

datEx_sub[["H1_CA1"]] <- Crop(datEx_sub[["H1"]],
                              x = c(8500, 9900), y = c(4900, 7600), coords = "plot")

datEx_sub[["I1_CA1"]] <- Crop(datEx_sub[["I1"]],
                              x = c(9200, 11100), y = c(4300, 6900), coords = "plot")

datEx_sub[["J1_CA1"]] <- Crop(datEx_sub[["J1"]],
                              x = c(600, 2500), y = c(4900, 7600), coords = "plot")

datEx_sub[["L10_CA1"]] <- Crop(datEx_sub[["L10"]],
                               x = c(300, 1800), y = c(2400, 5100), coords = "plot")

datEx_sub[["L20_CA1"]] <- Crop(datEx_sub[["L20"]],
                               x = c(6500, 8100), y = c(2400, 5100), coords = "plot")

datEx_sub[["L5_CA1"]] <- Crop(datEx_sub[["L5"]],
                              x = c(8700, 10800), y = c(4900, 7600), coords = "plot")

datEx_sub[["M70_CA1"]] <- Crop(datEx_sub[["M70"]],
                               x = c(0, 2600), y = c(3700, 6300), coords = "plot")

datEx_sub[["D8_CA1"]] <- Crop(datEx_sub[["D8"]],
                              x = c(7000, 9000), y = c(3500, 6700), coords = "plot")

datEx_sub[["P28_CA1"]] <- Crop(datEx_sub[["P28"]],
                               x = c(800, 3500), y = c(1200, 6500), coords = "plot")

datEx_sub[["P60_CA1"]] <- Crop(datEx_sub[["P60"]],
                               x = c(200, 2300), y = c(4500, 8600), coords = "plot")

datEx_sub[["P71_CA1"]] <- Crop(datEx_sub[["P71"]],
                               x = c(8000, 11500), y = c(4500, 8500), coords = "plot")

datEx_sub[["P51_CA1"]] <- Crop(datEx_sub[["P51"]],
                               x = c(4500, 9000), y = c(6000, 8500), coords = "plot")

datEx_sub[["M56_CA1"]] <- Crop(datEx_sub[["M56"]],
                               x = c(300, 2300), y = c(2500, 6000), coords = "plot")

datEx_sub[["D11_CA1"]] <- Crop(datEx_sub[["D11"]],
                               x = c(1500, 5000), y = c(400, 2500), coords = "plot")
# Save cropped images and metadata for reuse.
dir.create(cfg$output_dir_data, recursive = TRUE, showWarnings = FALSE)
saveRDS(datEx_sub, file.path(cfg$output_dir_data, "datEX_v8_CA1sub.rds"))

### ALIGN DEEP AND SUPERFICIAL CA1 POSITIONS ###################################

set.seed(42)

source_fovs <- sub("_CA1$", "", grep("_CA1$", names(datEx_sub@images), value = TRUE))

fov_crops <- paste0(source_fovs, "_CA1")

combined_list <- lapply(fov_crops, get_df2_for_fov, seurat_obj = datEx_sub,
                        depth_col = depth_col, groups_of_interest = groups_of_interest)

combined_list <- Filter(Negate(is.null), combined_list)

if (!length(combined_list)) stop("No deep/superficial CA1 cells found in cropped FOVs.")

combined <- bind_rows(combined_list)

# Both groups must be present to determine orientation and calculate a difference.
fov_counts <- combined %>%
        count(fov, .data[[depth_col]], name = "n_cells") %>%
        tidyr::pivot_wider(names_from = all_of(depth_col), values_from = n_cells, values_fill = 0)

for (group in groups_of_interest) if (!group %in% names(fov_counts)) fov_counts[[group]] <- 0L

invalid_fovs <- fov_counts %>% filter(.data[["dCA1"]] < 2 | .data[["sCA1"]] < 2) %>% pull(fov)

if (length(invalid_fovs)) {
        warning("Excluding FOVs with fewer than two cells in either group: ", paste(invalid_fovs, collapse = ", "))
        combined <- combined %>% filter(!fov %in% invalid_fovs)
}

if (!nrow(combined)) stop("No FOV has at least two cells of both CA1 subtypes.")

print(fov_counts)

# Flip the selected axis within each FOV so the deep CA1 median is higher.
# This is an orientation convention, not an independent test of directionality.
aligned <- combined %>%
        group_by(fov) %>%
        mutate(raw_metric = .data[[axis_used]],
               min_m = min(raw_metric), max_m = max(raw_metric),
               med_deep_raw = median(raw_metric[.data[[depth_col]] == "dCA1"]),
               med_super_raw = median(raw_metric[.data[[depth_col]] == "sCA1"]),
               flip_needed = med_deep_raw < med_super_raw,
               metric_aligned = if_else(flip_needed, (max_m + min_m) - raw_metric, raw_metric),
               x_aligned = if (axis_used == "x") if_else(flip_needed, max(x) + min(x) - x, x) else x,
               y_aligned = if (axis_used == "y") if_else(flip_needed, max(y) + min(y) - y, y) else y,
               x_norm01 = if (diff(range(x_aligned)) > 0) (x_aligned - min(x_aligned)) / diff(range(x_aligned)) else NA_real_,
               y_norm01 = if (diff(range(y_aligned)) > 0) (y_aligned - min(y_aligned)) / diff(range(y_aligned)) else NA_real_) %>%
        ungroup()

metric_col <- if (axis_used == "y") "y_norm01" else "x_norm01"

if (anyNA(aligned[[metric_col]])) stop("Some FOVs have zero coordinate range on the selected axis.")

### NORMALIZED SPATIAL SCATTER PLOTS ###########################################

p_scatter_facets_norm01 <- ggplot(aligned,
                                  aes(x = x_norm01, y = y_norm01, color = .data[[depth_col]])) +
        geom_point(size = 0.8) +
        scale_y_reverse() +  # Display y as in image coordinates.
        scale_color_manual(values = pal_exN_subtype_broad[groups_of_interest], drop = FALSE) +
        facet_wrap(~ fov, ncol = 4) +
        labs(x = "Scaled x-coordinate", y = "Scaled y-coordinate", color = NULL) +
        theme_classic(base_size = 10) + theme(legend.position = "right")

print(p_scatter_facets_norm01)

### PER-FOV MEDIAN DIFFERENCES AND BOOTSTRAP CIS ###############################

per_fov_eff <- aligned %>%
        group_by(fov) %>%
        summarise(flipped = first(flip_needed),
                  n_deep = sum(.data[[depth_col]] == "dCA1"),
                  n_super = sum(.data[[depth_col]] == "sCA1"),
                  med_deep = median(.data[[metric_col]][.data[[depth_col]] == "dCA1"]),
                  med_super = median(.data[[metric_col]][.data[[depth_col]] == "sCA1"]),
                  diff_med = med_deep - med_super,
                  ci = list(boot_diff(.data[[metric_col]][.data[[depth_col]] == "dCA1"],
                                      .data[[metric_col]][.data[[depth_col]] == "sCA1"], R = 1000)),
                  .groups = "drop") %>%
        unnest_wider(ci)

print(per_fov_eff, n = Inf, width = Inf)

write_csv(per_fov_eff, file.path(plot_dir, "CA1_depth_per_FOV_bootstrap.csv"))

p_forest <- ggplot(per_fov_eff, aes(y = reorder(fov, diff_med), x = diff_med)) +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_point() +
        geom_segment(aes(x = ci_low, xend = ci_high, y = reorder(fov, diff_med),
                         yend = reorder(fov, diff_med))) +
        labs(x = "Median difference (dCA1 - sCA1)", y = "FOV") +
        theme_classic(base_size = 12)

print(p_forest)

### POOLED VIOLIN PLOT #########################################################

# Rebase each FOV to its own superficial CA1 median before pooling cells.
aligned_rebased <- aligned %>%
        group_by(fov) %>%
        mutate(super_med = median(.data[[metric_col]][.data[[depth_col]] == "sCA1"]),
               metric_rebased = .data[[metric_col]] - super_med) %>%
        ungroup()

p_violin_rebased <- ggplot(aligned_rebased,
                           aes(x = .data[[depth_col]], y = metric_rebased, fill = .data[[depth_col]])) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_violin(trim = FALSE, alpha = 0.85) +
        geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
        scale_fill_manual(values = pal_exN_subtype_broad[groups_of_interest], drop = FALSE) +
        labs(x = NULL, y = "Normalized depth relative to sCA1 median") +
        theme_classic(base_size = 12) + theme(legend.position = "none")
print(p_violin_rebased)

### SAVE FIGURES ###############################################################
combined_plot <- (p_scatter_facets_norm01 / (p_violin_rebased | p_forest)) +
        plot_layout(heights = c(2.2, 1))
print(combined_plot)
ggsave(file.path(plot_dir, "CA1_deep_super_spatial_combined.pdf"), combined_plot,
       width = 15, height = 13)

### CLEAN UP ###################################################################
# Keep datEx, datEx_sub, aligned, aligned_rebased, per_fov_eff, and plots available.
rm(combined_list, combined, fov_counts, invalid_fovs, source_fovs)
