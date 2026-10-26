
################################################################################
#
# 05_microglia_spatial_expression.R
#
# Xenium analysis of microglia in human hippocampus samples
# Uses datM_clean_v8_annotated.rds generated in Script 01_microglia
#
# Analyses:
# 1. Spatial distribution of microglia cell states
# 2. Spatial expression of CD83 and CCL4
# 3. Regional expression of CD83 and CCL4 across ILAE groups
#
# Original analysis: Rachel J. Boyd (2026)
# Code cleanup and standardization: Larissa Kraus (October 2026)
#
### 0) Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
        library(conflicted)
        library(Seurat)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(patchwork)
        library(ggrastr)
        library(cowplot)
})

conflicted::conflicts_prefer(dplyr::filter, 
                             dplyr::select, 
                             dplyr::rename,
                             dplyr::lag, 
                             dplyr::count, 
                             dplyr::first)

### 1) Helper functions --------------------------------------------------------

# Reuse an existing object in the global environment; otherwise load its RDS file
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

# Build spatial cell-state maps with ILAE headers and a shared legend
fov_figure <- function(groups, n_col, show_label = FALSE, body_h = 10) {
        
        fill <- scale_fill_manual(values = pal_cluster_named[keep_states],
                                  limits = keep_states, na.value = NA, name = NULL)
        
        # Plot one FOV and rasterize the cell layer to reduce figure size
        panel <- function(fov) {
                p <- ImageDimPlot(FOV_plots, fov = fov, group.by = "cellstates",
                                  size = 0.5, border.size = NA, axes = FALSE,
                                  dark.background = FALSE) +
                        fill +
                        theme(legend.position = "none",
                              plot.margin = margin(2, 2, 2, 2))
                
                p$layers <- lapply(p$layers, ggrastr::rasterise, dpi = 300)
                
                if (show_label)
                        p <- p + annotate("text", x = -Inf, y = Inf, label = fov,
                                          hjust = -0.15, vjust = 1.5,
                                          size = 3.2, fontface = "bold")
                p
        }
        
        # Extract one shared cell-state legend
        legend <- cowplot::get_legend(
                ImageDimPlot(FOV_plots, fov = unlist(groups)[1],
                             group.by = "cellstates", size = 1.2,
                             border.size = NA, axes = FALSE,
                             dark.background = FALSE) +
                        fill + theme(legend.position = "right"))
        
        # Create one panel block per ILAE group
        blocks <- lapply(names(groups), function(score) {
                fovs <- groups[[score]]
                
                title <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
                        geom_tile(fill = pal_ILAE_score[score], width = 1, height = 0.5) +
                        annotate("text", x = 1, y = 1, label = group_titles[score],
                                 fontface = "bold", size = 4.5) +
                        theme_void() +
                        theme(plot.margin = margin(2, 2, 2, 2))
                
                title / wrap_plots(lapply(fovs, panel), ncol = n_col[score]) +
                        plot_layout(heights = c(1, body_h))
        })
        
        # Combine all ILAE blocks with one shared legend
        (wrap_plots(blocks, nrow = 1) | wrap_elements(legend)) +
                plot_layout(widths = c(length(blocks), 0.3))
}

# Make Imagefeatureplots where highest expressing cells are plotted on top. Made
# from Seurat function ImageFeaturePlot() implenting sorting as can be done for
# FeaturePlot()
ImageFeaturePlot_ordered <- function(object,features,fov = NULL,boundaries = NULL,
                                     size = 0.5,min.cutoff = NA, max.cutoff = NA,split.by = NULL,
                                     molecules = NULL,mols.size = 0.1,mols.cols = NULL,nmols = 1000,
                                     alpha = 1,border.color = "white",border.size = NULL,dark.background = TRUE,
                                     blend = FALSE,blend.threshold = 0.5,crop = FALSE,cells = NULL,
                                     scale = c("feature", "all", "none"),overlap = FALSE,axes = FALSE,
                                     combine = TRUE,coord.fixed = TRUE,
                                     order = FALSE,                   # NEW: turn ordering on/off
                                     order.by = c("value", "abs"),    # NEW: order by value or abs(value)
                                     cols = if (isTRUE(x = blend)) {
                                             c("lightgrey", "#ff0000", "#00ff00")
                                     } else {c("lightgrey", "firebrick1")}) {
        cells <- cells %||% Cells(x = object)
        scale <- scale[[1L]]
        scale <- match.arg(arg = scale)
        order.by <- match.arg(order.by)
        
        no.right <- theme(
                axis.line.y.right = element_blank(),
                axis.ticks.y.right = element_blank(),
                axis.text.y.right = element_blank(),
                axis.title.y.right = element_text(
                        face = "bold",
                        size = 14,
                        margin = margin(r = 7)
                )
        )
        
        fov <- fov %||% DefaultFOV(object = object)
        fov <- Filter(
                f = function(x) {
                        return(x %in% Images(object = object) && inherits(x = object[[x]], what = "FOV"))
                },
                x = fov
        )
        if (!length(x = fov)) {
                stop("No compatible spatial coordinates present")
        }
        
        boundaries <- boundaries %||% sapply(
                X = fov,
                FUN = function(x) {
                        return(DefaultBoundary(object = object[[x]]))
                },
                simplify = FALSE,
                USE.NAMES = TRUE
        )
        boundaries <- .BoundariesByImage(object = object, fov = fov, boundaries = boundaries)
        fov <- names(x = boundaries)
        
        if (isTRUE(x = blend) || !is.null(x = split.by)) {
                type <- ifelse(test = isTRUE(x = "blend"), yes = "Blended", no = "Split")
                if (length(x = fov) != 1L) {
                        fov <- fov[1L]
                        warning(
                                type,
                                " image feature plots can only be done on a single image, using \"",
                                fov, "\"",
                                call. = FALSE,
                                immediate. = TRUE
                        )
                }
                if (any(!overlap) && length(x = boundaries[[fov]]) > 1L) {
                        warning(
                                type,
                                " image feature plots require overlapped segmentations",
                                call. = FALSE,
                                immediate. = TRUE
                        )
                }
                overlap <- TRUE
        }
        overlap <- rep_len(x = overlap, length.out = length(x = fov))
        crop <- rep_len(x = crop, length.out = length(x = fov))
        names(x = crop) <- names(x = overlap) <- fov
        
        if (isTRUE(x = blend)) {
                if (length(x = features) != 2L) {
                        stop("Blended feature plots only works with two features")
                }
                default.colors <- eval(expr = formals(fun = ImageFeaturePlot)$cols)
                cols <- switch(
                        EXPR = as.character(x = length(x = cols)),
                        `0` = {
                                warning("No colors provided, using default colors", immediate. = TRUE)
                                default.colors
                        },
                        `1` = {
                                warning(
                                        "Only one color provided, assuming specified is double-negative and augmenting with default colors",
                                        immediate. = TRUE
                                )
                                c(cols, default.colors[2:3])
                        },
                        `2` = {
                                warning(
                                        "Only two colors provided, assuming specified are for features and augmenting with '",
                                        default.colors[1], "' for double-negatives",
                                        immediate. = TRUE
                                )
                                c(default.colors[1], cols)
                        },
                        `3` = cols,
                        {
                                warning("More than three colors provided, using only first three", immediate. = TRUE)
                                cols[1:3]
                        }
                )
        }
        
        md <- FetchData(object = object, vars = c(features, split.by[1L]), cells = cells)
        # Keep a copy of the original values (can have negatives)
        md.orig <- md
        
        split.by <- intersect(x = split.by, y = colnames(x = md))
        if (!length(x = split.by)) {
                split.by <- NULL
        }
        
        imax <- ifelse(
                test = is.null(x = split.by),
                yes = ncol(x = md),
                no = ncol(x = md) - length(x = split.by)
        )
        features <- colnames(x = md)[1:imax]
        
        min.cutoff <- mapply(
                FUN = function(cutoff, feature) {
                        return(ifelse(
                                test = is.na(x = cutoff),
                                yes  = min(md[[feature]]),
                                no   = cutoff
                        ))
                },
                cutoff = min.cutoff,
                feature = features
        )
        max.cutoff <- mapply(
                FUN = function(cutoff, feature) {
                        return(ifelse(
                                test = is.na(x = cutoff),
                                yes  = max(md[[feature]]),
                                no   = cutoff
                        ))
                },
                cutoff = max.cutoff,
                feature = features
        )
        
        check.lengths <- unique(x = vapply(
                X = list(features, min.cutoff, max.cutoff),
                FUN = length,
                FUN.VALUE = numeric(length = 1)
        ))
        if (length(x = check.lengths) != 1) {
                stop("There must be the same number of minimum and maximum cuttoffs as there are features")
        }
        
        brewer.gran <- ifelse(
                test = length(x = cols) == 1,
                yes  = brewer.pal.info[cols, ]$maxcolors,
                no   = length(x = cols)
        )
        for (i in seq_along(along.with = features)) {
                f <- features[[i]]
                data.feature <- md[[f]]
                min.use <- SetQuantile(cutoff = min.cutoff[i], data = data.feature)
                max.use <- SetQuantile(cutoff = max.cutoff[i], data = data.feature)
                data.feature[data.feature < min.use] <- min.use
                data.feature[data.feature > max.use] <- max.use
                if (brewer.gran != 2) {
                        data.feature <- if (all(data.feature == 0)) {
                                rep_len(x = 0, length.out = length(x = data.feature))
                        } else {
                                as.numeric(x = as.factor(x = cut(
                                        x = as.numeric(x = data.feature),
                                        breaks = brewer.gran
                                )))
                        }
                }
                md[[f]] <- data.feature
        }
        
        if (is.null(x = split.by)) {
                split.by <- RandomName()
                md[[split.by]] <- factor(x = split.by)
        }
        if (!is.factor(x = md[[split.by]])) {
                md[[split.by]] <- factor(x = md[[split.by]])
        }
        
        if (isTRUE(x = blend)) {
                md <- lapply(
                        X = levels(x = md[[split.by]]),
                        FUN = function(x) {
                                df <- md[as.character(x = md[[split.by]]) == x, , drop = FALSE]
                                no.expression <- features[colMeans(x = df[, features]) == 0]
                                if (length(x = no.expression)) {
                                        stop(
                                                "The following features have no value: ",
                                                paste(no.expression, collapse = ", ")
                                        )
                                }
                                return(cbind(
                                        df[, split.by, drop = FALSE],
                                        BlendExpression(data = df[, features])
                                ))
                        }
                )
                md <- do.call(what = "rbind", args = md)
                features <- setdiff(x = colnames(x = md), y = split.by)
        }
        
        pnames <- unlist(x = lapply(
                X = seq_along(along.with = fov),
                FUN = function(i) {
                        return(if (isTRUE(x = overlap[i])) {
                                fov[i]
                        } else {
                                paste(fov[i], boundaries[[i]], sep = "_")
                        })
                }
        ))
        pdata <- vector(mode = "list", length = length(x = pnames))
        names(x = pdata) <- pnames
        
        for (i in names(x = pdata)) {
                ul <- unlist(x = strsplit(x = i, split = "_"))
                if (length(ul) > 1) {
                        img <- paste(ul[1:length(ul) - 1], collapse = "_")
                        lyr <- ul[length(ul)]
                } else if (length(ul) == 1) {
                        img <- ul[1]
                        lyr <- "centroids"
                } else {
                        stop("the length of ul is 0. please check.")
                }
                if (is.na(x = lyr)) {
                        lyr <- boundaries[[img]]
                }
                pdata[[i]] <- lapply(
                        X = lyr,
                        FUN = function(l) {
                                df <- fortify(model = object[[img]][[l]])
                                df <- df[df$cell %in% cells, , drop = FALSE]
                                if (!is.null(x = md)) {
                                        df <- merge(x = df, y = md, by.x = "cell", by.y = 0, all.x = TRUE)
                                }
                                df$cell <- paste(l, df$cell, sep = "_")
                                df$boundary <- l
                                return(df)
                        }
                )
                pdata[[i]] <- if (!is_na(x = pdata[[i]])) {
                        do.call(what = "rbind", args = pdata[[i]])
                } else {
                        unlist(x = pdata[[i]])
                }
        }
        
        if (!is.null(x = molecules)) {
                molecules <- .MolsByFOV(object = object, fov = fov, molecules = molecules)
                mdata <- vector(mode = "list", length = length(x = fov))
                names(x = mdata) <- fov
                for (img in names(x = mdata)) {
                        idata <- object[[img]]
                        if (!img %in% names(x = molecules)) {
                                mdata[[img]] <- NULL
                                next
                        }
                        if (isTRUE(x = crop[img])) {
                                idata <- Overlay(x = idata, y = idata)
                        }
                        imols <- gsub(
                                pattern = paste0("^", Key(object = idata)),
                                replacement = "",
                                x = molecules[[img]]
                        )
                        mdata[[img]] <- FetchData(object = idata, vars = imols, nmols = nmols)
                }
        } else {
                mdata <- NULL
        }
        
        if (isTRUE(x = blend)) {
                ncol <- 4
                color.matrix <- BlendMatrix(
                        two.colors = cols[2:3],
                        col.threshold = blend.threshold,
                        negative.color = cols[1]
                )
                cols <- cols[2:3]
                colors <- list(
                        color.matrix[, 1],
                        color.matrix[1, ],
                        as.vector(x = color.matrix)
                )
                blend.legend <- BlendMap(color.matrix = color.matrix)
        }
        
        limits <- switch(
                EXPR = scale,
                all  = range(unlist(x = md[, features])),
                NULL
        )
        
        plots <- vector(
                mode = "list",
                length = length(x = levels(x = md[[split.by]]))
        )
        names(x = plots) <- levels(x = md[[split.by]])
        
        for (i in seq_along(along.with = levels(x = md[[split.by]]))) {
                ident <- levels(x = md[[split.by]])[i]
                plots[[ident]] <- vector(mode = "list", length = length(x = pdata))
                names(x = plots[[ident]]) <- names(x = pdata)
                
                if (isTRUE(x = blend)) {
                        blend.key <- suppressMessages(
                                expr = blend.legend +
                                        scale_y_continuous(
                                                sec.axis = dup_axis(
                                                        name = ifelse(
                                                                test = length(x = levels(x = md[[split.by]])) > 1,
                                                                yes  = ident,
                                                                no   = ""
                                                        )
                                                ),
                                                expand = c(0, 0)
                                        ) +
                                        labs(
                                                x = features[1L],
                                                y = features[2L],
                                                title = if (i == 1L) {
                                                        paste("Color threshold:", blend.threshold)
                                                } else {
                                                        NULL
                                                }
                                        ) +
                                        no.right
                        )
                }
                
                for (j in seq_along(along.with = pdata)) {
                        key <- names(x = pdata)[j]
                        img <- unlist(x = strsplit(x = key, split = "_"))[1L]
                        plots[[ident]][[key]] <- vector(
                                mode = "list",
                                length = length(x = features) + ifelse(
                                        test = isTRUE(x = blend),
                                        yes  = 1L,
                                        no   = 0L
                                )
                        )
                        
                        data.base <- pdata[[j]][
                                as.character(x = pdata[[j]][[split.by]]) == ident,
                                ,
                                drop = FALSE
                        ]
                        
                        for (y in seq_along(along.with = features)) {
                                feature <- features[y]
                                
                                # work on a copy for this feature
                                data.plot <- data.base
                                
                                # NEW: reorder rows so high (or |high|) values are drawn last / on top
                                if (isTRUE(order) &&
                                    !isTRUE(blend) &&
                                    feature %in% colnames(data.plot) &&
                                    feature %in% colnames(md.orig)) {
                                        
                                        # strip "boundary_" prefix from cell id to get original cell name
                                        cell_orig <- sub("^[^_]+_", "", data.plot$cell)
                                        
                                        ord_vals <- md.orig[cell_orig, feature, drop = TRUE]
                                        
                                        if (order.by == "abs") {
                                                ord_vals <- abs(ord_vals)
                                        }
                                        
                                        o <- order(ord_vals)  # low first, high last
                                        data.plot <- data.plot[o, , drop = FALSE]
                                }
                                
                                cols.use <- if (isTRUE(x = blend)) {
                                        cc <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
                                        colors[[y]][sort(unique(x = cc))]
                                } else {
                                        NULL
                                }
                                
                                colnames(data.plot) <- gsub("-", "_", colnames(data.plot))
                                p <- SingleImagePlot(
                                        data = data.plot,
                                        col.by = gsub("-", "_", feature),
                                        size = size,
                                        col.factor = blend,
                                        cols = cols.use,
                                        molecules = mdata[[img]],
                                        mols.size = mols.size,
                                        mols.cols = mols.cols,
                                        alpha = alpha,
                                        border.color = border.color,
                                        border.size = border.size,
                                        dark.background = dark.background
                                ) +
                                        CenterTitle() +
                                        labs(fill = feature)
                                
                                if (isTRUE(x = blend)) {
                                        p <- p + guides(fill = "none")
                                }
                                if (isTRUE(coord.fixed)) {
                                        p <- p + coord_fixed()
                                }
                                if (!isTRUE(x = axes)) {
                                        p <- p + NoAxes(panel.background = element_blank())
                                } else if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 1L) {
                                        if (y != 1L) {
                                                p <- p + theme(
                                                        axis.line.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.text.y = element_blank(),
                                                        axis.title.y.left = element_blank()
                                                )
                                        }
                                        if (i != length(x = levels(x = md[[split.by]]))) {
                                                p <- p + theme(
                                                        axis.line.x = element_blank(),
                                                        axis.ticks.x = element_blank(),
                                                        axis.text.x = element_blank(),
                                                        axis.title.x = element_blank()
                                                )
                                        }
                                }
                                
                                if (!isTRUE(x = blend)) {
                                        if (length(x = cols) == 1L) {
                                                p <- p + scale_fill_brewer(palette = cols)
                                        } else {
                                                cols.grad <- cols
                                                fexp <- data.base[data.base[[split.by]] == ident, feature, drop = TRUE]
                                                fexp <- unique(x = fexp)
                                                if (length(x = fexp) == 1L) {
                                                        warning(
                                                                "All cells have the same value (",
                                                                fexp,
                                                                ") of ",
                                                                feature,
                                                                call. = FALSE,
                                                                immediate. = TRUE
                                                        )
                                                        if (fexp == 0) {
                                                                cols.grad <- cols.grad[1L]
                                                        }
                                                }
                                                if (scale == "feature") {
                                                        limits <- range(pdata[[j]][[feature]])
                                                }
                                                p <- p + ggplot2::scale_fill_gradientn(
                                                        colors = cols.grad,
                                                        guide = "colorbar",
                                                        limits = limits
                                                )
                                        }
                                }
                                
                                p <- p + if (i == 1L) {
                                        ggplot2::labs(title = feature)
                                } else {
                                        ggplot2::labs(title = NULL)
                                }
                                
                                plots[[ident]][[key]][[y]] <- p
                        }
                        
                        if (isTRUE(x = blend)) {
                                plots[[ident]][[key]][[length(x = plots[[ident]][[key]])]] <- blend.key
                        } else if (length(x = levels(x = md[[split.by]])) > 1L) {
                                plots[[ident]][[key]][[length(features)]] <- suppressMessages(
                                        expr = plots[[ident]][[key]][[length(features)]] +
                                                scale_y_continuous(sec.axis = dup_axis(name = ident)) +
                                                no.right
                                )
                        }
                }
                
                plots[[ident]] <- unlist(x = plots[[ident]], recursive = FALSE, use.names = FALSE)
        }
        
        plots <- unlist(x = plots, recursive = FALSE, use.names = FALSE)
        
        if (isTRUE(x = combine)) {
                if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 1L) {
                        plots <- wrap_plots(
                                plots,
                                ncol = ifelse(test = isTRUE(x = blend), yes = 4L, no = length(x = features)),
                                nrow = length(x = levels(x = md[[split.by]])),
                                guides = "collect"
                        )
                } else {
                        plots <- wrap_plots(plots)
                }
        }
        return(plots)
}

# after defining ImageFeaturePlot_ordered very important to run this:
environment(ImageFeaturePlot_ordered) <- environment(Seurat::ImageFeaturePlot)


### 2) User-editable paths -----------------------------------------------------

cfg <- list(
        microglia_rds = NA_character_,
        output_dir_plots = NA_character_)

if (is.na(cfg$output_dir_plots) || !nzchar(cfg$output_dir_plots))
        stop("Set cfg$output_dir_plots")

dir.create(cfg$output_dir_plots, recursive = TRUE, showWarnings = FALSE)

### 3) Color palettes ----------------------------------------------------------

ilae_levels <- c("P", "0", "1", "2", "3")

pal_ILAE_score <- c("P" = "#C67BFF", "0" = "#F8766C",
                    "1" = "#7BAD00", "2" = "#00C0C5", "3" = "#E69F00")

pal_cluster_named <- c("Homeo" = "#ddcc77", "DIM" = "#117733",
                       "BAM" = "#88ccee", "DAM" = "#cc6677",
                       "Inflam.DAM" = "#aa4499", "T-cells" = "#332288",
                       "PeriNeurM" = "#44aa99")

group_titles <- c("P" = "Postmortem", "0" = "ILAE 0",
                  "1" = "ILAE 1", "2" = "ILAE 2", "3" = "ILAE 3")

### 4) Load and validate data --------------------------------------------------

# Load the annotated microglia Seurat object
datM_clean <- load_if_missing("datM_clean", cfg$microglia_rds, class = "Seurat")

# Check required metadata
required_meta <- c("cellstates", "ILAE_score", "region_md", "Sample")

if (!all(required_meta %in% colnames(datM_clean@meta.data)))
        stop("Missing metadata: ", paste(setdiff(required_meta, colnames(datM_clean@meta.data)), collapse = ", "))

if (anyNA(datM_clean$cellstates))
        stop("Missing cell-state annotations")

if (!all(unique(as.character(datM_clean$cellstates)) %in% names(pal_cluster_named)))
        stop("Cell-state palette does not cover all annotations")

if (!all(unique(as.character(datM_clean$ILAE_score)) %in% ilae_levels))
        stop("Unexpected ILAE score in microglia metadata")

# Apply consistent factor ordering
datM_clean$cellstates <- factor(as.character(datM_clean$cellstates),
                                levels = names(pal_cluster_named))

datM_clean$ILAE_score <- factor(as.character(datM_clean$ILAE_score),
                                levels = ilae_levels)

# Use normalized SCT expression for plotting gene-expression
DefaultAssay(datM_clean) <- "SCT"

### 5) Spatial cell-state maps -------------------------------------------------

# Purpose:
# Visualize the anatomical distribution of annotated microglia cell states.
#
# Only five states are displayed to match the original figure.

keep_states <- c("Homeo", "DIM", "BAM", "Inflam.DAM", "DAM")

FOV_plots <- subset(datM_clean, subset = cellstates %in% keep_states)

# Representative samples: one per ILAE group
representative_fovs <- list("P" = "M70", "0" = "I1",
                            "1" = "L15", "2" = "L1", "3" = "D9")

# Group all samples by ILAE score
sample_groups <- split(as.character(FOV_plots$Sample),
                       as.character(FOV_plots$ILAE_score))

sample_groups <- lapply(sample_groups, function(x) sort(unique(x)))
sample_groups <- sample_groups[ilae_levels[ilae_levels %in% names(sample_groups)]]

# Check that each sample belongs to one ILAE group and has a FOV
if (anyDuplicated(unlist(sample_groups)))
        stop("A sample occurs in multiple ILAE groups")

if (!all(unlist(sample_groups) %in% Images(FOV_plots)))
        stop("Some samples have no corresponding FOV image")

if (!all(unlist(representative_fovs) %in% Images(FOV_plots)))
        stop("Some representative FOVs are missing")

# Plot 1: Representative FOV for each ILAE group
spatial_representative <- fov_figure(
        representative_fovs,
        n_col = c("P" = 1, "0" = 1, "1" = 1, "2" = 1, "3" = 1),
        show_label = FALSE, body_h = 10)

print(spatial_representative)

# Plot 2: All samples grouped by ILAE score
spatial_all <- fov_figure(
        sample_groups,
        n_col = c("P" = 2, "0" = 2, "1" = 2, "2" = 2, "3" = 1),
        show_label = TRUE, body_h = 20)

print(spatial_all)

### 6) Spatial gene-expression maps -------------------------------------------

# Purpose:
# Visualize spatial expression of CD83 and CCL4 in selected samples.
#
# Expression is displayed using ImageFeaturePlot_ordered().
# Higher-expression cells are drawn on top of lower-expression cells.
#
# q2 and q95 specify the 2nd and 95th percentile expression cutoffs.
# Viridis is used for expression intensity.
#
# Each sample/gene combination is saved as a transparent PNG.

genes_of_interest <- c("CD83", "CCL4")
samples <- c("M70", "I1", "L15", "L1", "D9")

output_dir <- file.path(cfg$output_dir_plots, "ImageFeaturePlots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(samples %in% Images(datM_clean)))
        stop("Some requested spatial expression FOVs are missing")

if (!all(genes_of_interest %in% rownames(datM_clean[["SCT"]])))
        stop("Some requested genes are missing from the SCT assay")

# Plot and save spatial expression for every sample/gene combination
for (sample in samples) {
        for (gene in genes_of_interest) {
                
                p <- ImageFeaturePlot_ordered(
                        datM_clean, features = gene, fov = sample,
                        size = 1.5, border.size = NA,
                        axes = FALSE, dark.background = FALSE,
                        min.cutoff = "q2", max.cutoff = "q95",
                        order = TRUE) +
                        scale_fill_viridis_c(option = "viridis") +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position = "right")
                
                ggsave(file.path(output_dir, paste0(sample, "_", gene, "_viridis.png")),
                       plot = p, width = 6, height = 6, dpi = 300,
                       bg = "transparent")
        }
}

### 7) Regional gene expression ------------------------------------------------

# Purpose:
# Compare normalized CD83 and CCL4 expression across hippocampal
# anatomical regions and ILAE groups.
#
# Each tile represents mean normalized expression across all microglia
# cells belonging to one anatomical region and ILAE group.

# Extract SCT-normalized expression, retaining cell barcode alignment
expr_df <- FetchData(datM_clean, vars = genes_of_interest,
                     layer = "data", clean = FALSE)

meta <- datM_clean@meta.data
meta[rownames(expr_df), genes_of_interest] <- expr_df[, genes_of_interest]

# Calculate mean expression per anatomical region and ILAE group
tile_df_ilae <- meta %>%
        filter(!is.na(region_md), !is.na(ILAE_score)) %>%
        group_by(ILAE_score, region_md) %>%
        summarise(across(all_of(genes_of_interest),
                         ~ mean(.x, na.rm = TRUE),
                         .names = "mean_{.col}"),
                  .groups = "drop")

# Reshape to long format for plotting
tile_long_ilae <- tile_df_ilae %>%
        pivot_longer(cols = starts_with("mean_"),
                     names_to = "gene", values_to = "mean_expr") %>%
        mutate(gene = sub("^mean_", "", gene),
               gene = factor(gene, levels = genes_of_interest),
               ILAE_score = factor(ILAE_score, levels = ilae_levels))

# Preserve anatomical region ordering from the Seurat metadata
region_levels <- if (is.factor(datM_clean$region_md))
        levels(datM_clean$region_md) else unique(as.character(datM_clean$region_md))

tile_long_ilae$region_md <- factor(tile_long_ilae$region_md,
                                   levels = region_levels)

# Shared color scale across both genes
min_val <- min(tile_long_ilae$mean_expr, na.rm = TRUE)
max_val <- max(tile_long_ilae$mean_expr, na.rm = TRUE)

# Heatmap:
# X-axis = anatomical region
# Y-axis = ILAE group
# Fill = mean normalized expression across pooled microglia cells
# Separate panel for each gene

tileplot <- ggplot(tile_long_ilae,
                   aes(x = region_md, y = ILAE_score, fill = mean_expr)) +
        geom_tile(color = "grey90") +
        scale_fill_gradient(low = "white", high = "blue",
                            limits = c(min_val, max_val),
                            breaks = c(min_val, max_val),
                            labels = scales::number(c(min_val, max_val), accuracy = 0.1)) +
        facet_wrap(~ gene, ncol = 1, strip.position = "right") +
        labs(x = "Region", y = "ILAE score", fill = "Average expression") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              strip.text = element_text(face = "italic", size = 11))

print(tileplot)
