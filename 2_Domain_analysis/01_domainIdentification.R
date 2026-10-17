################################################################################
# 01_domainIdentification.R
#
# Xenium spatial domain atlasing (PRECAST)
#  - Bins raw Xenium transcripts into spatial bins (default 150x150um)
#  - Runs PRECAST across a range of k to find spatial domains
#  - Maps PRECAST clusters to anatomical regions (final choice: k = 23)
#  - Shiny app for interactive manual correction of region calls
#  - Joins final region calls back onto the cell-level MegaData object
#
# This script:
#   1) Loads MegaData and bins raw transcripts per sample into a Seurat list
#   2) Filters low-count bins, normalizes, builds a PRECASTObj
#   3) Runs PRECAST at k = 23; optional k = 1:30 model selection/AIC
#   4) Optionally plots PRECAST clusters per sample/k to compare models
#   5) Assigns clusters -> regions at k = 23, builds seumerged
#   6) Launches a Shiny app for manual region relabeling (incl. "hilus")
#   7) Joins region / region_md back onto MegaData (cell level)
#   8) Saves session info
#
# Original author: Mathias Delhaye, October 2026
# Code cleanup and standardization: Larissa Kraus (October 2026)

### 0) Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
        library(conflicted)
        library(tidyverse)
        library(Seurat)
        
        library(future)
        library(scCustomize)
        library(PRECAST)
        
        library(grid)
        library(plotly)
        library(DT)
        library(shiny)
        library(data.table)
        library(Matrix)
})

### PACKAGE CONFLICTS ---------------------------------------------------------
conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename, dplyr::lag)


### 1) User-editable config ----------------------------------------------------

cfg <- list(
        
        # Annotated MegaData object available from GEO
        megadata_file = NA_character_,
        rds_dir = NA_character_,
        
        # Plot outputs
        fig_dir = NA_character_,
        
        # Binning parameters - bins are bin_width x bin_height um squares
        bin_width  = 150,
        bin_height = 150,
        
        # QC filter: minimum transcript count per bin to keep
        min_nCount_bin = 10,
        
        # PRECAST k search range + shortlist compared visually + final choice
        k_range     = 1:30,
        k_shortlist = c(21, 22, 23),
        k_final     = 23,
        
        # PRECAST run settings
        precast_sigma_equal = FALSE,
        precast_coreNum     = 1,
        precast_maxIter     = 30,
        
        # Region palettes (used for region plots and the Shiny app)
        col_regions = setNames(
                c("#3D0F99", "#C7B8E6", "#0F8299", "#B8DEE6", "#54990F", "#CFE6B8",
                  "#333333", "#990F26", "#CCAA7A", "#99600F", "#FFD300"),
                c("GCL", "ML", "CA2_CA4", "CA1", "SL.SR.SLM", "GABA", "SUB", "WM", "pia", "vascular", "void")
        ),
        region_levels =  c("GCL","ML","CA2_CA4","CA1","SUB","GABA","SL.SR.SLM","WM","vascular","pia","void"),
        
        # smpl_ILAE factor level order (hand-set, matches manuscript sample ordering)
        smpl_ILAE_levels = c(
                "L5_P", "L10_P", "M56_P", "M70_P", "P28_P", "P51_P", "P60_P", "P71_P",
                "B1_0", "D8_0", "D11_0", "F1_0", "H1_0", "I1_0", "J1_0", "L20_0",
                "C1_1", "C2_1", "E008_1", "L14_1", "L15_1", "L18_1",
                "D12_2", "D13_2", "L1_2", "L3_2", "L21_2",
                "D9_3", "D14_3", "E015_3", "L2_3"
        ),
        
        # Output filenames
        f_precastobj              = "05_PRECASTobj.rds",
        f_precast_list            = "06_precast_list_k1to30.rds",
        f_seumerged               = "07_seumerged.rds",
        f_megadata_metadata       = "08_megadata_metadata_md.rds",
        f_fig_aic                 = "01_precast_AICvskcluster.png"  
)


# Create output dirs
dir.create(cfg$rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)

# Memory for large objects (adjust if needed)
options(future.globals.maxSize = 25000 * 1024^2)

# Fixed 30-color cluster palette (no palette-generating functions)
cluster_blues <- c("#F7FBFF", "#EDF4FA", "#E2EDF6", "#D7E7F2", "#CADFEB", "#BCD7E7", "#ACCFE3", "#9CC6DF", "#8CBEDB", "#7BB6D7", "#6AACD1", "#5BA3CB", "#4E9AC5", "#4292C6", "#3889BE", "#2F80B6", "#2677AD", "#1F6FA6", "#19669C", "#135D92", "#0F5588", "#0B4D7E", "#084675", "#083E6C", "#083663", "#082F5A", "#082751", "#081F48", "#07183F", "#061136")

### 2) Helper functions --------------------------------------------------------

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


# 2.1 Bin raw Xenium transcripts for one sample into a bin x gene count matrix,
#     wrapped as a Seurat object with bin-grid coordinates in meta.data
bin_one_sample <- function(im, megadata, genes_list, bin_width = 150, bin_height = 150) {
        fov <- megadata@images[[im]]@molecules$molecules
        
        # Pull all gene coordinates via lapply + rbindlist instead of growing a df in a loop
        transcript_coord <- rbindlist(lapply(genes_list, function(gene) {
                coords <- fov[[gene]]@coords
                data.table(gene = gene, x = coords[, 1], y = coords[, 2])
        }))
        
        # Vectorized binning
        transcript_coord[, `:=`(
                bin_x = ceiling(x / bin_width),
                bin_y = ceiling(y / bin_height)
        )]
        
        number_bin_x <- ceiling(max(transcript_coord$x) / bin_width)
        number_bin_y <- ceiling(max(transcript_coord$y) / bin_height)
        total_bins   <- number_bin_x * number_bin_y
        
        transcript_coord[, bin_number := (bin_y - 1L) * number_bin_x + bin_x]
        
        # Force factor levels to cover the FULL range of genes/bins so every
        # gene and every bin appears in the matrix, even with zero transcripts.
        # No need to manually add missing columns afterwards.
        gene_f <- factor(transcript_coord$gene, levels = sort(genes_list))
        bin_f  <- factor(transcript_coord$bin_number, levels = seq_len(total_bins))
        
        # sparseMatrix sums duplicate (i,j) triplets automatically -> this IS the
        # group_by + tally + pivot_wider step, done in one call.
        transcript_genexbin_mat <- sparseMatrix(
                i = as.integer(gene_f),
                j = as.integer(bin_f),
                x = 1,
                dims = c(nlevels(gene_f), nlevels(bin_f)),
                dimnames = list(levels(gene_f), paste0(im, "_", seq_len(total_bins)))
        )
        
        # Bin coordinate grid (col varies fastest, matching bin_number formula)
        bin_grid <- expand.grid(col = 1:number_bin_x, row = 1:number_bin_y)
        rownames(bin_grid) <- paste0(im, "_", seq_len(total_bins))
        bin_grid$x_window <- bin_grid$col * bin_width
        bin_grid$y_window <- bin_grid$row * bin_height
        
        seu <- CreateSeuratObject(
                counts = transcript_genexbin_mat,
                assay = "RNA",
                meta.data = bin_grid
        )
        
        levels(seu@meta.data$orig.ident) <- im
        seu
}

# 2.2 Ensure every factor level is present at least once (with NA coords, so it
#     is skipped by geom_point) so ggplot always draws a complete legend even
#     when a sample has zero bins in some level (e.g. a cluster/region absent
#     from that sample)
pad_missing_levels <- function(df, level_col, all_levels) {
        df[[level_col]] <- factor(df[[level_col]], levels = all_levels)
        
        present <- unique(df[[level_col]][!is.na(df[[level_col]])])
        missing_levels <- setdiff(all_levels, present)
        
        if (length(missing_levels) > 0) {
                dummy_rows <- df[rep(1, length(missing_levels)), ]
                dummy_rows[[level_col]] <- factor(missing_levels, levels = all_levels)
                dummy_rows$col <- NA
                dummy_rows$row <- NA
                df <- rbind(df, dummy_rows)
        }
        df
}

# 2.3 Shared spatial-domain plot: one sample, colored by a discrete metadata
#     column (cluster / region / region_md_hilus), saved to disk. Replaces the
#     three near-duplicate ggplot blocks from the original script.
plot_spatial_domains <- function(df,
                                 color_col,
                                 palette,
                                 out_file,
                                 title,
                                 legend_title    = color_col,
                                 highlight_level = NULL,
                                 highlight_color = "#d62728FF",
                                 width  = 6,
                                 height = 6) {
        
        # Optionally recolor one level (used to spotlight a single cluster at a time)
        if (!is.null(highlight_level)) {
                palette[highlight_level] <- highlight_color
        }
        
        max_dim <- max(c(df$col, df$row), na.rm = TRUE)
        point_size_pt <- (6 * 72) / max_dim
        
        p <- ggplot(df, aes(x = col, y = row, color = .data[[color_col]], fill = .data[[color_col]])) +
                geom_point(shape = 15, size = point_size_pt, na.rm = TRUE) +
                scale_color_manual(values = palette, drop = FALSE, guide = "none") +
                scale_fill_manual(
                        name = legend_title,
                        values = palette,
                        drop = FALSE,
                        guide = guide_legend(
                                ncol = 1,
                                override.aes = list(
                                        shape = 22,      # Square shape
                                        color = "black", # Black border outline
                                        stroke = 0.5,    # Outline thickness
                                        size = 4
                                )
                        )
                ) +
                scale_y_reverse(limits = c(max_dim, 1), expand = c(0, 0)) +
                scale_x_continuous(limits = c(1, max_dim), expand = c(0, 0)) +
                coord_fixed(expand = FALSE) +
                theme(
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        plot.background = element_rect(fill = "white", color = NA),
                        plot.title.position = "panel",
                        plot.title = element_text(
                                size = 20,
                                hjust = 0.01,
                                margin = margin(t = 0, b = -25),
                                face = "bold"
                        ),
                        legend.text = element_text(size = 8),
                        legend.title = element_text(size = 10, face = "bold"),
                        legend.box.margin = margin(l = -5),
                        legend.margin = margin(l = 0),
                        legend.key = element_blank(),
                        aspect.ratio = 1
                ) +
                labs(title = title, fill = legend_title) +
                annotation_custom(
                        grob = rectGrob(gp = gpar(col = "black", fill = NA, lwd = 2.5)),
                        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
                )
        
        ggsave(filename = out_file, plot = p, dpi = 300, width = width, height = height)
        invisible(p)
}

### 3) Load inputs -------------------------------------------------------------

# make sure to move megadata in the newly created rds folder
megadata <- load_if_missing("megadata", cfg$megadata_file, class = "Seurat")

# Names of every sample in megadata, in order assigned in megadata
samples_names <- names(megadata@images)

# Genes analyzed with Xenium (266 panel + 100 epilepsy-specific)
genes_list <- Features(megadata)

### 4) MegaData -> binned seulist ----------------------------------------------

message("Binning transcripts into ", cfg$bin_width, "x", cfg$bin_height,
        "um spatial bins for ", length(samples_names), " samples")

seulist <- purrr::map(
        samples_names, bin_one_sample,
        megadata = megadata, genes_list = genes_list,
        bin_width = cfg$bin_width, bin_height = cfg$bin_height
)
names(seulist) <- samples_names

# Filter bins by min transcript count + log-normalize
message("Filtering to bins with nCount_RNA > ", cfg$min_nCount_bin, " and log-normalizing")

seulist_log30 <- purrr::map(seulist, function(seu) {
        seu %>%
                subset(subset = nCount_RNA > cfg$min_nCount_bin) %>%
                NormalizeData(verbose = FALSE) %>%
                FindVariableFeatures(verbose = FALSE) %>%
                ScaleData(verbose = FALSE)
})

### 5) seulist -> PRECAST (default k = 23) -------------------------------------------

# Create PRECAST object
preobj <- CreatePRECASTObject(seuList = seulist_log30)

# Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
PRECASTObj <- AddAdjList(preobj, platform = "Other_SRT")

# Add a model setting in advance for a PRECASTObj object. verbose = TRUE helps
# output the information in the algorithm.
PRECASTObj <- AddParSetting(
        PRECASTObj,
        Sigma_equal = cfg$precast_sigma_equal,
        coreNum     = cfg$precast_coreNum,
        maxIter     = cfg$precast_maxIter,
        verbose     = TRUE
)
saveRDS(PRECASTObj, file.path(cfg$rds_dir, cfg$f_precastobj))

# DEFAULT: fit only k = 23 (selected in the original analysis).
# This is sufficient to test/run the full downstream pipeline.
message("Running PRECAST for selected k = ", cfg$k_final)
PRECASTObj_final <- PRECAST(PRECASTObj, K = cfg$k_final)
saveRDS(PRECASTObj_final, file.path(cfg$rds_dir, paste0("06_precast_k", cfg$k_final, ".rds")))

### 6) OPTIONAL: compare k = 1:30 and plot shortlisted models -----------------
# OPTIONAL MODEL SELECTION ONLY. Do not run this section for the standard
# pipeline/test run. Run this section independently AFTER Sections 1-5 if you
# want to reproduce model selection. It does not change cfg$k_final = 23.
# The original analysis compared k = 1:30, visually inspected 21, 22, 23,
# and selected k = 23 for the anatomical region mapping below.

if (FALSE) {  # Change to TRUE only when evaluating multiple k values
        precast_list_k1to30 <- list()
        for (k in cfg$k_range) {
                if (k == cfg$k_final) {
                        precast_list_k1to30[[paste0("k", k)]] <- PRECASTObj_final
                } else {
                        precast_list_k1to30[[paste0("k", k)]] <- PRECAST(PRECASTObj, K = k)
                }
        }
        saveRDS(precast_list_k1to30, file.path(cfg$rds_dir, cfg$f_precast_list))
        
        # Compare fit via AIC across k
        aic_list <- list()
        for (i in seq_along(precast_list_k1to30)) {
                PRECASTObj_k <- SelectModel(precast_list_k1to30[[i]], criteria = "AIC")
                aic_list[i] <- PRECASTObj_k@resList$icMat[2]
                message("k = ", cfg$k_range[i], " -> AIC = ", aic_list[i])
                rm(PRECASTObj_k)
                gc()
        }
        
        aic_df <- data.frame(Index = cfg$k_range, AIC = as.numeric(unlist(aic_list)))
        
        p_aic <- ggplot(aic_df, aes(x = Index, y = AIC)) +
                geom_line() +
                geom_point() +
                scale_y_continuous(labels = scales::scientific) +
                theme_minimal() +
                theme(
                        axis.title = element_text(size = 14),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 10, angle = 45)
                )
        
        out_dir <- file.path(cfg$fig_dir, "01_precast_AIC_plot")
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        ggsave(filename = file.path(out_dir, cfg$f_fig_aic), plot = p_aic, dpi = 300, width = 8, height = 6)
        
        # From the AIC curve, optimal k falls between 21 and 23 (see cfg$k_shortlist)
        
        
        # OPTIONAL: spatial plots for k = 21, 22, 23
        seumerged <- merge(seulist_log30[[1]], seulist_log30[-1])
        
        df_temp <- seumerged@meta.data %>%
                rownames_to_column(var = "id") %>%
                select(all_of(c("id", "col", "row")))
        
        # Loop through candidate k values, plot every cluster of every sample, to
        # visually pick the best k
        for (k in cfg$k_shortlist) {
                
                PRECASTObj_k <- precast_list_k1to30[[paste0("k", k)]]
                
                seuint <- PRECASTObj_k %>%
                        SelectModel() %>%
                        IntegrateSpaData(species = "unknown")
                
                metadata_seuint <- seuint@meta.data %>%
                        rownames_to_column(var = "id") %>%
                        left_join(df_temp, by = "id")
                
                all_clusters <- levels(factor(metadata_seuint$cluster))
                col_clusters_blues <- setNames(cluster_blues[seq_len(k)], all_clusters)
                
                for (smpl in samples_names) {
                        df <- metadata_seuint %>% filter(orig.ident == smpl)
                        df <- pad_missing_levels(df, level_col = "cluster", all_levels = all_clusters)
                        
                        out_dir <- file.path(cfg$fig_dir, "02_clusters_location", paste0("k", k), smpl)
                        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                        
                        for (current_cluster in all_clusters) {
                                plot_spatial_domains(
                                        df              = df,
                                        color_col       = "cluster",
                                        palette         = col_clusters_blues,
                                        highlight_level = current_cluster,
                                        legend_title    = "Cluster",
                                        title           = smpl,
                                        out_file        = file.path(out_dir, paste0(smpl, "_cluster", current_cluster, "_k", k, ".png"))
                                )
                        }
                }
        }
        
}

### 7) Finalize region assignment (k = cfg$k_final) ----------------------------

seuint <- PRECASTObj_final %>%
        SelectModel() %>%
        IntegrateSpaData(species = "unknown")

metadata_seuint <- seuint@meta.data %>%
        select(orig.ident, cluster) %>%
        rownames_to_column(var = "spot_id")

# Cells that survived PRECAST integration (i.e., got a cluster assignment)
kept_cells <- Cells(seuint)
message(length(kept_cells), " / ", ncol(seumerged), " bins kept after PRECAST integration at k = ", cfg$k_final)

# Filter each per-sample object to only these cells, BEFORE merging
seulist_log30_filtered <- lapply(seulist_log30, function(seu) {
        cells_to_keep <- intersect(Cells(seu), kept_cells)
        subset(seu, cells = cells_to_keep)
})

# Re-compute seumerged to omit the missing spots
seumerged <- merge(seulist_log30_filtered[[1]], seulist_log30_filtered[-1])

# Map PRECAST clusters -> anatomical regions
# IMPORTANT: this mapping is dataset/k-specific (k = 23). If you re-run PRECAST
# with different parameters, you must re-validate this mapping.
metadata_seumerged <- seumerged@meta.data %>%
        rownames_to_column(var = "spot_id") %>%
        left_join(metadata_seuint, by = c("spot_id", "orig.ident")) %>%
        filter(!is.na(cluster)) %>%
        mutate(region = case_when(
                cluster == 9  | cluster == 17                 ~ "GCL",
                cluster == 14                                   ~ "ML",
                cluster == 10  | cluster == 12 | cluster == 16 ~ "CA2_CA4",
                cluster == 2  | cluster == 19 ~ "CA1",
                cluster == 5  | cluster == 7  | cluster == 13 ~ "SL.SR.SLM",
                cluster == 3                                   ~ "GABA",
                cluster == 20                                  ~ "SUB",
                cluster == 15 | cluster == 22                 ~ "WM",
                cluster == 8                                  ~ "pia",
                cluster == 4 | cluster == 21 | cluster == 23 ~ "vascular",
                cluster == 1  | cluster == 6  | cluster == 11 | cluster == 18 ~ "void"
        ))

out_dir <- file.path(cfg$fig_dir, "03_regions_location", paste0("k_", cfg$k_final,"_regions"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (smpl in samples_names) {
        df <- metadata_seumerged %>% filter(orig.ident == smpl)
        plot_spatial_domains(
                df           = df,
                color_col    = "region",
                palette      = cfg$col_regions,
                legend_title = "Region",
                title        = smpl,
                out_file     = file.path(out_dir, paste0(smpl, "_k_", cfg$k_final, ".png"))
        )
}

### 8) Add ILAE score + interactive manual correction (Shiny) ------------------

# Generate table with ILAE score per sample
ILAE_table <- megadata@meta.data %>%
        group_by(orig.ident, ILAE_score) %>%
        tally() %>%
        select(orig.ident, ILAE_score) %>%
        mutate(ILAE_score = factor(ILAE_score, levels = c("P", "0", "1", "2", "3")))

# Add the ILAE score to the metadata, and set up the editable region_md column
metadata_seumerged <- metadata_seumerged %>%
        left_join(ILAE_table, by = "orig.ident") %>%
        rename(cluster_precast = cluster) %>%
        mutate(region_md = region)

# 8.1 Shiny app: interactively relabel spots (region_md) by lasso-selecting
#     points on the spatial scatter plot. The final in-session state is written
#     back to `metadata_seumerged` when the app is closed.
ui <- fluidPage(
        titlePanel("Interactive Spot Viewer"),
        tags$script(HTML("
    $(document).on('keydown', '#new_id', function(e) {
      if (e.key === 'Enter' || e.keyCode === 13) {
        e.preventDefault();
        $('#change_id').click();
      }
    });
  ")),
        sidebarLayout(
                sidebarPanel(
                        selectInput("case", "Select Case:", choices = unique(metadata_seumerged$orig.ident)),
                        actionButton("change_id", "Change Identity of Selected Points"),
                        textInput("new_id", "New Identity:", value = "")
                ),
                mainPanel(
                        plotlyOutput("scatterPlot", height = "85vh", width = "100%"),
                        DTOutput("selectedPoints")
                )
        )
)

server <- function(input, output, session) {
        metadata_reactive <- reactiveVal(metadata_seumerged)
        filteredData <- reactive({
                metadata_reactive() %>% filter(orig.ident == input$case)
        })
        
        output$scatterPlot <- renderPlotly({
                df <- filteredData()
                x_range <- range(df$col, na.rm = TRUE)
                
                plot_ly(
                        data = df,
                        x = ~col,
                        y = ~row,
                        type = 'scatter',
                        mode = 'markers',
                        color = ~region_md,
                        colors = cfg$col_regions,
                        text = ~paste("region:", region_md, "<br>Spot:", spot_id),
                        key = ~spot_id,
                        marker = list(size = 6, symbol = "square"),
                        source = "scatter"
                ) %>%
                        layout(
                                autosize = TRUE,
                                xaxis = list(
                                        scaleanchor = "y",
                                        scaleratio = 1,
                                        constrain = "domain",
                                        range = x_range
                                ),
                                yaxis = list(
                                        autorange = "reversed",
                                        constrain = "domain"
                                ),
                                dragmode = "pan"   # pan is now the default drag behavior
                        ) %>%
                        config(
                                scrollZoom = TRUE,
                                displaylogo = FALSE
                        ) %>%
                        htmlwidgets::onRender("
        function(el, x) {

          // ---- marker size scaling on zoom ----
          el.on('plotly_relayout', function(ev) {
            var xa = el.layout.xaxis;
            var rangeNow;
            if (ev['xaxis.range[0]'] !== undefined) {
              rangeNow = Math.abs(ev['xaxis.range[1]'] - ev['xaxis.range[0]']);
            } else if (ev['xaxis.autorange']) {
              rangeNow = Math.abs(xa.range[1] - xa.range[0]);
            } else {
              return;
            }
            el._baselineRange = el._baselineRange || rangeNow;
            var zoomFactor = el._baselineRange / rangeNow;
            var baseSize = 6;
            var newSize = Math.max(2, Math.min(40, baseSize * zoomFactor));
            Plotly.restyle(el, {'marker.size': newSize});
          });

          // ---- keyboard shortcut: space toggles lasso <-> pan ----
          var lassoActive = false;

          function isTyping(target) {
            var tag = target.tagName;
            return tag === 'INPUT' || tag === 'TEXTAREA' || target.isContentEditable;
          }

          document.addEventListener('keydown', function(e) {
            if (isTyping(e.target)) return;
            if (e.repeat) return; // ignore key-repeat

            if (e.code === 'Space') {
              e.preventDefault(); // stop page scroll
              lassoActive = !lassoActive;
              var newMode = lassoActive ? 'lasso' : 'pan';
              console.log('Space pressed -> switching to', newMode);
              Plotly.relayout(el, {dragmode: newMode});
            }
          });
        }
      ")
        })
        
        selected_spots <- reactive({
                event_data("plotly_selected", source = "scatter")
        })
        
        output$selectedPoints <- renderDT({
                if (is.null(selected_spots())) return(NULL)
                spots <- selected_spots()$key
                filteredData() %>% filter(spot_id %in% spots)
        })
        
        observeEvent(input$change_id, {
                sel <- selected_spots()
                if (!is.null(sel) && input$new_id != "") {
                        spots <- sel$key
                        updated <- metadata_reactive() %>%
                                mutate(region_md = ifelse(spot_id %in% spots, input$new_id, region_md))
                        metadata_reactive(updated)
                }
        })
        
        session$onSessionEnded(function() {
                final_metadata <- isolate(metadata_reactive())
                metadata_seumerged <<- final_metadata
        })
}

message("Launching Shiny app for manual region relabeling - close the app window to continue")
shinyApp(ui, server)

seumerged@meta.data <- metadata_seumerged %>% column_to_rownames(var = "spot_id")
rm(metadata_seumerged, metadata_seuint)

seumerged$region    <- factor(seumerged$region,    levels = cfg$region_levels)
seumerged$region_md <- factor(seumerged$region_md, levels = cfg$region_levels)

seumerged@meta.data <- seumerged@meta.data %>%
        unite(smpl_ILAE, c(orig.ident, ILAE_score), sep = "_", remove = FALSE)

seumerged$smpl_ILAE <- factor(seumerged$smpl_ILAE, levels = cfg$smpl_ILAE_levels)

saveRDS(seumerged, file.path(cfg$rds_dir, cfg$f_seumerged))
message("Saved seumerged to: ", file.path(cfg$rds_dir, cfg$f_seumerged))

# Plot final region_md calls
out_dir <- file.path(cfg$fig_dir, "03_regions_location", paste0("k_", cfg$k_final, "_regions_md"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (smpl in samples_names) {
        df <- seumerged@meta.data %>%
                filter(orig.ident == smpl)
        
        plot_spatial_domains(
                df           = df,
                color_col    = "region_md",
                palette      = cfg$col_regions,
                legend_title = "Region",
                title        = smpl,
                out_file     = file.path(out_dir, paste0(smpl, ".png"))
        )
}

### 9) Regions -> MegaData -----------------------------------------------------

# Cell centroid coordinates for every sample
cell_coord <- tibble(cell = character(), x = numeric(), y = numeric())

for (smpl in samples_names) {
        coord <- megadata@images[[smpl]]@boundaries[["centroids"]]@coords
        cell_names <- megadata@images[[smpl]]@boundaries[["centroids"]]@cells
        temp_df <- tibble(cell = cell_names, x = coord[, 1], y = coord[, 2])
        cell_coord <- rbind(cell_coord, temp_df)
        rm(temp_df)
}

# Add bin col/row for each cell, using the same bin size chosen in cfg
cell_coord <- cell_coord %>%
        mutate(col = ceiling(x / cfg$bin_width), row = ceiling(y / cfg$bin_height))

# Attach the sample name for each cell
metadata_megadata <- megadata@meta.data
cell_to_orig <- metadata_megadata %>%
        rownames_to_column(var = "cell") %>%
        select(cell, orig.ident)
cell_coord <- left_join(cell_coord, cell_to_orig, by = "cell")
rm(cell_to_orig)

# Join cell_coord and seumerged metadata by (orig.ident, col, row) to assign
# each cell to its bin's region / region_md / region_md_hilus
metadata_seumerged <- seumerged@meta.data
cell_coord <- cell_coord %>%
        left_join(
                metadata_seumerged %>% select(orig.ident, col, row, region, region_md),
                by = c("orig.ident", "col", "row")
        )

metadata_megadata <- metadata_megadata %>%
        rownames_to_column(var = "cell") %>%
        left_join(cell_coord %>% select(cell, region, region_md), by = "cell") %>%
        column_to_rownames(var = "cell")

saveRDS(metadata_megadata, file.path(cfg$rds_dir, cfg$f_megadata_metadata))
message("Saved MegaData metadata (with regions) to: ", file.path(cfg$rds_dir, cfg$f_megadata_metadata))
