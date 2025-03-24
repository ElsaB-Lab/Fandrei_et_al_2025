require(reshape2)
require(scDNA)

#-------------#
# utils
#-------------#
# pvalue to code
pval_to_signif_code <- function(p_vals) {
  return(ifelse(p_vals < 0.001, "***",
                ifelse(p_vals < 0.01, "**",
                       ifelse(p_vals < 0.05, "*", ""))))
}

#-------------#
# KM curve
#-------------#
print_and_flush <- function(string) {
    # Print a string in stdout and flush the result (useful to print in Jupyter Notebook in for loops)
    # → Arguments
    #     - string: string to print

    cat(string)
    flush.console()
}

print_p_value_from_survival_curve <- function(data, column_name, value_1_name, value_2_name) {
    # work in progress...
    stmp <- survdiff(as.formula(sprintf('Surv(os_diag_years, os_status) ~ %s', column_name)), data = data[data[,column_name] %in% c(value_1_name, value_2_name),])
    pval <- 1 - pchisq(stmp$chisq, length(stmp$n) - 1)
    
    print_and_flush(sprintf('%-21s and %21s: p = %.3f\n', value_1_name, value_2_name, pval))
}

print_p_value_from_boxplot <- function(data, column_name, value_1_name, value_2_name, clinical_name) {
    # work in progress...
    
    vec1 = data[ which(data[,column_name] %in% value_1_name), clinical_name]
    vec2 = data[ which(data[,column_name] %in% value_2_name), clinical_name]

    pval <- t.test(vec1,vec2)$p.val 
    
    print_and_flush(sprintf('%-21s and %21s: p = %.3f\n', value_1_name, value_2_name, pval))
}

km_plot <- function(
    x,
    time,
    status,
    groups,
    palette,
    title,
    pval_coord,
    limits=c(0,21),
    conf = F,
    break.by=10,
    plot_title = NULL
  ) {
  
  surv <- as.formula(paste0("Surv(", time, ",", status, ") ~", groups))
  surv_fit <- survfit(surv, data=x)
  surv_fit$call$formula <- surv
  
  print(surv_fit)
  
  g_surv <- ggsurvplot(surv_fit, data=x, 
                       color= groups, 
                       legend.title=title,
                       size=.75,
                       conf.int= conf,
                       xlim=limits,
                       break.time.by=break.by, 
                       pval=T,  pval.coord=pval_coord, pval.size=5,
                       surv.median.line="hv",
                       palette=palette,
                       risk.table=TRUE,
                       tables.col=groups, tables.y.text=FALSE,
                       risk.table.title="No. at risk",risk.table.fontsize=6,
                       ggtheme = theme_classic(),
  )
  
  g <- g_surv$plot +
    theme_classic() +
    gtheme() +
    theme(legend.key.width = unit(1.1,"cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_rect(linewidth=1, fill=NA),
          strip.text.x = element_text(size=10),
          plot.title = element_text(face="bold")
    ) + 
    ylab("OS Probability") + 
    guides(color = guide_legend(override.aes = list(shape = NA), nrow=2)) +
    xlab("Months after diagnosis") +
    scale_color_manual(values=palette, limits=levels(x[[groups]])) +
    ggtitle(plot_title)
  
  gt <- g_surv$table +
    theme(
      axis.line=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      legend.position="none",
      plot.title = element_text(size=18)
    ) +
    xlab("") +
    ylab("") +
    scale_color_manual(values=palette, limits=levels(x[[groups]]))
  
  g_g <- ggarrange(g + theme(legend.position="top"), gt +theme(plot.margin = unit(c(-0.5, 0, .5, .5), "cm")),
                   ncol=1, nrow=2, heights=c(5,1.5))
  
  return(g_g)
}

#-------------#
# Import per sample tapestri Seurat object
#-------------#
import_sample_norm_prot = function(sample) {
  
  # Using h5 object from single cell DNA pipeline 
  h5 = Sys.glob("~/Nextcloud/PPM1D/sc_tapestri/int_files/*.h5")
  names(h5) = basename(h5) %>% strsplit(".dna") %>% lapply(function(l) l[[1]]) %>% unlist()
  
  sce <- scDNA::tapestri_h5_to_sce(
    file=h5[[sample]],
    variant_set = final_detected[which(final_detected$tapestri_id == sample),],
    GT_cutoff = 0,
    DP_cutoff = 0,
    GQ_cutoff = 0,
    protein = T
  )
  
  sce<-scDNA::enumerate_clones(sce)
  droplet_metadata <- scDNA::extract_droplet_size(sce)
  
  exclude = which(rownames(altExp(sce, "Protein")) %in% "CD22")
  
  altExp(sce, "Protein") = altExp(sce, "Protein")[-exclude,]
 
  background_droplets<-droplet_metadata%>%
     dplyr::filter(Droplet_type=="Empty")%>%
     dplyr::filter(dna_size<1.5&dna_size>0.15)%>%
     pull(Cell)
   
  file <- sce@metadata$file
  protein_sce <- SingleCellExperiment::altExp(sce, "Protein")
  protein_mat <- protein_sce@assays@data$Protein
  print("DSB normalization")
  cells_of_interest <- colnames(protein_mat)
  all_protein_droplets <- rhdf5::h5read(file = file, name = "/all_barcodes/protein_read_counts/layers/read_counts")[-exclude,]
  colnames(all_protein_droplets) <- rhdf5::h5read(file = file, name = "/all_barcodes/protein_read_counts/ra/barcode")
  colnames(all_protein_droplets) <- gsub("-1", "", colnames(all_protein_droplets))
  colnames(all_protein_droplets) <- gsub("-1", "", colnames(all_protein_droplets))
  
  empty_drops_matrix_input <- data.frame(all_protein_droplets) %>%
    dplyr::select(all_of(background_droplets))

  rownames(empty_drops_matrix_input) <- rownames(protein_mat)
  isotype <- grep("IgG", rownames(protein_mat), value = TRUE)

  adt_norm <- dsb::DSBNormalizeProtein(
    cell_protein_matrix = protein_mat,
    empty_drop_matrix = empty_drops_matrix_input,
    denoise.counts = TRUE,
    use.isotype.control = TRUE,
    isotype.control.name.vec = isotype
  )

  assay(protein_sce, "DSB_norm") <- adt_norm

  protein_sce@metadata <- droplet_metadata

  cdata_prot = droplet_metadata %>% dplyr::filter(Cell %in% colnames(protein_mat)) %>% column_to_rownames("Cell") %>% S4Vectors::DataFrame()
  protein_sce@colData <- cdata_prot

  SingleCellExperiment::altExp(sce, "Protein") <- protein_sce

  # altExp(sce, "Protein") = altExp(sce, "Protein")[which(!rownames(altExp(sce, "Protein")) %in% c("IgG2a", "IgG2b")),]

  s<-Seurat::as.Seurat(x = altExp(sce),counts="Protein",data="DSB_norm")
  s<-SeuratObject::RenameAssays(object = s, originalexp = 'Protein')
  Seurat::DefaultAssay(s) <- "Protein"
  s[["Protein"]]@data[s[["Protein"]]@data<1] <- 0

  s <- Seurat::ScaleData(s, assay = "Protein")
  s <- Seurat::FindVariableFeatures(s,assay = "Protein")
  s <- Seurat::RunPCA(s, features = Seurat::VariableFeatures(object = s))

  # s <- Seurat::RunPCA(s, features = rownames(s)[1:16])
  adt.dist = dist(t(adt_norm))
  s[["UMAP"]] <- Seurat::RunUMAP(adt.dist, assay = "Protein", reduction.key = "UMAP_")
  s[["snn"]] <- Seurat::FindNeighbors(adt.dist)$snn
  s<- Seurat::FindClusters(s, resolution = 0.2, graph.name = "snn")

}

#-------------#
# Feature Plot
#-------------#
FeaturePlot_scCustom = function (seurat_object, features, colors_use = viridis_plasma_dark_high, 
                                 na_color = "lightgray", order = TRUE, pt.size = NULL, reduction = NULL, 
                                 na_cutoff = 1e-09, raster = NULL, raster.dpi = c(512, 512), 
                                 split.by = NULL, split_collect = NULL, aspect_ratio = NULL, 
                                 figure_plot = FALSE, num_columns = NULL,
                                 layer = "data", alpha_exp = NULL, alpha_na_exp = NULL, label = FALSE, 
                                 label_feature_yaxis = FALSE, combine = TRUE, ...) 
{
  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(object = seurat_object, meta_col_names = split.by, 
                             print_msg = FALSE, omit_warn = FALSE)[[1]]
  }
  if (is.null(x = split_collect)) {
    if (length(x = features) == 1) {
      split_collect <- TRUE
    }
    else {
      split_collect <- FALSE
    }
  }
  if (!is.null(x = split_collect)) {
    if (length(x = features) > 1 && isTRUE(x = split_collect)) {
      cli_abort(message = "{.code split_collect} cannot be set to {.field TRUE} if the number of features is greater than 1.")
    }
  }
  all_found_features <- features
  if (!is.null(x = split.by)) {
    split.by_length <- length(x = unique(x = seurat_object@meta.data[[split.by]]))
    if (!is.null(x = num_columns) && isTRUE(x = label_feature_yaxis)) {
      cli_warn(message = c("Setting number of columns is not permitted if {.code label_feature_yaxis = TRUE}", 
                           i = "Number of columns be automatically set to number of levels in `split.by` ({.field {split.by_length}})."))
      num_columns <- split.by_length
    }
    if (is.null(x = num_columns)) {
      num_columns <- split.by_length
    }
    num_rows <- ceiling(split.by_length/num_columns)
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.", 
                            `*` = "{.val {split.by}} only contains {.field {split.by_length}} variables.", 
                            i = "Please adjust {.code num_columns} to be less than or equal to {.field {split.by_length}}."))
    }
  }
  if (any(all_found_features %in% colnames(x = seurat_object@meta.data))) {
    cli_warn(message = c("Some of the plotted features are from meta.data slot.", 
                         `*` = "Please check that {.code na_cutoff} param is being set appropriately for those features."))
  }
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 
                           2e+05)
  if (is.null(x = pt.size)) {
    if (!is.null(x = split.by)) {
      cells_by_meta <- data.frame(table(seurat_object@meta.data[, 
                                                                split.by]))
      max_cells <- max(cells_by_meta$Freq)
      pt.size <- AutoPointSize_scCustom(data = max_cells, 
                                        raster = raster)
    }
    else {
      cells_total <- nrow(x = seurat_object@meta.data)
      pt.size <- AutoPointSize_scCustom(data = cells_total, 
                                        raster = raster)
    }
  }
  if (is.null(x = na_cutoff)) {
    na_cutoff <- NA
  }
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
  seurat_version <- packageVersion("Seurat")
  if (!is.null(x = alpha_exp) && seurat_version < "5") {
    colors_use <- alpha(colors_use, alpha_exp)
  }
  if (!is.null(x = alpha_na_exp) && seurat_version < "5") {
    na_color <- alpha(na_color, alpha_na_exp)
  }
  if (!is.null(x = alpha_na_exp) && seurat_version >= "5") {
    cli_warn(message = "{.code alpha_na_exp} is not currently supported for Seurat v5+")
  }
  if (is.null(x = alpha_exp) && seurat_version >= "5") {
    alpha_exp <- 1
  }
  if (is.null(x = split.by) && isTRUE(x = combine)) {
    if (seurat_version >= "5") {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, 
                                           features = all_found_features, order = order, 
                                           pt.size = pt.size, reduction = reduction, raster = raster, 
                                           split.by = split.by, ncol = num_columns, combine = combine, 
                                           raster.dpi = raster.dpi, label = label, alpha = alpha_exp, 
                                           ...) & scale_color_gradientn(colors = colors_use, 
                                                                        limits = c(na_cutoff, NA), na.value = na_color))
    }
    else {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, 
                                           features = all_found_features, order = order, 
                                           pt.size = pt.size, reduction = reduction, raster = raster, 
                                           split.by = split.by, ncol = num_columns, combine = combine, 
                                           raster.dpi = raster.dpi, label = label, ...) & 
                                 scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, 
                                                                                       NA), na.value = na_color))
    }
  }
  if (is.null(x = split.by) && isFALSE(x = combine)) {
    if (seurat_version >= "5") {
      plot_list <- suppressMessages(FeaturePlot(object = seurat_object, 
                                                features = all_found_features, order = order, 
                                                pt.size = pt.size, reduction = reduction, raster = raster, 
                                                split.by = split.by, ncol = num_columns, combine = combine, 
                                                raster.dpi = raster.dpi, label = label, alpha = alpha_exp, 
                                                ...))
      plot <- lapply(1:length(x = plot_list), function(i) {
        p[[i]] <- suppressMessages(p[[i]] + scale_color_gradientn(colors = colors_use, 
                                                                  limits = c(na_cutoff, NA), na.value = na_color))
      })
    }
    else {
      plot_list <- suppressMessages(FeaturePlot(object = seurat_object, 
                                                features = all_found_features, order = order, 
                                                pt.size = pt.size, reduction = reduction, raster = raster, 
                                                split.by = split.by, ncol = num_columns, combine = combine, 
                                                raster.dpi = raster.dpi, label = label, ...))
      plot <- lapply(1:length(x = plot_list), function(i) {
        p[[i]] <- suppressMessages(p[[i]] + scale_color_gradientn(colors = colors_use, 
                                                                  limits = c(na_cutoff, NA), na.value = na_color))
      })
    }
  }
  if (!is.null(x = split.by) && length(x = all_found_features) == 
      1) {
    feature_data <- FetchData(object = seurat_object, vars = all_found_features, 
                              layer = layer)
    max_exp_value <- max(feature_data)
    min_exp_value <- min(feature_data)
    if (seurat_version >= "5") {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, 
                                           features = all_found_features, order = order, 
                                           pt.size = pt.size, reduction = reduction, raster = raster, 
                                           split.by = split.by, raster.dpi = raster.dpi, 
                                           label = label, alpha = alpha_exp, ...) & scale_color_gradientn(colors = colors_use, 
                                                                                                          limits = c(na_cutoff, max_exp_value), na.value = na_color, 
                                                                                                          name = all_found_features)) & RestoreLegend() & 
        theme(axis.title.y.right = element_blank())
    }
    else {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, 
                                           features = all_found_features, order = order, 
                                           pt.size = pt.size, reduction = reduction, raster = raster, 
                                           split.by = split.by, raster.dpi = raster.dpi, 
                                           label = label, ...) & scale_color_gradientn(colors = colors_use, 
                                                                                       limits = c(na_cutoff, max_exp_value), na.value = na_color, 
                                                                                       name = all_found_features)) & RestoreLegend() & 
        theme(axis.title.y.right = element_blank())
    }
    if (isTRUE(x = label_feature_yaxis)) {
      plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns)
      plot <- plot & theme(legend.title = element_blank())
      plot <- suppressMessages(plot + scale_y_continuous(sec.axis = dup_axis(name = all_found_features))) + 
        No_Right()
    }
    else {
      if (isTRUE(x = split_collect)) {
        if (hasArg("keep.scale")) {
          cli_abort(message = "The parameter {.code keep.scale} cannot be set different from default if {.code split_collect - TRUE}.")
        }
        plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns, 
                                   guides = "collect")
      }
      else {
        plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns)
      }
    }
  }
  if (!is.null(x = split.by) && length(x = all_found_features) > 
      1) {
    plot_list <- lapply(1:length(x = all_found_features), 
                        function(i) {
                          feature_data <- FetchData(object = seurat_object, 
                                                    vars = all_found_features[i], layer = layer)
                          max_exp_value <- max(feature_data)
                          min_exp_value <- min(feature_data)
                          if (seurat_version >= "5") {
                            single_plot <- suppressMessages(FeaturePlot(object = seurat_object, 
                                                                        features = all_found_features[i], order = order, 
                                                                        pt.size = pt.size, reduction = reduction, 
                                                                        raster = raster, split.by = split.by, raster.dpi = raster.dpi, 
                                                                        label = label, alpha = alpha_exp, ...) & 
                                                              scale_color_gradientn(colors = colors_use, 
                                                                                    limits = c(na_cutoff, max_exp_value), na.value = na_color, 
                                                                                    name = all_found_features[i])) & RestoreLegend() & 
                              theme(axis.title.y.right = element_blank())
                          }
                          else {
                            single_plot <- suppressMessages(FeaturePlot(object = seurat_object, 
                                                                        features = all_found_features[i], order = order, 
                                                                        pt.size = pt.size, reduction = reduction, 
                                                                        raster = raster, split.by = split.by, raster.dpi = raster.dpi, 
                                                                        label = label, ...) & scale_color_gradientn(colors = colors_use, 
                                                                                                                    limits = c(na_cutoff, max_exp_value), na.value = na_color, 
                                                                                                                    name = features[i])) & RestoreLegend() & 
                              theme(axis.title.y.right = element_blank())
                          }
                          if (isTRUE(x = label_feature_yaxis)) {
                            single_plot <- single_plot + plot_layout(nrow = num_rows, 
                                                                     ncol = num_columns)
                            single_plot <- single_plot & theme(legend.title = element_blank())
                            single_plot <- suppressMessages(single_plot + 
                                                              scale_y_continuous(sec.axis = dup_axis(name = all_found_features[i]))) + 
                              No_Right()
                          }
                          else {
                            single_plot <- single_plot + plot_layout(nrow = num_rows, 
                                                                     ncol = num_columns)
                          }
                        })
    plot <- wrap_plots(plot_list) + plot_layout(ncol = 1)
  }
  if (getOption(x = "scCustomize_warn_na_cutoff", default = TRUE) && 
      !is.na(x = na_cutoff) && na_cutoff == 1e-09) {
    cli_inform(message = c("", "NOTE: {.field FeaturePlot_scCustom} uses a specified {.code na_cutoff} when plotting to", 
                           "color cells with no expression as background color separate from color scale.", 
                           "Please ensure {.code na_cutoff} value is appropriate for feature being plotted.", 
                           "Default setting is appropriate for use when plotting from 'RNA' assay.\n", 
                           "When {.code na_cutoff} not appropriate (e.g., module scores) set to NULL to \n", 
                           "plot all cells in gradient color palette.", "", 
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_na_cutoff = FALSE)
  }
  if (getOption(x = "scCustomize_warn_zero_na_cutoff", default = TRUE) && 
      !is.na(x = na_cutoff) && na_cutoff == 0) {
    cli_inform(message = c("", "NOTE: Specified {.code na_cutoff} is set to {.field zero (0)}. This means that only cells/nuclei", 
                           "with expression less than zero will be plotted with {.code na_color}: {.val {na_color}}.", 
                           "To plot cells with expression values of zero using {.code na_color} leave", 
                           "default {.code na_cutoff} value. If you want to plot full spectrum without", 
                           "{.code na_cutoff} (e.g., for module scores) then set {.code na_cutoff = NULL`}.", 
                           "", "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_na_cutoff = FALSE)
  }
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot <- plot & theme(aspect.ratio = aspect_ratio)
  }
  if (isTRUE(x = figure_plot)) {
    if (length(x = all_found_features) == 1) {
      plot <- Figure_Plot(plot = plot)
    }
    else {
      plot_list <- lapply(1:length(x = all_found_features), 
                          function(j) {
                            fig_plot <- Figure_Plot(plot = plot[[j]])
                          })
      plot <- wrap_plots(plot_list, ncol = num_columns)
    }
  }
  return(plot)
}
