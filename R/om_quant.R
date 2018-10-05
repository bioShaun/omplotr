#' plot expression heatmap
#' @param plot_data expression matrix
#' @param samples sample <-> group map, first column is group id
#' @param outdir out put directory, default is NULL, don't output file
#' @examples
#' om_heatmap(exp_test_data, test_sample_data)
om_heatmap <- function(exp_data, samples, outdir=NULL) {

  # normalize exp data
  plot_data <- norm_exp_data(exp_data)

  # heatmap outer bar color
  if (length(unique(samples$condition)) > length(heatmap_col)) {
    heatmap_col <- colorRampPalette(heatmap_col)(length(unique(samples$condition)))
  }

  # heatmap cell color
  cell_cols <- RColorBrewer::brewer.pal(10,"RdYlGn")
  gradiant_cell_cols <- rev(colorRampPalette(cell_cols)(100))

  # map group to bar color
  Group <- heatmap_col[1:length(unique(samples$condition))]
  names(Group) <- unique(samples$condition)
  ann_color = data.frame(group = samples$condition)
  rownames(ann_color) <- samples$sample
  ann_colors = list(group = Group)

  # adjust output height & width
  sample_num = length(colnames(plot_data))
  heatmap_width <- (sample_num - 5)/3 + 2
  heatmap_heigh <- (sample_num - 5)/3 + 4
  fontsize = (sample_num - 5)/10 + 4.5
  cellwidth <- (heatmap_width - 0.5) * 50/sample_num

  ## TODO legend size
  draw_heatmap <- function(){
    pheatmap::pheatmap(plot_data, show_rownames = F,
                       annotation_col = ann_color, annotation_colors = ann_colors,
                       annotation_legend = T, annotation_names_col = F,
                       color = gradiant_cell_cols,
                       treeheight_row = 0, scale = "row", fontsize = fontsize,
                       cellwidth = cellwidth, border_color = NA)
  }

  if (! is.null(outdir)) {
    out_prefix <- file.path(outdir, 'Diff.genes.heatmap')
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh)
  } else {
    return(draw_heatmap())
  }

}


#' plot expression box, violin, density plot or merged plots of all above
#' @param plot_data expression matrix
#' @param samples sample <-> group map, first column is group id
#' @param plot_type type of plot to choose: box, violin, density, all
#' @param outdir out put directory, default is NULL, don't output file
#' @examples
#' om_boxplot(exp_test_data, test_sample_data, 'box')
#' om_boxplot(exp_test_data, test_sample_data, 'violin')
#' om_boxplot(exp_test_data, test_sample_data, 'density')
#' om_boxplot(exp_test_data, test_sample_data, 'all')
om_boxplot <- function(exp_data, samples,
                       plot_type, outdir=NULL) {

  # normalize data
  plot_data <- norm_exp_data(exp_data)

  # reshape matrix and add group inf
  data.m <- reshape2::melt(as.matrix(plot_data))
  colnames(data.m) <- c('gene_id', 'sample_id', 'exp_level')
  row.names(samples) <- samples$sample
  data.m$sample_id <- as.character(data.m$sample_id)
  data.m$group <- samples[data.m$sample_id, 'condition']
  data.m$sample_id <- factor(data.m$sample_id, levels = samples$sample)
  data.m$exp_level <- as.numeric(data.m$exp_level)

  # get sample & group colors
  group_cols <- colorRampPalette(heatmap_col)(length(unique(data.m$group)))
  sample_cols <- colorRampPalette(sample_cols)(length(samples$sample))
  sample_num = length(unique(data.m$variable))

  theme_set(theme_onmath() +
              theme(axis.text.x = element_text(angle = -90, color = "black",
                                               vjust = 0.5, hjust = 0, size = rel(1.2)),
                    axis.text.y = element_text(size = rel(1.2)),
                    legend.text = element_text(size = rel(0.8)),
                    legend.key = element_blank(),
                    axis.title.x = element_blank()))


  boxplot <- ggplot(data.m, aes(x = sample_id, y = exp_level, fill = group)) +
    geom_boxplot(notch = T) +
    guides(fill = guide_legend(nrow = 8, title = "group")) +
    scale_fill_manual(values = group_cols) +
    ylab("Log2(TPM)")
  boxplot_no_guide <- boxplot + guides(fill = F)

  violin <- ggplot(data.m, aes(x = sample_id, y = exp_level, fill = group)) +
    geom_violin() +
    guides(fill = guide_legend(nrow = 8, title = "group")) +
    scale_fill_manual(values = group_cols) +
    ylab("Log2(TPM)")
  violin_no_guide <- violin + guides(fill = F)

  density_plot <- ggplot(data.m, aes(exp_level, color = sample_id, fill = group)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = group_cols) +
    scale_color_manual(values = sample_cols) +
    theme(axis.text.x = element_text(angle = 0)) +
    guides(fill = guide_legend(nrow = 8, title = "group"),
           color = guide_legend(nrow = 8, title = "sample")) +
    ylab('Log2(TPM)')

  merged_plot <- function() {
    gridExtra::grid.arrange(boxplot_no_guide, violin_no_guide, density_plot,
                              nrow = 2, layout_matrix = rbind(c(1,2), c(3, 3)))
  }

  ## save plots to list
  plots_list <- list(
    box=boxplot,
    violin=violin,
    density=density_plot
  )

  if (plot_type == 'all') {
    if (! is.null(outdir)) {
      out_prefix = file.path(outdir, 'Gene_expression')
      out_plot <- merged_plot()
      merge_width = 8 + sample_num/4
      merge_height = 12 + sample_num/8
      save_ggplot(out_plot, out_prefix,
                  width = merge_width,
                  height = merge_height)
    } else {
      return(merged_plot())
    }

  } else {
    out_plot <- plots_list[[plot_type]]
    if (! is.null(outdir)) {
      out_prefix = file.path(outdir,
                             paste('Gene_expression', plot_type, sep = '.'))
      plot_width = 6 + sample_num/4
      plot_height = 6 + sample_num/8
      save_ggplot(out_plot, out_prefix,
                  width = plot_width,
                  height = plot_height)
    } else {
      return(out_plot)
    }
  }

}

