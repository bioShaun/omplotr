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

