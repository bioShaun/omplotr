#' plot expression heatmap
#' @param plot_data expression matrix
#' @param samples sample <-> group map, first column is group id
#' @param outdir out put directory, default is NULL, don't output file
#' @examples
#' om_heatmap(exp_test_data, test_sample_data)
om_heatmap <- function(exp_data, samples,
                       outdir=NULL, out_prefix=NULL) {

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
                      height = heatmap_heigh,
                      plot_type='pdf')
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh,
                      plot_type='png')
  } else if (! is.null(out_prefix)) {
    out_dir <- dirname(out_prefix)
    path_name <- save_mkdir(out_dir)
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh,
                      plot_type='pdf')
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh,
                      plot_type='png')
  }
  return(draw_heatmap())

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
  sample_num = length(unique(data.m$sample_id))

  theme_set(theme_onmath() +
              theme(axis.text.x = element_text(angle = 90,
                                               vjust = 0.5, hjust = 0),
                    legend.text = element_text(size = rel(0.8)),
                    legend.key = element_blank()))


  boxplot <- ggplot(data.m, aes(x = sample_id, y = exp_level, fill = group)) +
    geom_boxplot(notch = T) +
    guides(fill = guide_legend(nrow = 8, title = "group")) +
    scale_fill_manual(values = group_cols) +
    ylab("Log2(TPM)") + xlab('')
  boxplot_no_guide <- boxplot + guides(fill = F)

  violin <- ggplot(data.m, aes(x = sample_id, y = exp_level, fill = group)) +
    geom_violin() +
    guides(fill = guide_legend(nrow = 8, title = "group")) +
    scale_fill_manual(values = group_cols) +
    ylab("Log2(TPM)") + xlab('')
  violin_no_guide <- violin + guides(fill = F)

  density_plot <- ggplot(data.m, aes(exp_level, color = sample_id, fill = group)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = group_cols) +
    scale_color_manual(values = sample_cols) +
    theme(axis.text.x = element_text(angle = 0)) +
    guides(fill = guide_legend(nrow = 8, title = "group"),
           color = guide_legend(nrow = 8, title = "sample")) +
    xlab('Log2(TPM)') +  ylab('Density')


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

#' plot diff expression volcano plot
#' @param dpa_results diff analysis results
#' @param compare_name compare name of diff analysis, 'ALL' for merged results
#' @param logfc log2foldchange cutoff for DEG analysis, default is 1
#' @param qvalue adjust pvalue cutoff for DEG analysis, default is 0.05
#' @param outdir out put directory, default is NULL, don't output file
#' @examples
#' om_volcano_plot(diff_test_data, 'Case_vs_Control')
#' om_volcano_plot(diff_test_data, 'ALL')
om_volcano_plot <- function(dpa_results, compare_name,
                            logfc=1, qvalue=0.05, outdir=NULL) {

  dpa_results$logFDR <- -log10(dpa_results$FDR)
  dpa_results$color <- "blue"
  up_name = unlist(strsplit(compare_name, split = "_vs_"))[1]
  down_name = unlist(strsplit(compare_name, split = "_vs_"))[2]
  for (i in 1:dim(dpa_results)[1]) {
    if (dpa_results$logFC[i] > logfc & dpa_results$FDR[i] < qvalue)
      dpa_results$color[i] <- "red" else if (dpa_results$logFC[i] < -(logfc) & dpa_results$FDR[i] < qvalue)
        dpa_results$color[i] <- "green"
  }

  # labels and titles
  dpa_results2 <- dpa_results
  count <- table(dpa_results2$color)
  count <- as.data.frame(count)

  red_count <- sum(count$Freq[which(count$Var1 == "red")])
  green_count <- sum(count$Freq[which(count$Var1 == "green")])
  all_count_number <- red_count + green_count
  # all_count <- paste('Differential Expressed Genes',all_count_number,sep = ':')

  if (all_count_number < 100) {
    dpa_results2$logFC <- ifelse(dpa_results2$logFC > 15,
                                 15 + (dpa_results2$logFC -15)/logFC_max,
                                 dpa_results2$logFC)
    dpa_results2$logFC <- ifelse(dpa_results2$logFC < -15,
                                 -15 - (dpa_results2$logFC +15)/logFC_min,
                                 dpa_results2$logFC)
    logFDR_max = max(dpa_results2$logFDR)
    dpa_results2$logFDR <- ifelse(dpa_results2$logFDR > 50,
                                  50 + (dpa_results2$logFDR -50)/logFDR_max,
                                  dpa_results2$logFDR)
  }

  logFC_max = max(dpa_results2$logFC)
  logFC_min = min(dpa_results2$logFC)
  logFC_limit <- ceiling(max(c(abs(logFC_min), logFC_max)))
  logFC_limit <- ifelse(logFC_limit < 8, 8, logFC_limit)
  logFC_limit <- ifelse(logFC_limit > 15, 15, logFC_limit)
  logFDR_limit <- ceiling(max(dpa_results2$logFDR))
  logFDR_limit <- ifelse(logFDR_limit > 50, 50, logFDR_limit)
  logFDR_limit <- ifelse(logFDR_limit < 35, 35, logFDR_limit)

  red_label <- paste("No.", up_name, "up-regulated genes:", red_count, sep = " ")
  green_label <- paste("No.", down_name, "up-regulated genes:", green_count, sep = " ")
  red <- RColorBrewer::brewer.pal(6, "Reds")[6]
  green <- RColorBrewer::brewer.pal(6, "Greens")[6]
  blue <- RColorBrewer::brewer.pal(4, "Blues")[4]

  theme_set(theme_onmath() +
              theme(legend.key = element_blank(),
                    panel.grid.major = element_blank(),
                    legend.position = "bottom",
                    legend.direction = "vertical"))


  y_line_pos = round(-log10(qvalue), 1)
  p <- ggplot(dpa_results2, aes(logFC, logFDR, colour = color)) +
    geom_point(size = 0.6) +
    geom_hline(yintercept = y_line_pos, lty = 4, size = 0.45) +
    geom_vline(xintercept = -(logfc),lty = 4, size = 0.45) +
    geom_vline(xintercept = logfc, lty = 4, size = 0.45) +
    xlab("logFC") + ylab("-log10(FDR)")

  compare_number = 0
  if (compare_name == 'ALL') {
    compare_number <- length(unique(dpa_results2$compare))
    facet_wrap_ncol = round(sqrt(compare_number))
    p <- p + guides(color = F) +
      scale_color_manual(values = c(red = red, green = green, blue = blue)) +
      facet_wrap(~compare, ncol = facet_wrap_ncol)
  } else {
    p <- p + guides(color = guide_legend(title = "")) +
      scale_color_manual(values = c(red = red,green = green, blue = blue),
                         breaks = c("green", "red"), labels = c(green_label,red_label)) +
      ggtitle(compare_name) + scale_y_continuous(breaks = c(0, 10, 20, 30),
                                                 limits = c(0, logFDR_limit)) +
      scale_x_continuous(breaks = c(-8,-4, -2, -1, 0, 1, 2, 4, 8),
                         limits = c(-logFC_limit, logFC_limit))

  }

  plot_height <- 8 + compare_number/4
  plot_width <- 6 + compare_number/4
  if (! is.null(outdir)) {
    out_prefix <- file.path(outdir, paste(compare_name, 'Volcano_plot', sep = '.'))
    save_ggplot(p, out_prefix, width = plot_width, height = plot_height)
  }
  return(p)

}

#' plot PCA point plot
#' @param exp_data expression matrix
#' @param samples sample <-> group map, first column is group id
#' @param outdir out put directory, default is NULL, don't output file
#' @examples
#' om_pca_plot(exp_test_data, test_sample_data)
om_pca_plot <- function(exp_data, samples, outdir=NULL,
                        out_prefix=NULL) {

  plot_data <- norm_exp_data(exp_data)
  PCA_data_mat <- t(apply(plot_data[, 1:dim(plot_data)[2]], 2, as.numeric))
  PCA <- prcomp(PCA_data_mat)
  Summary_PCA <- summary(PCA)

  PCA_data <- as.data.frame(PCA$x[, 1:dim(samples)[1]])
  sample_name <- rownames(PCA$x)
  match_index <- match(sample_name, samples$sample, nomatch = 0)
  group_name <- samples$condition[match_index]
  PCA_data$Sample <- sample_name
  PCA_data$Group <- group_name

  sample_types <- as.character(unique(samples$condition))
  nsamples <- length(sample_types)
  sample_colors <- colorRampPalette(heatmap_col)(nsamples)

  pca_plot <- ggplot(PCA_data, aes(PC1, PC2)) +
      geom_point(aes(colour = Group), size = rel(3)) +
      geom_text(aes(label = Sample), vjust = 0,
                hjust = 0.5, color = "black",
                size = rel(2)) +
      geom_vline(xintercept = 0, linetype = 2,
                 colour = "grey60", size = rel(0.5)) +
      geom_hline(yintercept = 0, linetype = 2,
                 colour = "grey60", size = rel(0.5)) +
      theme_onmath() +
      scale_color_manual(values = sample_colors) +
      labs(title = "PCA",
           x = paste0("PC1: ", Summary_PCA$importance[2, 1] * 100, "% variance"),
           y = paste0("PC2: ", Summary_PCA$importance[2, 2] * 100, "% variance"))

  if (! is.null(outdir)) {
    out_prefix <- file.path(outdir, "PCA_plot")
    save_ggplot(pca_plot, out_prefix)
  } else if (! is.null(out_prefix)) {
    out_dir <- dirname(out_prefix)
    path_name <- save_mkdir(out_dir)
    save_ggplot(ggplot_out = pca_plot,
                output = out_prefix)
  }

  return(pca_plot)
}


#' plot correlation heatmap
#' @param exp_data expression matrix
#' @param samples sample <-> group map, first column is group id
#' @param outdir out put directory, default is NULL, don't output file
#' @examples
#' om_correlation_plot(exp_test_data, test_sample_data)
om_correlation_plot <- function(exp_data, samples, outdir=NULL) {
  data <- norm_exp_data(exp_data)
  sample_types <- as.character(unique(samples$condition))
  rep_names <- as.character(samples$sample)
  data <- data[, colnames(data) %in% samples$sample, drop = F]
  nsamples <- length(sample_types)
  sample_colors <- colorRampPalette(heatmap_col)(nsamples)
  data <- as.matrix(data)
  sample_cor <- cor(data, method = "pearson", use = "pairwise.complete.obs")
  sample_cor_df <- as.data.frame(sample_cor)
  sample_cor_df <- cbind(Sample = rownames(sample_cor_df), sample_cor_df)

  Group <- sample_colors[1:length(unique(samples$condition))]
  names(Group) <- unique(samples$condition)
  ann_color = data.frame(group = samples$condition)
  rownames(ann_color) <- samples$sample
  ann_colors = list(group = Group)

  theme_set(theme_onmath() + theme(legend.position = c(0.5, 0.5)))
  sample_num = length(colnames(exp_data))
  heatmap_width <- (sample_num - 5)/5 + 7
  heatmap_heigh <- (sample_num - 5)/5 + 6
  fontsize = (sample_num - 5)/10 + 7
  cellwidth <- (heatmap_width - 1) * 50/sample_num

  draw_heatmap <- function(){
    pheatmap::  pheatmap(sample_cor, annotation_col = ann_color,
                         annotation_colors = ann_colors,
                         annotation_row = ann_color, annotation_names_row = F,
                         annotation_names_col = F,
                         color = rev(cor_plot_col), border_color = NA,
                         cellwidth = cellwidth, fontsize = fontsize)
  }

  if (! is.null(outdir)) {

    out_prefix <- file.path(outdir, 'Sample.correlation.heatmap')
    write.table(sample_cor_df, file = paste(out_prefix, "txt",sep = "."),
                quote = F, sep = "\t", row.names = F)
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh,
                      plot_type='pdf')
    save_general_plot(draw_heatmap(), out_prefix,
                      width = heatmap_width,
                      height = heatmap_heigh,
                      plot_type='png')

  }

  return(draw_heatmap())

}

#' plot correlation heatmap
#' @param plot_data cluster result
#' @param out_prefix out put prefix, default is NULL, don't output file
#' @examples
#' om_cluster_plot(cluster_test_data)
om_cluster_plot <- function(plot_data, out_prefix=NULL) {

  cluster_number <- length(unique(plot_data$cluster))
  col_theme <- colorRampPalette(heatmap_col)(cluster_number)
  theme_cluster <- theme_onmath() + theme(axis.text.x = element_text(vjust = -0.2,
                                                                     angle = 90,
                                                                     size = rel(0.8)))

  cluster_plot <- ggplot(plot_data, aes(x=variable, y=value, group = Gene_id, color=cluster)) +
    geom_line(alpha = 0.2) +
    scale_color_manual(guide=F, values = col_theme) +
    xlab("") + ylab("Scaled log10(tpm+1)") + theme_cluster

  if (cluster_number > 1) {
    facet_wrap_ncol = round(sqrt(cluster_number))
    cluster_plot <- cluster_plot + facet_wrap(~cluster, ncol = facet_wrap_ncol)
  }

  if (! is.null(out_prefix)) {
    plot_height <- 6 + cluster_number/4
    plot_width <- 8 + cluster_number/4
    save_ggplot(cluster_plot, out_prefix,
                width = plot_width, height = plot_height)
  }

  return(cluster_plot)

}
