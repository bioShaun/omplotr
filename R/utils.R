RANDOM_SEED <- 1

save_ggplot <- function(ggplot_out, output,
                      width=8, height=6) {
  ggsave(paste(output, "png", sep = "."),
         plot = ggplot_out,
         width = width,
         height = height,
         dpi = 300, type = "cairo")
  ggsave(paste(output, "pdf", sep = "."),
         plot = ggplot_out,
         width = width,
         height = height,
         device = cairo_pdf)
}

save_general_plot <- function(plot, output,
                              width=8, height=6) {

  pdf(paste(output, "pdf", sep = "."),
      width = width, height = height, onefile = F)
  tmp_plot <- dev.cur()
  png(paste(output, "png", sep = "."),
      width = width, height = height, units = "in", res = 300)
  dev.control("enable")
  plot
  dev.copy(which=tmp_plot)
  dev.off()
  dev.off()
}


norm_exp_data <- function(exp_data) {
  # filter all zero exp values
  plot_data <- exp_data[rowSums(exp_data) > 0, ]

  # log normalization
  plot_data <- log2(plot_data + 1)

  return(plot_data)
}

get_go_gene <- function(go_id, gene_go_df) {
  gene_list <- gene_go_df[gene_go_df[, 2] == go_id, 1]
  paste(gene_list, collapse = ",")
}

get_go_test_data <- function(go_anno_file, gene_length_file,
                             test_gene_num=1000) {
  go_anno_df <- data.table::fread(go_anno_file, sep = ',')
  colnames(go_anno_df) <- c('gene_id', 'go_id')
  go_anno_df <- go_anno_df[go_anno_df$go_id !='', ]
  test_genes <- unique(go_anno_df$gene_id)[1:test_gene_num]
  test_go_anno <- dplyr::filter(go_anno_df, gene_id %in% test_genes)
  set.seed(RANDOM_SEED)
  test_diff_genes <- sample(test_genes, test_gene_num/10)
  gene_len_df <- data.table::fread(gene_length_file, header = F)
  colnames(gene_len_df) <- c('gene_id', 'gene_len')
  test_gene_len <- dplyr::filter(gene_len_df, gene_id %in% test_genes)

  go_test_data_list <- list(
    test_go_anno=test_go_anno,
    test_diff_genes=test_diff_genes,
    test_gene_len=test_gene_len
  )

  return(go_test_data_list)
}

