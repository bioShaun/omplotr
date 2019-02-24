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
                              width=8, height=6,
                              plot_type='pdf') {
  if (plot_type == 'pdf') {
    pdf(paste(output, "pdf", sep = "."),
        width = width, height = height, onefile = F)
    plot
    dev.off()
  } else {
    png(paste(output, "png", sep = "."),
        width = width, height = height,
        units = "in", res = 300)
    plot
    dev.off()
  }
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


flat_gene_go_map <- function(each_row) {
  go_vector <- unlist(strsplit(as.character(each_row[2]), ","))
  gene_vector <- rep(as.character(each_row[1]), length(go_vector))
  go_df <- data.frame(gene_id=gene_vector, go_id=go_vector)
  return(go_df)
}

gene_map_to_go_anno <- function(gene_go_map) {
  gene_go_map_df <- data.table::fread(gene_go_map, header = F)
  gene_go_map_df <- data.frame(gene_go_map_df)
  go_anno_list <- lapply(data.frame(t(gene_go_map_df)), flat_gene_go_map)
  go_anno_df <- plyr::ldply(go_anno_list, data.frame)
  go_anno_df <- go_anno_df[, -1]
  return(go_anno_df)
}


save_mkdir <- function(dir_path) {
  if (! dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

palette_colors <- function(pal_name, sample_num) {
  col_pal_inf <- RColorBrewer::brewer.pal.info
  if (! pal_name %in% rownames(col_pal_inf)) {
    print('Palette for analysis:')
    print(rownames(col_pal_inf))
    stop('Wrong palette name!')
  }
  col_num <- col_pal_inf[pal_name, 'maxcolors']
  if (sample_num <= col_num) {
    return(RColorBrewer::brewer.pal(sample_num, pal_name))
  } else {
    return(colorRampPalette(RColorBrewer::brewer.pal(col_num, pal_name))(sample_num))
  }
}


str_percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

wrap_long_name <- function(name, width=30) {
  return(paste(strwrap(name, width = width), collapse="\n"))
}

clean_enrich_table <- function(enrich_file) {
  enrich_file_name <- basename(enrich_file)
  enrich_df <- read.delim(enrich_file)
  if (stringr::str_detect(enrich_file_name, 'go.enrichment')) {
    enrich_df <- enrich_df[, c('qvalue', 'term', 'ontology')]
    ylab_title <- '-log10(qvalue)'

  } else if (stringr::str_detect(enrich_file_name, 'kegg.enrichment')) {
    enrich_df <- enrich_df[, c('Corrected.P.Value', 'X.Term', 'Database')]
    colnames(enrich_df) <- c('qvalue', 'term', 'ontology')
    ylab_title <- '-log10(Corrected.P.Value)'
  }
  enrich_plot_data_list <- list(table=enrich_df, title=ylab_title)
  return(enrich_plot_data_list)
}
