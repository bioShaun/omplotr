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


# functions for rnaseq plot
check_input <- function(file_path, err_message=NULL) {
  if ( is.null(file_path) || ! file.exists(file_path)) {
    if (is.null(err_message)) {
      err_message <- paste(file_path, 'Not exist!')
    }
    stop(err_message)
  }
}


valid_input <- function() {
  print('Input file is valid.')
}


select_exp_data <- function(exp_obj, item_file, select='row') {
  if (! is.na(item_file)) {
    check_input(item_file)
    item_df <- read.delim(item_file, header = F)
    item_num <- dim(item_df)[1]
    if (select == 'row') {
      exp_obj <- dplyr::filter(exp_obj, target_id %in% item_df$V1)
      item_name <- 'targets'
    } else if (select == 'column') {
      exp_obj <- dplyr::select(exp_obj, c('target_id', item_df$V1))
      item_name <- 'samples'
    } else {
      stop('Wrong select parameter [row, column].')
    }
    print(paste('Select', item_num, item_name))
  }
  return(exp_obj)
}


load_exp_file <- function(exp_file, genes, samples) {
  exp_obj <- data.table::fread(exp_file, check.names=F)
  total_target <- dim(exp_obj)[1]
  total_sample <- dim(exp_obj)[2] - 1
  print(paste('Total', total_target, 'targets,',
              total_sample, 'samples.'))
  colnames(exp_obj)[1] <- 'target_id'
  exp_obj <- select_exp_data(exp_obj, genes, 'row')
  exp_obj <- select_exp_data(exp_obj, samples, 'column')
  valid_input()
  exp_df <- data.frame(exp_obj, check.names=F)
  rownames(exp_df) <- exp_df$target_id
  exp_df <- exp_df[, -1]
  if (dim(exp_df)[1] > DIFF_HEATMAP_GENE) {
    gene_mean_exp <- sort(rowMeans(exp_df),
                          decreasing = T)
    top_genes <- names(gene_mean_exp[1:DIFF_HEATMAP_GENE])
    exp_df <- exp_df[top_genes, ]
  }

  return(exp_df)
}


label_sample <- function(exp_df, group_vs_sample, sample_list, group_mean_exp=F) {
  if (group_mean_exp | is.na(group_vs_sample)) {
    sample_inf <- data.frame(condition=colnames(exp_df),
                             sample=colnames(exp_df))
  } else {
    sample_inf <- read.delim(group_vs_sample, header=F)
    colnames(sample_inf) <- c("condition", "sample")
    sample_df <- read.delim(sample_list, header = F)
    sample_inf <- dplyr::filter(sample_inf, condition %in% sample_df$V1)
  }
  return(sample_inf)
}


test_data <- function(exp_df, test) {
  if (test) {
    return(head(exp_df, 1000))
  } else {
    return(exp_df)
  }
}


exp_by_group <- function(exp_df, group_vs_sample) {
  sample_inf <- label_sample(exp_df,
                             group_vs_sample)
  group_exp_df <- merge(sample_inf, t(exp_df),
                        by.x='sample',
                        by.y=0)
  total_col <- dim(group_exp_df)[2]
  group_mean_df <- aggregate(group_exp_df[, 3:total_col],
                             list(group_exp_df$condition), mean)
  row.names(group_mean_df) <- group_mean_df$Group.1
  group_mean_df <- group_mean_df[,-1]
  t_group_mean_df <- t(group_mean_df)
  t_group_mean_df <- t_group_mean_df[,unique(sample_inf$condition)]
  return(t_group_mean_df)
}


cluster_plot <- function(exp_df, sample_inf, out_prefix) {
  # cluster plot
  cluster_prefix <- paste(out_prefix, '.cluster.P', cluster_cut_tree, sep = '')
  cluster_data_dir <- paste(out_prefix, '.cluster.P', cluster_cut_tree,'.all', sep = '')
  dir.create(cluster_data_dir, showWarnings = F)
  diff_matrix <- as.matrix(exp_df)
  diff_gene_count <- dim(diff_matrix)[1]
  log_diff_matrix <- log2(diff_matrix + 1)
  # center rows, mean substracted
  scale_log_diff_matrix = t(scale(t(log_diff_matrix), scale = F))

  # gene clustering according to centered distance values.
  gene_dist = dist(scale_log_diff_matrix, method = "euclidean")
  hc_genes = hclust(gene_dist, method = "complete")
  gene_partition_assignments <- cutree(as.hclust(hc_genes), h = cluster_cut_tree/100 * max(hc_genes$height))

  max_cluster_count = max(gene_partition_assignments)
  cluster_num_cutoff = max(c(MIN_CLUSTER_NUM, diff_gene_count * MIN_CLUSTER_POR))

  all_partition_list <- list()
  m = 1
  for (i in 1:max_cluster_count) {
    partition_i = (gene_partition_assignments == i)
    partition_data = scale_log_diff_matrix[partition_i, , drop = F]
    cluster_name <- paste("cluster", i, sep = "_")
    partition_data_df <- as.data.frame(partition_data)
    partition_data_df <- cbind(Gene_id = rownames(partition_data_df), partition_data_df)
    write.table(partition_data_df, file = paste(cluster_data_dir, "/", cluster_name,
                                                ".txt", sep = ""), quote = F, row.names = F, sep = "\t")
    partition_data_df$cluster <- cluster_name
    melt_partition_data_df <- melt(partition_data_df, id = c("cluster", "Gene_id"))
    melt_partition_data_df$variable <- factor(melt_partition_data_df$variable,
                                              levels = sample_inf$sample)
    out_prefix <- file.path(cluster_data_dir, cluster_name)
    om_cluster_plot(melt_partition_data_df, out_prefix = out_prefix)
    if (dim(partition_data)[1] > cluster_num_cutoff) {
      all_partition_list[[m]] <- partition_data_df
      m <- m + 1
    }
  }

  all_cluster_df <- ldply(all_partition_list, data.frame)
  colnames(all_cluster_df) <- colnames(partition_data_df)
  melt_all_cluster_df <- melt(all_cluster_df, id = c("cluster", "Gene_id"))
  melt_all_cluster_df$cluster <- factor(melt_all_cluster_df$cluster, levels = unique(melt_all_cluster_df$cluster))
  om_cluster_plot(melt_all_cluster_df, out_prefix = cluster_prefix)
}
