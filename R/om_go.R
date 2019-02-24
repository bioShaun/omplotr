#' goseq analysis
#' @param diff_genes diff gene vector
#' @param gene_length_df gene length dataframe
#' @param gene_go_map gene go map file
#' @param out_prefix out putfile prefix, default is NULL, don't output file
#' @examples
#' diff_genes <- go_test_data_list[['test_diff_genes']]
#' gene_length_df <- go_test_data_list[['test_gene_len']]
#' gene_go_map <- system.file("extdata", "topgo_test_data.txt", package = "omplotr")
#' goseq_output <- om_goseq(diff_genes, gene_length_df, gene_go_map)
#' head(goseq_output, 4)
om_goseq <- function(diff_genes, gene_length_df,
                     gene_go_map, out_prefix=NULL) {
  go_anno_df <- gene_map_to_go_anno(gene_go_map)
  all_id <- gene_length_df[, 1]
  gene.vector = as.integer(all_id %in% diff_genes)
  names(gene.vector) = all_id
  id_len <- gene_length_df[, 2]
  names(id_len) = all_id
  ## goseq
  pwf = goseq::nullp(gene.vector, bias.data = id_len, plot.fit = F)
  GO.wall = goseq::goseq(pwf, gene2cat = go_anno_df)
  GO.wall <- GO.wall[GO.wall$numDEInCat > 0, c(1, 2, 4, 5, 6, 7)]
  GO.wall$qvalue <- p.adjust(GO.wall$over_represented_pvalue,
                             method = "BH",
                             n = length(GO.wall$over_represented_pvalue))
  out_go <- GO.wall[, c(1, 2, 7, 3, 4, 5, 6)]
  out_go <- na.omit(out_go)
  ## add diff gene id to enrich table
  diff_go_anno_df <- go_anno_df[go_anno_df[, 1] %in% diff_genes, ]
  diff_go_anno_df <- diff_go_anno_df[diff_go_anno_df[, 1] != "", ]
  out_go_de_id <- unlist(lapply(out_go[, 1], get_go_gene,
                                gene_go_df = diff_go_anno_df))
  out_go$DE_id <- out_go_de_id

  if (! is.null(out_prefix)) {
    if (dim(out_go)[1] > 0) {
      write.table(out_go, file = paste(out_prefix, "txt", sep = "."),
                  quote = F, sep = "\t", row.names = F)
      # write.xlsx(out_go, file = paste(out_prefix, "xlsx", sep = "."),
      #            sheetName = "go.enrichment",
      #            append = FALSE, row.names = F)
    } else {
      print("No gene successfully annotated!")
    }
  }

  return(out_go)
}

#' plot enrichment barplot
#' @param enrich_df GO/KEGG enrichment result table
#' @param term_number term number to show, default is 30
#' @param ylab_title text of ylab, default is "-log10(qvalue)"
#' @param out_prefix output plot file prefix
#'
#' @examples
#' # kegg enrichment plot
#' kegg_enrich_file <- system.file("extdata", "enrichment", "example.kegg.enrichment.txt", package = "omplotr")
#' kegg_enrich_list <- clean_enrich_table(kegg_enrich_file)
#' head(kegg_enrich_list$table, 4)
#' om_enrich_bar_plot(kegg_enrich_list$table, ylab_title=kegg_enrich_list$title)
#' # go enrichment plot
#' go_enrich_file <- system.file("extdata", "enrichment", "example.go.enrichment.txt", package = "omplotr")
#' go_enrich_list <- clean_enrich_table(go_enrich_file)
#' om_enrich_bar_plot(go_enrich_list$table, ylab_title=go_enrich_list$title)
om_enrich_bar_plot <- function(enrich_df,
                               term_number=30,
                               ylab_title='-log10(qvalue)',
                               out_prefix=NULL) {
  enrich_df <- dplyr::filter(enrich_df, qvalue < 1)
  if (dim(enrich_df)[1] > 0) {
    enrich_df$log10qvalue <- -log10(enrich_df$qvalue)
    enrich_df$wrap_term <- sapply(enrich_df$term, wrap_long_name)
    enrich_df$wrap_term <- factor(enrich_df$wrap_term,
                                  levels = rev(enrich_df$wrap_term))
    top_df <- head(enrich_df, term_number)

    plot_witdh <- 6
    plot_height <- dim(top_df)[1] / 3

    p <- ggplot(top_df, aes(wrap_term, log10qvalue, fill = ontology, alpha = log10qvalue)) +
      geom_bar(stat = 'identity', width = 0.45) +
      coord_flip() +
      facet_grid(ontology~., scales = 'free_y',
                 space = 'free_y') +
      theme_onmath() +
      theme(axis.text.y = element_text(color = 'grey30', face = "plain",
                                       size = rel(0.75)),
            strip.text.y = element_text(angle = 0, colour = 'grey30')) +
      scale_fill_brewer(palette = 'Set1') +
      guides(fill = F, alpha = F) +
      xlab('') +
      ylab(ylab_title)
    if (! is.null(out_prefix)) {
      save_ggplot(p, out_prefix,
                  width = plot_witdh,
                  height = plot_height)
    }
    return(p)
  } else {
    return('Nothing to plot!')
  }
}

#' topGO plot
#' @param gene_go_map gene go map file
#' @param diff_genes diff gene vector
#' @param enrich_result_df goseq analysis result
#' @param name out put file name, default is NULL, don't output file
#' @param out_dir out put directory, default is NULL, don't output file
#' @examples
#' gene_go_map <- system.file("extdata", "topgo_test_data.txt", package = "omplotr")
#' diff_genes <- go_test_data_list[['test_diff_genes']]
#' gene_length_df <- go_test_data_list[['test_gene_len']]
#' enrich_result_df <- om_goseq(diff_genes, gene_length_df, gene_go_map)
#' om_topgo(gene_go_map, diff_genes, enrich_result_df)
om_topgo <- function(gene_go_map, diff_genes, enrich_result_df,
                      name=NULL, out_dir=NULL) {
  geneID2GO <- readMappings(file = gene_go_map)
  geneNames <- names(geneID2GO)
  geneList <- factor(as.integer(geneNames %in% diff_genes))
  names(geneList) <- geneNames
  go_catogary_vector <- c("MF", "CC", "BP")
  enrich_result_df <- dplyr::filter(enrich_result_df, numInCat >= 5)

  topgo_data <- list(
    plot_data <- list(),
    node_size <- list()
  )

  for (i in 1:length(go_catogary_vector)) {
    go_catogary <- go_catogary_vector[i]
    each_enrich_result <- dplyr::filter(enrich_result_df,
                                        ontology == go_catogary)
    each_go_qvalue <- each_enrich_result$qvalue
    each_go_qvalue[which(each_go_qvalue == 0)] <- 1e-100
    names(each_go_qvalue) <- each_enrich_result[, 1]
    if (dim(each_enrich_result)[1] < 2) {
      out_info <- paste("Too little gene annotated to ",
                        go_catogary, sep = "")
      print(out_info)
    } else {
      GOdata <- new("topGOdata", ontology = go_catogary,
                    allGenes = geneList,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO, nodeSize = 5)
      if (dim(each_enrich_result)[1] <= 10) {
        topgo_plot_fun <- function() {
          showSigOfNodes(GOdata, each_go_qvalue,
                         firstSigNodes = 1, useInfo = "all")
        }

      } else {
        topgo_plot_fun <- function() {
          showSigOfNodes(GOdata, each_go_qvalue,
                         firstSigNodes = 5, useInfo = "all")
        }
      }
      if (is.null(name) || is.null(out_dir)) {
        topgo_plot_fun()
      } else {
        out_prefix <- file.path(out_dir,
                                paste(name, go_catogary, 'GO.DAG', sep = '.'))
        save_general_plot(topgo_plot_fun(), out_prefix,
                          width = 8, height = 8,
                          plot_type='pdf')
        save_general_plot(topgo_plot_fun(), out_prefix,
                          width = 8, height = 8,
                          plot_type='png')
      }
    }
  }
}

