#' goseq analysis
#' @param diff_genes diff gene vector
#' @param gene_length_df gene length dataframe
#' @param go_anno_df go annotation file, seperated with ","
#' @param out_prefix out putfile prefix, default is NULL, don't output file
#' @examples
#' diff_genes <- go_test_data_list[['test_diff_genes']]
#' gene_length_df <- go_test_data_list[['test_gene_len']]
#' go_anno_df <- go_test_data_list[['test_go_anno']]
#' goseq_output <- run_goseq(diff_genes, gene_length_df, go_anno_df)
#' head(goseq_output, 4)
om_goseq <- function(diff_genes, gene_length_df, go_anno_df, out_prefix=NULL) {
  all_id <- gene_length_df[, 1]
  gene.vector = as.integer(all_id %in% diff_genes)
  names(gene.vector) = all_id
  id_len <- gene_length_df[, 2]
  names(id_len) = all_id
  ## goseq
  pwf = goseq::nullp(gene.vector, bias.data = id_len)
  GO.wall = goseq::goseq(pwf, gene2cat = go_anno_df)
  GO.wall <- GO.wall[GO.wall$numDEInCat > 0, c(1, 2, 4, 5, 6, 7)]
  GO.wall$qvalue <- p.adjust(GO.wall$over_represented_pvalue, method = "BH", n = length(GO.wall$over_represented_pvalue))
  out_go <- GO.wall[, c(1, 2, 7, 3, 4, 5, 6)]
  out_go <- na.omit(out_go)
  ## add diff gene id to enrich table
  diff_go_anno_df <- go_anno_df[go_anno_df[, 1] %in% diff_genes, ]
  diff_go_anno_df <- diff_go_anno_df[diff_go_anno_df[, 1] != "", ]
  out_go_de_id <- unlist(lapply(out_go[, 1], get_go_gene, gene_go_df = diff_go_anno_df))
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


# om_draw_topgo <- function(topgo_data) {
#
# }

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
#' go_anno_df <- go_test_data_list[['test_go_anno']]
#' enrich_result_df <- run_goseq(diff_genes, gene_length_df, go_anno_df)
#' run_topgo(gene_go_map, diff_genes, enrich_result_df)
run_topgo <- function(gene_go_map, diff_genes, enrich_result_df,
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
    each_enrich_result <- dplyr::filter(enrich_result_df, ontology == go_catogary)
    each_go_qvalue <- each_enrich_result$qvalue
    each_go_qvalue[which(each_go_qvalue == 0)] <- 1e-100
    names(each_go_qvalue) <- each_enrich_result[, 1]
    if (dim(each_enrich_result)[1] < 2) {
      out_info <- paste("Too little gene annotated to ", go_catogary, sep = "")
      print(out_info)
    } else {
      GOdata <- new("topGOdata", ontology = go_catogary, allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
      # topgo_data[['plot_data']][[i]] <- GOdata
      if (dim(each_enrich_result)[1] <= 10) {
        # topgo_plot_fun <- function() {
        #   showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = "all")
        # }
        # topgo_plots[[i]] <- showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = "all")
      #   pdf(file = paste(out_dir, "/", name, ".", go_catogary, ".GO.DAG.pdf",
      #                    sep = ""), width = 8, height = 8)
      #   topGO::showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = "all")
      #   dev.off()
      #   Cairo(file = paste(out_dir, "/", name, ".", go_catogary, ".GO.DAG.png",
      #                      sep = ""), type = "png", units = "in", width = 8, height = 8, pointsize = 12,
      #         dpi = 300, bg = "white")
      #   topGO::showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 1, useInfo = "all")
      #   dev.off()
      } else {
        showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 5, useInfo = "all")
        # pdf(file = paste(out_dir, "/", name, ".", go_catogary, ".GO.DAG.pdf",
        #                  sep = ""), width = 8, height = 8)
        # showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 5, useInfo = "all")
        # dev.off()
        # Cairo(file = paste(out_dir, "/", name, ".", go_catogary, ".GO.DAG.png",
        #                    sep = ""), type = "png", units = "in", width = 8, height = 8, pointsize = 12,
        #       dpi = 300, bg = "white")
        # showSigOfNodes(GOdata, each_go_qvalue, firstSigNodes = 5, useInfo = "all")
        # dev.off()
      }
    }
  }
}

