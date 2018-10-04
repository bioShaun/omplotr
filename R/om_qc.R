#' GC plot function
#' @param plot_data gc dataframe
#' @param output out put prefix, default is NULL, don't output file
#' @examples
#' gc_line_plot(gc_test_data)

gc_line_plot <- function(plot_data, output=NULL) {

  sample_number <- length(unique(plot_data$sample))
  seq_len <- round(max(plot_data[, 1])/2)
  max_gc <- max(plot_data$value) + 0.1

  gc_plot <- ggplot(plot_data, aes(x = X.Base, y = value, colour = variable)) +
    geom_line() +
    geom_vline(xintercept = seq_len, linetype = 2) +
    scale_x_continuous(breaks = seq(from = 0, to = 2 * seq_len, by = seq_len),
                       labels = seq(from = 0, to = 2 * seq_len, by = seq_len)) +
    scale_y_continuous(breaks = seq(0, max_gc, by = round((max_gc)/4, 1)),
                       labels = scales::percent(seq(0, max_gc, by = round((max_gc)/4, 1)))) +
    xlab("Postion") + ylab("Percent(%)") +
    guides(color = guide_legend(title = "")) +
    theme_onmath() +
    scale_color_brewer(palette = 'Set1')

  if (sample_number > 1) {
    facet_wrap_ncol = round(sqrt(sample_number))
    gc_plot <- gc_plot + facet_wrap(~sample, ncol = facet_wrap_ncol)
  }

  if (! is.null(output)) {
    plot_height <- 6 + sample_number/4
    plot_width <- 8 + sample_number/4
    save_ggplot(gc_plot, output,
              width=plot_width,
              height=plot_height)
  }

  return(gc_plot)

}

#' Reads Quality barplot function
#' @param plot_data reads quality dataframe
#' @param output out put prefix, default is NULL, don't output file
#' @examples
#' reads_quality_plot(rq_test_data)
reads_quality_plot <- function(plot_data, output=NULL) {

  col <- plot_data$color
  names(col) <- col
  sample_number <- length(unique(plot_data$sample))

  p <- ggplot(plot_data, aes(x = Quality, y = Proportion, fill = color)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col) +
    geom_segment(aes(x = 30, y = 0, xend = 30,yend = max(plot_data$Proportion)),
                 colour = "red", linetype = "dashed", size = 1) +
    theme_onmath() +
    guides(fill = F) +
    scale_y_continuous(breaks = seq(from = 0,to = max(plot_data$Proportion), by = 0.1),
                       labels = scales::percent(seq(from = 0,to = max(plot_data$Proportion), by = 0.1))) +
    xlab("Quality Score")

  if (sample_number > 1) {
    facet_wrap_ncol = round(sqrt(sample_number))
    p <- p + facet_wrap(~sample, ncol = facet_wrap_ncol)
  }



  if (! is.null(output)) {
    plot_height <- 6 + sample_number/4
    plot_width <- 6 + sample_number/4
    save_ggplot(p, output,
              width=plot_width,
              height=plot_height)
  }

  return(p)

}
