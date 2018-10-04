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
