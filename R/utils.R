save_plot <- function(ggplot_out, output,
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
