#' Principal component visualization
#'
#' @docType methods
#' @name pcView
#' @rdname pcView
#'
#' @param mat A data matrix.
#' @param color The column name specifying the color.
#' @param filename The file name of output figure.
#' @param width The width of the output figure.
#' @param height The height of the output figure.
#' @param ... parameters in ggsave.
#'
#' @return a ggplot instance.
#' @import ggplot2 ggpubr
#' @author Wubing Zhang
#' @export
#'
pcView <- function(mat, color = gsub(".*_", "", colnames(mat)),
                   filename = NULL, width = 5, height = 4, ...){
  requireNamespace("ggplot2")
  requireNamespace("ggpubr")
  mat = as.matrix(mat)
  tmp = prcomp(t(mat))
  gg = as.data.frame(tmp$x[,1:2], stringsAsFactors = FALSE)
  gg$Color = color
  p = ggplot(gg, aes_string("PC1", "PC2", color = "Color"))
  p = p + geom_point()
  p = p + labs(color = NULL)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))
  if(!is.null(filename))
    ggsave(filename, p, width = width, height = height, ...)
  return(p)
}
