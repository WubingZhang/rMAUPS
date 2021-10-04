#' Quality control of proteomic data
#'
#' @docType methods
#' @name countNA
#' @rdname countNA
#'
#' @param expr The expression profile.
#'
#' @return a ggplot object.
#' @import ggplot2 ggpubr
#'
#' @author Wubing Zhang
#' @export
#'
#'
countNA <- function(expr
                    # , filename = NULL, width = 4, height = 3, ...
                    ){
  requireNamespace("ggplot2")
  requireNamespace("ggpubr")
  if(!class(expr) %in% c("data.frame", "matrix")){
    exprSet <- read.table(expr, sep = "\t", header = TRUE, row.names = 1,
                          stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  }else{
    exprSet = expr
  }
  ##  Count NA
  genecount <- rowSums(!is.na(exprSet))
  gg = sapply(1:ncol(exprSet), function(x) sum(genecount>=x))
  gg = data.frame("Ratio" = (1:ncol(exprSet)),
                  "Count" = gg, stringsAsFactors = FALSE)
  p = ggplot(gg, aes_string("Ratio", "Count"))
  p = p + geom_point(color = "gray50")
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))
  p = p + labs(x = "Number of sample", y = "Number of gene")
  # if(!is.null(filename))
  #   ggsave(filename, p, width = width, height = height, ...)
  return(p)
}
