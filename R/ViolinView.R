#' Violin plot
#'
#' Violin plot for a given data frame
#'
#' @docType methods
#' @name ViolinView
#' @rdname ViolinView
#'
#' @param gg A data frame.
#' @param x Character, specifying the column name for x plotting.
#' @param y Character, specifying the column name for y plotting.
#' @param color Character, specifying the column name/color.
#' @param fill Character, specifying the column name/fill color.
#' @param group Character, specifying the column for grouping.
#' @param width Numeric, specifying the violin width.
#' @param size Numeric, specifying size of the violin.
#' @param add.box Boolean, whether add boxplot in the figure.
#' @param box.color Similar as `color`.
#' @param box.fill Similar as `fill`.
#' @param box.width Similar as `width`.
#' @param box.size Similar as `size`.
#' @param ... Other parameters in geom_violin.
#'
#' @return An object created by \code{ggplot}, which can be assigned
#' and further customized.
#'
#' @author Wubing Zhang
#'
#' @export

ViolinView <- function(gg, x, y,
                       color = NA,
                       fill = NA,
                       group = NA,
                       width = 0.6,
                       size = 1,
                       add.box = TRUE,
                       box.color = color,
                       box.fill = fill,
                       box.width = 0.3,
                       box.size = size,
                       ...){
  gg = as.data.frame(gg)
  p = ggplot(gg, aes_string(x, y))
  if(group%in%colnames(gg))
    p = ggplot(gg, aes_string(x, y, group=group))
  ## Check if color is valid color
  boo <- try(col2rgb(color), silent=TRUE)
  boo1 = "try-error" %in% class(boo)
  boo <- try(col2rgb(fill), silent=TRUE)
  boo2 = "try-error" %in% class(boo)
  ## Check if color is valid color
  boo <- try(col2rgb(color), silent=TRUE)
  boo1 = "try-error" %in% class(boo)
  boo <- try(col2rgb(fill), silent=TRUE)
  boo2 = "try-error" %in% class(boo)
  if(color %in% colnames(gg)){
    if(fill %in% colnames(gg))
      p = p + geom_violin(aes_string(color = color, fill = fill),
                           width = width, size = size, ...)
    else if(!boo2)
      p = p + geom_violin(aes_string(color = color), fill = fill,
                           width = width, size = size, ...)
    else
      p = p + geom_violin(aes_string(color = color), width = width, size = size, ...)
  }else if(!boo1){
    if(fill %in% colnames(gg))
      p = p + geom_violin(aes_string(fill = fill), color = color,
                           width = width, size = size, ...)
    else if(!boo2)
      p = p + geom_violin(color = color, fill = fill, width = width, size = size, ...)
    else
      p = p + geom_violin(color = color, width = width, size = size, ...)
  }else{
    if(fill %in% colnames(gg))
      p = p + geom_violin(aes_string(fill = fill), width = width, size = size, ...)
    else if(!boo2)
      p = p + geom_violin(fill = fill, width = width, size = size, ...)
    else
      p = p + geom_violin(width = width, size = size, ...)
  }
  ## Add box plot
  boo <- try(col2rgb(box.color), silent=TRUE)
  boo1 = "try-error" %in% class(boo)
  boo <- try(col2rgb(box.fill), silent=TRUE)
  boo2 = "try-error" %in% class(boo)
  if(box.color %in% colnames(gg)){
    if(box.fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(color = box.color, fill = box.fill),
                          width = box.width, size = box.size,
                          outlier.shape = NA)
    else if(!boo2)
      p = p + geom_boxplot(aes_string(color = box.color), fill = box.fill,
                          width = box.width, size = box.size,
                          outlier.shape = NA)
    else
      p = p + geom_boxplot(aes_string(color = box.color), width = box.width,
                           size = box.size, outlier.shape = NA)
  }else if(!boo1){
    if(box.fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(fill = box.fill), color = box.color,
                          width = box.width, size = box.size,
                          outlier.shape = NA)
    else if(!boo2)
      p = p + geom_boxplot(color = box.color, fill = box.fill, width = box.width,
                           size = box.size, outlier.shape = NA)
    else
      p = p + geom_boxplot(color = box.color, width = box.width, size = box.size,
                           outlier.shape = NA)
  }else{
    if(box.fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(fill = box.fill), width = box.width,
                           size = box.size, outlier.shape = NA)
    else if(!boo2)
      p = p + geom_boxplot(fill = box.fill, width = box.width, size = box.size,
                           outlier.shape = NA)
    else
      p = p + geom_boxplot(width = box.width, size = box.size,
                           outlier.shape = NA)
  }
  p = p + theme_bw(base_size = 12)
  return(p)
}
