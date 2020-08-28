#' Boxplot
#'
#' Boxplot for a given data frame
#'
#' @docType methods
#' @name BoxView
#' @rdname BoxView
#'
#' @param gg A data frame.
#' @param x Character, specifying the column name for x plotting.
#' @param y Character, specifying the column name for y plotting.
#' @param color Character, specifying the column name/color.
#' @param fill Character, specifying the column name/fill color.
#' @param group Character, specifying the column for grouping.
#' @param width Numeric, specifying the box width.
#' @param size Numeric, specifying size of the box.
#' @param add.jitter Boolean, whether add jitter into the plot.
#' @param jitter.color Character, specifying the column name/color of the points.
#' @param jitter.size Numeric, specifying size of the jitter.
#' @param ... Other parameters in geom_boxplot.
#'
#' @return An object created by \code{ggplot}, which can be assigned
#' and further customized.
#'
#' @author Wubing Zhang
#'
#' @export

BoxView <- function(gg, x, y,
                    color = NA,
                    fill = NA,
                    group = NA,
                    width = 0.6,
                    size = 1,
                    add.jitter = FALSE,
                    jitter.color = color,
                    jitter.size = size,
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
  if(color %in% colnames(gg)){
    if(fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(color = color, fill = fill),
                           width = width, size = size, ...)
    else if(!boo2)
      p = p + geom_boxplot(aes_string(color = color), fill = fill,
                           width = width, size = size, ...)
    else
      p = p + geom_boxplot(aes_string(color = color), width = width, size = size, ...)
  }else if(!boo1){
    if(fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(fill = fill), color = color,
                           width = width, size = size, ...)
    else if(!boo2)
      p = p + geom_boxplot(color = color, fill = fill, width = width, size = size, ...)
    else
      p = p + geom_boxplot(color = color, width = width, size = size, ...)
  }else{
    if(fill %in% colnames(gg))
      p = p + geom_boxplot(aes_string(fill = fill), width = width, size = size, ...)
    else if(!boo2)
      p = p + geom_boxplot(fill = fill, width = width, size = size, ...)
    else
      p = p + geom_boxplot(width = width, size = size, ...)
  }
  if(add.jitter){
    boo <- try(col2rgb(jitter.color), silent=TRUE)
    boo3 = "try-error" %in% class(boo)
    if(jitter.color %in% colnames(gg))
      p = p + geom_jitter(aes_string(color = jitter.color), size = jitter.size)
    else if(!boo3)
      p = p + geom_jitter(color = jitter.color, size = jitter.size)
    else
      p = p + geom_jitter(size = jitter.size)
  }
  p = p + theme_bw(base_size = 12)
  return(p)
}
