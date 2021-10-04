#' Density plot for gene beta scores in Control and Treatment
#'
#' Plot the density of gene beta scores in two samples.
#'
#' @docType methods
#' @name DensityView
#' @rdname DensityView
#'
#' @param beta Data frame, including \code{samples} as columns.
#' @param samples Character, specifying sample names in \code{beta}.
#' @param main As in 'plot'.
#' @param xlab As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{ViolinView}}
#'
#' @examples
#'
#' @importFrom data.table melt
#'
#' @export

DensityView <- function(beta, samples = NULL, main = NULL,xlab = "Beta Score",
                        filename = NULL, width = 5, height = 4, ...){
  if(!is.null(samples) && length(samples)>0){ beta = beta[, samples, drop = FALSE]}
  dd1 = data.table::melt(beta,id=NULL)
  if(!"variable" %in% colnames(dd1)){
    dd1$variable = colnames(beta)
  }
  #==========
  p=ggplot(data=dd1,aes_string(x="value",color="variable",group="variable"))
  p=p+geom_density()
  # p=p+facet_wrap(~variable,nrow=1)
  p=p+labs(color=NULL)
  p=p+theme(legend.justification = c(1, 1), legend.position = c(0.99, 0.99))
  # p=p+theme(legend.text = element_text(size=8))
  p=p+labs(x=xlab, y="Density", title=main)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}

#' Identical bar plot
#'
#' Identical bar plot
#'
#' @docType methods
#' @name IdentBarView
#' @rdname IdentBarView
#'
#' @param gg A data frame.
#' @param x A character, indicating column (in countSummary) of x-axis.
#' @param y A character, indicating column (in countSummary) of y-axis.
#' @param fill A character, indicating fill color of all bars.
#' @param main A charater, specifying the figure title.
#' @param xlab A character, specifying the title of x-axis.
#' @param ylab, A character, specifying the title of y-axis.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in ggsave.
#'
#' @author Wubing Zhang
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#'
#' @examples
#'
#' @import ggplot2
#' @export

IdentBarView <- function(gg, x = "x", y = "y", fill = c("#CF3C2B", "#394E80"),
                         main = NULL, xlab = NULL, ylab = NULL,
                         filename = NULL, width = 5, height = 4, ...){
  gg$x = gg[, x]
  gg$y = gg[, y]
  p <- ggplot(gg)
  p = p + geom_bar(aes(x, y), stat="identity", width=0.6, fill = fill[1], alpha=0.8)
  p = p + labs(x=xlab, y=ylab, title=main)
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + theme(text = element_text(colour="black",size = 14),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"),
                axis.text.x=element_text(angle = 45, hjust=1, vjust = 1))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())

  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}

#' Bar plot
#'
#' Bar plot
#'
#' @docType methods
#' @name BarView
#' @rdname BarView
#'
#' @param df A data frame.
#' @param x A character, specifying the x-axis.
#' @param y A character, specifying the x-axis.
#' @param fill A character, specifying the fill color.
#' @param bar.width A numeric, specifying the width of bar.
#' @param position "dodge" (default), "stack", "fill".
#' @param dodge.width A numeric, set the width in position_dodge.
#' @param main A charater, specifying the figure title.
#' @param xlab A character, specifying the title of x-axis.
#' @param ylab A character, specifying the title of y-axis.
#' @param ... Other parameters in geom_bar
#'
#' @author Wubing Zhang
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#'
#' @examples
#' mdata = data.frame(group=letters[1:5], count=sample(1:100,5))
#' BarView(mdata, x = "group", y = "count")
#' @import ggplot2 ggpubr
#' @importFrom grDevices col2rgb
#' @export

BarView <- function(df, x = "x", y = "y", fill = "#FC6665",
                    bar.width = 0.8, position = "dodge",
                    dodge.width = 0.8, main = NA,
                    xlab = NULL, ylab = NA, ...){

  ## Check if fill is valid color
  boo <- try(grDevices::col2rgb(fill), silent=TRUE)
  boo = "try-error" %in% class(boo)

  ## Use the order of x in the df
  df[,x] = factor(df[,x], levels = unique(df[,x]))

  if(boo){
    p <- ggplot(df, aes_string(x, y, fill=fill))
    if(position=="dodge"){
      p <- p + geom_bar(width = bar.width, stat="identity",
                        position=position_dodge(width = dodge.width,
                                                preserve = "single"), ...)
    }else{
      p <- p + geom_bar(width = bar.width, stat="identity", position=position, ...)
    }
  }else{
    p <- ggplot(df, aes_string(x, y))
    p = p + geom_bar(fill=fill, stat="identity", position=position, ...)
  }
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + labs(fill = NULL)
  if(!(length(xlab)==1 && is.na(xlab))) p = p + labs(x=xlab)
  if(!(length(ylab)==1 && is.na(ylab))) p = p + labs(y=ylab)
  if(!(length(main)==1 && is.na(main))) p = p + labs(title=main)
  p = p + theme_pubr()
  p = p + theme(plot.title = element_text(hjust = 0.5, size=16))
  return(p)
}

#' Violin plot
#'
#' Plots the violin of beta scores in Control and Treatment samples.
#'
#' @docType methods
#' @name ViolinView
#' @rdname ViolinView
#' @aliases violinview
#'
#' @param beta Data frame, , including \code{samples} as columns.
#' @param samples Character, specifying the name of samples to be compared.
#' @param main As in 'plot'.
#' @param ylab As in 'plot'.
#' @param filename Figure file name to create on disk. Default filename="NULL", which means
#' don't save the figure on disk.
#' @param width As in ggsave.
#' @param height As in ggsave.
#' @param ... Other available parameters in function 'ggsave'.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @seealso \code{\link{DensityView}}
#'
#' @examples
#'
#' @importFrom data.table melt
#'
#' @export
#'
ViolinView <- function(beta, samples=NULL, main=NULL, ylab="Beta Score",
                       filename=NULL, width=5, height=4, ...){
  requireNamespace("data.table", quietly=TRUE) || stop("need data.table package")

  if(!is.null(samples) && length(samples)>1){ beta = beta[, samples]}

  dd1 = melt(beta, id=NULL)
  if(!"variable" %in% colnames(dd1)){
    dd1$variable = colnames(beta)
  }
  #======
  p=ggplot(data=dd1,aes_string(x="variable",y="value",color="variable"))
  p=p+geom_violin()+geom_boxplot(width=.1, outlier.colour=NA)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(),
                legend.key = element_rect(fill = "transparent"))
  p=p+theme(legend.position = "none")
  p=p+labs(x=NULL,y=ylab,title=main)

  if(!is.null(filename)){
    ggsave(plot=p,filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}


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


#' Scatter plot
#'
#' Scatter plot supporting groups.
#'
#' @docType methods
#' @name ScatterView
#' @rdname ScatterView
#' @aliases ScatterView
#'
#' @param data Data frame.
#' @param x A character, specifying the x-axis.
#' @param y A character, specifying the y-axis.
#' @param label An integer or a character specifying the column used as the label, default value is 0 (row names).
#'
#' @param model One of "none" (default), "ninesquare", "volcano", and "rank".
#' @param x_cut An one or two-length numeric vector, specifying the cutoff used for x-axis.
#' @param y_cut An one or two-length numeric vector, specifying the cutoff used for y-axis.
#' @param slope A numberic value indicating slope of the diagonal cutoff.
#' @param intercept A numberic value indicating intercept of the diagonal cutoff.
#' @param auto_cut Boolean, take 1.5 fold standard deviation as cutoff.
#' @param auto_cut_x Boolean, take 1.5 fold standard deviation as cutoff on x-axis.
#' @param auto_cut_y Boolean, take 1.5 fold standard deviation as cutoff on y-axis.
#' @param auto_cut_diag Boolean, take 1.5 fold standard deviation as cutoff on diagonal.
#'
#' @param groups A character vector specifying groups. Optional groups include "top", "mid", "bottom",
#' "left", "center", "right", "topleft", "topcenter", "topright", "midleft", "midcenter",
#' "midright", "bottomleft", "bottomcenter", "bottomright".
#' @param group_col A vector of colors for specified groups.
#' @param groupnames A vector of group names to show on the legend.
#'
#' @param label.top Boolean, specifying whether label top hits.
#' @param top Integer, specifying the number of top terms in the groups to be labeled.
#' @param toplabels Character vector, specifying terms to be labeled.
#'
#' @param display_cut Boolean, indicating whether display the dashed line of cutoffs.
#'
#' @param color A character, specifying the column name of color in the data frame.
#' @param shape A character, specifying the column name of shape in the data frame.
#' @param size A character, specifying the column name of size in the data frame.
#'
#' @param main Title of the figure.
#' @param xlab Title of x-axis
#' @param ylab Title of y-axis.
#' @param legend One of "top", "bottom", "left", "right", and "none".
#' @param ... Other available parameters in function 'geom_text_repel'.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @examples
#'
#' @import ggplot2 ggpubr ggrepel
#' @export
#'
#'

ScatterView<-function(data, x = "x", y = "y", label = 0,
                      model = c("none", "ninesquare", "volcano", "rank")[1],
                      x_cut = NULL, y_cut = NULL, slope = 1, intercept = NULL,
                      auto_cut = FALSE, auto_cut_x = auto_cut,
                      auto_cut_y = auto_cut, auto_cut_diag = auto_cut,
                      groups = NULL, group_col = NULL, groupnames = NULL,
                      label.top = TRUE, top = 0, toplabels = NULL,
                      display_cut = FALSE, color = NULL, shape = 16, size = 1,
                      main = NULL, xlab = x, ylab = y, legend = "none", ...){
  requireNamespace("ggplot2", quietly=TRUE) || stop("need ggplot package")
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  requireNamespace("ggpubr", quietly=TRUE) || stop("need ggpubr package")
  data = as.data.frame(data, stringsAsFactors = FALSE)
  data = data[!(is.na(data[,x])|is.na(data[,y])), ]
  ## Add label column in the data frame.
  if(label==0) data$Label = rownames(data)
  else data$Label = as.character(data[, label])

  ## Compute the cutoff used for each dimension.
  model = tolower(model)
  if(model == "ninesquare"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
    if(length(y_cut)==0)
      y_cut = c(-CutoffCalling(data[,y], 1.5), CutoffCalling(data[,y], 1.5))
    if(length(intercept)==0)
      intercept = c(-CutoffCalling(data[,y]-data[,x], 1.5), CutoffCalling(data[,y]-data[,x], 1.5))
  }
  if(model == "volcano"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
    if(length(y_cut)==0) y_cut = -log10(0.05)
  }
  if(model == "rank"){
    if(length(x_cut)==0)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
  }
  if(model == "none"){
    if(auto_cut_x)
      x_cut = c(-CutoffCalling(data[,x], 1.5), CutoffCalling(data[,x], 1.5))
    if(auto_cut_y)
      y_cut = c(-CutoffCalling(data[,y], 1.5), CutoffCalling(data[,y], 1.5))
    if(auto_cut_diag)
      intercept = c(-CutoffCalling(data[,y]-data[,x], 1.5), CutoffCalling(data[,y]-data[,x], 1.5))
  }
  ## Decide the colored groups
  avail_groups = c("topleft", "topright", "bottomleft", "bottomright",
                   "midleft", "topcenter", "midright", "bottomcenter", "midcenter",
                   "top", "mid", "bottom", "left", "center", "right", "none")
  ## Select the colors
  mycolour = c("#1f78b4", "#fb8072", "#33a02c", "#ff7f00",
               "#bc80bd", "#66c2a5", "#6a3d9a", "#fdb462", "#ffed6f",
               "#e78ac3", "#fdb462", "#8da0cb", "#66c2a5", "#fccde5", "#fc8d62", "#d9d9d9")
  names(mycolour) = avail_groups

  if(model == "ninesquare") groups = c("midleft", "topcenter", "midright", "bottomcenter")
  if(model == "volcano") groups = c("topleft", "topright")
  if(model == "rank") groups = c("left", "right")
  groups = intersect(groups, avail_groups)

  ## Annotate the groups in the data frame
  if(length(x_cut)>0){
    idx1 = data[,x] < min(x_cut)
    idx2 = data[,x] > max(x_cut)
  }else{
    idx1 = NA
    idx2 = NA
  }
  if(length(y_cut)>0){
    idx3 = data[,y] < min(y_cut)
    idx4 = data[,y] > max(y_cut)
  }else{
    idx3 = NA
    idx4 = NA
  }
  if(length(intercept)>0){
    idx5 = data[,y]<slope*data[,x]+min(intercept)
    idx6 = data[,y]>slope*data[,x]+max(intercept)
  }else{
    idx5 = NA; idx6 = NA
  }
  data$group="none"
  for(gr in groups){
    if(gr=="topleft") idx = cbind(idx1, idx4, idx6)
    if(gr=="topcenter") idx = cbind(!idx1, !idx2, idx4, idx6)
    if(gr=="topright") idx = cbind(idx2, idx4, idx6)
    if(gr=="midleft") idx = cbind(idx1, idx6 , !idx3, !idx4)
    if(gr=="midcenter") idx = cbind(!idx1, !idx2, !idx3, !idx4, !idx5, !idx6)
    if(gr=="midright") idx = cbind(idx2, !idx3, !idx4, idx5)
    if(gr=="bottomleft") idx = cbind(idx1, idx3, idx5)
    if(gr=="bottomcenter") idx = cbind(!idx1, !idx2, idx3, idx5)
    if(gr=="bottomright") idx = cbind(idx2, idx3, idx5)
    if(gr=="top"){
      if(length(y_cut)>0 & length(intercept)>0)
        idx = idx4 & idx6
      else if(length(y_cut)>0)
        idx = idx4
      else idx = idx6
    }
    if(gr=="mid") idx = (!idx3) & (!idx4)
    if(gr=="bottom"){
      if(length(y_cut)>0 & length(intercept)>0)
        idx = idx3 & idx5
      else if(length(y_cut)>0)
        idx = idx3
      else idx = idx5
    }
    if(gr=="left"){
      if(length(x_cut)>0 & length(intercept)>0)
        if(slope>0) idx = idx1 & idx6 else idx = idx1 & idx5
        else if(length(x_cut)>0)
          idx = idx1
        else
          if(slope>0) idx = idx6 else idx = idx5
    }
    if(gr=="center") idx = (!idx1) & (!idx2)
    if(gr=="right"){
      if(length(x_cut)>0 & length(intercept)>0)
        if(slope>0) idx = idx2 & idx5 else idx = idx2 & idx6
        else if(length(x_cut)>0)
          idx = idx2
        else
          if(slope>0) idx = idx5 else idx = idx6
    }
    ## Assign groups
    if(is.null(ncol(idx))){
      if(sum(!is.na(idx))>0) data$group[idx] = gr
      else warning("No cutpoint for group:", gr)
    }else{
      idx = idx[, !is.na(idx[1,])]
      if(is.null(ncol(idx)))
        warning("No cutpoint for group:", gr)
      else if(ncol(idx)<4 & gr=="midcenter")
        warning("No cutpoint for group:", gr)
      else
        data$group[rowSums(idx)==ncol(idx)] = gr
    }
  }
  data$group=factor(data$group, levels = unique(c(groups, "none")))
  ## Group names
  if(length(groupnames)!=length(groups)) groupnames = groups
  if(length(groups)>0) names(groupnames) = groups
  if(length(group_col)==length(groups)) mycolour[groups] = group_col
  if(length(groups)==0) mycolour["none"] = "#FF6F61"

  ## Label top gene names ##
  data$rank = top + 1
  for(g in groups){
    idx1 = data$group==g
    x_symb = 0; y_symb = 0;
    if(g=="topleft"){ x_symb = 1; y_symb = -1 }
    if(g=="topcenter"){ x_symb = 0; y_symb = -1 }
    if(g=="topright"){ x_symb = -1; y_symb = -1 }
    if(g=="midleft"){ x_symb = 1; y_symb = 0 }
    if(g=="midright"){ x_symb = -1; y_symb = 0 }
    if(g=="bottomleft"){ x_symb = 1; y_symb = 1 }
    if(g=="bottomcenter"){ x_symb = 0; y_symb = 1 }
    if(g=="bottomright"){ x_symb = -1; y_symb = 1 }
    if(g=="top"){ x_symb = 0; y_symb = -1 }
    if(g=="bottom"){ x_symb = 0; y_symb = 1 }
    if(g=="left"){ x_symb = 1; y_symb = 0 }
    if(g=="right"){ x_symb = -1; y_symb = 0 }
    tmp = data[,c(x,y)]
    tmp[,x] = (tmp[,x]-min(tmp[,x])) / (max(tmp[,x])-min(tmp[,x]))
    tmp[,y] = (tmp[,y]-min(tmp[,y])) / (max(tmp[,y])-min(tmp[,y]))
    data$rank[idx1] = rank((x_symb*tmp[,x]+y_symb*tmp[,y])[idx1])
  }
  data$rank[data$rank==0] = Inf
  if(mode(toplabels)=="list"){
    data$Label[data$rank>top & !(data$Label %in% unlist(toplabels))] = ""
    data$group = data$Label;
    if(length(toplabels)>0){
      tmp = stack(toplabels)
      tmp = tmp[!duplicated(tmp[,1]), ]
      rownames(tmp) = tmp[,1]
      data$group[data$group%in%tmp[,1]] = as.character(tmp[data$group[data$group%in%tmp[,1]], 2])
      data$group[!(data$group%in%tmp[,2]) & data$group!=""] = "Top hits"
    }
  }else{
    data$Label[data$rank>top & !(data$Label %in% toplabels)] = ""
  }

  ## Color issue
  if(is.null(color)){
    color = "group"
  }else if(!color%in%colnames(gg)){
    color = "group"
    mycolour["none"] = color
  }
  ## Plot the scatter figure ##
  gg = data
  ## Plot the figure
  p = ggplot(gg, aes_string(x, y, label="Label", color = color))
  if(all(c(shape,size)%in%colnames(gg)))
    p = p + geom_point(aes(shape = shape, size = size), alpha = 0.6)
  else if(shape%in%colnames(gg))
    p = p + geom_point(aes(shape = shape), size = size, alpha = 0.6)
  else if(size%in%colnames(gg))
    p = p + geom_point(aes(size = size), shape = shape, alpha = 0.6)
  else
    p = p + geom_point(size = size, shape = shape, alpha = 0.6)

  ## Color
  if(color=="group"){
    if(mode(toplabels)!="list")
      p = p + scale_color_manual(values = mycolour, labels = groupnames)
    else
      p = p + scale_color_manual(values = c("#d9d9d9", "#fb8072", "#80b1d3", "#fdb462", "#bc80bd", "#b3de69", "#bebada", "#8dd3c7", "#ffffb3", "#fccde5", "#ccebc5", "#ffed6f"))
  }else if(color%in%colnames(gg)){
    if(mode(gg[,color])=="numeric")
      p = p + scale_color_gradient2(low = "#377eb8", high = "#e41a1c", midpoint = 0)
    else if(!"try-error"%in%class(try(col2rgb(gg[1,color]),silent=TRUE))){
      mycolour = unique(gg[,color]); names(mycolour) = mycolour
      p = p + scale_color_manual(values = mycolour)
    }
  }
  # else if(!"try-error"%in%class(try(col2rgb(x),silent=TRUE)))
  #   p = p + scale_color_manual(values = color)

  if(label.top)
    p = p + ggrepel::geom_text_repel(...)
  if(display_cut){
    if(length(x_cut)>0)
      p = p + geom_vline(xintercept = x_cut,linetype = "dotted")
    if(length(y_cut)>0)
      p = p + geom_hline(yintercept = y_cut,linetype = "dotted")
    if(length(intercept)>0)
      p = p + geom_abline(slope=slope, intercept=intercept, linetype = "dotted")
  }
  p = p + labs(x=xlab, y = ylab, title = main, color = NULL)
  p = p + ggpubr::theme_pubr(legend = legend) + theme(plot.title = element_text(hjust = 0.5))

  return(p)
}


#' Quantile of normal distribution.
#'
#' Compute cutoff from a normal-distributed vector.
#'
#' @docType methods
#' @name CutoffCalling
#' @rdname CutoffCalling
#'
#' @param d A numeric vector.
#' @param scale Boolean or numeric, specifying how many standard deviation will be used as cutoff.
#'
#' @return A numeric value.
#'
#' @import stats utils
#' @export
#' @examples
#' CutoffCalling(rnorm(10000))

CutoffCalling=function(d, scale=1){
  param=1
  if(is.logical(scale) & scale){
    param = round(length(d) / 20000, digits = 1)
  }else if(is.numeric(scale)){param = scale}

  Control_mean=0
  sorted_beta=sort(abs(d))
  temp=stats::quantile(sorted_beta,0.68)
  temp_2=stats::qnorm(0.84)
  cutoff=round(temp/temp_2,digits = 3)
  names(cutoff)=NULL
  cutoff=cutoff*param
  return(cutoff)
}

