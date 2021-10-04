#' Calculate score across genes and samples
#'
#' This wrapper function combines filtering out genes with low reads in a number of
#' samples (recommended for limma:voom) with normalization
#'
#' @docType methods
#' @name gsScore
#' @rdname gsScore
#'
#' @param dat normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples
#' @param gset Gene sets.
#' @param fun ("PC", default), Pearson, ssGSEA or mean (other value). fisher, stouffer
#' @return numeric vector or matrix.
#' @author Wubing Zhang
#' @importFrom GSVA gsva
#' @importFrom metap sumlog
#' @importFrom matrixStats colMedians
#' @export
gsScore <- function(dat, gset, fun="PC") {
  if(tolower(fun) == "ssgsea"){
    gss <- gsva(dat, gset, method="ssgsea")
    return(gss)
  }
  if(class(gset)!="list"){
    gset = list(gset = gset)
  }
  gss = t(sapply(gset, function(gs){
    gm = dat[rownames(dat)%in%gs, , drop = FALSE]
    if(nrow(gm)<0) return(rep(NA, ncol(dat)))
    if(nrow(gm)==1) return(gm[1,])
    if((nrow(gm)>1 & nrow(gm)<3) | tolower(fun) == "mean"){
      gss <- colMeans(gm)
    }else if(tolower(fun) == "median"){
      gss <- matrixStats::colMedians(as.matrix(gm))
    }else if (tolower(fun) == "pc") {
      pc <- prcomp(t(gm), retx=TRUE)
      gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
    } else if (tolower(fun) == "pearson"){
      corvalues = cor(t(gm))
      diag(corvalues) = 0
      corvalues[corvalues<0] = 0
      w = rowSums(corvalues) / sum(corvalues)
      gss <- (t(gm) %*% w)[,1]
    } else if (tolower(fun) == "fisher"){
      gss = metap::sumlog(gm[,1])$p
    }else if(tolower(fun) == "stouffer"){
      gss = colSums(gm) / sqrt(nrow(gm))
    }else{
      gss = rep(NA, ncol(dat))
    }
    return(gss)
  }))
  if(nrow(gss)==1){
    gss = t(gss)
    rownames(gss) = names(gset)
  }
  return(gss)
}
