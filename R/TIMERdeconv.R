#' Deconvolute the immune cell types.
#'
#' @param XX A matrix of immune expression data.
#' @param YY A matrix of cancer expression data.
#' @param w Weight.
#'
#' @return A list.
GetFractions.Abbas <- function(XX, YY, w=NA){
  ## XX is immune expression data
  ## YY is cancer expression data
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(setNames(rep(0, ncol(XX)), nm = colnames(XX)))
      tmp.XX=tmp.XX[, -ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0) return(setNames(rep(0, ncol(XX)), nm = colnames(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX, YY, intercept=F) else tmp=lsfit(tmp.XX, YY, w, intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0, ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}

#' Combine the cancer expression data and the immune expression data
#'
#' @param cancer.exp A data frame or matrix of raw read counts of bulk RNA-seq samples.
#' Column names correspond to sample names, row names to genes.
#' @param immune.exp A matrix of read counts representing the gene expression profiles of
#' all cell types that might contribute to the bulk samples.
#' @param immune.cellType Immune cell types.
#'
#' @importFrom sva ComBat
#' @return A list.
#'
RemoveBatchEffect <- function(cancer.exp, immune.exp, immune.cellType) {
  ## intersect the gene names of cancer.exp and immune.exp
  tmp.dd <- as.matrix(cancer.exp)
  tmp <- sapply(strsplit(rownames(cancer.exp), '\\|'),
                function(x) x[[1]])
  rownames(tmp.dd) <- tmp
  tmp.dd <- as.matrix(tmp.dd[which(nchar(tmp)>1), ])
  tmp.ss <- intersect(rownames(tmp.dd), rownames(immune.exp))
  ## bind cancer and immune expression data into one dataframe
  N1 <- ncol(tmp.dd)
  tmp.dd <- cbind(tmp.dd[tmp.ss, ], immune.exp[tmp.ss, ])
  tmp.dd <- as.matrix(tmp.dd)
  mode(tmp.dd) <- 'numeric'
  ## remove batch effects
  N2 <- ncol(immune.exp)
  tmp.batch <- c(rep(1, N1), rep(2, N2))
  tmp.dd0 <- sva::ComBat(tmp.dd, tmp.batch, c())
  ## separate cancer and immune expression data after batch effect removing
  dd.br <- tmp.dd0[, 1:N1]
  immune.exp.br <- tmp.dd0[, (N1+1):(N1+N2)]
  ## a immune category has multiple samples, use the median expression level for a gene
  tmp0 <- c()
  for(kk in unique(names(immune.cellType))){
    tmp.vv <- which(names(immune.cellType)==kk)
    tmp0 <- cbind(tmp0, apply(immune.exp.br[, tmp.vv], 1, median, na.rm=T))
  }
  immune.exp.agg.br <- tmp0
  colnames(immune.exp.agg.br) <- unique(names(immune.cellType))
  return(list(as.matrix(dd.br), immune.exp.br, immune.exp.agg.br))
}

#' Detection of the cancer type
#'
#' @param exprsn A data frame or matrix of gene expression profile.
#'
signatureDetect <- function(exprsn){
  avgTPM = file.path(system.file("extdata", package = "rMAUPS"), "deconv/avgTPM.Rdata")
  cancerTypeVector = get(load(avgTPM))
  ol_genes = intersect(rownames(cancerTypeVector), rownames(exprsn))
  sp_correlation = apply(cor(exprsn[ol_genes,],cancerTypeVector[ol_genes,],
                             method='spearman',use='complete.obs'),1,which.max)
  max_cnt = table(colnames(cancerTypeVector)[sp_correlation])
  names(max_cnt)[which.max(max_cnt)]
}

#'Runs TIMER for decomposing bulk RNA-seq samples
#'
#' @param e A data frame or matrix of raw read counts of bulk RNA-seq samples.
#' Column names correspond to sample names, row names to genes.
#' @param cancer A matrix of read counts representing the gene expression profiles of
#' all cell types that might contribute to the bulk samples.
#' @importFrom plyr revalue
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column column_to_rownames
#' @export
TIMERdeconv <- function(e, cancer = "AUTO") {
  pureImmune = file.path(system.file("extdata", package = "rMAUPS"),
                         "deconv/timer_pureImmune.Rdata")
  geneMarker = file.path(system.file("extdata", package = "rMAUPS"),
                         "deconv/timer_geneMarker.Rdata")
  immune = NULL
  load(pureImmune)
  load(geneMarker)
  cancer <- ifelse(cancer == 'AUTO', signatureDetect(e), cancer)
  immune.geneExpression <- immune$genes
  immune.cellTypes <- immune$celltypes
  cancer.expression <- e # replace file
  outlier.genes <- sort(unique(c(as.matrix(apply(cancer.expression, 2, function(x) {
    rownames(cancer.expression)[tail(order(x), 5)]})))))
  cancer.expression <- cancer.expression[!(rownames(cancer.expression) %in% outlier.genes), , drop=F]
  d.rmBatch <- RemoveBatchEffect(cancer.expression, immune.geneExpression, immune.cellTypes)
  cancer.expNorm <- d.rmBatch[[1]]
  immune.expNormMedian <- d.rmBatch[[3]]
  gene.selected.marker <- intersect(geneMarker[[cancer]], row.names(cancer.expNorm))
  X_immune = immune.expNormMedian[gene.selected.marker, c(-4)]
  Y_cancer = cancer.expNorm[gene.selected.marker, , drop=FALSE]
  transName <- c(B_cell = "B cell",
                 T_cell.CD4 = "T cell CD4+",
                 T_cell.CD8 = "T cell CD8+",
                 Neutrophil = "Neutrophil",
                 Macrophage = "Macrophage",
                 DC = "Myeloid dendritic cell")
  res = as.data.frame(apply(Y_cancer, 2, function(y){
    GetFractions.Abbas(X_immune, y)}))
  rownames(res) = transName[rownames(res)]
  res
}
