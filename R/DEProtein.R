#' Differential protein abundance analysis with adjustment of differential RNA expression.
#'
#' @docType methods
#' @name DEProtein
#' @rdname DEProtein
#'
#' @param pData A matrix of protein abundance profile.
#' @param design A vector, indicating the condition of samples in the `pData`.
#' @param rData A matrix of gene expression profile.
#' @param rData.type RNAseq or array.
#' @param toplabels Same as that in the ScatterView.
#'
#' @return A list including the differential analysis results and a scatter plot
#' indicating the adjustment of differential RNA expression from the protein abundance change.
#'
#' @author Wubing Zhang
#' @import ggplot2
#' @export

DEProtein <- function(pData, design, rData = NULL,
                      rData.type = "RNAseq",
                      toplabels = NULL){
  if(is.null(names(design)) & length(design)==ncol(pData)){
    names(design) = colnames(pData)
  }
  samples = intersect(colnames(pData), names(design))
  genes = rownames(pData)
  if(!is.null(rData)) {
    samples = intersect(samples, colnames(rData))
    genes = intersect(genes, rownames(rData))
  }
  if(length(samples)<2) stop("Too few samples ...")
  ann = data.frame(Condition = design[samples], row.names = samples)
  design = stats::model.matrix(~1+Condition, ann)
  rownames(design) = samples

  #"ls" for least squares or "robust" for robust regression
  pfit = limma::eBayes(limma::lmFit(pData[genes, samples], design))
  dep = limma::topTable(pfit, adjust.method="BH", coef=ncol(design), number = Inf)
  dep = dep[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  colnames(dep) = c("log2FC", "baseMean", "stat", "pvalue", "padj")

  # Create result list for return
  res = list(dep = dep, deg = NULL,
             protein.vs.rna.p = NULL,
             adj.vs.orig.p = NULL,
             dep.p = NULL,
             adj.dep.p = NULL)

  efit = NULL
  if(!is.null(rData)){
    if(tolower(rData.type)=="rnaseq"){
      efit = DESeq2::DESeqDataSetFromMatrix(rData[genes, samples],
                                            colData = ann, design = design)
      efit <- DESeq2::DESeq(efit)
      deg <- DESeq2::lfcShrink(efit, coef = ncol(design), quiet = TRUE)
      deg$padj[is.na(deg$padj)] = 1
      deg = deg[, c("log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
      colnames(deg) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(tolower(rData.type)=="array"){
      efit = limma::eBayes(limma::lmFit(rData[genes, samples], design, na.rm=TRUE))
      deg = limma::topTable(efit, adjust.method="BH", coef=ncol(design), number = Inf)
      deg = deg[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
      colnames(deg) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else{ stop("rData.type error ...") }
    deg = as.data.frame(deg)
    res$deg = deg
    ## Adjust the RNA differential expression
    dat = data.frame(dep.stat = dep[genes, "stat"], deg.stat = deg[genes, "stat"],
                     row.names = genes, stringsAsFactors = FALSE)
    dat[is.na(dat)] = 0
    mod = stats::lm(dep.stat ~ deg.stat, data = dat)
    dep[genes, "adj.stat"] = mod$residuals
    dep[genes, "adj.pvalue"] = pt(abs(dep[genes, "stat"]),
                              pfit$df.residual+pfit$df.prior,
                              lower.tail = FALSE)
    dep[, "adj.padj"] = p.adjust(dep[, "adj.pvalue"])
    res$dep = dep

    ## Protein VS RNA level statistics
    intercept = sd(abs(mod$residuals))
    p = ScatterView(dat, "deg.stat", "dep.stat",
                    slope = mod$coefficients[2],
                    intercept = c(-intercept, intercept),
                    groups = c("top", "bottom"),
                    display_cut = FALSE, top = 5,
                    toplabels = toplabels)
    p = p + labs(x = "RNA wald-statistic", y = "Proteomic t-statistic")
    res$protein.vs.rna.p = p

    ## Adjusted VS original results
    tmp = rownames(dep)[order(dep$adj.stat)]
    toplabels = c(tmp[c(1:6, length(tmp):(length(tmp)-5))], toplabels)
    intercept = sd(dep$adj.stat - dep$stat)
    p = ScatterView(dep, "stat", "adj.stat",
                    y_cut = c(-sd(dep$adj.stat), sd(dep$adj.stat)),
                    groups = c("top", "bottom"), display_cut = FALSE, top = 5,
                    toplabels = toplabels)
    p = p + labs(x = "Proteomic t-statistic", y = "Adjusted proteomic t-statistic")
    res$adj.vs.orig.p = p

    ## Volcano plot
    dat = dep
    dat$logP = -log10(dat$pvalue)
    p = ScatterView(dat, "log2FC", "logP", y_cut = 1, groups = c("top"),
                    top = 10, ylab = "-log10(p-value)", display_cut = TRUE)
    p = suppressMessages(p + ggplot2::scale_color_manual(values = c("#f03b20", "gray90")))
    res$dep.p = p

    dat$logP = -log10(dat$adj.pvalue)
    p = ScatterView(dat, "log2FC", "logP", y_cut = 1, groups = c("top"),
                    top = 10, ylab = "-log10(p-value)", display_cut = TRUE)
    p = suppressMessages(p + ggplot2::scale_color_manual(values = c("#f03b20", "gray90")))
    res$adj.dep.p = p

  }
  return(res)
}
