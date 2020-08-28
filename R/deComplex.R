#' Convert gene level differences to pathway level or complex level differences.
#'
#' @docType methods
#' @name DEComplex
#' @rdname DEComplex
#'
#' @param depres A data frame, with log2FC and pvalue in the column.
#' @param type Molecular signatures for testing, available datasets include Pathway
#' (KEGG, REACTOME, C2_CP), GO (GOBP, GOCC, GOMF), MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA),
#' C2_CGP), C3 (C3_MIR, C3_TFT), C4, C6, C7, HALLMARK) and Complex (CORUM). Any combination of them
#' are also accessible (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param lfc A character indicating the column of lfc in the depres.
#' @param pval A character indicating the column of pvalue in the depres.
#' @param method One of "gsea" (default), fisher.
#' @param limit A two-length vector, specifying the size of genesets for calculation.
#' @param ... Other available parameters in enrich.GSE.
#'
#' @return A list.
#' @author Wubing Zhang
#' @importFrom ggplot2 theme
#' @importFrom MAGeCKFlute enrich.GSE
#' @export
#'
DEComplex <- function(depres, type = "GOCC+CORUM",
                      lfc = "log2FC", pval = "pvalue",
                      method = "gsea", limit = c(2,200), ...){
  # requireNamespace("ggplot2")
  # requireNamespace("MAGeCKFlute")

  if(method=="fisher"){
    ## GOCC, CORUM, Pathway
    gsets = gsGetter(type = type, limit = limit)
    gsets$Symbol = TransGeneID(gsets$Gene, "Entrez", "Symbol")
    gsets = gsets[gsets$Symbol %in% rownames(depres), ]
    genesets = unstack(gsets[,c(4,2)])
    stoufferZ = gsScore(depres[,lfc, drop=FALSE], gset = genesets, fun = "stouffer")

    fisherp = gsScore(depres[,pval, drop=FALSE], gset = genesets, fun = "fisher")
    Gene = lapply(genesets, function(x){
      ogene = unique(intersect(x,rownames(depres)))
      Gene = paste0(ogene, collapse = ",")
      Count = length(ogene)
      c(Gene, Count)
    })
    Gene = unlist(Gene)
    merged_deres = data.frame(Gene = Gene[seq(1,length(Gene),2)],
                              Count = Gene[seq(2,length(Gene),2)],
                              stringsAsFactors = FALSE)
    merged_deres$Count = as.integer(merged_deres$Count)
    merged_deres = cbind.data.frame(merged_deres, Zscore = stoufferZ)
    merged_deres = cbind.data.frame(merged_deres, pvalue = fisherp)
    rownames(merged_deres) = names(genesets)
    tmp = gsets[!duplicated(gsets$PathwayID), ]
    rownames(tmp) = tmp$PathwayID
    merged_deres$Description = tmp[rownames(merged_deres), 3]
    merged_deres = merged_deres[, c(5,1:4)]
    merged_deres$p.adjust = p.adjust(merged_deres$pvalue, method = "fdr")
  }else if(method=="gsea"){
    gsets = gsGetter(type = type, limit = limit)
    gsets$Symbol = TransGeneID(gsets$Gene, "Entrez", "Symbol")
    gsets = gsets[gsets$Symbol %in% rownames(depres), ]
    genesets = unstack(gsets[,c(4,2)])
    stoufferZ = gsScore(depres[,lfc, drop=FALSE], gset = genesets, fun = "stouffer")

    genelist = -log10(depres[,pval])
    genelist[depres[,lfc]<0] = -genelist[depres[,lfc]<0]
    names(genelist) = rownames(depres)
    enrichRes = MAGeCKFlute::enrich.GSE(genelist, keytype = "symbol",
                                        type = type, pvalueCutoff = 1,
                                        limit = limit, ...)
    merged_deres = enrichRes@result
    merged_deres$Zscore = stoufferZ[rownames(merged_deres),1]
  }

  ## Visualize differential expressed complexes and pathways
  merged_deres$logFDR = -log10(merged_deres$p.adjust)
  res = list(deComplex = merged_deres)
  types = unlist(strsplit(type, split = "\\+"))
  types = gsub("BP|CC|MF", "", types)
  for(i in types){
    if(sum(grepl(i, rownames(merged_deres)))<2) next
    tmp = merged_deres[grepl(i, rownames(merged_deres)), ]
    p = ScatterView(tmp, x = "Zscore", y = "logFDR", label = "Description",
                     model = "volcano", auto_cut_x = TRUE, force = 5,
                     top = 5, main = i, ylab = "-log10(FDR)")
    p = p + ggplot2::theme(legend.position = "none")
    res[[paste0(tolower(i), ".p")]] = p
  }
  return(res)
}
