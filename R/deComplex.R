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
#' @param lfc Gene sets.
#' @param pval ("PC", default), Pearson, ssGSEA or mean (other value). fisher, stouffer.
#'
#' @return A list.
#' @author Wubing Zhang
#' @export
#'
DEComplex <- function(depres, type = "GOBP+GOCC+CORUM+REACTOME",
                      lfc = "log2FC", pval = "pvalue"){
  ## GOCC, CORUM, Pathway
  gsets = gsGetter(type = type, limit = c(0,200))
  gsets$Symbol = TransGeneID(gsets$Gene, "Entrez", "Symbol")
  gsets = gsets[gsets$Symbol %in% rownames(depres), ]
  genesets = unstack(gsets[,c(4,2)])

  fisherp = gsScore(depres[,pval, drop=FALSE], gset = genesets, fun = "fisher")
  stoufferZ = gsScore(depres[,lfc, drop=FALSE], gset = genesets, fun = "stouffer")
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

  ## Visualize differential expressed complexes and pathways
  merged_deres$logP = -log10(merged_deres$pvalue)
  tmp = merged_deres[grepl("C5_BP", rownames(merged_deres)), ]
  p1 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 3, main = "GOBP", ylab = "-log10(p.value)")
  p1 = p1 + theme(legend.position = "none")
  tmp = merged_deres[grepl("REACTOME", rownames(merged_deres)), ]
  p2 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 3, main = "REACTOME", ylab = "-log10(p.value)")
  p2 = p2 + theme(legend.position = "none")
  tmp = merged_deres[grepl("C5_CC", rownames(merged_deres)), ]
  tmp = tmp[tmp$Count<80,]
  p3 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 3, main = "GOCC", ylab = "-log10(p.value)")
  p3 = p3 + theme(legend.position = "none")
  tmp = merged_deres[grepl("CORUM", rownames(merged_deres)), ]
  tmp = tmp[tmp$Count<50,]
  p4 = ScatterView(tmp, x = "Zscore", y = "logP", label = "Description",
                   model = "volcano", auto_cut_x = TRUE, force = 5,
                   top = 3, main = "CORUM", ylab = "-log10(p.value)")
  p4 = p4 + theme(legend.position = "none")
  return(list(gocc.p = p3, corum.p = p4, gobp.p = p1, reactome.p = p2,
              deComplex = merged_deres))
}
