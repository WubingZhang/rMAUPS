#' Pairwise correlation of proteins within the same protein complexes
#'
#' @docType methods
#' @name corComplex
#' @rdname corComplex
#'
#' @param data Proteomics data matrix, with gene or protein as row names.
#' @param idType The input ID type, one of "symbol" (default), "entrez", "uniprot".
#'
#' @return A ggplot instance.
#'
#' @author Wubing Zhang
#'
#' @examples
#'
#' @export

corComplex <- function(data, idType = c("symbol", "entrez", "uniprot")[1]){
  gsets = gsGetter(type = "CORUM")
  if(tolower(idType)=="symbol"){
    gsets$Gene = TransGeneID(gsets$Gene, "Entrez", "Symbol")
    gsets = gsets[gsets$Gene %in% rownames(data), ]
    genesets = unstack(gsets[,c(1,2)])
    genesets = genesets[lengths(genesets)>1]
  }else if(tolower(idType)=="uniprot"){
    gsets$Gene = TransGeneID(gsets$Gene, "Entrez", "uniprot")
    gsets = gsets[gsets$Gene %in% rownames(data), ]
    genesets = unstack(gsets[,c(1,2)])
    genesets = genesets[lengths(genesets)>1]
  }else if(tolower(idType)!="entrez") stop("Unsupport id type ...")

  pcc = unlist(lapply(genesets, function(g) {
    tmp = cor(t(data[g,]), use = "pairwise.complete.obs")
    tmp[upper.tri(tmp)]
  }))
  pcc = pcc[!(is.na(pcc)|pcc==1)]
  gg = data.frame(PCC = pcc)
  p = DensityView(gg, main = "Protein complex")
  p = p + theme(legend.position = "none")
  p = p + labs(x = "Pairwise protein correlation")
  p = p + geom_vline(xintercept = median(gg$PCC), linetype = "dashed")
  p = p + annotate("text", label = round(median(gg$PCC),3), x = median(gg$PCC), y = 0,
                   color = "red", hjust = -0.1, vjust = -1)
  p
}
