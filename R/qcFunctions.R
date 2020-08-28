#' Get the number of detected proteins in the given data.
#'
#' @param dat A matrix-like object of proteomics data, in which rows are proteins,
#' and columns are samples.
#' @param samples A character vector, specifying the samples for the analysis.
#' @param condition A named character vector, specifying the condition of the samples.
#'
#' @return A data frame, in which the first column includes sample names,
#' and the second column includes the number of detected proteins.
#' @author Wubing Zhang
#' @export
#'
getDetection <- function(dat, samples = NA, condition = NA){
  detection = colSums(!is.na(dat))
  res = data.frame(Sample = colnames(dat), Detection = detection,
                   stringsAsFactors = FALSE)
  if(sum(res$Sample%in%samples)>1) res = res[res$Sample%in%samples, ]
  if(length(condition)>1) res$Condition = condition[res$Sample]
  return(res)
}

#' Get the number of missing proteins in the given data.
#'
#' @param dat A matrix-like object of proteomics data, in which rows are proteins,
#' and columns are samples.
#' @param samples A vector of samples, specifying the samples for the analysis.
#'
#' @return A data frame, in which the first column records the percentage of samples (X),
#' and the second column includes the number of missing proteins in more than X samples.
#' @author Wubing Zhang
#' @export
#'
getMissing <- function(dat, samples = NA,
                       percentage = seq(0,0.99,0.01)){
  if(sum(colnames(dat)%in%samples)>1) dat = dat[, colnames(dat)%in%samples]
  Missing = sapply(percentage, function(x)
    sum(rowSums(is.na(dat))>x*ncol(dat)))
  res = data.frame(Percentage = percentage*100, Missing = Missing,
                   stringsAsFactors = FALSE)
  return(res)
}


#' Pairwise correlation of proteins within the same protein complexes
#'
#' @docType methods
#' @name corComplex
#' @rdname corComplex
#'
#' @param dat A matrix-like object of proteomics data, with gene or protein as row names.
#' @param idType A character specifying the protein ID type, one of "symbol" (default), "entrez", "uniprot".
#'
#' @return A data frame including three columns, Protein1, Protein2, and PCC (Pearson Correlation Coefficient).
#'
#' @author Wubing Zhang
#'
#' @examples
#'
#' @export
corComplex <- function(dat, idType = c("symbol", "entrez", "uniprot")[1]){
  gsets = gsGetter(type = "CORUM")
  if(tolower(idType)=="symbol"){
    gsets$Gene = TransGeneID(gsets$Gene, "Entrez", "Symbol")
    gsets = gsets[gsets$Gene %in% rownames(dat), ]
    genesets = unstack(gsets[,c(1,2)])
    genesets = genesets[lengths(genesets)>1]
  }else if(tolower(idType)=="uniprot"){
    gsets$Gene = TransGeneID(gsets$Gene, "Entrez", "uniprot")
    gsets = gsets[gsets$Gene %in% rownames(dat), ]
    genesets = unstack(gsets[,c(1,2)])
    genesets = genesets[lengths(genesets)>1]
  }else if(tolower(idType)!="entrez") stop("Unsupport id type ...")

  pcc = unlist(lapply(genesets, function(g) {
    tmp = cor(t(dat[g,]), use = "pairwise.complete.obs")
    res = tmp[upper.tri(tmp)]
    pairs = matrix(paste0(rep(g, each = length(g)), "-", rep(g, length(g))), nrow = length(g))
    names(res) = pairs[upper.tri(pairs)]
    return(res)
  }))
  res = data.frame(Protein1 = gsub("-.*|.*\\.", "", names(pcc)),
                   Protein2 = gsub(".*-", "", names(pcc)),
                   PCC = pcc, stringsAsFactors = FALSE)
  res = res[order(-abs(res$PCC)), ]
  idx = is.na(res$PCC)|duplicated(paste0(res$Protein1, "-", res$Protein2))|res$Protein1==res$Protein2
  res = res[!idx, ]
  return(res)
}

#' Correlation between samples/genes in two different data sets.
#'
#' @docType methods
#' @name DataCorrelation
#' @rdname DataCorrelation
#'
#' @param dat1 A matrix-like object of proteomics data.
#' @param dat2 A matrix-like object of RNA-Seq data or RPPA data.
#' @param method One of "pearson", "kendall", or "spearman" (default).
#'
#' @return A two-column data frame, the first column is the `Correlation`, and the
#' second column is `Type`, either 'Sample' or 'Gene'.
#'
#' @import matrixStats
#' @author Wubing Zhang
#' @export

DataCorrelation <- function(dat1, dat2, method="spearman"){
  ## Get data matrix
  dat1 <- as.matrix(dat1)
  dat2 <- as.matrix(dat2)

  ## Overlap genes and samples
  overlap_gene <- intersect(rownames(dat1), rownames(dat2))
  overlap_sample <- intersect(colnames(dat1), colnames(dat2))
  if(com_gene<100 | com_sample<3){
    message("Number of overlap genes: ", com_gene,
            "\nNumber of overlap samples: ", com_sample)
    stop("limited number of samples or genes!")
  }

  ## Remove Zero Standard deviations
  idx1 = matrixStats::rowSds(dat1[overlap_gene, overlap_sample], na.rm = TRUE)!=0
  idx2 = matrixStats::rowSds(dat2[overlap_gene, overlap_sample], na.rm = TRUE)!=0
  overlap_gene = overlap_gene[idx1&idx2]

  ## Calculate correlation
  corvalues_g = unlist(lapply(overlap_gene, function(x){
    cor(dat1[x, overlap_sample], dat2[x, overlap_sample],
        method = method, use = "pairwise.complete.obs")
  }))
  corvalues_s = unlist(lapply(overlap_sample, function(y){
    cor(dat1[overlap_gene, y], dat2[overlap_gene, y],
        method = method, use = "pairwise.complete.obs")
  }))

  ## Sample correlation
  res = data.frame(Correlation = c(corvalues_s, corvalues_g),
                   Type = rep(c("Sample", "Gene"),
                   c(length(overlap_sample), length(overlap_gene))))
  return(res)
}
