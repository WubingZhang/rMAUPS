#' Normalization of ProteomeDiscoverer-exported proteomic data
#'
#' @docType methods
#' @name normalizeProteomeDiscoverer
#' @rdname normalizeProteomeDiscoverer
#'
#' @param prot A file (or folder) path of proteomics data. Or a data frame of proteomics data,
#' which includes columns of "Contaminant", "Number of unique peptides",
#' "Gene Symbol", and protein abundance of samples.
#' @param output The path to the output file or directory (when \code{prot} is a folder path).
#' @param norm method for normalization.
#' @param log2 Boolean, whether perform log2 transformation.
#' @param return Boolean, whether return the data object.
#'
#' @return A data frame.
#' @author Wubing Zhang
#' @export
#'
normalizeProteomeDiscoverer <- function(prot, output = NULL,
                                        norm = "median", log2 = FALSE,
                                        return = FALSE){
  if(class(prot)=="character" && dir.exists(prot)){# Process all export data in a folder
    files = list.files(prot, pattern = "export.*txt$", full.names = TRUE)
    for(f in files){
      outfile = basename(gsub("export.*", "normdata.csv", f))
      if(!is.null(output)) outfile = file.path(output, outfile)
      tmp = normalizeProteomeDiscoverer(f, output = outfile, norm = norm, log2 = log2)
    }
    return(TRUE)
  }else if(class(prot)=="character" && file.exists(prot)){# Process an export data
    message(format(Sys.time(), "%b-%d-%Y %X "), "Processing ", prot)
    dd = read.table(prot, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, comment.char = "")
  }else if(is.data.frame(prot)){
    dd = prot
  }else{return()}

  # Remove protein contaminants
  idx1 = toupper(dd$Contaminant)!="TRUE"
  # Remove data with unique peptides below threshold
  idx2 = dd[,grep("Unique.*Peptides", colnames(dd), value = TRUE)]>=2
  dd <- dd[idx1&idx2, ]

  # Identify and count reporter ion channels
  channels <- grep("^Abundance.*Sample", colnames(dd), value=TRUE)
  psm <- grep("PSMs", colnames(dd), value=TRUE)[1]
  df <- dd[, c(psm, channels)]

  genename = grep("symbol",colnames(dd), value = TRUE, ignore.case = TRUE)
  dupgenes = dd[, genename][duplicated(dd[, genename])]
  dupdf = t(sapply(unique(dupgenes), function(g){
    idx = dd[, genename]==g
    colMeans(df[idx, , drop = FALSE], na.rm = TRUE)
  }))
  idx = dd[, genename] %in% dupgenes
  df = rbind.data.frame(df[!idx, ], dupdf)
  allgenes <- c(dd[!idx, genename], rownames(dupdf))
  idx = is.na(allgenes) | duplicated(allgenes) | allgenes==""
  df = df[!idx, ]
  rownames(df) <- allgenes[!idx]
  colnames(df) <- gsub("Abundance.|Sample.", "", colnames(df))
  df[,1] = round(df[,1])
  colnames(df)[1] = "PSMs"
  # Remove data with NAs or low sum of reporter ion intensities
  idx1 <- rowSums(df[,-1], na.rm = TRUE)>0
  df <- df[idx1, ]

  normdf = df
  normdf[,-1] = normalizeProteomics(df[,-1], norm, log2)
  normdf = as.data.frame(normdf, stringsAsFactors = FALSE)
  if(!is.null(output)){
    write.csv(normdf, output, row.names = TRUE, quote = FALSE)
  }
  if(return) return(normdf) else return(TRUE)
}

#' Normalize Fisher's data
#'
#' @docType methods
#' @name normalizeProteomics
#' @rdname normalizeProteomics
#'
#' @param df A matrix-like object.
#' @param norm method for normalization, such as none, median, medianratio, scale, robustz, quantile, loess.
#' @param log2 Boolean, whether perform log2 transformation.
#'
#' @return A matrix.
#' @author Wubing Zhang
#' @export
#'
normalizeProteomics <- function(df, norm = "median", log2 = FALSE){
  normdf = df
  if(norm=="medianratio"){
    # Median ratio normalization
    geomeans <- exp(rowMeans(log(normdf)))
    Sizes <- apply(normdf, 2, function(cnts)
      median((cnts/geomeans)[geomeans > 0]))
    normdf = t(t(normdf)/Sizes)
    if(log2) normdf = log2(normdf)
  }else if(norm=="median"){
    if(log2){
      mid = apply(log2(normdf), 2, median, na.rm = TRUE)
      normdf = t(t(log2(normdf)) - mid)
    }else{
      mid = apply(normdf, 2, median, na.rm = TRUE)
      normdf = t(t(normdf) - mid)
    }
  }else if(norm=="scale"){
    if(log2){
      normdf = scale(log2(normdf))
    }else{
      normdf = scale(normdf)
    }
  }else if(norm=="robustz"){
    if(log2){
      M = apply(log2(normdf), 2, median, na.rm = TRUE)
      MAD = apply(log2(normdf), 2, mad, na.rm = TRUE)
      normdf = t((t(log2(normdf)) - M) / MAD)
    }else{
      M = apply(normdf, 2, median, na.rm = TRUE)
      MAD = apply(normdf, 2, mad, na.rm = TRUE)
      normdf = t((t(normdf) - M) / MAD)
    }
  }else if(norm=="quantile"){
    normdf = limma::normalizeQuantiles(normdf)
    if(log2) normdf = log2(normdf)
  }else if(norm=="loess"){
    normdf = limma::normalizeCyclicLoess(normdf)
  }
  return(normdf)
}
