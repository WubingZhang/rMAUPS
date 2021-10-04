#' Filter out rows with NA or low value
#'
#' @docType methods
#' @name filterN
#' @rdname filterN
#'
#' @param m A matrix-like object.
#' @param minS A numeric, specifying the minimum proportion/number of values should
#' be quantified for each row.
#' @param out A vector, specifying low-quality values (NA and 0).
#' @param impute Imputation method.
#' @return A matrix with the same columns as input matrix.
#' @author Wubing Zhang
#' @export
#'
filterN <- function(m, minS = 3, out = NA, impute = "none"){
  m = as.matrix(m)

  # if(is.null(design)){
  if(minS<1) minS = minS*ncol(m)
  idx = is.na(m)
  if(length(out[!is.na(out)])>0)
    idx = idx | (m==max(out[!is.na(out)]))
  sel = rowSums(!idx)>=minS
  m = imputeNA(m[sel,], method = impute)
  return(m)
}

#' Imputation
#'
#' @docType methods
#' @name imputeNA
#' @rdname imputeNA
#'
#' @param m A matrix-like object.
#' @param method method for imputation, such as knn, lowAbundanceResampling, ReplicateBasedResampling
#' @param k Integer, parameter for knn.
#' @param rowmax parameter for knn.
#' @param colmax parameter for knn.
#'
#' @return A matrix.
#' @author Wubing Zhang
#' @import impute
#' @export
#'
imputeNA <- function(m, method = "knn", k = 30, rowmax = 0.95, colmax = 0.95){
  requireNamespace("impute")
  imputed_m = m
  if(tolower(method) == "knn"){
    imputed_m = impute.knn(m, k = k, rowmax = rowmax, colmax = colmax,
                           maxp = floor(nrow(m)/1000)*1000)$data
  }else if(tolower(method) == "lowabundanceresampling"){
    imputed_m = lowAbundanceResampling(m)
  }else if(tolower(method) == "replicatebasedresampling"){
    imputed_m = replicateBasedResampling(m)
  }
  return(imputed_m)
}

#' low Abundance Resampling method from protein discover
#'
#' @docType methods
#' @name lowAbundanceResampling
#' @rdname lowAbundanceResampling
#'
#' @param df Matrix-like object.
#' @param percent cutoff for low abundance values.
lowAbundanceResampling <- function(df, percent = 0.05){
  distr = df[!is.na(df)]
  df[is.na(df)] = rnorm(sum(is.na(df)), mean(distr), sd(distr))
  return(df)
}

#' Replicate based resampling method from protein discover
#'
#' @docType methods
#' @name replicateBasedResampling
#' @rdname replicateBasedResampling
#'
#' @param df Matrix-like object.
#'
replicateBasedResampling <- function(df){
  replicate = rep(1, ncol(df))
  for(i in unique(replicate)){
    mid = rowMeans(df[, replicate==i], na.rm = TRUE)
    sd = matrixStats::rowSds(df[, replicate==i])
    mod = lm(sd~mid, data.frame(sd = sd[!is.na(sd)], mid = mid[!is.na(sd)]))
    pred = predict(mod, data.frame(mid = mid))
    for(j in which(is.na(sd)&(!is.na(mid)))){
      tmpj = is.na(df[j, replicate==i])
      df[j, which(replicate==i)[tmpj]] = rnorm(sum(tmpj), mid[j], pred[j])
    }
    return(lowAbundanceResampling(df))
  }
}
