#' Merge all the rMAUPS results
#'
#' @docType methods
#' @name arrangePngs
#' @rdname arrangePngs
#'
#' @param outdir Path to rMAUPS results.
#'
#' @return A list of merged figures.
#'
#' @author Wubing Zhang
#'
#' @import ggpubr
#' @export

arrangePngs <- function(outdir){
  allfolders = list.files(outdir)
  pnglist = list()

  ## QC folder
  if("qc" %in% allfolders){
    rds_list = list.files(paste0(outdir, "/qc"), ".rds", full.names = TRUE)
    qc = list()
    for(rds in rds_list){
      tmp = readRDS(rds)[c(1:6)]
      tmp = tmp[lengths(tmp)>0]
      qc[[gsub(".*\\/|.qc.*", "", rds)]] = ggarrange(plotlist = tmp, ncol = 3)
    }
    pnglist[["qc"]] = qc
  }

  ## Imputation folder
  if("imputation" %in% allfolders){
    rds_list = list.files(paste0(outdir, "/imputation"), ".rds", full.names = TRUE)
    impute = list()
    for(rds in rds_list){
      tmp = readRDS(rds)[c(1:6)]
      tmp = tmp[lengths(tmp)>0]
      impute[[gsub(".*\\/|.qc.*", "", rds)]] = ggarrange(plotlist = tmp, ncol = 3)
    }
    pnglist[["impute"]] = impute
  }

  ## Comparisons
  allfolders = setdiff(allfolders, c("qc", "imputations"))
  if(length(allfolders)>0){
    rds_list = list.files(outdir, ".rds", full.names = TRUE, recursive = TRUE)
    rds_list = rds_list[grepl("Comparison_.*/", rds_list)]
    merged_rds = rds_list[grepl("merged.rds$", rds_list)]
    rds_list = setdiff(rds_list, merged_rds)
    compare = list()
    for(rds in rds_list){
      tmp = readRDS(rds)
      tmp = tmp[-length(tmp)]
      tmp = tmp[lengths(tmp)>0]
      compare[[gsub(".*\\/|.rds", "", rds)]] = ggarrange(plotlist = tmp)
    }
    pnglist[["compare"]] = compare

    if(length(merged_rds)>0){
      integrate = list()
      for(rds in merged_rds){
        tmp = readRDS(rds)[1:4]
        tmp = tmp[lengths(tmp)>0]
        integrate[[gsub(".*\\/|_merged.rds", "", rds)]] = ggarrange(plotlist = tmp)
      }
      pnglist[["integrate"]] = integrate
    }
  }
  return(pnglist)
}
