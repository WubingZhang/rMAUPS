#' Standardize the Proteome_CDAP_Protein_Report data
#' @param CDAP A file path, specifying the COAP report data for processing.
#' You can download the CPTAC data (Proteome_CDAP_Protein_Report or Phosphoproteome_CDAP_Protein_Report)
#' from https://cptc-xfer.uis.georgetown.edu/publicData/.
#' @param biospecimens.csv Path to the biospecimens csv file, with "Aliquot" (corresponding to column names in CDAP)
#' and "Sample" (Real sample names) in the data.
#' @param summary Path to the summary file, which includes counts of the peptide spectrum match.
#' @param outdir Output directory.
#' @param prefix The prefix for the output files.
#'
#' @return A data frame, representing the stardardized proteome/phosphoproteome data,
#' in which the first column is median of PSM across multile parallels (if summary is available).
#'
#' @export
preprocess_cptac <- function(CDAP, biospecimens.csv = "", summary = "",
                             outdir = NULL, prefix = ""){
  options(stringsAsFactors = FALSE)
  pmat = read.table(CDAP, row.names = 1, sep = "\t",
                    header = TRUE, quote = "", check.names = FALSE)
  type = gsub(".*_", "", gsub("_CDAP.*", "", CDAP))
  if("NCBIGeneID" %in% colnames(pmat)){
    pmat = pmat[!is.na(pmat$NCBIGeneID), grepl("Unshared Log Ratio",colnames(pmat))]
    colnames(pmat) = gsub("\\..*| .*|_.*", "", colnames(pmat))
    dupsamples = unique(colnames(pmat)[duplicated(colnames(pmat))])
    if(length(dupsamples)>0){
      dupmat = sapply(dupsamples, function(x) rowMeans(pmat[,colnames(pmat)==x], na.rm = TRUE))
      pmat = cbind(pmat[, !(colnames(pmat)%in%dupsamples)], dupmat)
    }
  }else if("Gene" %in% colnames(pmat)){
    uniprotFile = file.path(system.file("extdata", package = "rMAUPS"),
                            "uniprot_protein_mapping_hsa.rds")
    uniprot_mapping = readRDS(uniprotFile)
    uniprot_mapping = uniprot_mapping[!(is.na(uniprot_mapping$refseq)|duplicated(uniprot_mapping$refseq)), ]
    rownames(uniprot_mapping) = uniprot_mapping$refseq
    pmat = read.table(CDAP, row.names = 1, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
    pmat = pmat[!(is.na(pmat$Gene)|pmat$Gene==""), ]
    pmat$uniprot = uniprot_mapping[gsub(":.*", "", rownames(pmat)), "uniprot"]
    pmat = pmat[!is.na(pmat$uniprot), ]
    tmp = pmat[pmat$uniprot%in%uniprot_mapping$Canonical, ]
    tmp$uniprot = gsub("-.*", "", tmp$uniprot)
    pmat = rbind(pmat[!(gsub("-.*","",pmat$uniprot)%in%gsub("-.*","",tmp$uniprot)), ], tmp)
    sites = paste0(pmat$Gene, ":", pmat$uniprot, ":", gsub(".*:","",rownames(pmat)))
    pmat = pmat[, grepl("Log Ratio",colnames(pmat))]
    colnames(pmat) = gsub("\\..*| .*|_.*", "", colnames(pmat))
    dupsamples = unique(colnames(pmat)[duplicated(colnames(pmat))])
    if(length(dupsamples)>0){
      dupmat = sapply(dupsamples, function(x) rowMeans(pmat[,colnames(pmat)==x], na.rm = TRUE))
      pmat = cbind(pmat[, !(colnames(pmat)%in%dupsamples)], dupmat)
    }
    dupsites = unique(sites[duplicated(sites)])
    if(length(dupsites)>0){
      dupmat = sapply(dupsites, function(x) colMeans(pmat[rownames(pmat)==x, ], na.rm = TRUE))
      pmat = pmat[!(sites%in%dupsites), ]
      rownames(pmat) = sites[!(sites%in%dupsites)]
      pmat = rbind(pmat, t(dupmat))
    }else{
      rownames(pmat) = sites
    }
  }else{
    stop("Incorrect CDAP protein report data ...")
  }
  if(file.exists(biospecimens.csv)){
    sampleann = read.csv(biospecimens.csv, header = TRUE)
    sampleann$Aliquot = gsub("\\..*| .*|_.*", "", sampleann$Aliquot)
    sampleann = sampleann[!duplicated(sampleann$Aliquot), ]
    rownames(sampleann) = sampleann$Aliquot
    idx = colnames(pmat) %in% rownames(sampleann)
    colnames(pmat)[idx] = sampleann[colnames(pmat)[idx], "Sample"]
    if(!is.null(outdir)){
      write.csv(sampleann, file.path(outdir, paste0(prefix, "_Biospecimens.csv")),
                row.names = FALSE, quote = FALSE)
    }
  }
  if(file.exists(summary)){
    summary = read.table(summary, sep = "\t", row.names = 1, header = TRUE,
                         quote = "",check.names = FALSE)
    summary = summary[!is.na(summary$NCBIGeneID), grepl("Unshared",colnames(summary))]
    PSMs = floor(Biobase::rowMedians(as.matrix(summary)))
    names(PSMs) = rownames(summary)
    pmat = cbind.data.frame(PSM = PSMs[gsub(":.*", "", rownames(pmat))], pmat)
  }
  if(grepl("TCGA", CDAP)){
    idx = nchar(colnames(pmat))==11
    colnames(pmat)[idx] = paste0("TCGA-", colnames(pmat)[idx])
    colnames(pmat)[idx] = gsub(".$", "", colnames(pmat)[idx])
  }
  if(!is.null(outdir)){
    write.csv(pmat, file.path(outdir, paste0(prefix, "_", type, ".csv")), quote = FALSE)
  }
  return(pmat)
}
