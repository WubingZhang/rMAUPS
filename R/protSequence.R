#' Search protein sequences using a regular expression
#'
#' @docType methods
#' @name searchProtSeq
#' @rdname searchProtSeq
#'
#' @param prot_seq_df A data.frame containing uniprot IDs under the "ID" column and protein sequence under the "protein_sequence" column
#' @param regex a string containing the regular expression to search against the protein sequences
#'
#' @return A data frame object
#'
#' @author Collin Tokheim
#'
#' @importFrom purrr map imap reduce compact
#' @importFrom dplyr mutate
#' @importFrom stringr str_locate_all
#' @export
searchProtSeq <- function(prot_seq_df, regex){
  # regex search all protein sequences
  motif_hits <- stringr::str_locate_all(prot_seq_df$protein_sequence, regex)

  # format motif hits
  motif_hits <- purrr::map(motif_hits, as.data.frame)
  names(motif_hits) <- prot_seq_df$ID
  motif_hits <- purrr::imap(motif_hits, ~.x %>% dplyr::mutate(UniprotId = .y))
  motif_hits_df <- purrr::reduce(purrr::compact(motif_hits), rbind)

  myorder <- c('UniprotId', 'start', 'end')
  return(motif_hits_df[,myorder])
}


#' Reads in the ubiquitination site data from Phosphosite Plus database
#'
#' @docType methods
#' @name readPSPUbiquitin
#' @rdname readPSPUbiquitin
#'
#' @param path File path to Ubiquitination_site_dataset from phosphositeplus db
#' @param minStudies Minimum number of studies reporting UB site to keep it
#'
#' @return A data frame object
#'
#' @author Collin Tokheim
#'
#' @importFrom stringr str_sub str_length
#' @export
readPSPUbiquitin <- function(path, minStudies=1) {
  # read in ub site data
  ub <- read.delim(path, sep='\t', stringsAsFactors=F)

  # only keep human cases
  ub = ub[ub['ORGANISM']=='human',]

  # extract out the position of the UB
  ub['position'] <- stringr::str_sub(ub$MOD_RSD, 2, stringr::str_length(ub$MOD_RSD)-3)
  ub$position <- as.numeric(ub$position)

  # only keep sites that have been reported at least minStudies times
  numReports <- rowSums(ub[,c('LT_LIT', 'MS_LIT', 'MS_CST')], na.rm=T)
  ub <- ub[numReports>=minStudies,]

  # rename columns
  names(ub)[names(ub)=='ACC_ID'] <- 'ID'
  names(ub)[names(ub)=='GENE'] <- 'gene'

  return(ub)
}


#' View protein sequence positions on protein structure
#'
#' @docType methods
#' @name browseProtStructure
#' @rdname browseProtStructure
#'
#' @param protId a uniprot id string
#' @param start a vector of start positions
#' @param end a vector of end positions
#' @param doBrowse a boolean indicator on whether to open browser to view protein structure
#' @param baseUrl string for location of mupit service
#' @param checkBaseUrl string for location of mupit service to check availability of protein structure
#'
#' @return None
#'
#' @author Collin Tokheim
#'
#' @importFrom httr GET content
#' @export
browseProtStructure <- function(protId, start, end,
                                doBrowse=TRUE,
                                baseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/?gm=',
                                checkBaseUrl='https://mupit.icm.jhu.edu/MuPIT_Interactive/rest/showstructure/check?pos='){

  stopifnot(length(start)==length(end))

  index = lapply(1:length(start), function(i) seq(start[i], end[i]))
  uniProtSeq = paste(paste0(protId, ":", unlist(index)), collapse = ",")

  # check if prot structure available
  jsonResponse <- httr::GET(paste0(checkBaseUrl, uniProtSeq, '&protquery=y'))
  jsonResponseParsed <- httr::content(jsonResponse, as="parsed")

  # construct url
  fullUrl <- paste0(baseUrl, uniProtSeq, '&protquery=y')
  # browse url if there is available prot structure
  if (jsonResponseParsed$hit){
    cat(paste0(fullUrl, '\n'))
    if (!doBrowse) {
      # pass, do nothing
    } else if (Sys.getenv('R_BROWSER')!="") {
      browseURL(fullUrl)
    } else {
      cat('\nBrowser not set!\n')
    }
  } else {
    cat('\nNo protein structure available!\n')
  }
}
