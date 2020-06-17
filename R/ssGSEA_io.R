########################################
### GCT class and method definitions ###
########################################

#' An S4 class to represent a GCT object
#'
#' @slot mat a numeric matrix
#' @slot rid a character vector of row ids
#' @slot cid a character vector of column ids
#' @slot rdesc a \code{data.frame} of row descriptors
#' @slot rdesc a \code{data.frame} of column descriptors
#' @slot src a character indicating the source (usually file path) of the data
#'
#' @description The GCT class serves to represent annotated
#'   matrices. The \code{mat} slot contains said data and the
#'   \code{rdesc} and \code{cdesc} slots contain data frames with
#'   annotations about the rows and columns, respectively
#'
methods::setClass("GCT",
                  methods::representation(
                    mat = "matrix",
                    rid = "character",
                    cid = "character",
                    rdesc = "data.frame",
                    cdesc = "data.frame",
                    version = "character",
                    src = "character"
                  )
)


# set up methods for checking GCT validity
methods::setValidity("GCT",
                     function(object) {
                       # check whether dimensions of various
                       # slots are in sync
                       nrows <- nrow(object@mat)
                       ncols <- ncol(object@mat)
                       if (nrows != length(object@rid)) {
                         return("rid must be the same length as number of matrix rows")
                       }
                       if (ncols != length(object@cid)) {
                         return("cid must be the same length as number of matrix columns")
                       }
                       if (length(object@cid) > length(unique(object@cid))) {
                         return("cid must be unique")
                       }
                       if (length(object@rid) > length(unique(object@rid))) {
                         return("rid must be unique")
                       }
                       if (nrow(object@cdesc) != ncols & nrow(object@cdesc) != 0) {
                         return("cdesc must either have 0 rows or the same number of rows as matrix has columns")
                       }
                       if (nrow(object@rdesc) != nrows & nrow(object@rdesc) != 0) {
                         return("rdesc must either have 0 rows or the same number of rows as matrix has rows")
                       }
                       else {
                         return(T)
                       }
                     }
)

# suppressMessages({
#   # set method for displaying a GCT object
#   # just use the 'str' function to show its structure
#   setMethod("show", methods::signature("GCT"), function(object) {
#     utils::str(object)
#   })
#
#   # dim, nrow and ncol to display the # of rows and columns
#   # for a GCT object's matrix
#   setMethod("ncol", methods::signature("GCT"), function(x) {
#     ncol(x@mat)
#   })
#   setMethod("nrow", methods::signature("GCT"), function(x) {
#     nrow(x@mat)
#   })
#   setMethod("dim", methods::signature("GCT"), function(x) {
#     dim(x@mat)
#   })
#   setMethod("range", methods::signature("GCT"), function(x, na.rm=F, finite=F) {
#     range(x@mat, na.rm=na.rm, finite=finite)
#   })
#   setMethod("max", methods::signature("GCT"), function(x, na.rm=F) {
#     max(x@mat, na.rm=na.rm)
#   })
#   setMethod("min", methods::signature("GCT"), function(x, na.rm=F) {
#     min(x@mat, na.rm=na.rm)
#   })
#   setMethod("diag", methods::signature("GCT"), function(x) {
#     diag(x@mat)
#   })
# })


# # define the initialization method for the GCT class
# methods::setMethod("initialize",
#                    signature = "GCT",
#                    definition = function(.Object, mat=NULL, rdesc=NULL, cdesc=NULL, src=NULL, rid=NULL, cid=NULL,
#                                          matrix_only=FALSE) {
#                      # if we were supplied a matrix and annotations, use them
#                      if (!is.null(mat)) {
#                        .Object@mat <- mat
#                        # if given rid and cid, use those as well
#                        if (!is.null(rid)) {
#                          .Object@rid <- rid
#                        } else {
#                          .Object@rid <- rownames(mat)
#                        }
#                        if (!is.null(cid)) {
#                          .Object@cid <- cid
#                        } else {
#                          .Object@cid <- colnames(mat)
#                        }
#                      }
#                      if (!is.null(rdesc)) {
#                        .Object@rdesc <- rdesc
#                      }
#                      if (!is.null(cdesc)) {
#                        .Object@cdesc <- cdesc
#                      } else if (!is.null(src)) {
#                        # we were not given a matrix, were we given a src file?
#                        # check to make sure it's either .gct or .gctx
#                        if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
#                          stop("Either a .gct or .gctx file must be given")
#                        if (grepl(".gct$", src)) {
#                          if ( ! is.null(rid) || !is.null(cid) )
#                            warning(paste("rid and cid values may only be given for .gctx files, not .gct files\n",
#                                          "ignoring"))
#                          # parse the .gct
#                          .Object@src <- src
#                          # get the .gct version by reading first line
#                          .Object@version <- scan(src, what = "", nlines = 1, sep = "\t", quiet = TRUE)[1]
#                          # get matrix dimensions by reading second line
#                          dimensions <- scan(src, what = double(0), nlines = 1, skip = 1, sep = "\t", quiet = TRUE)
#                          nrmat <- dimensions[1]
#                          ncmat <- dimensions[2]
#                          if (length(dimensions)==4) {
#                            # a #1.3 file
#                            message("parsing as GCT v1.3")
#                            nrhd <- dimensions[3]
#                            nchd <- dimensions[4]
#                          } else {
#                            # a #1.2 file
#                            message("parsing as GCT v1.2")
#                            nrhd <- 0
#                            nchd <- 0
#                          }
#                          message(paste(src, nrmat, "rows,", ncmat, "cols,", nrhd, "row descriptors,", nchd, "col descriptors"))
#                          # read in header line
#                          header <- scan(src, what = "", nlines = 1, skip = 2, sep = "\t", quote = NULL, quiet = TRUE)
#                          # construct row header and column id's from the header line
#                          if ( nrhd > 0 ) {
#                            rhd <- header[2:(nrhd+1)]
#                            cid <- header[-(nrhd+1):-1]
#                            col_offset <- 1
#                          }
#                          else {
#                            if (any(grepl("description", header, ignore.case=T))) {
#                              # check for presence of description column in v1.2 files
#                              col_offset <- 2
#                            } else {
#                              col_offset <- 1
#                            }
#                            rhd <- NULL
#                            cid <- header[(1+col_offset):length(header)]
#                          }
#                          # read in the next set of headers (column annotations) and shape into a matrix
#                          if ( nchd > 0 ) {
#                            header <- scan(src, what = "", nlines = nchd, skip = 3, sep = "\t",
#                                           quote = NULL, quiet = TRUE)
#                            header <- matrix(header, nrow = nchd,
#                                             ncol = ncmat + nrhd + 1, byrow = TRUE)
#                            # extract the column header and column descriptions
#                            chd <- header[,1]
#                            cdesc <- header[,-(nrhd+1):-1]
#                            # need to transpose in the case where there's only one column annotation
#                            if ( nchd == 1 )
#                              cdesc <- t(cdesc)
#                          }
#                          else {
#                            chd = NULL
#                            cdesc <- data.frame(id=cid)
#                          }
#                          # read in the data matrix and row descriptions, shape into a matrix
#                          mat <- scan(src, what = "", nlines = nrmat,
#                                      skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
#                          mat <- matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + col_offset,
#                                        byrow = TRUE)
#                          # message(paste(dim(mat), collapse="\t"))
#                          # Extract the row id's row descriptions, and the data matrix
#                          rid <- mat[,1]
#                          if ( nrhd > 0 ) {
#                            # need as.matrix for the case where there's only one row annotation
#                            rdesc <- as.matrix(mat[,2:(nrhd + 1)])
#                            mat <- matrix(as.numeric(mat[,-(nrhd + 1):-1]),
#                                          nrow = nrmat, ncol = ncmat)
#                          }
#                          else {
#                            rdesc <- data.frame(id=rid)
#                            mat <- matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]), nrow = nrmat, ncol = ncmat)
#                          }
#                          # assign names to the data matrix and the row and column descriptions
#                          # message(paste(dim(mat), collapse="\t"))
#                          dimnames(mat) <- list(rid, cid)
#                          if ( nrhd > 0 ) {
#                            dimnames(rdesc) <- list(rid, rhd)
#                            rdesc <- as.data.frame(rdesc, stringsAsFactors = FALSE)
#                          }
#                          if ( nchd > 0 ) {
#                            cdesc <- t(cdesc)
#                            dimnames(cdesc) <- list(cid,chd)
#                            cdesc <- as.data.frame(cdesc, stringsAsFactors = FALSE)
#                          }
#                          # assign to the GCT slots
#                          .Object@mat <- mat
#                          .Object@rid <- rownames(mat)
#                          .Object@cid <- colnames(mat)
#                          if (!matrix_only) {
#                            # return annotations as well as matrix
#                            .Object@rdesc <- fix.datatypes(rdesc)
#                            .Object@cdesc <- fix.datatypes(cdesc)
#                            # add id columns to rdesc and cdesc
#                            .Object@rdesc$id <- rownames(.Object@rdesc)
#                            .Object@cdesc$id <- rownames(.Object@cdesc)
#                          }
#                        }
#                        else {
#                          # parse the .gctx
#                          message(paste("reading", src))
#                          .Object@src <- src
#                          # if the rid's or column id's are .grp files, read them in
#                          if ( length(rid) == 1 && grepl(".grp$", rid) )
#                            rid <- parse.grp(rid)
#                          if ( length(cid) == 1 && grepl(".grp$", cid) )
#                            cid <- parse.grp(cid)
#                          # get all the row and column ids
#                          all_rid <- read.gctx.ids(src, dimension="row")
#                          all_cid <- read.gctx.ids(src, dimension="col")
#                          # if rid or cid specified, read only those rows/columns
#                          # if already numeric, use as is
#                          # else convert to numeric indices
#                          processed_rids <- process_ids(rid, all_rid, type="rid")
#                          processed_cids <- process_ids(cid, all_cid, type="cid")
#                          # read the data matrix
#                          .Object@mat <- rhdf5::h5read(src, name="0/DATA/0/matrix",
#                                                       index=list(processed_rids$idx, processed_cids$idx))
#                          # set the row and column ids, casting as characters
#                          .Object@rid <- processed_rids$ids
#                          .Object@cid <- processed_cids$ids
#                          rownames(.Object@mat) <- processed_rids$ids
#                          colnames(.Object@mat) <- processed_cids$ids
#                          # get the meta data
#                          if (!matrix_only) {
#                            .Object@rdesc <- read.gctx.meta(src, dimension="row", ids=processed_rids$ids)
#                            .Object@cdesc <- read.gctx.meta(src, dimension="col", ids=processed_cids$ids)
#                          }
#                          else {
#                            .Object@rdesc <- data.frame(id=.Object@rid, stringsAsFactors = F)
#                            .Object@cdesc <- data.frame(id=.Object@cid, stringsAsFactors = F)
#                          }
#                          # close any open handles and return the object
#                          if(utils::packageVersion('rhdf5') < "2.23.0") {
#                            rhdf5::H5close()
#                          } else {
#                            rhdf5::h5closeAll()
#                          }
#                          message("done")
#                        }
#                      }
#                      # finally, make sure object is valid before returning
#                      ok <- methods::validObject(.Object)
#                      return(.Object)
#                    }
# )

#' Read a GMT file and return a list
#' @param fname the file path to be parsed
#'
#' @return a list of the contents of \code{fname}. See details.
#'
#' @details \code{parse.gmt} returns a nested list object. The top
#'   level contains one list per row in \code{fname}. Each of
#'   these is itself a list with the following fields:
#'   - \code{head}: the name of the data (row in \code{fname})
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#'
#' @examples
#' gmt_path <- system.file("extdata", "query_up.gmt", package="cmapR")
#' gmt <- parse.gmt(gmt_path)
#' str(gmt)
#'
#' @family CMap parsing functions
#' @export
parse.gmt <- function(fname) {
  gmt.lines <- scan(fname, what = "", sep = "\n",
                    quote = NULL, quiet = TRUE)
  tmp <- lapply(gmt.lines, function(x) unlist(strsplit(x, "\t")))
  mk.gmt.entry <- function(x) {
    L <- list()
    L[["head"]] <- x[1]
    L[["desc"]] <- x[2]
    l.entry <- x[-c(1:2)]
    idx <- l.entry != ""
    L[["entry"]] <- l.entry[idx]
    L[["len"]] <- length(L[["entry"]])
    return(L)
  }
  L <- lapply(tmp, function(x) mk.gmt.entry(x))
  names(L) <- unlist(lapply(L, function(x) x$head))
  return(L)
}


#' Write a nested list to a GMT file
#'
#' @param lst the nested list to write. See \code{details}.
#' @param fname the desired file name
#'
#' @details \code{lst} needs to be a nested list where each
#'   sub-list is itself a list with the following fields:
#'   - \code{head}: the name of the data
#'   - \code{desc}: description of the corresponding data
#'   - \code{len}: the number of data items
#'   - \code{entry}: a vector of the data items
#'
#' @examples
#' \dontrun{
#' write.gmt(gene_set, "gene_set.gmt")
#' }
#'
#' @family CMap parsing functions
#' @export
write.gmt <- function(lst, fname) {
  # assumes that each element of the list will have the fields
  # head, desc, entry
  if (file.exists(fname)) {
    message(paste(fname, "exists, deleting..."))
    file.remove(fname)
  }
  for (i in 1:length(lst)) {
    el <- lst[[i]]
    ncolumns <- 2 + length(el$entry)
    write(c(el$head, el$desc, el$entry), file=fname, sep="\t", append=T, ncolumns=ncolumns)
  }
}
