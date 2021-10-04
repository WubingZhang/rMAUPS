#' Quality control of proteomic data
#'
#' @param data A data frame of the proteomics data.
#' @param condition A vector of conditions corresponding to
#' the samples in the proteomics data.
#' @param proj.name A character as the prefix of output files.
#' @param outdir The path to a local directory.
#'
#' @return a list of ggplot objects.
#' @import ggplot2
#' @author Wubing Zhang
#' @export
#'
ProteomicsQC <- function(data, condition = NULL, proj.name = NA, outdir = NULL){
  ## QCs associated with missing values
  missflag = FALSE
  if(sum(is.na(data))>0){
    missflag = TRUE
    gg = data.frame(gene = rownames(data), NAs = rowSums(is.na(data)))
    p4 = DensityView(gg[,2,drop=FALSE], xlab = "The number of missing value")
    p4 = p4 + theme(legend.position = "none") + labs(title = proj.name)
    gg = data.frame(sample = colnames(data), Detection = colSums(!is.na(data)))
    p5 = DensityView(gg[,2,drop=FALSE], xlab = "The number of detected gene")
    p5 = p5 + theme(legend.position = "none") + labs(title = proj.name)
    p6 = BarView(gg, "sample", "Detection", fill = "#8da0cb",
                 ylab = "The number of detected gene", main = proj.name)
    p6 = p6 + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
    p7 = countNA(data) + labs(title = proj.name)

    if(!is.null(outdir)){
      ggsave(paste0(outdir, "/", proj.name, "_density_NA_acrossGene.png"), p4, width = 4, height = 3.5)
      ggsave(paste0(outdir, "/", proj.name, "_density_detection.png"), p5, width = 4, height = 3.5)
      ggsave(paste0(outdir, "/", proj.name, "_bar_detection.png"), p6, width = 4, height = 3.5)
      ggsave(paste0(outdir, "/", proj.name, "_count_NA.png"), p7, width = 4, height = 3.5)
    }
  }else{
    p4 = p5 = p6 = p7 = NULL
  }
  if(missflag) {
    data = filterN(data, minS = ncol(data)*0.5)
    data[is.na(data)] = median(data, na.rm = TRUE)
  }
  ## Consistency between samples
  p1 = ViolinView(data, ylab = "Logarithmic abundance", main = proj.name)
  p1 = p1 + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
  p2 = pcView(data, color = condition) + labs(title = proj.name)
  ## Protein correlation within protein complexes
  p3 = corComplex(data) + labs(title = proj.name)
  if(!is.null(outdir)){
    ggsave(paste0(outdir, "/", proj.name, "_violin_abundance.png"),
           p1, width = 4, height = 3.5)
    ggsave(paste0(outdir, "/", proj.name, "_pcview_samples.png"),
           p2, width = 4, height = 3.5)
    ggsave(paste0(outdir, "/", proj.name, "_correlation_complexes.png"),
           p3, width = 4, height = 3.5)
  }
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, p7 = p7))
}
