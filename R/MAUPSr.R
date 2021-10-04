#' rMAUPS pipeline - QC, differential analysis, integrative analysis.
#'
#' @docType methods
#' @name MAUPSr
#' @rdname MAUPSr
#'
#' @param metadata File path or data frame of the meta data, including columns of
#' Experiment, Sample, Condition, and multiple comparisons.
#' @param outdir Output directory.
#' @param qc Boolean, specifying whether perform the quanlity control.
#' @param type Could be "msms", "RNASeq", or "Arrary".
#' @param method A character, specifying the method for differential analysis.
#' Optional methods include DESeq2 (RNASeq), GFOLD (RNASeq), limma, glm.pois, glm.qlll, and glm.nb.
#'
#' @author Wubing Zhang
#' @return No return value. Output multiple files into local folder.
#' @export
#'
MAUPSr <- function(metadata, outdir = "./", qc = TRUE,
                   type = "msms", method = "limma"){
  options(stringsAsFactors = FALSE)
  if(!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  message(Sys.time(), " Reading the metadata ...")
  if(length(metadata)==1 && file.exists(metadata)){# Read meta file
    if(grepl(".csv$", metadata))
      metadata = read.csv(metadata, header = TRUE, quote="")
    else
      metadata = read.table(metadata, sep = "\t", header = TRUE,
                            comment.char = "", quote="")
  }

  message(Sys.time(), " Analyzing experiments one by one ...")
  ## Analyze each data separately
  for(experiment in unique(metadata[,1])){
    message("\t", Sys.time(), " ", experiment, " ...")
    rawdata = read.csv(experiment, header = TRUE, row.names = 1)
    meta = metadata[metadata[,1]==experiment, ]
    rownames(meta) = meta[,2]
    meta = meta[, colSums(is.na(meta))<nrow(meta)]
    data = rawdata[, rownames(meta)]
    proj.name = gsub(".*\\/|_normdata.*|\\..*", "", experiment)
    dir.create(file.path(outdir, "qc"), showWarnings = FALSE)
    if(qc){
      plist = ProteomicsQC(data, condition = meta[colnames(data), 3],
                           proj.name = proj.name,
                           outdir = file.path(outdir, "qc"))
      saveRDS(plist, paste0(outdir, "/qc/", proj.name, ".qc.rds"))
    }

    if(sum(is.na(data))>0){
      rowmax = round(max(rowSums(is.na(data))) / ncol(data),2)+0.01
      colmax = round(max(colSums(is.na(data))) / nrow(data),2)+0.01
      data = filterN(data)
      imputeddata = imputeNA(as.matrix(data), rowmax = rowmax,
                            colmax = colmax, k = 5)
      dir.create(file.path(outdir, "imputation"), showWarnings = FALSE)
      tmpfile = paste0(outdir,"/imputation/", basename(gsub("\\..*",
                                        "_imputed.csv", experiment)))
      write.csv(imputeddata, tmpfile, row.names = TRUE, quote = FALSE)
      data = imputeddata
      if(qc){
        plist = ProteomicsQC(data, condition = meta[colnames(data), 3],
                           proj.name = proj.name,
                           outdir = file.path(outdir, "imputation"))
        saveRDS(plist, paste0(outdir, "/imputation/", proj.name, ".qc.rds"))
      }
    }
    comparisons = grep("comparison", colnames(meta),
                       ignore.case = TRUE, value = TRUE)
    for(comp in comparisons){
      dir.create(file.path(outdir, comp), showWarnings = FALSE)
      prefix = paste0(proj.name, ".", comp)
      SA = meta[!is.na(meta[,comp]), comp, drop = FALSE]
      message("\t", Sys.time(), " ", comp, " ...")

      psm.p = NULL
      if(sum(grepl("PSMs", colnames(rawdata)))>0){
        colnames(SA) = "Condition"
        design = model.matrix(~1+Condition, SA)
        rownames(design) = rownames(SA)
        fit = eBayes(lmFit(data[,rownames(design)], design, na.rm=TRUE))
        fit$count = rawdata[rownames(fit$coefficients),
                            grep("PSMs", colnames(rawdata), value = TRUE)[1]]
        fit2 = SpectraCounteBayes(fit)
        psm.p = VarianceBoxplot(fit2,n=30, main=prefix)
        ggsave(paste0(outdir,"/",comp,"/",prefix,"_var_boxplot.png"),
               psm.p, width = 6, height = 5)
      }
      deres_p = DEAnalyze(data, SA, type = type, method = method)
      write.csv(deres_p, paste0(outdir,"/",comp,"/",prefix,"_dep.csv"),
                row.names = TRUE, quote = FALSE)
      deres_p$logP = -log(deres_p$pvalue)
      p1 = ScatterView(deres_p, x = "log2FC", y = "logP",
                       model = "volcano", x_cut = c(-0.2,0.2), force = 5,
                       top = 5, ylab = "-log10(p.value)",
                       main = paste0(prefix, "(Gene)"))
      p1 = p1 + theme(legend.position = "none")
      ggsave(paste0(outdir,"/",comp,"/",prefix, "_dep_volcano.png"),
             p1, width = 6, height = 5)
      res = list()
      res$psm.p = psm.p
      res$dep.p = p1
      res = c(res, DeComplex(deres_p))
      write.csv(res$deComplex, paste0(outdir,"/",comp,"/",prefix, "_dePathway.csv"),
                row.names = TRUE, quote = FALSE)
      res$gobp.p = res$gobp.p + labs(title = paste0(prefix, "(BP)"))
      res$reactome.p = res$reactome.p + labs(title = paste0(prefix, "(REACTOME)"))
      res$gocc.p = res$gocc.p + labs(title = paste0(prefix, "(CC)"))
      res$corum.p = res$corum.p + labs(title = paste0(prefix, "(CORUM)"))
      ggsave(paste0(outdir,"/",comp,"/",prefix, "_deBP_volcano.png"),
             res$gobp.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/",comp,"/",prefix, "_deREACTOME_volcano.png"),
             res$reactome.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/",comp,"/",prefix, "_deCC_volcano.png"),
             res$gocc.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/",comp,"/",prefix, "_deCORUM_volcano.png"),
             res$corum.p, width = 6, height = 5)
      saveRDS(res, paste0(outdir,"/",comp,"/",prefix, ".rds"))
    }
  }
  ## Merge the same comparisons in different experiments
  message(Sys.time(), " Merge the same comparisons in different experiments ...")
  comparisons = grep("comparison", colnames(metadata),
                     ignore.case = TRUE, value = TRUE)
  for(comp in comparisons){
    experiments = unique(metadata[!is.na(metadata[,comp]),1])
    if(length(experiments)>1){
      message("\t", Sys.time(), " ", comp, " ...")
      proj.names = gsub(".*\\/|_normdata.*|\\..*", "", experiments)
      prefix = paste0(proj.names, ".", comp)
      DEPs = paste0(outdir,"/",comp,"/",prefix, "_dep.csv")
      summary = data.frame()
      for(r in DEPs){
        tmp = read.csv(r, row.names = 1, header = TRUE)
        proteins = unique(c(rownames(summary), rownames(tmp)))
        summary = cbind(summary[proteins,], tmp[proteins, c(1,4)])
        rownames(summary) = proteins
      }
      colnames(summary) = paste0(rep(proj.names,each=2), ".",
                                 rep(c("log2FC", "pvalue"),length(DEPs)))
      mergedDep = t(apply(summary, 1, function(x){
        lfc = x[seq(1, length(x),2)]; pval = x[seq(2, length(x),2)]
        lfc = lfc[!is.na(lfc)]; pval = pval[!is.na(pval)]
        if(length(lfc)==0) return(c(NA, NA))
        if(length(lfc)==1) return(c(lfc, pval))
        if(length(lfc)>1) c(sum(lfc)/sqrt(length(lfc)), metap::sumlog(pval)$p)
      }))
      mergedDep = as.data.frame(mergedDep, stringsAsFactors = FALSE)
      colnames(mergedDep) = c("Zscore", "pvalue")
      write.csv(mergedDep, paste0(outdir,"/",comp,"/",comp, "_merged.dep.csv"), quote = FALSE)
      mergedDep$logP = -log10(mergedDep$pvalue)
      p1 = ScatterView(mergedDep, x = "Zscore", y = "logP",
                       model = "volcano", x_cut = c(-0.2,0.2), force = 5,
                       top = 5, ylab = "-log10(p.value)",
                       main = paste0(prefix, "(Gene)"))
      p1 = p1 + theme(legend.position = "none")
      ggsave(paste0(outdir,"/", comp,"/", comp, "_merged_dep_volcano.png"),
             p1, width = 6, height = 5)

      res = list()
      res$dep.p = p1
      res = c(res, DeComplex(mergedDep, lfc = "Zscore"))
      res$gobp.p = res$gobp.p + labs(title = paste0(prefix, "(BP)"))
      res$reactome.p = res$reactome.p + labs(title = paste0(prefix, "(REACTOME)"))
      res$gocc.p = res$gocc.p + labs(title = paste0(prefix, "(CC)"))
      res$corum.p = res$corum.p + labs(title = paste0(prefix, "(CORUM)"))
      write.csv(res$deComplex, paste0(outdir,"/",comp,"/", comp, "_merged.dePathway.csv"),
                row.names = TRUE, quote = FALSE)
      ggsave(paste0(outdir,"/", comp,"/", comp, "_merged_deBP_volcano.png"),
             res$gobp.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/", comp,"/", comp, "_merged_deREACTOME_volcano.png"),
             res$reactome.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/", comp,"/", comp, "_merged_deCC_volcano.png"),
             res$gocc.p, width = 6, height = 5)
      ggsave(paste0(outdir,"/", comp,"/", comp, "_merged_deCORUM_volcano.png"),
             res$corum.p, width = 6, height = 5)
      saveRDS(res, paste0(outdir,"/", comp,"/", comp, "_merged.rds"))
    }
  }
  message(Sys.time(), " Arrange figures for visualization ...")
  pnglist = arrangePngs(outdir)
  saveRDS(pnglist, paste0(outdir,"/", "summary_list.rds"))
}
