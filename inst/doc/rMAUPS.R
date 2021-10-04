## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval = FALSE---------------------------------------------------------
#  $ conda install -c r r r-markdown r-recommended
#  $ conda install -c anaconda libgfortran libnetcdf libxml2
#  $ conda install -c conda-forge pandoc r-data.table r-rcpparmadillo

## ----installr, eval=FALSE--------------------------------------------------
#  install.packages(c("devtools", "BiocManager"), repos = "https://cloud.r-project.org")
#  # Install dependencies
#  BiocManager::install(c("ggpubr", "metap", "ggrepel", "GSVA", "DESeq2", "limma", "impute", "biomaRt", "msigdbr", "BiocStyle", "msmsTests"))
#  # Install rMAUPS from github
#  devtools::install_github("WubingZhang/rMAUPS")

## ----libs, warning=FALSE---------------------------------------------------
library(ggplot2)
library(rMAUPS)

## ----listfile--------------------------------------------------------------
datapath = system.file("extdata", package = "rMAUPS")
list.files(datapath, pattern = "export_proteins")

## ----preprocess------------------------------------------------------------
normalizeProteomeDiscoverer(datapath, output = "./", log2 = TRUE)

## ----onedata---------------------------------------------------------------
normdata = normalizeProteomeDiscoverer(file.path(datapath, "experiment1_export_proteins.txt"), log2 = TRUE, return = TRUE)
head(normdata)

## ----readData--------------------------------------------------------------
metadata = read.csv(file.path(datapath, "metadata.csv"))
head(metadata)

## ----pipeline, eval=FALSE--------------------------------------------------
#  MAUPSr(metadata, outdir = "analysis/")
#  ## Or
#  MAUPSr(system.file("extdata", "metadata.csv", package = "rMAUPS"), outdir = "analysis/")

## ----visualize, eval=FALSE-------------------------------------------------
#  view(outdir = "analysis/")

## ----simulatedata----------------------------------------------------------
data = as.matrix(normdata[,-1])
meta = metadata[metadata$Experiment=="experiment1_normdata.csv", -1]
rownames(meta) = meta[,1]
simulated = data
idx = sample(1:length(simulated), round(0.1*length(simulated)))
simulated[idx] = NA

## ----qc--------------------------------------------------------------------
qc = ProteomicsQC(simulated, condition = meta[colnames(data), 2], proj.name = "TestQC")
qc$p1
qc$p2
qc$p4
qc$p5
qc$p6
qc$p7

## ----normalize-------------------------------------------------------------
normalized = normalizeProteomics(simulated, norm = "median", log2 = FALSE)

## ----knn-------------------------------------------------------------------
imputed = imputeNA(normalized)
plot(imputed[idx], data[idx])

## ----dep-------------------------------------------------------------------
deres = DEAnalyze(data, meta[,-1], type = "msms", method = "limma")
## Visualize the results
deres$logFDR = log10(deres$padj)
ScatterView(deres, x = "log2FC", y = "logFDR", 
            x_cut = c(-0.5,0.5), y_cut = -2, 
            groups = c("bottomleft", "bottomright"), top = 5)

## ----decomplex-------------------------------------------------------------
res = DeComplex(deres)
head(res$deComplex)
res$gobp.p
res$reactome.p
res$gocc.p
res$corum.p

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

