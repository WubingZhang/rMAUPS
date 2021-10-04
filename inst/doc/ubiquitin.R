## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----libs, warning=FALSE---------------------------------------------------
library(rMAUPS)

## ----listfile--------------------------------------------------------------
uniprot_path <- system.file("extdata", "human_uniprot_seq.txt", package = "rMAUPS")
uniprot <- read.delim(uniprot_path, sep='\t')
head(uniprot)

## ----search----------------------------------------------------------------
lysine_pos <- searchProtSeq(head(uniprot), 'K')
head(lysine_pos)

## ----viewProt--------------------------------------------------------------
# get hits for P01116
id <- 'P01116'
start <- lysine_pos[lysine_pos['UniprotId']==id, 'start']
end <- lysine_pos[lysine_pos['UniprotId']==id, 'end']

# view on protein structure
browseProtStructure(id, start, end)

## ----viewProtUb------------------------------------------------------------
# Here, we manually create a dataframe for the KRAS ubiqutination sites.
# In the realisitic scenario, you should load the "Ubiquitination_site_dataset" file
# from PhosphositePlus DB using the readPSPUbiquitin function. You'll need to download
# it your self due to copy protection.
id <- 'P01116'
ub <- data.frame(gene='KRAS', ID=id, position=c(117, 128, 147))

# get ub sites for P01116
start <- ub[ub['ID']==id, 'position']
end <- ub[ub['ID']==id, 'position']

# view on protein structure
browseProtStructure(id, start, end)

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

