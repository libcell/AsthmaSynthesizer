
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 01: Obtain the informations of datasets on asthma, from GEO.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-01. Preparing the GEO accession number for GEO Series.

# load("Rimage_6.11.RData")

# My Mac OS. 
# setwd("/Users/libo/Test/Asthma/AsthmaData/")

# My windows.
setwd("F:/Asthma/AsthmaData/test")

source("http://bioconductor.org/biocLite.R")

gse.asthma <- c("GSE470", "GSE4302", "GSE18965", "GSE41861", "GSE44037", 
                "GSE64913", "GSE67472", "GSE89809", "GSE104468", "GSE63142")

library(GEOquery)

GSEInfo <- list()

p <- 0

for (s in gse.asthma) {
  
  gpl <- getGEO(s)
  
  # GSM samples. 
  
  gsm <- gpl[[1]]$geo_accession
  
  p <- p + 1
  
  GSEInfo[[p]] <- gsm
  
}

names(GSEInfo) <- gse.asthma

save(GSEInfo, file = "GSEInfo.RData")

# GSEInfo <- get(load("GSEInfo.RData"))
# 
# The End.