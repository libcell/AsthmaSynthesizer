
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 02: Download the datasets, asthma.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-01. Download these ten datasets on asthma from GEO.

library(GEOquery)

gse.download <- gse.asthma

gse_non <- NULL

for (i in gse.download) {
  
  test <- try(getGEOSuppFiles(i, makeDirectory = TRUE, baseDir = getwd(),
                              fetch_files = TRUE, filter_regex = NULL), 
              silent=TRUE)
  
  if (class(test) == "NULL") {
    
    gse_non <- c(gse_non, i)
    
    next
    
  }
  
}

gse_downloaded <- dir()

gse_non

save(gse_non, file = "gse_non.RData")

### End of Step-01.
### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###