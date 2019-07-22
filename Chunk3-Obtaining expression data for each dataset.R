
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 03: Obtaining expression data for each dataset.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-01. Counting the sample data, with GSM accession number.

library(openxlsx)

select.asthma <- read.xlsx("Supplementary Material S1-All datasets and samples primarily selected in this study.xlsx", 
                           sheet = 1, 
                           rowNames = FALSE)

cel.list <- select.asthma$Sample_ID

### End of Step-01. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. Reading, and normalizing all raw data for each dataset, firtstly.  

library(GEOquery)
library(affy)
library(limma)

file.dir <- dir()

dir.pos <- grep("^GSE[0-9]", file.dir)

file.dir <- file.dir[dir.pos]

file.dir <- file.dir[order(as.numeric(gsub("GSE", "", file.dir)))]


# Unzip all compressed files into the current folder

f <- file.dir

for (d in f[-1]) {
  
  setwd(d)
  
  rd <- dir()
  
  untar(rd)
  
  file.remove(rd)
  
  setwd("..")
  
}

# The decompression process is complete

eset.asthma <- list()

obj.asthma <- list()

p <- 0

for (f in file.dir) {
  
  setwd(f)
  
  # rd <- dir()
  # untar(rd)
  # file.remove(rd)
  
  p <- p + 1
  
  if (f != "GSE104468"  & f != "GSE63142") {
    
    dat <- ReadAffy()
    
    eset <- rma(dat)
    
    eset2 <- exprs(eset)
    
    obj.asthma[[p]] <- eset
    
  } else {
    
    #. dat <- read.maimages(files = dir(), columns = list(
    #.   G = "gMedianSignal", Gb = "gBGMedianSignal", 
    #.   R = "gMedianSignal", Rb = "gBGMedianSignal", 
    #.   annotation = c("ProbeName")
    #. ))
    
    dat <- read.maimages(files = dir(),
                         source="agilent", 
                         columns = list(G = "gMedianSignal", 
                                        Gb = "gBGMedianSignal", 
                                        R = "gMedianSignal", 
                                        Rb = "gBGMedianSignal"), 
                         annotation = c("Row", 
                                        "Col",
                                        "FeatureNum",
                                        "ControlType",
                                        "ProbeName",
                                        "SystematicName"))
    
    dat2 <- backgroundCorrect(dat, "normexp", offset = 50)
    eset <- normalizeBetweenArrays(dat2$R, method = "quantile")
    eset2 <- log(eset)
    
    rownames(eset2) <- dat$genes$ProbeName
    
    asdf <- ExpressionSet(assayData = eset)
    
    obj.asthma[[p]] <- asdf
    
  }
  
  eset.asthma[[p]] <- eset2
  
  setwd("..")
  
}

names(eset.asthma) <- names(obj.asthma) <- file.dir

save(eset.asthma, file = "eset.asthma-10.RData")
save(obj.asthma, file = "obj.asthma-10.RData")

### End of Step-02.
### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###

