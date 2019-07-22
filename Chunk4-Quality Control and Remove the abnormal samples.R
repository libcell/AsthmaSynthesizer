
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 04: Quality Control and Remove the abnormal samples.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-01. Remove the unwanted cel-format files, according to object 'cel.list'.

for (i in file.dir) {
  
  setwd(i)
  
  pos.seq <- NULL
  
  for (s in cel.list) {
    
    pos <- grep(s, dir())
    
    pos.seq <- c(pos.seq, pos)
    
  }
  
  file.remove(dir()[-pos.seq])
  
  setwd("..")
  
}


# eset.asthma <- get(load("eset.asthma.RData"))

for (i in 1:length(eset.asthma)) {
  
  print(names(eset.asthma)[i])
  
  print(dim(eset.asthma[[i]]))
  
}

eset.asthma

### End of Step-01. 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-02. Carry out the QC check, and remove the incorrect samples.

# -- QC --

# obj.asthma <- get(load("obj.asthma.RData"))

length(obj.asthma)

err.samples <- NULL

for (d in 1:length(obj.asthma)) {
  
  obj <- obj.asthma[[d]]
  
  dir.nam <- paste("QC_report_for", names(obj.asthma)[d], sep = "_")
  
  library(arrayQualityMetrics)
  err.pos <- arrayQualityMetrics(expressionset = obj, 
                      outdir = dir.nam, 
                      force = TRUE)
  
  err.cel <- which(err.pos$arrayTable == "x", arr.ind = TRUE)[, 1]
  
  err.sam <- err.pos$arrayTable$sampleNames[as.numeric(names(table(err.cel))[table(err.cel) > 0])]
  
  err.samples <- c(err.samples, err.sam)
  
}

# -- End --

# save(err.samples, file = "err.samples.RData")

# -- Remove --

for (e in file.dir) {
  
  setwd(e)
  
  rem.pos <- na.omit(match(err.samples, dir()))
  
  rem.file <- dir()[rem.pos]
  
  print(rem.file)
  
  file.remove(rem.file)
  
  setwd("..")
  
}

### End of Step-02. 
### ------------------------------------------------------------------------ ###

### Step-03. Re-reading and normalizing the raw data after QC check.

exprs.asthma <- list()

p <- 0

for (f in file.dir) {
  
  setwd(f)
  
  # rd <- dir()
  # untar(rd)
  # file.remove(rd)
  
  p <- p + 1
  
  if (f != "GSE104468" & f != "GSE63142") {
    
    dat <- ReadAffy()
    
    eset <- rma(dat)
    
    eset2 <- exprs(eset)
    
    # obj.asthma[[p]] <- eset
    
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
    
    # obj.asthma[[p]] <- asdf
    
  }
  
  exprs.asthma[[p]] <- eset2
  
  setwd("..")
  
}

names(exprs.asthma) <- file.dir

# ncol(exprs.asthma$GSE104468)


# save(exprs.asthma, file = "exprs.asthma-10_sets-June_8.RData")

# save(obj.asthma, file = "obj.asthma.RData")

# exprs.asthma <- get(load("exprs.asthma.RData"))

for (s in 1:10) {
  
  print(dim(exprs.asthma[[s]]))
  
}

### End of Step-04. 
### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###
