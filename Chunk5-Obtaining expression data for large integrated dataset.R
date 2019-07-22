
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 05: Obtaining expression data for large dataset.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-5.1: Prepare the annotation informations for all datasets.

### !!! Warnning!!! -------------------------------------------------------- ### 
### Analyzing data from line 384. 

err.samples <- get(load("err.samples.RData"))

exprs.asthma <- get(load("exprs.asthma.RData"))

### Prepareing the annotation informations for five types of microarrays. 

hgu95av2 <- read.csv("GPL8300-hgu95av2.txt", sep = "\t", header = TRUE)
hgu133a <- read.csv("GPL96-hgu133a.txt", sep = "\t", header = TRUE)
hgu133plus2 <- read.csv("GPL570-hgu133plus2.txt", sep = "\t", header = TRUE)
hthgu133pm <- read.csv("GPL13158-hthgu133pm.txt", sep = "\t", header = TRUE)
SurePrintG3 <- read.csv("GPL21185-Agilent-SurePrintG3.txt", sep = "\t", header = TRUE)
G4112F <- read.csv("GPL6480-Agilent-G4112F.txt", sep = "\t", header = TRUE)

dim(hgu133a)
dim(hgu95av2)
dim(hgu133plus2)
dim(hthgu133pm)
dim(SurePrintG3)
dim(G4112F)

head(hthgu133pm)


### End of Step-5.1 
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-5.2: Prepare gene expression levels for large dataset before integration.

eset.tidy <- list()

### ------------------------------------------------------------------------ ###
### Step-5.2.1: For hgu95av2, annotation informations, GSE470.

### Step-5.2.1.1
### Find out the probes which match multiple genes per probe, or are unknown. 

o2m.p.95av2 <- NULL

len.stat <- NULL

for (i in 1:nrow(hgu95av2)) {
  
  a <- as.character(hgu95av2$ENTREZ_GENE_ID[i])
  
  tmp <- strsplit(a, " /// ")[[1]]
  
  # adding a condition, "| is.na(tmp)"
  
  len.stat <- c(len.stat, length(tmp))
  
  if (length(tmp) > 1 | length(tmp) == 0) {
    
    o2m.p.95av2 <- c(o2m.p.95av2, i)
    
  }   
  
}

o2m.p.95av2
length(o2m.p.95av2)

# sum(table(len.stat))-11433

### Step-5.2.1.2
### Remove the o2m probes. 

anno.95av2 <- hgu95av2[-o2m.p.95av2, c(1, 11, 12)]

table(is.na(anno.95av2$ENTREZ_GENE_ID))

# DT::datatable(anno.95av2)

table(unique(anno.95av2$ID) == anno.95av2$ID)

# one probe, one gene. 

o2o.95av2 <-  names(table(anno.95av2$ENTREZ_GENE_ID))[table(anno.95av2$ENTREZ_GENE_ID) == 1]
length(o2o.95av2)

# more probes, one gene. 

m2o.95av2 <-  names(table(anno.95av2$ENTREZ_GENE_ID))[table(anno.95av2$ENTREZ_GENE_ID) > 1]
length(m2o.95av2)


### Step-5.2.1.3
### Extracting the gene expression values from GSE470 dataset.  

### 1) Extracting one2one gene expression levels. 

mat.470 <- NULL

probe.470 <- rownames(exprs.asthma$GSE470)

for (s in o2o.95av2) {
  
  s.pos <- which(anno.95av2$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE470[which(probe.470 == anno.95av2$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.470 <- rbind(mat.470, s.tmp)
  
}

### 2) Extracting more2one gene expression levels. 

for (m in m2o.95av2) {
  
  m.pos <- which(anno.95av2$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE470[match(anno.95av2$ID[m.pos], probe.470), ] 
  
  if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.470 <- rbind(mat.470, tmp)
  
}

rownames(mat.470) <- c(o2o.95av2, m2o.95av2)

head(mat.470)
tail(mat.470)
dim(mat.470)

eset.tidy$mat.470 <- mat.470

### End of Step-5.2.1
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### Step-5.2.2 For hgu133plus2, annotation informations, GSE4302, GSE41861, GSE64913, GSE67472.

### Srep-5.2.2.1
### Find out the probes which match multiple genes per probe, or are unknown.  

o2m.p.133p2 <- NULL

len.stat <- NULL
# b <- 0

for (i in 1:nrow(hgu133plus2)) {
  
  a <- as.character(hgu133plus2$ENTREZ_GENE_ID[i])
  
  tmp <- strsplit(a, " /// ")[[1]]
  
  len.stat <- c(len.stat, length(tmp))
  
  if (length(tmp) > 1 | length(tmp) == 0) {
    
    o2m.p.133p2 <- c(o2m.p.133p2, i)
    
    # b <- b + 1
    
  }   
  
}

length(o2m.p.133p2)
# b
# sum(table(len.stat)) - 41834

table(is.na(hgu133plus2$ENTREZ_GENE_ID))

# which(o2m.p.95av2 == "4721")

### Srep-5.2.2.1
### Remove the o2m probes. 

anno.133p2 <- hgu133plus2[-o2m.p.133p2, c(1, 11, 12)]

# DT::datatable(anno.133p2)

table(unique(anno.133p2$ID) == anno.133p2$ID)

# one probe, one gene. 

o2o.133p2 <-  names(table(anno.133p2$ENTREZ_GENE_ID))[table(anno.133p2$ENTREZ_GENE_ID) == 1]
length(o2o.133p2)

# more probes, one gene. 

m2o.133p2 <-  names(table(anno.133p2$ENTREZ_GENE_ID))[table(anno.133p2$ENTREZ_GENE_ID) > 1]
length(m2o.133p2)

### Srep-5.2.2.3
### Extracting the gene expression values from GSE4302, GSE41861, GSE64913, GSE67472.  

### =================================== GSE4302 ================================= ###

### 1) Extracting one2one gene expression levels. 

mat.4302 <- NULL

probe.4302 <- rownames(exprs.asthma$GSE4302)

for (s in o2o.133p2) {
  
  s.pos <- which(anno.133p2$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE4302[which(probe.4302 == anno.133p2$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.4302 <- rbind(mat.4302, s.tmp)
  
}

### 2) Extracting more2one gene expression levels. 

for (m in m2o.133p2) {
  
  m.pos <- which(anno.133p2$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE4302[match(anno.133p2$ID[m.pos], probe.4302), ] 
  
  # if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.4302 <- rbind(mat.4302, tmp)
  
}

rownames(mat.4302) <- c(o2o.133p2, m2o.133p2)

head(mat.4302)
tail(mat.4302)
dim(mat.4302)

eset.tidy$mat.4302 <- mat.4302

### ======================================================================== ###

### ================================ GSE41861 ============================== ###

### 1) Extracting one2one gene expression levels. 

mat.41861 <- NULL

probe.41861 <- rownames(exprs.asthma$GSE41861)

for (s in o2o.133p2) {
  
  s.pos <- which(anno.133p2$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE41861[which(probe.41861 == anno.133p2$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.41861 <- rbind(mat.41861, s.tmp)
  
}

### 2) Extracting more2one gene expression levels. 

for (m in m2o.133p2) {
  
  m.pos <- which(anno.133p2$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE41861[match(anno.133p2$ID[m.pos], probe.41861), ] 
  
  # if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.41861 <- rbind(mat.41861, tmp)
  
}

rownames(mat.41861) <- c(o2o.133p2, m2o.133p2)

head(mat.41861)
tail(mat.41861)
dim(mat.41861)

eset.tidy$mat.41861 <- mat.41861

### ======================================================================== ###

### ================================ GSE64913 ============================== ###

### 1) Extracting one2one gene expression levels. 

mat.64913 <- NULL

probe.64913 <- rownames(exprs.asthma$GSE64913)

for (s in o2o.133p2) {
  
  s.pos <- which(anno.133p2$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE64913[which(probe.64913 == anno.133p2$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.64913 <- rbind(mat.64913, s.tmp)
  
}

### 2) Extracting more2one gene expression levels. 

for (m in m2o.133p2) {
  
  m.pos <- which(anno.133p2$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE64913[match(anno.133p2$ID[m.pos], probe.64913), ] 
  
  # if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.64913 <- rbind(mat.64913, tmp)
  
}

rownames(mat.64913) <- c(o2o.133p2, m2o.133p2)

head(mat.64913)
tail(mat.64913)
dim(mat.64913)

eset.tidy$mat.64913 <- mat.64913

### ======================================================================== ###

### ================================ GSE67472 ============================== ###

### 1) Extracting one2one gene expression levels. 

mat.67472 <- NULL

probe.67472 <- rownames(exprs.asthma$GSE67472)

for (s in o2o.133p2) {
  
  s.pos <- which(anno.133p2$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE67472[which(probe.67472 == anno.133p2$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.67472 <- rbind(mat.67472, s.tmp)
  
}

### 2) Extracting more2one gene expression levels. 

for (m in m2o.133p2) {
  
  m.pos <- which(anno.133p2$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE67472[match(anno.133p2$ID[m.pos], probe.67472), ] 
  
  # if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.67472 <- rbind(mat.67472, tmp)
  
}

rownames(mat.67472) <- c(o2o.133p2, m2o.133p2)

head(mat.67472)
tail(mat.67472)
dim(mat.67472)

eset.tidy$mat.67472 <- mat.67472

### ======================================================================== ###


### ------------------------------------------------------------------------ ###
### Step-5.2.3 For hgu133a, annotation informations, GSE18965.

### Step-5.2.3.1
### Find out the probes which match multiple genes per probe, or are unknown.  

o2m.p.133a <- NULL

# b <- 0
len.stat <- NULL

for (i in 1:nrow(hgu133a)) {
  
  a <- as.character(hgu133a$ENTREZ_GENE_ID[i])
  
  tmp <- strsplit(a, " /// ")[[1]]
  
  len.stat <- c(len.stat, length(tmp))
  
  if (length(tmp) > 1 | length(tmp) == 0) {
    
    o2m.p.133a <- c(o2m.p.133a, i)
    
    # b <- b + 1
    
  }   
  
}

length(o2m.p.133a)
# b

# sum(table(len.stat)) - 19726

# which(o2m.p.95av2 == "4721")

### Step-5.2.3.2 
### Remove the o2m probes.

anno.133a <- hgu133a[-o2m.p.133a, c(1, 11, 12)]

# DT::datatable(anno.95av2)

table(unique(anno.133a$ID) == anno.133a$ID)

# one probe, one gene. 

o2o.133a <-  names(table(anno.133a$ENTREZ_GENE_ID))[table(anno.133a$ENTREZ_GENE_ID) == 1]
length(o2o.133a)

# more probes, one gene. 

m2o.133a <-  names(table(anno.133a$ENTREZ_GENE_ID))[table(anno.133a$ENTREZ_GENE_ID) > 1]
length(m2o.133a)

### Step-5.2.3.3
### Extracting the gene expression values from GSE470 dataset.  

### 1) Extracting one2one gene expression levels. 

mat.18965 <- NULL

probe.18965 <- rownames(exprs.asthma$GSE18965)

for (s in o2o.133a) {
  
  s.pos <- which(anno.133a$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE18965[which(probe.18965 == anno.133a$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.18965 <- rbind(mat.18965, s.tmp)
  
}

### 2) Extracting more2one gene expression levels. 

for (m in m2o.133a) {
  
  m.pos <- which(anno.133a$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE18965[match(anno.133a$ID[m.pos], probe.18965), ] 
  
  # if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.18965 <- rbind(mat.18965, tmp)
  
}

rownames(mat.18965) <- c(o2o.133a, m2o.133a)

head(mat.18965)
tail(mat.18965)
dim(mat.18965)

eset.tidy$mat.18965 <- mat.18965

### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-5.2.4 hthgu133pm, annotation informations, GSE44037, GSE51392, GSE89809.

# Find out the probes which match multiple genes per probe, or are unknown.  
o2m.p.133pm <- NULL

# b <- 0

len.stat <- NULL

for (i in 1:nrow(hthgu133pm)) {
  
  a <- as.character(hthgu133pm$ENTREZ_GENE_ID[i])
  
  tmp <- strsplit(a, " /// ")[[1]]
  
  len.stat <- c(len.stat, length(tmp))
  
  if (length(tmp) > 1 | length(tmp) == 0) {
    
    o2m.p.133pm <- c(o2m.p.133pm, i)
    
    # b <- b + 1
    
  }   
  
}

length(o2m.p.133pm)
# b

# sum(table(len.stat)) - 40740

# which(o2m.p.95av2 == "4721")

# Remove the o2m probes. 

anno.133pm <- hthgu133pm[-o2m.p.133pm, c(1, 11, 12)]

# DT::datatable(anno.133p2)

table(unique(anno.133pm$ID) == anno.133pm$ID)

# one probe, one gene. 

o2o.133pm <-  names(table(anno.133pm$ENTREZ_GENE_ID))[table(anno.133pm$ENTREZ_GENE_ID) == 1]
length(o2o.133pm)

# more probes, one gene. 

m2o.133pm <-  names(table(anno.133pm$ENTREZ_GENE_ID))[table(anno.133pm$ENTREZ_GENE_ID) > 1]
length(m2o.133pm)

### Extracting the gene expression values from GSE44037, GSE89809.  
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

### ================================ GSE44037 ============================== ###

# GSE44037, GSE89809

# 1) Extracting one2one gene expression levels. 

mat.44037 <- NULL

probe.44037 <- rownames(exprs.asthma$GSE44037)

for (s in o2o.133pm) {
  
  s.pos <- which(anno.133pm$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE44037[which(probe.44037 == anno.133pm$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.44037 <- rbind(mat.44037, s.tmp)
  
}

# 2) Extracting more2one gene expression levels. 

for (m in m2o.133pm) {
  
  m.pos <- which(anno.133pm$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE44037[match(anno.133pm$ID[m.pos], probe.44037), ] 
  
  # if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.44037 <- rbind(mat.44037, tmp)
  
}

rownames(mat.44037) <- c(o2o.133pm, m2o.133pm)

head(mat.44037)
tail(mat.44037)
dim(mat.44037)

eset.tidy$mat.44037 <- mat.44037

### ======================================================================== ###


### =================================== GSE89809 ================================= ###

# GSE44037, GSE51392, GSE89809

# 1) Extracting one2one gene expression levels. 

mat.89809 <- NULL

probe.89809 <- rownames(exprs.asthma$GSE89809)

for (s in o2o.133pm) {
  
  s.pos <- which(anno.133pm$ENTREZ_GENE_ID == s) 
  
  if (length(s.pos) > 1) print(s.pos)
  
  s.tmp <- exprs.asthma$GSE89809[which(probe.89809 == anno.133pm$ID[s.pos]), ] 
  
  # rownames(mat.tmp) <- s
  
  mat.89809 <- rbind(mat.89809, s.tmp)
  
}

# 2) Extracting more2one gene expression levels. 

for (m in m2o.133pm) {
  
  m.pos <- which(anno.133pm$ENTREZ_GENE_ID == m) 
  
  m.tmp <- exprs.asthma$GSE89809[match(anno.133pm$ID[m.pos], probe.89809), ] 
  
  # if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.89809 <- rbind(mat.89809, tmp)
  
}

rownames(mat.89809) <- c(o2o.133pm, m2o.133pm)

head(mat.89809)
tail(mat.89809)
dim(mat.89809)

eset.tidy$mat.89809 <- mat.89809



### ------------------------------------------------------------------------ ###
### Step-5.2.5 For Agilent G4112F, annotation informations, GSE63142.

### Step-5.2.5.1
### Find out the probes which match multiple genes per probe, or are unknown. 

o2m.p.4112f <- NULL

# b <- 0

len.stat <- NULL

for (i in 1:nrow(G4112F)) {
  
  a <- as.character(G4112F$GENE[i])
  
  tmp <- strsplit(a, " /// ")[[1]]
  
  len.stat <- c(len.stat, length(tmp))
  
  if (length(tmp) > 1 | length(tmp) == 0 | is.na(tmp)) {
    
    o2m.p.4112f <- c(o2m.p.4112f, i)
    
    # b <- b + 1
    
  }   
  
  # print(i); print(tmp); print(o2m.p.PrintG3)
  
}

length(o2m.p.4112f)

table(len.stat)

# b

# which(o2m.p.95av2 == "4721")

### Step-5.2.5.2
### Remove the o2m probes. 

###!!!!!!!!!!!!!!!!!!!!!!
anno.4112f <- G4112F[-o2m.p.4112f, c(1, 6, 7)]
length(is.na(anno.4112f$GENE))
table(is.na(anno.4112f$GENE))
# DT::datatable(anno.95av2)

table(unique(anno.4112f$ID) == anno.4112f$ID)

# one probe, one gene. 

o2o.4112f <-  names(table(anno.4112f$GENE))[table(anno.4112f$GENE) == 1]
length(o2o.4112f)

# more probes, one gene. 

m2o.4112f <-  names(table(anno.4112f$GENE))[table(anno.4112f$GENE) > 1]
length(m2o.4112f)

### Step-5.2.5.3
### Extracting the gene expression values from GSE470 dataset.  

### 1) Extracting one2one gene expression levels. 

mat.63142 <- NULL

probe.63142 <- rownames(exprs.asthma$GSE63142)

count.o2o <- NULL

for (s in o2o.4112f) { # 16513 genes. 
  
  s.pos <- which(anno.4112f$GENE == s) 
  
  # if (length(s.pos) > 1) print(s.pos)
  
  probe.pos <- which(probe.63142 == anno.4112f$ID[s.pos])
  
  # print(probe.pos)
  
  s.tmp <- exprs.asthma$GSE63142[probe.pos, ] 
  
  if (length(probe.pos) > 1) {
    
    s.tmp <- apply(s.tmp, 2, mean)
    
  }
  
  # rownames(mat.tmp) <- s
  
  # if (length(which(probe.104468 == anno.PrintG3$ID[s.pos])) > 1) {print(s.tmp)}
  
  mat.63142 <- rbind(mat.63142, s.tmp)
  
}

dim(mat.63142)

### 2) Extracting more2one gene expression levels. 

for (m in m2o.4112f) { # 8531 genes. 
  
  m.pos <- which(anno.4112f$GENE == m) 
  
  m.tmp <- exprs.asthma$GSE63142[match(anno.4112f$ID[m.pos], probe.63142), ] 
  
  if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.63142 <- rbind(mat.63142, tmp)
  
}

rownames(mat.63142) <- c(o2o.4112f, m2o.4112f)

head(mat.63142)
tail(mat.63142)
dim(mat.63142)

eset.tidy$mat.63142 <- mat.63142



### ------------------------------------------------------------------------ ###
### Step-5.2.6 For Agilent SurePrint3G, annotation informations, GSE104468.

### Step-5.2.6.1
### Find out the probes which match multiple genes per probe, or are unknown. 

o2m.p.PrintG3 <- NULL

# b <- 0

len.stat <- NULL

for (i in 1:nrow(SurePrintG3)) {
  
  a <- as.character(SurePrintG3$LOCUSLINK_ID[i])
  
  tmp <- strsplit(a, " /// ")[[1]]
  
  len.stat <- c(len.stat, length(tmp))
  
  if (length(tmp) > 1 | length(tmp) == 0 | is.na(tmp)) {
    
    o2m.p.PrintG3 <- c(o2m.p.PrintG3, i)
    
    # b <- b + 1
    
  }   
  
  # print(i); print(tmp); print(o2m.p.PrintG3)
  
}

length(o2m.p.PrintG3)

table(len.stat)

# b

# which(o2m.p.95av2 == "4721")

### Step-5.2.6.2
### Remove the o2m probes. 

###!!!!!!!!!!!!!!!!!!!!!!
anno.PrintG3 <- SurePrintG3[-o2m.p.PrintG3, c(1, 5, 6)]
length(is.na(anno.PrintG3$LOCUSLINK_ID))
table(is.na(anno.PrintG3$LOCUSLINK_ID))
# DT::datatable(anno.95av2)

table(unique(anno.PrintG3$ID) == anno.PrintG3$ID)

# one probe, one gene. 

o2o.PrintG3 <-  names(table(anno.PrintG3$LOCUSLINK_ID))[table(anno.PrintG3$LOCUSLINK_ID) == 1]
length(o2o.PrintG3)

# more probes, one gene. 

m2o.PrintG3 <-  names(table(anno.PrintG3$LOCUSLINK_ID))[table(anno.PrintG3$LOCUSLINK_ID) > 1]
length(m2o.PrintG3)

### Step-5.2.6.3
### Extracting the gene expression values from GSE470 dataset.  

### 1) Extracting one2one gene expression levels. 

mat.104468 <- NULL

probe.104468 <- rownames(exprs.asthma$GSE104468)

count.o2o <- NULL

for (s in o2o.PrintG3) { # 16513 genes. 
  
  s.pos <- which(anno.PrintG3$LOCUSLINK_ID == s) 
  
  # if (length(s.pos) > 1) print(s.pos)
  
  probe.pos <- which(probe.104468 == anno.PrintG3$ID[s.pos])
  
  # print(probe.pos)
  
  s.tmp <- exprs.asthma$GSE104468[probe.pos, ] 
  
  if (length(probe.pos) > 1) {
    
    s.tmp <- apply(s.tmp, 2, mean)
    
  }
  
  # rownames(mat.tmp) <- s
  
  # if (length(which(probe.104468 == anno.PrintG3$ID[s.pos])) > 1) {print(s.tmp)}
  
  mat.104468 <- rbind(mat.104468, s.tmp)
  
}

dim(mat.104468)

### 2) Extracting more2one gene expression levels. 

for (m in m2o.PrintG3) { # 8531 genes. 
  
  m.pos <- which(anno.PrintG3$LOCUSLINK_ID == m) 
  
  m.tmp <- exprs.asthma$GSE104468[match(anno.PrintG3$ID[m.pos], probe.104468), ] 
  
  if (nrow(m.tmp) > 1) print(m.tmp)
  
  tmp <- m.tmp[which.max(apply(m.tmp, 1, IQR)), ]
  
  mat.104468 <- rbind(mat.104468, tmp)
  
}

rownames(mat.104468) <- c(o2o.PrintG3, m2o.PrintG3)

head(mat.104468)
tail(mat.104468)
dim(mat.104468)

eset.tidy$mat.104468 <- mat.104468

# save(eset.tidy, file = "eset.tidy_10_sets.RData")
### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###