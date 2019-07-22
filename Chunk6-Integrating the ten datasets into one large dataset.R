
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 06: Integrating the ten datasets into one large dataset.
### ****************************************************************************

### ------------------------------------------------------------------------ ###
### Step-6.1: Determination of the common genes among five microarray platforms. 

exprs.asthma <- get(load("exprs.asthma.RData"))

exprs2.asthma <- get(load("exprs.asthma2.RData"))

### Step-6.1.1
### Extracting the gene list from five platforms. 

gene.95av2 <- rownames(mat.470)
gene.133a <- rownames(mat.18965)
gene.133p2 <- rownames(mat.4302)
gene.133pm <- rownames(mat.44037)
gene.print3G <- rownames(mat.104468)
gene.4112f <- rownames(mat.63142)


### Step-6.1.2
### Identification of common genes ammong five platforms. 

common.gene <- Reduce(intersect, list(gene.95av2, 
                                      gene.133a, 
                                      gene.133p2, 
                                      gene.133pm, 
                                      gene.print3G, 
                                      gene.4112f))

length(common.gene)

### Step-6.1.3
### Visulization for common genes among five platforms. 

gene.collect <- list(hgu95av2 = gene.95av2, hgu133plus2 = gene.133p2, 
                     hgu133a = gene.133a, hthgu133pm = gene.133pm, 
                     SurePrint3G = gene.print3G, G4112F = gene.4112f)

library(venn)

venn(gene.collect, ilab=TRUE, zcolor = "style")


### ------------------------------------------------------------------------ ###


### ------------------------------------------------------------------------ ###
### Step-6.2: Extracting the common genes and expression matrix. 

### The transfromation from entrez ID to gene symbol.
### please input a numeric vector or matrix of entrez ID.
### The output is a vector with names(gene ID).

library(annotate)
library(org.Hs.eg.db)
id2symbol <- function(x) {
  if (class(x)=="data.frame") print("Error! Please input a numeric vector") 
  else 
  {x <- as.character(x)
  unlist(lookUp(x,'org.Hs.eg','SYMBOL'))
  }
}

com.gene <- id2symbol(common.gene)

# write.csv(com.gene, file = "common_genes.csv")

# com.gene[429]

### Determination of repeat genes in common genes. 

### Find the missing values in common genes. 
com.gene[is.na(com.gene)]
com.gene["1607"]
com.gene["9142"]
com.gene["26148"]

###............................................................................#

#. names(table(com.gene))[table(com.gene) > 1]
#. table(table(com.gene) > 1)
#. length(com.gene) == length(unique(com.gene))

### End.

### ------------------------------------------------------------------------ ###
### Step-6.3: Integrating the ten datasets into one big dataset. 

merge.DS <- cbind(
  
  mat.470[common.gene, ], 
  mat.4302[common.gene, ], 
  mat.18965[common.gene, ], 
  mat.41861[common.gene, ], 
  mat.44037[common.gene, ], 
  mat.63142[common.gene, ], 
  mat.64913[common.gene, ], 
  mat.67472[common.gene, ], 
  mat.89809[common.gene, ], 
  mat.104468[common.gene, ]
  
)

dim(merge.DS)

save(merge.DS, file = "merge.DS.RData")

# merge.DS <- get(load("merge.DS.RData"))

library(ggsci)
mycol <- pal_npg("nrc", alpha = 1)(10)

col.seq <- rep(mycol, time = c(ncol(mat.470), 
                                     ncol(mat.4302), 
                                     ncol(mat.18965), 
                                     ncol(mat.41861), 
                                     ncol(mat.44037), 
                                     ncol(mat.63142), 
                                     ncol(mat.64913), 
                                     ncol(mat.67472), 
                                     ncol(mat.89809), 
                                     ncol(mat.104468)))

batch.seq <- rep(1:10, time = c(ncol(mat.470), 
                                     ncol(mat.4302), 
                                     ncol(mat.18965), 
                                     ncol(mat.41861), 
                                     ncol(mat.44037), 
                                     ncol(mat.63142), 
                                     ncol(mat.64913), 
                                     ncol(mat.67472), 
                                     ncol(mat.89809), 
                                     ncol(mat.104468)))

### ------------------------------------------------------------------------ ###
### Step-6.4: Remove the batch effect from merged dataset. 

### Batch effect removal. 
library(sva)
asthma.DS <- ComBat(dat = merge.DS, batch = batch.seq, mod = NULL, 
                    par.prior = TRUE, prior.plots = FALSE)

# save(asthma.DS.10sets, file = "asthma.DS.10sets.RData")

### Unified the names of all .cel files. 
cel.name <- colnames(asthma.DS)
cel.name <- gsub(".CEL.gz", "", cel.name, ignore.case = TRUE)
# strsplit(cel.name, "_")
cel <- NULL
for (i in cel.name) {
  tmp <- gsub(".CEL.gz", "", i, ignore.case = TRUE)
  cel <- c(cel, strsplit(tmp, "_")[[1]][1])
}
cel

colnames(merge.DS) <- colnames(asthma.DS) <- cel

### Visualization for dataset after batch effect removal.  
### 1) Using boxplot. 

op <- par(mfrow = c(2, 1))
boxplot(merge.DS, col = col.seq)
boxplot(asthma.DS, col = col.seq)
par(op)

### ------------------------------------------------------------------------ ###
### 2) Using tree plot, with the samples names. 

library(gridExtra)
library(grid)
library(ape)
library(ggplot2)
library(ggtree)

groupInfo <- list(
  GSE470 = colnames(mat.470), 
  GSE4302 = colnames(mat.4302),
  GSE18965 = colnames(mat.18965),
  GSE41861 = colnames(mat.41861),
  GSE44037 = colnames(mat.44037),
  GSE63142 = colnames(mat.63142),
  GSE64913 = colnames(mat.64913),
  GSE67472 = colnames(mat.67472),
  GSE89809 = colnames(mat.89809),
  GSE104468 = colnames(mat.104468)
) 
grp <- groupInfo
for (i in 1:length(groupInfo)) {
  x <- gsub(".CEL.gz", "", groupInfo[[i]], ignore.case = TRUE)
  x <- sapply(strsplit(x, "_"), function(x) x[1])
  grp[[i]] <- x
} 
grp

save(groupInfo, file = "groupInfo.RData")

iris.before <- t(merge.DS)
d.before <- dist(iris.before, method = "euclidean")
h.before <- hclust(d.before)
tree.before1 <- as.phylo(h.before)
tree.before2 <- groupOTU(tree.before1, groupInfo)

iris.after <- t(asthma.DS)
d.after <- dist(iris.after, method = "euclidean")
h.after <- hclust(d.after)
tree.after1 <- as.phylo(h.after)
tree.after2 <- groupOTU(tree.after1, groupInfo)

library(ggsci)
mycol <- pal_npg("nrc", alpha = 1)(10)
mycol <- pal_aaas("default", alpha = 1)(10)

p3b <- ggtree(tree.before2, aes(), alpha=.2, layout='circular') + 
  scale_colour_manual(values = mycol) +
  geom_tiplab(size=1, aes(angle=angle, color=group)) + 
  geom_tippoint(aes(color=group), size=2, alpha=.3) + 
  theme(legend.position="right") 

p3a <- ggtree(tree.after2, aes(), alpha=.2, layout='circular') + 
  scale_colour_manual(values = mycol) +
  geom_tiplab(size=1, aes(angle=angle, color=group)) + 
  # theme(legend.position="right") +
  geom_tippoint(aes(color=group), size=2, alpha=.3) 

grid.arrange(p3b, p3a, ncol=2)

### ------------------------------------------------------------------------ ###

### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ###
### 3) Using tree plot, without the samples names. 

library(ggplot2)
library(ggtree)
library(ggnewscale)

mat.label <- select.asthma[, c(1, 3, 4)]
rownames(mat.label) <- mat.label$Sample_ID
mat.label <- mat.label[, c(2, 3)]
mat.label <- mat.label[cel, ]

library(ape)
rownames(iris.before) <- rownames(iris.after) <- cel
tree.before <- hclust(dist(iris.before))
tree.after <- hclust(dist(iris.after))
tree.before <- as.phylo(tree.before)
tree.after <- as.phylo(tree.after)
circ.b <- ggtree(tree.before, branch.length='none', layout = "circular")
circ.a <- ggtree(tree.after, branch.length='none', layout = "circular")

df <- as.data.frame(as.matrix(mat.label)[, 2]); colnames(df) <- "Status"

p.b <- gheatmap(circ.b, df, offset=.1, width=.05,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

p.a <- gheatmap(circ.a, df, offset=.1, width=.05,
                colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")

grid.arrange(p.b, p.a, ncol=2)

### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ###
### ------------------------------------------------------------------------ ###

save(merge.DS, file = "merge.DS_6.12.RData")
save(asthma.DS, file = "asthma.DS_6.12.RData")

### Extract the labels for each samples. 

asthma.label <- select.asthma$Label[match(cel, select.asthma$Sample_ID)]

asthma.label <- as.factor(asthma.label)

### ------------------------------------------------------------------------ ###
### ------------------------------------------------------------------------ ###

