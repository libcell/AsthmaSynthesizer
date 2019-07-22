
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 08: Visualization of differentially expressed genes. 
### ****************************************************************************

DEGs.1.5 <- DEGs[DEGs$FC_grade == ">1.5-fold", c(1, 8)]

#. DEGs.1.5 <- DEGs[order(DEGs$logFC)[c(1:50, 104:153)], c(1, 8)]
# DEGs.1.5 <- DEGs[, c(1, 8)]

degs1.5.DS <- asthma.DS[rownames(DEGs.1.5), ]

degs1.5.DS <- t(degs1.5.DS)

tmp <- gsub(".CEL.GZ", "", rownames(degs1.5.DS), ignore.case = TRUE)

tmp <- strsplit(tmp, "_")

tmp <- unlist(sapply(tmp, function(x) x[1]))

rownames(degs1.5.DS) <- tmp

#. heatmap(degs1.5.DS)

library(ggplot2)
library(ggtree)
library(ggnewscale)

# Create a matrix, with samples in rows and features in columns.

iris.mat <- degs1.5.DS

tree.iris <- hclust(dist(iris.mat, method = "euclidean"))

library(ape)

tree <- as.phylo(tree.iris)

circ <- ggtree(tree, branch.length='none', layout = "circular")

circ <- groupOTU(circ, grp, 'Dataset') + 
   aes(color=Dataset) + 
   theme(legend.position="right")

circ

df <- as.data.frame(as.matrix(mat.label)); colnames(df) <- c("Status", "Dataset")

df2 <- iris.mat[, abs(DEGs.1.5$logFC) > log2(1.5)]

# abs(DEGs.1.5$logFC) > log2(1.5) 

dim(df2)

p1 <- gheatmap(circ, df, offset=0, width=0.3,
               colnames_angle = 95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option = "D", name = "discrete\nvalue")

p1

p2 <- p1 + new_scale_fill()

gheatmap(p2, df2, offset = 2, width=1, low = "gray",
         high = "black", color = "white", 
         colnames_angle = 90, colnames_offset_y = .25) + 
  scale_fill_viridis_c(option="A", name="continuous\nvalue")

### ************************************************************************ ###
### ************************************************************************ ###


### Construct the correlation network!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

### Positional gene sets. 

library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
hasProteinData(edb)

#. select(edb, keys = "ZBTB16", keytype = "GENENAME", columns = "UNIPROTID")
#. select(edb, keys = "ZBTB16", keytype = "GENENAME", columns = "ENTREZID")

#. AllY.granges.tx <- genes(edb,
#.                          filter = SeqNameFilter("Y"),
#.                          columns = c("gene_id", "seq_name",
#.                                      "seq_strand", "tx_id", "tx_biotype",
#.                                      "tx_seq_start", "tx_seq_end"),
#.                          order.by = "tx_seq_start")

All.granges <- genes(edb, return.type = "DataFrame")

rownames(All.granges) <- All.granges$entrezid

head(All.granges)

All.granges[All.granges$entrezid == "100033816", ]


# degs <- NULL

for (g in rownames(DEGs)) {
  
  print("#####################################################################")
  
  x <- unlist(xx[g])
  
  print(x)
  
  Sys.sleep(2)
  # namex(x) <- NULL
  
  # degs <- c(degs, x)
  
}

degs <- DEGs$Symbol

degs.posi <- NULL

for (i in degs) {
  
  x <- All.granges[All.granges$symbol == i, c(2, 4, 5)][1, ]
  
  degs.posi <- rbind(degs.posi, x)

}

head(degs.posi)

degs.posi[degs.posi$gene_name == "GATA2", ] <- All.granges[All.granges$symbol == "GATA2", c(2, 4, 5)][2, ]
degs.posi[degs.posi$gene_name == "SF1", ] <- All.granges[All.granges$symbol == "SF1", c(2, 4, 5)][2, ]

degs.posi[degs.posi$gene_name == "PROS1", 2] <- 93873033
degs.posi[degs.posi$gene_name == "PROS1", 3] <- 93974089

degs.posi[degs.posi$gene_name == "LRIG1", 2] <- 66378797
degs.posi[degs.posi$gene_name == "LRIG1", 3] <- 66500932

degs.posi[degs.posi$gene_name == "MUC2", 2] <- 1102455
degs.posi[degs.posi$gene_name == "MUC2", 3] <- 1103456

degs.posi[degs.posi$gene_name == "ADAM9", 2] <- 38996986
degs.posi[degs.posi$gene_name == "ADAM9", 3] <- 39105001

dim(degs.posi)
asthma_degs.pos <- degs.posi
#. asthma_degs.pos <- data.frame(chr = paste0("chr", unlist(xx[rownames(DEGs)])), 
#.                         start = degs.posi$gene_seq_start, 
#.                         end = degs.posi$gene_seq_end) 

#. rownames(asthma_degs.pos) <- rownames(DEGs)

#. BiocManager::install("karyoploteR")

library(karyoploteR)
kp <- plotKaryotype(genome="hg38", plot.type=2, chromosomes=c("chr3", "chr4", "chr18"))

poi.data <- openxlsx::read.xlsx("Position_enrichment-6.4.xlsx", 1)

poi.data <- poi.data[-c(1:11), ]

colnames(poi.data) <- poi.data[1, ]

poi.data <- poi.data[-1, ]

DT::datatable(poi.data)

position.data <- poi.data[1:22, ]


kpDataBackground(kp, r1=0.45, data.panel=1)

kpAddBaseNumbers(kp)

kpPoints(kp, chr = c(rep("chr3", 8), rep("chr4", 6), rep("chr18", 5), rep("chr3", 3)), 
         x = as.numeric(position.data$Start), 
         y = runif(22, 0.05, 0.40), 
         col = position.data$Regulation)


kpText(kp, chr = c(rep("chr3", 8), rep("chr4", 6), rep("chr18", 5), rep("chr3", 3)), 
       x = as.numeric(position.data$Start), 
       y = runif(22, 0.5, 1.5), 
       labels = position.data$`Gene Symbol`)

asthma_hyper <- asthma_degs.pos[DEGs$Status == "Up", ]
asthma_hypo <- asthma_degs.pos[DEGs$Status == "Down", ]

library(circlize)
circos.initializeWithIdeogram()
bed_list = list(asthma_hyper, asthma_hypo)
head(asthma_hyper$chr)
length(asthma_hypo$chr)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))

circos.genomicDensity(asthma_hyper, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(asthma_hypo, col = c("#0000FF80"), track.height = 0.1)

circos.genomicDensity(asthma_degs.pos, col = c("green"), track.height = 0.1)

degs.posi["13", ]

################################################################################
save.image("Rimage_6.15.RData")
################################################################################


### ****************************************************************************
### code chunk number 08: Identification of asthma-related drugs or compounds. 
### ****************************************************************************

#. biocLite("PharmacoGx")
#. biocLite("hgu133a.db")

## download and process the HDAC signature

library(PharmacoGx)

HDAC.asthma <- DEGs[DEGs$FC_grade == ">1.5-fold", c(1, 8)]

# Top10 and bottom10. 
#. HDAC.asthma <- DEGs[order(DEGs$logFC), c(1, 8)]
#. HDAC.asthma <- HDAC.asthma[c(1:10, (nrow(HDAC.asthma)-9):nrow(HDAC.asthma)), ]


### ---------------------- Convert entrez gene IDs to ENSEMBL IDs. --------- ###

library(org.Hs.eg.db)
#. help(package = "org.Hs.eg.db")

## Bimap interface:
x <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the Ensembl gene IDs for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

en.seq <- NULL
for (i in rownames(HDAC.asthma)) {
  # print(xx[[i]])
  en.seq <- c(en.seq, xx[[i]][1])
  # Sys.sleep(1)
  
}

rownames(HDAC.asthma) <- en.seq
HDAC.asthma <- HDAC.asthma[, 2]
names(HDAC.asthma) <- en.seq
HDAC.asthma[HDAC.asthma == "Up"] <- 1
HDAC.asthma[HDAC.asthma == "Down"] <- -1
HDAC.asthma <- as.numeric(HDAC.asthma)
names(HDAC.asthma) <- en.seq

### ------------------------------------------------------------------------ ###

#. drug.perturbation <- PharmacoGx::downloadPertSig("CMAP")

# drug.perturbation <- get(load("drug.perturbation.RData"))

message("Be aware that computing sensitivity will take some time...")

f.cmap <- function(x, HDAC) { 
  return(PharmacoGx::connectivityScore(x=x, y=HDAC, method="gsea", nperm=1000))
}

res <- apply(drug.perturbation[ , , c("tstat", "fdr")], 2, f.cmap, HDAC = HDAC.asthma)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
DT::datatable(res)
res <- as.data.frame(res)
save(res, file = "result.RData")

drug <- res[res$`P Value` < 0.05, ]

drug <- drug[order(drug$`P Value`), ]

drug$`P Value` <- round(drug$`P Value`, digits = 4)

#### GSEA by using ClusterProfiler

library(clusterProfiler)
library(DOSE)
asthma.geneList <- DEGs$logFC
names(asthma.geneList) <- rownames(DEGs)

asthma.geneList <- sort(asthma.geneList, decreasing = TRUE)

library(DOSE)

asthma.geneList

gene <- names(asthma.geneList)

library(enrichplot)
edo <- enrichDGN(gene)

# barplot(edo, showCategory=5)

## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

heatplot(edox, foldChange=asthma.geneList)

### End of all codes. 
