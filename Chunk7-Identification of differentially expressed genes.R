
################################################################################
### Meta-analysis of asthma-related microarray data in GEO, Feb 26th, 2019.
### Author: Bo Li, Xiner Nie.
################################################################################

### ****************************************************************************
### code chunk number 07: Identification of differentially expressed genes. 
### ****************************************************************************

library(limma)

#  Moderated t-statistic

asthma.group <- abs(as.numeric(asthma.label) - 2)
asthma.degisn <- cbind(Intercept = 1, Group = asthma.group)
fit <- lmFit(asthma.DS, asthma.degisn)
fit <- eBayes(fit)
Gene.detect <-topTable(fit, coef = 2, n = nrow(asthma.DS))

dim(Gene.detect)

### Vocalno plot 

p.col <- NULL
p.pch <- NULL

for (i in 1:nrow(Gene.detect)) {
  
  if (abs(Gene.detect$logFC[i]) < log2(1.2)) {p.col[i] <- "black"; p.pch[i] <- 1} else 
    if (Gene.detect$logFC[i] < -log2(1.2)) {p.col[i] <- "green"; p.pch[i] <- 16} else 
      {p.col[i] <- "red"; p.pch[i] <- 16}
  
} 

x <- Gene.detect$logFC; y <- -1*log10(Gene.detect$adj.P.Val)

plot(x, y, 
     xlab = "Log2(Fold Change)", ylab = "-log10(adj.P.Val)", 
     col = p.col, pch = p.pch)

abline(v = c(-log2(1.2), log2(1.2)), lty = 3, col = 4, lwd = 1.5)
abline(h = c(-log10(0.05)), lty = 3, col = 4, lwd = 1.5)

for (i in 1:nrow(Gene.detect)) {
  
  if (abs(x[i]) > log2(1.5)) {
    
    text(x[i] + 0.05, y[i] + 0.5, rownames(Gene.detect)[i])
    
    Sys.sleep(0.1)
    
    } else next

}


# DT::datatable(Gene.detect) # FC more than 1.2

DEGs <- Gene.detect[abs(Gene.detect$logFC) > log2(1.2) & Gene.detect$adj.P.Val < 0.05, ]

status <- rep(NA, nrow(DEGs))
status[DEGs$logFC > 0] <- "Up"
status[DEGs$logFC < 0] <- "Down"

FC_grade <- rep(NA, nrow(DEGs))
FC_grade[abs(DEGs$logFC) > log2(1.5)] <- ">1.5-fold"
FC_grade[abs(DEGs$logFC) < log2(1.5)] <- ">1.2-fold"

DEGs$Status <- status

DEGs$FC_grade <- FC_grade

DEGs$Symbol <- com.gene[rownames(DEGs)]

DT::datatable(DEGs)

library(openxlsx)
# write.xlsx(DEGs, file = "DEGs_limma_Combat_10_1.2fold.xlsx", asTable = TRUE, row.names = TRUE)

common_gene <- data.frame(Entrez_ID = common.gene, Gene_Symbol = com.gene)

length(common.gene)

# write.xlsx(common_gene, file = "common_genes_among_10_datasets.xlsx", row.names = TRUE)

length(com.gene)
gene.deg <- rownames(DEGs)
gene.up <- rownames(DEGs)[DEGs$logFC > 0]
gene.down <- rownames(DEGs)[DEGs$logFC < 0]
# write.csv(gene.deg, "deg.csv")

#. id2symbol(gene.up)
#. id2symbol(gene.down)

# save.image("matrix_practice.RData")

