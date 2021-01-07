

rm(list = ls())
library(limma)
library(tidyverse)
library(ggplot2)
library(ReactomePA)
library(org.Hs.eg.db)
library(readxl)
## setwd
setwd("C:/Users/fbinsatter/Box/Brain paper 2")
df = read_xlsx("C:\\Users\\fbinsatter\\Box\\Brain paper 2\\Brain2_Supp_tables_v1.xlsx",
               sheet = 5, skip = 1)
df = as.data.frame(df)
colnames(df) ; row.names(df)

rownames(df) = df[ , 1]
df = df[ , - 1]

## genelist
glist = read_xlsx("C:\\Users\\fbinsatter\\Box\\Brain paper 2\\Brain2_Supp_tables_v1.xlsx",
                  sheet = 2, skip = 1)
glist2 = glist$`Gene Symbol`

## pheno
pheno = read_xlsx("C:\\Users\\fbinsatter\\Box\\Brain paper 2\\Brain2_Supp_tables_v1.xlsx",
                  sheet = 4, skip = 1)
pheno = pheno[match(rownames(df), pheno$Nanostring.or.submitted.ID), ]
match(rownames(df), pheno$Nanostring.or.submitted.ID)

# 168 gene expression
df2 = df[ , colnames(df) %in% glist2]
match(rownames(df2), pheno$Nanostring.or.submitted.ID)
df3 = data.frame(pheno$Sup168, df2)
df3 = df3 %>% filter(pheno.Sup168 != "Ambi")

### limma
design = model.matrix( ~0 + df3$pheno.Sup168)
colnames(design) = c("TP1", "TP2", "TP3", "TP4")
fit = lmFit(t(df3[ , -1]), design)
contrast.matrix = makeContrasts(
  "TP1" = TP1 - (TP2 + TP3 + TP4)/3,
  "TP2" = TP2 - (TP1 + TP3 + TP4)/3,
  "TP3" = TP3 - (TP1 + TP2 + TP4)/3,
  "TP4" = TP4 - (TP1 + TP2 + TP3)/3, levels = design
)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)
allToptable = topTable(fit2, n = dim(df3)[2 -1], adjust.method = "fdr")
TP1.tt = topTable(fit2, coef = "TP1", n = dim(df3)[2]-1, adjust.method = "fdr")
TP2.tt = topTable(fit2, coef = "TP2", n = dim(df3)[2]-1, adjust.method = "fdr")
TP3.tt = topTable(fit2, coef = "TP3", n = dim(df3)[2]-1, adjust.method = "fdr")
TP4.tt = topTable(fit2, coef = "TP4", n = dim(df3)[2]-1, adjust.method = "fdr")


## id conversion
library(org.Hs.eg.db)
library(ReactomePA)
geneID = mapIds(org.Hs.eg.db, rownames(TP1.tt), "ENTREZID", "SYMBOL")
fchange = TP1.tt$logFC
x = enrichPathway(geneID, organism = "human", pAdjustMethod = "BH", minGSSize = 10, readable = T)
p1 = cnetplot(x, foldChange = abs(fchange), categorySize = "pvalue", colorEdge = T) 


### GSEA
library(fgsea)
pathways = gmtPathways("C:\\Users\\fbinsatter\\Box\\R\\Renal cancer\\Surv_analysis\\c2.cp.reactome.v7.2.entrez.gmt")
fgene = function(x){
  genes = rownames(x)
  entrezID = mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
  x1 = data.frame(entrezID, x)
  x1 = x1[!duplicated(x1$entrezID), ]
  x1 = x1[!is.na(x1$entrezID), ]
  ranks = x1$t
  names(ranks) = x1$entrezID
  fgsea.x = fgsea(pathways = pathways,
                  stats = ranks,
                  minSize = 15,
                  eps = 0,
                  maxSize = 500)
  fgsea.x = fgsea.x[order(fgsea.x$padj, decreasing = F), ]
  return(fgsea.x)
}
franks = function(x){
  genes = rownames(x)
  entrezID = mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
  x1 = data.frame(entrezID, x)
  x1 = x1[!duplicated(x1$entrezID), ]
  x1 = x1[!is.na(x1$entrezID), ]
  ranks = x1$t
  names(ranks) = x1$entrezID
  return(ranks)
}
fgsea.TP1 = fgene(TP1.tt)
fgsea.TP2 = fgene(TP2.tt)
fgsea.TP3 = fgene(TP3.tt)
fgsea.TP4 = fgene(TP4.tt)


### barplot

pdf("GSEA_by_Group.pdf", width = 12, height = 6)
ggplot(fgsea.TP1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() + labs(x = "Pathway", y = "Normalized Enrichment Score",
                      title = "Reactome TP1") +
  theme_minimal()

ggplot(fgsea.TP2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() + labs(x = "Pathway", y = "Normalized Enrichment Score",
                      title = "Reactome TP2a") +
  theme_minimal()
ggplot(fgsea.TP3, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() + labs(x = "Pathway", y = "Normalized Enrichment Score",
                      title = "Reactome TP2b") +
  theme_minimal()
ggplot(fgsea.TP4, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() + labs(x = "Pathway", y = "Normalized Enrichment Score",
                      title = "Reactome TP4") +
  theme_minimal()
dev.off()
