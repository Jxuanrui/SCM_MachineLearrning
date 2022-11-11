#Figure7_homology mapping_Enrichment Analysis
####1.Preparation####
library(tidyverse)
library(ggplot2)
library(magrittr)
library(ggrepel)
library(janitor)
library(cowplot)
library(clusterProfiler)
library(enrichplot)
library(GEOquery)
library(limma)
library(AnnoProbe)
library("pathview")
library(biomaRt)
library(randomForest)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)
####2.Preparing data for input####
#Diagnostic markers obtained from machine learning screening
marker
####3.Homology mapping####
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
Rat <- useMart('ensembl',dataset = "rnorvegicus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")

#h2m.g <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
#                values =marker ,mart = human,
#                attributesL = c("rgd_symbol"),
#                martL = Rat,uniqueRows = T)
h2m.g <- getLDS(attributes = c("rgd_symbol"),filters = "rgd_symbol",
                values =marker ,mart = Rat,
                attributesL = c("hgnc_symbol"),
                martL = human,uniqueRows = T)

h2m <- h2m.g[!duplicated(h2m.g$RGD.symbol),]
x <- data.frame(symbol = h2m$HGNC.symbol)
finall <- merge(x,Model_LV_DEG,by = "symbol")

####4.
####4.Enrichment Analysis####
load("SCM_ML/0.Preparing/expression.Rdata")
df <- SCM_control_exp
data <- as.data.frame(t(df))

y <- as.numeric(data[,"SEMA6B"]) 
colnames <- colnames(data)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  print(i)
  test <- cor.test(as.numeric(data[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}

names(cor_data_df) <- c("Symbol","correlation","pvalue")

cor_data_sig <- cor_data_df %>% 
  filter(pvalue < 0.05) %>% 
  filter(abs(correlation)>0.6) %>% 
  arrange(desc(abs(correlation)))

gene <- cor_data_sig$Symbol
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene          = gene$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable= TRUE)
barplot(ego, showCategory=10) 

eKEGG <- enrichKEGG(gene    = gene$ENTREZID,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.1)
barplot(eKEGG, showCategory=15) 

df <- SCM_control_exp
data <- as.data.frame(t(df))

y <- as.numeric(data[,"PLAT"]) 
colnames <- colnames(data)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  print(i)
  test <- cor.test(as.numeric(data[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}

names(cor_data_df) <- c("Symbol","correlation","pvalue")

cor_data_sig <- cor_data_df %>% 
  filter(pvalue < 0.05) %>% 
  filter(abs(correlation)>0.6) %>% 
  arrange(desc(abs(correlation)))

gene <- cor_data_sig$Symbol
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene <- cbind(cor_data_sig,gene)
identical(gene$Symbol,gene$SYMBOL)

gene$correlation <- sort(gene$correlation,decreasing = T)
geneList <- gene[,2]
names(geneList) = as.character(gene[,1])

go_result <-ego@result 

go_interesting <- go_result[c(12,19,),]
KEGG_finally <-rbind(y[20:22,],KEGG_interesting) 
ekegg@result <- KEGG_finally

eKEGG <- enrichKEGG(gene    = gene$ENTREZID,
                    pAdjustMethod = "BH")
barplot(eKEGG) 
#Showing one of the groups