#Figure6_Enrichment Analysis
####1.Preparatio####
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(clusterProfiler)
library(GOplot)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
####2.Preparing data for input####
rm(list = ls())
load("SCM_ML/0.Preparing/rat_DEG.Rdata")
diffgene <-  Model_LV_DEG %>%
  filter(FDR < 0.05) %>%
  filter(abs(logFC) >2)

####3.Enrichment Analysis####
#GO
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Rn.eg.db")
ego <- enrichGO(gene          = gene$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable= TRUE)

GO_result <- ego@result
GO_interesting <- GO_result[c(13,102,117,121,122,124,140,146),]
GO_finally <-rbind(GO_result[10:12,],GO_interesting) 

ego@result <- GO_finally
if(T){
  x = ego
  dd =x@result
  dd$richFactor =dd$Count / as.numeric(sub("/\\d+", "", dd$BgRatio))
  dd <- dd[dd$p.adjust < 0.05,]
  
  library(ggplot2)
  ggplot(dd,aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_viridis_c(begin = 0.3, end = 1) +
    scale_size_continuous(range=c(2, 10)) +
    theme_bw() + 
    xlab("Rich factor") +
    ylab(NULL) + 
    ggtitle("")
}

#KEGG
#https://www.genome.jp/kegg/catalog/org_list.html
ekegg <- enrichKEGG(gene          = gene$ENTREZID,
                    organism =    "rno",
                    pAdjustMethod = "BH")
y <- ekegg@result

KEGG_interesting <- y[c(1,11,13,15,162,168,172),]
KEGG_finally <-rbind(y[20:22,],KEGG_interesting) 
ekegg@result <- KEGG_finally

if(T){
  x = ekegg
  dd =x@result
  dd$richFactor =dd$Count / as.numeric(sub("/\\d+", "", dd$BgRatio))
  library(ggplot2)
  ggplot(dd,aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_viridis_c(begin = 0.3, end = 1) +
    scale_size_continuous(range=c(2, 10)) +
    theme_bw() + 
    xlab("Rich factor") +
    ylab(NULL) + 
    ggtitle("")
}


ego@result <- GO_interesting
p3 <- cnetplot(ego, foldChange=gene$ENTREZID, circular = TRUE, colorEdge = TRUE)
