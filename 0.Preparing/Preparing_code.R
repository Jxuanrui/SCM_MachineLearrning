#SCM_Machine learning
####1.Preparation####
library(tidyverse)
library(ggplot2)
library(magrittr)
library(ggrepel)
library(janitor)
library(cowplot)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(GEOquery)
library(limma)
library(AnnoProbe)
library(biomaRt)
library(glmnet)
library(VennDiagram)
library(sigFeature)
library(e1071)
library(randomForest)
library(caret)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)
####2.Preparing data for input####
gset = getGEO('GSE95368', destdir=".",getGPL = F)

gset=gset[[1]]

pdata=pData(gset)

SCM_control_pData <- pdata[c(1:6,19:33),]
SCM_control_pData$characteristics_ch1.2[4] <- "group: Normal Control"
group_list <- SCM_control_pData$characteristics_ch1.2
group_list <- unlist(strsplit(group_list,": "))[seq(2,42,2)]
group_list=factor(group_list,levels = c("Normal Control","Acute SCM"))
group_list

exprSet=exprs(gset)
SCM_control_exp <- exprSet[,match(rownames(SCM_control_pData),colnames(exprSet))]
boxplot(SCM_control_exp,outline=FALSE, notch=T,col=group_list, las=2)
SCM_control_exp=normalizeBetweenArrays(SCM_control_exp)

ex <- SCM_control_exp
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
SCM_control_exp <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

probe2symbol <- data.table::fread("SCM_ML/0.Preparing/platform.txt",data.table = FALSE)
probe2symbol <- probe2symbol[,c(1,5)]
colnames(probe2symbol) <- c("probeset","symbol")

SCM_control_exp <- as.data.frame(SCM_control_exp)

SCM_control_exp <- SCM_control_exp %>% 
  rownames_to_column(var="probeset") %>% 
  inner_join(probe2symbol,by="probeset") %>% 
  dplyr::select(-probeset) %>% 
  dplyr::select(symbol,everything()) %>% 
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>% 
  filter(symbol != "NA") %>% 
  arrange(desc(rowMean)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-rowMean) %>% 
  column_to_rownames(var = "symbol")

SCM_control_exp <- SCM_control_exp[-grep(" ",rownames(SCM_control_exp)),]

#save(SCM_control_exp,file = "SCM_ML/0.Preparing/expression.Rdata")

design=model.matrix(~ group_list)
colnames(design) <- levels(group_list)
design
fit=lmFit(SCM_control_exp,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf)

diffgene <- allDiff %>%
  filter(adj.P.Val < 0.05) %>%
  filter(abs(logFC) >1)

#save(allDiff,diffgene,group_list,file = "SCM_ML/0.Preparing/DEG.Rdata")
