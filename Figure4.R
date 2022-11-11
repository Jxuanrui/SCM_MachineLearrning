#Figure4_Immune infiltration
####1.Preparation####
library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(curl)
library(devtools)
library(MCPcounter)
library(pheatmap)
library(vioplot)  
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
source("SCM_ML/0.Preparing/Fig4/CIBERSORT.R")
####2.CIBERSORT####
load("SCM_ML/0.Preparing/expression.Rdata")
#write.table(cbind(rownames(exprSet), exprSet),"exprSet.txt",quote = F, sep = "\t", row.names=FALSE)
sig_matrix <- "SCM_ML/0.Preparing/Fig4/LM22.txt"
mixture_file = 'SCM_ML/0.Preparing/Fig4/exprSet.txt'
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=10, QN=TRUE)
res_cibersort <- data.frame(res_cibersort)
res_cibersort <- res_cibersort[,-c(23,24,25)]
res_cibersort <- res_cibersort[,apply(res_cibersort,2,sum) > 0]
colnames(res_cibersort)<-paste(colnames(res_cibersort),"CIBERSORT",sep="_")
####3.MCPcounter####
genes <- data.table::fread("SCM_ML/0.Preparing/Fig4/genes.txt",data.table = F)
probesets <- data.table::fread("SCM_ML/0.Preparing/Fig4/probesets.txt",data.table = F,header = F)
example <- read.table("SCM_ML/0.Preparing/Fig4/exprSet.txt", header=T, sep="\t", check.names=F,row.names = 1)
results<- MCPcounter.estimate(example,
                              featuresType= "HUGO_symbols",
                              probesets=probesets,
                              genes=genes)
results <- data.frame(t(results))
colnames(results)<-paste(colnames(results),"MCPCOUNTER",sep="_")

####4.Figure4A####
hmdat <- cbind(res_cibersort,results)
methods.col <- brewer.pal(n = 2,name = "Paired")
annCol <- data.frame(group = group_list,
                     row.names = rownames(results),
                     stringsAsFactors = F)
annRow <- data.frame(Methods = factor(c(rep("CIBERSORT",length(grep("_CIBERSORT",colnames(hmdat)))),rep("MCPCOUNTER",length(grep("_MCPCOUNTER",colnames(hmdat))))),levels = c("CIBERSORT","MCPCOUNTER")),
                     row.names = colnames(hmdat),
                     stringsAsFactors = F)
annColors <- list(Methods = c("CIBERSORT" = methods.col[1], 
                              "MCPCOUNTER" = methods.col[2]))
indata <- t(hmdat)
indata <- indata[,colSums(indata) > 0]
plotdata <- standarize.fun(indata,halfwidth = 2)

pheatmap::pheatmap(mat = as.matrix(plotdata), 
                   border_color = NA, 
                   color = bluered(64), 
                   cluster_rows = F, 
                   cluster_cols = F, 
                   show_rownames = T, 
                   show_colnames = F, 
                   annotation_col = annCol, 
                   annotation_row = annRow, 
                   annotation_colors = annColors, 
                   gaps_col = table(annCol$group)[1], 
                   gaps_row = cumsum(table(annRow$Methods)), 
                   cellwidth = 5, 
                   cellheight = 10, 
                   filename = "immune heatmap by pheatmap.pdf")

####5.Figure4B####
normal=6                                                          
tumor=15
res_cibersort
tmp<-strsplit(colnames(res_cibersort),split="_",fixed=TRUE)
colnames(res_cibersort) <- unlist(lapply(tmp,head,1)) 
rt <- res_cibersort

pdf("vioplot.pdf",height=8,width=13)             
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,58),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
for(i in 1:ncol(rt)){
  normalData=rt[1:normal,i]
  tumorData=rt[(normal+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02,labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)), cex = 0.8)
  text(seq(1,60,3),-0.03,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
dev.off()
