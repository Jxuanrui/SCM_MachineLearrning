#Figure5_MachineLearning
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
library(pheatmap)
library(ReactomePA)
library("pathview")
library(biomaRt)
library(glmnet)
library(VennDiagram)
library(sigFeature)
library(e1071)
library(randomForest)
library(caret)
source('SCM_ML/0.Preparing/msvmRFE.R')
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
####2.Preparing data for input####
rm(list = ls())
load("SCM_ML/0.Preparing/expression.Rdata")
####3.Random Forest####
set.seed(123)
ranger_data <- SCM_control_exp
ranger_data <- data.frame(t(ranger_data))
ranger_data$group <- rownames(ranger_data)
rownames(ranger_data) <- NULL
ranger_data <- ranger_data %>% 
  dplyr::select(group,everything())

y <- ifelse(group_list == "Normal Control","normal","trt")
ranger_data$group <- y
ranger_data$group <- as.factor(ranger_data$group)

rf_ntree<- randomForest(group ~ ., data=ranger_data,  
                        ntree=2000,important=TRUE,proximity=TRUE)

importance_value <- data.frame(importance(rf_ntree))
importance_value$symbol <- rownames(importance_value)
importance_value <- importance_value[order(importance_value$MeanDecreaseGini,decreasing = TRUE),]
head(importance_value,20)

x <- importance_value$MeanDecreaseGini[1:20]
names(x) <- importance_value$symbol[1:20]


xrange <- range(pretty(range(x))) 
yrange <- c(1,length(x)) 
par(bty = "o", mgp = c(1.5,.33,0), mar = c(3,7,1,2),las = 1, tcl = -.25)
plot(NULL,NULL,
     xlim = xrange,
     ylim = yrange,
     xlab = "Variable Importance",
     ylab = "",
     yaxt = "n",
     las = 1)
axis(side = 2,at = 1:length(x),rev(names(x))) 

for (i in 1:length(x)) { 
  lines(c(xrange[1],rev(x)[i]),
        c(i,i),
        lwd = 2.5,
        col = "steelblue")
}

####4.Friends analysis####
rm(list = ls())
load("SCM_ML/0.Preparing/DEG.Rdata")

allDiff <- allDiff %>%
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC)>1)

gene <- rownames(allDiff)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene)

mf <- godata('org.Hs.eg.db', ont="MF", computeIC = FALSE)
cc <- godata('org.Hs.eg.db', ont="CC", computeIC = FALSE)

simmf <- mgeneSim(gene$ENTREZID, semData = mf, measure = "Wang", drop = NULL, combine = "BMA")
simcc <- mgeneSim(gene$ENTREZID, semData = cc, measure = "Wang", drop = NULL, combine = "BMA")

fsim <- sqrt(simmf * simcc)

colnames(fsim) = gene$SYMBOL
rownames(fsim) = gene$SYMBOL

for (i in 1:ncol(fsim)){
  fsim[i,i] <- NA
}

y <- melt(fsim) 
y <- y[!is.na(y$value),] 
y <- y[,c(1,3)]

y.mean <- aggregate(.~Var1,y,mean) 
m <- y.mean$value
names(m) <- y.mean$Var1
y$Var1 <- factor(y$Var1, levels=names(sort(m)))

f <- function(y) {
  r <- quantile(y, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  r[3] <- mean(y)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

d_palettes <- palettes_d_names
mycol<-paletteer_d("ggsci::category20_d3",n=20)

p1 <- ggplot(y, aes(Var1, value, fill = factor(Var1))) + 
  scale_fill_manual(values=mycol)+
  guides(fill=FALSE) + 
  
  stat_summary(fun.data= f, geom='boxplot',alpha = 0.8) + 
  geom_hline(aes(yintercept=0.75), linetype="dashed") + 
  
  coord_flip() + 
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) 
####5.SVM-rfe and Lasso####
rm(list = ls())
load("SCM_ML/0.Preparing/expression.Rdata")

#Lasso
data <- SCM_control_exp

data <- data.frame(t(data))

y <- ifelse(group_list == "Normal Control","normal","trt")

x <- data %>% 
  dplyr::mutate(group = y) %>% 
  dplyr::select(group,everything())

y <- ifelse(x$group == "normal", 0,1) 
x <- as.matrix(x[,-1]) 

fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
plot(fit, xvar = "dev", label = TRUE)

cvfit = cv.glmnet(x, y, 
                  nfold=10, 
                  family = "binomial", type.measure = "class")
plot(cvfit)

cvfit$lambda.min 

myCoefs <- coef(cvfit, s="lambda.min")
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea <- lasso_fea[-1]
lasso_fea
#[1] "LTA4H"

#SVM-rfm
data <- SCM_control_exp

data <- data.frame(t(data))

y <- ifelse(group_list == "Normal Control","normal","trt")

x <- data %>% 
  dplyr::mutate(group = y) %>% 
  dplyr::select(group,everything())

x$group <-  factor(x$group)


input <- x

svmRFE(input, k = 10, halve.above = 100) 

nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) 

top.features = WriteFeatures(results, input, save=F) 
head(top.features)

#featsweep = lapply(1:300, FeatSweep.wrap, results, input)
#save(featsweep,file = "featsweep.RData")                       

(load("SCM_ML/0.Preparing/featsweep.RData"))

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
PlotErrors(errors, no.info=no.info)
Plotaccuracy(1-errors,no.info=no.info) 
which.min(errors) 
top.features[1:which.min(errors), "FeatureName"]
#[1] "LTA4H"

par(bty = "o", mgp = c(1.5,.33,0), mar = c(3,7,1,2),las = 1, tcl = -.25)
head(top.features,2)
df <- top.features[1:20,c(1,3)]
df

x <- df$AvgRank[1:20]
names(x) <- df$FeatureName[1:20]
x
xrange <- range(pretty(range(x))) 
yrange <- c(1,length(x)) 

plot(NULL,
     xlim = xrange,
     ylim = yrange,
     xlab = "AvgRank",
     ylab = "",
     yaxt = "n",
     las = 1)
axis(side = 2,at = 1:length(x),rev(names(x))) 

for (i in 1:length(x)) { 
  lines(c(xrange[1],rev(x)[i]),
        c(i,i),
        lwd = 2.5,
        col = "steelblue")
}