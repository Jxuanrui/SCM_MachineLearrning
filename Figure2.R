#Figure2A_Volcano plot
####1.Preparation####
library(ggplot2)
library(tidyverse)
library(ggrepel)
####2.Preparing data for input####
rm(list = ls())
load("SCM_ML/0.Preparing/DEG.Rdata")

####3.Figure2A####
x <- allDiff
x$label <- rownames(allDiff)

plot_mode <- "advanced" 

logFCcut <- 1 #log2-foldchange
pvalCut <- 0.05 #P.value
FDR <- 0.05 #adj.P.value

logFCcut2 <- 2
logFCcut3 <- 5
pvalCut2 <- 0.0001
pvalCut3 <- 0.00001

xmin <- (range(x$logFC)[1]- (range(x$logFC)[1]+ 10))
xmax <- (range(x$logFC)[1]+ (10-range(x$logFC)[1]))
ymin <- 0
ymax <- -log10(x$P.Value)[3] * 1.1

mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

if (plot_mode == "classic"){
  x$color_transparent <- ifelse((x$PValue < pvalCut & x$logFC > logFCcut), "red", ifelse((x$PValue < pvalCut & x$logFC < -logFCcut), "blue","grey30"))
  size <- ifelse((x$PValue < pvalCut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  n1 <- length(x[, 1])
  cols <- rep("grey30", n1)
  names(cols)<- rownames(x)
  
  cols[x$P.Value< pvalCut & x$logFC >logFCcut]<- "#FB9A99"
  cols[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- "#ED4F4F"
  cols[x$P.Value < pvalCut & x$logFC < -logFCcut]<- "#B2DF8A"
  cols[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- "#329E3F"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  size[x$P.Value < pvalCut & x$logFC > logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC > logFCcut3]<- 6
  size[x$P.Value < pvalCut & x$logFC < -logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC < -logFCcut3]<- 6
  
} else {
  stop("Unsupport mode")
}

p1 <- ggplot(data=x, aes(logFC, -log10(P.Value), label = label)) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(-10, -5, -logFCcut, 0, logFCcut, 5, 10), 
    labels = c(-10, -5, -logFCcut, 0, logFCcut, 5, 10),
    limits = c(-11, 11) 
  ) +
  #xlim(c(xmin, xmax)) + 
  
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + 
  geom_hline(yintercept = -log10(pvalCut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12#, base_family = "Times" 
  ) +
  theme(panel.grid=element_blank())

if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut2), color="grey40", 
               linetype="longdash", lwd = 0.5)
}
p1

n = 2
p1 + geom_text_repel(aes(x = logFC, y = -log10(P.Value), 
                         label = ifelse(abs(logFC) > n, rownames(x),"")),
                     colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines"))

selectgenes <- x %>%
  filter(P.Value < 0.05) %>%
  filter(abs(logFC) >2)


p2 <- p1 + 
  geom_point(data = selectgenes, alpha = 1, size = 3, shape = 1, 
             stroke = 1, 
             color = "black") +
  
  scale_color_manual(values = mycol) + 
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, 
                  size = 5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines"),
                  direction = "y") +
  guides(color=guide_legend(title = NULL)) 