#Figure3_Enrichment Analysis
####1.Preparatio####
library(clusterProfiler)
library(GOplot)
library(tidyverse)
library(data.table)
library(ggraph)
library(tidygraph)
source(file = "SCM_ML/0.Preparing/gather_graph_node.R")
source(file = "SCM_ML/0.Preparing/gather_graph_edge.R")
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
####2.Preparing data for input####
load("SCM_ML/0.Preparing/DEG.Rdata")
####3.Figure3A&B####
#GO
rm(list=ls())
load("SCM_ML/0.Preparing/DEG.Rdata")

diffgene <- allDiff %>%
  filter(adj.P.Val < 0.05) %>%
  filter(abs(logFC) >0)
gene <- rownames(diffgene)

gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ego <- enrichGO(gene          = gene$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable= TRUE)

GO_result <- ego@result

GO_interesting <- GO_result[c(58,59,68,99,98,117),]
GO_finally <-rbind(GO_result[13:16,],GO_interesting) 

mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

colnames(GO_finally) 
p3 <- ggplot(data = GO_finally,
             aes(x = Count,
                 y = reorder(Description,Count)))+ 
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = 1) +
  labs(x = "",
       y = "",
       title = "",
       size = "Count") +
  mytheme
p3

#KEGG
rm(list=ls())
load("SCM_ML/0.Preparing/DEG.Rdata")

gene <- rownames(allDiff)

gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eKEGG <- enrichKEGG(gene          = gene$ENTREZID,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2)

KEGG_result <- eKEGG @result

KEGG_interesting <- KEGG_result[c(10,39,100,139,142,143),]
KEGG_finally <-rbind(KEGG_result[1:4,],KEGG_interesting) 

mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11))

colnames(KEGG_finally) 
p3 <- ggplot(data = KEGG_finally,
             aes(x = Count,
                 y = reorder(Description,Count)))+ 
  geom_point(aes(size = Count,
                 color = -log10(pvalue)))+
  theme_bw()+
  scale_color_distiller(palette = "Spectral",direction = 1) +
  labs(x = "",
       y = "",
       title = "",
       size = "Count") +
  mytheme
####4.Figure3C####
df <- allDiff
df$SYMBOL <- rownames(df)
gsym.fc <- df %>% 
  select(SYMBOL,logFC)

gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)

gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]

id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID

kk <- gseKEGG(id.fc, organism = "hsa",pvalueCutoff = 1)
dim(kk)

kk.gsym <- setReadable(kk, 'org.Hs.eg.db', 
                       'ENTREZID')
sortkk <- kk.gsym[kk.gsym@result$Description %like% "Platelet activation" | 
                    kk.gsym@result$Description %like% "Neuroactive ligand-receptor interaction" | 
                    kk.gsym@result$Description %like% "Lipid and atherosclerosis",]

go <- data.frame(Category = "KEGG",
                 ID = sortkk$ID,
                 Term = sortkk$Description, 
                 Genes = gsub("/", ", ", sortkk$core_enrichment), 
                 adj_pval = sortkk$p.adjust)

genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)

circ <- circle_dat(go, genelist)
head(circ)
#write.csv(circ[,c(3,5,6)],"SCM_ML/0.Preparing/very_easy_input.csv", quote = F, row.names = F)
df <- read.csv("SCM_ML/0.Preparing/very_easy_input.csv")
nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root="all")
edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
nodes <- nodes %>% mutate_at(c("node.level","node.branch"),as.character)

graph <- tbl_graph(nodes, edges)

gc1 <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 

  geom_edge_diagonal(aes(color = node1.node.branch,
                         filter=node1.node.level!="all"), 
                     alpha = 1/3,edge_width=1)+ 
  scale_edge_color_manual(values = c("#61C3ED","red","purple","darkgreen")) + 
  

  geom_node_point(aes(size = node.size,
                      filter=node.level!="all"), 
                  #alpha = 1/3,
                  color = "#61C3ED") + 
  #scale_size(range = c(0.5,80)) + 
  theme(legend.position = "none") + 
  
  geom_node_text(
    aes(
      x = 1.05 * x, 
      y = 1.05 * y, 
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = leaf
    ),
    color="black", 
    size = 6, hjust = 'outward') +
  
  geom_node_text(
    aes(label=node.short_name,
        filter = !leaf & (node.level != "all")
    ),
    color="black", 
    fontface="bold",
    size=6,
    family="sans"
  ) + 
  theme(panel.background = element_rect(fill = NA)) + 
  coord_cartesian(xlim=c(-1.3,1.3),ylim = c(-1.3,1.3)) 

ggsave("gc1.pdf",width = 14,height = 14)
####5.Figure3D####
#GSVA requires an additional full process

#Preparatio
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

#Preparing data for input
msigdbr_species()
h <- msigdbr(species = "Homo sapiens", 
             category = "H") 
#https://blog.csdn.net/weixin_43569478/article/details/83744521
h <- dplyr::select(h, gs_name, gene_symbol) %>% 
  as.data.frame %>% 
  split(., .$gs_name) %>% 
  lapply(., function(x)(x$gene_symbol)) 

gs <- lapply(h, unique)

count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 2))
gs <- lapply(gs, function(x) intersect(keep, x))
gs <- gs[lapply(gs, length) > 0]
#save(gs, file = "hallmark.gs.RData")

load("SCM_ML/0.Preparing/expression.Rdata")

#Analyze and visualize
gsym.expr <- SCM_control_exp
gsva_es <- gsva(as.matrix(gsym.expr), gs)
head(gsva_es)

group_list <- data.frame(sample = colnames(gsva_es), group = group_list)
head(group_list)

group_list$group <- c(rep("control",6),rep("SCM",15))

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design


contrast.matrix <- makeContrasts(SCM-control, levels = design)


fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)

cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))


sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) + 
  

  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2,
             size = 0.3) + 
  

  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),
            size = 3, 
            hjust = "outward" ) +  
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.03, label=ID, color = group),
            size = 3, hjust = "inward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score")+
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴

#ggsave("gsva.pdf", width = 8, height = 8)