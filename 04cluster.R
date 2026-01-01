setwd("~/wsx/projects/jlx/20220928/4-cluster/")
library(tidyverse)
library(data.table)
library(ConsensusClusterPlus)
load("../1-dat/01geo31312.cel.RData")
load("../5-model/05model.RData")

# 一致性聚类
cc.data <- as.data.frame(expr31312) %>%
  dplyr::filter(rownames(.)%in%unique(stepCox.genes))
cc.data <- sweep(cc.data,1, apply(cc.data,1,median,na.rm=T))
# res <- ConsensusClusterPlus(as.matrix(cc.data),maxK = 10,reps = 1000,pItem = 0.8,
#                             pFeature = 1,clusterAlg = "pam",
#                             distance = "pearson",innerLinkage="complete",seed = 123456,
#                             plot = "png",title = "./consensusOut/")
# anno <- as.data.frame(res[[2]]$consensusClass) %>%
#   dplyr::rename(cluster="res[[2]]$consensusClass") %>%
#   rownames_to_column(var = "gsm") %>%
#   inner_join(clin31312,by="gsm") %>%
#   dplyr::arrange(cluster) %>%
#   dplyr::mutate(cluster=str_c("MC",cluster))
# save(res,anno,file = "04cluster.res.RData")
# write.table(anno,file = "04cluster.xlsx",sep = "\t",row.names = T,quote = F,col.names = T)

load("04cluster.res.RData")

# PCA分析
library(ggplot2)
sr <- rowSums(expr31312)
pca.dat <- expr31312 %>% dplyr::filter(sr>0)
dat.pca <- prcomp(t(pca.dat),center = T,scale. = T)
plot.pca <- predict(dat.pca)[,1:2]
plot.pca <- data.frame(plot.pca) %>%
  rownames_to_column(var = "gsm") %>%
  inner_join(anno[,c("gsm","cluster","immu.subtype")],by="gsm")
pdf("../figure/04pca.gse31312.pdf",height = 6,width = 6,onefile = FALSE)
ggplot(plot.pca,aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size=3) +
  stat_ellipse()+
  theme_bw()+
  # scale_x_continuous(limits=c(-150,150))+
  theme(text = element_text(size=20),
        legend.position = c(0.8,0.85),
        legend.background = element_blank())
dev.off()

# 生存分析
library(survival)
library(survivalROC)
require(survminer)
fit <- survfit(Surv(os.time,os)~cluster,data=anno)
pdf("../figure/04cluster.surv.pdf",onefile = F)
ggsurvplot(fit, data = anno,
           conf.int = TRUE,
           pval = TRUE,
           pval.coord = c(1,1),
           fun = "pct",
           risk.table = TRUE,
           risk.table.height = 0.2,
           size = 1,
           linetype = "strata",
           ##调色板
           legend = c(0.13,0.16),
           # legend.labs = c("high-score","low-score"),
           legend.title = "",
           ggtheme = theme_survminer()+theme(legend.background = element_rect(fill = "transparent",colour = NA))
) 
dev.off()

# 热图
library(pheatmap)
anno <- anno %>% dplyr::arrange(cluster)
pheat.dat <- expr31312 %>% 
  dplyr::select(anno$gsm) %>% 
  dplyr::filter(rownames(.)%in%stepCox.genes)
pheat.anno <- anno %>% 
  column_to_rownames(var = "gsm") %>% 
  dplyr::select(cluster)
pdf("../figure/04heatmap.geo31312.pdf",height = 6,width = 12,onefile = FALSE)
pheatmap(as.matrix(pheat.dat),show_rownames = T, show_colnames = F,
         # color = colorRampPalette(colors = c("#3581b9", "white","#CF5246"))(100),
         # color = colorRampPalette(colors = c(sc.cell.type[5],"lightyellow",sc.cell.type[1]))(100),
         # annotation_colors=list(group=c(MC1="#CF5246", MC2="#3581b9")),
         cluster_rows = T,cluster_cols = F,scale = "row",
         annotation_col = pheat.anno, annotation_names_col = T,
         use_raster=T,fontsize = 20)
dev.off()
