setwd("~/wsx/projects/jlx/20220928/2-exosome/")
library(tidyverse)
library(data.table)
load("../1-dat/01geo31312.cel.RData")
load("../1-dat/01tcga.DLBC.RData")
# load("../1-dat/01geo10846.RData")
# load("../1-dat/01geo31312.RData")
# load("../1-dat/01geo83632.RData")
# load("../1-dat/01geo31312.RData")
# load("../1-dat/01geo9327.RData")

# 读取exosome
exosome <- openxlsx::read.xlsx("34496872.xlsx",colNames = F)
# exosome <- read.csv("Gene general (mRNA).csv")
exosome <- list(exosome=exosome$X1)

# ssGSEA
require(GSVA)
require(GSEABase)
#log2/FPKM/TPM用Gaussian算法，原始count数据用Possion算法，单样本使用ssgsea方法
Es.exosome = gsva(expr=as.matrix(expr31312),gset.idx.list = exosome,kcdf="Gaussian",parallel.sz=10,method="ssgsea")
Es.exosome <- as.data.frame(t(Es.exosome)) %>% rownames_to_column(var = "gsm")

# 生存分析
library(ggplot2)
library(survival)
library(survivalROC)
require(survminer)
surv.dat <- clin31312 %>% 
  dplyr::inner_join(Es.exosome,by="gsm")
# cut.point <- surv_cutpoint(surv.dat,time = "os.time",event = "os",variables = "exosome")
# cut.point <- as.vector(cut.point$cutpoint[1,1])
# surv.dat <- surv.dat %>% dplyr::mutate(group=ifelse(exosome>cut.point,"high","low"))
surv.dat <- surv.dat %>% dplyr::mutate(group=ifelse(exosome>median(exosome),"high","low"))
fit <- survfit(Surv(os.time,os)~group,data=surv.dat)
pdf("../figure/02gse31312.surv.pdf",onefile = F)
ggsurvplot(fit, data = surv.dat,
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
           legend.labs = c("high-score","low-score"),
           legend.title = "",
           ggtheme = theme_survminer()+theme(legend.background = element_rect(fill = "transparent",colour = NA))
) 
dev.off()

# 临床相关性
library(ggpubr)
surv.dat <- surv.dat %>% 
  dplyr::mutate(age.group=ifelse(age>=40,"age>=40","age<40"))
pdf("../figure/02half.violin.age.pdf",height=4,width=2,onefile = FALSE)
ggviolin(surv.dat %>% dplyr::filter(!is.na(age.group)),
         x="age.group",
         y="exosome",
         fill = "age.group",
         palette = c(),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = list(c("age>=40","age<40")), label = "p.signif",
                     label.y = 1.1,
                     method = "wilcox.test",paired = F)+
  theme_bw()+
  # scale_y_log10() +
  # ylim(0.4,1.4)+
  ylab("Enrichment Score of Exosome") + xlab("")+
  theme(axis.text.x = element_text(angle = 0,size = 8),
        legend.position = "none")
dev.off()

surv.dat <- surv.dat %>% dplyr::mutate(gender=ifelse(gender=="F","Female","Male"))
table(surv.dat$gender)
pdf("../figure/02half.violin.gender.pdf",height=4,width=2,onefile = FALSE)
ggviolin(surv.dat,
         x="gender",
         y="exosome",
         fill = "gender",
         palette = c(),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = list(c("Female","Male")), label = "p.signif",
                     label.y = 1.1,
                     method = "wilcox.test",paired = F)+
  theme_bw()+
  # scale_y_log10() +
  # ylim(0.4,1.4)+
  ylab("Enrichment Score of Exosome") + xlab("")+
  theme(axis.text.x = element_text(angle = 0,size = 8),
        legend.position = "none")
dev.off()

table(surv.dat$immu.subtype)
pdf("../figure/02half.violin.immu.subtype.pdf",height=4,width=2,onefile = FALSE)
ggviolin(surv.dat %>% dplyr::filter(!is.na(immu.subtype)),
         x="immu.subtype",
         y="exosome",
         fill = "immu.subtype",
         palette = c(),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = list(c("ABC","GCB")), label = "p.signif",
                     label.y = 1.1,
                     method = "wilcox.test",paired = F)+
  theme_bw()+
  # scale_y_log10() +
  # ylim(0.4,1.4)+
  ylab("Enrichment Score of Exosome") + xlab("")+
  theme(axis.text.x = element_text(angle = 0,size = 8),
        legend.position = "none")
dev.off()

table(surv.dat$stage)
surv.dat <- surv.dat %>% dplyr::mutate(stage=factor(stage,levels = c("stage I","stage II","stage III","stage IV")))
pdf("../figure/02half.violin.stage.pdf",height=4,width=4,onefile = FALSE)
ggviolin(surv.dat %>% dplyr::filter(!is.na(stage)),
         x="stage",
         y="exosome",
         fill = "stage",
         palette = c(),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = list(c("stage I","stage II"),c("stage I","stage III"),
                                        c("stage I","stage IV"),c("stage II","stage III"),
                                        c("stage II","stage IV"),c("stage III","stage IV")), 
                     label = "p.signif",
                     # label.y = 1.1,
                     method = "wilcox.test",paired = F)+
  theme_bw()+
  # scale_y_log10() +
  # ylim(0.4,1.4)+
  ylab("Enrichment Score of Exosome") + xlab("")+
  theme(axis.text.x = element_text(angle = 0,size = 8),
        legend.position = "none")
dev.off()

table(surv.dat$ecog)
surv.dat <- surv.dat %>% 
  dplyr::mutate(ecog=ifelse(ecog%in%c(0,1,2),"ECOG 0/1/2","ECOG 3/4"))
pdf("../figure/02half.violin.ecog.pdf",height=4,width=2,onefile = FALSE)
ggviolin(surv.dat %>% dplyr::filter(!is.na(ecog)),
         x="ecog",
         y="exosome",
         fill = "ecog",
         palette = c(),
         add = "boxplot",
         add.params = list(fill="white"))+
  stat_compare_means(comparisons = list(c("ECOG 0/1/2","ECOG 3/4")), label = "p.signif",
                     label.y = 1.1,
                     method = "wilcox.test",paired = F)+
  theme_bw()+
  # scale_y_log10() +
  # ylim(0.4,1.4)+
  ylab("Enrichment Score of Exosome") + xlab("")+
  theme(axis.text.x = element_text(angle = 0,size = 8),
        legend.position = "none")
dev.off()

# 差异分析
library(limma)
##分组信息
index <- match(colnames(expr31312),surv.dat$gsm)
surv.dat <- surv.dat[index,]
group_list=factor(surv.dat$group)
table(group_list)
## 限定对比
group_list <- relevel(group_list, ref="low")
#无配对的
design=model.matrix(~ group_list)
fit=lmFit(expr31312,design)
fit=eBayes(fit) 
allDiff.group=topTable(fit,adjust='fdr',number=Inf) 
allDiff.group <- na.omit(allDiff.group)
diff.group=allDiff.group %>% dplyr::filter(abs(logFC)>0.5,adj.P.Val<0.05)
write.csv(diff.group,file = "02diff.csv")
## 火山图
library(ggpubr)
colors <-  c("#3581b9", "grey","#CF5246")
volca.dat <- allDiff.group %>% 
  dplyr::mutate(yplot=-log10(adj.P.Val),
                state=ifelse(logFC>0.5 & adj.P.Val<0.05,"up",ifelse(logFC<(-0.5) & adj.P.Val<0.05,"down","unchange")),
                logFC=ifelse(logFC>10,10,ifelse(logFC<(-10),-10,logFC))
  )
pdf("../figure/02volcano.geo.pdf",height = 6,width = 7,onefile = FALSE)
ggscatter(data=volca.dat,x="logFC",y="yplot",color = "state",alpha = 0.3,
          size = 2,palette = colors,xlab = "logFC",ylab = "-log10(Pvalue)")+
  grids(linetype = "solid")+
  geom_vline(xintercept = c(-0.5,0.5),lty=4,col="black",lwd=1)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=1)+
  theme_bw()+
  theme(text = element_text(size=14),
        legend.position = c(0.2,0.85),
        legend.background = element_blank())
dev.off()

## 热图
library(pheatmap)
surv.dat <- surv.dat %>% dplyr::arrange(group)
pheat.dat <- expr31312 %>% 
  dplyr::select(surv.dat$gsm) %>% 
  dplyr::filter(rownames(.)%in%rownames(diff.group))
pheat.anno <- surv.dat %>% 
  column_to_rownames(var = "gsm") %>% 
  dplyr::select(group)
pdf("../figure/02heatmap.geo31312.pdf",height = 6,width = 8,onefile = FALSE)
pheatmap(as.matrix(pheat.dat),show_rownames = F, show_colnames = F,
         color = colorRampPalette(colors = c("#3581b9", "white","#CF5246"))(100),
         # color = colorRampPalette(colors = c(sc.cell.type[5],"lightyellow",sc.cell.type[1]))(100),
         annotation_colors=list(group=c(low="#3581b9", high="#CF5246")),
         cluster_rows = T,cluster_cols = F,scale = "row",
         annotation_col = pheat.anno, annotation_names_col = T,
         use_raster=T,fontsize = 20)
dev.off()

save(surv.dat,allDiff.group,diff.group,file = "02EXs.RData")

