rm(list = ls())
setwd("~/wsx/projects/jlx/20220928/6-diff/")
library(tidyverse)
library(data.table)
library(limma)
load("../1-dat/01geo31312.cel.RData")
load("../5-model/05model.RData")

##分组信息
index <- match(colnames(expr31312),stepCox.dat$gsm)
stepCox.dat <- stepCox.dat[index,]
group_list=factor(stepCox.dat$group)
## 限定对比
group_list <- relevel(group_list, ref="low")
#无配对的
design=model.matrix(~ group_list)
fit=lmFit(expr31312,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',number=Inf) 
allDiff <- na.omit(allDiff)
diff=allDiff %>% dplyr::filter(abs(logFC)>0.5,adj.P.Val<0.05)
write.csv(diff,file = "06diff.csv")

## 火山图
library(ggpubr)
colors <-  c("#3581b9", "grey","#CF5246")
volca.dat <- allDiff %>% 
  dplyr::mutate(yplot=-log10(adj.P.Val),
                state=ifelse(logFC>0.5 & adj.P.Val<0.05,"up",ifelse(logFC<(-0.5) & adj.P.Val<0.05,"down","unchange")),
                logFC=ifelse(logFC>1,1,ifelse(logFC<(-1),-1,logFC)))
pdf("../figure/06volcano.geo.pdf",height = 6,width = 7,onefile = FALSE)
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
stepCox.dat <- stepCox.dat %>% dplyr::arrange(group)
pheat.dat <- as.data.frame(expr31312) %>% 
  dplyr::filter(rownames(.)%in%rownames(diff)) %>% 
  dplyr::select(stepCox.dat$gsm)
anno <- data.frame(group=stepCox.dat$group,row.names = stepCox.dat$gsm)
pdf("../figure/06heatmap.geo.pdf",height = 6,width = 8,onefile = FALSE)
pheatmap(as.matrix(pheat.dat),show_rownames = F, show_colnames = F,
         # color = colorRampPalette(colors = c("#3581b9", "white","#CF5246"))(100),
         # color = colorRampPalette(colors = c(sc.cell.type[5],"lightyellow",sc.cell.type[1]))(100),
         annotation_colors=list(group=c(high="#CF5246", low="#3581b9")),
         cluster_rows = T,cluster_cols = F,scale = "row",
         annotation_col = anno, annotation_names_col = T,
         use_raster=T,fontsize = 20)
dev.off()
dev.off()

# GO和KEGG
# 颜色
color2 <- c("#6F80be","#f47e62")
color3 <- c("#6F80be","#f47e62","#BD8EC0")
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(GOplot)#弦图，弦表图，系统聚类图
genes <- rownames(diff)
diff$gene <- rownames(diff)
eg <- bitr(genes, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db") %>% 
  dplyr::rename(gene=SYMBOL) %>% 
  inner_join(diff[,c("gene","logFC")],by="gene")
ego <- enrichGO(gene          = unique(eg$ENTREZID),
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                qvalueCutoff  = 0.05)
pdf("../figure/06diff.go.pdf",height = 8,width = 6,onefile = FALSE)
enrichplot::dotplot(ego,showCategory = 6,split="ONTOLOGY")+
  scale_color_continuous(low = color2[1],high = color2[2])+
  facet_grid(ONTOLOGY~.,scale="free")
dev.off()

kk <- enrichKEGG(gene         = eg$ENTREZID,
                 organism     = 'hsa',
                 qvalueCutoff = 0.05)
goPlotKegg <- kk@result[,c(1,2,6,8)] %>% 
  dplyr::mutate(geneID=str_replace_all(geneID,"/",",")) %>% 
  dplyr::mutate(Category="KEGG")
names(goPlotKegg)<-c('ID','Term','adj_pval','Genes','Category')
genedata <- diff %>% 
  inner_join(eg[,1:2],by="gene")%>% 
  dplyr::select(ID=ENTREZID,logFC=logFC)
circ_kegg<-GOplot::circle_dat(goPlotKegg,genedata) 
pdf("../figure/06diff.kegg.pdf",height = 8,width = 6,onefile = FALSE)
GOCircle(circ_kegg,table.legend=F,
         zsc.col = c(color2[1],"white",color2[2]),
         lfc.col = c(color2[1],color2[2]))+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom" )
dev.off()
ego.d <- ego@result
kk.d <- kk@result
d <- ego.d %>% 
  bind_rows(kk.d) 
write.table(d,file = "06go.kegg.xlsx",sep = "\t",row.names = T,quote = F,col.names = T)


# GSEA分析
library(msigdbr)
library(future.apply)
geneList <- eg$logFC
names(geneList) <- eg$gene
geneList <- sort(geneList,decreasing = T)
msigdb.c2 <- read.gmt("c2.cp.v7.5.1.symbols.gmt")
set.seed(123456)
egsea <- GSEA(geneList,TERM2GENE=msigdb.c2, pvalueCutoff = 1)
# egsea1 <- egsea %>%
#   dplyr::arrange(qvalues)
# egsea1@result <- egsea1@result[c(1:8),]
egsea1 <- egsea
rownames(egsea1@result)
egsea1@result <- egsea@result[c(4,5,9,21,34,39,41),]
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    yulab.utils::str_wrap(str, n)
  }
}
ridgeplot.gseaResult <- function(x, showCategory=30, fill="p.adjust",
                                 core_enrichment = TRUE, label_format = 30,
                                 orderBy = "NES", decreasing = FALSE) {
  if (!is(x, "gseaResult"))
    stop("currently only support gseaResult")

  ## fill <- match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
  if (fill == "qvalue") {
    fill <- "qvalues"
  }
  if (!fill %in% colnames(x@result)) {
    stop("'fill' variable not available ...")
  }

  ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')
  if (orderBy !=  'NES' && !orderBy %in% colnames(x@result)) {
    message('wrong orderBy parameter; set to default `orderBy = "NES"`')
    orderBy <- "NES"
  }
  n <- showCategory
  if (core_enrichment) {
    gs2id <- geneInCategory(x)[seq_len(n)]
  } else {
    gs2id <- x@geneSets[x$ID[seq_len(n)]]
  }

  if (x@readable) {
    id <- match(names(x@geneList), names(x@gene2Symbol))
    names(x@geneList) <- x@gene2Symbol[id]
  }

  gs2val <- lapply(gs2id, function(id) {
    res <- x@geneList[id]
    res <- res[!is.na(res)]
  })

  nn <- names(gs2val)
  i <- match(nn, x$ID)
  nn <- x$Description[i]

  # j <- order(x$NES[i], decreasing=FALSE)
  j <- order(x@result[[orderBy]][i], decreasing = decreasing)
  len <- sapply(gs2val, length)
  gs2val.df <- data.frame(category = rep(nn, times=len),
                          color = rep(x[i, fill], times=len),
                          value = unlist(gs2val))

  colnames(gs2val.df)[2] <- fill
  gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])

  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }

  ggplot(gs2val.df, aes_string(x="value", y="category", fill=fill)) +
    ggridges::geom_density_ridges() +
    ## scale_x_reverse() +
    scale_fill_continuous(low = "#CF5246",high = "#3581b9", name = fill,
                          guide=guide_colorbar(reverse=TRUE)) +
    scale_y_discrete(labels = label_func) +
    ## scale_fill_gradientn(name = fill, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ## geom_vline(xintercept=0, color='firebrick', linetype='dashed') +
    xlab(NULL) + ylab(NULL) +  theme_dose()+
    theme(text = element_text(size=20))
}
pdf("../figure/06gsea.ridge.pdf",height = 4,width = 10,onefile = FALSE)
ridgeplot.gseaResult(egsea1)
dev.off()

gseaScores <- getFromNamespace("gseaScores", "DOSE")
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
gseaplot2.self <- function(x, geneSetID, title = "", color="green", base_size = 11,
                      rel_heights=c(1.5, .5, 1), subplots = 1:3,
                      pvalue_table = FALSE, ES_geom="line",v=0,h=0,num=100) {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
  geneList <- position <- NULL ## to satisfy codetool
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
    anno <- x[geneSetID, c("NES", "pvalue", "p.adjust")]
    lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                          size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                           size=1, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer +
    theme(legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))+
    annotate("text", num, x[geneSetID, "enrichmentScore"] * .75,label = lab,
             hjust=h, vjust=v)
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  if (length(geneSetID) == 1) {
    ## geneList <- gsdata$geneList
    ## j <- which.min(abs(geneList))
    ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
    ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]
    
    ## v <- sort(c(v1, v2))
    ## inv <- findInterval(geneList, v)
    
    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      alpha=.9,
      inherit.aes=FALSE)
  }
  
  ## p2 <- p2 +
  ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
  ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
  ## theme(legend.position="none") +
  ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                            color="grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +
    xlab("Rank in Ordered Dataset") +
    theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    rownames(pd) <- pd$Description
    
    pd <- pd[,-1]
    # pd <- round(pd, 4)
    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 4)
    }
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())
  
  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                    l=.2, unit="cm")))
  
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  # aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
  aplot::gglist(gglist = plotlist, ncol=1, heights=rel_heights) 
}

pdf(str_c("../figure/06gsea.",1,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 1, title = egsea$Description[1],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",2,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 2, title = egsea$Description[2],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",3,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 3, title = egsea$Description[3],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",4,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 4, title = egsea$Description[4],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",5,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 5, title = egsea$Description[5],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",6,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 6, title = egsea$Description[6],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",7,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 7, title = egsea$Description[7],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",8,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 8, title = egsea$Description[8],v=0.5,h=1,num=100)
dev.off()
pdf(str_c("../figure/06gsea.",9,".pdf"),height = 4,width = 8,onefile = FALSE)
gseaplot2.self(egsea, geneSetID = 9, title = egsea$Description[9],v=0.5,h=1,num=100)
dev.off()

# PPI
ppi.net <- read.csv("miRWalk_miRNA_Targets.csv") %>% 
  dplyr::filter(validated!="") %>% 
  dplyr::mutate(select=str_c(mirnaid,genesymbol)) %>% 
  dplyr::distinct(select,.keep_all = T)
write.csv(ppi.net,file = "miRWalk.valid.csv")

# 相关性分析
library(corrplot)
ppi.hub <- read.csv("miRWalk_miRNA_Targets.csv")
ppi.hub <- unique(ppi.hub$genesymbol)
load("../2-exosome/02EXs.RData")
cor.dat <- expr31312 %>% 
  dplyr::filter(rownames(.)%in%ppi.hub) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "gsm") %>% 
  inner_join(surv.dat[,c("gsm","exosome")],by="gsm") %>% 
  dplyr::select(-gsm)
cor.res <- cor(as.matrix(cor.dat), method = c("spearman"))
p.res <- cor.mtest(as.matrix(cor.dat))
pdf("../figure/06exosome.corrplot.pdf",onefile = FALSE)
corrplot.mixed(cor.res,upper = "ellipse",lower = "square",tl.pos = "lt",
               number.cex=0.3,tl.col="black",
               p.mat=p.res$p,sig.level = 0.05, insig = "label_sig",pch.cex=0.8)
dev.off()


save(ppi.hub,file = "06ppi.hub.RData")
