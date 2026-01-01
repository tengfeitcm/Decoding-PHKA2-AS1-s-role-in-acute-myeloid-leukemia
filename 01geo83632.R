setwd("~/wsx/projects/jlx/20220928/1-dat/")
library(tidyverse)
library(data.table)
library(GEOquery)
library(AnnoProbe)

if (!file.exists("./GSE83632")) {dir.create("./GSE83632")}
# gset = getGEO(filename = "./GSE83632/GSE83632_series_matrix.txt.gz",getGPL = F)
gset = getGEO("GSE83632",destdir="./GSE83632",getGPL = F)
## 获取分组信息,点开查阅信息
# pdata=pData(gset)
pdata=pData(gset[[1]])
##获取矩阵并校正
# exprSet=exprs(gset)
exprSet=exprs(gset[[1]])
exprSet[1:4,1:4]
# boxplot(exprSet,outline=FALSE, notch=T,col=group, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
# boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
##判断取log2
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

##制作转换文件（有些平台是没有对应的R包的，可以自己下载平台信息来做）
index = gset[[1]]@annotation
# probe2gene=idmap(index,destdir = "./GSE83632")
gpl <- getGEO(index, destdir="./GSE83632/")
probe2gene <- Table(gpl) %>% 
  dplyr::select(probe_id=ID,symbol=gene_assignment) %>% 
  dplyr::mutate(symbol=str_split(symbol," // ",simplify = T)[,2],
                probe_id=as.character(probe_id)) %>% 
  dplyr::filter(symbol!="---",!is.na(symbol),symbol!="",)

exprSet <- as.data.frame(exprSet) %>% 
  rownames_to_column(var="probe_id") %>% 
  #合并探针的信息
  inner_join(probe2gene,by="probe_id") %>% 
  #去掉多余信息
  dplyr::select(-probe_id) %>% 
  #重新排列
  dplyr::select(symbol,everything()) %>% 
  #去除symbol中的NA
  dplyr::filter(symbol != "NA") %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise_all(mean) %>% 
  # 列名变成行名
  column_to_rownames(var = "symbol")
##临床数据处理
# x <- clin$sample.type
# sum(is.na(x))
# table(x)
# max(x)
clin <- pdata %>% 
  dplyr::select(gsm=geo_accession,age=`age:ch1`,gender=`gender:ch1`,
                sample.type=source_name_ch1) %>% 
  # inner_join(paper.clin,by="title") %>% 
  # dplyr::filter(os.time>=1) %>% 
  dplyr::mutate(location="blood")

clin83632 <- clin
expr83632 <- exprSet %>% 
  dplyr::select(clin83632$gsm)

save(clin83632,expr83632,file = "01geo83632.RData")
