setwd("~/wsx/projects/jlx/20220928/1-dat/")
library(tidyverse)
library(data.table)
library(GEOquery)
library(AnnoProbe)

if (!file.exists("./GSE56315")) {dir.create("./GSE56315")}
# gset = getGEO(filename = "./GSE56315/GSE56315_series_matrix.txt.gz",getGPL = F)
gset = getGEO("GSE56315",destdir="./GSE56315",getGPL = F)
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
probe2gene=idmap(index,destdir = "./GSE56315")

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
x <- clin$immu.subtype
sum(is.na(x))
table(x)
# max(x)
clin <- pdata %>% 
  dplyr::select(gsm=geo_accession,
                sample.type=characteristics_ch1,
                immu.subtype=characteristics_ch1.5) %>% 
  dplyr::mutate(sample.type=ifelse(sample.type=="tissue: human healthy tonsils","health","DLBCL"),
                immu.subtype=str_replace_all(immu.subtype,"abc/gcb/nc subclass: ",""),
                immu.subtype=ifelse(immu.subtype%in%c("","Unclassified"),NA,immu.subtype)) 

clin56315 <- clin
expr56315 <- exprSet %>% 
  dplyr::select(clin56315$gsm)

save(clin56315,expr56315,file = "01geo56315.RData")
