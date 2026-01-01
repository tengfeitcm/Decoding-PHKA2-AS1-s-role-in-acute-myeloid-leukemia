setwd("~/wsx/projects/jlx/20220928/1-dat/")
library(tidyverse)
library(data.table)
library(GEOquery)
library(AnnoProbe)

if (!file.exists("./GSE87371")) {dir.create("./GSE87371")}
# gset = getGEO(filename = "./GSE87371/GSE87371_series_matrix.txt.gz",getGPL = F)
gset = getGEO("GSE87371",destdir="./GSE87371",getGPL = F)
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
# probe2gene=idmap(index,destdir = "./GSE87371")
gpl <- getGEO(index, destdir="./GSE87371/")
probe2gene <- Table(gpl) %>%
  dplyr::select(probe_id=ID,symbol=`Gene Symbol`) %>%
  dplyr::filter(!is.na(symbol),symbol!="",) %>% 
  dplyr::mutate(probe_id=as.character(probe_id))

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
x <- clin$os.time
sum(is.na(x))
table(x)
# max(x)
clin <- pdata %>% 
  dplyr::select(gsm=geo_accession,
                age=`age:ch1`,gender=`Sex:ch1`,
                sample.type=`diagnosis:ch1`,
                immu.subtype=`coo:ch1`,
                chemo=`treatment:ch1`,
                pfs=`cens_pfs:ch1`,
                pfs.time=`pfs_time:ch1`,
                os=`cens_os:ch1`,
                os.time=`os_time:ch1`) %>% 
  dplyr::mutate(age=as.numeric(age),
                pfs=as.numeric(pfs),
                pfs.time=as.numeric(pfs.time),
                os=as.numeric(os),
                os.time=as.numeric(os.time)) %>% 
  dplyr::filter(!is.na(os),
                os.time>=1)

clin87371 <- clin
expr87371 <- exprSet %>% 
  dplyr::select(clin87371$gsm)

save(clin87371,expr87371,file = "01geo87371.RData")
