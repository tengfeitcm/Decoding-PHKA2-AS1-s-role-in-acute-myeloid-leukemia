setwd("~/wsx/projects/jlx/20220928/1-dat/")
library(tidyverse)
library(data.table)
library(GEOquery)
library(AnnoProbe)

if (!file.exists("./GSE31312")) {dir.create("./GSE31312")}
# gset = getGEO(filename = "./GSE31312/GSE31312_series_matrix.txt.gz",getGPL = F)
gset = getGEO("GSE31312",destdir="./GSE31312",getGPL = F)
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
probe2gene=idmap(index,destdir = "./GSE31312")

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
paper.clin1 <- openxlsx::read.xlsx("./GSE31312/GSE31312.xlsx",sheet = "Table 1")
paper.clin2 <- openxlsx::read.xlsx("./GSE31312/GSE31312.xlsx",sheet = "Table 2")
paper.clin <- paper.clin2 %>% 
  inner_join(paper.clin1,by="Case") %>% 
  dplyr::select(title=GEO.Depository,age=Age,
                gender=Gender,stage=AAS,immu.subtype=GEP,
                os=OScensor,os.time=OS,
                ecog=ECOG) %>% 
  dplyr::mutate(stage=ifelse(stage=="1","stage I",
                             ifelse(stage=="2","stage II",
                                    ifelse(stage=="3","stage III",
                                           ifelse(stage=="4","stage IV", NA)))))
# x <- clin$ecog
# sum(is.na(x))
# table(x)
clin <- pdata %>% 
  dplyr::select(title,gsm=geo_accession) %>% 
  inner_join(paper.clin,by="title") %>% 
  dplyr::mutate(immu.subtype=ifelse(immu.subtype=="UC",NA,immu.subtype)) %>% 
  dplyr::filter(os.time>=1)

clin31312 <- clin
expr31312 <- exprSet %>% 
  dplyr::select(clin31312$gsm)

save(clin31312,expr31312,file = "01geo31312.RData")
