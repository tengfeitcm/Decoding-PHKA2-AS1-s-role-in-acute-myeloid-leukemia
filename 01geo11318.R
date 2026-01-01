setwd("~/wsx/projects/jlx/20220928/1-dat/")
library(tidyverse)
library(data.table)
library(GEOquery)
library(AnnoProbe)

if (!file.exists("./GSE11318")) {dir.create("./GSE11318")}
# gset = getGEO(filename = "./GSE11318/GSE11318_series_matrix.txt.gz",getGPL = F)
gset = getGEO("GSE11318",destdir="./GSE11318",getGPL = F)
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
probe2gene=idmap(index,destdir = "./GSE11318")
# gpl <- getGEO(index, destdir="./GSE11318/")
# probe2gene <- Table(gpl) %>% 
#   dplyr::select(probe_id=ID,symbol=gene_assignment) %>% 
#   dplyr::mutate(symbol=str_split(symbol," // ",simplify = T)[,2],
#                 probe_id=as.character(probe_id)) %>% 
#   dplyr::filter(symbol!="---",!is.na(symbol),symbol!="",)

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
  dplyr::select(gsm=geo_accession,age=`Age:ch1`,gender=`Gender:ch1`,
                sample.type=`Disease state:ch1`,location=`Tissue:ch1`,
                stage=characteristics_ch1.11,
                immu.subtype=characteristics_ch1.6,
                os=characteristics_ch1.7,
                os.time=characteristics_ch1.8,
                chemo=characteristics_ch1.9,
                ecog=characteristics_ch1.10,
                ldh=characteristics_ch1.12) %>% 
  dplyr::mutate(stage=str_replace_all(stage,"Clinical info: Stage: ",""),
                immu.subtype=str_replace_all(immu.subtype,"Clinical info: Final microarray diagnosis: ",""),
                immu.subtype=str_replace_all(immu.subtype," DLBCL",""),
                os=ifelse(os=="Clinical info: Follow up status: DEAD",1,0),
                os.time=str_replace_all(os.time,"Clinical info: Follow up years: ",""),
                chemo=str_replace_all(chemo,"Clinical info: Chemotherapy: ",""),
                ecog=str_replace_all(ecog,"Clinical info: ECOG performance status: ",""),
                ldh=str_replace_all(ldh,"Clinical info: LDH ratio: ",""),
                age=as.numeric(age),
                stage=as.numeric(stage),
                os.time=as.numeric(os.time),
                ecog=as.numeric(ecog),
                ldh=as.numeric(ldh),
                os.time=os.time*12) %>% 
  dplyr::filter(!is.na(os),!is.na(os.time)) %>% 
  # inner_join(paper.clin,by="title") %>% 
  dplyr::filter(os.time>=1)

clin11318 <- clin
expr11318 <- exprSet %>% 
  dplyr::select(clin11318$gsm)

save(clin11318,expr11318,file = "01geo11318.RData")
