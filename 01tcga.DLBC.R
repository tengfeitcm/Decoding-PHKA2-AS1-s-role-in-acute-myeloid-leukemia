setwd("~/wsx/projects/jlx/20220928/1-dat/")
library(tidyverse)
library(data.table)
library(TCGAbiolinks)
library(SummarizedExperiment)

projects <- getGDCprojects()
query.exp <- GDCquery(
  project = "TCGA-DLBC", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query.exp)
TCGA_DLBC_Exp <- GDCprepare(query = query.exp)
# 临床数据整理
tcga.DLBC.clin <-SummarizedExperiment::colData(TCGA_DLBC_Exp) %>% 
  as.data.frame(.) %>% 
  dplyr::select(barcode,sample.type=definition,age=age_at_diagnosis,
                gender=gender,race=race,
                pathology=primary_diagnosis,
                location=site_of_resection_or_biopsy,
                stage=ann_arbor_clinical_stage,
                os=vital_status,
                day.to.last.follow.up=days_to_last_follow_up,
                day.to.death=days_to_death) %>% 
  dplyr::filter(!is.na(age),
                # !is.na(stage),
                ) %>% 
  dplyr::mutate(age=floor(age/365),
                barcode=str_sub(barcode,1,15),
                race=ifelse(race=="not reported","unkown",race),
                stage=str_replace_all(stage,"A|B",""),
                os=ifelse(os=="Alive",0,1),
                os.time=day.to.last.follow.up,
                os.time=ifelse(is.na(os.time),day.to.death,os.time),
                day.to.death=ifelse(is.na(day.to.death),day.to.last.follow.up,day.to.death),
                os.time=ifelse(os.time<day.to.death,day.to.death,os.time)) %>% 
  # dplyr::filter(os.time>30) %>% 
  dplyr::select(-day.to.death,-day.to.last.follow.up) %>% 
  dplyr::distinct(barcode,.keep_all = T)
table(tcga.DLBC.clin$sample.type)
# 表达矩阵整理
anno <- as.data.frame(TCGA_DLBC_Exp@rowRanges@elementMetadata@listData)
tcga.DLBC.exp <- as.data.frame(SummarizedExperiment::assay(TCGA_DLBC_Exp,1)) %>% 
  rownames_to_column(var = "gene_id") %>% 
  mutate(rowMean =rowMeans(.[,grep("TCGA", names(.))])) %>% 
  left_join(anno[,c("gene_id","gene_name")],by="gene_id") %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-c(rowMean,gene_id)) %>% 
  dplyr::filter(!is.na(gene_name)) %>%
  column_to_rownames(var = "gene_name") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "id") %>% 
  dplyr::mutate(id=str_sub(id,1,15)) %>% 
  group_by(id) %>% 
  summarise_all(mean) %>% 
  column_to_rownames(var = "id") %>% 
  t() %>% as.data.frame()
tcga.DLBC.tpm <- as.data.frame(SummarizedExperiment::assay(TCGA_DLBC_Exp,4)) %>% 
  rownames_to_column(var = "gene_id") %>% 
  mutate(rowMean =rowMeans(.[,grep("TCGA", names(.))])) %>% 
  left_join(anno[,c("gene_id","gene_name")],by="gene_id") %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-c(rowMean,gene_id)) %>% 
  dplyr::filter(!is.na(gene_name)) %>%
  column_to_rownames(var = "gene_name")%>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "id") %>% 
  dplyr::mutate(id=str_sub(id,1,15)) %>% 
  group_by(id) %>% 
  summarise_all(mean) %>% 
  column_to_rownames(var = "id") %>% 
  t() %>% as.data.frame()
tcga.DLBC.fpkm <- as.data.frame(SummarizedExperiment::assay(TCGA_DLBC_Exp,5)) %>% 
  rownames_to_column(var = "gene_id") %>% 
  mutate(rowMean =rowMeans(.[,grep("TCGA", names(.))])) %>% 
  left_join(anno[,c("gene_id","gene_name")],by="gene_id") %>% 
  arrange(desc(rowMean)) %>% 
  distinct(gene_name,.keep_all = T) %>% 
  dplyr::select(-c(rowMean,gene_id)) %>% 
  dplyr::filter(!is.na(gene_name)) %>%
  column_to_rownames(var = "gene_name")%>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "id") %>% 
  dplyr::mutate(id=str_sub(id,1,15)) %>% 
  group_by(id) %>% 
  summarise_all(mean) %>% 
  column_to_rownames(var = "id") %>% 
  t() %>% as.data.frame()
# 取交集
index <- intersect(colnames(tcga.DLBC.exp),tcga.DLBC.clin$barcode)
tcga.DLBC.clin <- tcga.DLBC.clin %>% dplyr::filter(barcode%in%index)
tcga.DLBC.exp <- tcga.DLBC.exp %>% dplyr::select(tcga.DLBC.clin$barcode)
tcga.DLBC.tpm <- tcga.DLBC.tpm %>% dplyr::select(tcga.DLBC.clin$barcode)
tcga.DLBC.tpm <- log2(tcga.DLBC.tpm+1)
tcga.DLBC.fpkm <- tcga.DLBC.fpkm %>% dplyr::select(tcga.DLBC.clin$barcode)
tcga.DLBC.fpkm <- log2(tcga.DLBC.fpkm+1)
save(tcga.DLBC.clin,tcga.DLBC.exp,tcga.DLBC.tpm,tcga.DLBC.fpkm,file = "01tcga.DLBC.RData")

# SNV数据下载
query.maf <- GDCquery(
  project = "TCGA-DLBC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  # legacy = FALSE,
  # data.type = "Masked Somatic Mutation",
)
GDCdownload(query.maf)
maf.DLBC <- GDCprepare(query = query.maf)
save(maf.DLBC,file = "01tcga.DLBC.snv.RData")

# CNV数据下载
query.CNV <- GDCquery(
  project = "TCGA-DLBC",
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",
  sample.type = c("Primary Tumor")
)
GDCdownload(query.CNV)
cnv.DLBC <- GDCprepare(query = query.CNV)
save(cnv.DLBC,file = "01tcga.DLBC.cnv.RData")
# 
# # 甲基化数据
# query.met.hg38 <- GDCquery(
#   project = "TCGA-DLBC",
#   data.category = "DNA Methylation",
#   platform = "Illumina Human Methylation 450",
#   data.type = "Methylation Beta Value"
# )
# GDCdownload(query.met.hg38)
# met.DLBC <- GDCprepare(query = query.met.hg38)
# save(met.DLBC,file = "01tcga.DLBC.met.RData")

load("01tcga.DLBC.RData")
table(tcga.DLBC.clin$sample.type)
