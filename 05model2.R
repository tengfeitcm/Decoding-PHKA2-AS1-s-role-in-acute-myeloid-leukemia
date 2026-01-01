rm(list = ls())
setwd("~/wsx/projects/jlx/20220928/5-model/")
library(tidyverse)
library(data.table)
library(glmnet)
library(future.apply)
library(rlist)
library(survival)
library(survivalROC)
require(survminer)
load("../1-dat/01geo31312.cel.RData")
load("../1-dat/01tcga.DLBC.RData")
load("../3-wgcna/03mut.genes.RData")


# 数据处理
tcga.DLBC.clin <- tcga.DLBC.clin %>% dplyr::mutate(os.time=os.time/30)
rownames(tcga.DLBC.tpm) <- str_replace_all(rownames(tcga.DLBC.tpm),"-","_")
train.dat <- expr31312 %>% 
  dplyr::filter(rownames(.)%in%mut.genes) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "gsm") %>% 
  inner_join(clin31312[,c("gsm","os","os.time")],by="gsm")

# 设置参数
seed=123456

# 设置函数
surv <- function(yT,xVar,time_inc,year){
  cutoff <- time_inc*year
  sur.fit= survivalROC(Stime=yT$os.time,##生存时间
                       status=yT$os,## 终止事件    
                       marker = yT[,xVar], ## marker value    
                       predict.time = cutoff,## 预测时间截点
                       method="KM")##span,NNE法的namda
  return(sur.fit$AUC)
}
unicox<-function(yT,xVar){
  FML<-as.formula(paste0("BaSurv~",xVar))
  Gcox<-coxph(FML,data=yT)
  Gsum<-summary(Gcox)
  HR<-round(Gsum$coefficients[,2],2)
  Pvalue<-round(Gsum$coefficients[,5],3)
  CI<-data.frame(lower=round(Gsum$conf.int[,3],2),upper=round(Gsum$conf.int[,4],2)) %>% 
    dplyr::mutate(ci=str_c(lower,"-",upper)) %>% 
    .[,"ci"]
  lowCI<-round(Gsum$conf.int[,3],2)
  highCI<-round(Gsum$conf.int[,4],2)
  coef<-round(Gsum$coefficients[,1],2)
  lowCoef<-round(log(Gsum$conf.int[,3]),2)
  highCoef<-round(log(Gsum$conf.int[,4]),2)
  unicox<-data.frame("Characteristics"=xVar,
                     "group"=rownames(Gsum$coefficients),
                     "coef"=coef,
                     "lowCoef"=lowCoef,
                     "highCoef"=highCoef,
                     "Hazard Ratio"=HR,
                     "lowCI"=lowCI,
                     "highCI"=highCI,
                     "CI95"=CI,
                     "P Value"=Pvalue)
  return(unicox)
}

# 单因素分析
colnames(train.dat) <- str_replace_all(colnames(train.dat),"-","_")
##准备BaSurv
BaSurv<-Surv(time = train.dat$os.time,event = train.dat$os)
##提取变量名
names(train.dat)
varnames<-colnames(train.dat)[2:(length(mut.genes)+1)]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(train.dat,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
length(univar$Characteristics[univar$P.Value<0.05])
cox.genes <- univar$Characteristics[univar$P.Value<0.05]
cox.genes
write.csv(univar,file = "05univar.result.csv")

##lasso回归
lasso.dat <- train.dat %>% dplyr::select(gsm,os,os.time,all_of(cox.genes))
xT <- as.matrix(lasso.dat[,c(4:ncol(lasso.dat))])
yT <- Surv(lasso.dat$os.time,lasso.dat$os)
set.seed(123456)
cv.fit<-cv.glmnet(xT,yT,type.measure = "deviance",family = "cox")
pdf("../figure/05lasso1.pdf",height = 6,width = 8,onefile = FALSE)
plot(cv.fit,cex.lab=2,cex.axis=2)
dev.off()
fit<-glmnet(xT,yT,alpha=1,family='cox')
pdf("../figure/05lasso2.pdf",height = 6,width = 8,onefile = FALSE)
plot(fit, xvar = "lambda", label = TRUE,cex.lab=2,cex.axis=2)
abline(v=log(cv.fit$lambda.1se),lty=2)
abline(v=log(cv.fit$lambda.min),lty=2)
dev.off()
lasso.coef <- as.matrix(coef(cv.fit,s=cv.fit$lambda.min))
lasso.coef <- as.data.frame(lasso.coef) %>% 
  rownames_to_column(var = "gene") %>% 
  dplyr::rename(coef="1") %>% 
  dplyr::filter(coef!=0) %>% 
  .[-1,]
write.csv(lasso.coef,file = "05lasso.coef.csv")

# stepCox
stepCox.dat <- lasso.dat %>% dplyr::select(os,os.time,all_of(lasso.coef$gene))
fit.backward <- step(coxph(Surv(os.time,os)~.,stepCox.dat),direction = "backward")
stepCox.dat$pred=predict(fit.backward,type = 'risk',newdata = stepCox.dat)
stepCox.dat <- stepCox.dat %>% 
  dplyr::mutate(pred=(pred-min(pred))/(max(pred)-min(pred)))
roc.backward <- future_lapply(1:5,function(x){surv(stepCox.dat,"pred",12,x)},future.seed = NULL)
stepCox.genes <- rownames(summary(fit.backward)$coefficients)
fit.backward

# 绘图
# 生存曲线图
library(ggplot2)
# cut.point <- surv_cutpoint(stepCox.dat,time = "os.time",event = "os",variables = "pred")
# cut.point <- as.vector(cut.point$cutpoint[1,1])
# stepCox.dat <- stepCox.dat %>% dplyr::mutate(group=ifelse(pred>cut.point,"high","low"))
stepCox.dat <- stepCox.dat %>% dplyr::mutate(group=ifelse(pred>median(pred),"high","low"))
stepCox.dat <- cbind(gsm=lasso.dat[,1],stepCox.dat)
cut.point <- median(stepCox.dat$pred)
fit <- survfit(Surv(os.time,os)~group,data=stepCox.dat)
pdf("../figure/05train.surv.pdf",onefile = F)
ggsurvplot(fit, data = stepCox.dat,
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
# roc曲线图
sur.fit1 <- survivalROC(Stime=stepCox.dat$os.time,##生存时间
                        status=stepCox.dat$os,## 终止事件    
                        marker = stepCox.dat[,"pred"], ## marker value    
                        predict.time = 12*3,## 预测时间截点
                        method="KM")##span,NNE法的namda
sur.fit2 <- survivalROC(Stime=stepCox.dat$os.time,##生存时间
                        status=stepCox.dat$os,## 终止事件    
                        marker = stepCox.dat[,"pred"], ## marker value    
                        predict.time = 12*5,## 预测时间截点
                        method="KM")##span,NNE法的namda
pdf("../figure/05train.AUC.pdf",onefile = F)
par(mar=c(5,5,4,2))
plot(sur.fit1$FP, sur.fit1$TP, ## x=FP,y=TP
     type="l",lwd=4,col="#F7903D", ##线条设置
     xlim=c(0,1), ylim=c(0,1),
     cex.axis=2,cex.lab=2,
     xlab="False Positive",
     ylab="True Positive",
     main="The AUC of FARS")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色
lines(sur.fit2$FP, sur.fit2$TP, type="l",lwd=4,
      col="#4D85BD",
      xlim=c(0,1), ylim=c(0,1))
legend(0.3,0.3,c(paste("AUC of 3 YEAR =",sprintf("%0.3f",sur.fit1$AUC)),
                 paste("AUC of 5 YEAR =",sprintf("%0.3f",sur.fit2$AUC))),
       x.intersp=1, y.intersp=1,
       lty= 1 ,lwd= 4,col=c("#F7903D","#4D85BD","#59A95A"),
       bty = "n",# bty框的类型
       seg.len=1,cex=1.5)# 
dev.off()
# 风险因子三联图
library(ggrisk)
library(rms)
risk.dat <- train.dat %>% 
  dplyr::select(os,os.time,all_of(stepCox.genes)) %>% 
  dplyr::rename(time=os.time,status=os)
dd <- datadist(risk.dat)
options(datadist="dd")
fit <- cph(Surv(time,status)~CHSY1+LAMB1+PTTG1IP+RGS2+SPP1,risk.dat)
pdf("../figure/05train.ggrisk.pdf",height = 6,width = 6,onefile = F)
ggrisk(fit,
       color.A=c(low='#3581b9',high='#CF5246'),
       color.B=c(code.0='#3581b9',code.1='#CF5246'),
       color.C=c(low='#3581b9',median='white',high='#CF5246'))
dev.off()

# tcga
table(tcga.DLBC.clin$sample.type)
valid.tcga <- as.data.frame(tcga.DLBC.tpm) %>% 
  dplyr::filter(rownames(.)%in%lasso.coef$gene) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "barcode") %>% 
  inner_join(tcga.DLBC.clin[,c("barcode","os","os.time")],by="barcode") 
valid.tcga$pred=predict(fit.backward,type = 'risk',newdata = valid.tcga)
valid.tcga <- valid.tcga %>% 
  dplyr::mutate(pred=(pred-min(pred))/(max(pred)-min(pred)))
v.roc.backward <- future_lapply(1:5,function(x){surv(valid.tcga,"pred",12,x)},future.seed = NULL)
# 绘图
# 生存曲线图
library(ggplot2)
valid.tcga <- valid.tcga %>% dplyr::mutate(group=ifelse(pred>cut.point,"high","low"))
fit <- survfit(Surv(os.time,os)~group,data=valid.tcga)
pdf("../figure/05tcga.surv.pdf",onefile = F)
ggsurvplot(fit, data = valid.tcga,
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
# roc曲线图
sur.fit1 <- survivalROC(Stime=valid.tcga$os.time,##生存时间
                        status=valid.tcga$os,## 终止事件    
                        marker = valid.tcga[,"pred"], ## marker value    
                        predict.time = 12*3,## 预测时间截点
                        method="KM")##span,NNE法的namda
sur.fit2 <- survivalROC(Stime=valid.tcga$os.time,##生存时间
                        status=valid.tcga$os,## 终止事件    
                        marker = valid.tcga[,"pred"], ## marker value    
                        predict.time = 12*5,## 预测时间截点
                        method="KM")##span,NNE法的namda
pdf("../figure/05tcga.AUC.pdf",onefile = F)
par(mar=c(5,5,4,2))
plot(sur.fit1$FP, sur.fit1$TP, ## x=FP,y=TP
     type="l",lwd=4,col="#F7903D", ##线条设置
     xlim=c(0,1), ylim=c(0,1),
     cex.axis=2,cex.lab=2,
     xlab="False Positive",
     ylab="True Positive",
     main="The AUC of FARS")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色
lines(sur.fit2$FP, sur.fit2$TP, type="l",lwd=4,
      col="#4D85BD",
      xlim=c(0,1), ylim=c(0,1))
legend(0.3,0.3,c(paste("AUC of 3 YEAR =",sprintf("%0.3f",sur.fit1$AUC)),
                 paste("AUC of 5 YEAR =",sprintf("%0.3f",sur.fit2$AUC))),
       x.intersp=1, y.intersp=1,
       lty= 1 ,lwd= 4,col=c("#F7903D","#4D85BD","#59A95A"),
       bty = "n",# bty框的类型
       seg.len=1,cex=1.5)# 
dev.off()
# 风险因子三联图
library(ggrisk)
library(rms)
risk.dat <- valid.tcga %>% 
  dplyr::select(os,os.time,all_of(stepCox.genes)) %>% 
  dplyr::rename(time=os.time,status=os)
dd <- datadist(risk.dat)
options(datadist="dd")
fit <- cph(Surv(time,status)~CHSY1+LAMB1+PTTG1IP+RGS2+SPP1,risk.dat)
pdf("../figure/05tcga.ggrisk.pdf",height = 6,width = 6,onefile = F)
ggrisk(fit,
       color.A=c(low='#3581b9',high='#CF5246'),
       color.B=c(code.0='#3581b9',code.1='#CF5246'),
       color.C=c(low='#3581b9',median='white',high='#CF5246'))
dev.off()

save(stepCox.dat,valid.tcga,stepCox.genes,fit.backward,file = "05model.RData")
