rm(list = ls())
setwd("~/wsx/projects/jlx/20220928/9-nomogram/")
library(tidyverse)
library(data.table)
library(survival)
library(survivalROC)
require(survminer)
load("../1-dat/01geo31312.cel.RData")
load("../5-model/05model.RData")

##单因素函数
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
## nomogram总
nomo.clin <- clin31312 %>% 
  inner_join(stepCox.dat[,c("gsm","pred","group")],by="gsm") %>% 
  dplyr::mutate(gender=ifelse(gender=="M","male","female"),
                # ecog=ifelse(ecog%in%c("0","1"),"0/1","2/3/4"),
                ecog=str_c("ECOG ",ecog),
                # stage=ifelse(stage%in%c("stage I","stage II"),"stage I/II","stage III/IV"),
                ) %>% 
  dplyr::rename(subtype=immu.subtype) %>% 
  dplyr::select(3:10)

##准备BaSurv
BaSurv<-Surv(time = nomo.clin$os.time,event = nomo.clin$os)
##提取变量名
names(nomo.clin)
varnames<-colnames(nomo.clin)[c(1:4,7:8)]
varnames
##循环单因素
univar<-lapply(varnames,function(x){unicox(nomo.clin,x)})
univar<-do.call(rbind,univar)
univar<-na.omit(univar)
length(univar$Characteristics[univar$P.Value<0.05])

##多因素分析
factors <- paste0(unique(univar$Characteristics[univar$P.Value<0.05]),collapse = "+")
fml<-as.formula(paste0("BaSurv~",factors))
multicox<-coxph(fml,data = nomo.clin)
multisum<-summary(multicox)
multisum
multiname<-rownames(multisum$coefficients)
MHR<-round(multisum$coefficient[,2],2)
MPV<-round(multisum$coefficient[,5],3)
MCIL<-round(multisum$conf.int[,3],2)
MCIH<-round(multisum$conf.int[,4],2)
MCI<-paste0(MCIL,"-",MCIH)
Mcoef<-round(multisum$coefficients[,1],2)
MlowCoef<-round(log(multisum$conf.int[,3]),2)
MhighCoef<-round(log(multisum$conf.int[,4]),2)
mulcox<-data.frame("Characteristics"=multiname,
                   "coef"=Mcoef,
                   "lowCoef"=MlowCoef,
                   "highCoef"=MhighCoef,
                   "Hazard Ratio"=MHR,
                   "lowCI"=MCIL,
                   "highCI"=MCIH,
                   "CI95"=MCI,
                   "P Value"=MPV)
##绘制森林图
pdf("../figure/09forest.pdf",height=4,width=6,onefile = FALSE)
ggforest(multicox,
         # main = "Hazard ratio", # 设置标题
         cpositions = c(0.02, 0.12, 0.35), # 设置前三列的相对距离
         fontsize = 0.8, # 设置字体大小
         refLabel = "reference",
         noDigits = 2)
dev.off()
# nomogram绘图
##绘制列线图
library(rms)
dd <- datadist(nomo.clin)
options(datadist="dd")
BaSurv<-Surv(nomo.clin$os.time, nomo.clin$os)
# factors <- paste0(mulcox$Characteristics[mulcox$P.Value<0.05],collapse = "+")
# factors <- paste0(c("pred","pNstage"),collapse = "+")
# fml<-as.formula(paste0("BaSurv~",factors))
coxm <- cph(BaSurv~pred+ecog,x=T,y=T,data = nomo.clin,surv = T)
coxm
summary(coxm)

surv <- Survival(coxm)
surv1 <- function(x)surv(12*1,lp=x)
surv2 <- function(x)surv(12*2,lp=x)
surv3 <- function(x)surv(12*3,lp=x)

norm1 <- nomogram(coxm,fun = list(surv1,surv2,surv3),lp=F,
                  funlabel=c("1-Year Survival probability",
                             "2-Year Survival probability",
                             "3-Year Survival probability"),
                  maxscale = 100, fun.at = c("0.9","0.85","0.8","0.7","0.6","0.5","0.4","0.3","0.2","0.1"))
pdf("../figure/09nomogram.pdf",onefile = F,width = 8,height = 5)
plot(norm1,cex.var = 0.8)
dev.off()
##校正曲线
coxm <- cph(fml,x=T,y=T,data = nomo.clin,surv = T,time.inc = 12)
cal <- calibrate(coxm, cmethod = "KM", method = "boot", u=12*1, m=155, B=466)
pdf("../figure/09nomogram2.pdf",onefile = F)
plot(cal, lwd=2,lty=1,subtitles = F,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim=c(0.9,1),ylim=c(0.9,1),xlab = "Nomogram-Predicted Probability of 1 year",
     ylab = "Actual 1 year (Probability)",col=c(rgb(192,98,83,maxColorValue = 255)))
# lines(cal[,c("mean.predicted","KM")],type = "b", lwd=2,col=c(rgb(192,98,83,maxColorValue = 255)),pch=16)
abline(0,1,lty=3,lwd=2,col="grey")
dev.off()

coxm <- cph(fml,x=T,y=T,data = nomo.clin,surv = T,time.inc = 12*3)
cal <- calibrate(coxm, cmethod = "KM", method = "boot", u=12*3, m=155, B=466)
pdf("../figure/09nomogram3.pdf",onefile = F)
plot(cal, lwd=2,lty=1,subtitles = F,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim=c(0.4,1),ylim=c(0.4,1),xlab = "Nomogram-Predicted Probability of 3 year",
     ylab = "Actual 3 year (Probability)",col=c(rgb(192,98,83,maxColorValue = 255)))
# lines(cal[,c("mean.predicted","KM")],type = "b", lwd=2,col=c(rgb(192,98,83,maxColorValue = 255)),pch=16)
abline(0,1,lty=3,lwd=2,col="grey")
dev.off()

coxm <- cph(fml,x=T,y=T,data = nomo.clin,surv = T,time.inc = 12*5)
cal <- calibrate(coxm, cmethod = "KM", method = "boot", u=12*5, m=155, B=466)
pdf("../figure/09nomogram4.pdf",onefile = F)
plot(cal, lwd=2,lty=1,subtitles = F,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim=c(0.1,0.7),ylim=c(0.1,0.7),xlab = "Nomogram-Predicted Probability of 5 year",
     ylab = "Actual 5 year (Probability)",col=c(rgb(192,98,83,maxColorValue = 255)))
# lines(cal[,c("mean.predicted","KM")],type = "b", lwd=2,col=c(rgb(192,98,83,maxColorValue = 255)),pch=16)
abline(0,1,lty=3,lwd=2,col="grey")
dev.off()
