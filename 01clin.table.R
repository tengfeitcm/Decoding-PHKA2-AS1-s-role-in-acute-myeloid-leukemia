setwd("~/wsx/projects/jlx/20220928/1-dat/")
require(tidyverse)
library(tableone)
load("01geo31312.cel.RData")
load("01tcga.DLBC.RData")

clin31312 <- clin31312 %>% 
  dplyr::mutate(gender=ifelse(gender=="F","female","male"),
                data="gse31312") %>% 
  dplyr::select(age,gender,stage,ecog,data)
tcga.DLBC.clin <- tcga.DLBC.clin %>% 
  dplyr::mutate(data="tcga") %>% 
  dplyr::select(age,gender,race,stage,data)
clin <- clin31312 %>% 
  bind_rows(tcga.DLBC.clin) %>% 
  mutate(stage=tolower(stage))

names(clin)
vars <- colnames(clin)
table1 <- CreateTableOne(vars = vars, strata = c("data"), data = clin)
table1 <- print(table1,
                showAllLevels = TRUE)
write.csv(table1,file="01clin.table.csv")
