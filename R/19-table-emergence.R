rm(list = ls())
gc()
# Notes on indexing the state
# index 1 is GI3 infection
# index 2 is other Gi infection
# index 3 is GII4 infection
# index 4 is other GII infection
# For states E, I and A the first number is the active infecyive strain
# the following numbers are the carrying immunity
# R compartments show carrying immunity in ascending order
# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root             <- here::here()

source(file.path(root, "R", "modify_attach.R"))
# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(abind, include.only = c("abind"))
modify_attach(tidyr, include.only = c("gather"))


library(ggplot2)
library(dplyr)
library(foreach)
library(matlib)
library(matrixStats)
library(tidyr)


# Load Packages -----------------------------------------------------------


library(here)
library(ggmatplot)
library(lubridate)
library(mcstate)
library(matrixStats)

# Select model ------------------------------------------------------------


# Cross protection 5%

#model       <- "_0" # Simple SEIAR
#model       <- "_1" # Simple SEIAR no reinf no cross-prot
model       <- "_2" # Full version with reinf and cross-protection
#model      <- "_3" # Full version with reinf and cross-protection and drop immunity


# Cross protection 25%

#model      <- "_4" # Simple SEIAR
#model      <- "_5" # Simple SEIAR no reinf no cross-prot
#model      <- "_6" # Full version with reinf and cross-protection
#model      <- "_7" # Full version with reinf and cross-protection and drop immunity


# Cross protection 50%

#model       <- "_8" # Simple SEIAR
#model      <- "_9" # Simple SEIAR no reinf no cross-prot
#model      <- "_10" # Full version with reinf and cross-protection
#model      <- "_11" # Full version with reinf and cross-protection and drop immunity

# Source Odin model

infile  <- file.path(root, "output", paste0("emergence",model,".qs2"))
outfile <- file.path(root, "output", paste0("table_emergence",model,".csv"))


# Set outputs path



df<-qs_read(infile)

tt<-5*365


base_all     <-as.matrix(df$scen0$runs_gii4)
tmp_all      <-as.matrix(df$scen0$runs_gii4)
tmp_gii4     <-df$scen0$gii4
tmp<-(colCumsums(tmp_all)-colCumsums(base_all))/colCumsums(tmp_all)
rel<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))

scen0  <-c(cumsum(tmp_gii4$`50%`)[tt],
           cumsum(tmp_gii4$`2.5%`)[tt],
           cumsum(tmp_gii4$`97.5%`)[tt], 
           NA,
           NA,
           NA)

#Scen1
tmp_all      <-as.matrix(df$scen1$runs_gii4)
tmp_gii4     <-df$scen1$gii4
tmp<-(colCumsums(tmp_all)-colCumsums(base_all))/colCumsums(tmp_all)
rel<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))

scen1  <-c(cumsum(tmp_gii4$`50%`)[tt],
           cumsum(tmp_gii4$`2.5%`)[tt],
           cumsum(tmp_gii4$`97.5%`)[tt], 
           rel[tt,2]*100,
           rel[tt,1]*100,
           rel[tt,3]*100)

#Scen2
tmp_all      <-as.matrix(df$scen2$runs_gii4)
tmp_gii4     <-df$scen2$gii4
tmp<-(colCumsums(tmp_all)-colCumsums(base_all))/colCumsums(tmp_all)
rel<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))

scen2  <-c(cumsum(tmp_gii4$`50%`)[tt],
           cumsum(tmp_gii4$`2.5%`)[tt],
           cumsum(tmp_gii4$`97.5%`)[tt], 
           rel[tt,2]*100,
           rel[tt,1]*100,
           rel[tt,3]*100)

#Scen3
tmp_all      <-as.matrix(df$scen3$runs_gii4)
tmp_gii4     <-df$scen3$gii4
tmp<-(colCumsums(tmp_all)-colCumsums(base_all))/colCumsums(tmp_all)
rel<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))

scen3  <-c(cumsum(tmp_gii4$`50%`)[tt],
           cumsum(tmp_gii4$`2.5%`)[tt],
           cumsum(tmp_gii4$`97.5%`)[tt], 
           rel[tt,2]*100,
           rel[tt,1]*100,
           rel[tt,3]*100)

#Scen4
tmp_all      <-as.matrix(df$scen4$runs_gii4)
tmp_gii4     <-df$scen4$gii4
tmp<-(colCumsums(tmp_all)-colCumsums(base_all))/colCumsums(tmp_all)
rel<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))

scen4  <-c(cumsum(tmp_gii4$`50%`)[tt],
           cumsum(tmp_gii4$`2.5%`)[tt],
           cumsum(tmp_gii4$`97.5%`)[tt], 
           rel[tt,2]*100,
           rel[tt,1]*100,
           rel[tt,3]*100)


#Scen5
tmp_all      <-as.matrix(df$scen5$runs_gii4)
tmp_gii4     <-df$scen5$gii4
tmp<-(colCumsums(tmp_all)-colCumsums(base_all))/colCumsums(tmp_all)
rel<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))

scen5  <-c(cumsum(tmp_gii4$`50%`)[tt],
           cumsum(tmp_gii4$`2.5%`)[tt],
           cumsum(tmp_gii4$`97.5%`)[tt], 
           rel[tt,2]*100,
           rel[tt,1]*100,
           rel[tt,3]*100)


tab<-rbind(scen0,scen1,scen2,scen3,scen4,scen5)




# by season realtive ------------------------------------------------------
ids<-which(epiweek(df$scen0$gii4$date)==27)
start<-ids[1]
df$scen0$gii4$date[ids]


get_table<-function(df_base,df_scen,wks=ids){
  
  
  start<-wks[1]
  season1<-c(start:wks[8])
  season2<-c(wks[8]:wks[8+7])
  season3<-c(wks[8+7]:wks[8+14])
  
  base_all     <-as.matrix(df_base$runs_gi3+df_base$runs_gi+df_base$runs_gii4+df_base$runs_gii)
  tmp_all      <-as.matrix(df_scen$runs_gi3+df_scen$runs_gi+df_scen$runs_gii4+df_scen$runs_gii)
  tmp_gii4     <-as.matrix(df_scen$runs_gii4)
  
  tmp<-(colCumsums(tmp_gii4[season1,]))/colCumsums(tmp_all[season1,])
  frac_gii4_s1<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))
  
  tmp<-(colCumsums(tmp_gii4[season2,]))/colCumsums(tmp_all[season2,])
  frac_gii4_s2<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))
  
  tmp<-(colCumsums(tmp_gii4[season3,]))/colCumsums(tmp_all[season3,])
  frac_gii4_s3<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))
  
  tmp<-(colCumsums(tmp_all[season1,])-colCumsums(base_all[season1,]))/colCumsums(tmp_all[season1,])
  rel_s1<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))
  
  tmp<-(colCumsums(tmp_all[season2,])-colCumsums(base_all[season2,]))/colCumsums(tmp_all[season2,])
  rel_s2<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))
  
  tmp<-(colCumsums(tmp_all[season3,])-colCumsums(base_all[season3,]))/colCumsums(tmp_all[season3,])
  rel_s3<-as.data.frame(rowQuantiles((tmp),probs = c(0.025, 0.5, 0.975)))
  
  
  out  <-c(rel_s1[365,2],
           rel_s1[365,1],
           rel_s1[365,3],
           frac_gii4_s1[365,2],
           frac_gii4_s1[365,1],
           frac_gii4_s1[365,3],
           rel_s2[365,2],
           rel_s2[365,1],
           rel_s2[365,3],
           frac_gii4_s2[365,2],
           frac_gii4_s2[365,1],
           frac_gii4_s2[365,3],
           rel_s3[365,2],
           rel_s3[365,1],
           rel_s3[365,3],
           frac_gii4_s3[365,2],
           frac_gii4_s3[365,1],
           frac_gii4_s3[365,3]
  )
  
  
  return(out)
}

scen0<-get_table(df$scen0,df$scen0)
scen1<-get_table(df$scen0,df$scen1)
scen2<-get_table(df$scen0,df$scen2)
scen3<-get_table(df$scen0,df$scen3)
scen4<-get_table(df$scen0,df$scen4)
scen5<-get_table(df$scen0,df$scen5)



tab2<-rbind(scen0,scen1,scen2,scen3,scen4,scen5)
write.csv(tab2,outfile)
