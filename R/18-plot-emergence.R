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




res_emergence<-qs_read(infile)






sce<-c(
  ' \nBaseline: no emergence',
  '0% cross-immunity,\n5% less transmissible',
  'Same cross-immunity,\n5% less transmissible',
  '0% cross-immunity,\nequally transmissible',
  '0% cross-immunity,\n10% more transmissible',
  '0% cross-immunity,\n25% more transmissible'
)




#############  by strain

df<-data.frame(c(
  res_emergence$scen0$gi3$`50%`,
  res_emergence$scen0$gi$`50%`,
  res_emergence$scen0$gii4$`50%`,
  res_emergence$scen0$gii$`50%`
))



get_scen_cases<-function(scen){
  
  df<-data.frame(c(
    scen$gii$`50%`, # GI3
    scen$gii$`50%`,  # Other Gi
    scen$gii4$`50%`,
    scen$gii$`50%`
  ))
  
  df$strain<-factor(c(
    rep("GI3",length(scen$gi$`50%`)),
    rep("GI",length(scen$gi$`50%`)),
    rep("GII4",length(scen$gi$`50%`)),
    rep("GII",length(scen$gi$`50%`))
  ), levels =c("GI3","GI","GII4","GII"))
  
  df$date<-c(
    rep(scen$gi3$date,4))
  
  names(df)<-paste(c("cases","strain","date"))
  
  return(df)
}

df0<-get_scen_cases( res_emergence$scen0)
df1<-get_scen_cases(res_emergence$scen1)
df2<-get_scen_cases(res_emergence$scen2)
df3<-get_scen_cases(res_emergence$scen3)
df4<-get_scen_cases(res_emergence$scen4)
df5<-get_scen_cases(res_emergence$scen5)


get_scen_plot<-function(df,tag,ii){  

  if (ii>4){
    y_up<- 6e4
  } else{
    
    y_up<- 6e4
    
  }
  
  gp<-ggplot(df, aes(x = date, y = cases, fill = strain, width=1)) + 
    geom_bar(stat = "identity")+
    scale_fill_manual(values=c('skyblue3','firebrick2','yellow2','grey28'))+
    scale_x_date(date_breaks = "4 months", date_labels = "%b",
                 limits = as.Date(c('2025-12-01','2029-08-31'))) +
    labs(y="Cases",x="",title=tag)+
    
    ylim(c(0,y_up))+
    
    theme_minimal() +  
    theme(
      #legend.position = "none",
      title.text = element_text(size = 10),
      panel.background = element_blank(),
      axis.text = element_text(colour = "black", size = 10),
      axis.title.y = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9), legend.key = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, size=8),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0.1, 0.95),
      panel.border = element_rect(colour = "black", fill=NA)
    )
  
  return(gp)
  
}

p0<-get_scen_plot(df0,sce[1],1)
p1<-get_scen_plot(df1,sce[2],2)
p2<-get_scen_plot(df2,sce[3],3)
p3<-get_scen_plot(df3,sce[4],4)
p4<-get_scen_plot(df4,sce[5],5)
p5<-get_scen_plot(df5,sce[6],6)

windows()
gridExtra::grid.arrange( p0,p1,p2,p3,p4,p5,
                         layout_matrix = rbind(c(1, 2),
                                               c(3, 4),
                                               c(5, 6)))


## By age

get_scen_cases_age<-function(scen){
  
  df<-data.frame(c(
    scen$a1$`50%`, # GI3
    scen$a2$`50%`,  # Other Gi
    scen$a3$`50%`,
    scen$a4$`50%`
  ))
  
  df$age<-factor(c(
    rep("0_4",length(scen$a1$`50%`)),
    rep("5_15",length(scen$a1$`50%`)),
    rep("16_65",length(scen$a1$`50%`)),
    rep("65p",length(scen$a1$`50%`))
  ), levels =c("0_4","5_15","16_65","65p"))
  
  df$date<-c(
    rep(scen$a1$date,4))
  
  names(df)<-paste(c("cases","age","date"))
  
  return(df)
}

df0<-get_scen_cases_age(res_emergence$scen0)
df1<-get_scen_cases_age(res_emergence$scen1)
df2<-get_scen_cases_age(res_emergence$scen2)
df3<-get_scen_cases_age(res_emergence$scen3)
df4<-get_scen_cases_age(res_emergence$scen4)
df5<-get_scen_cases_age(res_emergence$scen5)


get_scen_plot_age<-function(df,tag,ii){  
  if (ii>4){
    y_up<- 16e4
  } else{
    
    y_up<- 10e4
    
  }
  gp<-ggplot(df, aes(x = date, y = cases, fill = age, width=1)) + 
    geom_bar(stat = "identity")+
    scale_fill_manual(values=c('blue3','aquamarine3','dodgerblue3','deepskyblue3'))+
    scale_x_date(date_breaks = "4 months", date_labels = "%b",
                 limits = as.Date(c('2025-12-01','2029-08-31'))) +
    labs(y="Cases",x="",title=tag)+
    ylim(c(0,y_up))+
    theme_minimal() +  
    theme(
      #legend.position = "none",
      title.text = element_text(size = 10),
      panel.background = element_blank(),
      axis.text = element_text(colour = "black", size = 10),
      axis.title.y = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9), legend.key = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, size=8),
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0.1, 0.95),
      panel.border = element_rect(colour = "black", fill=NA)
    )
  
  return(gp)
  
}

p0<-get_scen_plot_age(df0,sce[1],1)
p1<-get_scen_plot_age(df1,sce[2],2)
p2<-get_scen_plot_age(df2,sce[3],3)
p3<-get_scen_plot_age(df3,sce[4],4)
p4<-get_scen_plot_age(df4,sce[5],5)
p5<-get_scen_plot_age(df5,sce[6],6)

windows()
gridExtra::grid.arrange( p0,p1,p2,p3,p4,p5,
                         layout_matrix = rbind(c(1, 2),
                                               c(3, 4),
                                               c(5, 6)))


