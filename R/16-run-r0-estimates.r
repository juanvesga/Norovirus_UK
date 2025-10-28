rm(list = ls())
gc()  # garbage collection

# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root             <- here::here()
infile0          <- file.path(root,"output", "parameters_short.qs2")
infile_dataS     <- file.path(root,"output", "data_short2.qs2")
infile_datafit     <- file.path(root,"output", "data_short.qs2")
infile_dataL     <- file.path(root,"output", "data_long.qs2")
infile_dataPlots <- file.path(root, "output", "data_for_plots.qs2")
infile_input     <- file.path(root,"output", "params_list_short.qs2")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))
source(file.path(root, "R", "collect_function.R"))
source(file.path(root, "R", "r0_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(abind, include.only = c("abind"))
modify_attach(tidyr, include.only = c("gather"))

library(odin2)
library(dust2)
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

infile_runs      <- file.path(root,"output", paste0("mcmc_fits",model,".qs2"))


# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------

extract_iteration <- function(pars, i) {
  lapply(pars, function(param) param[[i]])
}



data_all   <- qs_read(infile_dataS)
data_fit   <- qs_read(infile_datafit)
data_plots <- qs_read(infile_dataPlots)
mcmc_runs  <- qs_read(infile_runs)
sims       <- mcmc_runs$runs
pars       <- mcmc_runs$parameters
nsamps     <- dim(sims)[2]
p          <- qs_read(infile0)

endsim<-24109
times<-seq(0,endsim)


nt<-24109#length(pars_list$fixed_pars$school_time)
# Rs_gi3_re<- array(0, c(nsamps,nt))  
# Rs_gi_re<- array(0, c(nsamps, nt))
# Rs_gii4_re<- array(0, c(nsamps,nt))
# Rs_gii_re<- array(0, c(nsamps,nt))
# 
# Rs_gi3_symp<- array(0, c(nsamps,nt))  
# Rs_gi_symp<- array(0, c(nsamps, nt))
# Rs_gii4_symp<- array(0, c(nsamps,nt))
# Rs_gii_symp<- array(0, c(nsamps,nt))
# 
# Rs_gi3_cross<- array(0, c(nsamps,nt))  
# Rs_gi_cross<- array(0, c(nsamps, nt))
# Rs_gii4_cross<- array(0, c(nsamps,nt))
# Rs_gii_cross<- array(0, c(nsamps,nt))
# season<- array(0, c(nsamps,nt))
# 
# 
# R0_gi3<- array(0, c(nsamps))  
# R0_gi<- array(0, c(nsamps))
# R0_gii4<- array(0, c(nsamps))
# R0_gii<- array(0, c(nsamps))




foreach_fun<-function(pars,
                      get_R0,
                      get_addresses,
                      extract_iteration,
                      .combine='comb', 
                      .multicombine=TRUE){
  
  
  foreach(iii = 1:nsamps) %dopar% {
    
   
    # Use it for a specific iteration
    pars_i <- extract_iteration(pars, iii)
    
  
    nt<-24109
    
   
    
    res1<-get_R0(pars_i,nt,get_addresses,"GI3")
    res2<-get_R0(pars_i,nt,get_addresses,"GI")
    res3<-get_R0(pars_i,nt,get_addresses,"GII4")
    res4<-get_R0(pars_i,nt,get_addresses,"GII")
    
  
    
    list(
      R0_gi3=res1$R0_symp , 
      R0_gi= res2$R0_symp,
      R0_gii4=res3$R0_symp,#*(1-0.2)
      R0_gii= res4$R0_symp,
      
      R_full_gi3=res1$R0_full , 
      R_full_gi= res2$R0_full,
      R_full_gii4=res3$R0_full,#*(1-0.2)
      R_full_gii= res4$R0_full,
      

      Rs_gi3_re=res1$Rs_re , 
      Rs_gi_re= res2$Rs_re,
      Rs_gii4_re=res3$Rs_re,
      Rs_gii_re= res4$Rs_re,
      
      Rs_gi3_symp=res1$Rs_symp,  
      Rs_gi_symp= res2$Rs_symp,
      Rs_gii4_symp=res3$Rs_symp,
      Rs_gii_symp= res4$Rs_symp,
      
      Rs_gi3_cross=res1$Rs_cross,  
      Rs_gi_cross= res2$Rs_cross,
      Rs_gii4_cross=res3$Rs_cross,
      Rs_gii_cross= res4$Rs_cross,
      
      season=as.numeric(res1$season)
    )
    
    
    #print(iii/nsamples)
    
  }
}


# Call cores and register clusters

ncores<-parallel::detectCores()

cl <- parallel::makeCluster(ncores-1)
doParallel::registerDoParallel(cl)


results<-foreach_fun(pars,
                     get_R0,
                     get_addresses,
                     extract_iteration)

# Stop cluster
parallel::stopCluster(cl)

# POst processing



R0_gi3 <- do.call(abind,c(lapply(results,function(x){x[[1]]}),along=1) )
R0_gi <- do.call(abind,c(lapply(results,function(x){x[[2]]}),along=1) )
R0_gii4 <- do.call(abind,c(lapply(results,function(x){x[[3]]}),along=1) )
R0_gii <- do.call(abind,c(lapply(results,function(x){x[[4]]}),along=1) )

Rs_gi3 <- do.call(abind,c(lapply(results,function(x){x[[5]]}),along=1) )
Rs_gi <- do.call(abind,c(lapply(results,function(x){x[[6]]}),along=1) )
Rs_gii4 <- do.call(abind,c(lapply(results,function(x){x[[7]]}),along=1) )
Rs_gii <- do.call(abind,c(lapply(results,function(x){x[[8]]}),along=1) )

Rs_re_gi3 <- do.call(abind,c(lapply(results,function(x){x[[9]]}),along=0) )
Rs_re_gi <- do.call(abind,c(lapply(results,function(x){x[[10]]}),along=0) )
Rs_re_gii4 <- do.call(abind,c(lapply(results,function(x){x[[11]]}),along=0) )
Rs_re_gii <- do.call(abind,c(lapply(results,function(x){x[[12]]}),along=0) )

Rs_symp_gi3 <- do.call(abind,c(lapply(results,function(x){x[[13]]}),along=0) )
Rs_symp_gi <- do.call(abind,c(lapply(results,function(x){x[[14]]}),along=0) )
Rs_symp_gii4 <- do.call(abind,c(lapply(results,function(x){x[[15]]}),along=0) )
Rs_symp_gii <- do.call(abind,c(lapply(results,function(x){x[[16]]}),along=0) )

Rs_cross_gi3 <- do.call(abind,c(lapply(results,function(x){x[[17]]}),along=0) )
Rs_cross_gi <- do.call(abind,c(lapply(results,function(x){x[[18]]}),along=0) )
Rs_cross_gii4 <- do.call(abind,c(lapply(results,function(x){x[[19]]}),along=0) )
Rs_cross_gii <- do.call(abind,c(lapply(results,function(x){x[[20]]}),along=0) )



season = do.call(abind,c(lapply(results,function(x){x[[21]]}),along=0) )



r0<-data.frame(rbind(
  quantile(R0_gi3,c(0.025,0.5,0.975)),
  quantile(R0_gi,c(0.025,0.5,0.975)),
  quantile(R0_gii4,c(0.025,0.5,0.975)),
  quantile(R0_gii,c(0.025,0.5,0.975))))



r0$Strain<-factor(c("GI.3","Other GI","GII.4","Other GII"),
                  levels=c("GI.3","Other GI","GII.4","Other GII"))

singleR0<- ggplot(r0,                                      # Grouped barplot using ggplot2
                  aes(x = Strain,
                      y = X50.,
                      fill = Strain)) +
  geom_bar(stat = "identity",position = "dodge", alpha=0.7)+
  geom_errorbar(data = r0, aes(x = Strain, ymin = X2.5., ymax = X97.5.),
                width = 0.2, position = position_dodge(.9)) +
  ylim(0,10)+
  labs(tag = "A",
       x = " ", y = "R0 (Basic Reproduction Number)") +
  scale_fill_manual(values=c('skyblue3','firebrick2','yellow3','grey28'))+
  theme_light() +  
  theme(
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.1, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )


#info<-list(index=idx_all)
#id<-index(info)
#labs<-names(id$state)
#idx<-seq(1,length(labs))

#names(idx)<-labs
#labs<-names(idx)
#dimnames(sims)[1]<-list(c(labs))

pop_all<-
  (sims["pop_all1",,]+
  sims["pop_all2",,]+
  sims["pop_all3",,]+
  sims["pop_all4",,]+
  sims["pop_all5",,]+
  sims["pop_all6",,]+
  sims["pop_all7",,]+
  sims["pop_all8",,]+
  sims["pop_all9",,]+
  sims["pop_all10",,])

sus_gi3<-  sims["sus_gi3",,]/pop_all
sus_gi<-   sims["sus_gi",,]/pop_all
sus_gii4<- sims["sus_gii4",,]/pop_all#(1-0.2)*sims["sus_gii4",,]/pop_all
sus_gii<-  sims["sus_gii",,]/pop_all

sus_gi3_re<-  sims["sus_gi3_re",,]/pop_all
sus_gi_re<-   sims["sus_gi_re",,]/pop_all
sus_gii4_re<- sims["sus_gii4_re",,]/pop_all#(1-0.2)*
sus_gii_re<-  sims["sus_gii_re",,]/pop_all

sus_gi3_cross<-  sims["sus_gi3_cross",,]/pop_all
sus_gi_cross<-   sims["sus_gi_cross",,]/pop_all
sus_gii4_cross<- sims["sus_gii4_cross",,]/pop_all
sus_gii_cross<-  sims["sus_gii_cross",,]/pop_all


sus_gi3  <- sus_gi3[,-1]
sus_gi   <- sus_gi[,-1]
sus_gii4 <- sus_gii4[,-1]
sus_gii  <- sus_gii[,-1]

sus_gi3_re  <- sus_gi3_re[,-1]
sus_gi_re   <- sus_gi_re[,-1]
sus_gii4_re <- sus_gii4_re[,-1]
sus_gii_re  <- sus_gii_re[,-1]

sus_gi3_cross  <- sus_gi3_cross[,-1]
sus_gi_cross   <- sus_gi_cross[,-1]
sus_gii4_cross <- sus_gii4_cross[,-1]
sus_gii_cross  <- sus_gii_cross[,-1]

# sus_gi3[,1] <-1
# sus_gi[,1]  <-1
# sus_gii4[,1]<-1
# sus_gii[,1] <-1
# 
# sus_gi3_re[,1] <-0
# sus_gi_re[,1]  <-0
# sus_gii4_re[,1]<-0
# sus_gii_re[,1] <-0
# 
# sus_gi3_cross[,1] <-0
# sus_gi_cross[,1]  <-0
# sus_gii4_cross[,1]<-0
# sus_gii_cross[,1] <-0

season<-season[1,]


S_gi3 <-(sus_gi3 +sus_gi3_re +sus_gi3_cross) 
S_gi  <-(sus_gi  +sus_gi_re  +sus_gi_cross) 
S_gii4<-(sus_gii4+sus_gii4_re+sus_gii4_cross) 
S_gii <-(sus_gii +sus_gii_re +sus_gii_cross) 

Reff_gi3 <-S_gi3  * as.numeric(Rs_gi3) 
Reff_gi  <-S_gi   * as.numeric(Rs_gi) 
Reff_gii4<-S_gii4 * as.numeric(Rs_gii4)
Reff_gii <-S_gii * as.numeric(Rs_gii)

Reff_gi3_s <-  t(Reff_gi3) * season                    
Reff_gi_s  <-  t(Reff_gi) * season                    
Reff_gii4_s<-  t(Reff_gii4) * season                    
Reff_gii_s <-  t(Reff_gii) * season                    

## aletarntive season R0

# Rs_gi3 <- as.data.frame(R0_gi3) 
# Rs_gi  <- as.data.frame(R0_gi) 
# Rs_gii4<- as.data.frame(R0_gii4)
# Rs_gii <- as.data.frame(R0_gii)
# 
# Reff_gi3_s <-  t(outer(Rs_gi3$R0_gi3, season))                    
# Reff_gi_s  <-  t(outer(Rs_gi$R0_gi, season))                   
# Reff_gii4_s<-  t(outer(Rs_gii4$R0_gii4, season))                    
# Reff_gii_s <-  t(outer(Rs_gii$R0_gii, season))   



Reff_gi3_s   <-  t(Rs_symp_gi3 * sus_gi3 + Rs_re_gi3 * sus_gi3_re  + Rs_cross_gi3 * sus_gi3_cross)                   
Reff_gi_s    <-  t(Rs_symp_gi * sus_gi + Rs_re_gi * sus_gi_re  + Rs_cross_gi * sus_gi_cross)                   
Reff_gii4_s  <-  t(Rs_symp_gii4 * sus_gii4 + Rs_re_gii4 * sus_gii4_re  + Rs_cross_gii4 * sus_gii4_cross)                   
Reff_gii_s   <-  t(Rs_symp_gii * sus_gii + Rs_re_gii * sus_gii_re  + Rs_cross_gii * sus_gii_cross)                   


probs <- c(0.025, 0.5, 0.975)

# Row quantiles
q_gi3  <- (rowQuantiles((Reff_gi3_s), probs = probs))* 1.4                   
q_gi   <- (rowQuantiles((Reff_gi_s), probs = probs)  ) *1.4
q_gii4 <- (rowQuantiles((Reff_gii4_s), probs = probs) ) 
q_gii  <- (rowQuantiles((Reff_gii_s), probs = probs)  )* 1.2

cols=c('skyblue3','firebrick2','yellow3','grey28')

plot(q_gi3[,2],type = "l",col=cols[1],ylim = c(0,6),xlim = c(20000,25000))
lines(q_gi[,2],col=cols[2])
lines(q_gii4[,2],col=cols[3])
lines(q_gii[,2],col=cols[4])
lines(c(0,length(season)),c(1,1),type = "l")

dates<-c(p$time_vec,p$time_vec[length(p$time_vec)]+1)-190
dates<-dates[-1]


df1<-data.frame(dates,(q_gi3))
df1$strain<-"GI.3"
names(df1)<-paste(c("date","low","mid","up","Strain"))

df2<-data.frame(dates,(q_gi))
df2$strain<-"Other GI"
names(df2)<-paste(c("date","low","mid","up","Strain"))

df3<-data.frame(dates,(q_gii4))
df3$strain<-"GII.4"
names(df3)<-paste(c("date","low","mid","up","Strain"))

df4<-data.frame(dates,(q_gii))
df4$strain<-"Other GII"
names(df4)<-paste(c("date","low","mid","up","Strain"))


df<-rbind(df1,df2,df3,df4)

df$Strain<-factor(df$Strain,levels = c("GI.3","Other GI","GII.4","Other GII")) 

Rt_long <- ggplot(df, aes(x=date, color=Strain) ) +
  geom_ribbon(aes(ymin = low, ymax = up, fill=Strain),alpha = 0.3,colour = NA)+
  geom_line(aes(y=mid),  lwd = 1) +
  geom_hline(yintercept=1,linetype="dashed", color = "navyblue",lwd = 1) +
  labs(tag = "B",x="", y=bquote(Seasonal-R[eff]))+
  ylim(0,2)  +
  scale_color_manual(values=c('skyblue3','firebrick2','yellow3','grey28'))+
  scale_fill_manual(values=c('skyblue3','firebrick2','yellow3','grey28'))+
  scale_x_date(date_breaks = "6 months", date_labels = "%b-%Y",
               limits = as.Date(c('2018-01-12','2025-05-31'))) +
  # theme_classic() +  
  # theme(
  #   legend.position = "none",
  #   panel.background = element_blank(),
  #   axis.text = element_text(colour = "black", size = 10),
  #   axis.title.y = element_text(size = 10, face = "bold"),
  #   legend.text = element_text(size = 9), legend.key = element_blank(),
  #   axis.text.x = element_text(angle = 60, hjust = 1),
  #   plot.tag = element_text(size = 12, face = "bold"),
  #   plot.tag.position = c(0.1, 0.95),
  #   panel.border = element_rect(colour = "black", fill=NA)
  # )
  theme_light() +  
  theme(
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.1, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )





gridExtra::grid.arrange(singleR0,Rt_long )


# Covid phase plot --------------------------------------------------------

covid_labs<-c(
  "Lockdown 1",
  "Lockdown 1 \n easing",
  "Reduced \n restrictions",
  "Schools \n open",
  "Lockdown 2",
  "Lockdown 2 \n easing",
  "Christmas \n easing",
  "Lockdown 3",
  "Lockdown 3 \n with schools open",
  "Lockdown 3 end"
)

# Lockdown 1 = 23rd March -  3rd June 2020
lock_start1<-epi_week("2020-03-23")
lock_end1  <-epi_week("2020-06-03")

# Lockdown 1 easing =  4th June - 29th July 2020
lock_start2<-epi_week("2020-06-04")
lock_end2  <-epi_week("2020-07-29")

# Reduced restrictions = 30th July - 3rd Sep 2020
lock_start3<-epi_week("2020-07-30")
lock_end3  <-epi_week("2020-09-03")

# Schools open = 4th Sept - 26th October 2020
lock_start4<-epi_week("2020-09-04")
lock_end4  <-epi_week("2020-11-04")

# Lockdown 2 = 5th November - 2nd December 2020
lock_start5<-epi_week("2020-11-05") #
lock_end5  <-epi_week("2020-12-02")

# Lockdown 2 easing = 3rd December - 19th December 2020
lock_start6<-epi_week("2020-12-03") 
lock_end6  <-epi_week("2020-12-19")

# Christmas = 20 December 2020 - 2nd January 2021
lock_start7<-epi_week("2020-12-20") 
lock_end7  <-epi_week("2021-01-04")#<- changed for continuity

# Lockdown 3 = 5th January - 8th March 2021
lock_start8<-epi_week("2021-01-05") 
lock_end8  <-epi_week("2021-03-08")

# Lockdown 3 with schools open = 8th March - 16th March 2021
lock_start9<-epi_week("2021-03-09") 
# lock_end9  <-epi_week("2021-03-16")
lock_end9  <-epi_week("2021-06-01")





Rt_covid <- ggplot(df, aes(x=x, y=value, color=strain) ) +
  geom_line(alpha = 0.6, lwd = 0.001) +
  geom_hline(yintercept=1,linetype="dashed", color = "red") +
  labs( x = "time (days)")+
  
  geom_vline(xintercept=lock_start1,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start2,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start3,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start4,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start5,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start6,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start7,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start8,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start9,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_end9,linetype="dashed", color = "purple") +
  
  ylab(bquote(Seasonal-R[eff]))+ 
  ylim(0,5)  +
  scale_color_manual(values=c('skyblue3','firebrick2','yellow2','grey28'))+
  scale_x_date(date_breaks = "1 months", date_labels = "%b-%Y",
               limits = as.Date(c('2019-12-01','2021-12-31'))) +
  theme_classic() +  
  theme(
    #legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.2, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )



inc_all<-
  sims["inc_day_gii4_1",,]+
  sims["inc_day_gii4_2",,]+
  sims["inc_day_gii4_3",,]+
  sims["inc_day_gii4_4",,]+
  sims["inc_day_gii4_5",,]+
  sims["inc_day_gii_1",,]+
  sims["inc_day_gii_2",,]+      
  sims["inc_day_gii_3",,]+      
  sims["inc_day_gii_4",,]+ 
  sims["inc_day_gii_5",,]+ 
  sims["inc_day_gi3_1",,]+ 
  sims["inc_day_gi3_2",,]+
  sims["inc_day_gi3_3",,]+
  sims["inc_day_gi3_4",,]+
  sims["inc_day_gi3_5",,]+ 
  sims["inc_day_gi_1",,]+
  sims["inc_day_gi_2",,]+
  sims["inc_day_gi_3",,]+
  sims["inc_day_gi_4",,]+
  sims["inc_day_gi_5",,] 


qtls_inc <- as.data.frame(
  rowQuantiles(t(inc_all),
               probs = c(0.025, 0.5, 0.975)
  )
)
qtls_inc$x<-c(p$time_vec[1]-1,p$time_vec) 

df<-qtls_inc



inc_covid <- ggplot(df, aes(x=x) ) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
  geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
  ylim(0,200000)  +
  labs(x = "Date", y="New daily symptomatic cases")+
  geom_vline(xintercept=lock_start1,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start2,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start3,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start4,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start5,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start6,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start7,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start8,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start9,linetype="dashed", color = "purple") +
  
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y",
               limits = as.Date(c('2017-12-01','2024-12-31'))) +
  theme_classic() +  
  theme(
    #legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0.2, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )

size_tag=3

inc_covid_short <- ggplot(df, aes(x=x) ) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
  geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
  ylim(0,200000)  +
  labs( x = "Date", y="New daily symptomatic cases")+
  geom_vline(xintercept=lock_start1,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start2,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start3,linetype="dashed", color = "purple") +
  #geom_vline(xintercept=lock_start4,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start5,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start6,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start7,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start8,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_start9,linetype="dashed", color = "purple") +
  geom_vline(xintercept=lock_end9,linetype="dashed", color = "purple") +
  annotate("text", x = lock_start1-30, y = 2e5, label = covid_labs[1],size = size_tag)+
  annotate("text", x = lock_start2-30, y = 1.75e5, label = covid_labs[2],size = size_tag)+
  annotate("text", x = lock_start3-30, y = 1.6e5, label = covid_labs[3],size = size_tag)+
  #annotate("text", x = lock_start4-30, y = 1.35e5, label = covid_labs[4],size = size_tag)+
  annotate("text", x = lock_start5-30, y = 2e5, label = covid_labs[5],size = size_tag)+
  annotate("text", x = lock_start6-30, y = 1.7e5, label = covid_labs[6],size = size_tag)+
  annotate("text", x = lock_start7-30, y = 1.25e5, label = covid_labs[7],size = size_tag)+
  annotate("text", x = lock_start8-30, y = 0.9e5, label = covid_labs[8],size = size_tag)+
  annotate("text", x = lock_start9-30, y = 0.7e5, label = covid_labs[9],size = size_tag)+
  annotate("text", x = lock_end9-30, y = 1e5, label = covid_labs[10],size = size_tag)+
  
  scale_x_date(date_breaks = "1 months", date_labels = "%b-%Y",
               limits = as.Date(c('2019-01-01','2021-07-30'))) +
  theme_classic() +  
  theme(
    #legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.2, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )

windows()
gridExtra::grid.arrange(inc_covid_short  )











