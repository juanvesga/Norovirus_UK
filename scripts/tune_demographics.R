
# Call SEIAR Noro model using an 
rm(list = ls()) 
#remotes::install_github("mrc-ide/odin.dust", upgrade = FALSE)
library(odin.dust)
library(lubridate)
library(here)
library(socialmixr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggmatplot)
#################################
# Source scripts
vanHoek<-FALSE

source(here("scripts","setup.model.R"))
source(here("src","model_functions.R"))
seiar <- odin.dust::odin_dust("src/seiar.age.R")

########################################
########################################
##Fine tune age distribution through mortality 
mu0<-approx   (x=c(0,length(ages)), y = c(0.0005,0.003), xout=seq(1,length(ages),1), method="linear")
mu<-mu0$y*10000
sim<-run_demog_model(mu)

ntimes<-365*5
p<-sim$age.distr.sim
plot_age_distr(p,params)
idday<-seq(1,ntimes, 1/params$dt)
states<-sim$states
idx<-sim$idx

sample<- 
  states[idx$M,]+
  states[idx$G,]+
  states[idx$S,]+
  states[idx$E,]+
  states[idx$I,]+
  states[idx$A,]+
  states[idx$R,]
n_age<-as.data.frame(t(sample))
names(n_age) <-contact$demography$age.group
n_age$days<-idday
n_age.m<-reshape2::melt(n_age, id="days")
head(n_age.m)

n_all<- colSums(sample)

df<-data.frame(x=idday,y=n_all)

# Plot all population
p1<-ggplot(df, aes(x=x,y=y))+
  geom_line(col="purple",linewidth=2)+
  ylim(0,70e6)+
  labs(y="Total population")

gridExtra::grid.arrange(p1)


### Plot population by age group
p2<-ggplot(data=n_age.m, aes(x=days, y= value, col=variable))+
  geom_line()
gridExtra::grid.arrange(p2)
########################## 
# ###############


loglikelihood <- function(p){
  llk<-0
  noise<-rexp(length(p),1e6)
  if ( any(p < 0)){
    print(c("bad",t(p)))
    lsq<- Inf

  }else{
    sim<-run_demog_model(p)
    # llk<-sum(dbinom(round(contact$demography$proportion*1000),
    #                 1000,
    #                 noise+(1000*sim$age.distr.sim), log = TRUE),
    #          na.rm=T)*-1
    # 
    lsq<-sum(sqrt(round(contact$demography$proportion*1000)-noise+(1000*sim$age.distr.sim)))

  }
  return(lsq)

}
# 


fminsearch <- neldermead::fminsearch
opt <- neldermead::optimset(TolX = 1.e-2,  Display = "iter", MaxFunEvals =500)
x1 <- fminsearch(fun = loglikelihood,
                 x0 = mu,
                 # xmin =c(1,0.0001,2),
                 # xmax = c(10,1,5),
                 options = opt
)

x  <- c(x1$optbase$xopt)
mu1<-x/10000
write.csv(mu1, here("data",paste("mortality_MLE",".csv",sep = "")))

ssi<-run_demog_model(x)
plot_age_distr(ssi$age.distr.sim,params)
