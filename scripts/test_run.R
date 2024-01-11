rm(list = ls()) 

#dev.off()

library(odin.dust)
library(here)
library(socialmixr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggmatplot)
library(profvis)
library(lubridate)

#################################
# 1 Source scripts
source(here("src","utility_functions.R"))
source(here("scripts","load_data.R"))
source(here("src","model_functions.R"))
source(here("src","make_transform_imm1par.R"))
source(here("scripts","setup.model.R"))
source(here("src","plot_single_fits.R"))



stochastic <- 0

#############
# 2 Create dust object 
#model_path<-here("src","seiar.age.2strain_alternative.R")
#model_path<-here("src","seiar.age.imm2par_noreinfection.R")
model_path<-here("src","seiar.age.immdrop_noreinfection.R")
seiar <- odin.dust::odin_dust(model_path)

#############
# 3 Create filter
if (stochastic==1){
  filter <- mcstate::particle_filter$new(data_all, 
                                         model = seiar, 
                                         n_particles = 10,
                                         compare = compare, 
                                         index = index)
}else{
  
  filter <- mcstate::particle_deterministic$new(data_all, 
                                                model = seiar, 
                                                compare = compare, 
                                                index = index)
  
  
  filter2 <- mcstate::particle_deterministic$new(data=data_empty,
                                                 model = seiar,
                                                 compare=compare_empty,
                                                 index = index)
}


fac<-1.5
  

theta=c(
    beta_1 =0.20901773*fac,    
    beta_2 =   0.22070250*fac ,
    beta_3 =   0.75916798*fac,   
    beta_4 =  0.46708090*fac,
    aduRR =      0.08377272,
    maternalAB =   81.59086447,  
    imm_yr = 19.61791224  ,
    imm_fac =  2.11062235,
    w1_1 =    0.31567480 ,
    repfac_0=159.98810378 ,   
    repfac_5= 510, 
    repfac_15=480, 
    repfac_65p= 40, 
    crossp_12=0.04876855,   
    crossp_21=0.05012397,    
    crossp_34=0.05404153,     
    crossp_43=0.05126302)  


# xvals<-(read.csv( here("output",paste("simplex_set",".csv",sep = ""))))
# xx<-xvals$x/unlist(params$scaling_fac)
# theta<-xx

pars<-list(
  beta_1  = theta[['beta_1']],   # transm coefficient
  beta_2  = theta[['beta_2']],   # transm coefficient
  beta_3  = theta[['beta_3']],   # transm coefficient
  beta_4  = theta[['beta_4']],   # transm coefficient
  maternalAB  = theta[['maternalAB']],   
  imm_yr = theta[['imm_yr']],
  imm_fac = theta[['imm_fac']],
  w1_1 = theta[['w1_1']],
  repfac_0=theta[['repfac_0']],
  repfac_5=theta[['repfac_5']],
  repfac_15=theta[['repfac_15']],
  repfac_65p=theta[['repfac_65p']],
  aduRR=theta[["aduRR"]],
  crossp_12=theta[['crossp_12']],
  crossp_21=theta[['crossp_21']],
  crossp_34=theta[['crossp_34']],
  crossp_43=theta[['crossp_43']],
  pop  = pop,
  init  = init,
  mu    = params$mu,
  m=params$transmission,
  m_holi=params$transmission_holi,
  cmx_1=params$cmx_1,
  cmx_2=params$cmx_2,
  cmx_3=params$cmx_3,
  cmx_4=params$cmx_4,
  cmx_5=params$cmx_5,
  cmx_6=params$cmx_6,
  cmx_7=params$cmx_7,
  cmx_8=params$cmx_8,
  cmx_9=params$cmx_9,
  aging_vec =params$aging_vec,
  school_step= as.double(params$school_uk),
  n_school_steps=params$n_school_steps,
  covid_step= as.double(params$covid_sche),
  n_covid_steps=params$n_covid_steps,
  N_age = params$N_age,
  vac_camp_cov=params$vac_camp_cov,
  vac_sche_cov=params$vac_sche_cov,
  dt=1
) # number of age groups


#pars2<-mapply(c, pars, pars, SIMPLIFY=FALSE)


## Run Filter
# 
filter2$run(pars, save_history = TRUE)
sims<-filter2$history()
data<-data_all
plot_single_fits(sims,data)



#### Other simulations(Stochastic)


reps<-2
model <- seiar$new(pars, 0, reps)
n_times<-length(days_vec)
tt<-seq(1,n_times,1)
sim<- model$simulate(tt)
sim<-drop(sim)

## Seasonality


idx<-model$info()$index
inci<-sim[idx$seasonality,1,]

plot(days_vec,inci*2, type = "l",
     xlim = c(as.Date("2014-01-01"),as.Date("2020-01-01")),
     ylim = c(0,10))
lines(total_cases$date,log(total_cases$cases),col="firebrick")



# 


## Cases week
idx<-model$info()$index
inci<-sim[idx$new_cases_week,1,]
plot(days_vec,inci, type = "l",
     xlim = c(as.Date("2012-01-01"),as.Date("2023-01-01")),
     ylim = c(0,500000)
     )



incgi<-sim[idx$new_cases_week_gi,1,]
incgi3<-sim[idx$new_cases_week_gi3,1,]
incgii<-sim[idx$new_cases_week_gii,1,]
incgii4<-sim[idx$new_cases_week_gii4,1,]

plot(days_vec,incgii4, type = "l",col="red",
     xlim = c(as.Date("2002-01-01"),as.Date("2023-01-01")))
lines(days_vec,incgii,col="black")
lines(days_vec,incgi3,col="purple")
lines(days_vec,incgi,col="darkgreen")



p_31<-1 - exp(- (1-(incgi3[991]/(incgi[991]+incgi3[991]))))

incgi<-colSums(sim[idx$inc_day_gi,1,])
incgi3<-colSums(sim[idx$inc_day_gi3,1,])
incgii<-colSums(sim[idx$inc_day_gii,1,])
incgii4<-colSums(sim[idx$inc_day_gii4,1,])

plot(days_vec,incgii4, type = "l",col="red",
     xlim = c(as.Date("2002-01-01"),as.Date("2023-01-01")))
lines(days_vec,incgii,col="black")
lines(days_vec,incgi3,col="purple")
lines(days_vec,incgi,col="darkgreen")





comix_line<-sim[idx$covid_timeline,1,]
inci<-sim[idx$new_cases,1,]
plot(days_vec,inci, type = "l", 
     xlim = c(as.Date("2018-01-01"),as.Date("2023-01-01")),
     ylim = c(0,10000))
lines(days_vec,comix_line*1000,col="firebrick")


cases_day<-sim[idx$new_cases,1,]
cases_week<-sim[idx$new_cases_week,1,]
plot(days_vec,cases_day, type = "l", 
     xlim = c(as.Date("2018-11-01"),as.Date("2020-12-01")),
     ylim = c(0,100000))
lines(days_vec,cases_week,col="firebrick")


idd<-total_cases$day

cases_day<-sim[idx$new_cases,1,]

cases_week2<-sim[idx$new_cases_week,1,idd]
plot(total_cases$date,cases_week2, type = "l", 
     xlim = c(as.Date("2018-11-01"),as.Date("2019-03-01")),
     ylim = c(0,100000))
lines(days_vec,cases_week,col="firebrick")





prev<-sim[idx$seroprev_num[2:7],1,]/sim[idx$seroprev_den[2:7],1,]
t<-sim[idx$time,1,]
plot(t,prev[1,], type = "l", col = "blue",  
        xlab = "day", ylab = "prevalence", las = 1,
        ylim=c(0,1))
lines(t,prev[2,], col = "red")
lines(t,prev[3,], col = "orange")
lines(t,prev[4,], col = "green")
lines(t,prev[5,], col = "black")
lines(t,prev[6,], col = "purple")
lines(c(7852,7852),c(0,1))


idx<-model$info()$index
prev<-sim[idx$seroprev[2:7],1,]
t<-sim[idx$time,1,]
plot(t,prev[1,], type = "l", col = "blue",  
     xlab = "day", ylab = "prevalence", las = 1,
     ylim=c(0,1))
lines(t,prev[2,], col = "red")
lines(t,prev[3,], col = "orange")
lines(t,prev[4,], col = "green")
lines(t,prev[5,], col = "black")
lines(t,prev[6,], col = "purple")
lines(c(6935,6935),c(0,1))





# Incidence
## Weekly cass reported by UKHSA
cols <- c("#8c8cd9", "#e67300", "#d279a6", "#ff4d4d", "#999966",
          "#660000")

idx<-model$info()$index
cases<-(sim[idx$infections_day_gii,,])
t<-sim[idx$time,1,]
plot(t,cases[1,1,], type = "l", col = cols[1],  
        xlab = "day", ylab = "incidence", las = 1,ylim = c(0,250000))
lines(t,cases[2,1,], type = "l", col = cols[2])
lines(t,cases[3,1,], type = "l", col = cols[3])
lines(t,cases[4,1,], type = "l", col = cols[4])


idx<-model$info()$index
cases<- cbind(
  (sim[idx$inc_year_gi3[1],1,]+
          sim[idx$inc_year_gi[1],1,]+
          sim[idx$inc_year_gii4[1],1,]+
          sim[idx$inc_year_gii[1],1,])/sim[idx$pop_by4age[1],1,] ,
(sim[idx$inc_year_gi3[2],1,]+
  sim[idx$inc_year_gi[2],1,]+
  sim[idx$inc_year_gii4[2],1,]+
  sim[idx$inc_year_gii[2],1,])/sim[idx$pop_by4age[2],1,],
(sim[idx$inc_year_gi3[3],1,]+
  sim[idx$inc_year_gi[3],1,]+
  sim[idx$inc_year_gii4[3],1,]+
  sim[idx$inc_year_gii[3],1,])/sim[idx$pop_by4age[3],1,],
(sim[idx$inc_year_gi3[4],1,]+
  sim[idx$inc_year_gi[4],1,]+
  sim[idx$inc_year_gii4[4],1,]+
  sim[idx$inc_year_gii[4],1,])/sim[idx$pop_by4age[4],1,],
(sim[idx$inc_year_gi3[5],1,]+
  sim[idx$inc_year_gi[5],1,]+
  sim[idx$inc_year_gii4[5],1,]+
  sim[idx$inc_year_gii[5],1,])/sim[idx$pop_by4age[5],1,])


cases<-t(cases*1000)  

t<-sim[idx$time,1,]
plot(t,cases[1,], type = "l", col = cols[1],  
     xlab = "day", ylab = "incidence", las = 1)
lines(t,cases[2,], type = "l", col = cols[2])
lines(t,cases[3,], type = "l", col = cols[3])
lines(t,cases[4,], type = "l", col = cols[4])
lines(t,cases[5,], type = "l", col = cols[5])
lines(c(6935,6935),c(0,2000))
lines(c(0,12000),c(310,310))

cases[,6935]

##Inc calc2

idx<-model$info()$index
cases<-(sim[idx$inc_day_all,,])*365
t<-sim[idx$time,1,]
plot(t,cases[1,1,], type = "l", col = cols[1],  
     xlab = "day", ylab = "incidence", las = 1, ylim = c(0,1000))
lines(t,cases[2,1,], type = "l", col = cols[2])
lines(t,cases[3,1,], type = "l", col = cols[3])
lines(t,cases[4,1,], type = "l", col = cols[4])
lines(t,cases[5,1,], type = "l", col = cols[5])





idx<-model$info()$index
cases<-sim[idx$infections_day_k,,]
t<-sim[idx$time,1,]
matplot(t,t(cases), type = "l", col = "#ff4d4d",  
        xlab = "day", ylab = "incidence", las = 1,
        xlim = c(0,365*2))



#cases year
idx<-model$info()$index
cases_yr<-sim[idx$cases_year[1],,]
t<-sim[idx$time,1,]
matplot(t,t(cases_yr), type = "l", col = "#e673001A",  
        xlab = "day", ylab = "incidence", las = 1,
        xlim = c(0,365*3))
lines(c(1095,1095),c(0,1e6))


# SIR
sir_col <- c("#8c8cd9", "#ff4d4d", "#999966","#e67300")
sir_col_transp <- paste0(sir_col, "1A")

M<-data.frame(M= t(colSums(sim[idx$M,,])))
G<-data.frame(G= t(colSums(sim[idx$G,,])))
S<-data.frame(S= t(colSums(sim[idx$S,,])))
E<-data.frame(I= t(colSums(sim[idx$E,,])))
I<-data.frame(I= t(colSums(sim[idx$I,,])))
A<-data.frame(R= t(colSums(sim[idx$A,,])))
R<-data.frame(R= t(colSums(sim[idx$R,,])))

sims<-(cbind(S,I,A,R))

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(t, sims, xlab = "Time", ylab = "Number of individuals",
        type = "l", 
        col = c(rep(sir_col_transp[1], reps),
                rep(sir_col_transp[2], reps),
                rep(sir_col_transp[3], reps),
                rep(sir_col_transp[4], reps)),
        lty = 1,
        xlim = c(0,365*2))
legend("right", lwd = 1, col = sir_col,
       legend = c("S", "I", "A",  "R"), bty = "n")

## Reproductive number
m <- params$contact$matrix # age-structured contact matrix
ngm<-c_mat*0
beta<- theta[['beta']]


# Next Generation matrix
for (i in 1:params$N_age){
  for (j in 1:params$N_age){
    ngm[i,j] <- beta *  m[i, j] * (params$theta + ((params$sigma+params$epsilon) * params$rho))
    
  }
}

r0<-eigen(ngm)
R0<-r0$values[1]

Rt<- R0*(S/(M+G+S+E+I+A+R))

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(t, Rt, xlab = "Time", ylab = "Effective Reproductive number",
        type = "l", 
        col = "#6600001A" ,
        lty = 1,
        ylim = c(0,R0))


## Seaonlaity
# School vs data 
x<-seq(1,365,7.2)
length(x)
length(sgss$cases)
aux<-sim[idx$aux,,]*1800
t<-sim[idx$time,1,]
matplot(t,t(aux), type = "l", col = "#e673001A",  
        xlab = "day", 
        ylab = "incidence", 
        las = 1,
        xlim = c(0,365),
        ylim=c(0,260)
)
points(x, sgss$cases)
lines(school_sche$day.no.,school_sche$School*200, col="red")


#
idx<-model$info()$index
cases<-sim[idx$reported_wk,,]
t<-sim[idx$time,1,]
matplot(t,t(cases), type = "l", col = "#e673001A",  
        xlab = "day", 
        ylab = "incidence", 
        las = 1,
        xlim = c(sgss$day[1],sgss$day[length(sgss$day)]),
        ylim = c(0,max(sgss$cases))
)
points(sgss$day, sgss$cases)




idx<-model$info()$index
cases<-sim[idx$reported_wk0_4,,]
t<-sim[idx$time,1,]
matplot(t,t(cases), type = "l", col = "red",  
        xlab = "day", 
        ylab = "Cases under4", 
        las = 1)

idx<-model$info()$index
cases<-sim[idx$reported_wk5_65,,]
t<-sim[idx$time,1,]
matplot(t,t(cases), type = "l", col = "red",  
        xlab = "day", 
        ylab = "Cases 5-65", 
        las = 1)

idx<-model$info()$index
cases<-sim[idx$reported_wk65_p,,]
t<-sim[idx$time,1,]
matplot(t,t(cases), type = "l", col = "red",  
        xlab = "day", 
        ylab = "Cases 65p", 
        las = 1)









