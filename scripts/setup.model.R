
library(qs) # comix data
#######################################################################
## Set up necessary model structures and load data and parameters
############### Demographic model #######################################
vanHoek<- TRUE

# Mort rates (from ONS life tables)
mort_rates<-c(   #read.csv(here("data",paste("mortality_MLE",".csv",sep = "")))
  0.0039640, #1
  0.0002937, #2
  0.0001537, #3
  0.0001155, #4
  0.0000924, #5
  0.0000924, #6
  0.0000874, #7
  0.0000891, #7-15
  0.0002980, #16-25
  0.0005477, #26-35
  0.0011747, #36-45
  0.0026180, #46-55
  0.0064792, #56-65
  0.0163398) #66-75



# Load polymod

# we use 10 age bands 
#            1 2 3 4 5 6 7 8  9  10 11 12 13 14     
ages   <-  c(1,2,3,4,5,6,7,15,25,35,45,55,65,75) # Upper end of age bands



#ages   <-  seq(1,75,1) # Upper end of age bands
age_sq <- ages-1
adults <-  seq(tail(which(ages<=5),1),length(ages),1) 
infa_id<- which(ages<5)
adult_id<-which(ages>=15)
da <- diff(c(0,ages))

aging_vec = c(1/diff(ages), 0)
# find diagonal matrix for aging transitions
aging_mat <- diag(-1/da)
aging_mat[row(aging_mat)-col(aging_mat)==1] <- 1/head(da,-1)
aging_mat[length(ages),length(ages)]<-0 # Last age group aging accounted for with mortality
age.categories <- as.factor(ages)

#Load age contact matrix for Infants in UK according to vanHoek supplement





data(polymod, package = "socialmixr") # POLYMOD for all other contacts

contact = contact_matrix(
  polymod, countries = "United Kingdom", 
  age.limits = c(0,as.numeric(as.character(age.categories[1:(length(ages)-1)]))),
  symmetric = TRUE)

contact_holi = contact_matrix(
  polymod, countries = "United Kingdom", 
  age.limits = c(0,as.numeric(as.character(age.categories[1:(length(ages)-1)]))),
  filter = list(cnt_school=0),
  symmetric = TRUE)




# replace with vanHoek infants matrix

if(vanHoek==TRUE){
  vh<-read.csv(here("data","vanHoek.csv"), header=TRUE)#, sep=,)
  
  vhoek<-contact$matrix*0
  vhoek.upper.age<-c(1, 4,9,14,19,24,29,34,39,44,49,54,59,64,69,100)
  cmat.upper.age<-ages
  
  for (ii in 1:length(ages)){
    
    ii.cmx<- min(which( vhoek.upper.age>=cmat.upper.age[ii]))
    
    for (jj in 1:length(ages)){
      
      jj.cmx<- min(which( vhoek.upper.age>=cmat.upper.age[jj]))
      
      
      
      vhoek[ii,jj]<-vh[ii.cmx,jj.cmx]
      
      
    }
  }
  
  
  # Make matrix symmetric
  x1<-contact$matrix
  x1[1,]<-vhoek[1,]
  x1[,1]<-vhoek[,1]
  x<-((x1+t(x1))/2)
  contact$matrix<-x
  
  x2<-contact_holi$matrix
  x2[1,]<-vhoek[1,]
  x2[,1]<-vhoek[,1]
  x<-((x2+t(x2))/2)
  contact_holi$matrix<-x
  
}

pop = rep(contact$demography$population) # Population numbers

# Matrix to input into transmission formula, note is corrected for pop size
transmission <- contact$matrix /
  rep(contact$demography$population, each = ncol(contact$matrix))

transmission_holi <- contact_holi$matrix /
  rep(contact$demography$population, each = ncol(contact$matrix))
# UK holiday schedule 
school_sche<- read.csv(here("data","uk_holidays.csv"), header=TRUE)#, sep=,)

id<-which(school_sche$Date=="28-Feb")
school_year<-school_sche$School
val<-school_year[id]
school_year_leap <-c(school_year[1:id],            # Applying rbind function
                     val,
                     school_year[- (1:id)])

# School year vector accounting for leap years (2004 & 2008)
school_uk<-c(
  school_year, school_year, # 2002, 2003
  school_year_leap,         # 2004
  school_year, school_year, school_year, # 2005-2007
  school_year_leap,         # 2008
  school_year, school_year, school_year, # 2009-2011
  school_year_leap, # 2012
  school_year, school_year, school_year, # 2013-2015
  school_year_leap, # 2016
  school_year, school_year, school_year, # 2017-2019
  school_year_leap, # 2020
  school_year, school_year, school_year, # 2021-2023
  school_year_leap, # 2024
  school_year, school_year, school_year, # 2025-2027
  school_year_leap, # 2028
  school_year, school_year, school_year, # 2029-2031
  school_year_leap, # 2032
  school_year, school_year, school_year, # 2033-2035
  school_year_leap)


######### COMIX 9 periods of lockdown






comix<- read.csv(here("data","contact_matrices_9_periods.csv"), header=TRUE)#, sep=,)
comix[,1]<-NULL
comix.period <- names(table(comix$period))
names(comix) <- c("contactee","contact","vals","period")  # so can be used in the function
head(comix)
p1<-t(matrix(comix[comix$period =="1. Lockdown 1",3], ncol = length(unique(comix$contactee))))
p2<-t(matrix(comix[comix$period =="2. Lockdown 1 easing",3], ncol = length(unique(comix$contactee))))
p3<-t(matrix(comix[comix$period =="3. Relaxed restrictions",3], ncol = length(unique(comix$contactee))))
p4<-t(matrix(comix[comix$period =="4. School reopening",3], ncol = length(unique(comix$contactee))))
p5<-t(matrix(comix[comix$period =="5. Lockdown 2",3], ncol = length(unique(comix$contactee))))
p6<-t(matrix(comix[comix$period =="6. Lockdown 2 easing",3], ncol = length(unique(comix$contactee))))
p7<-t(matrix(comix[comix$period =="7. Christmas",3], ncol = length(unique(comix$contactee))))
p8<-t(matrix(comix[comix$period =="8. Lockdown 3",3], ncol = length(unique(comix$contactee))))
p9<-t(matrix(comix[comix$period =="9. Lockdown 3 + schools",3], ncol = length(unique(comix$contactee))))

cmx_1<-transmission*0
cmx_2<-transmission*0
cmx_3<-transmission*0
cmx_4<-transmission*0
cmx_5<-transmission*0
cmx_6<-transmission*0
cmx_7<-transmission*0
cmx_8<-transmission*0
cmx_9<-transmission*0
comix.upper.age<-c(4,11,17,29,39,49,59,69,100)
cmat.upper.age<-ages

for (ii in 1:length(ages)){
  
  ii.cmx<- min(which(comix.upper.age>=cmat.upper.age[ii]))
  
  for (jj in 1:length(ages)){
    
    jj.cmx<- min(which(comix.upper.age>=cmat.upper.age[jj]))
    
    
    
    cmx_1[ii,jj]<-p1[ii.cmx,jj.cmx]
    cmx_2[ii,jj]<-p2[ii.cmx,jj.cmx]
    cmx_3[ii,jj]<-p3[ii.cmx,jj.cmx]
    cmx_4[ii,jj]<-p4[ii.cmx,jj.cmx]
    cmx_5[ii,jj]<-p5[ii.cmx,jj.cmx]
    cmx_6[ii,jj]<-p6[ii.cmx,jj.cmx]
    cmx_7[ii,jj]<-p7[ii.cmx,jj.cmx]
    cmx_8[ii,jj]<-p8[ii.cmx,jj.cmx]
    cmx_9[ii,jj]<-p9[ii.cmx,jj.cmx]
    
    
  }
}

# tmp <- qread(here("data","4_Schools_return_cms.qs"))
# tmp1 <- apply(tmp,1,median)
# comix[comix$period=="4. School reopening",]$vals/tmp1


#### Check contact matrices
# age_labs<-c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75")
# c1<-contact$matrix
# rownames(c1)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(c1)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_1)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_1)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_2)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_2)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_3)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_3)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_4)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_4)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_5)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_5)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_6)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_6)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_7)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_7)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_8)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_8)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# rownames(cmx_9)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# colnames(cmx_9)<-paste(c("a1","a2","a3","a4","a5","a6","a7","a15","a25","a35","a45","a55","a65","a75"))
# 
# library(reshape2)
# melted_c1 <- melt(c1)
# 
# comix1<-cmx_1
# melted_comix1 <- melt(comix1)
# 
# comix2<-cmx_2
# melted_comix2 <- melt(comix2)
# 
# comix3<-cmx_3
# melted_comix3 <- melt(comix3)
# 
# comix4<-cmx_4
# melted_comix4 <- melt(comix4)
# 
# comix5<-cmx_5
# melted_comix5 <- melt(comix5)
# 
# comix6<-cmx_6
# melted_comix6 <- melt(comix6)
# 
# comix7<-cmx_7
# melted_comix7 <- melt(comix7)
# 
# comix8<-cmx_8
# melted_comix8 <- melt(comix8)
# 
# comix9<-cmx_9
# melted_comix9 <- melt(comix9)
# 
# library(ggplot2)
# ggplot(data = melted_c1, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()
# 
# 
# ggplot(data = melted_comix1, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix1")
# ggplot(data = melted_comix2, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix2")
# ggplot(data = melted_comix3, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix3")
# ggplot(data = melted_comix4, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix4")
# ggplot(data = melted_comix5, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix5")
# ggplot(data = melted_comix6, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix6")
# ggplot(data = melted_comix7, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix7")
# ggplot(data = melted_comix8, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix8")
# ggplot(data = melted_comix9, aes(x=Var1, y=contact.age.group, fill=value)) + 
#   geom_tile()+
#   ggtitle("comix9")
# 
# 
# sum(c1)
# sum(cmx_1)
# sum(cmx_2)
# sum(cmx_3)
# sum(cmx_4)
# sum(cmx_5)
# sum(cmx_6)
# sum(cmx_7)
# sum(cmx_8)
# sum(cmx_9)


################################
# 
# devtools::install_github('Mikata-Project/ggthemr')
# library(patchwork)
# library(ggthemr)
# library(socialmixr)
# library(data.table)
# 
# ggthemr("fresh")
# ## Create Figure of Contact Matrices
# 
# source('r/functions/sm_to_gg_matrix.R')
# source('r/functions/utility_functions.R')
# 
# h2020_cm <- readRDS('data/contact_matrices/h2020_cm.rds')
# h2020_cm_imputed <- readRDS('data/contact_matrices/h2020_cm_imputed.rds')
# polymod_cm <- readRDS('data/contact_matrices/polymod_cm.rds')
# 
# 
# rowSums(polymod_cm)/rowSums(h2020_cm)
# colnames(polymod_cm)
# 
# cmatrices <- list(
#   "POLYMOD" = c1,
#   "CoMix" = cmx_1
# )
# 
# ## Transform data
# cm_dt <- sm_to_gg_matrix(cmatrices)
# 
# cm_dt$contact_age
# 
# ## Create the plot
# matrix_plot <- gg_matrix(
#   cm_dt, 
#   breaks = c(0, 2,4,6, 8),
#   age_lab = age_labs
# ) + 
#   ggtitle("A") 
# 
# matrix_plot







###################################

## make symetric matrices 
# make_symme<-function(C){
#   mat<-((C+t(C))/2)
# }
#   
cmx_1<-symm_mat(cmx_1)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_2<-symm_mat(cmx_2)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_3<-symm_mat(cmx_3)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_4<-symm_mat(cmx_4)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_5<-symm_mat(cmx_5)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_6<-symm_mat(cmx_6)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_7<-symm_mat(cmx_7)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_8<-symm_mat(cmx_8)/rep(contact$demography$population, each = ncol(contact$matrix))
cmx_9<-symm_mat(cmx_9)/rep(contact$demography$population, each = ncol(contact$matrix))


## Lockdown dates
covid_start<-as.Date("2020-03-23")
covid_end  <-as.Date("2021-03-16")
covid_vec<-seq(covid_start,covid_end,"day")
id_day_covid<-which(days_vec%in%covid_vec)
covid_sche<-school_uk*0
covid_sche[id_day_covid]<-1
# plot(days_vec,school_uk,type = "l")
# lines(days_vec,covid_sche, col="red")

# Lockdown 1 = 23rd March -  3rd June 2020
lock_start1<-as.Date("2020-03-23")
lock_end1  <-as.Date("2020-06-03")
lock_vec<-seq(lock_start1,lock_end1,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-1

# Lockdown 1 easing =  4th June - 29th July 2020
lock_start2<-as.Date("2020-06-04")
lock_end2  <-as.Date("2020-07-29")
lock_vec<-seq(lock_start2,lock_end2,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-2

# Reduced restrictions = 30th July - 3rd Sep 2020
lock_start3<-as.Date("2020-07-30")
lock_end3  <-as.Date("2020-09-03")
lock_vec<-seq(lock_start3,lock_end3,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-3

# Schools open = 4th Sept - 26th October 2020
lock_start4<-as.Date("2020-09-04")
lock_end4  <-as.Date("2020-11-04")
lock_vec<-seq(lock_start4,lock_end4,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-4

# Lockdown 2 = 5th November - 2nd December 2020
lock_start5<-as.Date("2020-11-05") #
lock_end5  <-as.Date("2020-12-02")
lock_vec<-seq(lock_start5,lock_end5,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-5

# Lockdown 2 easing = 3rd December - 19th December 2020
lock_start6<-as.Date("2020-12-03") 
lock_end6  <-as.Date("2020-12-19")
lock_vec<-seq(lock_start6,lock_end6,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-6

# Christmas = 20 December 2020 - 2nd January 2021
lock_start7<-as.Date("2020-12-20") 
lock_end7  <-as.Date("2021-01-04")#<- changed for continuity
lock_vec<-seq(lock_start7,lock_end7,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-7

# Lockdown 3 = 5th January - 8th March 2021
lock_start8<-as.Date("2021-01-05") 
lock_end8  <-as.Date("2021-03-08")
lock_vec<-seq(lock_start8,lock_end8,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-8

# Lockdown 3 with schools open = 8th March - 16th March 2021
lock_start9<-as.Date("2021-03-09") 
# lock_end9  <-as.Date("2021-03-16")
lock_end9  <-as.Date("2021-06-01")
lock_vec<-seq(lock_start9,lock_end9,"day")
id_day_lock<-which(days_vec%in%lock_vec)
covid_sche[id_day_lock]<-9




# 
# sgssdata<-data.frame(cases=total_cases$cases, day=dates2)
# 
# plot(days_vec,school_uk,type = "l", xlim = as.Date(c("2012-01-01","2023-12-31")),
#      ylim = c(0,200))
# lines(days_vec,covid_sche*10,col='red')
# lines(sgssdata$day, sgssdata$cases)


########## Model parameters ##############################################
params<-list(
  
  age_sq=age_sq,  
  aging_mat=aging_mat/365, # get the day step aging
  aging_vec=aging_vec/365,# yearly aging rates 
  age.categories = age.categories,
  contact = contact,
  contact_holi=contact_holi,
  pop = pop,
  transmission =transmission,
  transmission_holi=transmission_holi,
  cmx_1=cmx_1,
  cmx_2=cmx_2,
  cmx_3=cmx_3,
  cmx_4=cmx_4,
  cmx_5=cmx_5,
  cmx_6=cmx_6,
  cmx_7=cmx_7,
  cmx_8=cmx_8,
  cmx_9=cmx_9,
  school_uk=school_uk,
  n_school_steps=length(school_uk),
  covid_sche=covid_sche,
  n_covid_steps=length(covid_sche),
  N_age = length(ages),
  age_select= c(1,seq(2,length(ages),1)*0),
  maternalAB = 25,  # maternal Ab decay (days)
  epsilon = 1,   # incubation
  theta = 2,   # duration symptoms
  sigma = 15, # duration asymp shedding (days)
  imm_yr   = 5.1 ,#(365*5.1),    # duration immunity
  rr_inf_asymp   = 0.05, # rel infect asymptomatic 
  p_nonsecretor=0.2, # Fraction immune genetically
  mu    = mort_rates/365,
  age_beta = 1+(seq(1,length(ages),1)*0),
  adult_id = adult_id,# indices for adult groups 
  infa_id=infa_id,
  aduRR=0.1, # cofactor of infectiousness for adults
  repfac= 287,# factor to amplify reported to community
  repfac_65= 287,# factor to amplify reported to community
  adult_beta=1, # adult with adult transmission scalar
  w1_1 = 0.8, # sesonality (if 0 not seasonal)
  w2 = -2, # 
  alpha = 1, # relative suscept in R compartment 
  # simulation
  dt=0.5,
  index_idd2=index_idd2,
  scaling_fac=list(
    beta_1 = 1000,
    beta_2 = 1000,
    beta_3 = 1000,
    beta_4 = 1000,
    aduRR=100,
    maternalAB = 1,
    imm_yr = 10,
    imm_fac = 100,
    repfac_0=0.25,
    repfac_5=0.25,
    repfac_15=0.25,
    repfac_65p=0.25,
    crossp_GI = 1000,
    crossp_GII = 1000
  ),
  
  vac_camp_cov= c(seq(1,length(ages),1)*0),
  vac_sche_cov= c(seq(1,length(ages),1)*0)
)

# Initial conditions of the model (M,G,S,E,I,A,R) x age cats X 4

#1:3 M,G,S,
#4:7 Ej,Ij,Aj,Rj
#8:10 Ekj,Ikj,Akj
#11:14 Ejk,Ijk,Ajk,Rjk
#15:18 Ek,Ik,Ak,Rk
#19:22 El,Il,Al,Rl
#23:25 Eml,Iml,Aml
#26:29 Elm,Ilm,Alm,Rlm
#30:33 Em,Im,Am,Rm

init<-matrix(0, nrow = 114, ncol = length(ages))
init[2,]<-round(params$pop*params$p_nonsecretor) # G
init[3,]<-params$pop - round(params$pop*params$p_nonsecretor)
init[5,5]<-1 # Seed in 15yrs old in Ij
init[9,5]<-1 # Seed in 15yrs old in Ik
init[13,5]<-1 # Seed in 15yrs old in Il
init[17,5]<-1 # Seed in 15yrs old in Im
#init[7,1]<-500


## transform params function


c1<-params$transmission
c2<-params$transmission_holi 
aduid<-params$adult_id 
cmx_1=params$cmx_1
cmx_2=params$cmx_2
cmx_3=params$cmx_3
cmx_4=params$cmx_4
cmx_5=params$cmx_5
cmx_6=params$cmx_6
cmx_7=params$cmx_7
cmx_8=params$cmx_8
cmx_9=params$cmx_9
covid_step= as.double(params$covid_sche)
n_covid_steps=params$n_covid_steps
mu<-params$mu 
school<-params$school_uk
n_school_steps<-params$n_school_steps
n_age<-params$N_age
aging_vec<-params$aging_vec
scalefc<-params$scaling_fac
pop<-params$pop
vac_camp_cov<-params$vac_camp_cov
vac_sche_cov<-params$vac_sche_cov


footransform <-  make_transform(c1,
                               c2, 
                               cmx_1,
                               cmx_2,
                               cmx_3,
                               cmx_4,
                               cmx_5,
                               cmx_6,
                               cmx_7,
                               cmx_8,
                               cmx_9,
                               aduid, 
                               mu, 
                               school,
                               n_school_steps,
                               covid_step,
                               n_covid_steps,
                               n_age,
                               aging_vec,
                               init,
                               pop,
                               scalefc,
                               vac_camp_cov,
                               vac_sche_cov)



## Clear space
rm(contact, 
   contact_holi, 
   vh,
   transmission,
   transmission_holi,
   fulldata_iid2,
   #   data_iid2.c4,
   #   sgss,
   #   sero,
   asymp,
   synth,
   agg_synth,
   d,
   newd,
   mort_rates,
   polymod,
   aging,
   p1,p2,p3,p4,p5,p6,p7,p8,p9)#,
#school_sche)





