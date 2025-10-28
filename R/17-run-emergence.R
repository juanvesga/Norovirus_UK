rm(list = ls())
gc()  # garbage collection
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
infile0          <- file.path(root,"output", "parameters_short.qs2")
infile_dataS     <- file.path(root,"output", "data_short2.qs2")
infile_datafit   <- file.path(root,"output", "data_short.qs2")
infile_dataL     <- file.path(root,"output", "data_long.qs2")
infile_dataPlots <- file.path(root, "output", "data_for_plots.qs2")
infile_input     <- file.path(root,"output", "params_list_short.qs2")
infile_model0    <- file.path(root,"models", "model_simple.R")
infile_model2    <- file.path(root,"models", "model_2pars_emerg.R")
infile_model3    <- file.path(root,"models", "model_drop.R")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))
source(file.path(root, "R", "collect_function.R"))
source(file.path(root, "R", "r0_functions.R"))
source(file.path(root, "R", "run_emergence_scenarios.R"))


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
if (model=="_0"||model=="_4"||model=="_8"){
  
  source(infile_model0)
  
  
}else if((model=="_1"||model=="_5"||model=="_9")){
  
  source(infile_model1)
  
  
}else if((model=="_2"||model=="_6"||model=="_10")){
  
  source(infile_model2)
  
  
}else if((model=="_3"||model=="_7"||model=="_11")){
  
  source(infile_model3)
  
}

# Source Odin model

infile_runs      <- file.path(root,"output", paste0("mcmc_fits",model,".qs2"))
outfile         <- file.path(root, "output", paste0("emergence",model,".qs2"))
# Load Packages -------------------------------------------------------------

extract_iteration <- function(pars, i) {
  lapply(pars, function(param) param[[i]])
}



data_all   <- qs_read(infile_dataS)
data_fit   <- qs_read(infile_datafit)
data_plots <- qs_read(infile_dataPlots)
mcmc_runs  <- qs_read(infile_runs)
sims       <- mcmc_runs$runs
pars       <- mcmc_runs$parameters
init       <- mcmc_runs$initState
id        <- mcmc_runs$index
nsim       <- dim(sims)[2]
p          <- qs_read(infile0)
rm(mcmc_runs)

endsim<-24109
times<-seq(0,endsim)
step_yr<-365
nt<-24109#length(pars_list$fixed_pars$school_time)
tt<-seq(0,endsim)
times<-seq(tt[length(t)],tt[length(tt)]+(365*5),1)



# Pass GII4 states -----------------------------------------------------------



foreach_fun<-function(cross_p,
                      transm_r,
                      params,
                      id,
                      init,
                      nsamps,
                      noro_model,
                      times,
                      update_parameters=update_parameters,
                      collect_fn=collect_fn,
                      extract_fn=extract_iteration,
                      .combine='comb', 
                      .multicombine=TRUE){
  
  foreach(jj = 1:nsamps) %dopar% {
    
    
    pars        <- extract_fn(params,jj)
    
    init_state0 <- unlist(extract_fn(init,jj))
    init_state  <- init_state0
    
    age_cats<-10
    i_id<-function(ind,st){
      
      m2d<-matrix(ind,nrow=age_cats)
      return(m2d[,st])
      
    }
    
    gii4s<-c( i_id(id$E,3), i_id(id$I,3), i_id(id$A,3), i_id(id$R,3),
              i_id(id$E1,2), i_id(id$I1,2), i_id(id$A1,2), i_id(id$R1,2),
              i_id(id$E2,2), i_id(id$I2,2), i_id(id$A2,2), i_id(id$R2,2),
              i_id(id$E4,3), i_id(id$I4,3), i_id(id$A4,3), i_id(id$R4,3),
              i_id(id$E12,1), i_id(id$I12,1), i_id(id$A12,1), i_id(id$R12,1),
              i_id(id$E14,2), i_id(id$I14,2), i_id(id$A14,2), i_id(id$R14,2),
              i_id(id$E24,2), i_id(id$I24,2), i_id(id$A24,2), i_id(id$R24,2),
              i_id(id$E4rd,3), i_id(id$I4rd,3), i_id(id$A4rd,3), i_id(id$R4rd,3))
    
    
    
    init_state[c(i_id(id$E,4), i_id(id$I,4), i_id(id$A,4), i_id(id$R,4))]<-
      init_state[c(i_id(id$E,4), i_id(id$I,4), i_id(id$A,4), i_id(id$R,4))]+
      init_state[c(i_id(id$E,3), i_id(id$I,3), i_id(id$A,3), i_id(id$R,3))]
    
    init_state[c(i_id(id$E1,3), i_id(id$I1,3), i_id(id$A1,3), i_id(id$R1,3))] <-
      init_state[c(i_id(id$E1,3), i_id(id$I1,3), i_id(id$A1,3), i_id(id$R1,3))] +
      init_state[c(i_id(id$E1,2), i_id(id$I1,2), i_id(id$A1,2), i_id(id$R1,2))]
    
    init_state[c(i_id(id$E2,3), i_id(id$I2,3), i_id(id$A2,3), i_id(id$R2,3))] <-
      init_state[c(i_id(id$E2,3), i_id(id$I2,3), i_id(id$A2,3), i_id(id$R2,3))] +
      init_state[c(i_id(id$E2,2), i_id(id$I2,2), i_id(id$A2,2), i_id(id$R2,2))]
    
    init_state[c(i_id(id$E3,3), i_id(id$I3,3), i_id(id$A3,3), i_id(id$R3,3))] <-
      init_state[c(i_id(id$E3,3), i_id(id$I3,3), i_id(id$A3,3), i_id(id$R3,3))] +
      init_state[c(i_id(id$E4,3), i_id(id$I4,3), i_id(id$A4,3), i_id(id$R4,3))]
    
    init_state[c(i_id(id$E12,2), i_id(id$I12,2), i_id(id$A12,2), i_id(id$R12,2))] <-
      init_state[c(i_id(id$E12,2), i_id(id$I12,2), i_id(id$A12,2), i_id(id$R12,2))] +
      init_state[c(i_id(id$E12,1), i_id(id$I12,1), i_id(id$A12,1), i_id(id$R12,1))]
    
    init_state[c(i_id(id$E13,2), i_id(id$I13,2), i_id(id$A13,2), i_id(id$R13,2))] <-
      init_state[c(i_id(id$E13,2), i_id(id$I13,2), i_id(id$A13,2), i_id(id$R13,2))] +
      init_state[c(i_id(id$E14,2), i_id(id$I14,2), i_id(id$A14,2), i_id(id$R14,2))]
    
    init_state[c(i_id(id$E23,2), i_id(id$I23,2), i_id(id$A23,2), i_id(id$R23,2))] <-
      init_state[c(i_id(id$E23,2), i_id(id$I23,2), i_id(id$A23,2), i_id(id$R23,2))] +
      init_state[c(i_id(id$E24,2), i_id(id$I24,2), i_id(id$A24,2), i_id(id$R24,2))]
    
    init_state[c(i_id(id$E4rd,4), i_id(id$I4rd,4), i_id(id$A4rd,4), i_id(id$R4rd,4))]<-
      init_state[c(i_id(id$E4rd,4), i_id(id$I4rd,4), i_id(id$A4rd,4), i_id(id$R4rd,4))]+
      init_state[c(i_id(id$E4rd,3), i_id(id$I4rd,3), i_id(id$A4rd,3), i_id(id$R4rd,3))]
    
    
    init_state[gii4s]<-0
    
    #seed
    init_state[i_id(id$E,3)]<-1
    
    
    
    if(transm_r==1 && cross_p==1){
      
      initial=init_state0
    }else{
      
      initial=init_state
      
    }
    
    
    
    pars$crossp_GII   <- pars$crossp_GII*cross_p
    #pars$log_beta_3<- log( exp(pars$log_beta_3) * transm_r)
    pars$transm_emerg <- transm_r
    
    # Determinsitic object
    sys  <- dust2::dust_system_create(noro_model, pars, deterministic = FALSE, n_particles = 1)
    
    
    dust2::dust_system_set_state(sys, initial)
    
    runs1<- dust2::dust_system_simulate(sys, times)
    
    
    
    
    
    res<-data.frame(
      runs_gi3=colSums(runs1[id$inc_day_gi3,]),
      runs_gi =colSums(runs1[id$inc_day_gi,]),
      runs_gii4=colSums(runs1[id$inc_day_gii4,]),
      runs_gii=colSums(runs1[id$inc_day_gii,]),
      reported=colSums(runs1[id$reported_wk,]),
      
      
      runs_age1  =
        runs1[id$inc_day_gi3[1],] +
        runs1[id$inc_day_gi3[2],]+
        runs1[id$inc_day_gi[1],] +
        runs1[id$inc_day_gi[2],]+
        runs1[id$inc_day_gii4[1],]+
        runs1[id$inc_day_gii4[2],]+
        runs1[id$inc_day_gii[1],] +
        runs1[id$inc_day_gii[2],],
      
      runs_age2  =
        runs1[id$inc_day_gi3[3],] +
        runs1[id$inc_day_gi[3],] +
        runs1[id$inc_day_gii4[3],]+
        runs1[id$inc_day_gii[3],],
      
      runs_age3=
        runs1[id$inc_day_gi3[4],] +
        runs1[id$inc_day_gi[4],] +
        runs1[id$inc_day_gii4[4],]+
        runs1[id$inc_day_gii[4],],
      
      runs_age4=
        runs1[id$inc_day_gi3[5],] +
        runs1[id$inc_day_gi[5],] +
        runs1[id$inc_day_gii4[5],]+
        runs1[id$inc_day_gii[5],]
    )
    
    
    
  }
}

scenario_pipeline<-function(cross_p,
                            transm_r,
                            pars,
                            id,
                            init,
                            nsim,
                            noro_model,
                            times,
                            par_fun,
                            time_vec){
  
  
  # Call cores and register clusters
  
  ncores<-parallel::detectCores()
  
  cl <- parallel::makeCluster(ncores-1)
  doParallel::registerDoParallel(cl)
  
  
  results    <-par_fun(cross_p,
                       transm_r,
                       pars,
                       id,
                       init,
                       nsim,
                       noro_model,
                       times)
  
  # Stop cluster
  parallel::stopCluster(cl)   
  
  
  # Get all column names
  col_names <- colnames(results[[1]])
  
  # Create a list where each element is a dataframe for one column
  combined_list <- lapply(col_names, function(col) {
    data.frame(sapply(results, function(df) df[[col]]))
  })
  
  # Name the list elements
  names(combined_list) <- col_names
  
  
  gi3 <- as.data.frame(
    rowQuantiles(as.matrix(combined_list$runs_gi3),
                 probs = c(0.025, 0.5, 0.975))
  )
  gi <- as.data.frame(
    rowQuantiles(as.matrix(combined_list$runs_gi),
                 probs = c(0.025, 0.5, 0.975))
  )
  gii4 <- as.data.frame(
    rowQuantiles(as.matrix(combined_list$runs_gii4),
                 probs = c(0.025, 0.5, 0.975))
  )
  gii <- as.data.frame(
    rowQuantiles(as.matrix(combined_list$runs_gii),
                 probs = c(0.025, 0.5, 0.975))
  )
  
  
  gi3$day<-times
  gi$day<-times
  gii4$day<-times
  gii$day<-times
  
  
  gii4$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  gi3$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  gi$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  gii$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  
  
  
  ## By age
  
  a1 <- as.data.frame( 
    rowQuantiles(as.matrix(combined_list$runs_age1),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  a2 <- as.data.frame( 
    rowQuantiles(as.matrix(combined_list$runs_age2),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  a3 <- as.data.frame( 
    rowQuantiles(as.matrix(combined_list$runs_age3),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  a4 <- as.data.frame( 
    rowQuantiles(as.matrix(combined_list$runs_age4),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  
  a1$day<-times
  a2$day<-times
  a3$day<-times
  a4$day<-times
  
  
  a1$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  a2$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  a3$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  a4$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  
  
  
  res<-list(
    runs_gi3=combined_list$runs_gi3,
    runs_gi=combined_list$runs_gi,
    runs_gii4=combined_list$runs_gii4,
    runs_gii=combined_list$runs_gii,
    reported=combined_list$reported , 
    gi3=gi3,
    gi=gi,
    gii4=gii4,
    gii=gii,
    a1=a1,
    a2=a2,
    a3=a3,
    a4=a4
  )
}






scen0<-scenario_pipeline(1,1,pars,id, init, nsim, noro_model(),times, foreach_fun, p$time_vec)

scen1<-scenario_pipeline(0,0.95,pars,id, init, nsim, noro_model(),times, foreach_fun, p$time_vec)

scen2<-scenario_pipeline(1,0.95,pars,id, init, nsim, noro_model(),times, foreach_fun, p$time_vec)

scen3<-scenario_pipeline(0,1,pars,id, init, nsim, noro_model(),times, foreach_fun, p$time_vec)

scen4<-scenario_pipeline(1,1.1,pars,id, init, nsim, noro_model(),times, foreach_fun, p$time_vec)

scen5<-scenario_pipeline(0,1.25,pars,id, init, nsim, noro_model(),times, foreach_fun, p$time_vec)



res_emergence<-list(
  scen0=scen0,
  scen1=scen1,
  scen2=scen2,
  scen3=scen3,
  scen4=scen4,
  scen5=scen5
  
)

qs_save(res_emergence,outfile)




























