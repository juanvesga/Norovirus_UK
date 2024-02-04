#############
## PMCMC
#############




# Define  priors 
scaling_fac<-params$scaling_fac

if (start_from_best){
  
  
  # Load previous chains
  if (stochastic==1){
    load(here("output","processed_samples.RData")) 
  }else{
    
    
    load(here("output",paste0("chains_",scenario,".RData"))) 
    
    # xvals<-read.csv( here("output",paste("simplex_set",".csv",sep = "")))
    
  }
  
  id<- which(processed_chains$probabilities[,"log_posterior"] == 
               max(processed_chains$probabilities[,"log_posterior"]))[1]
  thetas<-processed_chains$pars
  
  mean_hpd <- apply(thetas, 2, mean)
  last_best <- thetas[id,]
  rm(processed_chains)
  # last_best2<-xvals$x
  # names(last_best2)<-names(scaling_fac)
  # # 
  # last_best<-last_best2
  
  priors <- list(
    mcstate::pmcmc_parameter("beta_1",
                             initial= last_best[["beta_1"]],#0.4*scaling_fac[["beta_1"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_1"]]),
    mcstate::pmcmc_parameter("beta_2",
                             initial= last_best[["beta_2"]],#0.1*scaling_fac[["beta_2"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_2"]]),
    
    mcstate::pmcmc_parameter("beta_3",
                             initial= last_best[["beta_3"]],#0.4*scaling_fac[["beta_1"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_3"]]),
    mcstate::pmcmc_parameter("beta_4",
                             initial= last_best[["beta_4"]],#0.1*scaling_fac[["beta_2"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_4"]]),
    
    mcstate::pmcmc_parameter("aduRR",
                             initial= last_best[["aduRR"]],#0.3*scaling_fac[["aduRR"]], # Excess infectiousness in Und 5
                             min = 0*scaling_fac[["aduRR"]], max = 1*scaling_fac[["aduRR"]]),
    
    mcstate::pmcmc_parameter("maternalAB",
                             initial= last_best[["maternalAB"]],#170*scaling_fac[["maternalAB"]], # Duration of maternal Abs
                             min = 1*scaling_fac[["maternalAB"]], max = 365*scaling_fac[["maternalAB"]] ),
    
    mcstate::pmcmc_parameter("imm_yr",
                             initial=last_best[["imm_yr"]],#2*scaling_fac[["imm_yr"]], # Duration of immunity
                             min = 0.5*scaling_fac[["imm_yr"]], max = 100*scaling_fac[["imm_yr"]] ),
    
    mcstate::pmcmc_parameter("imm_fac",
                             initial=last_best[["imm_fac"]],#2*scaling_fac[["imm_yr"]], # Duration of immunity
                             min = 1*scaling_fac[["imm_fac"]] ),
   
    
    mcstate::pmcmc_parameter("imm_fac",
                             initial=last_best[["imm_fac"]],#2*scaling_fac[["imm_yr"]], # Duration of immunity
                             min = 1*scaling_fac[["imm_fac"]] ),
    
    
    mcstate::pmcmc_parameter("repfac_0",
                             initial= last_best[["repfac_0"]],#287*scaling_fac[["repfac"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_0"]], max = 3000*scaling_fac[["repfac_0"]]),
    
    mcstate::pmcmc_parameter("repfac_5",
                             initial= last_best[["repfac_5"]],#287*scaling_fac[["repfac"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_5"]], max = 3000*scaling_fac[["repfac_5"]]),
    
    mcstate::pmcmc_parameter("repfac_15",
                             initial= last_best[["repfac_15"]],#287*scaling_fac[["repfac"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_15"]], max = 3000*scaling_fac[["repfac_15"]]),
    
    
    mcstate::pmcmc_parameter("repfac_65p",
                             initial= last_best[["repfac_65p"]],#287*scaling_fac[["repfac"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_65p"]], max = 3000*scaling_fac[["repfac_65p"]]),
    
    mcstate::pmcmc_parameter("crossp_GI",
                             initial= last_best[["crossp_GI"]],#0.05*scaling_fac[["crossp_GI"]], # Duration of maternal Abs
                             min = 0*scaling_fac[["crossp_GI"]],max = 1*scaling_fac[["crossp_GI"]]),
    
    mcstate::pmcmc_parameter("crossp_GII",
                             initial= last_best[["crossp_GII"]],#0.05*scaling_fac[["crossp_GI"]], # Duration of maternal Abs
                             min = 0*scaling_fac[["crossp_GII"]],max = 1*scaling_fac[["crossp_GII"]]))
  
  # Covariance matrix
  vcv <-cov(thetas)
  
  
}else {
  
  priors <- list(
    mcstate::pmcmc_parameter("beta_1",
                             initial=   0.06*scaling_fac[["beta_1"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_1"]]),
    mcstate::pmcmc_parameter("beta_2",
                             initial=  0.09*scaling_fac[["beta_2"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_2"]]),
    
    mcstate::pmcmc_parameter("beta_3",
                             initial=  0.18*scaling_fac[["beta_3"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_3"]]),
    mcstate::pmcmc_parameter("beta_4",
                             initial= 0.1*scaling_fac[["beta_4"]], # Transmission per capita
                             min = 0*scaling_fac[["beta_4"]]),
    
    mcstate::pmcmc_parameter("aduRR",
                             initial=  0.2*scaling_fac[["aduRR"]], # Excess infectiousness in Und 5
                             min = 0*scaling_fac[["aduRR"]], max = 1*scaling_fac[["aduRR"]]),
    
    mcstate::pmcmc_parameter("maternalAB",
                             initial=  249 *scaling_fac[["maternalAB"]], # Duration of maternal Abs
                             min = 1*scaling_fac[["maternalAB"]], max = 365*scaling_fac[["maternalAB"]] ),
    
    
    mcstate::pmcmc_parameter("imm_yr",
                             initial= 30*scaling_fac[["imm_yr"]], # Duration of immunity
                             min = 0.5*scaling_fac[["imm_yr"]], max = 100*scaling_fac[["imm_yr"]] ),
    
    mcstate::pmcmc_parameter("imm_fac",
                             initial=  3*scaling_fac[["imm_fac"]], # Duration of immunity
                             min = 1*scaling_fac[["imm_fac"]] ),
    
    
    mcstate::pmcmc_parameter("repfac_0",
                             initial= 1000*scaling_fac[["repfac_0"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_0"]], max = 3000*scaling_fac[["repfac_0"]]),
    
    mcstate::pmcmc_parameter("repfac_5",
                             initial= 670*scaling_fac[["repfac_5"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_5"]], max = 3000*scaling_fac[["repfac_5"]]),
    
    mcstate::pmcmc_parameter("repfac_15",
                             initial= 480*scaling_fac[["repfac_15"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_15"]], max = 3000*scaling_fac[["repfac_15"]]),
    
    mcstate::pmcmc_parameter("repfac_65p",
                             initial= 40*scaling_fac[["repfac_65p"]], # Seasonality amplitude
                             min = 1*scaling_fac[["repfac_65p"]], max = 3000*scaling_fac[["repfac_65p"]]),
    
    
    mcstate::pmcmc_parameter("crossp_GI",
                             initial= 0.2*scaling_fac[["crossp_GI"]], # Duration of maternal Abs
                             min = 0*scaling_fac[["crossp_GI"]],max = 1*scaling_fac[["crossp_GI"]]),
    
    
    mcstate::pmcmc_parameter("crossp_GII",
                             initial= 0.2*scaling_fac[["crossp_GII"]], # Duration of maternal Abs
                             min = 0*scaling_fac[["crossp_GII"]],max = 1*scaling_fac[["crossp_GII"]]))
  
  # Define covariance matrix at start
  ini<-c(priors[[1]]$mean,
         priors[[2]]$mean,
         priors[[3]]$mean,
         priors[[4]]$mean,
         priors[[5]]$mean,
         priors[[6]]$mean,
         priors[[7]]$mean,
         priors[[8]]$mean,
         priors[[9]]$mean,
         priors[[10]]$mean,
         priors[[11]]$mean,
         priors[[12]]$mean,
         priors[[13]]$mean,
         priors[[14]]$mean
  )*0.2
  
  vcv <- diag(ini,length(ini))
  
}




# Create pmcmc parameters
mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform = footransform)

if (stochastic==1){
  control <- mcstate::pmcmc_control(
    n_steps = n_steps,
    n_chains = 1,
    n_threads_total = 2,
    n_workers = 1,
    save_state = FALSE,
    save_trajectories = FALSE,
    progress = TRUE)
}else{
  control <- mcstate::pmcmc_control(
    n_steps = n_steps,
    n_chains = chainsn,
    n_threads_total = 16,
    n_workers = chainsn,
    save_state = FALSE,
    save_trajectories = FALSE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_burnin = n_burnin,
    n_steps_retain = n_out)
  
}
control_premulti<- mcstate::pmcmc_control(
  n_steps = 20,
  n_chains = 1,
  n_threads_total = 1,
  n_workers = 1,
  save_state = FALSE,
  save_trajectories = FALSE,
  progress = TRUE,
  adaptive_proposal = FALSE)

samples_temp <- mcstate::pmcmc(mcmc_pars, filter, control = control_premulti)
samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)


# Plot initial 
plot(samples$probabilities[, "log_posterior"], type = "l",
     xlab = "Sample", ylab = "Log posterior")


# Process samples
matplot(samples$pars,type = "l" , xlim = c(0,n_steps))


processed_chains <- samples#mcstate::pmcmc_thin(samples, 
#   burnin = n_burnin, 
#    thin = n_thin)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)

print(parameter_mean_hpd/unlist(scaling_fac))

plot(processed_chains$probabilities[, "log_posterior"], type = "l",
     xlab = "Sample", ylab = "Log posterior")



if (stochastic==1){
  save(processed_chains,file=here("output","processed_samples.RData"))
  
}else{
  
  save(processed_chains,file=here("output",paste0("chains_",scenario,".RData")))
  
  
}



