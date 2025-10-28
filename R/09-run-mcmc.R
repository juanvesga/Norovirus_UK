rm(list = ls())
gc()  # garbage collection
# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters_short.qs2")
infile_dataS   <- file.path(root,"output", "data_short.qs2")
infile_dataL   <- file.path(root,"output", "data_long.qs2")
infile_input   <- file.path(root,"output", "params_list_short.qs2")
infile_input2   <- file.path(root,"output", "params_list2_short.qs2")
infile_prior   <- file.path(root,"output", "priors.qs2")
infile_limits  <- file.path(root, "output", "limits.qs2")
infile_prior2  <- file.path(root, "output", "priors2.qs2")
infile_limits2 <- file.path(root, "output", "limits2.qs2")
infile_model0  <- file.path(root,"models", "model_simple_fits_monty.R")
infile_model1  <- file.path(root,"models", "model_simple_no_reinf.R")
infile_model2  <- file.path(root,"models", "model_2pars_fits_monty.R")
infile_model3  <- file.path(root,"models", "model_drop_fits_monty.R")


#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))


# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(odin2)
library(dust2)
library(monty)
library(posterior)
library(bayesplot)
library(lattice)

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------


## Cross protection 5%

#model       <- "_0" # Simple SEIAR
#model       <- "_1" # Simple SEIAR no reinf no cross-prot
model       <- "_2" # Full version with reinf and cross-protection
#model        <- "_3" # Full

## Cross protection 25%

#model      <- "_4" # Simple SEIAR
#model      <- "_5" # Simple SEIAR no reinf no cross-prot
#model      <- "_6" # Full version with reinf and cross-protection
#model      <- "_7" # Full version with reinf and cross-protection and drop immunity


## Cross protection 50%

#model       <- "_8" # Simple SEIAR
#model      <- "_9" # Simple SEIAR no reinf no cross-prot
#model      <- "_10" # Full version with reinf and cross-protection
#model      <- "_11" # Full version with reinf and cross-protection and drop immunity

outfile     <- file.path(root, "output", paste0("monty_chains",model,".rds"))
outfile2    <- file.path(root, "output", paste0("sample_coda",model,".qs2"))
outfile3    <- file.path(root, "output", paste0("sample_df",model,".qs2"))
infile_imis <- file.path(root, "output", paste0("imis",model,".qs2"))

# MCMC settings -----------------------------------------------------------
# MCMC conditions
nsamples    = 2
0000
use_imis    = FALSE # uses imis output as starting values
use_last    = TRUE # uses last sampling object
use_mle     = FALSE # uses last MLE run


nout    = nsamples 
nchains = 3 
to_burn = nsamples*0
each_n = round((nsamples-to_burn)/nout) 
set.seed(42)


# Source Odin model
if (model=="_0"||model=="_4"||model=="_8"){
  
  source(infile_model0)
  
  prior      <- qs_read(infile_prior2)
  limits     <- qs_read(infile_limits2)
  
  if (model=="_0"){
    crossp <-0.05
  }else if (model=="_4"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
  parameters <- qs_read(infile0)
  data       <- qs_read(infile_dataS)
  pars_list  <- qs_read(infile_input2)
  pars_list$pars_list$crossp_GI=crossp
  pars_list$pars_list$crossp_GII=crossp
  pars<-pars_list$pars_list
  packer<-pars_list$packer
  start_list<-pars_list$start_pars
  
 
  
}else if((model=="_1"||model=="_5"||model=="_9")){
  
  source(infile_model1)
  
  prior      <- qs_read(infile_prior2)
  limits     <- qs_read(infile_limits2)
  
  if (model=="_1"){
    crossp <-0.05
  }else if (model=="_5"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
  parameters <- qs_read(infile0)
  data       <- qs_read(infile_dataS)
  pars_list  <- qs_read(infile_input2)
  pars_list$pars_list$crossp_GI=crossp
  pars_list$pars_list$crossp_GII=crossp
  pars<-pars_list$pars_list
  packer<-pars_list$packer
  start_list<-pars_list$start_pars
  
}else if((model=="_2"||model=="_6"||model=="_10")){
  source(infile_model2)
  
  prior      <- qs_read(infile_prior)
  limits     <- qs_read(infile_limits)
  
  if (model=="_2"){
    crossp <-0.05
  }else if (model=="_6"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
  parameters <- qs_read(infile0)
  data       <- qs_read(infile_dataS)
  pars_list  <- qs_read(infile_input)
  pars_list$pars_list$crossp_GI=crossp
  pars_list$pars_list$crossp_GII=crossp
  pars<-pars_list$pars_list
  packer<-pars_list$packer
  start_list<-pars_list$start_pars

  
}else if((model=="_3"||model=="_7"||model=="_11")){
  source(infile_model3)
  
  prior      <- qs_read(infile_prior2)
  limits     <- qs_read(infile_limits2)
  
  if (model=="_3"){
    crossp <-0.05
  }else if (model=="_7"){
    crossp <-0.25
  } else {
    crossp <-0.5
  }
  
  parameters <- qs_read(infile0)
  data       <- qs_read(infile_dataS)
  pars_list  <- qs_read(infile_input2)
  pars_list$pars_list$crossp_GI=crossp
  pars_list$pars_list$crossp_GII=crossp
  pars<-pars_list$pars_list
  packer<-pars_list$packer
  start_list<-pars_list$start_pars
  
}


if (use_imis){
 

  
  imis<-qs_read(infile_imis)
  
  params<-imis$best_parameters
  
  
  fac<-0.1
  
  vcv<-diag(abs(params)*fac)

  # starting_points <- imis_out$mcmc_initialization$starting_points
  # proposal_cov    <- imis_out$mcmc_initialization$proposal_cov
  # 
  # rm(imis_out)
  # 
  # params<-starting_points[c(1:3)]
  # 
  # epsilon = 1e-6
  # 
  # d = nrow(proposal_cov)
  # 
  # # Create the stability matrix: epsilon * I
  # stability_matrix = epsilon * diag(d)
  # 
  # vcv <- (2.38^2 / d) * proposal_cov +   stability_matrix
  




  
}else if (use_last){

  tmp                 <- readRDS(outfile) 
  chains              <- coda::as.mcmc.list(tmp)
  longchain           <- do.call(rbind, chains)
  colnames(longchain) <- colnames(chains[[1]])
  covar               <- cov(longchain)
  density             <- matrix(as.vector(tmp$density), ncol = 1)

  mle                 <- which(density==max(density))[1]
  p_mle               <- longchain[mle,]
  
  if (use_mle){
  
  params<-p_mle # Use MLE as starting values

  }else{

  params<-tmp  # Start MCMC from exact point in last run

  }
  
  rm(tmp)
  
  
  epsilon = 1e-6
  
  d = nrow(covar)
  
  # Create the stability matrix: epsilon * I
  stability_matrix = epsilon * diag(d)
  
  vcv <- (2.38^2 / d) * covar +   stability_matrix
  
} else {
  
  if((model=="_0"||model=="_1"||model=="_4"||model=="_5"||model=="_8"||model=="_9")){
    
    # start of model immunity drop                                
    params<-c(
      factor_beta_1   =  0.0005947734 ,          
      factor_beta_2   =    0.0012281536,         
      log_beta_3      =   -1.5688842005 ,          
      factor_beta_4   =   0.3352941290   ,       
      log_maternalAB  =   5.3486274395  ,      
      log_aduRR       =   -0.0370342793 ,  
      imm_yr          =   15.6276032587 ,    
      log_ratio_0   =  1.4720076837 , 
      log_ratio_5   =   3.7545344692 , 
      log_ratio_15   =  2.8019324301,  
      log_repfac_65p=   log(150)) 
    
    

    rng <- monty::monty_rng_create()
    theta <- monty::monty_model_direct_sample(prior, rng)
    params<-theta
    
    
  }else if((model=="_2"||model=="_6"||model=="_10")){
    
    params<-c(
      factor_beta_1   =  0.200000000 ,              
      factor_beta_2   =   0.200000000 ,        
      log_beta_3      =   -2.059477382 ,          
      factor_beta_4   =   0.700000000   ,       
      log_maternalAB  =   4.122434293    , 
      log_aduRR       =    -0.001791197  ,  
      imm_yr          = 31.614787219  , 
      imm_fac        =  1.033024271   , 
      log_ratio_0   =  1.540533057  ,  
      log_ratio_5   =  3.547943527  ,  
      log_ratio_15   = 3.129645174  , 
      log_repfac_65p=   5.029942097 ,
      geno_frac       = 0.2) 
    
    # rng <- monty::monty_rng_create()
    # theta <- monty::monty_model_direct_sample(prior, rng)
    # params<-theta
    
  } else if((model=="_3"||model=="_7"||model=="_11")){
    
    # start of model immunity drop                                
    params<-c(
      factor_beta_1   =  0.0005947734 ,          
      factor_beta_2   =    0.0012281536,         
      log_beta_3      =   -1.5688842005 ,          
      factor_beta_4   =   0.3352941290   ,       
      log_maternalAB  =   5.3486274395  ,      
      log_aduRR       =   -0.0370342793 ,  
      imm_yr          =   15.6276032587 ,    
      log_ratio_0   =  1.4720076837 , 
      log_ratio_5   =   3.7545344692 , 
      log_ratio_15   =  2.8019324301,  
      log_repfac_65p=   3.9120230054) 
    
    # rng <- monty::monty_rng_create()
    # theta <- monty::monty_model_direct_sample(prior, rng)
    # params<-theta
    # 
    
    
  }
  

 
    fac<-0.25
  #browser()
  vcv<-diag(abs(params))*fac
  #vcv<-diag(length(params))*fac

  
 

}


# Create adaptive sampler
sampler <- monty::monty_sampler_adaptive(vcv)

  # MCMC function  
  run_mcmc <- function(sampler_obj,
                       nsamples, 
                       pars0,
                       nchains,
                       to_burn,
                       each_n) {
    
    # Create deterministic filter
    # Recreate all objects inside the function
    filter_det <- dust_unfilter_create(noro_model(), 
                                       time_start = 0, 
                                       data = data)
    
    likelihood <- dust_likelihood_monty(filter_det, 
                                        packer,
                                        save_trajectories = FALSE)
    
    

    posterior <- likelihood + prior
    
    # Use callr runner
    runner <- monty::monty_runner_callr(3)
    
    samples <- monty::monty_sample(model = posterior, 
                                   sampler = sampler_obj, 
                                   n_steps = nsamples,
                                   initial = pars0,
                                   n_chains = nchains, 
                                   runner = runner,
                                   burnin = to_burn,
                                   thinning_factor = each_n,
                                   restartable = TRUE)
    return(samples)
  }
  
  

  

  result<- run_mcmc(sampler,
                    nsamples, 
                    params,
                    nchains,
                    to_burn,
                    each_n)
  


# Save
#qs_save(result, outfile)


saveRDS(result, outfile)



# Quick exploration
samples_coda <- coda::as.mcmc.list(result)
samples_df <- posterior::as_draws_df(result)

qs_save(samples_coda, outfile2)
qs_save(samples_df, outfile3)



matplot(result$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")


plot1<-xyplot(samples_coda,layout = c(4, 4))
gridExtra::grid.arrange(plot1)

bayesplot::mcmc_trace(samples_df)

