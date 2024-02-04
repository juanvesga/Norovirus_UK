rm(list = ls(all.names = TRUE))
gc()


# Load Packages -----------------------------------------------------------

library(odin.dust)
library(here)
library(socialmixr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggmatplot)
library(profvis)
library(lubridate)
library(mcstate)


# Set outputs path

output_dir<- paste(here(),"_output", sep="")

# Set up model parameters -------------------------------------------------
actions<-c("MCMC","RUN_SIMS")
#actions<-"MCMC"
#actions<-"MCMC_DX"
#actions<- "RUN_SIMS"
#actions<-"PLOT_FITS"
#actions<- "DYNAMICS"
#actions<- "DIC"

# Set model scenario ------------------------------------------------------
scenarios<- c("imm1par","imm2par","imm4par","immdrop","immdrop_noreinfection","imm2par_noreinfection")
#scenarios<-"imm1par"
#scenarios<-"imm2par"
#scenarios<-"imm4par"
#scenarios<- "imm2par_noreinfection"
#scenarios<-"immdrop"
#scenarios<-"immdrop_noreinfection"
#scenarios<-"full"


# MCMC parameters 
nsamples<-500 # from the posterior
n_steps  <- 50000
n_out <- 5000
start_previous<-FALSE
burnin_prop<-0.5
chain_selection<-c(1,2,3,4)

for (aa in 1:length(actions)){
  
  action<-actions[aa]
  
  for (bb in 1:length(scenarios)){
    
    scenario<-scenarios[bb]
    
    # Source general files ------------------------------------------------------------
    source(here("src","utility_functions.R"))
    source(here("scripts","load_data.R"))
    source(here("src","model_functions.R"))
    source(here("src","run_simulations.R"))
    
    
    if(scenario == "imm4par"){
      
      source(here("src","make_transform_imm4par.R"))
      source(here("scripts","setup.model_imm4par.R"))
      
      
    }else {
      
      source(here("src","make_transform_imm1par.R"))
      source(here("scripts","setup.model.R"))
      
      
    }
    
    
    # Create model objects ----------------------------------------------------
    
    
    model_path<-here("src",paste0("seiar.age.",scenario,".R"))
    seiar <- odin.dust::odin_dust(model_path)
    
    #seiar_odin<-odin::odin(model_path)
    
    stochastic <- 0
    if (stochastic==1){
      filter <- mcstate::particle_filter$new(data_all, 
                                             model = seiar, 
                                             n_particles = 11,
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
    
    
    # MCMC --------------------------------------------------------------------
    
    if (action=="MCMC"){
      
      
      
      # Controls 
      start_from_best <-start_previous
      chainsn  <- 4
      n_steps  <- n_steps
      n_burnin <- round(n_steps*burnin_prop)
      n_out    <-  n_out
      n_thin   <- round((n_steps-n_burnin)/n_out)
      
      if(scenario == "imm4par"){
        
        source(here("scripts","pMCMC_imm4par.R"))
        
      }else if (scenario == "imm2par_noreinfection"){
        
        
        source(here("scripts","pMCMC_imm2par_noreinfection.R"))
        
        
      }else if (scenario == "immdrop_noreinfection"){
        
        
        source(here("scripts","pMCMC_immdrop_noreinfection.R"))
        
      }else{
        
        
        source(here("scripts","pMCMC_1.R"))
        
        
      }
      
      rm(processed_chains)
      # MCMC Dx -----------------------------------------------------------------
      
    } else if(action=="MCMC_DX"){
      
      
      load(paste0(output_dir,"/chains_",scenario,".RData")) 
      
      
      source(here("src","parameter_plot.R"))
      
      figs<-parameter_plot(processed_chains,nsamples,scalefc,scenario)
      
      
      gridExtra::grid.arrange(figs$plot1)
      
      gridExtra::grid.arrange(figs$immplot)
      
      gridExtra::grid.arrange(figs$immbox)
      
   
      
      #   
      # mcmc1 <- coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars))
      # summary(mcmc1)
      # windows()
      # plot(mcmc1)
      # 
      # 
      # 
      # coda::effectiveSize(mcmc1)
      # 1 - coda::rejectionRate(mcmc1)
      
      
      
      # Run Sims ----------------------------------------------------------------
      
      
    }else if(action=="RUN_SIMS"){
      
      
      load(paste0(output_dir,"/chains_",scenario,".RData"))
      
      
    

      
      id<-which(processed_chains$chain ==chain_selection)
      par<-processed_chains$pars[id,]
      M <- par[sample(nrow(par), nsamples), ]
      
    
      out<-run_simulations(nsamples,M,filter2,footransform)
      
      save(out,file=paste0(output_dir,"/runs_mcmc_",scenario,".Rdata"))
      
      
      rm(processed_chains)
      rm(out)
      
      # Plot MCMC fits ----------------------------------------------------------
      
      
    } else if(action=="PLOT_FITS"){
      
      
      load(paste0(output_dir,"/runs_mcmc_",scenario,".RData")) 
      
      source(here("scripts","plot_model_fits.R"))
      
      
      
      # Plot Dynamics -----------------------------------------------------------
      
    }else if(action=="DYNAMICS"){
      
      source(here("scripts","plot_Rt.R"))
   
      
      source(here("src","parameter_plot.R"))
      browser()
      source(here("scripts","plot_dynamics.R"))
     
      # Run DIC selections script -----------------------------------------------
    }else if(action=="DIC"){
     
      
      source(here("scripts","compareR0.R"))
       
      source(here("scripts","DICselection.R"))
      
    
      
       
    }
  }
}


