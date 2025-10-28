rm(list = ls())
gc()  # garbage collection
# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")

outfile1        <- file.path(root, "output", "mcmc_fits.qs2")
outfile2        <- file.path(root, "output", "mcmc_inits.qs2")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))
source(file.path(root, "R", "collect_function.R"))
# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(abind, include.only = c("abind"))

library(ggplot2)
library(dplyr)
library(lattice)  ## for the 'xyplot' command
library(fitR)
library(bayesplot)
library(coda)
library(viridisLite)

# -------------------------------------------------------------------------
#  -------------------------------------------------------------------
# -------------------------------------------------------------------------

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

burn_in    <- 0
thin_in    <- 10


infile_chains  <- file.path(root, "output", paste0("monty_chains_weight",model,".rds"))
result         <- readRDS(infile_chains)

samples_coda <- coda::as.mcmc.list(result)
samples_df <- posterior::as_draws_df(result)


# 
# infile_chains2 <- file.path(root, "output", paste0("sample_coda",model,".qs2"))
# infile_chains3 <- file.path(root, "output", paste0("sample_df",model,".qs2"))
# 
# samples_coda   <- qs_read(infile_chains2)
# samples_df     <- qs_read(infile_chains2)



windows()
bp<-bayesplot::mcmc_trace(samples_df)
gridExtra::grid.arrange(bp)

summary<-posterior::summarise_draws(samples_df)

print(summary)

summary$mean
# browser()
# 
# d<- ncol(trace[[1]])
# trace[[1]][,1:(d-1)] <-(trace[[1]][,1:(d-1)])
# trace[[2]][,1:(d-1)] <-(trace[[2]][,1:(d-1)])
# trace[[3]][,1:(d-1)] <-(trace[[3]][,1:(d-1)])
# trace[[4]][,1:(d-1)] <-(trace[[4]][,1:(d-1)])


# burn and thin
traceBurn <- burnAndThin(samples_coda, burn = burn_in)
chains <- burnAndThin(traceBurn, thin = thin_in)



plot1<-xyplot(chains,layout = c(4, 4))
windows()
gridExtra::grid.arrange(plot1)


# removing the burn-in increases the ESS
print(coda::effectiveSize(chains))


##           R_0         D_lat         D_inf         alpha         D_imm 
##      2495.223      3024.481      3132.444      3201.725      3188.532 
##           rho      logPrior logLikelihood    logDensity 
##      3350.921         0.000      2681.293      2681.293

# autocorrelation


plot3<-coda::acfplot(chains, lag.max = nrow(chains[[1]]))
windows()
gridExtra::grid.arrange(plot3)



# Note that plotPosteriorDensity can take a list of mcmc.list It will plot the
# different mcmc.list by combining their elements Let's plot the combined
# unthinned trace vs the combined thinned trace.
plot4<-plotPosteriorDensity(list(full =traceBurn , thinned_burn = chains))
windows()
gridExtra::grid.arrange(plot4)

plot5<-levelplot(chains[[1]], col.regions = heat.colors(100))
windows()
gridExtra::grid.arrange(plot5)


windows()
mcmc_pairs(samples_df,pars=c("log_beta_3","factor_beta_1","factor_beta_2","factor_beta_4"),
           off_diag_args = list(size = 0.75))



thin_burned<- seq(max(burn_in,1),nrow(result$density),by=max(thin_in,1))

windows()

n_chains <- ncol(result$density[thin_burned,])
matplot(result$density[thin_burned,], type = "l", lty = 1,
        col = viridis(n_chains, alpha = 0.7),
        xlab = "Sample", ylab = "Log posterior probability density")


