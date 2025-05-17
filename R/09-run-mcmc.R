# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()
infile0        <- file.path(root,"output", "parameters.qs2")
infile_dataS   <- file.path(root,"output", "data_short.qs2")
infile_dataL   <- file.path(root,"output", "data_long.qs2")
infile_input   <- file.path(root,"output", "params_list.qs2")
infile_prior   <- file.path(root,"output", "priors.qs2")
infile_model2  <- file.path(root,"models", "model_2pars.R")
infile_compare <- file.path(root,"models", "compare_model.R") 

outfile <- file.path(root, "output", "chains.qs")
outfile_index <- file.path(root, "output", "index.qs")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "index_function.R"))
source(infile_compare)

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))

library(odin)
library(dust)
library(mcstate)
library(posterior)
library(bayesplot)

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------
# MCMC conditions
nsamples  = 100
n_out     = 100 
nchains  = 3 
to_burn  = nsamples*0


set.seed(42)

# Source Odin model
parameters <- qs_read(infile0)
data       <- qs_read(infile_dataS)
pars_list  <- qs_read(infile_input)
prior      <- qs_read(infile_prior)

#Parameters object
pars<-pars_list$pars_list

start_list<-pars_list$start_pars

footransform<-pars_list$footransform

vcv <- diag(abs(start_list),length(start_list))*0.01


# Create pmcmc parameters
mcmc_pars <- mcstate::pmcmc_parameters$new(prior, vcv, transform = footransform)

control <- mcstate::pmcmc_control(
  n_steps = nsamples,
  n_chains = nchains,
  n_threads_total = nchains*4,
  n_workers = nchains,
  save_state = FALSE,
  save_trajectories = FALSE,
  progress = TRUE,
  adaptive_proposal = TRUE,
  n_burnin = to_burn,
  n_steps_retain = n_out)

control_premulti<- mcstate::pmcmc_control(
  n_steps = 1,
  n_chains = 1,
  n_threads_total = 1,
  n_workers = 1,
  save_state = FALSE,
  save_trajectories = FALSE,
  progress = TRUE,
  adaptive_proposal = FALSE)

# call model object
seiar <- odin.dust::odin_dust(infile_model2)

data_compare<-mcstate::particle_filter_data(data, tstep, 1, 0)



filter <- mcstate::particle_deterministic$new(data_compare,
                                              model = seiar,
                                              compare = compare,
                                              index = index,
                                              initial = basic_initial)





#samples_temp <- mcstate::pmcmc(mcmc_pars, filter, control = control_premulti)


samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)


chain1<-which(samples$chain==1)
chain2<-which(samples$chain==2)
chain3<-which(samples$chain==3)

samples_df1 <- posterior::as_draws_df(exp(samples$pars[chain1,]))
samples_df2 <- posterior::as_draws_df(exp(samples$pars[chain2,]))
samples_df3 <- posterior::as_draws_df(exp(samples$pars[chain3,]))


samples_array  <- posterior::bind_draws(
  samples_df1, 
  samples_df2,
  samples_df3,
  along = "chain")

print(samples_array)
samples_df <- posterior::as_draws_array(exp(samples$pars), .nchains=3)

summary<-posterior::summarise_draws(samples_df)


matplot(samples$probabilities[, "log_posterior"] , type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")




plot(samples$probabilities[chain1, "log_posterior"],
     type="l",xlab = "iteration",ylab = "log-Posterior", col="grey")
lines(samples$probabilities[chain2, "log_posterior"], col="firebrick")
lines(samples$probabilities[chain3, "log_posterior"], col="navy")



root<-here::here()

bayesplot::mcmc_trace(samples_array)


bayesplot::mcmc_dens_overlay(samples_array)

bayesplot::mcmc_hist(samples_array)

bayesplot::mcmc_violin(samples_array)

bayesplot::mcmc_intervals(samples_array)

bayesplot::mcmc_acf(samples_array)




# Update covariance matrix ------------------------------------------------


vcv <-cov(samples$pars)


# Create pmcmc parameters
mcmc_pars <- mcstate::pmcmc_parameters$new(prior, vcv, transform = footransform)


samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)


chain1<-which(samples$chain==1)
chain2<-which(samples$chain==2)
chain3<-which(samples$chain==3)

samples_df1 <- posterior::as_draws_df(exp(samples$pars[chain1,]))
samples_df2 <- posterior::as_draws_df(exp(samples$pars[chain2,]))
samples_df3 <- posterior::as_draws_df(exp(samples$pars[chain3,]))


samples_array  <- posterior::bind_draws(
  samples_df1, 
  samples_df2,
  samples_df3,
  along = "chain")

print(samples_array)
samples_df <- posterior::as_draws_array(exp(samples$pars), .nchains=3)

summary<-posterior::summarise_draws(samples_df)


matplot(samples$probabilities[, "log_posterior"] , type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")




plot(samples$probabilities[chain1, "log_posterior"],
     type="l",xlab = "iteration",ylab = "log-Posterior", col="grey")
lines(samples$probabilities[chain2, "log_posterior"], col="firebrick")
lines(samples$probabilities[chain3, "log_posterior"], col="navy")



root<-here::here()

bayesplot::mcmc_trace(samples_array)


bayesplot::mcmc_dens_overlay(samples_array)

bayesplot::mcmc_hist(samples_array)

bayesplot::mcmc_violin(samples_array)

bayesplot::mcmc_intervals(samples_array)

bayesplot::mcmc_acf(samples_array)




qs_save(samples, outfile)
qs_save(index,outfile_index)


