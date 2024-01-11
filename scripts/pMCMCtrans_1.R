#############
## PMCMC
#############
# Controls 
n_steps  <- 1000
n_burnin <- round(n_steps*0.1)
n_out    <- round(n_steps*0.50)
n_thin   <- round((n_steps-n_burnin)/n_out)

## Transform parameters 
# Define minimum values & space range 

theta_min<-list(
  beta = 0.01,
  und5inf = 1,
  delta = 1,
  rho = 0.001,
  tau = 0.1,
  w1 = 0.01
)

theta_range<-list(
  beta = 1-0.01,
  und5inf = 10-1,
  delta = 365*2-1,
  rho = 1-0.001,
  tau = 25-0.1,
  w1 = 1-0.01
)

trans_par<-function(x,min,range,n){
  
  y <- n * (x-min)/range 
}
backtrans_par<-function(y,min,range,n){
  x <- (min*n + y*range)/n 
}


# Load previous chains
load(here("output","processed_samples.RData"))
id<- which(processed_chains$probabilities[,"log_posterior"] == 
             max(processed_chains$probabilities[,"log_posterior"]))[1]
# 
thetas<-processed_chains$pars
mean_hpd <- apply(thetas, 2, mean)
last_best <- thetas[id,]
# Define  priors 

priors <- list(
  mcstate::pmcmc_parameter("beta", 
                           initial= trans_par(0.05, 
                                              theta_min[["beta"]],
                                              theta_range[["beta"]],100), # Transmission per capita
                           min = 0),
  
  mcstate::pmcmc_parameter("und5inf",
                           initial= trans_par(5, 
                                              theta_min[["und5inf"]],
                                              theta_range[["und5inf"]],100), # Transmission per capita
                           min = 0),
  
  mcstate::pmcmc_parameter("delta", 
                           initial= trans_par(79, 
                                              theta_min[["delta"]],
                                              theta_range[["delta"]],100), # Transmission per capita
                           min = 0),
  
  mcstate::pmcmc_parameter("rho", 
                           initial= trans_par(0.05, 
                                              theta_min[["rho"]],
                                              theta_range[["rho"]],100), # Transmission per capita
                           min = 0),
  
  mcstate::pmcmc_parameter("tau", 
                           initial= trans_par(2, 
                                              theta_min[["tau"]],
                                              theta_range[["tau"]],100), # Transmission per capita
                           min = 0),
  
  mcstate::pmcmc_parameter("w1", 
                           initial= trans_par(0.15, 
                                              theta_min[["w1"]],
                                              theta_range[["w1"]],100), # Transmission per capita
                           min = 0)
)



# Create params transformation function
make_transform <- function(c_mat,
                           c_mat2, 
                           infa_id, 
                           mu, 
                           school,
                           n_school_steps,
                           n_age,
                           aging_mat,
                           ini,
                           pop,
                           theta_min,
                           theta_range,
                           backfoo) {
  
  function(theta) {
    
    theta_back<-c(
      beta = backfoo(theta[["beta"]],
                     theta_min[["beta"]], 
                     theta_range[["beta"]],100),
      und5inf = backfoo(theta[["und5inf"]],
                        theta_min[["und5inf"]],
                        theta_range[["und5inf"]],100),
      delta = backfoo(theta[["delta"]],
                      theta_min[["delta"]],
                      theta_range[["delta"]],100),
      rho = backfoo(theta[["rho"]],
                    theta_min[["rho"]],
                    theta_range[["rho"]],100),
      tau = backfoo(theta[["tau"]],
                    theta_min[["tau"]],
                    theta_range[["tau"]],100),
      w1 = backfoo(theta[["w1"]],
                   theta_min[["w1"]],
                   theta_range[["w1"]],100)
    )
    
    # print(theta)
    # print(theta_back)
    c_mat[infa_id,infa_id]<-c_mat[infa_id,infa_id]* theta_back[["und5inf"]]
    c_mat2[infa_id,infa_id]<-c_mat2[infa_id,infa_id]* theta_back[["und5inf"]]
    
    c(list(
      pop  = pop,
      init  = ini,
      mu    = mu,
      m     = c_mat,
      m_holi= c_mat2,
      school_step= as.double(school),
      n_school_steps=n_school_steps,
      N_age = n_age,
      aging_mat = aging_mat ),
      as.list(theta_back))
  }
}
# make_transform <- function(par,ini) {
#   function(theta) {
#     
#     c_mat<-par$transmission
#     c_mat2<-par$transmission_holi
#     
#     c_mat[par$infa_id,par$infa_id]<-c_mat[par$infa_id,par$infa_id]*theta[["und5inf"]]
#     c_mat2[par$infa_id,par$infa_id]<-c_mat2[par$infa_id,par$infa_id]*theta[["und5inf"]]
#     
#     list(
#       beta  = theta[["beta"]] ,   # transm coefficient
#       repfac= theta[["repfac"]],
#       delta = theta[["delta"]],
#       tau   = theta[["tau"]],
#       init  = ini,
#       mu    = par$mu,
#       m     = c_mat,
#       m_holi= c_mat2,
#       school= as.double(par$school_uk),
#       aging_mat= par$aging_mat, 
#       N_age = par$N_age)
#   }
# }

c1<-params$transmission
c2<-params$transmission_holi 
id<-params$infa_id 
mu<-params$mu 
school<-params$school_uk
n_school_steps<-params$n_school_steps
n_age<-params$N_age
aging_mat<-params$aging_mat
transform <- make_transform(c1,
                            c2, 
                            id, 
                            mu, 
                            school,
                            n_school_steps,
                            n_age,
                            aging_mat,
                            init,
                            pop=params$pop,
                            theta_min=theta_min,
                            theta_range=theta_min,
                            backfoo = backtrans_par )


# Define covariance matrix at start
ini<-c(priors[[1]]$mean,
       priors[[2]]$mean,
       priors[[3]]$mean,
       priors[[4]]$mean,
       priors[[5]]$mean,
       priors[[6]]$mean
)*0.15


vcv <- diag(ini, 6)
vcv<-cov(thetas) 

# Create pmcmc parameters
mcmc_pars <- mcstate::pmcmc_parameters$new(priors, vcv, transform)

# bb<-mcmc_pars$propose(ini,vcv)
# v<-backtrans_par(bb[1],theta_min[["beta"]],theta_range[["beta"]],100)
# v<-backtrans_par(bb[2],theta_min[["und5inf"]],theta_range[["und5inf"]],100)
# v<-backtrans_par(bb[3],theta_min[["delta"]],theta_range[["delta"]],100)
# v<-backtrans_par(0,theta_min[["rho"]],theta_range[["rho"]],100)
# v<-backtrans_par(0,theta_min[["tau"]],theta_range[["tau"]],100)
# v<-backtrans_par(0,theta_min[["w1"]],theta_range[["w1"]],100)

control <- mcstate::pmcmc_control(
  n_steps = n_steps,
  n_chains = 1,
  n_threads_total = 1,
  n_workers = 1,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE)


samples <- mcstate::pmcmc(mcmc_pars, filter, control = control)


# Plot initial 
plot(samples$probabilities[, "log_posterior"], type = "s",
     xlab = "Sample", ylab = "Log posterior",
     xlim = c(0,n_burnin))


# Process samples

processed_chains <- mcstate::pmcmc_thin(samples, burnin = n_burnin, thin = n_thin)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
parameter_mean_hpd

save(processed_chains,file=here("output","processed_samples.RData"))

mcmc1 <- coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars))
summary(mcmc1)
windows()
plot(mcmc1)


coda::effectiveSize(mcmc1)
1 - coda::rejectionRate(mcmc1)


