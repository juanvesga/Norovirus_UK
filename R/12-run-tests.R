# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile0  <- file.path(root,"output", "parameters.qs2")
infile_dataS  <- file.path(root,"output", "data_short.qs2")
infile_dataL  <- file.path(root,"output", "data_long.qs2")
infile_input  <- file.path(root,"output", "params_list.qs2")
infile_prior  <- file.path(root,"output", "priors.qs2")
infile_model2  <- file.path(root,"models", "model_2pars.R")

outfile <- file.path(root, "output", "chains.qs")
outfile_index <- file.path(root, "output", "index.qs")

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

# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------


# Source Odin model
source(infile_model2)
parameters <- qs_read(infile0)
data       <- qs_read(infile_dataS)
pars_list  <- qs_read(infile_input)
prior      <- qs_read(infile_prior)
#Parameters object
pars<-pars_list$pars_list
packer<-pars_list$packer
start_list<-pars_list$start_pars


#dust object
sys  <- dust_system_create(noro_model(), pars, n_particles = 1, deterministic = TRUE)
index<- dust_unpack_index(sys)
state<-dust_system_state(sys)
named_state<-dust_unpack_state(sys,state)

# Generate Initial conditions vector
seed<-1# Number infected per strain
#init<-generate_initial_state(index,state,parameters,seed,2)


# check initial conditions in current object
dust_system_state(sys)

# Pass initial conditions to model
dust_system_set_state_initial(sys)

# Check again initial conditions in object
dust_system_state(sys)

#Simulate system up to a point 
t <- seq(0, 20000)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)


matplot(t(y$pop_all),type='l')

plot(y$pop_all[8,])

# plot(t, y$new_cases_week, type = "l", xlab = "Time", ylab = "Infected population",
#      xlim = c(1e4,2e4),ylim=c(0,5e4))
# 
# plot(t[t %% 7 == 0], y$new_cases_week[t %% 7 == 0], type = "o", pch = 19,
#      ylab = "Infection incidence", xlab = "Time",xlim = c(1e4,10600),ylim=c(0,15e4))
# lines(t, y$new_cases_week, col = "red")