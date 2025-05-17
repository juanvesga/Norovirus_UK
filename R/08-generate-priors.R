# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile_input  <- file.path(root,"output", "params_list.qs2")
outfile <- file.path(root, "output", "priors.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(mcstate, include.only = c("pmcmc_parameter"))

# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------

pars_list  <- qs_read(infile_input)

last_best<-pars_list$start_pars

## Priors
prior <- list(
  pmcmc_parameter("beta_1", last_best[["beta_1"]], min = log(0),prior = function(p) {
    dgamma(exp(p), shape = 2, scale = 0.15/2, log = TRUE) }), # Transmission per capita
  
  pmcmc_parameter("beta_2", last_best[["beta_2"]],min = log(0), prior = function(p) {
    dgamma(exp(p), shape = 2, scale = 0.15/2, log = TRUE) }), # Transmission per capita
  
  pmcmc_parameter("beta_3", last_best[["beta_3"]], min = log(0), prior = function(p) {
    dgamma(exp(p), shape = 3, scale = 0.2/3, log = TRUE) }), # Transmission per capita
  
  pmcmc_parameter("beta_4", last_best[["beta_4"]],min = log(0), prior = function(p) {
    dgamma(exp(p), shape = 3, scale = 0.2/3, log = TRUE) }), # Transmission per capita
  
  pmcmc_parameter("aduRR", last_best[["aduRR"]], min = log(0), max = log(1)),# prior = function(p) {
  #  dbeta(exp(p), shape1=3, shape2 =  3, log = TRUE) }),             # Excess infectiousness in Und 5
  
  pmcmc_parameter("maternalAB", last_best[["maternalAB"]], min = log(0), max = log(1000)),# min = log(0),prior = function(p) {
    # dgamma(exp(p), shape = 3, scale = 150/3, log = TRUE) }), 
  
  pmcmc_parameter("imm_yr", last_best[["imm_yr"]], min = log(0), prior = function(p) {
    dgamma(exp(p), shape = 4, scale = 20/4, log = TRUE) }),
  
  pmcmc_parameter("imm_fac", last_best[["imm_fac"]], min = log(0), max = log(5)), # Imm factor infection>1
  
  pmcmc_parameter("repfac_0",last_best[["repfac_0"]],min = log(0),  max = log(900)), #prior = function(p) {
  #  dgamma(exp(p), shape = 6, scale = 128/6, log = TRUE) }),
  
  pmcmc_parameter("repfac_5",last_best[["repfac_5"]], min = log(0),  max = log(1500)), #prior = function(p) {
  #  dgamma(exp(p), shape = 12, scale = 523/12, log = TRUE) }),
  
  pmcmc_parameter("repfac_15",last_best[["repfac_15"]],min = log(0), max = log(900)), #prior = function(p) {
  #  dgamma(exp(p), shape = 12, scale = 500/12, log = TRUE) }),
  
  pmcmc_parameter("repfac_65p",last_best[["repfac_65p"]],min = log(0), max = log(900)), #prior = function(p) {
  #  dgamma(exp(p), shape = 6, scale = 98/6, log = TRUE) }),
  
  pmcmc_parameter("crossp_GI", last_best[["crossp_GI"]],min = log(0), max = log(1)), #prior = function(p) {
  #  dbeta(exp(p), shape1=2, shape2 =  6, log = TRUE) }), 
  
  pmcmc_parameter("crossp_GII", last_best[["crossp_GII"]], min = log(0), max = log(1)), #prior = function(p) {
  #  dbeta(exp(p), shape1=2, shape2 =  6, log = TRUE) }), 
  
  pmcmc_parameter("reported_var", last_best[["reported_var"]], min = log(0), max = log(1))
)



# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(prior, outfile)
