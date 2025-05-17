
# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
infile1  <- file.path(root,"output", "parameters.qs2")
infile2  <- file.path(root,"output", "polymod.qs2")
infile3  <- file.path(root,"output", "comix.qs2")
infile4  <- file.path(root,"output", "school_uk.qs2")
infile5  <- file.path(root,"output", "covid_sche.qs2")

outfile <- file.path(root, "output", "params_list.qs2")
#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))


# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------

parameters <- qs_read(infile1)
cmat       <- qs_read(infile2)
comix      <- qs_read(infile3)
school_uk  <- qs_read(infile4)
covid_sche <- qs_read(infile5)

# Create parameters transform object
make_transform <- function(
    c_mat,
    c_mat2, 
    cmx_1,
    cmx_2,
    cmx_3,
    cmx_4,
    cmx_5,
    cmx_6,
    cmx_7,
    cmx_8,
    cmx_9,
    mu, 
    school,
    n_school_steps,
    covid_step,
    n_covid_steps,
    n_age,
    aging_vec,
    pop,
    dt) {
  
  function(theta) {
    
    theta_back<-exp(c(
      beta_1    = theta[["beta_1"]],
      beta_2    = theta[["beta_2"]],
      beta_3    = theta[["beta_3"]],
      beta_4    = theta[["beta_4"]],
      maternalAB = theta[["maternalAB"]],
      aduRR      = theta[["aduRR"]],
      imm_yr     = theta[["imm_yr"]], 
      imm_fac    = theta[["imm_fac"]], 
      repfac_0   = theta[["repfac_0"]],
      repfac_5   = theta[["repfac_5"]],
      repfac_15  = theta[["repfac_15"]],
      repfac_65p = theta[["repfac_65p"]],
      crossp_GI  = theta[["crossp_GI"]],
      crossp_GII = theta[["crossp_GII"]],
      reported_var=theta[["reported_var"]]
    ))
  
    c(list(
      pop  = pop,
      mu    = mu,
      m     = c_mat,
      m_holi= c_mat2,
      cmx_1=cmx_1,
      cmx_2=cmx_2,
      cmx_3=cmx_3,
      cmx_4=cmx_4,
      cmx_5=cmx_5,
      cmx_6=cmx_6,
      cmx_7=cmx_7,
      cmx_8=cmx_8,
      cmx_9=cmx_9,
      covid_step= as.double(covid_step),
      n_covid_steps= n_covid_steps,
      school_step= as.double(school),
      n_school_steps=n_school_steps,
      N_age = n_age,
      aging_vec = aging_vec,
      dt=dt
      ),
      as.list(theta_back))
  }
}

# Create parameters objects

footransform <-  make_transform(
  c_mat  = cmat$transmission,
  c_mat2 = cmat$transmission_holi, 
  cmx_1  = comix$cmx_1,
  cmx_2  = comix$cmx_2,
  cmx_3  = comix$cmx_3,
  cmx_4  = comix$cmx_4,
  cmx_5  = comix$cmx_5,
  cmx_6  = comix$cmx_6,
  cmx_7  = comix$cmx_7,
  cmx_8  = comix$cmx_8,
  cmx_9  = comix$cmx_9,
   mu    = parameters$mort_rates/(365),
  school = school_uk,
  n_school_steps = length(school_uk),
  covid_step     = covid_sche,
  n_covid_steps  = length(covid_sche),
  n_age          = length(parameters$ages),
  aging_vec      = c(1/head(diff(c(0,parameters$ages)),-1),0)/(365),
  pop            = parameters$pop,
  dt             = 1 )



start_pars<- log(c(
  beta_1       = 0.13953720,
  beta_2       = 0.12680114,
  beta_3       = 0.20410557,
  beta_4       = 0.27306245,
  aduRR        = 0.69372532  ,
  maternalAB   = 520.61197156,
  imm_yr       = 17.18686356,
  imm_fac      = 1.68360017 ,
  repfac_0     = 105.07189709,
  repfac_5     = 954.17803464 ,
  repfac_15    = 581.94243972 ,
  repfac_65p   = 95.41216406  ,
  crossp_GI    = 0.04103662   ,
  crossp_GII   = 0.01982703,
  reported_var =  0.35656712))


# Create filter
input_object<-list(
start_pars=start_pars,
footransform=footransform)

# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(input_object, outfile)



