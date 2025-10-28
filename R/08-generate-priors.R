# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root    <- here::here()
outfile <- file.path(root, "output", "priors.qs2")
outfile2 <- file.path(root, "output", "limits.qs2")

outfile3 <- file.path(root, "output", "priors2.qs2")
outfile4 <- file.path(root, "output", "limits2.qs2")

# Packages
source(file.path(root, "R", "modify_attach.R"))
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(monty, include.only = c("monty_dsl"))

# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------

# mean
last_best=
  c(
    log_beta_1      =     0.13410824,
    log_beta_2      =     0.12708480,
    log_beta_3      =     0.20300807,
    log_beta_4      =     0.27545582,
    aduRR       =     0.58903565,
    maternalAB  =     367.16148967,
    imm_yr      =     16.15967444 ,
    imm_fac     =     1.70697752,
    repfac_0    =     64.26091774, 
    repfac_5    =     450.23810289,
    repfac_15   =     345.14095794,  
    repfac_65p  =    59.89184305#,
  #  reported_var=      0.5,
    #season_lag  =      30,
    #season_amp  =      0.15
  )


prior <- monty_dsl({
  factor_beta_1    ~ Uniform(0.15,1)#Beta(2,2)
  factor_beta_2    ~ Uniform(0.15,1)#Beta(2,2)
  log_beta_3       ~ Normal(log(0.18),0.5)
  factor_beta_4    ~ Uniform(0.15,1)#Beta(2,2)
  log_maternalAB  ~ Uniform(log(1),log(365))#  ~ Gamma(shape = 4, scale = 150/3)
  log_aduRR       ~ Uniform(log(0.01),log(1))
  imm_yr          ~ Gamma(shape = 4, scale = 10/4)
  imm_fac         ~ Uniform(1,5)
  log_ratio_0     ~ Normal(log(2.7), 0.5)        # repfac_0/repfac_65p ≈ 136/51 ≈ 2.7
  log_ratio_5     ~ Normal(log(28), 0.5)         # repfac_5/repfac_65p ≈ 1444/51 ≈ 28  
  log_ratio_15    ~ Normal(log(2.9), 0.5)        # repfac_15/repfac_65p ≈ 150/51 ≈ 2.9
  log_repfac_65p  ~ Normal(log(90), 0.3)         # reference: lowest underre0porting
  geno_frac       ~ Beta(3, 17)

})

prior2 <- monty_dsl({
  factor_beta_1    ~ Uniform(0.15,1)#Beta(2,2)
  factor_beta_2    ~ Uniform(0.15,1)#Beta(2,2)
  log_beta_3       ~ Normal(log(0.18),0.5)
  factor_beta_4    ~ Uniform(0.15,1)#Beta(2,2)
  log_maternalAB  ~ Uniform(log(1),log(365))
  log_aduRR       ~ Uniform(log(0.01),log(1))
  imm_yr          ~ Gamma(shape = 4, scale = 10/4)
  log_ratio_0     ~ Normal(log(2.7), 0.5)        # repfac_0/repfac_65p ≈ 136/51 ≈ 2.7
  log_ratio_5     ~ Normal(log(28), 0.5)         # repfac_5/repfac_65p ≈ 1444/51 ≈ 28  
  log_ratio_15    ~ Normal(log(2.9), 0.5)        # repfac_15/repfac_65p ≈ 150/51 ≈ 2.9
  log_repfac_65p  ~ Normal(log(90), 0.3)         # reference: lowest underre0porting
  geno_frac       ~ Beta(3, 17)
  # season_lag      ~ Uniform(0,120)
  # season_amp      ~ Uniform(0,1)
})


# prior2 <- monty_dsl({
#   beta_1 ~     Uniform(0,0.25)#Gamma(shape = 2, scale = 0.15/2)
#   beta_2 ~     Uniform(0,0.25)#Gamma(shape = 2, scale = 0.15/2)
#   beta_3 ~     Uniform(0,0.25)#Gamma(shape = 3, scale = 0.2/3)
#   beta_4 ~     Uniform(0,0.25)#Gamma(shape = 3, scale = 0.2/3)
#   aduRR  ~     Uniform(0,1)
#   maternalAB ~ Uniform(1,365*2)#Gamma(shape = 3, scale = 150/3)
#   imm_yr~      Uniform(0,100)#Gamma(shape = 4, scale = 20/4)
#   repfac_0 ~   Uniform(0,900)#Gamma(shape = 6, scale = 128/6)
#   repfac_5 ~   Uniform(0,1500)#Gamma(shape = 12, scale = 523/12)
#   repfac_15~   Uniform(0,900)#Gamma(shape = 12, scale = 500/12)
#   repfac_65p ~ Uniform(0,900)#Gamma(shape = 6, scale = 98/6)
#   season_lag  ~ Uniform(-10+20,15+20)
#   season_amp  ~ Uniform(0,1)
# })



# mean
limits<-list(
lower=
  (c(
    log_beta_1      = log(0.01),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
    log_beta_2      = log(0.01),        # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
    log_beta_3      = log(0.01),        # Normal(log(0.155), 0.5) ≈ Normal(-1.86, 0.5)
    log_beta_4      = log(0.01),        # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
    maternalAB      = 1,              # Gamma(4, 50) covers ~30-250 days
    log_aduRR       = log(0.01),
    imm_yr          = 5,
    imm_fac         = 1,
    log_ratio_0     = log(1),
    log_ratio_5     = log(1),
    log_ratio_15    = log(1),
    log_repfac_65p  = log(1)
    # log_season_lag      = -10+20,
    # log_season_amp      =      0        # Normal(log(0.09), 0.3) ≈ Normal(-2.41, 0.3)
    
  )),
upper=
  (c(
    log_beta_1      =  log(0.5),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
    log_beta_2      =  log(0.5),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
    log_beta_3      =  log(0.5),       # Normal(log(0.155), 0.5) ≈ Normal(-1.86, 0.5)
    log_beta_4     =  log(0.5),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
    maternalAB      =  700,          # Gamma(4, 50) covers ~30-250 days
    log_aduRR       =  log(1),       # Normal(log(0.5), 0.3) ≈ Normal(-0.69, 0.3)
    imm_yr          =  20,            # Gamma(4, 2.5) covers ~5-20 years
    imm_fac         =  5,             # Uniform(1, 5)
    log_ratio_0     =  log(10),         # Normal(log(2.7), 0.5) ≈ Normal(0.99, 0.5)
    log_ratio_5     =  log(40),         # Normal(log(28), 0.5) ≈ Normal(3.33, 0.5)
    log_ratio_15    =  log(10),         # Normal(log(2.9), 0.5) ≈ Normal(1.06, 0.5)
    log_repfac_65p  =  log(600)         # Normal(log(51), 0.3) ≈ Normal(3.93, 0.3)
    # log_season_lag      = 15+20,         # Normal(log(11), 0.3) ≈ Normal(2.40, 0.3)
    # log_season_amp      =  1        # Normal(log(0.09), 0.3) ≈ Normal(-2.41, 0.3)
  ))
)


limits2<-list(
  lower=
    (c(
      log_beta_1      = log(0.01),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
      log_beta_2      = log(0.01),        # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
      log_beta_3      = log(0.01),        # Normal(log(0.155), 0.5) ≈ Normal(-1.86, 0.5)
      log_beta_4      = log(0.01),        # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
      maternalAB      = 1,              # Gamma(4, 50) covers ~30-250 days
      log_aduRR       = log(0.01),
      imm_yr          = 5,
      log_ratio_0     = log(1),
      log_ratio_5     = log(1),
      log_ratio_15    = log(1),
      log_repfac_65p  = log(1)
      # log_season_lag      = -10+20,
      # log_season_amp      =      0        # Normal(log(0.09), 0.3) ≈ Normal(-2.41, 0.3)
      # 
      
    )),
  upper=
    (c(
      log_beta_1      =  log(0.5),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
      log_beta_2      =  log(0.5),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
      log_beta_3      =  log(0.5),       # Normal(log(0.155), 0.5) ≈ Normal(-1.86, 0.5)
      log_beta_4     =  log(0.5),       # Normal(log(0.048), 0.5) ≈ Normal(-3.04, 0.5)
      maternalAB      =  700,          # Gamma(4, 50) covers ~30-250 days
      log_aduRR       =  log(1),       # Normal(log(0.5), 0.3) ≈ Normal(-0.69, 0.3)
      imm_yr          =  20,            # Gamma(4, 2.5) covers ~5-20 years
      log_ratio_0     =  log(10),         # Normal(log(2.7), 0.5) ≈ Normal(0.99, 0.5)
      log_ratio_5     =  log(40) ,        # Normal(log(28), 0.5) ≈ Normal(3.33, 0.5)
      log_ratio_15    =  log(10),         # Normal(log(2.9), 0.5) ≈ Normal(1.06, 0.5)
      log_repfac_65p  =  log(600)         # Normal(log(51), 0.3) ≈ Normal(3.93, 0.3)
      # log_season_lag      = 15+20,         # Normal(log(11), 0.3) ≈ Normal(2.40, 0.3)
      # log_season_amp      =  1        # Normal(log(0.09), 0.3) ≈ Normal(-2.41, 0.3)
    ))
)


# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(prior, outfile)
qs_save(limits, outfile2)
qs_save(prior2, outfile3)
qs_save(limits2, outfile4)
