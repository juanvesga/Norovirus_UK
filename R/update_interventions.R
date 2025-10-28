# Process model output function -------------------------------------------
update_parameters<-function(ins,params){

  
  if("imm_fac" %in% names(ins)){
    
    params$factor_beta_1       = ins["factor_beta_1"] 
    params$factor_beta_2       = ins["factor_beta_2"] 
    params$log_beta_3       = ins["log_beta_3"] 
    params$factor_beta_4       = ins["factor_beta_4"] 
    params$log_maternalAB      = ins["log_maternalAB"]
    params$log_aduRR       = ins["log_aduRR"]
    params$imm_yr          = ins["imm_yr"]
    params$imm_fac         = ins["imm_fac"]
    params$log_ratio_0     = ins["log_ratio_0"]
    params$log_ratio_5     = ins["log_ratio_5"]
    params$log_ratio_15    = ins["log_ratio_15"]
    params$log_repfac_65p  = ins["log_repfac_65p"]
    params$geno_frac       = ins["geno_frac"]

    
  } else {
    params$factor_beta_1   = ins["factor_beta_1"] 
    params$factor_beta_2   = ins["factor_beta_2"] 
    params$log_beta_3      = ins["log_beta_3"] 
    params$factor_beta_4   = ins["factor_beta_4"] 
    params$log_maternalAB      = ins["log_maternalAB"]
    params$log_aduRR       = ins["log_aduRR"]
    params$imm_yr          = ins["imm_yr"]
    params$log_ratio_0     = ins["log_ratio_0"]
    params$log_ratio_5     = ins["log_ratio_5"]
    params$log_ratio_15    = ins["log_ratio_15"]
    params$log_repfac_65p  = ins["log_repfac_65p"]
    params$geno_frac       = ins["geno_frac"]
  }
  
  
  
  return(params)
}


update_interventions<-function(mod,start_t,state0,run_t,param,label){
  
  # Set time   
  dust2::dust_system_set_time(mod, start_t)
  
  # Set model state  
  dust2::dust_system_set_state(mod,state0)
  
  # Set parameters
  dust2::dust_system_update_pars(mod, pars=param)
  
  # run simulation
  y0<- dust2::dust_system_simulate(mod, seq(start_t,start_t+run_t))
  
  # unpack results
  y0 <- dust2::dust_unpack_state(mod, y0)
  
  
  

  # Create output dataframe

  out<-data.frame(
    days                   = seq(0,run_t),  
    new_yearly_hosp        = y0$new_hosp_elder_yr+y0$new_hosp_adult_yr,
    new_weekly_cases       = y0$new_cases_week,
    new_weekly_deaths     =  y0$death_reported_wk,
    new_daily_cases        = y0$new_cases,
    cum_new_daily_cases    = cumsum(y0$new_cases),
    cum_new_weekly_cases   = cumsum(y0$new_cases_week),
    vacc_doses             = y0$n_vacc_doses,
    cum_vacc_doses         = cumsum(y0$n_vacc_doses),
    new_weekly_cases_gi3    = y0$new_cases_week_gi3,
    new_weekly_cases_gi     = y0$new_cases_week_gi,
    new_weekly_cases_gii4   = y0$new_cases_week_gii4,
    new_weekly_cases_gii    = y0$new_cases_week_gii,
    new_weekly_reported     = colSums(y0$reported_wk)
    
  )
  
  out$scenario<-label
  
  return(out)
  
}