
compare_empty <- function(state, observed, pars = NULL) {
  
  llk<- -100
}


# ----------------------- Compare 
compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6
  
  llk_funs<-pars$llk_funs
  reported_t<-pars$reported_t
  
  noise<-rexp(n = length(state['inc_day_gii4_1',]), rate = exp_noise)
  
  

# IID2 data incidence -----------------------------------------------------

  
  modelled_irate <-rbind(
    ( (  state['inc_year_gi3_1',]+
                    state['inc_year_gi_1',]+
                    state['inc_year_gii4_1',]+
                    state['inc_year_gii_1',])/(state['pop_by4age1',])) + noise,

    ( (  state['inc_year_gi3_2',]+
                    state['inc_year_gi_2',]+
                    state['inc_year_gii4_2',]+
                    state['inc_year_gii_2',])/(state['pop_by4age2',])) + noise,

    ( (  state['inc_year_gi3_3',]+
                    state['inc_year_gi_3',]+
                    state['inc_year_gii4_3',]+
                    state['inc_year_gii_3',])/(state['pop_by4age3',])) + noise,

    ( (  state['inc_year_gi3_4',]+
                    state['inc_year_gi_4',]+
                    state['inc_year_gii4_4',]+
                    state['inc_year_gii_4',])/(state['pop_by4age4',])) + noise,

    ( (  state['inc_year_gi3_5',]+
                    state['inc_year_gi_5',]+
                    state['inc_year_gii4_5',]+
                    state['inc_year_gii_5',])/(state['pop_by4age5',])) + noise
    )*1000


  
  

  observed_size<-c(
    26.9,   # Person-years in 0 to 1
    190.8,  # Person-years in 1 to 4
    424.1,  # Person-years in 5 to 14
    2647.8, # Person-years in 15 to 65
    1369.1  # Person-years in 65plus
  )
  
  
  
  observations_irate <-rbind(
    observed$cases_a1,
    observed$cases_a2,
    observed$cases_a3,
    observed$cases_a4,
    observed$cases_a5
  )
  
  # llk_irate<-colSums(dbinom(x =  observations_irate,
  #                               size = round(observed_size),
  #                               prob =  modelled_irate,
  #                               log = TRUE),na.rm=TRUE)
  

  llk_irate<-colSums(dpois(x =  round(observations_irate),
                               lambda = modelled_irate,
                                log = TRUE),na.rm=TRUE)
  
  

  
  # SGSS weekly reported series all -------------------------------------------------------
  modelled_report<-
    state['reported_wk_1',] +
    state['reported_wk_2',] +
    state['reported_wk_3',] +
    state['reported_wk_4',] +
    state['reported_wk_5',] +
    state['reported_wk_6',] +
    state['reported_wk_7',] +
    state['reported_wk_8',] +
    state['reported_wk_9',] +
    state['reported_wk_10',] +
    state['reported_wk_11',] +
    state['reported_wk_12',] +
    state['reported_wk_13',] +
    state['reported_wk_14',] + noise
  
  observations_reported<-observed$reported
  
  
  negbin_dispersion<-pars$reported_var
  # # #
  # llk_reported<-dnbinom(x = observations_reported,
  #                       size = modelled_report,prob= negbin_dispersion , log = TRUE)/313

  llk_reported<-dpois(x = observations_reported,
                          lambda= modelled_report, log = TRUE)/313


  
  # SGSS weekly reported series by strain-------------------------------------------------------
  
  
  modelled_report_gi3<- state['reported_wk_gi3',] + noise
  modelled_report_gi<- state['reported_wk_gi',] + noise
  modelled_report_gii4<- state['reported_wk_gii4',] + noise
  modelled_report_gii<- state['reported_wk_gii',] + noise
  
  observations_reported_gi3<-observed$reported_gi3
  observations_reported_gi<-observed$reported_gi
  observations_reported_gii4<-observed$reported_gii4
  observations_reported_gii<-observed$reported_gii
  
  
  # negbin_dispersion<-pars$reported_var  
  # # # #  
  # llk_reported_gi3<-dnbinom(x = observations_reported_gi3,
  #                       size = modelled_report_gi3,prob= negbin_dispersion , log = TRUE)/52
  # 
  # llk_reported_gi<-dnbinom(x = observations_reported_gi,
  #                       size = modelled_report_gi,prob= negbin_dispersion , log = TRUE)/52
  # 
  # llk_reported_gii4<-dnbinom(x = observations_reported_gii4,
  #                       size = modelled_report_gii4,prob= negbin_dispersion , log = TRUE)/52
  # 
  # llk_reported_gii<-dnbinom(x = observations_reported_gii,
  #                       size = modelled_report_gii,prob= negbin_dispersion , log = TRUE)/52

  

  llk_reported_gi3<-dpois(x = observations_reported_gi3,
                            lambda= modelled_report_gi3, log = TRUE)/52

  llk_reported_gi<-dpois(x = observations_reported_gi,
                           lambda = modelled_report_gi, log = TRUE)/52

  llk_reported_gii4<-dpois(x = observations_reported_gii4,
                             lambda = modelled_report_gii4, log = TRUE)/52

  llk_reported_gii<-dpois(x = observations_reported_gii,
                            lambda = modelled_report_gii, log = TRUE)/52

  
  
  
  # Seroprev GII4 children 1 to 7 -------------------------------------------------------
  
  
  modelled_sero<-rbind(
    state['seroprev1.2',]  + noise ,
    state['seroprev2.3',]  + noise,
    state['seroprev3.4',]  + noise,
    state['seroprev4.5',]  + noise,
    state['seroprev5.6',]  + noise,
    state['seroprev6.7',] + noise)
  
  
  
  
  observed_event<-rbind(
    observed$sero1,
    observed$sero2,
    observed$sero3,
    observed$sero4,
    observed$sero5,
    observed$sero6
  )
  
  observed_size<-rbind(
    103,
    107,
    121,
    124,
    122,
    109
  )
  
  llk_sero<-colSums(dbinom(x=round(observed_event),
                           size = observed_size,
                           prob = modelled_sero,
                           log = TRUE),na.rm = TRUE)
  
  
  # SGSS age reporting fractions -------------------------------------------------------
  
  # # Weekly cases reported 0-4
  modelled_04<-rbind(
    ((state['cumm_incday_gii4_1',]+
        state['cumm_incday_gii_1',]+
        state['cumm_incday_gi3_1',]+
        state['cumm_incday_gi_1',]) /
       (    state['cumm_incday_gii4_1',]+
              state['cumm_incday_gii4_2',]+
              state['cumm_incday_gii4_3',]+
              state['cumm_incday_gii4_4',]+
              state['cumm_incday_gii_1',]+
              state['cumm_incday_gii_2',]+
              state['cumm_incday_gii_3',]+
              state['cumm_incday_gii_4',]+
              state['cumm_incday_gi3_1',]+
              state['cumm_incday_gi3_2',]+
              state['cumm_incday_gi3_3',]+
              state['cumm_incday_gi3_4',]+
              state['cumm_incday_gi_1',]+
              state['cumm_incday_gi_2',]+
              state['cumm_incday_gi_3',]+
              state['cumm_incday_gi_4',])) + noise)
  
  
  
  # 
  # # Weekly cases reported 5-14
  modelled_05<-rbind(
    ((state['cumm_incday_gii4_2',]+
        state['cumm_incday_gii_2',]+
        state['cumm_incday_gi3_2',]+
        state['cumm_incday_gi_2',]) /
       (    state['cumm_incday_gii4_1',]+
              state['cumm_incday_gii4_2',]+
              state['cumm_incday_gii4_3',]+
              state['cumm_incday_gii4_4',]+
              state['cumm_incday_gii_1',]+
              state['cumm_incday_gii_2',]+
              state['cumm_incday_gii_3',]+
              state['cumm_incday_gii_4',]+
              state['cumm_incday_gi3_1',]+
              state['cumm_incday_gi3_2',]+
              state['cumm_incday_gi3_3',]+
              state['cumm_incday_gi3_4',]+
              state['cumm_incday_gi_1',]+
              state['cumm_incday_gi_2',]+
              state['cumm_incday_gi_3',]+
              state['cumm_incday_gi_4',])) + noise)
  
  
  
  # # Weekly cases reported 65+
  modelled_15<-rbind(
    ((state['cumm_incday_gii4_3',]+
        state['cumm_incday_gii_3',]+
        state['cumm_incday_gi3_3',]+
        state['cumm_incday_gi_3',]) /
       (  state['cumm_incday_gii4_1',]+
            state['cumm_incday_gii4_2',]+
            state['cumm_incday_gii4_3',]+
            state['cumm_incday_gii4_4',]+
            state['cumm_incday_gii_1',]+
            state['cumm_incday_gii_2',]+
            state['cumm_incday_gii_3',]+
            state['cumm_incday_gii_4',]+
            state['cumm_incday_gi3_1',]+
            state['cumm_incday_gi3_2',]+
            state['cumm_incday_gi3_3',]+
            state['cumm_incday_gi3_4',]+
            state['cumm_incday_gi_1',]+
            state['cumm_incday_gi_2',]+
            state['cumm_incday_gi_3',]+
            state['cumm_incday_gi_4',])) + noise)
  
  
  
  # Weekly cases reported 65+
  modelled_65<-rbind(
    ((  state['cumm_incday_gii4_4',]+
          state['cumm_incday_gii_4',]+
          state['cumm_incday_gi3_4',]+
          state['cumm_incday_gi_4',]) /
       (    state['cumm_incday_gii4_1',]+
              state['cumm_incday_gii4_2',]+
              state['cumm_incday_gii4_3',]+
              state['cumm_incday_gii4_4',]+
              state['cumm_incday_gii_1',]+
              state['cumm_incday_gii_2',]+
              state['cumm_incday_gii_3',]+
              state['cumm_incday_gii_4',]+
              state['cumm_incday_gi3_1',]+
              state['cumm_incday_gi3_2',]+
              state['cumm_incday_gi3_3',]+
              state['cumm_incday_gi3_4',]+
              state['cumm_incday_gi_1',]+
              state['cumm_incday_gi_2',]+
              state['cumm_incday_gi_3',]+
              state['cumm_incday_gi_4',])) + noise)
  
  
  sample_size_fac<-12
  observed_event<-10152/sample_size_fac
  
  observed_size<-round(56308/sample_size_fac)
  
  
  
  llk_reported_04<-dbinom(x    =  round(observed$a0_event),
                          size = observed_size,
                          prob =  modelled_04,
                          log = TRUE)
  
  
  observed_event<-1612/sample_size_fac
  
  
  llk_reported_05<-dbinom(x    =  round(observed$a5_event),
                          size = observed_size,
                          prob =  modelled_05,
                          log = TRUE)
  
  
  observed_event<-11305/sample_size_fac
  
  
  llk_reported_15<-dbinom(x    = round(observed$a15_event),
                          size = observed_size,
                          prob =  modelled_15,
                          log = TRUE)
  
  
  observed_event<- 33239/sample_size_fac
  
  
  
  llk_reported_65<-dbinom(x    = round(observed$a65_event),
                          size = observed_size,
                          prob =  modelled_65,
                          log = TRUE)
  
  
  llk_reported_byage<-sum(c(llk_reported_04, llk_reported_05,
                            llk_reported_15,llk_reported_65)
                          , na.rm = TRUE)
  
  
  
  # Posterior -------------------------------------------------------
  posterior<-colSums(
    rbind(                  llk_irate,
                            llk_reported,
                            llk_reported_gi3,
                            llk_reported_gi,
                            llk_reported_gii4,
                            llk_reported_gii,
                            llk_reported_byage,
                            llk_sero
    )
    ,na.rm=T)
  


  return(posterior)
  
}
