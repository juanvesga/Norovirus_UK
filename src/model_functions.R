

# ----------------------- Compare 
compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6
  
  noise<-rexp(n = length(state['inc_year_gii4_1',]), rate = exp_noise)

  # Incidence rates per 1000 py
  modelled_irate <-rbind(
    1000*( (state['inc_year_gi3_1',]+
              state['inc_year_gi_1',]+
              state['inc_year_gii4_1',]+
              state['inc_year_gii_1',])/(state['pop_by4age1',])) + noise,
    
    1000*((state['inc_year_gi3_2',]+
             state['inc_year_gi_2',]+
             state['inc_year_gii4_2',]+
             state['inc_year_gii_2',])/(state['pop_by4age2',])) + noise,
    
    1000*((state['inc_year_gi3_3',]+
             state['inc_year_gi_3',]+
             state['inc_year_gii4_3',]+
             state['inc_year_gii_3',])/(state['pop_by4age3',])) + noise,
    
    1000*((state['inc_year_gi3_4',]+
             state['inc_year_gi_4',]+
             state['inc_year_gii4_4',]+
             state['inc_year_gii_4',])/(state['pop_by4age4',])) + noise, 
    
    1000*((state['inc_year_gi3_5',]+
             state['inc_year_gi_5',]+
             state['inc_year_gii4_5',]+
             state['inc_year_gii_5',])/(state['pop_by4age5',])) + noise)
  

  observations_irate <-rbind(
    observed$cases_a1,
    observed$cases_a2,
    observed$cases_a3,
    observed$cases_a4,
    observed$cases_a5
  ) 
  llk_irate<-colSums(dpois(x = observations_irate, lambda = modelled_irate, log = TRUE),na.rm=TRUE)
  
  
  # Strain distribution by Age in IID2 incidence - Harris
  
  denom_1<-(state['inc_year_gi3_1',]+
              state['inc_year_gi_1',]+
              state['inc_year_gii4_1',]+
              state['inc_year_gii_1',])+
    (state['inc_year_gi3_2',]+
       state['inc_year_gi_2',]+
       state['inc_year_gii4_2',]+
       state['inc_year_gii_2',])
  
  denom_2<-(state['inc_year_gi3_3',]+
              state['inc_year_gi_3',]+
              state['inc_year_gii4_3',]+
              state['inc_year_gii_3',])
  
  denom_3<-(state['inc_year_gi3_4',]+
              state['inc_year_gi_4',]+
              state['inc_year_gii4_4',]+
              state['inc_year_gii_4',])
  
  denom_4<-(state['inc_year_gi3_5',]+
              state['inc_year_gi_5',]+
              state['inc_year_gii4_5',]+
              state['inc_year_gii_5',])
    
  modelled_gi3_1 <- ((state['inc_year_gi3_1',]+ state['inc_year_gi3_2',])/denom_1) +noise
  modelled_gi_1  <- ((state['inc_year_gi_1',]+ state['inc_year_gi_2',])/denom_1) +noise
  modelled_gii4_1<- ((state['inc_year_gii4_1',]+ state['inc_year_gii4_2',])/denom_1) +noise
  modelled_gii_1 <- ((state['inc_year_gii_1',]+ state['inc_year_gii_2',])/denom_1) +noise
  
  modelled_gi3_2 <- ((state['inc_year_gi3_3',])/denom_2) +noise
  modelled_gi_2  <- ((state['inc_year_gi_3',])/denom_2) +noise
  modelled_gii4_2<- ((state['inc_year_gii4_3',])/denom_2) +noise
  modelled_gii_2 <- ((state['inc_year_gii_3',])/denom_2) +noise
  
  modelled_gi3_3 <- ((state['inc_year_gi3_4',])/denom_3) +noise
  modelled_gi_3  <- ((state['inc_year_gi_4',])/denom_3) +noise
  modelled_gii4_3<- ((state['inc_year_gii4_4',])/denom_3) +noise
  modelled_gii_3 <- ((state['inc_year_gii_4',])/denom_3) +noise
  
  modelled_gi3_4 <- ((state['inc_year_gi3_5',])/denom_4) +noise
  modelled_gi_4  <- ((state['inc_year_gi_5',])/denom_4) +noise
  modelled_gii4_4<- ((state['inc_year_gii4_5',])/denom_4) +noise
  modelled_gii_4 <- ((state['inc_year_gii_5',])/denom_4) +noise
  
  observed_size_1<-68
  observed_size_2<-41
  observed_size_3<-153
  observed_size_4<-48
  
  llk_gi3_1 <-dbinom(x=round(observed$gi3_prop_1*observed_size_1), size = observed_size_1, prob = modelled_gi3_1, log = TRUE )
  llk_gi_1  <-dbinom(x=round(observed$gi_prop_1*observed_size_1),  size = observed_size_1, prob = modelled_gi_1, log = TRUE )
  llk_gii4_1<-dbinom(x=round(observed$gii4_prop_1*observed_size_1),size = observed_size_1, prob = modelled_gii4_1, log = TRUE )
  llk_gii_1 <-dbinom(x=round(observed$gii_prop_1*observed_size_1), size = observed_size_1, prob = modelled_gii_1, log = TRUE )

  llk_gi3_2 <-dbinom(x=round(observed$gi3_prop_2*observed_size_2), size = observed_size_2, prob = modelled_gi3_2, log = TRUE )
  llk_gi_2  <-dbinom(x=round(observed$gi_prop_2*observed_size_2),  size = observed_size_2, prob = modelled_gi_2, log = TRUE )
  llk_gii4_2<-dbinom(x=round(observed$gii4_prop_2*observed_size_2),size = observed_size_2, prob = modelled_gii4_2, log = TRUE )
  llk_gii_2 <-dbinom(x=round(observed$gii_prop_2*observed_size_2), size = observed_size_2, prob = modelled_gii_2, log = TRUE )

  llk_gi3_3 <-dbinom(x=round(observed$gi3_prop_3*observed_size_3), size = observed_size_3, prob = modelled_gi3_3, log = TRUE )
  llk_gi_3  <-dbinom(x=round(observed$gi_prop_3*observed_size_3),  size = observed_size_3, prob = modelled_gi_3, log = TRUE )
  llk_gii4_3<-dbinom(x=round(observed$gii4_prop_3*observed_size_3),size = observed_size_3, prob = modelled_gii4_3, log = TRUE )
  llk_gii_3 <-dbinom(x=round(observed$gii_prop_3*observed_size_3), size = observed_size_3, prob = modelled_gii_3, log = TRUE )

  llk_gi3_4 <-dbinom(x=round(observed$gi3_prop_4*observed_size_4), size = observed_size_4, prob = modelled_gi3_4, log = TRUE )
  llk_gi_4  <-dbinom(x=round(observed$gi_prop_4*observed_size_4),  size = observed_size_4, prob = modelled_gi_4, log = TRUE )
  llk_gii4_4<-dbinom(x=round(observed$gii4_prop_4*observed_size_4),size = observed_size_4, prob = modelled_gii4_4, log = TRUE )
  llk_gii_4 <-dbinom(x=round(observed$gii_prop_4*observed_size_4), size = observed_size_4, prob = modelled_gii_4, log = TRUE )
  

  llk_genotype_byage<-sum(c(llk_gi3_1, llk_gi3_2, llk_gi3_3, llk_gi3_4,
                            llk_gi_1,  llk_gi_2,  llk_gi_3,  llk_gi_4,
                            llk_gii4_1,llk_gii4_2,llk_gii4_3,llk_gii4_4,
                            llk_gii_1, llk_gii_2, llk_gii_3, llk_gii_4), na.rm = TRUE)

  # # Strain prevalence GI.3 in SGSS
  # modelled_gi3<-rbind(
  #   ((state['cumm_incday_gi3_1',]+
  #       state['cumm_incday_gi3_2',]+
  #       state['cumm_incday_gi3_3',]+
  #       state['cumm_incday_gi3_4',]) /
  #      
  #      (    state['cumm_incday_gii4_1',]+
  #           state['cumm_incday_gii4_2',]+
  #           state['cumm_incday_gii4_3',]+
  #           state['cumm_incday_gii4_4',]+
  #           state['cumm_incday_gii_1',]+
  #           state['cumm_incday_gii_2',]+
  #           state['cumm_incday_gii_3',]+
  #           state['cumm_incday_gii_4',]+
  #           state['cumm_incday_gi3_1',]+
  #           state['cumm_incday_gi3_2',]+
  #           state['cumm_incday_gi3_3',]+
  #           state['cumm_incday_gi3_4',]+
  #           state['cumm_incday_gi_1',]+
  #           state['cumm_incday_gi_2',]+
  #           state['cumm_incday_gi_3',]+
  #           state['cumm_incday_gi_4',])
  #    ) + noise)
  # 
  # 
  # observed_event<-164
  # 
  # observed_size<-2710
  # 
  # 
  # llk_gi3<-dbinom(x=round( observed$gi3_prop*observed_size),
  #                 size = observed_size,
  #                 prob =  modelled_gi3,
  #                 log = TRUE)
  # 
  # 
  # 
  # # Strain prevalence GI inSGSS
  # modelled_gi<-rbind(
  #   ((  state['cumm_incday_gi_1',]+
  #         state['cumm_incday_gi_2',]+
  #         state['cumm_incday_gi_3',]+
  #         state['cumm_incday_gi_4',]) /
  #      (    
  #           state['cumm_incday_gii4_1',]+
  #           state['cumm_incday_gii4_2',]+
  #           state['cumm_incday_gii4_3',]+
  #           state['cumm_incday_gii4_4',]+
  #           state['cumm_incday_gii_1',]+
  #           state['cumm_incday_gii_2',]+
  #           state['cumm_incday_gii_3',]+
  #           state['cumm_incday_gii_4',]+
  #           state['cumm_incday_gi3_1',]+
  #           state['cumm_incday_gi3_2',]+
  #           state['cumm_incday_gi3_3',]+
  #           state['cumm_incday_gi3_4',]+
  #           state['cumm_incday_gi_1',]+
  #           state['cumm_incday_gi_2',]+
  #           state['cumm_incday_gi_3',]+
  #           state['cumm_incday_gi_4',])) + noise)
  # 
  # observed_event<-210
  # 
  # observed_size<-2710
  # 
  # llk_gi<-dbinom(x=round( observed$gi_prop*observed_size),
  #                size = observed_size,
  #                prob =  modelled_gi,
  #                log = TRUE)
  # 
  # 
  # # Strain prevalence GII.4 in SGSS
  # modelled_gii4<-rbind(
  #   ((state['cumm_incday_gii4_1',]+
  #       state['cumm_incday_gii4_2',]+
  #       state['cumm_incday_gii4_3',]+
  #       state['cumm_incday_gii4_4',]) /
  #      (  
  #           state['cumm_incday_gii4_1',]+
  #           state['cumm_incday_gii4_2',]+
  #           state['cumm_incday_gii4_3',]+
  #           state['cumm_incday_gii4_4',]+
  #           state['cumm_incday_gii_1',]+
  #           state['cumm_incday_gii_2',]+
  #           state['cumm_incday_gii_3',]+
  #           state['cumm_incday_gii_4',]+
  #           state['cumm_incday_gi3_1',]+
  #           state['cumm_incday_gi3_2',]+
  #           state['cumm_incday_gi3_3',]+
  #           state['cumm_incday_gi3_4',]+
  #           state['cumm_incday_gi_1',]+
  #           state['cumm_incday_gi_2',]+
  #           state['cumm_incday_gi_3',]+
  #           state['cumm_incday_gi_4',])) + noise)
  # 
  # 
  # observed_event<-1532
  # 
  # observed_size<-2710
  # 
  # llk_gii4<-dbinom(x=round( observed$gii4_prop*observed_size),
  #                size = observed_size,
  #                prob =  modelled_gii4,
  #                log = TRUE)
  # 
  # 
  # # Strain prevalence GII SGSS
  # modelled_gii<-rbind(
  #   ((state['cumm_incday_gii_1',]+
  #       state['cumm_incday_gii_2',]+
  #       state['cumm_incday_gii_3',]+
  #       state['cumm_incday_gii_4',]) /
  #      (   
  #             state['cumm_incday_gii4_1',]+
  #             state['cumm_incday_gii4_2',]+
  #             state['cumm_incday_gii4_3',]+
  #             state['cumm_incday_gii4_4',]+
  #             state['cumm_incday_gii_1',]+
  #             state['cumm_incday_gii_2',]+
  #             state['cumm_incday_gii_3',]+
  #             state['cumm_incday_gii_4',]+
  #             state['cumm_incday_gi3_1',]+
  #             state['cumm_incday_gi3_2',]+
  #             state['cumm_incday_gi3_3',]+
  #             state['cumm_incday_gi3_4',]+
  #             state['cumm_incday_gi_1',]+
  #             state['cumm_incday_gi_2',]+
  #             state['cumm_incday_gi_3',]+
  #             state['cumm_incday_gi_4',])) + noise)
  # 
  # 
  # observed_event<-210
  # 
  # observed_size<-2710
  # 
  # llk_gii<-dbinom(x=round( observed$gii_prop*observed_size),
  #                  size = observed_size,
  #                  prob =  modelled_gii,
  #                  log = TRUE)
  
  
  
  ##### By age
  
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
  
  
  observed_event<-10152
  
  observed_size<-56308
  llk_reported_04<-dbinom(x    = round(observed$a0_prop*observed_size),
                          size = observed_size,
                          prob =  modelled_04,
                          log = TRUE)
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
  
  
  observed_event<-1612
  
  observed_size<-56308
  llk_reported_05<-dbinom(x    = round(observed$a5_prop*observed_size),
                          size = observed_size,
                          prob =  modelled_05,
                          log = TRUE)
  # 
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
  
  
  observed_event<-11305
  
  observed_size<-56308
  llk_reported_15<-dbinom(x    = round(observed$a15_prop*observed_size),
                          size = observed_size,
                          prob =  modelled_15,
                          log = TRUE)
  # 
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
  
  
  observed_event<- 33239
  
  observed_size<-56308
  
  llk_reported_65<-dbinom(x    = round(observed$a65_prop*observed_size),
                          size = observed_size,
                          prob =  modelled_65,
                          log = TRUE)
  
  
  ## Reported weekly
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
  llk_reported<-dpois(x = observations_reported, 
                      lambda = modelled_report, log = TRUE)/496
  
  # Seroprevalence in children 1 to 7
  # modelled_sero<-rbind(
  #   state['seroprev_num1.2',]/state['seroprev_den1.2',] + noise,
  #   state['seroprev_num2.3',]/state['seroprev_den2.3',] + noise,
  #   state['seroprev_num3.4',]/state['seroprev_den3.4',] + noise,
  #   state['seroprev_num4.5',]/state['seroprev_den4.5',] + noise,
  #   state['seroprev_num5.6',]/state['seroprev_den5.6',] + noise,
  #   state['seroprev_num6.7',]/state['seroprev_den6.7',] + noise
  # )
  modelled_sero<-rbind(
    state['seroprev1.2',]  + noise ,
    state['seroprev2.3',]  + noise,
    state['seroprev3.4',]  + noise,
    state['seroprev4.5',]  + noise,
    state['seroprev5.6',]  + noise,
    state['seroprev6.7',] + noise)
  
  
  observed_event<-rbind(
    observed$sero1*103,
    observed$sero2*107,
    observed$sero3*121,
    observed$sero4*124,
    observed$sero5*122,
    observed$sero6*109
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
  
  
  posterior<-colSums(rbind(llk_irate,
                           # llk_gi3,
                           # llk_gi,
                           # llk_gii4,
                           # llk_gii, 
                           llk_genotype_byage,
                           llk_reported_04,
                           llk_reported_05,
                           llk_reported_15,
                           llk_reported_65,
                           llk_reported,
                           llk_sero
  ),na.rm=T)

  # print(posterior)
 
  return(posterior)
  
}

index <- function(info) {
  list(run = c(inc_year_gii4_1 = info$index$inc_year_gii4[1],
               inc_year_gii4_2 = info$index$inc_year_gii4[2],
               inc_year_gii4_3 = info$index$inc_year_gii4[3],
               inc_year_gii4_4 = info$index$inc_year_gii4[4],
               inc_year_gii4_5 = info$index$inc_year_gii4[5],
               
               inc_year_gii_1 = info$index$inc_year_gii[1],
               inc_year_gii_2 = info$index$inc_year_gii[2],
               inc_year_gii_3 = info$index$inc_year_gii[3],
               inc_year_gii_4 = info$index$inc_year_gii[4],
               inc_year_gii_5 = info$index$inc_year_gii[5],
               
               inc_year_gi3_1 = info$index$inc_year_gi3[1],
               inc_year_gi3_2 = info$index$inc_year_gi3[2],
               inc_year_gi3_3 = info$index$inc_year_gi3[3],
               inc_year_gi3_4 = info$index$inc_year_gi3[4],
               inc_year_gi3_5 = info$index$inc_year_gi3[5],
               
               inc_year_gi_1 = info$index$inc_year_gi[1],
               inc_year_gi_2 = info$index$inc_year_gi[2],
               inc_year_gi_3 = info$index$inc_year_gi[3],
               inc_year_gi_4 = info$index$inc_year_gi[4],
               inc_year_gi_5 = info$index$inc_year_gi[5],
               
               inc_day_gii4_1 = info$index$inc_day_gii4[1],
               inc_day_gii4_2 = info$index$inc_day_gii4[2],
               inc_day_gii4_3 = info$index$inc_day_gii4[3],
               inc_day_gii4_4 = info$index$inc_day_gii4[4],
               inc_day_gii4_5 = info$index$inc_day_gii4[5],
               
               inc_day_gii_1 = info$index$inc_day_gii[1],
               inc_day_gii_2 = info$index$inc_day_gii[2],
               inc_day_gii_3 = info$index$inc_day_gii[3],
               inc_day_gii_4 = info$index$inc_day_gii[4],
               inc_day_gii_5 = info$index$inc_day_gii[5],
               
               inc_day_gi3_1 = info$index$inc_day_gi3[1],
               inc_day_gi3_2 = info$index$inc_day_gi3[2],
               inc_day_gi3_3 = info$index$inc_day_gi3[3],
               inc_day_gi3_4 = info$index$inc_day_gi3[4],
               inc_day_gi3_5 = info$index$inc_day_gi3[5],
               
               inc_day_gi_1 = info$index$inc_day_gi[1],
               inc_day_gi_2 = info$index$inc_day_gi[2],
               inc_day_gi_3 = info$index$inc_day_gi[3],
               inc_day_gi_4 = info$index$inc_day_gi[4],
               inc_day_gi_5 = info$index$inc_day_gi[5],
               
               inc_rate_1 = info$index$inc_day_all[1],
               inc_rate_2 = info$index$inc_day_all[2],
               inc_rate_3 = info$index$inc_day_all[3],
               inc_rate_4 = info$index$inc_day_all[4],
               inc_rate_5 = info$index$inc_day_all[5],
               
               pop_by4age1 = info$index$pop_by4age[1],
               pop_by4age2 = info$index$pop_by4age[2],
               pop_by4age3 = info$index$pop_by4age[3],
               pop_by4age4 = info$index$pop_by4age[4],
               pop_by4age5 = info$index$pop_by4age[5],
               
               cumm_incday_gii4_1 = info$index$infections_day_gii4[1],
               cumm_incday_gii4_2 = info$index$infections_day_gii4[2],
               cumm_incday_gii4_3 = info$index$infections_day_gii4[3],
               cumm_incday_gii4_4 = info$index$infections_day_gii4[4],
               
               cumm_incday_gii_1 = info$index$infections_day_gii[1],
               cumm_incday_gii_2 = info$index$infections_day_gii[2],
               cumm_incday_gii_3 = info$index$infections_day_gii[3],
               cumm_incday_gii_4 = info$index$infections_day_gii[4],
               
               cumm_incday_gi3_1 = info$index$infections_day_gi3[1],
               cumm_incday_gi3_2 = info$index$infections_day_gi3[2],
               cumm_incday_gi3_3 = info$index$infections_day_gi3[3],
               cumm_incday_gi3_4 = info$index$infections_day_gi3[4],
               
               cumm_incday_gi_1 = info$index$infections_day_gi[1],
               cumm_incday_gi_2 = info$index$infections_day_gi[2],
               cumm_incday_gi_3 = info$index$infections_day_gi[3],
               cumm_incday_gi_4 = info$index$infections_day_gi[4],
               
               reported_wk_1 = info$index$reported_wk[1],
               reported_wk_2 = info$index$reported_wk[2],
               reported_wk_3 = info$index$reported_wk[3],
               reported_wk_4 = info$index$reported_wk[4],
               reported_wk_5 = info$index$reported_wk[5],
               reported_wk_6 = info$index$reported_wk[6],
               reported_wk_7 = info$index$reported_wk[7],
               reported_wk_8 = info$index$reported_wk[8],
               reported_wk_9 = info$index$reported_wk[9],
               reported_wk_10 = info$index$reported_wk[10],
               reported_wk_11 = info$index$reported_wk[11],
               reported_wk_12 = info$index$reported_wk[12],
               reported_wk_13 = info$index$reported_wk[13],
               reported_wk_14 = info$index$reported_wk[14],
               
               seroprev_num1.2 = info$index$seroprev_num[2],
               seroprev_num2.3 = info$index$seroprev_num[3],
               seroprev_num3.4 = info$index$seroprev_num[4],
               seroprev_num4.5 = info$index$seroprev_num[5],
               seroprev_num5.6 = info$index$seroprev_num[6],
               seroprev_num6.7 = info$index$seroprev_num[7],
               seroprev_den1.2 = info$index$seroprev_den[2],
               seroprev_den2.3 = info$index$seroprev_den[3],
               seroprev_den3.4 = info$index$seroprev_den[4],
               seroprev_den4.5 = info$index$seroprev_den[5],
               seroprev_den5.6 = info$index$seroprev_den[6],
               seroprev_den6.7 = info$index$seroprev_den[7],
               
               seroprev1.2 = info$index$seroprev[2],
               seroprev2.3 = info$index$seroprev[3],
               seroprev3.4 = info$index$seroprev[4],
               seroprev4.5 = info$index$seroprev[5],
               seroprev5.6 = info$index$seroprev[6],
               seroprev6.7 = info$index$seroprev[7],
               
               sus_gi3 = info$index$sus_gi3,
               sus_gi = info$index$sus_gi,
               sus_gii4 = info$index$sus_gii4,
               sus_gii = info$index$sus_gii

               ),
       
       state = c(
         t = info$index$time,
         inc_day_gi3 = info$index$infections_day_gi3,
         inc_day_gi = info$index$infections_day_gi,
         inc_day_gii4 = info$index$infections_day_gii4,
         inc_day_gii = info$index$infections_day_gii,
         
         inc_year_gii4_1 = info$index$inc_year_gii4[1],
         inc_year_gii4_2 = info$index$inc_year_gii4[2],
         inc_year_gii4_3 = info$index$inc_year_gii4[3],
         inc_year_gii4_4 = info$index$inc_year_gii4[4],
         inc_year_gii4_5 = info$index$inc_year_gii4[5],
         
         inc_year_gii_1 = info$index$inc_year_gii[1],
         inc_year_gii_2 = info$index$inc_year_gii[2],
         inc_year_gii_3 = info$index$inc_year_gii[3],
         inc_year_gii_4 = info$index$inc_year_gii[4],
         inc_year_gii_5 = info$index$inc_year_gii[5],
         
         inc_year_gi3_1 = info$index$inc_year_gi3[1],
         inc_year_gi3_2 = info$index$inc_year_gi3[2],
         inc_year_gi3_3 = info$index$inc_year_gi3[3],
         inc_year_gi3_4 = info$index$inc_year_gi3[4],
         inc_year_gi3_5 = info$index$inc_year_gi3[5],
         
         inc_year_gi_1 = info$index$inc_year_gi[1],
         inc_year_gi_2 = info$index$inc_year_gi[2],
         inc_year_gi_3 = info$index$inc_year_gi[3],
         inc_year_gi_4 = info$index$inc_year_gi[4],
         inc_year_gi_5 = info$index$inc_year_gi[5],
         
         inc_day_gii4_1 = info$index$inc_day_gii4[1],
         inc_day_gii4_2 = info$index$inc_day_gii4[2],
         inc_day_gii4_3 = info$index$inc_day_gii4[3],
         inc_day_gii4_4 = info$index$inc_day_gii4[4],
         inc_day_gii4_5 = info$index$inc_day_gii4[5],
         
         inc_day_gii_1 = info$index$inc_day_gii[1],
         inc_day_gii_2 = info$index$inc_day_gii[2],
         inc_day_gii_3 = info$index$inc_day_gii[3],
         inc_day_gii_4 = info$index$inc_day_gii[4],
         inc_day_gii_5 = info$index$inc_day_gii[5],
         
         inc_day_gi3_1 = info$index$inc_day_gi3[1],
         inc_day_gi3_2 = info$index$inc_day_gi3[2],
         inc_day_gi3_3 = info$index$inc_day_gi3[3],
         inc_day_gi3_4 = info$index$inc_day_gi3[4],
         inc_day_gi3_5 = info$index$inc_day_gi3[5],
         
         inc_day_gi_1 = info$index$inc_day_gi[1],
         inc_day_gi_2 = info$index$inc_day_gi[2],
         inc_day_gi_3 = info$index$inc_day_gi[3],
         inc_day_gi_4 = info$index$inc_day_gi[4],
         inc_day_gi_5 = info$index$inc_day_gi[5],
         
         inc_rate_1 = info$index$inc_day_all[1],
         inc_rate_2 = info$index$inc_day_all[2],
         inc_rate_3 = info$index$inc_day_all[3],
         inc_rate_4 = info$index$inc_day_all[4],
         inc_rate_5 = info$index$inc_day_all[5],
         
         pop_by4age1 = info$index$pop_by4age[1],
         pop_by4age2 = info$index$pop_by4age[2],
         pop_by4age3 = info$index$pop_by4age[3],
         pop_by4age4 = info$index$pop_by4age[4],
         pop_by4age5 = info$index$pop_by4age[5],
         cumm_incday_gii4_1 = info$index$infections_day_gii4[1],
         cumm_incday_gii4_2 = info$index$infections_day_gii4[2],
         cumm_incday_gii4_3 = info$index$infections_day_gii4[3],
         cumm_incday_gii4_4 = info$index$infections_day_gii4[4],
         
         cumm_incday_gii_1 = info$index$infections_day_gii[1],
         cumm_incday_gii_2 = info$index$infections_day_gii[2],
         cumm_incday_gii_3 = info$index$infections_day_gii[3],
         cumm_incday_gii_4 = info$index$infections_day_gii[4],
         
         cumm_incday_gi3_1 = info$index$infections_day_gi3[1],
         cumm_incday_gi3_2 = info$index$infections_day_gi3[2],
         cumm_incday_gi3_3 = info$index$infections_day_gi3[3],
         cumm_incday_gi3_4 = info$index$infections_day_gi3[4],
         
         cumm_incday_gi_1 = info$index$infections_day_gi[1],
         cumm_incday_gi_2 = info$index$infections_day_gi[2],
         cumm_incday_gi_3 = info$index$infections_day_gi[3],
         cumm_incday_gi_4 = info$index$infections_day_gi[4],
         
         reported_wk_1 = info$index$reported_wk[1],
         reported_wk_2 = info$index$reported_wk[2],
         reported_wk_3 = info$index$reported_wk[3],
         reported_wk_4 = info$index$reported_wk[4],
         reported_wk_5 = info$index$reported_wk[5],
         reported_wk_6 = info$index$reported_wk[6],
         reported_wk_7 = info$index$reported_wk[7],
         reported_wk_8 = info$index$reported_wk[8],
         reported_wk_9 = info$index$reported_wk[9],
         reported_wk_10 = info$index$reported_wk[10],
         reported_wk_11 = info$index$reported_wk[11],
         reported_wk_12 = info$index$reported_wk[12],
         reported_wk_13 = info$index$reported_wk[13],
         reported_wk_14 = info$index$reported_wk[14],
         seroprev_num1.2 = info$index$seroprev_num[2],
         seroprev_num2.3 = info$index$seroprev_num[3],
         seroprev_num3.4 = info$index$seroprev_num[4],
         seroprev_num4.5 = info$index$seroprev_num[5],
         seroprev_num5.6 = info$index$seroprev_num[6],
         seroprev_num6.7 = info$index$seroprev_num[7],
         seroprev_den1.2 = info$index$seroprev_den[2],
         seroprev_den2.3 = info$index$seroprev_den[3],
         seroprev_den3.4 = info$index$seroprev_den[4],
         seroprev_den4.5 = info$index$seroprev_den[5],
         seroprev_den5.6 = info$index$seroprev_den[6],
         seroprev_den6.7 = info$index$seroprev_den[7],
         
         seroprev1.2 = info$index$seroprev[2],
         seroprev2.3 = info$index$seroprev[3],
         seroprev3.4 = info$index$seroprev[4],
         seroprev4.5 = info$index$seroprev[5],
         seroprev5.6 = info$index$seroprev[6],
         seroprev6.7 = info$index$seroprev[7],
         
         sus_gi3 = info$index$sus_gi3,
         sus_gi = info$index$sus_gi,
         sus_gii4 = info$index$sus_gii4,
         sus_gii = info$index$sus_gii
         
       )
  )
}


compare_empty <- function(state, observed, pars = NULL) {
  
  llk<- -100
}



run_vaccine_scenario<-function(
    parameters, 
    states,
    footransform
){
  
  # Run simulation
  n_times<- (10 * 365)  
  tt<-seq(1,n_times,1)
  reps<-20
  sims<-dim(parameters)[1]
  
  all_inf_day<-matrix(0,sims*reps,n_times)
  gi3_inf_day<-matrix(0,sims*reps,n_times)
  gi_inf_day<-matrix(0,sims*reps,n_times)
  gii4_inf_day<-matrix(0,sims*reps,n_times)
  gii_inf_day<-matrix(0,sims*reps,n_times)
  vaccines_day<-matrix(0,sims*reps,n_times)
  
  for (ii in 1:sims){
    
    print(100*(ii/sims))
    
    input<-footransform(parameters[ii,]) 
    
    model <- seiar$new(input, 0, reps)
    model$update_state(state = states[,ii])
    sim<- model$simulate(tt)
    id<-model$info()$index
    span<- (1:reps)+ reps*(ii-1)
    gi3_inf_day[span,]<-t(drop(sim[id$infections_day_gi3,,]))
    gi_inf_day[span,]<-t(drop(sim[id$infections_day_gi,,]))
    gii4_inf_day[span,]<-t(drop(sim[id$infections_day_gii4,,]))
    gii_inf_day[span,]<-t(drop(sim[id$infections_day_gii,,]))
    vaccines_day[span,]<-t(drop(sim[id$vaccines_perday,,]))
    
    all_inf_day[span,]<- gi3_inf_day[span,]+gi_inf_day[span,]+
      gii4_inf_day[span,]+gii_inf_day[span,]
    
    
    
  }
  
  res<-list(
    time = tt,
    all_inf_day=all_inf_day,
    gi3_inf_day=gi3_inf_day,
    gi_inf_day=gi_inf_day,
    gii4_inf_day=gii4_inf_day,
    gii_inf_day=gii_inf_day,
    vaccines_day=vaccines_day
  )
  
  return(res)
  
}


















## Plot fits 
plot_fits<-function(sims,data){
  
  ## Community incidence (IID2)
  t<-which(sims["t",1,]%in%
             data$time_end[which(!is.na(data$cases_a1))])
  
  irates<-1000*(cbind(sims["cumu_inc1",,t],
                      sims["cumu_inc2",,t],
                      sims["cumu_inc3",,t],
                      sims["cumu_inc4",,t])/
                  (cbind(sims["n_age1",,t],
                         sims["n_age2",,t],
                         sims["n_age3",,t],
                         sims["n_age4",,t])/365))
  
  irate_obs<-c(data$cases_a1[1],data$cases_a2[1],
               data$cases_a3[1],data$cases_a4[1])
  matplot(c(1,2,3,4),t(irates), type = "p", col = "#00000011", 
          xlab = "Age", ylab = "Incidence per 1000", las = 1,ylim=c(0,250),xaxt="n")
  
  xtick<-seq(1, 4, by=1)
  axis(side=1, at=xtick, labels = c("0_4","5_14","15_64","65+"))
  arrows(x0=c(1,2,3,4), y0=data_iid2.c4$CI_lower, 
         x1=c(1,2,3,4), y1=data_iid2.c4$CI_upper,
         code=3, angle=90, length=0.1)
  points(irate_obs , pch = 19, col = "red")
  
  ## Weekly cass reported by UKHSA
  id<-which(sims["t",1,]%in%data_all$time_end[which(!is.na(data_all$reported))])
  reported<-sims["reported_wk",,id]
  reported_obs<-data$reported[which(!is.na(data$reported))]
  t<-seq(1,length(id),by=1)#sims["t",1,id]
  matplot(t,t(reported), type = "l", col = "#00000011", 
          xlab = "week", ylab = "Weekly reported cases", las = 1,ylim = c(0,200))
  points(t,reported_obs, pch = 19, col = "red")
  
  ## Seroprevalence 
  
  id<-which(sims["t",1,]%in%
              data$time_end[which(!is.na(data$sero1))])
  sero_model<-rbind(
    sims['seroprev1.2',,id],
    sims['seroprev2.3',,id],
    sims['seroprev3.4',,id],
    sims['seroprev4.5',,id],
    sims['seroprev5.6',,id],
    sims['seroprev6.7',,id]
  )
  id<-which(!is.na(data$sero1))
  sero_obs<-c(data$sero1[id],
              data$sero2[id],
              data$sero3[id],
              data$sero4[id],
              data$sero5[id],
              data$sero6[id])
  
  matplot(c(1,2,3,4,5,6),sero_model, type = "p", col = "#00000011", 
          xlab = "Age", ylab = "Seropositivity", las = 1,ylim=c(0,1),xaxt="n")
  xtick<-seq(1, 6, by=1)
  axis(side=1, at=xtick, labels = c("0_1","1_2","2_3","3_4","5_6","6_7"))
  arrows(x0=c(1,2,3,4,5,6), y0=sero$V2, 
         x1=c(1,2,3,4,5,6), y1=sero$V3,
         code=3, angle=90, length=0.1)
  points(sero_obs , pch = 19, col = "red")
  
}



# Demographicmodel
run_demog_model<-function(mu, ini=init, p=params, times=365*5, seiar_inst=seiar){
  
  
  
  n_particles <- 1L
  c_mat<-p$transmission
  c_mat[p$infa_id,p$infa_id]<-c_mat[p$infa_id,p$infa_id]*p$aduRR
  
  c_mat2<-p$transmission_holi
  c_mat2[p$infa_id,p$infa_id]<-c_mat2[p$infa_id,p$infa_id]*p$aduRR
  
  
  mort<-c(mu)/10000
  
  pars = list(
    beta = 0.0 ,   # transm coefficient
    repfac = 287,
    rho   = p$rho, # rel infect asymptomatic 
    init  = ini,
    mu    = mort/365,
    m     = c_mat,
    m_holi= c_mat2,
    aging_mat= p$aging_mat, 
    N_age = p$N_age,
    w1 = p$w1,
    school= as.double(p$school_uk))
  
  
  model <- seiar_inst$new(pars, 0, 1)
  n_times<-times
  tt<-seq(1,n_times,1)
  x<- model$simulate(tt)
  x<-drop(x)
  
  time <- x[1, ]
  
  # Plotting the trajectories
  idx<-model$info()$index
  
  # Check age distribution
  age.init<-rbind(
    x[idx$M,dim(x)[2]],
    x[idx$G,dim(x)[2]],
    x[idx$S,dim(x)[2]],
    x[idx$E,dim(x)[2]],
    x[idx$I,dim(x)[2]],
    x[idx$A,dim(x)[2]],
    x[idx$R,dim(x)[2]])
  
  age.distr.sim=colSums(age.init)/sum(age.init)
  idx<-model$info()$index
  results<-list(
    idx = idx,
    states=x,
    age.distr.sim=age.distr.sim
    
  )
  
  return(results)
}

########################
## Run model

run_model<-function(pars, times, seiar=model.seiar){
  
  n_particles <- 3L
  dt <- pars$dt
  seed<-1
  g.ini<-round(pars$pop*pars$p_nonsecretor)
  s.ini<-pars$pop-g.ini
  i.ini<-pars$contact$demography$population*0
  i.ini[2]<-seed
  mort.rates<- pars$mu
  c_mat<-pars$transmission
  c_mat[pars$infa_id,pars$infa_id]<-c_mat[pars$infa_id,pars$infa_id]*pars$aduRR
  
  model <- seiar$new(pars = list(dt = dt,
                                 M_ini= pars$contact$demography$population*0,
                                 G_ini= g.ini,
                                 S_ini=s.ini,
                                 E_ini=pars$contact$demography$population*0,
                                 I_ini=i.ini,
                                 A_ini=pars$contact$demography$population*0,
                                 R_ini=pars$contact$demography$population*0,
                                 age_select=pars$age_select,
                                 beta = pars$beta,   # transm coefficient
                                 rho   = pars$rho, # rel infect asymptomatic 
                                 mu    = mort.rates,
                                 m = c_mat,
                                 aging_mat= pars$aging_mat, 
                                 N_age = pars$N_age,
                                 w1 = pars$w1),
                     time = 0,
                     n_particles = n_particles,
                     n_threads = 4L,
                     seed = 1L)
  
  # Define how long the model runs for, number of time steps
  n_times <- times
  x <- array(NA, dim = c(model$info()$len,n_particles, n_times))
  
  # For loop to run the model iteratively
  
  # for (t in seq_len(n_times))
  #   {
  # x[ ,  ,t] <- model$run(t)
  # }
  
  tt<-seq(1/pars$dt,n_times,1/pars$dt)# output only the day step 
  x<- model$simulate(tt)
  
  
  time <- x[1, 1, ]
  # Plotting the trajectories
  idx<-model$info()$index
  
  
  # Check age distribution
  age.init<-rbind(
    x[idx$M,1 ,dim(x)[3]],
    x[idx$G,1 ,dim(x)[3]],
    x[idx$S,1 ,dim(x)[3]],
    x[idx$E,1 ,dim(x)[3]],
    x[idx$I,1 ,dim(x)[3]],
    x[idx$A,1 ,dim(x)[3]],
    x[idx$R,1 ,dim(x)[3]])
  
  age.distr.sim=colSums(age.init)/sum(age.init)
  results<-list(
    idx=idx,
    states=x,
    age.distr.sim=age.distr.sim
  )
  
  return(results)
}


get_output<-function(theta,pars, ntimes, seiar=model.seiar){
  
  
  nruns <- nrow(theta)
  
  
  # Allocate memory
  irate <- matrix(NA, nrow = nruns, ncol = 4)
  
  
  for (jj in 1:nruns) {
    
    pars$beta<-theta[['beta']][jj]
    pars$aduRR<-theta[['aduRR']][jj]
    
    sim<-run_model(pars,ntimes, seiar)
    
    states<-sim$states
    idday<-seq(1,ntimes, 1/pars$dt)
    idx<-sim$idx
    
    N<- states[idx$M,1,]+
      states[idx$G,1,]+
      states[idx$S,1,]+
      states[idx$R,1,]
    
    st<-which(days_vec==idd2_startdate)* 1/pars$dt
    ed<-which(days_vec==idd2_enddate) * 1/pars$dt
    span<-seq(st,ed,1)
    cases<-states[idx$cumu_inc,1,span]
    cases<-tail(t(cases),1)-head(t(cases),1)
    cases.4<-c(sum(cases[c(1,2)]),cases[3],sum(cases[c(4,5,6,7,8)]),sum(cases[c(9,10)]))
    
    PY<-(states[idx$M,1,span]+
           states[idx$G,1,span]+
           states[idx$S,1,span]+
           states[idx$R,1,span])
    
    PY<-apply(PY,1,cumsum)
    PY<-tail(PY,1)-head(PY,1)
    PY<-PY/(365/pars$dt)
    PY.4<-c(sum(PY[c(1,2)]),PY[3],sum(PY[c(4,5,6,7,8)]),sum(PY[c(9,10)]))
    
    irate[jj,]<-1000*cases.4/PY.4
    
  }
  
  
  out<-list(
    irate_pyear=irate#irate_pyear
  )
  
  return(out)
}


################################
# Plot age dsitributions
plot_age_distr<-function(sim_age, pars=pars){
  df<-data.frame(
    age=factor(pars$contact$demography$age.group),
    age.distr.sim=sim_age,
    age.distr.ons=pars$contact$demography$proportion
  )
  
  df$age <- factor(df$age, levels = df$age)
  
  df.melt<- reshape2::melt(df) 
  
  p1<-ggplot(data = df.melt, aes(x=age,y=value, fill=variable) )+
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_minimal()
  p1 + scale_fill_manual(values=c('#999999','#E69F00'))
}


