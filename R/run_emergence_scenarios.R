
run_emergence_scenarios <- function(cross_p,
                                    transm_r,
                                    parameters,
                                    model,
                                    init_state0,
                                    init_state,
                                    times,
                                    time_vec,
                                    nsamples,
                                    footr) {
  
  runs_gi3<- matrix(0, nrow=length(times),ncol=nsamples);  
  runs_gi<- matrix(0, nrow=length(times),ncol=nsamples);
  runs_gii4<- matrix(0, nrow=length(times),ncol=nsamples);
  runs_gii<- matrix(0, nrow=length(times),ncol=nsamples);
  reported<- matrix(0, nrow=length(times),ncol=nsamples);
  
  runs_age1<- matrix(0, nrow=length(times),ncol=nsamples);   
  runs_age2<- matrix(0, nrow=length(times),ncol=nsamples);  
  runs_age3<- matrix(0, nrow=length(times),ncol=nsamples);  
  runs_age4<- matrix(0, nrow=length(times),ncol=nsamples);  
  
  
  if(transm_r==1 && cross_p==1){
    
    init=init_state0
  }else{
    
    init=init_state
    
  }
  
  
  
  for (jj in 1:nsamples){
    
    
    
    pa<-footr(parameters[jj,])
    pa$crossp_GII<-pa$crossp_GII*cross_p
    pa$beta_3<-pa$beta_3*transm_r
    model$update_state(
      pars=pa,
      state =  init[,jj],
      time = times[1]
    )
    
    runs1<-model$simulate(times)
    
    runs_gi3[,jj]<-colSums(runs1[id$inc_day_gi3,,])
    runs_gi[,jj]<-colSums(runs1[id$inc_day_gi,,])
    runs_gii4[,jj]<-colSums(runs1[id$inc_day_gii4,,])
    runs_gii[,jj]<-colSums(runs1[id$inc_day_gii,,])
    reported[,jj]<-colSums(runs1[id$reported_wk,,])
    
    
    runs_age1[,jj]  <-
      runs1[id$inc_day_gi3[1],,] + 
      runs1[id$inc_day_gi3[2],,]+
      runs1[id$inc_day_gi[1],,] + 
      runs1[id$inc_day_gi[2],,]+
      runs1[id$inc_day_gii4[1],,]+ 
      runs1[id$inc_day_gii4[2],,]+
      runs1[id$inc_day_gii[1],,] + 
      runs1[id$inc_day_gii[2],,]
    
    runs_age2[,jj]  <-
      runs1[id$inc_day_gi3[3],,] + 
      runs1[id$inc_day_gi[3],,] + 
      runs1[id$inc_day_gii4[3],,]+ 
      runs1[id$inc_day_gii[3],,]
    
    runs_age3[,jj]  <-
      runs1[id$inc_day_gi3[4],,] + 
      runs1[id$inc_day_gi[4],,] + 
      runs1[id$inc_day_gii4[4],,]+ 
      runs1[id$inc_day_gii[4],,]
    
    runs_age4[,jj]  <-
      runs1[id$inc_day_gi3[5],,] + 
      runs1[id$inc_day_gi[5],,] + 
      runs1[id$inc_day_gii4[5],,]+ 
      runs1[id$inc_day_gii[5],,]
    
    
    print(1e2*(jj/nsamples))
  }
  
  
  gi3 <- as.data.frame( 
    rowQuantiles((runs_gi3),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  gi <- as.data.frame( 
    rowQuantiles((runs_gi),
                 probs = c(0.025, 0.5, 0.975)
    ))
  gii4 <- as.data.frame( 
    rowQuantiles((runs_gii4),
                 probs = c(0.025, 0.5, 0.975)
    ))
  gii <- as.data.frame( 
    rowQuantiles((runs_gii),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  
  gi3$day<-times
  gi$day<-times
  gii4$day<-times
  gii$day<-times
  
  
  gii4$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  gi3$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  gi$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  gii$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  
  
  
  ## By age
  
  a1 <- as.data.frame( 
    rowQuantiles((runs_age1),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  a2 <- as.data.frame( 
    rowQuantiles((runs_age2),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  a3 <- as.data.frame( 
    rowQuantiles((runs_age3),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  a4 <- as.data.frame( 
    rowQuantiles((runs_age4),
                 probs = c(0.025, 0.5, 0.975)
    ))
  
  
  a1$day<-times
  a2$day<-times
  a3$day<-times
  a4$day<-times
  
  
  a1$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  a2$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  a3$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  a4$date<-seq(time_vec[length(time_vec)],length.out=length(times), by =1)
  
  
  
  res<-list(
    runs_gi3=runs_gi3,
    runs_gi=runs_gi,
    runs_gii4=runs_gii4,
    runs_gii=runs_gii,
    reported=reported , 
    gi3=gi3,
    gi=gi,
    gii4=gii4,
    gii=gii,
    a1=a1,
    a2=a2,
    a3=a3,
    a4=a4
  )
  return(res)
  
}