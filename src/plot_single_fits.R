plot_single_fits<-function(sims,data){
  
 
  data_col="dodgerblue"
  ## Community incidence (IID2)
  t<-which(sims["t",1,]%in%
             data_all$time_end[which(!is.na(data_all$cases_a1))])


  # rate1<-mean(sims['inc_rate_1',,(t-365):t])
  # rate2<-mean(sims['inc_rate_2',,(t-365):t])
  # rate3<-mean(sims['inc_rate_3',,(t-365):t])
  # rate4<-mean(sims['inc_rate_4',,(t-365):t])
  # rate5<-mean(sims['inc_rate_5',,(t-365):t])
  # 

  irates <-(cbind(
    (  sims['inc_year_gi3_1',,]+
       sims['inc_year_gi_1',,]+
       sims['inc_year_gii4_1',,]+
       sims['inc_year_gii_1',,])/(sims['pop_by4age1',,]),
    
    (  sims['inc_year_gi3_2',,]+
       sims['inc_year_gi_2',,]+
       sims['inc_year_gii4_2',,]+
       sims['inc_year_gii_2',,])/(sims['pop_by4age2',,]),
    
    (  sims['inc_year_gi3_3',,]+
       sims['inc_year_gi_3',,]+
       sims['inc_year_gii4_3',,]+
       sims['inc_year_gii_3',,])/(sims['pop_by4age3',,]) ,
    
    (  sims['inc_year_gi3_4',,]+
       sims['inc_year_gi_4',,]+
       sims['inc_year_gii4_4',,]+
       sims['inc_year_gii_4',,])/(sims['pop_by4age4',,]) ,
    
    (  sims['inc_year_gi3_5',,]+
       sims['inc_year_gi_5',,]+
       sims['inc_year_gii4_5',,]+
       sims['inc_year_gii_5',,])/(sims['pop_by4age5',,])))
  
  irates<-1000*t(irates)

  irates<-irates[,t-1]
  
  irate_obs<- data_iid2.c4$per1000personyears
               # c(data$cases_a1[t-1],
               # data$cases_a2[t-1],
               # data$cases_a3[t-1],
               # data$cases_a4[t-1],
               # data$cases_a5[t-1])
 
  np<-length(irate_obs)
  xtick<-seq(1, np, by=1)

  matplot(xtick,(irates), pch=18, col =  data_col, cex=2,
          xlab = "Age", ylab = "Incidence per 1000", las = 1,
          ylim=c(0,max(max(data_iid2.c4$CI_uppe),max(irates))),xaxt="n")
  
  axis(side=1, at=xtick, labels = c("0_1","1_4","5_14","15_64","65+"))
  arrows(x0=xtick, y0=data_iid2.c4$CI_lower, 
         x1=xtick, y1=data_iid2.c4$CI_upper,
         code=3, angle=90, length=0.1)
  points(xtick,irate_obs , pch = 19, col = "black")
  


  ## sero 
  id<-which(sims["t",1,]%in%
              data_all$time_end[which(!is.na(data_all$sero1))])
  # sero_model<-rbind(
  #   sims['seroprev_num1.2',,id]/sims['seroprev_den1.2',,id] ,
  #   sims['seroprev_num2.3',,id]/sims['seroprev_den2.3',,id] ,
  #   sims['seroprev_num3.4',,id]/sims['seroprev_den3.4',,id] ,
  #   sims['seroprev_num4.5',,id]/sims['seroprev_den4.5',,id] ,
  #   sims['seroprev_num5.6',,id]/sims['seroprev_den5.6',,id] ,
  #   sims['seroprev_num6.7',,id]/sims['seroprev_den6.7',,id] )
  
  sero_model<-rbind(
    sims['seroprev1.2',,id] ,
    sims['seroprev2.3',,id] ,
    sims['seroprev3.4',,id] ,
    sims['seroprev4.5',,id] ,
    sims['seroprev5.6',,id] ,
    sims['seroprev6.7',,id])
  
  
  
  id<-which(!is.na(data_all$sero1))
  sero_obs<-c(data_all$sero1[id],
              data_all$sero2[id],
              data_all$sero3[id],
              data_all$sero4[id],
              data_all$sero5[id],
              data_all$sero6[id])
  
  matplot(c(1,2,3,4,5,6),sero_model, pch=18, col =  data_col,cex=2, 
          xlab = "Age", ylab = "Seropositivity", las = 1,ylim=c(0,1),xaxt="n")
  xtick<-seq(1, 6, by=1)
  axis(side=1, at=xtick, labels = c("1_2","2_3","3_4","4_5","5_6","6_7"))
  arrows(x0=c(1,2,3,4,5,6), y0=sero$V2, 
         x1=c(1,2,3,4,5,6), y1=sero$V3,
         code=3, angle=90, length=0.1)
  points(sero_obs , pch = 19, col = "black")
  
  
  

  #
  reported_sim<-drop(
    sims["reported_wk_1",,] +
    sims["reported_wk_2",,] +
    sims["reported_wk_3",,] +
    sims["reported_wk_4",,] +
    sims["reported_wk_5",,] +
    sims["reported_wk_6",,] +
    sims["reported_wk_7",,] +
    sims["reported_wk_8",,] +
    sims["reported_wk_9",,] +
    sims["reported_wk_10",,] +
    sims["reported_wk_11",,] +
    sims["reported_wk_12",,] +
    sims["reported_wk_13",,] +
    sims["reported_wk_14",,])
  
  reported_obs_y<-total_cases$cases
  reported_obs_x<-total_cases$date
  

  matplot(days_vec,reported_sim[2:dim(sims)[3]], type = "l", col = data_col,  
          xlab = "week", ylab = "Weekly reported cases", las = 1,
          xlim = c(dates2[1],dates2[length(dates2)]),
          ylim = c(0 , max(reported_obs_y)))
  points(reported_obs_x, reported_obs_y, pch = 19, col = "black")
  
  
  
  
  ## Weekly cass reported by UKHSA
  # id<-which(sims["t",1,]%in%data_all$time_end[which(!is.na(data_all$reported))])
  # 
  # days<-days_vec[id]
  # 
  # reported<-
  #   sims["reported_wk_1",,id] +
  #   sims["reported_wk_2",,id] +
  #   sims["reported_wk_3",,id] +
  #   sims["reported_wk_4",,id] +
  #   sims["reported_wk_5",,id] +
  #   sims["reported_wk_6",,id] +
  #   sims["reported_wk_7",,id] +
  #   sims["reported_wk_8",,id] +
  #   sims["reported_wk_9",,id] +
  #   sims["reported_wk_10",,id] +
  #   sims["reported_wk_11",,id] +
  #   sims["reported_wk_12",,id] +
  #   sims["reported_wk_13",,id] +
  #   sims["reported_wk_14",,id]
  # 
  # reported_obs<-data_all$reported[which(!is.na(data_all$reported))]
  # t<-sims["t",1,id]
  # matplot(days,reported, type = "l", col = data_col,  
  #         xlab = "week", ylab = "Weekly reported cases", las = 1,
  #         ylim = c(0 , max(reported_obs)))
  # points(reported_obs_x, reported_obs_y, pch = 19, col = "black")
  
  
  ##covid
  # matplot(t,reported, type = "l", col = data_col,  
  #         xlab = "week", ylab = "Weekly reported cases", las = 1,
  #         ylim = c(0 , 600),xlim=c(which(days_vec%in%covid_start)-365, which(days_vec%in%covid_end)+365))
  # points(t, reported_obs, pch = 19, col = "black")
  # lines(c(which(days_vec%in%covid_start),which(days_vec%in%covid_start)),
  #       c(0,500), col="red")
  # lines(c(which(days_vec%in%covid_end),which(days_vec%in%covid_end)),
  #       c(0,500), col="red")
  # 
  # lines(c(which(days_vec%in%lock_start1),which(days_vec%in%lock_start1)),
  #       c(0,500), col="blue")
  # lines(c(which(days_vec%in%lock_end1),which(days_vec%in%lock_end1)),
  #       c(0,500), col="blue")
  # lines(c(which(days_vec%in%lock_start1),which(days_vec%in%lock_end1)),
  #       c(500,500), col="blue")
  # 
  # lines(c(which(days_vec%in%lock_start2),which(days_vec%in%lock_start2)),
  #       c(0,500), col="black")
  # lines(c(which(days_vec%in%lock_end2),which(days_vec%in%lock_end2)),
  #       c(0,500), col="black")
  # lines(c(which(days_vec%in%lock_start2),which(days_vec%in%lock_end2)),
  #       c(500,500), col="black")
  # 
  # lines(c(which(days_vec%in%lock_start3),which(days_vec%in%lock_start3)),
  #       c(0,500), col="purple")
  # lines(c(which(days_vec%in%lock_end3),which(days_vec%in%lock_end3)),
  #       c(0,500), col="purple")
  # lines(c(which(days_vec%in%lock_start3),which(days_vec%in%lock_end3)),
  #       c(500,500), col="purple")
  # 
  # lines(c(which(days_vec%in%lock_start4),which(days_vec%in%lock_start4)),
  #       c(0,500), col="orange")
  # lines(c(which(days_vec%in%lock_end4),which(days_vec%in%lock_end4)),
  #       c(0,500), col="orange")
  # lines(c(which(days_vec%in%lock_start4),which(days_vec%in%lock_end4)),
  #       c(500,500), col="orange")
  # 
  # lines(c(which(days_vec%in%lock_start5),which(days_vec%in%lock_start5)),
  #       c(0,500), col="forestgreen")
  # lines(c(which(days_vec%in%lock_end5),which(days_vec%in%lock_end5)),
  #       c(0,500), col="forestgreen")
  # lines(c(which(days_vec%in%lock_start5),which(days_vec%in%lock_end5)),
  #       c(500,500), col="forestgreen")
  # 
  # lines(c(which(days_vec%in%lock_start6),which(days_vec%in%lock_start6)),
  #       c(0,500), col="firebrick")
  # lines(c(which(days_vec%in%lock_end6),which(days_vec%in%lock_end6)),
  #       c(0,500), col="firebrick")
  # lines(c(which(days_vec%in%lock_start6),which(days_vec%in%lock_end6)),
  #       c(500,500), col="firebrick")
  # 
  # lines(c(which(days_vec%in%lock_start7),which(days_vec%in%lock_start7)),
  #       c(0,500), col="navy")
  # lines(c(which(days_vec%in%lock_end7),which(days_vec%in%lock_end7)),
  #       c(0,500), col="navy")
  # lines(c(which(days_vec%in%lock_start7),which(days_vec%in%lock_end7)),
  #       c(500,500), col="navy")
  # 
  # lines(c(which(days_vec%in%lock_start8),which(days_vec%in%lock_start8)),
  #       c(0,500), col="brown")
  # lines(c(which(days_vec%in%lock_end8),which(days_vec%in%lock_end8)),
  #       c(0,500), col="brown")
  # lines(c(which(days_vec%in%lock_start8),which(days_vec%in%lock_end8)),
  #       c(500,500), col="brown")
  # 
  # lines(c(which(days_vec%in%lock_start9),which(days_vec%in%lock_start9)),
  #       c(0,500), col="violet")
  # lines(c(which(days_vec%in%lock_end9),which(days_vec%in%lock_end9)),
  #       c(0,500), col="violet")
  # lines(c(which(days_vec%in%lock_start9),which(days_vec%in%lock_end9)),
  #       c(500,500), col="violet")
  
  

  
  # Strains 
  id<-which(sims["t",1,]%in%
              data$time_end[which(!is.na(data$gi_prop))])
  
  
  gi3_model<-((sims['cumm_incday_gi3_1',,id]+
                 sims['cumm_incday_gi3_2',,id]+
                 sims['cumm_incday_gi3_3',,id]+
                 sims['cumm_incday_gi3_4',,id]) /
                (    
                     sims['cumm_incday_gii4_1',,id]+
                     sims['cumm_incday_gii4_2',,id]+
                     sims['cumm_incday_gii4_3',,id]+
                     sims['cumm_incday_gii4_4',,id]+
                     sims['cumm_incday_gii_1',,id]+
                     sims['cumm_incday_gii_2',,id]+
                     sims['cumm_incday_gii_3',,id]+
                     sims['cumm_incday_gii_4',,id]+
                     sims['cumm_incday_gi3_1',,id]+
                     sims['cumm_incday_gi3_2',,id]+
                     sims['cumm_incday_gi3_3',,id]+
                     sims['cumm_incday_gi3_4',,id]+
                     sims['cumm_incday_gi_1',,id]+
                     sims['cumm_incday_gi_2',,id]+
                     sims['cumm_incday_gi_3',,id]+
                     sims['cumm_incday_gi_4',,id]))
  
  gi_model<-((sims['cumm_incday_gi_1',,id]+
                sims['cumm_incday_gi_2',,id]+
                sims['cumm_incday_gi_3',,id]+
                sims['cumm_incday_gi_4',,id]) /
               (    
                 sims['cumm_incday_gii4_1',,id]+
                   sims['cumm_incday_gii4_2',,id]+
                   sims['cumm_incday_gii4_3',,id]+
                   sims['cumm_incday_gii4_4',,id]+
                   sims['cumm_incday_gii_1',,id]+
                   sims['cumm_incday_gii_2',,id]+
                   sims['cumm_incday_gii_3',,id]+
                   sims['cumm_incday_gii_4',,id]+
                   sims['cumm_incday_gi3_1',,id]+
                   sims['cumm_incday_gi3_2',,id]+
                   sims['cumm_incday_gi3_3',,id]+
                   sims['cumm_incday_gi3_4',,id]+
                   sims['cumm_incday_gi_1',,id]+
                   sims['cumm_incday_gi_2',,id]+
                   sims['cumm_incday_gi_3',,id]+
                   sims['cumm_incday_gi_4',,id]))
  
  gii4_model<-((sims['cumm_incday_gii4_1',,id]+
                  sims['cumm_incday_gii4_2',,id]+
                  sims['cumm_incday_gii4_3',,id]+
                  sims['cumm_incday_gii4_4',,id]) /
                 (    
                   sims['cumm_incday_gii4_1',,id]+
                     sims['cumm_incday_gii4_2',,id]+
                     sims['cumm_incday_gii4_3',,id]+
                     sims['cumm_incday_gii4_4',,id]+
                     sims['cumm_incday_gii_1',,id]+
                     sims['cumm_incday_gii_2',,id]+
                     sims['cumm_incday_gii_3',,id]+
                     sims['cumm_incday_gii_4',,id]+
                     sims['cumm_incday_gi3_1',,id]+
                     sims['cumm_incday_gi3_2',,id]+
                     sims['cumm_incday_gi3_3',,id]+
                     sims['cumm_incday_gi3_4',,id]+
                     sims['cumm_incday_gi_1',,id]+
                     sims['cumm_incday_gi_2',,id]+
                     sims['cumm_incday_gi_3',,id]+
                     sims['cumm_incday_gi_4',,id]))
  
  gii_model<-((sims['cumm_incday_gii_1',,id]+
                 sims['cumm_incday_gii_2',,id]+
                 sims['cumm_incday_gii_3',,id]+
                 sims['cumm_incday_gii_4',,id]) /
                (    
                  sims['cumm_incday_gii4_1',,id]+
                    sims['cumm_incday_gii4_2',,id]+
                    sims['cumm_incday_gii4_3',,id]+
                    sims['cumm_incday_gii4_4',,id]+
                    sims['cumm_incday_gii_1',,id]+
                    sims['cumm_incday_gii_2',,id]+
                    sims['cumm_incday_gii_3',,id]+
                    sims['cumm_incday_gii_4',,id]+
                    sims['cumm_incday_gi3_1',,id]+
                    sims['cumm_incday_gi3_2',,id]+
                    sims['cumm_incday_gi3_3',,id]+
                    sims['cumm_incday_gi3_4',,id]+
                    sims['cumm_incday_gi_1',,id]+
                    sims['cumm_incday_gi_2',,id]+
                    sims['cumm_incday_gi_3',,id]+
                    sims['cumm_incday_gi_4',,id]))
  
  strain_model<-cbind(gi3_model,gi_model, gii4_model,gii_model )*100
  

  x_d <- c(1,3,5,7) # bin x axis positions
  
  id<-which(!is.na(data$gi_prop))
  strain_obs<-c(data$gi3_prop[id],
                data$gi_prop[id],
                data$gii4_prop[id],
                data$gii_prop[id])*100
  
  plot(x_d, strain_model, pch = 18, col = "darkcyan", 
       xlab = "Strain", ylab = "%",
       ylim=c(0,100),
       xaxt="n",
       cex=2)
  axis(side=1, at= x_d, labels = c("GI3","Other GI","GII4","Other GII"))
  points(x_d,strain_obs , pch=19,col = "black", cex=1)
  
  
  
  ### Age inc
  
  id<-which(sims["t",1,]%in%
              data$time_end[which(!is.na(data$a0_prop))])
  
  
  a0_model<-((sims['cumm_incday_gi3_1',,id]+
                sims['cumm_incday_gi_1',,id]+
                sims['cumm_incday_gii4_1',,id]+
                sims['cumm_incday_gii_1',,id]) /
               (    
                 sims['cumm_incday_gii4_1',,id]+
                   sims['cumm_incday_gii4_2',,id]+
                   sims['cumm_incday_gii4_3',,id]+
                   sims['cumm_incday_gii4_4',,id]+
                   sims['cumm_incday_gii_1',,id]+
                   sims['cumm_incday_gii_2',,id]+
                   sims['cumm_incday_gii_3',,id]+
                   sims['cumm_incday_gii_4',,id]+
                   sims['cumm_incday_gi3_1',,id]+
                   sims['cumm_incday_gi3_2',,id]+
                   sims['cumm_incday_gi3_3',,id]+
                   sims['cumm_incday_gi3_4',,id]+
                   sims['cumm_incday_gi_1',,id]+
                   sims['cumm_incday_gi_2',,id]+
                   sims['cumm_incday_gi_3',,id]+
                   sims['cumm_incday_gi_4',,id]))
  
  
  a5_model<-((sims['cumm_incday_gi3_2',,id]+
                sims['cumm_incday_gi_2',,id]+
                sims['cumm_incday_gii4_2',,id]+
                sims['cumm_incday_gii_2',,id]) /
               (    
                 sims['cumm_incday_gii4_1',,id]+
                   sims['cumm_incday_gii4_2',,id]+
                   sims['cumm_incday_gii4_3',,id]+
                   sims['cumm_incday_gii4_4',,id]+
                   sims['cumm_incday_gii_1',,id]+
                   sims['cumm_incday_gii_2',,id]+
                   sims['cumm_incday_gii_3',,id]+
                   sims['cumm_incday_gii_4',,id]+
                   sims['cumm_incday_gi3_1',,id]+
                   sims['cumm_incday_gi3_2',,id]+
                   sims['cumm_incday_gi3_3',,id]+
                   sims['cumm_incday_gi3_4',,id]+
                   sims['cumm_incday_gi_1',,id]+
                   sims['cumm_incday_gi_2',,id]+
                   sims['cumm_incday_gi_3',,id]+
                   sims['cumm_incday_gi_4',,id]))
  
  a15_model<-((sims['cumm_incday_gi3_3',,id]+
                 sims['cumm_incday_gi_3',,id]+
                 sims['cumm_incday_gii4_3',,id]+
                 sims['cumm_incday_gii_3',,id]) /
                (    
                  sims['cumm_incday_gii4_1',,id]+
                    sims['cumm_incday_gii4_2',,id]+
                    sims['cumm_incday_gii4_3',,id]+
                    sims['cumm_incday_gii4_4',,id]+
                    sims['cumm_incday_gii_1',,id]+
                    sims['cumm_incday_gii_2',,id]+
                    sims['cumm_incday_gii_3',,id]+
                    sims['cumm_incday_gii_4',,id]+
                    sims['cumm_incday_gi3_1',,id]+
                    sims['cumm_incday_gi3_2',,id]+
                    sims['cumm_incday_gi3_3',,id]+
                    sims['cumm_incday_gi3_4',,id]+
                    sims['cumm_incday_gi_1',,id]+
                    sims['cumm_incday_gi_2',,id]+
                    sims['cumm_incday_gi_3',,id]+
                    sims['cumm_incday_gi_4',,id]))
  
  a65_model<-((  sims['cumm_incday_gi3_4',,id]+
                 sims['cumm_incday_gi_4',,id]+
                 sims['cumm_incday_gii4_4',,id]+
                 sims['cumm_incday_gii_4',,id]) /
                (    
                    sims['cumm_incday_gii4_1',,id]+
                    sims['cumm_incday_gii4_2',,id]+
                    sims['cumm_incday_gii4_3',,id]+
                    sims['cumm_incday_gii4_4',,id]+
                    sims['cumm_incday_gii_1',,id]+
                    sims['cumm_incday_gii_2',,id]+
                    sims['cumm_incday_gii_3',,id]+
                    sims['cumm_incday_gii_4',,id]+
                    sims['cumm_incday_gi3_1',,id]+
                    sims['cumm_incday_gi3_2',,id]+
                    sims['cumm_incday_gi3_3',,id]+
                    sims['cumm_incday_gi3_4',,id]+
                    sims['cumm_incday_gi_1',,id]+
                    sims['cumm_incday_gi_2',,id]+
                    sims['cumm_incday_gi_3',,id]+
                    sims['cumm_incday_gi_4',,id]))
  
  
 
  
  age_model<-cbind(a0_model,a5_model, a15_model,a65_model )*100
  
  
  id<-which(!is.na(data$gi_prop))
  age_obs<-c(data$a0_prop[id],
             data$a5_prop[id],
             data$a15_prop[id],
             data$a65_prop[id])*100
  
  x_d <- c(1,3,5,7) # bin x axis positions
  
  plot(x_d, age_obs, pch = 19, col = "black", 
       xlab = "Age", ylab = "%",
       ylim=c(0,100),
       xaxt="n",
       cex=1)
  axis(side=1, at= x_d, labels = c("0_4","5_14","15_64","65p"))
  points(x_d,age_model , pch=18,col = data_col,cex=2)
  
  
  

  
  
  
  ### Strains long
  # 
  # 
  # id<-seq(2,length(days_vec)+1)
  # 
  # days<-days_vec
  # 
  # gi3<-
  #   sims["inc_day_gi3_1",,id] +
  #   sims["inc_day_gi3_2",,id] +
  #   sims["inc_day_gi3_3",,id] +
  #   sims["inc_day_gi3_4",,id] +
  #   sims["inc_day_gi3_5",,id]
  # 
  # 
  # gi<-
  #   sims["inc_day_gi_1",,id] +
  #   sims["inc_day_gi_2",,id] +
  #   sims["inc_day_gi_3",,id] +
  #   sims["inc_day_gi_4",,id] +
  #   sims["inc_day_gi_5",,id]
  # 
  # 
  # 
  # gii<-
  #   sims["inc_day_gii_1",,id] +
  #   sims["inc_day_gii_2",,id] +
  #   sims["inc_day_gii_3",,id] +
  #   sims["inc_day_gii_4",,id] +
  #   sims["inc_day_gii_5",,id]
  # 
  #   
  # gii4<-
  #   sims["inc_day_gii4_1",,id] +
  #   sims["inc_day_gii4_2",,id] +
  #   sims["inc_day_gii4_3",,id] +
  #   sims["inc_day_gii4_4",,id] +
  #   sims["inc_day_gii4_5",,id]
  # 
  # 
  # plot(days_vec,gii4, type = "l", col = "firebrick",  
  #         xlab = "date", ylab = "Daily reported cases", las = 1,
  #      xlim = c(as.Date("2012-01-01"),as.Date("2023-01-01")),
  #         ylim = c(0 , 1e5))
  # lines(days_vec,gii, type = "l", col = "navyblue")
  # lines(days_vec,gi3, type = "l", col = "darkgreen")
  # lines(days_vec,gi, type = "l", col = "purple")
  # 
  
  
  
  
  
}