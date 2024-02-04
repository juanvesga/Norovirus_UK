library(matrixStats)
## Plot model fits
save_plot=0
fake=0
nsim <- nsamples
#ids<- which(processed_chains$probabilities[, "log_posterior"] > -9500) 
#posteriors<-processed_chains$trajectories$state#[,ids,]
# Sample
sims <- out$runs# posteriors[,sample(ncol(posteriors), nsim), ]

data<-data_all

######## IID2 Community incidence
x_d <- c(1, 3, 5, 7, 9) # bin x axis positions

df_d <- data.frame(
  x = x_d,
  inc = data_iid2.c4$per1000personyears,
  low = data_iid2.c4$CI_lower,
  up = data_iid2.c4$CI_upper
)

idx<-out$idx

t<-which(sims[idx$t,1,]%in%
           data$time_end[which(!is.na(data$cases_a1))])
irates <-(cbind(
  (  sims[idx$inc_year_gi3_1,,t]+
       sims[idx$inc_year_gi_1,,t]+
       sims[idx$inc_year_gii4_1,,t]+
       sims[idx$inc_year_gii_1,,t])/(sims[idx$pop_by4age1,,t]),
  
  (  sims[idx$inc_year_gi3_2,,t]+
       sims[idx$inc_year_gi_2,,t]+
       sims[idx$inc_year_gii4_2,,t]+
       sims[idx$inc_year_gii_2,,t])/(sims[idx$pop_by4age2,,t]),
  
  (  sims[idx$inc_year_gi3_3,,t]+
       sims[idx$inc_year_gi_3,,t]+
       sims[idx$inc_year_gii4_3,,t]+
       sims[idx$inc_year_gii_3,,t])/(sims[idx$pop_by4age3,,t]) ,
  
  (  sims[idx$inc_year_gi3_4,,t]+
       sims[idx$inc_year_gi_4,,t]+
       sims[idx$inc_year_gii4_4,,t]+
       sims[idx$inc_year_gii_4,,t])/(sims[idx$pop_by4age4,,t]) ,
  
  (  sims[idx$inc_year_gi3_5,,t]+
       sims[idx$inc_year_gi_5,,t]+
       sims[idx$inc_year_gii4_5,,t]+
       sims[idx$inc_year_gii_5,,t])/(sims[idx$pop_by4age5,,t])))

irates<-1000*t(irates)

df_qtls <- as.data.frame(rowQuantiles((irates),
                                      probs = c(0.025, 0.5, 0.975)))

df1 <- data.frame(t(irates)) 
colnames(df1) <- paste(c("0_1","1_4","5_14","15_64","65+"))
df_m <- reshape2::melt(df1)
df_m$variable <- as.factor(df_m$variable)
df_d$x <- factor(c("0_1","1_4","5_14","15_64","65+"))
df_qtls$x <- factor(c("0_1","1_4","5_14","15_64","65+"))


viol_col <- "#ff4d4d"# "chartreuse3" # "yellow3"
err_col <- "black"
data_col <- "black"

fits_iid2 <- ggplot() +
  geom_violin(
    data = df_m,
    aes(x = variable, y = value, fill = "Posterior Density"),
    draw_quantiles = c(0.5),
    width = 1,
    linetype = 1,
    trim = FALSE,
    color = "white",
    alpha = 0.7
  ) +
  geom_point(data = df_d, mapping = aes(x = x, y = inc, color = "Data (95% CI)"), size = 2, shape = 15) +
  geom_errorbar(
    mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
    width = .2, position = position_dodge(.9)
  ) +
  labs(title = "Model vs IID2 community incidence", x = "Age group (years)", y = "Incidence per 1000 \n person-year") +
  theme_minimal() +
  ylim(0, 600) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


if (save_plot==1){
  save(fits_iid2, file = here("output","fits_iid2.rdata"))
}


######## SGSS cases reported 

id<- seq(2,length(sims[idx$reported_wk_1,1,]))
reported<-
  sims[idx$reported_wk_1,,id] +
  sims[idx$reported_wk_2,,id] +
  sims[idx$reported_wk_3,,id] +
  sims[idx$reported_wk_4,,id] +
  sims[idx$reported_wk_5,,id] +
  sims[idx$reported_wk_6,,id] +
  sims[idx$reported_wk_7,,id] +
  sims[idx$reported_wk_8,,id] +
  sims[idx$reported_wk_9,,id] +
  sims[idx$reported_wk_10,,id] +
  sims[idx$reported_wk_11,,id] +
  sims[idx$reported_wk_12,,id] +
  sims[idx$reported_wk_13,,id] +
  sims[idx$reported_wk_14,,id]
  

reported_obs<-(data$reported[which(!is.na(data$reported))])


mo <- as.Date(seq(as.Date(sgss_start),
                  by = "week",
                  length.out = dim(reported)[2]
))

df_s <- as.data.frame(
  rowQuantiles(t(reported),
               probs = c(0.025, 0.5, 0.975)
  )
)
df_s$x <- days_vec

df_d <- data.frame(
  x = total_cases$date,
  cases = reported_obs
)

df_sim <- data.frame(
  x = days_vec,
  t(reported)
)

dat_sim <- reshape2::melt(df_sim, id = "x")


fits_sgss <- ggplot(data = df_s, aes(x = x)) +
  # geom_line(
  #   data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
  #   alpha = 0.2, lwd = 0.4
  # ) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
  geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
  geom_point(data = df_d, aes(x = x, y = cases)) +
  labs(title = "Model vs SGSS reported cases", x = " ", y = "Cases reported\n (per week)") +
  theme_minimal() +
  scale_x_date(date_breaks = "4 month", date_labels = "%b-%Y",
               limits = as.Date(c('2014-01-01','2020-12-01'))) +
  ylim(c(0,5000))+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


if (save_plot==1){
  save(fits_sgss, file = here("output", "fits_sgss.rdata"))
}

######## Lindesmith seroprevallence children 
x_d <- c(1, 3, 5, 7, 9, 11) # bin x axis positions

id<-which(!is.na(data$sero1))

sero_obs<-c(data$sero1[id], data$sero2[id], data$sero3[id], data$sero4[id],
            data$sero5[id], data$sero6[id])*100

df_d <- data.frame(
  x = x_d,
  sero = sero_obs,
  low = sero_obs,
  up = sero_obs
)

fac<-c(2,3,4,5,6,7)

id<-which(sims[idx$t,1,]%in%
            data$time_end[which(!is.na(data$sero1))])

sero_model<-rbind(
  sims[idx$seroprev1.2,,id] ,
  sims[idx$seroprev2.3,,id] ,
  sims[idx$seroprev3.4,,id] ,
  sims[idx$seroprev4.5,,id] ,
  sims[idx$seroprev5.6,,id] ,
  sims[idx$seroprev6.7,,id])*100



df_qtls <- as.data.frame(rowQuantiles((sero_model),
                                      probs = c(0.025, 0.5, 0.975)))

df1 <- data.frame(t(sero_model)) 
colnames(df1) <- paste(c("0_1","1_2","2_3","3_4","5_6","6_7"))
df_m <- reshape2::melt(df1)
df_m$variable <- as.factor(df_m$variable)
df_d$x <- factor(c("0_1","1_2","2_3","3_4","5_6","6_7"))
df_qtls$x <- factor(c("0_1","1_2","2_3","3_4","5_6","6_7"))


viol_col <-  "yellow3"
err_col <- "black"
data_col <- "black"

fits_sero <- ggplot() +
  geom_violin(
    data = df_m,
    aes(x = variable, y = value, fill = "Posterior Density"),
    draw_quantiles = c(0.5),
    width = 0.8,
    linetype = 1,
    trim = FALSE,
    color = viol_col,
    alpha = 0.5
  ) +
  geom_point(data = df_d, mapping = aes(x = x, y = sero, 
                                        color = "Data (95% CI)"), 
             size = 2, shape = 15) +
  # geom_errorbar(
  #   mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
  #   width = .2, position = position_dodge(.9)
  # ) +
  labs(title = "Model vs Seroprevalence of GII.4", x = "Age group (years)", y = "Seroprevalence (%)") +
  theme_minimal() +
  ylim(0, 100) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


if (save_plot==1){
  save(fits_sero, file = here("output","fits_sero.rdata"))
}




# Genotype by age  --------------------------------------------------------


id<-which(sims[idx$t,1,]%in%
            data$time_end[which(!is.na(data$gi_prop_1))])


denom_1<-( sims[idx$inc_year_gi3_1,,id] +
             sims[idx$inc_year_gi_1,,id] +
             sims[idx$inc_year_gii4_1,,id] +
             sims[idx$inc_year_gii_1,,id] )+
  ( sims[idx$inc_year_gi3_2,,id] +
      sims[idx$inc_year_gi_2,,id] +
      sims[idx$inc_year_gii4_2,,id] +
      sims[idx$inc_year_gii_2,,id] )

denom_2<-( sims[idx$inc_year_gi3_3,,id] +
             sims[idx$inc_year_gi_3,,id] +
             sims[idx$inc_year_gii4_3,,id] +
             sims[idx$inc_year_gii_3,,id] )

denom_3<-( sims[idx$inc_year_gi3_4,,id] +
             sims[idx$inc_year_gi_4,,id] +
             sims[idx$inc_year_gii4_4,,id] +
             sims[idx$inc_year_gii_4,,id] )

denom_4<-( sims[idx$inc_year_gi3_5,,id] +
             sims[idx$inc_year_gi_5,,id] +
             sims[idx$inc_year_gii4_5,,id] +
             sims[idx$inc_year_gii_5,,id] )

modelled_gi3_1 <- (( sims[idx$inc_year_gi3_1,,id] +  sims[idx$inc_year_gi3_2,,id] )/denom_1)
modelled_gi_1  <- (( sims[idx$inc_year_gi_1,,id] +  sims[idx$inc_year_gi_2,,id] )/denom_1)
modelled_gii4_1<- (( sims[idx$inc_year_gii4_1,,id] +  sims[idx$inc_year_gii4_2,,id] )/denom_1)
modelled_gii_1 <- (( sims[idx$inc_year_gii_1,,id] +  sims[idx$inc_year_gii_2,,id] )/denom_1)

modelled_gi3_2 <- (( sims[idx$inc_year_gi3_3,,id] )/denom_2)
modelled_gi_2  <- (( sims[idx$inc_year_gi_3,,id] )/denom_2)
modelled_gii4_2<- (( sims[idx$inc_year_gii4_3,,id] )/denom_2)
modelled_gii_2 <- (( sims[idx$inc_year_gii_3,,id] )/denom_2)

modelled_gi3_3 <- (( sims[idx$inc_year_gi3_4,,id] )/denom_3)
modelled_gi_3  <- (( sims[idx$inc_year_gi_4,,id] )/denom_3)
modelled_gii4_3<- (( sims[idx$inc_year_gii4_4,,id] )/denom_3)
modelled_gii_3 <- (( sims[idx$inc_year_gii_4,,id] )/denom_3)

modelled_gi3_4 <- (( sims[idx$inc_year_gi3_5,,id] )/denom_4)
modelled_gi_4  <- (( sims[idx$inc_year_gi_5,,id] )/denom_4)
modelled_gii4_4<- (( sims[idx$inc_year_gii4_5,,id] )/denom_4)
modelled_gii_4 <- (( sims[idx$inc_year_gii_5,,id] )/denom_4)


strain_model<-cbind(
  modelled_gi3_1,
  modelled_gi_1  ,
  modelled_gii4_1,
  modelled_gii_1 ,
  modelled_gi3_2 ,
  modelled_gi_2  ,
  modelled_gii4_2,
  modelled_gii_2 ,
  
  modelled_gi3_3 ,
  modelled_gi_3  ,
  modelled_gii4_3,
  modelled_gii_3 ,
  
  modelled_gi3_4 ,
  modelled_gi_4 ,
  modelled_gii4_4,
  modelled_gii_4 
)*100

data_df<- read.csv(here("data","IID2genotype_by_age.csv"), header=TRUE)#, sep=,)
head(data_df)
data_df$prop[1]<-0.01

data_df$Age<-factor(data_df$Age, levels = c("0-4", "5 to 15", "15 to 64", "65+"))
data_df$Genotype<-factor(data_df$genotype, levels = c("GI3", "GI ", "GII4", "GII "))
data_df$y<-data_df$prop*100

df_qtls <- as.data.frame(rowQuantiles(t(strain_model),
                                      probs = c(0.025, 0.5, 0.975)))

df_qtls$Age<-c("0-4","0-4","0-4","0-4",
               "5 to 15","5 to 15","5 to 15","5 to 15",
               "15 to 64","15 to 64","15 to 64","15 to 64",
               "65+","65+","65+","65+")

df_qtls$Age<-factor(df_qtls$Age, levels = c("0-4", "5 to 15", "15 to 64", "65+"))

df_qtls$Genotype<- c("GI3", "GI ", "GII4", "GII ", "GI3", "GI ", "GII4", "GII ","GI3", "GI ", "GII4", "GII ","GI3", "GI ", "GII4", "GII ") 
df_qtls$Genotype<- factor(df_qtls$Genotype, levels = c("GI3", "GI ", "GII4", "GII "))


df<- data_df %>% 
  left_join(df_qtls,by = c('Age','Genotype'))


fitagegeno<-ggplot(df,aes(x=Age, y=y, fill=Genotype ))+
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar( mapping = aes (x=Age, ymin=`2.5%`, ymax=`97.5%`), width=.2, position=position_dodge(.9))+
  geom_point(aes(x = Age, y = `50%`, color = "Mean"), width=.2,position=position_dodge(.9))+
  labs(title = "Model fit to Genotype Distribution by Age group in IID2", y='%')+
  theme_classic() +
  ylim(0, 100) +
  scale_fill_manual(values=c('skyblue3','firebrick2','yellow3','grey40'))+
  scale_color_manual(name = "Model", values = "black") +
  guides(fill = guide_legend(override.aes = list(shape = NA)))+
  theme(
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0)
  )



# 
# gi3_model<-((sims[idx$cumm_incday_gi3_1,,id]+
#                sims[idx$cumm_incday_gi3_2,,id]+
#                sims[idx$cumm_incday_gi3_3,,id]+
#                sims[idx$cumm_incday_gi3_4,,id]) /
#               (    
#                 sims[idx$cumm_incday_gii4_1,,id]+
#                   sims[idx$cumm_incday_gii4_2,,id]+
#                   sims[idx$cumm_incday_gii4_3,,id]+
#                   sims[idx$cumm_incday_gii4_4,,id]+
#                   sims[idx$cumm_incday_gii_1,,id]+
#                   sims[idx$cumm_incday_gii_2,,id]+
#                   sims[idx$cumm_incday_gii_3,,id]+
#                   sims[idx$cumm_incday_gii_4,,id]+
#                   sims[idx$cumm_incday_gi3_1,,id]+
#                   sims[idx$cumm_incday_gi3_2,,id]+
#                   sims[idx$cumm_incday_gi3_3,,id]+
#                   sims[idx$cumm_incday_gi3_4,,id]+
#                   sims[idx$cumm_incday_gi_1,,id]+
#                   sims[idx$cumm_incday_gi_2,,id]+
#                   sims[idx$cumm_incday_gi_3,,id]+
#                   sims[idx$cumm_incday_gi_4,,id]))
# 
# gi_model<-((sims[idx$cumm_incday_gi_1,,id]+
#               sims[idx$cumm_incday_gi_2,,id]+
#               sims[idx$cumm_incday_gi_3,,id]+
#               sims[idx$cumm_incday_gi_4,,id]) /
#              (    
#                sims[idx$cumm_incday_gii4_1,,id]+
#                  sims[idx$cumm_incday_gii4_2,,id]+
#                  sims[idx$cumm_incday_gii4_3,,id]+
#                  sims[idx$cumm_incday_gii4_4,,id]+
#                  sims[idx$cumm_incday_gii_1,,id]+
#                  sims[idx$cumm_incday_gii_2,,id]+
#                  sims[idx$cumm_incday_gii_3,,id]+
#                  sims[idx$cumm_incday_gii_4,,id]+
#                  sims[idx$cumm_incday_gi3_1,,id]+
#                  sims[idx$cumm_incday_gi3_2,,id]+
#                  sims[idx$cumm_incday_gi3_3,,id]+
#                  sims[idx$cumm_incday_gi3_4,,id]+
#                  sims[idx$cumm_incday_gi_1,,id]+
#                  sims[idx$cumm_incday_gi_2,,id]+
#                  sims[idx$cumm_incday_gi_3,,id]+
#                  sims[idx$cumm_incday_gi_4,,id]))
# 
# gii4_model<-((sims[idx$cumm_incday_gii4_1,,id]+
#                 sims[idx$cumm_incday_gii4_2,,id]+
#                 sims[idx$cumm_incday_gii4_3,,id]+
#                 sims[idx$cumm_incday_gii4_4,,id]) /
#                (    
#                  sims[idx$cumm_incday_gii4_1,,id]+
#                    sims[idx$cumm_incday_gii4_2,,id]+
#                    sims[idx$cumm_incday_gii4_3,,id]+
#                    sims[idx$cumm_incday_gii4_4,,id]+
#                    sims[idx$cumm_incday_gii_1,,id]+
#                    sims[idx$cumm_incday_gii_2,,id]+
#                    sims[idx$cumm_incday_gii_3,,id]+
#                    sims[idx$cumm_incday_gii_4,,id]+
#                    sims[idx$cumm_incday_gi3_1,,id]+
#                    sims[idx$cumm_incday_gi3_2,,id]+
#                    sims[idx$cumm_incday_gi3_3,,id]+
#                    sims[idx$cumm_incday_gi3_4,,id]+
#                    sims[idx$cumm_incday_gi_1,,id]+
#                    sims[idx$cumm_incday_gi_2,,id]+
#                    sims[idx$cumm_incday_gi_3,,id]+
#                    sims[idx$cumm_incday_gi_4,,id]))
# 
# gii_model<-((sims[idx$cumm_incday_gii_1,,id]+
#                sims[idx$cumm_incday_gii_2,,id]+
#                sims[idx$cumm_incday_gii_3,,id]+
#                sims[idx$cumm_incday_gii_4,,id]) /
#               (    
#                 sims[idx$cumm_incday_gii4_1,,id]+
#                   sims[idx$cumm_incday_gii4_2,,id]+
#                   sims[idx$cumm_incday_gii4_3,,id]+
#                   sims[idx$cumm_incday_gii4_4,,id]+
#                   sims[idx$cumm_incday_gii_1,,id]+
#                   sims[idx$cumm_incday_gii_2,,id]+
#                   sims[idx$cumm_incday_gii_3,,id]+
#                   sims[idx$cumm_incday_gii_4,,id]+
#                   sims[idx$cumm_incday_gi3_1,,id]+
#                   sims[idx$cumm_incday_gi3_2,,id]+
#                   sims[idx$cumm_incday_gi3_3,,id]+
#                   sims[idx$cumm_incday_gi3_4,,id]+
#                   sims[idx$cumm_incday_gi_1,,id]+
#                   sims[idx$cumm_incday_gi_2,,id]+
#                   sims[idx$cumm_incday_gi_3,,id]+
#                   sims[idx$cumm_incday_gi_4,,id]))
# 
# 
# 
# 
# 
# 
# if (fake==1){
#   
#   g2<-rnorm(n = nsim, mean = 88, sd = 3 )
#   g2[g2>100]<-100
#   g1<-rnorm(n = nsim, mean = 100-88, sd = 3 )
#   g1[g1<0]<-0
#   
#   strain_model<-cbind(g1, g2)
# } else{
#   strain_model<-cbind(gi3_model,gi_model, gii4_model,gii_model )*100
# }
# 
# 
# id<-which(!is.na(data$gi_prop))
# strain_obs<-c(data$gi3_prop[id],
#               data$gi_prop[id],
#               data$gii4_prop[id],
#               data$gii_prop[id])*100
# 
# x_d <- c(1,3,5,7) # bin x axis positions
# 
# 
# df_d <- data.frame(
#   x = x_d,
#   strain = strain_obs,
#   low = strain_obs,
#   up = strain_obs
# )
# 
# df_qtls <- as.data.frame(rowQuantiles(t(strain_model),
#                                       probs = c(0.025, 0.5, 0.975)))
# 
# df1 <- data.frame((strain_model)) 
# colnames(df1) <- paste(c("GI.3","Other GI","GII.4", "Other GII"))
# df_m <- reshape2::melt(df1)
# df_m$variable <- as.factor(df_m$variable)
# df_d$x <- factor(c("GI.3","Other GI","GII.4", "Other GII"))
# df_qtls$x <- factor(c("GI.3","Other GI","GII.4", "Other GII"))
# 
# 
# viol_col <-  "orange"
# err_col <- "black"
# data_col <- "black"
# 
# fits_strain <- ggplot() +
#   geom_violin(
#     data = df_m,
#     aes(x = variable, y = value, fill = variable),
#     draw_quantiles = c(0.5),
#     width = 0.8,
#     linetype = 1,
#     color = NA,
#     trim = FALSE,
#     alpha = 0.5
#   ) +
#   geom_point(data = df_d, mapping = aes(x = x, y = strain), 
#                                         color = 'black', 
#              size = 2, shape = 15) +
#   geom_errorbar(
#     mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
#     width = .2, position = position_dodge(.9)
#   ) +
#   labs(title = "Model vs Genotype type", x = "Genotype", y = "Proportion (%)") +
#   theme_minimal() +
#   ylim(0, 100) +
#   scale_fill_manual(values=c('skyblue3','firebrick2','yellow3','grey28'))+
#  
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),
#     axis.text = element_text(colour = "black", size = 12, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11), legend.key = element_blank(),
#     axis.text.x = element_text(angle = 60, hjust = 1)
#   )
# 
# 
# if (save_plot==1){
#   save(fits_strain, file = here("output","fits_strain.rdata"))
# }
# 

## Age 
#################

id<-which(sims[idx$t,1,]%in%
            data$time_end[which(!is.na(data$a0_prop))])


a0_model<-((sims[idx$cumm_incday_gi3_1,,id]+
              sims[idx$cumm_incday_gi_1,,id]+
              sims[idx$cumm_incday_gii4_1,,id]+
              sims[idx$cumm_incday_gii_1,,id]) /
             (    
               sims[idx$cumm_incday_gii4_1,,id]+
                 sims[idx$cumm_incday_gii4_2,,id]+
                 sims[idx$cumm_incday_gii4_3,,id]+
                 sims[idx$cumm_incday_gii4_4,,id]+
                 sims[idx$cumm_incday_gii_1,,id]+
                 sims[idx$cumm_incday_gii_2,,id]+
                 sims[idx$cumm_incday_gii_3,,id]+
                 sims[idx$cumm_incday_gii_4,,id]+
                 sims[idx$cumm_incday_gi3_1,,id]+
                 sims[idx$cumm_incday_gi3_2,,id]+
                 sims[idx$cumm_incday_gi3_3,,id]+
                 sims[idx$cumm_incday_gi3_4,,id]+
                 sims[idx$cumm_incday_gi_1,,id]+
                 sims[idx$cumm_incday_gi_2,,id]+
                 sims[idx$cumm_incday_gi_3,,id]+
                 sims[idx$cumm_incday_gi_4,,id]))


a5_model<-((sims[idx$cumm_incday_gi3_2,,id]+
              sims[idx$cumm_incday_gi_2,,id]+
              sims[idx$cumm_incday_gii4_2,,id]+
              sims[idx$cumm_incday_gii_2,,id]) /
             (    
               sims[idx$cumm_incday_gii4_1,,id]+
                 sims[idx$cumm_incday_gii4_2,,id]+
                 sims[idx$cumm_incday_gii4_3,,id]+
                 sims[idx$cumm_incday_gii4_4,,id]+
                 sims[idx$cumm_incday_gii_1,,id]+
                 sims[idx$cumm_incday_gii_2,,id]+
                 sims[idx$cumm_incday_gii_3,,id]+
                 sims[idx$cumm_incday_gii_4,,id]+
                 sims[idx$cumm_incday_gi3_1,,id]+
                 sims[idx$cumm_incday_gi3_2,,id]+
                 sims[idx$cumm_incday_gi3_3,,id]+
                 sims[idx$cumm_incday_gi3_4,,id]+
                 sims[idx$cumm_incday_gi_1,,id]+
                 sims[idx$cumm_incday_gi_2,,id]+
                 sims[idx$cumm_incday_gi_3,,id]+
                 sims[idx$cumm_incday_gi_4,,id]))

a15_model<-((sims[idx$cumm_incday_gi3_3,,id]+
               sims[idx$cumm_incday_gi_3,,id]+
               sims[idx$cumm_incday_gii4_3,,id]+
               sims[idx$cumm_incday_gii_3,,id]) /
              (    
                sims[idx$cumm_incday_gii4_1,,id]+
                  sims[idx$cumm_incday_gii4_2,,id]+
                  sims[idx$cumm_incday_gii4_3,,id]+
                  sims[idx$cumm_incday_gii4_4,,id]+
                  sims[idx$cumm_incday_gii_1,,id]+
                  sims[idx$cumm_incday_gii_2,,id]+
                  sims[idx$cumm_incday_gii_3,,id]+
                  sims[idx$cumm_incday_gii_4,,id]+
                  sims[idx$cumm_incday_gi3_1,,id]+
                  sims[idx$cumm_incday_gi3_2,,id]+
                  sims[idx$cumm_incday_gi3_3,,id]+
                  sims[idx$cumm_incday_gi3_4,,id]+
                  sims[idx$cumm_incday_gi_1,,id]+
                  sims[idx$cumm_incday_gi_2,,id]+
                  sims[idx$cumm_incday_gi_3,,id]+
                  sims[idx$cumm_incday_gi_4,,id]))

a65_model<-((  sims[idx$cumm_incday_gi3_4,,id]+
                 sims[idx$cumm_incday_gi_4,,id]+
                 sims[idx$cumm_incday_gii4_4,,id]+
                 sims[idx$cumm_incday_gii_4,,id]) /
              (    
                sims[idx$cumm_incday_gii4_1,,id]+
                  sims[idx$cumm_incday_gii4_2,,id]+
                  sims[idx$cumm_incday_gii4_3,,id]+
                  sims[idx$cumm_incday_gii4_4,,id]+
                  sims[idx$cumm_incday_gii_1,,id]+
                  sims[idx$cumm_incday_gii_2,,id]+
                  sims[idx$cumm_incday_gii_3,,id]+
                  sims[idx$cumm_incday_gii_4,,id]+
                  sims[idx$cumm_incday_gi3_1,,id]+
                  sims[idx$cumm_incday_gi3_2,,id]+
                  sims[idx$cumm_incday_gi3_3,,id]+
                  sims[idx$cumm_incday_gi3_4,,id]+
                  sims[idx$cumm_incday_gi_1,,id]+
                  sims[idx$cumm_incday_gi_2,,id]+
                  sims[idx$cumm_incday_gi_3,,id]+
                  sims[idx$cumm_incday_gi_4,,id]))






if (fake==1){
  
  g2<-rnorm(n = nsim, mean = 88, sd = 3 )
  g2[g2>100]<-100
  g1<-rnorm(n = nsim, mean = 100-88, sd = 3 )
  g1[g1<0]<-0
  
  age_model<-cbind(g1, g2)
} else{
  age_model<-cbind(a0_model,a5_model, a15_model,a65_model )*100
}


id<-which(!is.na(data$a0_prop))
age_obs<-c(data$a0_prop[id],
              data$a5_prop[id],
              data$a15_prop[id],
              data$a65_prop[id])*100

x_d <- c(1,3,5,7) # bin x axis positions


df_d <- data.frame(
  x = x_d,
  age = age_obs,
  low = age_obs,
  up = age_obs
)

df_qtls <- as.data.frame(rowQuantiles(t(age_model),
                                      probs = c(0.025, 0.5, 0.975)))

df1 <- data.frame((age_model)) 
colnames(df1) <- paste(c("0-4","5-14","15-64", "65+"))
df_m <- reshape2::melt(df1)
df_m$variable <- as.factor(df_m$variable)
df_d$x <- factor(c("0-4","5-14","15-64", "65+"))
df_qtls$x <- factor(c("0-4","5-14","15-64", "65+"))


viol_col <-  "orange"
err_col <- "black"
data_col <- "black"

fits_age <- ggplot() +
  geom_violin(
    data = df_m,
    aes(x = variable, y = value, fill = "Posterior Density"),
    draw_quantiles = c(0.5),
    width = 0.8,
    linetype = 1,
    trim = FALSE,
    color = viol_col,
    alpha = 0.5
  ) +
  geom_point(data = df_d, mapping = aes(x = x, y = age, 
                                        color = "Data (95% CI)"), 
             size = 2, shape = 15) +
  # geom_errorbar(
  #   mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
  #   width = .2, position = position_dodge(.9)
  # ) +
  labs(title = "Reported by Age (2014-2023)", x = "Age", y = "Proportion (%)") +
  theme_minimal() +
  ylim(0, 100) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


if (save_plot==1){
  save(fits_age, file = here("output","fits_age.rdata"))
}






##
gridExtra::grid.arrange(fits_iid2)
gridExtra::grid.arrange(fits_sero)
gridExtra::grid.arrange(fits_sgss)
gridExtra::grid.arrange(fitagegeno)
gridExtra::grid.arrange(fits_age)

