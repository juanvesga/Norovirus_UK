library(matrixStats)
## Plot model fits
save_plot=0
nsim <- 100
posteriors<-processed_chains$pars
# Sample
pars <- posteriors[sample(nrow(posteriors), nsim), ]

inc_GI<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
inc_GII<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
inc_GI3<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
inc_OGI<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
inc_GII4<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
inc_OGII<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)

cases1_1<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases4_1<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases14_1<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases64_1<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases65p_1<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases1_2<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases4_2<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases14_2<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases64_2<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)
cases65p_2<-matrix(0,ncol = nsim, nrow = length(days_vec)+1)




for (h in 1:nsim){
print(100*(h/nsim))
theta=c(pars[h,]/unlist(params$scaling_fac))



p<-list(
  beta_1  = theta[['beta_1']],   # transm coefficient
  beta_2  = theta[['beta_2']],   # transm coefficient
  beta_3  = theta[['beta_3']],   # transm coefficient
  beta_4  = theta[['beta_4']],   # transm coefficient
  maternalAB  = theta[['maternalAB']],   
  imm_yr = theta[['imm_yr']],
  imm_fac = theta[['imm_fac']],
  w1_1 = theta[['w1_1']],
  w2 = theta[['w2']],
  repfac=theta[['repfac']],
  aduRR=theta[["aduRR"]],
  repfac_65p=theta[['repfac_65p']],
  crossp_12=theta[['crossp_12']],
  crossp_21=theta[['crossp_21']],
  crossp_34=theta[['crossp_34']],
  crossp_43=theta[['crossp_43']],
  pop  = pop,
  init  = init,
  mu    = params$mu,
  m=params$transmission,
  m_holi=params$transmission_holi,
  cmx_1=params$cmx_1,
  cmx_2=params$cmx_2,
  cmx_3=params$cmx_3,
  cmx_4=params$cmx_4,
  cmx_5=params$cmx_5,
  cmx_6=params$cmx_6,
  cmx_7=params$cmx_7,
  cmx_8=params$cmx_8,
  cmx_9=params$cmx_9,
  aging_vec =params$aging_vec,
  school_step= as.double(params$school_uk),
  n_school_steps=params$n_school_steps,
  covid_step= as.double(params$covid_sche),
  n_covid_steps=params$n_covid_steps,
  N_age = params$N_age,
  vac_camp_cov=params$vac_camp_cov,
  vac_sche_cov=params$vac_sche_cov,
  dt=1
) # number of age groups


filter2$run(p, save_history = TRUE)
sims<-drop(filter2$history())
state<-filter2$state()

inc_GI[,h]<- (sims['inc_day_gi3',] + sims['inc_day_gi',])
inc_GII[,h]<-(sims['inc_day_gii4',] + sims['inc_day_gii',])
inc_GI[,h]<- (sims['inc_day_gi3',] + sims['inc_day_gi',])
inc_GII[,h]<-(sims['inc_day_gii4',] + sims['inc_day_gii',])

cases1_1[,h]<-sims['inc_day_gi3_1',] + sims['inc_day_gi_1',]
cases4_1[,h]<-sims['inc_day_gi3_2',] + sims['inc_day_gi_2',]
cases14_1[,h]<-sims['inc_day_gi3_3',] + sims['inc_day_gi_3',]
cases64_1[,h]<-sims['inc_day_gi3_4',] + sims['inc_day_gi_4',]
cases65p_1[,h]<-sims['inc_day_gi3_5',] + sims['inc_day_gi_5',]
cases1_2[,h] <-sims['inc_day_gii4_1',] + sims['inc_day_gii_1',]
cases4_2[,h] <-sims['inc_day_gii4_2',] + sims['inc_day_gii_2',]
cases14_2[,h]<-sims['inc_day_gii4_3',] + sims['inc_day_gii_3',]
cases64_2[,h]<-sims['inc_day_gii4_4',] + sims['inc_day_gii_4',]
cases65p_2[,h]<-sims['inc_day_gii4_5',] + sims['inc_day_gii_5',]


}



########  GI
GIcol<-"#69b3a2"

x <- days_vec
y<-inc_GI
y<- y[-1,]

df_s <- as.data.frame(
  rowQuantiles(y,
               probs = c(0.25, 0.5, 0.75)
  )
)
df_s$x <- x


df_sim <- data.frame(
  x = x,
  y=y
)

dat_sim <- reshape2::melt(df_sim, id = "x")


incGI <- ggplot(data = df_s, aes(x = x)) +
  geom_line(
    data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
    alpha = 0.2, lwd = 0.4
  ) +
  geom_ribbon(aes(ymin = `25%`, ymax = `75%`), fill = GIcol, alpha = 0.2) +
  geom_line(aes(y = `50%`), col = GIcol, lwd = 1) +
  labs(title = "GI incidence", x = " ", y = "Incidence (symptomatic)\n (per wk)") +
  theme_minimal() +
  scale_x_date(limits = c(as.Date("2012-01-01"),as.Date("2020-01-01") )) +
  ylim(0,1e5)+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

if (save_plot==1){
  save(
    incGI, file = here("output", "inc_GI.rdata"))
}






########  GII

GIIcol <- "#ff4d4d"

x <- days_vec
y<-inc_GII
y<- y[-1,]

df_s <- as.data.frame(
  rowQuantiles(y,
               probs = c(0.25, 0.5, 0.75)
  )
)
df_s$x <- x


df_sim <- data.frame(
  x = x,
  y=y
)

dat_sim <- reshape2::melt(df_sim, id = "x")


incGII <- ggplot(data = df_s, aes(x = x)) +
  geom_line(
    data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
    alpha = 0.2, lwd = 0.4
  ) +
  geom_ribbon(aes(ymin = `25%`, ymax = `75%`), fill = GIIcol, alpha = 0.2) +
  geom_line(aes(y = `50%`), col = GIIcol, lwd = 1) +
  labs(title = "GII incidence", x = " ", y = "Incidence (symptomatic)\n (per wk)") +
  theme_minimal() +
  scale_x_date(limits = c(as.Date("2010-01-01"),as.Date("2020-01-01") )) +
  ylim(0,1e5)+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )
windows()
gridExtra::grid.arrange(incGI,incGII)
if (save_plot==1){
  save(
    incGII, file = here("output", "inc_GII.rdata"))
}


# Stacked incidence
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
mo <- as.Date(seq(as.Date(sim_startdate),as.Date(sim_enddate),
                  by = "year")) 

id<-which(days_vec %in% mo)
mo<-days_vec

c1<-cases1_1[-1,]
c2<-cases4_1[-1,]
c3<-cases14_1[-1,]
c4<-cases64_1[-1,]
c5<-cases65p_1[-1,]


df_1 <- as.data.frame(rowQuantiles(c1,probs = c(0.025, 0.5, 0.975)))
df_1$x <- days_vec

df_2 <- as.data.frame(rowQuantiles(c2,probs = c(0.025, 0.5, 0.975)))
df_2$x <- days_vec

df_3 <- as.data.frame(rowQuantiles(c3, probs = c(0.025, 0.5, 0.975)))
df_3$x <- days_vec

df_4 <- as.data.frame(rowQuantiles(c4,probs = c(0.025, 0.5, 0.975)))
df_4$x <- days_vec

df_5 <- as.data.frame(rowQuantiles(c5,probs = c(0.025, 0.5, 0.975)))
df_5$x <- days_vec


tmp1 <- data.frame(date_index = mo, count = round(df_1$`50%`), incident = "0_1", id = 1)
tmp2 <- data.frame(date_index = mo, count = round(df_2$`50%`), incident = "1_4", id = 2)
tmp3 <- data.frame(date_index = mo, count = round(df_3$`50%`), incident = "5_14", id = 3)
tmp4 <- data.frame(date_index = mo, count = round(df_4$`50%`), incident = "15_65", id = 4)
tmp5 <- data.frame(date_index = mo, count = round(df_5$`50%`), incident = "65p", id = 5)

linelist <- rbind(tmp1, tmp2, tmp3, tmp4,tmp5)
linelist$Outcome <- factor(linelist$incident, levels = c("0_1", "1_4", "5_14", "15_65", "65p"))
colors <- c("0_1" = "brown1", "1_4" = "gold", "5_14" = "lightseagreen",
            "15_65" = "grey43", "65p"="orange")

linelistg2<-linelist

stack1 <- ggplot(linelist, aes(
  fill = Outcome,
  y = count, x = date_index
)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Incidence GI", x = " ", y = "Incidence\n (symptomatic) per wk") +
  theme_classic() +
  scale_x_date(date_breaks = "3 month", 
               date_labels = "%b-%Y", 
               limits = as.Date(c("2017-06-01", "2018-05-30"))) +
  theme( # legend.position = "bottom",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  guides(col = guide_legend(title = NULL))

windows()
stack
save(stack, file = here("output",  "stack_incidenceGII.rdata"))


# Stacked incidence
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
mo <- as.Date(seq(as.Date(sim_startdate),as.Date(sim_enddate),
                  by = "year")) 

id<-which(days_vec %in% mo)
mo<-days_vec

c1<-cases1_2[-1,]
c2<-cases4_2[-1,]
c3<-cases14_2[-1,]
c4<-cases64_2[-1,]
c5<-cases65p_2[-1,]


df_1 <- as.data.frame(rowQuantiles(c1,probs = c(0.025, 0.5, 0.975)))
df_1$x <- days_vec

df_2 <- as.data.frame(rowQuantiles(c2,probs = c(0.025, 0.5, 0.975)))
df_2$x <- days_vec

df_3 <- as.data.frame(rowQuantiles(c3, probs = c(0.025, 0.5, 0.975)))
df_3$x <- days_vec

df_4 <- as.data.frame(rowQuantiles(c4,probs = c(0.025, 0.5, 0.975)))
df_4$x <- days_vec

df_5 <- as.data.frame(rowQuantiles(c5,probs = c(0.025, 0.5, 0.975)))
df_5$x <- days_vec


tmp1 <- data.frame(date_index = mo, count = round(df_1$`50%`), incident = "0_1", id = 1)
tmp2 <- data.frame(date_index = mo, count = round(df_2$`50%`), incident = "1_4", id = 2)
tmp3 <- data.frame(date_index = mo, count = round(df_3$`50%`), incident = "5_14", id = 3)
tmp4 <- data.frame(date_index = mo, count = round(df_4$`50%`), incident = "15_65", id = 4)
tmp5 <- data.frame(date_index = mo, count = round(df_5$`50%`), incident = "65p", id = 5)

linelist <- rbind(tmp1, tmp2, tmp3, tmp4,tmp5)
linelist$Outcome <- factor(linelist$incident, levels = c("0_1", "1_4", "5_14", "15_65", "65p"))
colors <- c("0_1" = "brown1", "1_4" = "gold", "5_14" = "lightseagreen",
            "15_65" = "grey43", "65p"="orange")

linelistg1<-linelist

stack2 <- ggplot(linelist, aes(
  fill = Outcome,
  y = count, x = date_index
)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Incidence GII", x = " ", y = "Incidence\n (symptomatic) per wk") +
  theme_classic() +
  scale_x_date(date_breaks = "3 month", 
               date_labels = "%b-%Y", 
               limits = as.Date(c("2017-06-01", "2018-05-30"))) +
  theme( # legend.position = "bottom",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  guides(col = guide_legend(title = NULL))

windows()
stack2
save(stack, file = here("output",  "stack_incidenceGII.rdata"))

windows()
gridExtra::grid.arrange(stack1,stack2)



## All Inc stacked

mo <- as.Date(seq(as.Date(sim_startdate),as.Date(sim_enddate),
                  by = "year")) 

id<-which(days_vec %in% mo)
mo<-days_vec

c1<-cases1_1[-1,]+cases1_2[-1,]
c2<-cases4_1[-1,]+cases4_2[-1,]
c3<-cases14_1[-1,]+cases14_2[-1,]
c4<-cases64_1[-1,]+cases64_2[-1,]
c5<-cases65p_1[-1,]+cases65p_2[-1,]


df_1 <- as.data.frame(rowQuantiles(c1,probs = c(0.025, 0.5, 0.975)))
df_1$x <- days_vec

df_2 <- as.data.frame(rowQuantiles(c2,probs = c(0.025, 0.5, 0.975)))
df_2$x <- days_vec

df_3 <- as.data.frame(rowQuantiles(c3, probs = c(0.025, 0.5, 0.975)))
df_3$x <- days_vec

df_4 <- as.data.frame(rowQuantiles(c4,probs = c(0.025, 0.5, 0.975)))
df_4$x <- days_vec

df_5 <- as.data.frame(rowQuantiles(c5,probs = c(0.025, 0.5, 0.975)))
df_5$x <- days_vec


tmp1 <- data.frame(date_index = mo, count = round(df_1$`50%`), incident = "0_1", id = 1)
tmp2 <- data.frame(date_index = mo, count = round(df_2$`50%`), incident = "1_4", id = 2)
tmp3 <- data.frame(date_index = mo, count = round(df_3$`50%`), incident = "5_14", id = 3)
tmp4 <- data.frame(date_index = mo, count = round(df_4$`50%`), incident = "15_65", id = 4)
tmp5 <- data.frame(date_index = mo, count = round(df_5$`50%`), incident = "65p", id = 5)

linelist <- rbind(tmp1, tmp2, tmp3, tmp4,tmp5)
linelist$Outcome <- factor(linelist$incident, levels = c("0_1", "1_4", "5_14", "15_65", "65p"))
colors <- c("0_1" = "brown1", "1_4" = "gold", "5_14" = "lightseagreen",
            "15_65" = "grey43", "65p"="orange")

linelistg2<-linelist

stackall <- ggplot(linelist, aes(
  fill = Outcome,
  y = count, x = date_index
)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Incidence by age", x = " ", y = "Incidence\n (symptomatic) per wk") +
  theme_classic() +
  scale_x_date(date_breaks = "3 month", 
               date_labels = "%b-%Y", 
               limits = as.Date(c("2012-06-01", "2018-05-30"))) +
  theme( # legend.position = "bottom",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  ) +
  guides(col = guide_legend(title = NULL))

windows()
stackall
save(stack, file = here("output",  "stack_incidenceGII.rdata"))


########


incig1<-ggplot(linelistg1, aes(x=date_index, y=count)) +
  geom_line() +
  facet_wrap(~Outcome)+
  labs(title = "Incidence GI", x = " ", y = "Incidence\n (symptomatic) per year") +
  #theme_classic() +
  scale_x_date(date_breaks = "6 month", 
               date_labels = "%b-%Y", 
               limits = as.Date(c("2014-06-01", "2018-05-30"))) +
  ylim(0,10000)+
  theme( # legend.position = "bottom",
    #panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )#, scales="free")
  windows()
gridExtra::grid.arrange(incig1)


incig2<-ggplot(linelistg2, aes(x=date_index, y=count)) +
  geom_line() +
  facet_wrap(~Outcome)+
  labs(title = "Incidence GII", x = " ", y = "Incidence\n (symptomatic) per year") +
  #theme_classic() +
  scale_x_date(date_breaks = "6 month", 
               date_labels = "%b-%Y", 
               limits = as.Date(c("2014-06-01", "2018-05-30"))) +
  ylim(0,10000)+
  theme( # legend.position = "bottom",
    #panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )#, scales="free")
windows()
gridExtra::grid.arrange(incig2)


