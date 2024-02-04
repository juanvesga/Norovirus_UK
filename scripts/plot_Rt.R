library(matrixStats)
# Load runs
load(here(paste0(output_dir,"/runs_mcmc_",scenario,".RData"))) 
load(here(paste0(output_dir,"/chains_",scenario,".RData"))) 



# Workout R0 and Rt
nsamples<-dim(out$pars)[1]
nt<-dim(out$runs)[3]

R0_gi3<- array(0, c(nsamples, 1));  
R0_gi<- array(0, c(nsamples, 1));  
R0_gii4<- array(0, c(nsamples, 1));  
R0_gii<- array(0, c(nsamples, 1));  

Rt_gi3<- array(0, c(nsamples, nt));  
Rt_gi<- array(0, c(nsamples, nt));  
Rt_gii4<- array(0, c(nsamples, nt));  
Rt_gii<- array(0, c(nsamples, nt));  


for (ii in 1:nsamples){

  pars<-out$pars[ii,]/unlist(params$scaling_fac)
  rr_inf_asymp <- 0.05 
  D5  <- 2.5
  D5p <- 1.5
  Dasy<- 15 * rr_inf_asymp # duration asymp shedding
  w2  <- -3
  
  aduRR    <-pars["aduRR"]
  beta_gi3 <-pars["beta_1"]    
  beta_gi  <-pars["beta_2"] 
  beta_gii4<-pars["beta_3"]   
  beta_gii <-pars["beta_4"]
  w1       <-pars["w1_1"]
 
  # R0
  
  S<-rep(init[3,],times = length(init[3,]))
  C<-params$transmission 
  BC_gi3 <- beta_gi3*C
  BC_gi3[8:14,] <- BC_gi3[8:14,]*aduRR
  BC_gi <- beta_gi*C
  BC_gi[8:14,] <- BC_gi[8:14,]*aduRR
  BC_gii4 <- beta_gii4*C
  BC_gii4[8:14,] <- BC_gii4[8:14,]*aduRR
  BC_gii <- beta_gii*C
  BC_gii[8:14,] <- BC_gii[8:14,]*aduRR
  
  D<-C*0
  D[,1:5]<-D5 + Dasy
  D[,6:14]<-D5p + Dasy
  
  R0_gi3[ii] <- eigen(S*BC_gi3*D)$values[1] 
  R0_gi[ii] <- eigen(S*BC_gi*D)$values[1] 
  R0_gii4[ii] <- eigen(S*BC_gii4*D)$values[1] 
  R0_gii[ii] <- eigen(S*BC_gii*D)$values[1] 

  # Rt   
  tt<-seq(1,nt,1)
  seasonality<- (1 + w1*cos((2*pi*tt)/364 + (w2/12)*pi))
  
  id<-out$idx
  sus_gi3<-out$runs[id$sus_gi3,ii,]/out$runs[id$sus_gi3,ii,2]
  sus_gi <-out$runs[id$sus_gi,ii,]/out$runs[id$sus_gi,ii,2]
  sus_gii4<-out$runs[id$sus_gii4,ii,]/out$runs[id$sus_gii4,ii,2] 
  sus_gii<-out$runs[id$sus_gii,ii,]/out$runs[id$sus_gii,ii,2]
  
  Rt_gi3[ii,]<-R0_gi3[ii] * seasonality*sus_gi3
  Rt_gi[ii,]<-R0_gi[ii] * seasonality*sus_gi
  Rt_gii4[ii,]<-R0_gii4[ii] * seasonality*sus_gii4
  Rt_gii[ii,]<-R0_gii[ii] * seasonality*sus_gii
  
}



# Plot R0

df_qtls <- as.data.frame(rowQuantiles(t(cbind(R0_gi3, R0_gi, R0_gii4, R0_gii)),
                                      probs = c(0.025, 0.5, 0.975)))

names(df_qtls)<-paste(c("low","mid","up"))

df_qtls$x <- factor(c("GI.3","Other GI","GII.4","Other GII"))


viol_col <- "#ff4d4d"# "chartreuse3" # "yellow3"
err_col <- "black"
data_col <- "black"

figR0 <- ggplot() +
  geom_point(data = df_qtls, mapping = aes(x = x, y = mid), size = 2, shape = 15) +
  geom_errorbar(data = df_qtls, aes(x = x, ymin = low, ymax = up),
                width = .2, position = position_dodge(.9)) +
  labs(title = "Modeled R0",
       x = "Strain", y = "Basic Reproduction Number ") +
  theme_minimal() +
  ylim(0, 7) +
  geom_hline(yintercept=1,linetype="dashed", color = "red")+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )

## Rt

qtlsgi3 <- as.data.frame(rowQuantiles(t(Rt_gi3),
                                      probs = c(0.025, 0.5, 0.975)))
qtlsgi <- as.data.frame(rowQuantiles(t(Rt_gi),
                        probs = c(0.025, 0.5, 0.975)))
qtlsgii4 <- as.data.frame(rowQuantiles(t(Rt_gii4),
                          probs = c(0.025, 0.5, 0.975)))
qtlsgii <- as.data.frame(rowQuantiles(t(Rt_gii),
                         probs = c(0.025, 0.5, 0.975)))

qtlsgi3<-qtlsgi3[-c(1),]
qtlsgi<-qtlsgi[-c(1),]
qtlsgii4<-qtlsgii4[-c(1),]
qtlsgii<-qtlsgii[-c(1),]


Rt_gi3<-as.data.frame((qtlsgi3))
Rt_gi<-as.data.frame((qtlsgi))
Rt_gii4<-as.data.frame((qtlsgii4))
Rt_gii<-as.data.frame((qtlsgii))

Rt_gi$x<-days_vec
Rt_gi3$x<-days_vec
Rt_gii4$x<-days_vec
Rt_gii$x<-days_vec

# Rt_gi3<-reshape2::melt(Rt_gi3, id.vars = "x")
# Rt_gi<-reshape2::melt(Rt_gi, id.vars = "x")
# Rt_gii4<-reshape2::melt(Rt_gii4, id.vars = "x")
# Rt_gii<-reshape2::melt(Rt_gii, id.vars = "x")

Rt_gi3$group<-"GI.3"
Rt_gi$group<-"Other GI"
Rt_gii4$group<-"GII.4"
Rt_gii$group<-"Other GII"

Rt<- rbind(Rt_gi3, Rt_gi, Rt_gii4, Rt_gii)
head(Rt)


figRt<- ggplot(data=Rt, aes(x=x, group=group)) +
  geom_line(aes(y = `50%`, color=group), lwd = 1) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill=group), alpha = 0.5) +
  geom_hline(yintercept=1,linetype="dashed", color = "red")+
  labs(title = "Modeled Rt",
       x = "Time (days)", y = "Effective Reproduction Number ") +
  scale_fill_manual(values=c('skyblue3','firebrick2','yellow2','grey28'))+
  theme_minimal() +
  ylim(0, max(Rt$`97.5%`)*1.2) +
  geom_hline(yintercept=1,linetype="dashed", color = "red")+
  theme(
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11), legend.key = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1)
  )


##
gridExtra::grid.arrange(figR0)
gridExtra::grid.arrange(figRt)

