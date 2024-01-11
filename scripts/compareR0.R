library(matrixStats)

get_R0<-function(pars,scaling_fac,params,init,nt){
    
  rr_inf_asymp <- 0.05 
  D5  <- 2.5
  D5p <- 1.5
  Dasy<- 15 * rr_inf_asymp # duration asymp shedding
  w2  <- -3
  
  
  aduRR    <-pars["aduRR"]/scaling_fac["aduRR"]
  beta_gi3 <-pars["beta_1"]/scaling_fac["beta_1"]    
  beta_gi  <-pars["beta_2"]/scaling_fac["beta_2"] 
  beta_gii4<-pars["beta_3"]/scaling_fac["beta_3"]   
  beta_gii <-pars["beta_4"]/scaling_fac["beta_4"]
  w1       <-pars["w1_1"]/scaling_fac["w1_1"]


  tt<-seq(1,nt,1)
  seasonality<- (1 + w1*cos((2*pi*tt)/364 + (w2/12)*pi))
  # R0
  
  
  S<-rep(init[3,],times = length(init[3,]))
  C<-params$transmission 
  BC_gi3 <- beta_gi3*C*seasonality[1]
  BC_gi3[8:14,] <- BC_gi3[8:14,]*aduRR
  BC_gi <- beta_gi*C*seasonality[1]
  BC_gi[8:14,] <- BC_gi[8:14,]*aduRR
  BC_gii4 <- beta_gii4*C*seasonality[1]
  BC_gii4[8:14,] <- BC_gii4[8:14,]*aduRR
  BC_gii <- beta_gii*C*seasonality[1]
  BC_gii[8:14,] <- BC_gii[8:14,]*aduRR
  
  D<-C*0
  D[,1:5]<-D5 + Dasy
  D[,6:14]<-D5p + Dasy
  
  res<-list()
  res$R0_gi3 <- eigen(S*BC_gi3*D)$values[1] 
  res$R0_gi <- eigen(S*BC_gi*D)$values[1] 
  res$R0_gii4 <- eigen(S*BC_gii4*D)$values[1] 
  res$R0_gii <- eigen(S*BC_gii*D)$values[1] 

    return(res)
}

scenarios<-c("imm1par","imm2par","imm4par","imm2par_noreinfection","immdrop","immdrop_noreinfection")
samples<-nsamples



R0_gi3<- array(0, c(nsamples, length(scenarios)));  
R0_gi<- array(0, c(nsamples, length(scenarios)));
R0_gii4<- array(0, c(nsamples, length(scenarios)));
R0_gii<- array(0, c(nsamples, length(scenarios)));
wane<- array(0, c(nsamples, length(scenarios)));

for (ii in 1:length(scenarios)){


  scenario<-scenarios[ii]
  load(here("output",paste0("runs_mcmc_",scenario,".RData"))) 

  nt<-dim(out$runs)[3]
  
    for (iii in 1:nsamples){

    pars<-out$pars[iii,]

    res<-get_R0(pars,unlist(params$scaling_fac),params,init,nt)


    R0_gi3[iii,ii]<-res$R0_gi3  
    R0_gi[iii,ii]<- res$R0_gi
    R0_gii4[iii,ii]<-res$R0_gii4
    R0_gii[iii,ii]<- res$R0_gii
    wane[iii,ii]<-pars["imm_yr"]/unlist(params$scaling_fac)["imm_yr"]
    }

  
}



# Plot error bars of R0 for each strain and grouped by scenario

qtls_R0_gi3<-as.data.frame(rowQuantiles(t(R0_gi3), probs = c(0.025, 0.5, 0.975)))
qtls_R0_gi<-as.data.frame(rowQuantiles(t(R0_gi), probs = c(0.025, 0.5, 0.975)))
qtls_R0_gii4<-as.data.frame(rowQuantiles(t(R0_gii4), probs = c(0.025, 0.5, 0.975)))
qtls_R0_gii<-as.data.frame(rowQuantiles(t(R0_gii), probs = c(0.025, 0.5, 0.975)))

names(qtls_R0_gi3)<-paste(c("low","mid","up"))
names(qtls_R0_gi)<-paste(c("low","mid","up"))
names(qtls_R0_gii4)<-paste(c("low","mid","up"))
names(qtls_R0_gii)<-paste(c("low","mid","up"))

qtls_R0_gi3$strain<-rep("GI.3",times = nrow(qtls_R0_gi3))
qtls_R0_gi$strain<-rep("Other GI",times = nrow(qtls_R0_gi))
qtls_R0_gii4$strain<-rep("GII.4",times = nrow(qtls_R0_gii4))
qtls_R0_gii$strain<-rep("Other GII",times = nrow(qtls_R0_gii))

qtls_R0_gi3$scenario<-rep(scenarios,times = nrow(qtls_R0_gi3)/length(scenarios))
qtls_R0_gi$scenario<-rep(scenarios,times = nrow(qtls_R0_gi)/length(scenarios))
qtls_R0_gii4$scenario<-rep(scenarios,times = nrow(qtls_R0_gii4)/length(scenarios))
qtls_R0_gii$scenario<-rep(scenarios,times = nrow(qtls_R0_gii)/length(scenarios))

df_qtls<-rbind(qtls_R0_gi3,qtls_R0_gi,qtls_R0_gii4,qtls_R0_gii)
df_qtls$strain<-factor(df_qtls$strain,levels = c("GI.3","Other GI","GII.4","Other GII"))
# melt using reshape2 with scenario as id
df_qtlsm<-reshape2::melt(df_qtls,id.vars = c("strain","scenario"))
head(df_qtlsm)
head(df_qtls)


# use ggplot to plot errorbars and points for each strain grouped by scenario

 compR0<- ggplot(df_qtls,                                      # Grouped barplot using ggplot2
         aes(x = scenario,
             y = mid,
             fill = strain)) +
  geom_bar(stat = "identity",
           position = "dodge")+
  geom_errorbar(data = df_qtls, aes(x = scenario, ymin = low, ymax = up),
                width = .2, position = position_dodge(.9)) +
  labs(title = "Modeled R0 by genotype and immunity assumption",
       x = "Model Assumption", y = "R0 (Basic Reproduction Number) ") +
  geom_hline(yintercept=1,linetype="dashed", color = "red") +
  scale_fill_manual(values=c('skyblue3','firebrick2','yellow2','grey28'))+
  theme_minimal() +  
  theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.background = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 12, face = "bold")) 
 
 gridExtra::grid.arrange(compR0, ncol = 1, nrow = 1)
 
###
 
 qtl_wane<-as.data.frame(rowQuantiles(t(wane), probs = c(0.025, 0.5, 0.975)))
 
 names( qtl_wane)<-paste(c("low","mid","up"))

 qtl_wane$scenario<-rep(scenarios,times = nrow(qtls_R0_gi3)/length(scenarios))
 
 head(qtl_wane)
 
 qtl_wane <-qtl_wane[order(qtl_wane$mid),]
 qtl_wane$scenario<-factor(qtl_wane$scenario,levels = qtl_wane$scenario)
 
 # use ggplot to plot errorbars and points for each strain grouped by scenario
 
 compWane<- ggplot(qtl_wane,                                      # Grouped barplot using ggplot2
                 aes(x = scenario,
                     y = mid,
                     fill = scenario)) +
   geom_bar(stat = "identity",
            position = "dodge")+
   geom_errorbar(data = qtl_wane, aes(x = scenario, ymin = low, ymax = up),
                 width = .2, position = position_dodge(.9)) +
   labs(title = "Modeled mean waning time by immunity assumption",
        x = "Model Assumption", y = "Waning (year)") +
   scale_fill_manual(values=c('slateblue2','seagreen','hotpink3','orange2','dodgerblue2','brown2','darkslategrey'))+
   theme_minimal() +  
   theme(
     legend.position = "none",
     legend.title = element_blank(),
     legend.text = element_text(size = 12),
     panel.background = element_blank(),
     axis.text = element_text(colour = "black"),
     axis.text.y = element_text(size=12),
     axis.title.y = element_text(size = 12, face = "bold")) 
 
 gridExtra::grid.arrange(compWane, ncol = 1, nrow = 1)
 
 
 
 