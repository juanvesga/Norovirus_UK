
parameter_plot<-function(mcmc,nsim,scale,scenario){
  
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lemon)
  
  
  posteriors<-mcmc$pars
  
  pars <- posteriors[sample(nrow(posteriors), nsim), ]
  pars<-as.data.frame(t(t(pars) / unlist(scale)))
  
  pars_m<-reshape2::melt(pars)
  colnames(pars_m)<-c("parameter","value")
  pars_m$parameter<-as.character(pars_m$parameter)
  
  parplot<-ggplot(pars_m, aes(x=value, fill=parameter)) +
    geom_histogram() +
    facet_rep_wrap(~parameter, scales="free")+
    theme(
      axis.text = element_text(colour = "black", size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11), legend.key = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1)
    )
  
  
  
  #####################
  c1<-which(mcmc$chain==1)
  c2<-which(mcmc$chain==2)
  c3<-which(mcmc$chain==3)
  c4<-which(mcmc$chain==4)
  
  chain1<-(mcmc$probabilities[c1,3])
  chain2<-(mcmc$probabilities[c2,3])
  chain3<-(mcmc$probabilities[c3,3])
  chain4<-(mcmc$probabilities[c4,3])
  
  lims<-c(max(mcmc$probabilities[,3]), min(mcmc$probabilities[,3]) )
  
  llkplot=plot(chain1,type = "l", col="limegreen",ylab = "Log-Posterior",
               ylim = lims, xlab = "Iteration")
  lines(chain2,type = "l", col="yellow3")
  lines(chain3,type = "l", col="royalblue1")
  lines(chain4,type = "l", col="violetred")
  legend(x = "topright",  box.lwd = 2 , title="",  
         legend=c("Chain 1", "Chain 2","Chain 3","Chain 4"),  
         fill = c("limegreen","yellow3","royalblue1","violetred"))
  
  # Immunity ----------------------------------------------------------------
  # Do a barchart with errorbars with the immunity parameters using ggplot2
  
  if (scenario=="imm1par"){
    
    df<-data.frame(infection1=pars$imm_yr,
                   infection2=pars$imm_yr,
                   infection3=pars$imm_yr,
                   infection4=pars$imm_yr)
    title_text<-"Immunity: One tier waning"
    
  }else if (scenario=="imm2par" || scenario=="imm2par_noreinfection"){
    
    df<-data.frame(infection1=pars$imm_yr,
                   infection2=pars$imm_yr*pars$imm_fac,
                   infection3=pars$imm_yr*pars$imm_fac,
                   infection4=pars$imm_yr*pars$imm_fac)
    
    title_text<-"Immunity: Two tier waning No reinfection"
    
    
  }else if (scenario=="imm4par"){
    
    df<-data.frame(infection1=pars$imm_yr,
                   infection2=pars$imm_yr*pars$imm_fac1,
                   infection3=pars$imm_yr*pars$imm_fac2,
                   infection4=pars$imm_yr*pars$imm_fac3)
    
    title_text<-"Immunity: Four tier waning"
    
  } else if (scenario=="immdrop"){
    
    df<-data.frame(infection1=pars$imm_yr,
                   infection2=pars$imm_yr,
                   infection3=pars$imm_yr,
                   infection4=pars$imm_yr)
    title_text<-"Immunity: Drop one infection"
  } else if (scenario=="immdrop_noreinfection"){
    
    df<-data.frame(infection1=pars$imm_yr,
                   infection2=pars$imm_yr,
                   infection3=pars$imm_yr,
                   infection4=pars$imm_yr)
    title_text<-"Immunity: Drop one infection, No reinfection"}
  

  # use ggplot2 to plot a boxplot for each infection event
  
  df_m<-reshape2::melt(df)
  
  box<-ggplot(df_m, aes(x=variable, y=value)) +
    geom_boxplot(fill="indianred", outlier.colour = NA) +
    ylim(0, 100)+
    theme_minimal()+
    theme(
      axis.text = element_text(colour = "black", size = 10, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11), legend.key = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1)
    )+
    labs(title = title_text, x = "Infection event", y = "Average waning (Years)")
  
  
  
  df_qtl<-as.data.frame(t(data.frame(
    infection1=quantile(df$infection1,c(0.025,0.5,0.975)),
    infection2=quantile(df$infection2,c(0.025,0.5,0.975)),
    infection3=quantile(df$infection3,c(0.025,0.5,0.975)),
    infection4=quantile(df$infection4,c(0.025,0.5,0.975)))))
  
  df_qtl$ninfection<- row.names(df_qtl)
  
  df_m<-reshape2::melt(df_qtl)
  
  colnames(df_m)<-c("ninfection","qtl","value")
  
  df_m$ninfection<-as.character(df_m$ninfection)
  
  immplot<-ggplot(df_m[df_m$qtl=='50%',], aes(x=ninfection, y=value)) +
    geom_bar(stat="identity", position=position_dodge(), fill="indianred")+
    theme_minimal()+
    theme(
      axis.text = element_text(colour = "black", size = 10, face = "bold"),
      axis.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11), legend.key = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1)
    )+
    labs(title = title_text, x = "Infection event", y = "Average waning (Years)")+
    # add the value as text at the top of each bar
    geom_text(aes(label=round(value,2)), vjust=-0.5, size=3.5)
  
  
  return(
    list(
      plot1=parplot,
      llkplot,
      immplot=immplot,
      immbox=box
    ))
  
  
}

