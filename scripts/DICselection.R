rm(list = ls()) 
# DIC as 0.5 * var(-2*logLike) + mean(-2*logLike)


get_dic<-function(posterior){
  
  D<- -2*posterior
  DIC<- 0.5* var(D) + mean(D)
  return(DIC)
}

samples<-500
scenarios<-c("imm1par","imm2par","imm4par","imm2par_noreinfection","immdrop","immdrop_noreinfection")

posterior<-matrix(NA,nrow =samples, ncol=length(scenarios))
DIC<-seq(1,length(scenarios))*0
# Load posterior
for (i in 1:length(scenarios)){
  
  scenario<-scenarios[i]

load(here("output",paste0("chains_",scenario,".RData"))) 
postD<-processed_chains$probabilities[,3]
M <- postD[sample(length(postD), samples)]
posterior[,i]<-M
DIC[i]<-get_dic(M)
}


# Barchart of DIC values by scenario
df<-data.frame(scenarios,DIC)
df<-df[order(df$DIC),]
df$scenarios<-factor(df$scenarios,levels = df$scenarios)
dic_plot<-ggplot(df,aes(x=scenarios,y=DIC))+
  geom_bar(stat = "identity",fill="grey")+
  geom_text(aes(label=round(DIC,1)),vjust=-0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Scenario",y="DIC")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gridExtra::grid.arrange(dic_plot)


