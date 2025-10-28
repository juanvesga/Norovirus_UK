rm(list = ls())
gc()  # garbage collection

# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root             <- here::here()
infile0          <- file.path(root,"output", "parameters_short.qs2")
infile_dataS     <- file.path(root,"output", "data_short2.qs2")
infile_datafit     <- file.path(root,"output", "data_short.qs2")
infile_dataL     <- file.path(root,"output", "data_long.qs2")
infile_dataPlots <- file.path(root, "output", "data_for_plots.qs2")
infile_input     <- file.path(root,"output", "params_list_short.qs2")

#Functions
source(file.path(root, "R", "modify_attach.R"))
source(file.path(root, "R", "utility_functions.R"))
source(file.path(root, "R", "update_interventions.R"))
source(file.path(root, "R", "collect_function.R"))
source(file.path(root, "R", "r0_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(abind, include.only = c("abind"))
modify_attach(tidyr, include.only = c("gather"))

library(odin2)
library(dust2)
library(ggplot2)
library(dplyr)
library(foreach)
library(matlib)
library(matrixStats)
library(tidyr)

# Select model ------------------------------------------------------------


# Cross protection 5%

#model       <- "_0" # Simple SEIAR
#model       <- "_1" # Simple SEIAR no reinf no cross-prot
model       <- "_2" # Full version with reinf and cross-protection
#model      <- "_3" # Full version with reinf and cross-protection and drop immunity


# Cross protection 25%

#model      <- "_4" # Simple SEIAR
#model      <- "_5" # Simple SEIAR no reinf no cross-prot
#model      <- "_6" # Full version with reinf and cross-protection
#model      <- "_7" # Full version with reinf and cross-protection and drop immunity


# Cross protection 50%

#model       <- "_8" # Simple SEIAR
#model      <- "_9" # Simple SEIAR no reinf no cross-prot
#model      <- "_10" # Full version with reinf and cross-protection
#model      <- "_11" # Full version with reinf and cross-protection and drop immunity




# Source Odin model

infile_runs      <- file.path(root,"output", paste0("mcmc_fits",model,".qs2"))


# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------

extract_iteration <- function(pars, i) {
  lapply(pars, function(param) param[[i]])
}



data_all   <- qs_read(infile_dataS)
data_fit   <- qs_read(infile_datafit)
data_plots <- qs_read(infile_dataPlots)
mcmc_runs  <- qs_read(infile_runs)
sims       <- mcmc_runs$runs
pars       <- mcmc_runs$parameters
nsim     <- dim(sims)[2]
p          <- qs_read(infile0)

endsim<-24109
times<-seq(0,endsim)
step_yr<-365


nt<-24109#length(pars_list$fixed_pars$school_time)

idx<- seq(1,dim(sims)[1],1) 
names(idx)<-paste(rownames(sims))



agg_ages<-function(index,sim){
  
  
  id<-c(idx[paste("inc_gi3",as.character(index[1]),sep="_")],
        idx[paste("inc_gi",as.character(index[1]),sep="_")],
        idx[paste("inc_gii4",as.character(index[1]),sep="_")],
        idx[paste("inc_gii",as.character(index[1]),sep="_")])
  
  id2<- idx[paste("pop_all",as.character(index[1]),sep="")]
  times<- (step_yr*2):dim(sim)[3]
  
  inc<-sim[id,,times]
  pop<-sim[id2,,times]
  
  
  if (length(index)>1){
    for (ii in 2:length(index)){
      
      id<-c(idx[paste("inc_gi3",as.character(index[ii]),sep="_")],
            idx[paste("inc_gi",as.character(index[ii]),sep="_")],
            idx[paste("inc_gii4",as.character(index[ii]),sep="_")],
            idx[paste("inc_gii",as.character(index[ii]),sep="_")])
      
      id2<- idx[paste("pop_all",as.character(index[1]),sep="")]
      
      inc<- inc + sim[id,,times]
      pop<- pop + sim[id2,,times] 
      
    }
  }
  
  tmp1<-colSums(inc,dims = 1 )
  
  out<-list(
    inc=t(apply(tmp1, 1, cumsum)),
    pop=pop
  )
  
  return(out)
  
}

agg_ages_strain<-function(index,sim,strain){
  
  id<-idx[paste("inc",strain,as.character(index[1]),sep="_")]
  
  id2<- idx[paste("pop_all",as.character(index[1]),sep="")]
  times<- (step_yr*2):dim(sim)[3]
  
  inc<-sim[id,,times]
  pop<-sim[id2,,times]
  
  
  if (length(index)>1){
    for (ii in 2:length(index)){
      id<-idx[paste("inc",strain,as.character(index[ii]),sep="_")]
      
      id2<- idx[paste("pop_all",as.character(index[ii]),sep="")]
      
      inc<- inc + sim[id,,times]
      pop<- pop + sim[id2,,times] 
      
    }
  }
  
  tmp1<-inc
  
  out<-list(
    inc=t(apply(tmp1, 1, cumsum)),
    pop=pop
  )
  
  return(out)
  
}



index<-1:4
inc_a0_4<-agg_ages(index,sims)$inc
inc_a0_4_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a0_4_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a0_4_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a0_4_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a0_4<-agg_ages(index,sims)$pop

le<-2 # time point 
infpc_0_4<-inc_a0_4[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gi3<-inc_a0_4_gi3[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gi<-inc_a0_4_gi[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gii4<-inc_a0_4_gii4[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])
infpc_0_4_gii<-inc_a0_4_gii[,step_yr*le]/rowMeans(pop_a0_4[,1:step_yr*le])

##
index<-5:8
inc_a5_15<-agg_ages(index,sims)$inc
inc_a5_15_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a5_15_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a5_15_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a5_15_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a5_15<-agg_ages(index,sims)$pop
le<-5
infpc_5_15<-inc_a5_15[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gi3<-inc_a5_15_gi3[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gi<-inc_a5_15_gi[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gii4<-inc_a5_15_gii4[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])
infpc_5_15_gii<-inc_a5_15_gii[,step_yr*le]/rowMeans(pop_a5_15[,1:step_yr*le])


index<-9
inc_a16_65<-agg_ages(index,sims)$inc
inc_a16_65_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a16_65_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a16_65_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a16_65_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a16_65<-agg_ages(index,sims)$pop

le<-24
infpc_16_65<-inc_a16_65[,step_yr*le]/rowMeans(pop_a16_65[,1:step_yr*le])
infpc_16_65_gi3<-inc_a16_65_gi3[,step_yr*le]/rowMeans(pop_a16_65[,1:step_yr*le])
infpc_16_65_gi<-inc_a16_65_gi[,step_yr*le]/rowMeans(pop_a16_65[,1:step_yr*le])
infpc_16_65_gii4<-inc_a16_65_gii4[,step_yr*le]/rowMeans(pop_a16_65[,1:step_yr*le])
infpc_16_65_gii<-inc_a16_65_gii[,step_yr*le]/rowMeans(pop_a16_65[,1:step_yr*le])


index<-10
inc_a65p<-agg_ages(index,sims)$inc
inc_a65p_gi3<-agg_ages_strain(index,sims,'gi3')$inc
inc_a65p_gi<-agg_ages_strain(index,sims,'gi')$inc
inc_a65p_gii4<-agg_ages_strain(index,sims,'gii4')$inc
inc_a65p_gii<-agg_ages_strain(index,sims,'gii')$inc

pop_a65p<-agg_ages(index,sims)$pop

le<-7
infpc_65p<-inc_a65p[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gi3<-inc_a65p_gi3[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gi<-inc_a65p_gi[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gii4<-inc_a65p_gii4[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])
infpc_65p_gii<-inc_a65p_gii[,step_yr*le]/rowMeans(pop_a65p[,1:step_yr*le])

len=length(infpc_5_15)
age_cats<-4
df1<-data.frame(
  y = c(infpc_0_4_gi3,
        infpc_0_4_gi3+infpc_5_15_gi3,
        infpc_0_4_gi3+infpc_5_15_gi3+infpc_16_65_gi3,
        infpc_0_4_gi3+infpc_5_15_gi3+infpc_16_65_gi3+infpc_65p_gi3),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("GI.3",len*age_cats)
  
)

df2<-data.frame(
  y = c(infpc_0_4_gi,
        infpc_0_4_gi+infpc_5_15_gi,
        infpc_0_4_gi+infpc_5_15_gi+infpc_16_65_gi,
        infpc_0_4_gi+infpc_5_15_gi+infpc_16_65_gi+infpc_65p_gi),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("Other GI",len*age_cats)
  
)

df3<-data.frame(
  y = c(infpc_0_4_gii4,
        infpc_0_4_gii4+infpc_5_15_gii4,
        infpc_0_4_gii4+infpc_5_15_gii4+infpc_16_65_gii4,
        infpc_0_4_gii4+infpc_5_15_gii4+infpc_16_65_gii4+infpc_65p_gii4),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("GII.4",len*age_cats)
  
)

df4<-data.frame(
  y = c(infpc_0_4_gii,
        infpc_0_4_gii+infpc_5_15_gii,
        infpc_0_4_gii+infpc_5_15_gii+infpc_16_65_gii,
        infpc_0_4_gii+infpc_5_15_gii+infpc_16_65_gii+infpc_65p_gii),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("Other GII",len*age_cats)
  
)

df5<-data.frame(
  y = c(infpc_0_4,
        infpc_0_4+infpc_5_15,
        infpc_0_4+infpc_5_15+infpc_16_65,
        infpc_0_4+infpc_5_15+infpc_16_65+infpc_65p),
  Age=c(rep("at 4",len),
        rep("at 15",len),
        rep("at 65",len),
        rep("at 80",len)),
  Strain=rep("All",len*age_cats)
  
)


df<-rbind(df1,df2,df3,df4,df5)
df$Age<-factor(df$Age,
               levels=c("at 4","at 15","at 65","at 80"))
df$Strain<-factor(df$Strain,
                  levels=c("GI.3","Other GI","GII.4","Other GII","All"))

library(dplyr)
library(viridis)
plot_strains<-df %>%
  ggplot( aes(x=Age, y=y, fill=Strain)) +
  #geom_jitter(aes(color=Strain), size=0.4, alpha=0.5) +
  geom_boxplot(outlier.shape = NA)  +
  geom_boxplot(aes(color = Strain),
               fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0,
               show.legend = F)+
  theme_light() +
  scale_fill_manual(values = c('skyblue3','firebrick2','yellow3','grey40','slateblue2'))+

  #ggtitle("Accrued Norovirus infection") +
  ylim(0,8)+
  xlab("Age") + ylab("AGE episodes")  +
  theme(
    #legend.position = "none",
    panel.background = element_blank(),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10), legend.key = element_blank(),
    # plot.tag = element_text(size = 12, face = "bold"),
    #  plot.tag.position = c(0.1, 0.95),
    panel.border = element_rect(colour = "black", fill=NA)
  )


windows()
plot_strains


len=length(infpc_5_15)
df1<- data.frame(y=c(infpc_0_4,
                     infpc_0_4+infpc_5_15, 
                     infpc_0_4+infpc_5_15+infpc_16_65,
                     infpc_0_4+infpc_5_15+infpc_16_65+infpc_65p),
                 age=c(rep("at 4",len),
                       rep("at 15",len),
                       rep("at 65",len),
                       rep("at 80",len)))

df1$age<-factor(df1$age,
                levels=c("at 4","at 15","at 65","at 80")) 


plot_all<-df1 %>%
  ggplot( aes(x=age, y=y, fill=age)) +
  geom_jitter(color="grey", size=0.4, alpha=0.5) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  theme_minimal() +
  ggtitle("Accrued Norovirus infection") +
  xlab("Age") + ylab("Mean number of symptomatic events")                 

windows()
gridExtra::grid.arrange(plot_strains)
windows()
gridExtra::grid.arrange(plot_all)
######

library(dplyr)

# Calculate mean, standard error, and 95% CI for each age group
summary_stats <- df1 %>%
  group_by(age) %>%
  summarise(
    n = n(),
    mean_y = mean(y),
    sd_y = sd(y),
    se_y = sd_y / sqrt(n),
    ci_lower = mean_y - qt(0.975, df = n - 1) * se_y,
    ci_upper = mean_y + qt(0.975, df = n - 1) * se_y
  )

print(summary_stats)



rm(sims,out)




