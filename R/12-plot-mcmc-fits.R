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
source(file.path(root, "R", "plot_functions.R"))

# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))
modify_attach(tidyr, include.only = c("gather"))


library(posterior)
library(bayesplot)
library(matrixStats)
library(ggplot2)
#library(introdataviz)
library(see)
library(dplyr)
library(lubridate)
library(magrittr)
library(tidyr)
# Select model ------------------------------------------------------------


# Cross protection 5%

#model       <- "_0" # Simple SEIAR
#model       <- "_1" # Simple SEIAR no reinf no cross-prot
model       <- "_2" # Full version with reinf and cross-protection
#model      <- "_3" # Full version with reinf and cross-protection and drop immunity


# Cross protection 25%

#model      <- "_4" # Simple SEIAR
#model      <- "_5" # Simple SEIAR no reinf no cross-prot#
#model      <- "_6" # Full version with reinf and cross-protection
#model      <- "_7" # Full version with reinf and cross-protection and drop immunity


# Cross protection 50%

#model       <- "_8" # Simple SEIAR
#model      <- "_9" # Simple SEIAR no reinf no cross-prot
#model      <- "_10" # Full version with reinf and cross-protection
#model      <- "_11" # Full version with reinf and cross-protection and drop immunity




# Source Odin model

infile_runs      <- file.path(root,"output", paste0("mcmc_fits",model,".qs2"))
fitsfile         <- file.path(root,"output", paste0("model_fits",model,".png"))


# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------

data_all   <- qs_read(infile_dataS)
data_fit   <- qs_read(infile_datafit)
data_plots <- qs_read(infile_dataPlots)
mcmc_runs  <- qs_read(infile_runs)
state      <- mcmc_runs$runs
nsamp <- dim(state)[2]

endsim<-24109
times<-seq(0,endsim)

av_cases<-if (length(which(!is.na(data_all$reported_cases_week)))>53) FALSE else TRUE

tt<-which(times%in% data_all$time)
#state<-state[,,tt]

pars<-qs_read(infile0)


# SGSS reported time series -----------------------------------------------

nice_cols = c( "#007A87","#FF5A5F","#FFB400", "#007A87", 
               "#8CE071",  "#00D1C1", "#FFAA91", "#B4A76C", 
               "#9CA299", "#565A5C", "#00A04B", "#E54C20")



reported<-(colSums(state[c("reported_wk_1","reported_wk_2","reported_wk_3","reported_wk_4",
                           "reported_wk_5","reported_wk_6","reported_wk_7","reported_wk_8",
                           "reported_wk_9","reported_wk_10"), , ]))



ii<-which(!is.na(data_all$reported_cases_week))
model_t<-data_all$time[which(!is.na(data_all$reported_cases_week))]+1 # model strats from 0 and to match by week output



df_s <- as.data.frame(
  rowQuantiles(t(reported[,model_t]),
               probs = c(0.025, 0.5, 0.975)
  )
)

df<-data.frame(t(reported[,model_t]))

if(av_cases){

  df_s$x<-data_plots$total_cases_avg$date
  
df$x<-data_plots$total_cases_avg$date# total_cases$date 

dat_sim <- reshape2::melt(df, id = "x")

df_d <- data.frame(
  x = data_plots$total_cases_avg$date,
  y = data_plots$total_cases_avg$avg_cases
)

start_date <- data_plots$total_cases_avg$date[1]# ymd("2014-01-01")
end_date <- data_plots$total_cases_avg$date[length(data_plots$total_cases_avg$date)]# ymd("2020-31-12")

steps<-"1 months"

} else{

  df_s$x<-data_plots$total_cases$date
  
df$x<-data_plots$total_cases$date 

dat_sim <- reshape2::melt(df, id = "x")


df_d <- data.frame(
  x = data_plots$total_cases$date,
  y = data_plots$total_cases$cases
)

start_date <- data_plots$total_cases$date[1]# ymd("2014-01-01")
end_date <- data_plots$total_cases$date[length(data_plots$total_cases$date)]# ymd("2020-31-12")

steps<-"4 months"


}



col_0<-"#cc0044"
col_1<-"#E54C20"
col_2<-"#8c8cd9"
col_3<-"#00A04B"
col_4<-"grey40"


fits_sgss <- ggplot() +
  geom_col(data = df_d, aes(x = x, y = y), width = 7, 
           fill="grey18",alpha=0.5) +
  geom_line(data = dat_sim, aes(x = x, y=value, group=variable),
            col = nice_cols[4] , alpha = 0.05, lwd = 0.001) +
  geom_line(data = df_s, aes(x = x, y=`50%`), col = "black",lwd = 0.5) +
  labs( x = "", y = "Weekly cases\n reported") +
  theme_minimal() +
  ylim(c(0,max(df_d$y)*1.4))+
  scale_x_date(
    date_labels = "%b-%Y",     # Sept-2014
    date_breaks = steps,  # Every 4 months
    date_minor_breaks = "1 week",
    limits = c( start_date, end_date), expand=c(0,0)
  )+
  theme(
    legend.position = "none",
    panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.title = element_text(size = 10)
  )

#}
print(fits_sgss)

# SGSS by Strain ----------------------------------------------------------



df_data<-data.frame(
  date          = data_plots$reported_monthly_strain$date,
  reported_gi3  = data_plots$reported_monthly_strain$reported_gi3,
  reported_gi   = data_plots$reported_monthly_strain$reported_gi,
  reported_gii4 = data_plots$reported_monthly_strain$reported_gii4,
  reported_gii  = data_plots$reported_monthly_strain$reported_gii)



df_long <-df_data  %>%
  tidyr::pivot_longer(
    cols    = 2:5,
    names_to = "variable",
    values_to = "value"
  ) %>%
  # Extract strain and type (observations vs modelled)
  mutate(
    strain = case_when(
      grepl("gi3", variable) ~ "GI3",
      grepl("gi[^i]|gi$", variable) ~ "GI",
      grepl("gii4", variable) ~ "GII4",
      grepl("gii[^4]|gii$", variable) ~ "GII"
    ),
    type = case_when(
      grepl("gi", variable) ~ "Observations",
    )
  )



# Model
ii<-which(!is.na(data_all$reported_gi3))
model_t<-data_all$time[which(!is.na(data_all$reported_gi3))]+1





df_gi3<- data.frame(
  date=data_plots$reported_monthly_strain$date,
  t(state["reported_gi3",,model_t]))
df_gi3$strain <- "GI3"

df_gi<- data.frame(
  date=data_plots$reported_monthly_strain$date,
  t(state["reported_gi",,model_t]))
df_gi$strain <- "GI"

df_gii4<- data.frame(
  date=data_plots$reported_monthly_strain$date,
  t(state["reported_gii4",,model_t]))
df_gii4$strain <- "GII4"

df_gii<- data.frame(
  date=data_plots$reported_monthly_strain$date,
  t(state["reported_gii",,model_t]))
df_gii$strain <- "GII"

df_model<-rbind(
  df_gi3,df_gi,df_gii4,df_gii)

# 
sim_df_monthly <- df_model %>%
  # Convert to long format (melt the simulation columns)
  tidyr::pivot_longer(
    cols = 2:(nsamp+1),  # Columns 2 to 401 (the simulation runs)
    names_to = "simulation_id",
    values_to = "cases"
  )
# # 
# 
# fits_bystrain<-ggplot() +
#   # Dark grey columns for observed data
#   geom_col(data = filter(df_long, type == "Observations"),
#            aes(x = date, y = value),
#            fill = "darkgrey",
#            width = 25) +
# 
#   # Individual simulation points as colored dots
#   geom_point(data = sim_df_monthly,
#              aes(x = date, y = cases, color = as.factor(strain)),
#              alpha = 0.5,  # semi-transparent
#              size = 1.5,
#              position = position_jitter(width = 2, height = 0)) +  # slight jitter for visibility
# 
#   # Transparent boxplot for simulations
#   geom_boxplot(data = sim_df_monthly,
#                aes(x = date, y = cases, group = date),
#                fill = NA,
#                color = "black",
#                alpha = 0.7,
#                outlier.shape = NA,
#                width = 15,
#                linewidth = 0.7)+  # make lines thicker (default is usually 0.5)
#   labs(
#          x = "",
#          y = "Monthly reported\n Cases"
#        ) +
# 
# # #   # Color scale for the dots
#   scale_color_manual( values = cols<-c("GI3"="dodgerblue",
#                                                     "GI" ="#cc0044",
#                                                     "GII4"="#FFB400",
#                                                     "GII"="#00A04B")) +
# 
#   facet_wrap(~strain)+#, scales = "free_y") +
#   theme_minimal()+
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
#     panel.border = element_rect(colour = "black", fill=NA),
#     axis.title = element_text(size = 10)
#   )
# 
# print(fits_bystrain)

fits_bystrain <- ggplot() +
  # Dark grey columns for observed data
  geom_col(data = filter(df_long, type == "Observations"),
           aes(x = date, y = value),
           fill = "darkgrey",
           width = 25) +
  # Individual simulation points as colored dots
  geom_point(data = sim_df_monthly,
             aes(x = date, y = cases, color = as.factor(strain)),
             alpha = 0.5,  # semi-transparent
             size = 1.5,
             position = position_jitter(width = 2, height = 0)) +  # slight jitter for visibility
  # Transparent boxplot for simulations
  geom_boxplot(data = sim_df_monthly,
               aes(x = date, y = cases, group = date),
               fill = NA,
               color = "black",
               alpha = 0.7,
               outlier.shape = NA,
               width = 15,
               linewidth = 0.7) +  # make lines thicker
  # Add strain labels inside panels at top-left
  geom_text(data = sim_df_monthly %>% 
              group_by(strain) %>% 
              summarise(min_date = min(date), 
                        max_cases = max(cases, na.rm = TRUE)),
            aes(x = min_date, y = Inf, label = strain),
            hjust = 0.1, vjust = 1.5, 
            size = 4, fontface = "bold") +
  labs(x = "",
       y = "Monthly reported\n Cases") +
  # Color scale for the dots
  scale_color_manual(values = cols <- c("GI3" = "dodgerblue",
                                        "GI" = "#cc0044",
                                        "GII4" = "#FFB400",
                                        "GII" = "#00A04B")) +
  facet_wrap(~strain) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.title = element_text(size = 10),
    strip.text = element_blank()  # Hide the default facet labels
  )

print(fits_bystrain)


# IID incidence by strain and age -----------------------------------------
# Prepare data for incidence plot
# Prepare data for incidence plot
ii <- data_all$time[which(!is.na(data_all$cases_a1))] + 1

modelled_irate <- rbind(
  ((state['inc_year_gi3_1',,ii] +
      state['inc_year_gi_1',,ii] +
      state['inc_year_gii4_1',,ii] +
      state['inc_year_gii_1',,ii]) / (state['pop_by4age1',,ii])),
  
  ((state['inc_year_gi3_2',,ii] +
      state['inc_year_gi_2',,ii] +
      state['inc_year_gii4_2',,ii] +
      state['inc_year_gii_2',,ii]) / (state['pop_by4age2',,ii])),
  
  ((state['inc_year_gi3_3',,ii] +
      state['inc_year_gi_3',,ii] +
      state['inc_year_gii4_3',,ii] +
      state['inc_year_gii_3',,ii]) / (state['pop_by4age3',,ii])),
  
  ((state['inc_year_gi3_4',,ii] +
      state['inc_year_gi_4',,ii] +
      state['inc_year_gii4_4',,ii] +
      state['inc_year_gii_4',,ii]) / (state['pop_by4age4',,ii])),
  
  ((state['inc_year_gi3_5',,ii] +
      state['inc_year_gi_5',,ii] +
      state['inc_year_gii4_5',,ii] +
      state['inc_year_gii_5',,ii]) / (state['pop_by4age5',,ii]))
) * 1000

# Data frame for observed data (shifted left)
df_d <- data.frame(
  x = c(0.75, 1.75, 2.75, 3.75, 4.75),  # Numeric positions shifted left
  x_label = factor(c("[0 1)", "[1 4)", "[5 14)", "[15 64)", "65+"),
                   levels = c("[0 1)", "[1 4)", "[5 14)", "[15 64)", "65+")),
  inc = data_plots$data_iid2.c4$per1000personyears,
  low = data_plots$data_iid2.c4$CI_lower,
  up = data_plots$data_iid2.c4$CI_upper
)

# Model simulations data frame (shifted right)
df1 <- data.frame(t(modelled_irate)) 
colnames(df1) <- c("[0 1)", "[1 4)", "[5 14)", "[15 64)", "65+")
df_m <- reshape2::melt(df1)
df_m$variable <- factor(df_m$variable, 
                        levels = c("[0 1)", "[1 4)", "[5 14)", "[15 64)", "65+"))
# Add numeric x position for simulations (shifted right)
df_m$x_num <- as.numeric(df_m$variable) + 0.25

# Color settings
viol_col <- nice_cols[4]
data_col <- "black"

# Create plot
iid2_all <- ggplot() +
  # Data error bars with dots (on the left)
  geom_errorbar(
    data = df_d,
    aes(x = x, ymin = low, ymax = up), 
    width = 0.15
  ) +
  geom_point(
    data = df_d, 
    aes(x = x, y = inc, color = "Data (95% CI)"), 
    size = 2.5, 
    shape = 16
  ) +
  # Simulation colored dots with jitter (on the right)
  geom_point(
    data = df_m,
    aes(x = x_num, y = value),
    color = viol_col,
    position = position_jitter(width = 0.1, seed = 123), 
    size = 1.5,
    alpha = 0.3
  ) +
  # Transparent boxplot overlay (on the right)
  geom_boxplot(
    data = df_m,
    aes(x = x_num, y = value, group = variable, fill = "Posterior Density"),
    width = 0.35,
    outlier.shape = NA,
    alpha = 0.3,
    color = "black"
  ) +
  # Set x-axis breaks and labels
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5),
    labels = c("[0 1)", "[1 4)", "[5 14)", "[15 64)", "65+")
  ) +
  labs( 
    x = "Age", 
    y = "Incidence per 1k\nperson-year"
  ) +
  theme_minimal() +
  ylim(0, max(df_d$up) * 1.2) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10)
  )






# GII4 prevalence in children  --------------------------------------------

# Prepare data
ii <- data_all$time[which(!is.na(data_all$sero1))] + 1

observed_size <- c(
  103*10,
  107*10,
  121*10,
  124*10,
  122*10,
  109*10
)

# Data frame for observed data (shifted left)
df_d <- data.frame(
  x = c(0.75, 1.75, 2.75, 3.75, 4.75, 5.75),  # Numeric positions shifted left
  x_label = factor(c("[0 1)", "[1 2)", "[2 3)", "[3 4)", "[4 5)", "[5 6)"),
                   levels = c("[0 1)", "[1 2)", "[2 3)", "[3 4)", "[4 5)", "[5 6)")),
  sero = data_plots$dfsero$mean * 100,
  low = data_plots$dfsero$low * 100,
  up = data_plots$dfsero$up * 100
)

# Model simulations
sero_model <- rbind(
  state['seroprev1.2',,ii],
  state['seroprev2.3',,ii],
  state['seroprev3.4',,ii],
  state['seroprev4.5',,ii],
  state['seroprev5.6',,ii],
  state['seroprev6.7',,ii]
) * 100

# Model simulations data frame (shifted right)
df1 <- data.frame(t(sero_model)) 
colnames(df1) <- c("[0 1)", "[1 2)", "[2 3)", "[3 4)", "[4 5)", "[5 6)")
df_m <- reshape2::melt(df1)
df_m$variable <- factor(df_m$variable, 
                        levels = c("[0 1)", "[1 2)", "[2 3)", "[3 4)", "[4 5)", "[5 6)"))
# Add numeric x position for simulations (shifted right)
df_m$x_num <- as.numeric(df_m$variable) + 0.25

# Color settings
viol_col <- nice_cols[3]
data_col <- "black"

# Create plot
fits_sero <- ggplot() +
  # Data error bars with dots (on the left)
  geom_errorbar(
    data = df_d,
    aes(x = x, ymin = low, ymax = up), 
    width = 0.15
  ) +
  geom_point(
    data = df_d, 
    aes(x = x, y = sero, color = "Data (95% CI)"), 
    size = 2.5, 
    shape = 16  # Circle/dot instead of square
  ) +
  # Simulation colored dots with jitter (on the right)
  geom_point(
    data = df_m,
    aes(x = x_num, y = value),
    color = viol_col,
    position = position_jitter(width = 0.1, seed = 123), 
    size = 1.5,
    alpha = 0.3
  ) +
  # Transparent boxplot overlay (on the right)
  geom_boxplot(
    data = df_m,
    aes(x = x_num, y = value, group = variable, fill = "Posterior Density"),
    width = 0.35,
    outlier.shape = NA,
    alpha = 0.3,
    color = "black"
  ) +
  # Set x-axis breaks and labels
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 6),
    labels = c("[0 1)", "[1 2)", "[2 3)", "[3 4)", "[4 5)", "[5 6)")
  ) +
  labs(
    x = "Age", 
    y = "GII.4\nprevalence (%)"
  ) +
  theme_minimal() +
  ylim(0, max(df_d$up) * 1.2) +
  scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
  scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 10)
  )






# SGSS age reporting fractions --------------------------------------------


#ii<-which(!is.na(data_all$a0_event))
ii<-data_all$time[which(!is.na(data_all$a0_event))]+1 # model strats from 0 and to match by week output



a0_model<-((state['cumm_incday_gi3_1',,ii]+
              state['cumm_incday_gi_1',,ii]+
              state['cumm_incday_gii4_1',,ii]+
              state['cumm_incday_gii_1',,ii]) /
             (    
               state['cumm_incday_gii4_1',,ii]+
                 state['cumm_incday_gii4_2',,ii]+
                 state['cumm_incday_gii4_3',,ii]+
                 state['cumm_incday_gii4_4',,ii]+
                 state['cumm_incday_gii_1',,ii]+
                 state['cumm_incday_gii_2',,ii]+
                 state['cumm_incday_gii_3',,ii]+
                 state['cumm_incday_gii_4',,ii]+
                 state['cumm_incday_gi3_1',,ii]+
                 state['cumm_incday_gi3_2',,ii]+
                 state['cumm_incday_gi3_3',,ii]+
                 state['cumm_incday_gi3_4',,ii]+
                 state['cumm_incday_gi_1',,ii]+
                 state['cumm_incday_gi_2',,ii]+
                 state['cumm_incday_gi_3',,ii]+
                 state['cumm_incday_gi_4',,ii]))


a5_model<-((state['cumm_incday_gi3_2',,ii]+
              state['cumm_incday_gi_2',,ii]+
              state['cumm_incday_gii4_2',,ii]+
              state['cumm_incday_gii_2',,ii]) /
             (    
               state['cumm_incday_gii4_1',,ii]+
                 state['cumm_incday_gii4_2',,ii]+
                 state['cumm_incday_gii4_3',,ii]+
                 state['cumm_incday_gii4_4',,ii]+
                 state['cumm_incday_gii_1',,ii]+
                 state['cumm_incday_gii_2',,ii]+
                 state['cumm_incday_gii_3',,ii]+
                 state['cumm_incday_gii_4',,ii]+
                 state['cumm_incday_gi3_1',,ii]+
                 state['cumm_incday_gi3_2',,ii]+
                 state['cumm_incday_gi3_3',,ii]+
                 state['cumm_incday_gi3_4',,ii]+
                 state['cumm_incday_gi_1',,ii]+
                 state['cumm_incday_gi_2',,ii]+
                 state['cumm_incday_gi_3',,ii]+
                 state['cumm_incday_gi_4',,ii]))

a15_model<-((state['cumm_incday_gi3_3',,ii]+
               state['cumm_incday_gi_3',,ii]+
               state['cumm_incday_gii4_3',,ii]+
               state['cumm_incday_gii_3',,ii]) /
              (    
                state['cumm_incday_gii4_1',,ii]+
                  state['cumm_incday_gii4_2',,ii]+
                  state['cumm_incday_gii4_3',,ii]+
                  state['cumm_incday_gii4_4',,ii]+
                  state['cumm_incday_gii_1',,ii]+
                  state['cumm_incday_gii_2',,ii]+
                  state['cumm_incday_gii_3',,ii]+
                  state['cumm_incday_gii_4',,ii]+
                  state['cumm_incday_gi3_1',,ii]+
                  state['cumm_incday_gi3_2',,ii]+
                  state['cumm_incday_gi3_3',,ii]+
                  state['cumm_incday_gi3_4',,ii]+
                  state['cumm_incday_gi_1',,ii]+
                  state['cumm_incday_gi_2',,ii]+
                  state['cumm_incday_gi_3',,ii]+
                  state['cumm_incday_gi_4',,ii]))

a65_model<-((  state['cumm_incday_gi3_4',,ii]+
                 state['cumm_incday_gi_4',,ii]+
                 state['cumm_incday_gii4_4',,ii]+
                 state['cumm_incday_gii_4',,ii]) /
              (    
                state['cumm_incday_gii4_1',,ii]+
                  state['cumm_incday_gii4_2',,ii]+
                  state['cumm_incday_gii4_3',,ii]+
                  state['cumm_incday_gii4_4',,ii]+
                  state['cumm_incday_gii_1',,ii]+
                  state['cumm_incday_gii_2',,ii]+
                  state['cumm_incday_gii_3',,ii]+
                  state['cumm_incday_gii_4',,ii]+
                  state['cumm_incday_gi3_1',,ii]+
                  state['cumm_incday_gi3_2',,ii]+
                  state['cumm_incday_gi3_3',,ii]+
                  state['cumm_incday_gi3_4',,ii]+
                  state['cumm_incday_gi_1',,ii]+
                  state['cumm_incday_gi_2',,ii]+
                  state['cumm_incday_gi_3',,ii]+
                  state['cumm_incday_gi_4',,ii]))




age_model <- cbind(a0_model, a5_model, a15_model, a65_model) * 100

# Data frame for column bars (showing only mean of data)
df_d <- data.frame(
  x = factor(c("[0 4)", "[5 14)", "[15 64)", "65+"), 
             levels = c("[0 4)", "[5 14)", "[15 64)", "65+")),
  age_mean = data_plots$agg_age[, 2] * 100  # Only the mean
)

# Model simulations data frame
df1 <- data.frame(age_model) 
colnames(df1) <- c("[0 4)", "[5 14)", "[15 64)", "65+")
df_m <- reshape2::melt(df1)
df_m$variable <- factor(df_m$variable, 
                        levels = c("[0 4)", "[5 14)", "[15 64)", "65+"))

# Color settings
viol_col <- nice_cols[9]
data_col <- "darkgrey"  # Column color for data
point_col <- nice_cols[4]  # Color for simulation points

fits_age <- ggplot() +
  # Data as columns (mean only)
  geom_col(
    data = df_d,
    aes(x = x, y = age_mean),
    fill = data_col,
    alpha = 0.7,
    width = 0.6
  ) +
  # Model simulation points (colored dots with jitter)
  geom_point(
    data = df_m,
    aes(x = variable, y = value),
    color = point_col,
    position = position_jitter(width = 0.15, seed = 123), 
    size = 1.5,
    alpha = 0.3
  ) +
  # Transparent boxplot over the points
  geom_boxplot(
    data = df_m,
    aes(x = variable, y = value),
    width = 0.4,
    outlier.shape = NA,
    alpha = 0.3,
    fill = viol_col,
    color = "black"
  ) +
  labs(x = "Age", y = "Reported (%)") +
  theme_minimal() +
  ylim(0, max(df_d$age_mean) * 1.2) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text = element_text(colour = "black", size = 10),
    axis.title = element_text(size = 11)
  )





# Get plots together  -----------------------------------------------------
library(grid)
blank<-grid.rect(gp=gpar(col="white"))



p_fits<-gridExtra::grid.arrange(
  fits_sgss,
  fits_bystrain,
  iid2_all,
  fits_sero,
  fits_age,
  layout_matrix = rbind(c(1 ,1),
                        c(2, 3),
                        c(4, 5)))
  



windows()
gridExtra::grid.arrange(
  fits_sgss,
  fits_bystrain,
  iid2_all,
  fits_sero,
  fits_age,
  layout_matrix = rbind(c(1 ,1),
                        c(2, 3),
                        c(4, 5)))


ggsave(p_fits,filename=fitsfile)



