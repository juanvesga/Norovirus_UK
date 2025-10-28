# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
root           <- here::here()



#Functions
source(file.path(root, "R", "modify_attach.R"))


# Packages
modify_attach(qs2, include.only = c("qs_save","qs_read"))


# -------------------------------------------------------------------------
# Load inputs ---------------------------------------------------------------
# -------------------------------------------------------------------------


## Cross protection 5%

#model       <- "_0" # Simple SEIAR
#model       <- "_1" # Simple SEIAR no reinf no cross-prot
#model       <- "_2" # Full version with reinf and cross-protection
#model        <- "_3" # Full

## Cross protection 25%

#model      <- "_4" # Simple SEIAR
#model      <- "_5" # Simple SEIAR no reinf no cross-prot
#model      <- "_6" # Full version with reinf and cross-protection
#model      <- "_7" # Full version with reinf and cross-protection and drop immunity


## Cross protection 50%

#model       <- "_8" # Simple SEIAR
#model      <- "_9" # Simple SEIAR no reinf no cross-prot
#model      <- "_10" # Full version with reinf and cross-protection
model      <- "_11" # Full version with reinf and cross-protection and drop immunity


if (model=="_0"||model=="_4"||model=="_8"){
  
  Model<-"Model_A_"
  
  if (model=="_0"){
    
    cross<-"5pct"
    
  }  else if (model=="_4"){
    
    cross<-"25pct"
    
    
  }else{
    cross<-"50pct"
  }
  
}else if((model=="_2"||model=="_6"||model=="_10")){
  
  Model<-"Model_B_"
  
  if (model=="_2"){
    
    cross<-"5pct"
    
  }  else if (model=="_6"){
    
    cross<-"25pct"
    
  }else{
    cross<-"50pct"
  }
  
}else if((model=="_3"||model=="_7"||model=="_11")){
  
  Model<-"Model_C_"
  
  if (model=="_3"){
    
    cross<-"5pct"
    
  }  else if (model=="_7"){
    
    cross<-"25pct"
    
  }else{
    cross<-"50pct"
  }
  
}




results_file  <- file.path(root, "output", paste0("monty_chains",model,".rds"))
sample_df_file     <- file.path(root, "output", paste0("sample_df",model,".qs2"))

out_results   <- file.path(root, "output","final", paste0("results_",Model,cross,".rds"))
out_sample_df <- file.path(root, "output","final", paste0("sample_df_",Model,cross,".qs2"))


###### Save

results <- readRDS(results_file)

saveRDS(results, out_results)

samples_df<- qs_read(sample_df_file )

qs_save(samples_df, out_sample_df)
