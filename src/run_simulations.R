run_simulations<-function(nsamples,M,filter_obj,transformfunc){
  
  
  single_run <- function(single_par,func,footransform) {
    
    p<-footransform(single_par)   
    func$run(p, save_history = TRUE)
    
    return(sims<-func$history())
    
  }
  
  
  
  
  # run the model for each set of parameters in the MCMC
  tmp<-single_run(M[1,],filter_obj,transformfunc)
  nst<-dim(tmp)[1]
  nt<-dim(tmp)[3]
  idx<- seq(1, nst, 1)
  names(idx)<- paste(row.names(tmp))
  idx<-as.list(idx)
  
  runs<- array(0, c(nst, nsamples, nt));  
  
  pb = txtProgressBar(min = 0, max = (dim(M)[1]), initial = 0) 
  
  for (ii in 1:(dim(M)[1])){
    setTxtProgressBar(pb,ii)
    runs[,ii,]<-single_run(M[ii,],filter2,footransform)
    
  }
  
  
  out<-list(
    idx=idx,
    runs=runs,
    pars=M)
  
  return(out)
  
  
}