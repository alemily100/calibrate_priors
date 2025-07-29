library(nnet)
library(binom)
library(dfcrm)
library(tidyverse)
library(lattice)

#################

###########################time_to_dlt#######################################
#Description: Simulates correlated C-DLT and P-DLT as per the Clayton model

#Input: 
#dose - dose for patient i
#clin_tox - vector of true probability of C-DLT for each of the 5 dosages 
#pat_tox - vector of true probability of P-DLT for each of the 5 dosages 
#phi - correlation parameter

#Output: Vector; time to C-DLT and P-DLT for patient i 

time_to_dlt<- function(dose,clin_tox, pat_tox, phi){
  lambda_c<- -log(1-clin_tox[dose])
  lambda_p<- -log(1-pat_tox[dose])
  u1<-runif(1)
  time_c<- (-log(u1))/lambda_c
  u2<-runif(1)
  a<- u2/(u1^(-(phi+1)/phi))
  b<- u1^(-1/phi)-1
  time_p<- (phi/lambda_p)*log(a^(1/(-phi-1))-b)
  times<- c(time_c, time_p)
  return(times)
}

########################### pro_crm #######################################
#Description: Run a single (original) PRO-CRM trial with log-normal priors

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_p - vector of true P-DLT rates for doses
#skeletonc - vector of C-DLT skeleton 
#skeletonp - vector of P-DLT skeleton 
#sc - numeric variance for C-DLT exponent 
#sp - numeric variance for P-DLT exponent 
#targetc - numeric target C-DLT rate 
#targetp - numeric target P-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#n.stop - Number of patients needed on one combination to stop the trial
#start - starting dose 
#cl - confidence level for stopping rules
#phi - correlation paramter for clayton copula between C-DLT and P-DLT observations 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of P-DLTs at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 

procrm<-function(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl, phi){
  
  bcrmh<-function(a,p,y,n,s2){
    lik=exp(-0.5*a*a/s2)
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  bcrmht<-function(a,p,y,n,s2){
    lik=a*exp(-0.5*a*a/s2)
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  ###simulate a trial       
  ndose = length(skeletonc);   #number of combos
  yc=yp=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
  curr = start;  # current dose level           
  ptox.hatc = ptox.hatp = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  stopc=stopp=0; #indicate if trial stops early
  i=1      
  while(((i-1)*cohortsize+1) <= ncohort)
  {
    # generate data for a new cohort of patients
    clayton<- sapply(rep(curr, times =cohortsize), function (k) time_to_dlt(k, truthc, truthp, phi));
    yc[curr] = yc[curr] + sum(as.numeric(clayton[1,]<=1));
    yp[curr] = yp[curr] + sum(as.numeric(clayton[2,]<=1));
    n[curr] = n[curr] + cohortsize;
    
    if(any(n>n.stop)){
      stop<-0
      break
    }
    ##Model-based estimation of DLT probabilities
    marginalc = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonc, y=yc,n=n, s2=sc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=-10,upper=10, skeletonc, yc, n, sc,abs.tol = 0)$value/marginalc
    #ptox.hatc=sapply(1:5, function (k) integrate(bcrmht,lower=-10,upper=10,i=k, skeletonc, yc, n, sc,abs.tol = 0)$value)/marginalc
    
    marginalp = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonp, y=yp,n=n,s2=sp,abs.tol = 0)$value;
    estp=integrate(bcrmht,lower=-10,upper=10, skeletonp, yp, n, sp, abs.tol = 0)$value/marginalp
    #ptox.hatp=sapply(1:5, function (k) integrate(bcrmht,lower=-10,upper=10, i=k,skeletonp, yp, n, sp, abs.tol = 0)$value)/marginalp
    
    ptox.hatc=skeletonc**exp(estc)
    ptox.hatp=skeletonp**exp(estp)
    
    #########stopping rules
    safetyc=binom.confint(yc[1],n[1],conf.level=cl,methods="agresti-coull")$lower
    if(safetyc>targetc){
      stopc=1
      break
    }
    
    safetyp=binom.confint(yp[1],n[1],conf.level=cl,methods="agresti-coull")$lower
    if(safetyp>targetp){
      stopp=1
      break
    }
    
    ##Allocation algorithm
    distancec=abs(ptox.hatc-targetc)
    distancep=abs(ptox.hatp-targetp)
    
    bestc=which.is.max(-distancec)
    bestp=which.is.max(-distancep)
    best = min(bestc,bestp)
    curr=min(best,curr+1)
    
    i=i+1
  }
  if(stopc==0 & stopp==0){
    dose.select[curr]=dose.select[curr]+1;
  }
  return(list(dose.select=dose.select,tox.datac=yc,tox.datap=yp,pt.allocation=n,stopc=stopc,stopp=stopp))
}

########################### procrm.sim #######################################
#Description: Run a simulation study for (original) PRO-CRM trial design with log-normal priors

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_p - vector of true P-DLT rates for doses
#skeletonc - vector of C-DLT skeleton 
#skeletonp - vector of P-DLT skeleton 
#sc - numeric variance for C-DLT exponent 
#sp - numeric variance for P-DLT exponent 
#targetc - numeric target C-DLT rate 
#targetp - numeric target P-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#n.stop - Number of patients needed on one combination to stop the trial
#start - starting dose 
#cl - confidence level for stopping rules
#ntrial - number of trial simulations to run 
#phi - correlation paramter for clayton copula between C-DLT and P-DLT observations 

#Output: list; final dose selections, matrix of C-DLTs at each dose for each trial , number of P-DLTs at each dose for each trial, patient allocation at each dose for each trial. 


procrm.sim<-function(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl,ntrial, phi){
  ndose=length(truthc)
  
  dose.select<-n<-matrix(nrow=ntrial,ncol=ndose)
  yc<-yp<- c()
  nstopc=nstopp=0
  n<- matrix(nrow=ntrial, ncol=ndose)
  for(i in 1:ntrial){
    result<-procrm(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl, phi)
    dose.select[i,]=result$dose.select
    yc[i]=sum(result$tox.datac)
    yp[i]=sum(result$tox.datap)
    n[i,]<- c(result$pt.allocation)
  }
  
  return(list(colMeans(dose.select), yc,yp, n))
}

########################### mu_k #######################################
#Description: estimated probability a dose is the MTD for joint PRO-CRM dose recommendation 

#Input: 
#mtd - dose of interest 
#beta_shape - shape hyperparameter for C-DLT parameter 
#gamma_shape - shape hyperparameter for P-DLT parameter 
#beta_rate - rate hyperparameter for C-DLT parameter 
#gamma_rate - rate hyperparameter for P-DLT parameter 
#no.dose - number of doses under investigation 
#Hset_beta - matrix of the H-sets for C-DLT parameter
#Hset_gamma - matrix of the H-sets for P-DLT parameter

#Output: numeric probability specified dose is the MTD

mu_k<-function(mtd, beta_shape, gamma_shape, beta_rate, gamma_rate, no.doses, Hset_beta, Hset_gamma){
  if(mtd==1){
    return(1-(1-pgamma(Hset_beta[mtd+1,1], shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd+1,1], shape=gamma_shape, rate=gamma_rate)))
  }
  if(mtd==no.doses){
    return((1-pgamma(Hset_beta[mtd,1],shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd,1], shape=gamma_shape, rate=gamma_rate)))
  }
  if(mtd>1 && mtd < no.doses){
    return((1-pgamma(Hset_beta[mtd,1],shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd,1], shape=gamma_shape, rate=gamma_rate)) -
             (1-pgamma(Hset_beta[mtd+1,1], shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd+1,1], shape=gamma_shape, rate=gamma_rate)))
  }
}


########################### KL #######################################
#Description: KL divergence of a priori MTD probabilities and discrete uniform distribution for joint PRO-CRM dose recommendation 

#Input: 
#vec - vector containing C-DLT and P-DLT shape parameters and C-DLT and P-DLT rate parameters respectively. 
#no.dose - number of doses under investigation 
#Hset_b - matrix of the H-sets for C-DLT parameter
#Hset_g - matrix of the H-sets for P-DLT parameter

#Output: numeric KL divergence 

KL<- function(vec, no.d, Hset_b, Hset_g){
  b_shape<-vec[1]
  g_shape<- vec[2]
  b_rate<- vec[3]
  g_rate<- vec[4]
  kl<-sum(sapply(1:5, function(k) mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g)*log(no.d*mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g))))
  return(kl)
}

########################### mu_marginal #######################################
#Description: estimated probability a dose is the MTD for marginal PRO-CRM dose recommendation (ie C-MTD or P-MTD)

#Input: 
#mtd - dose of interest 
#shape - shape hyperparameter for C-DLT/P-DLT parameter 
#rate - rate hyperparameter for C-DLT/P-DLT parameter 
#no.dose - number of doses under investigation 
#Hset - matrix of the H-sets for C-DLT/P-DLT parameter

#Output: numeric probability specified dose is the MTD

mu_marginal<-function(mtd, shape, rate, no.d, Hset){
  if(mtd==1){
    return(pgamma(Hset[mtd+1,1], shape=shape, rate=rate))
  }
  if(mtd==no.d){
    return(1-pgamma(Hset[mtd,1], shape=shape, rate=rate))
  }
  if(mtd>1 && mtd < no.d){
    return(pgamma(Hset[mtd,2], shape=shape, rate=rate)-pgamma(Hset[mtd, 1], shape=shape, rate=rate))
  }
}

########################### KL_marginal #######################################
#Description: KL divergence of a priori MTD probabilities and discrete uniform distribution for marginal PRO-CRM dose recommendation (ie C-MTD or P-MTD)

#Input: 
#shape - shape hyperparameter for C-DLT/P-DLT parameter 
#rate - rate hyperparameter for C-DLT/P-DLT parameter 
#no.d - number of doses under investigation 
#Hset - matrix of the H-sets for C-DLT/P-DLT parameter

#Output: numeric KL divergence 

KL_marginal<- function(shape, rate, no.d, Hset){
  kl<-sum(sapply(1:5, function(k) mu_marginal(k,shape, rate, no.d, Hset)*log(no.d*mu_marginal(k,shape, rate, no.d, Hset))))
  return(kl)
}

########################### procrm_gamma #######################################
#Description: Run a single modified PRO-CRM trial with gamma priors 

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_p - vector of true P-DLT rates for doses
#skeletonc - vector of C-DLT skeleton 
#skeletonp - vector of P-DLT skeleton 
#ac - numeric shape for C-DLT parameter  
#bc - numeric rate for C-DLT parameter
#ap - numeric shape for P-DLT parameter
#bp - numeric rate for P-DLT parameter
#targetc - numeric target C-DLT rate 
#targetp - numeric target P-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#n.stop - Number of patients needed on one combination to stop the trial
#start - starting dose 
#cl - confidence level for stopping rules
#phi - correlation paramter for clayton copula between C-DLT and P-DLT observations 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of P-DLTs at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 


procrm_gamma<-function(truthc,truthp,skeletonc,skeletonp,ac, bc, ap, bp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl, phi){
  bcrmh<-function(x,p,y,n,a,b){
    lik=(b^a)*(x^(a-1))*exp(-b*x)
    for(j in 1:length(p)){
      pj=p[j]**x
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  bcrmht<-function(x,p,y,n,a,b){
    lik=x*(b^a)*(x^(a-1))*exp(-b*x)
    for(j in 1:length(p)){
      pj=p[j]**x
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  ###simulate a trial
  ndose = length(skeletonc);   #number of combos
  yc=yp=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
  curr = start;  # current dose level
  ptox.hatc = ptox.hatp = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  stopc=stopp=0; #indicate if trial stops early
  i=1  
  while(((i-1)*cohortsize+1) <= ncohort)
  {
    # generate data for a new cohort of patients
    clayton<- sapply(rep(curr, times =cohortsize), function (k) time_to_dlt(k, truthc, truthp, phi));
    yc[curr] = yc[curr] + sum(as.numeric(clayton[1,]<=1));
    yp[curr] = yp[curr] + sum(as.numeric(clayton[2,]<=1));
    n[curr] = n[curr] + cohortsize;
    
    if(any(n>n.stop)){
      stop<-0
      break
    }
    ##Model-based estimation of DLT probabilities
    marginalc = integrate(bcrmh,lower=0,upper=Inf, p=skeletonc, y=yc,n=n, a=ac, b=bc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=0,upper=10, p=skeletonc, y=yc, n=n, a=ac, b=bc,abs.tol = 0)$value/marginalc
    #ptox.hatc=sapply(1:5, function (k) integrate(bcrmht,lower=0,upper=10, i=k, p=skeletonc, y=yc, n=n, a=ac, b=bc,abs.tol = 0)$value)/marginalc
    
    marginalp = integrate(bcrmh,lower=0,upper=Inf, p=skeletonp, y=yp,n=n,a=ap, b=bp,abs.tol = 0)$value;
    estp=integrate(bcrmht,lower=0,upper=10, p=skeletonp, y=yp, n=n, a=ap, b=bp, abs.tol = 0)$value/marginalp
    #ptox.hatp=sapply(1:5, function (k) integrate(bcrmht,lower=0,upper=10, i=k, p=skeletonp, y=yp, n=n, a=ap, b=bp, abs.tol = 0)$value)/marginalp
    
    ptox.hatc=skeletonc**estc
    ptox.hatp=skeletonp**estp
    #########stopping rules
    safetyc=binom.confint(yc[1],n[1],conf.level=cl,methods="agresti-coull")$lower
    if(safetyc>targetc){
      stopc=1
      break
    }
    
    safetyp=binom.confint(yp[1],n[1],conf.level=cl,methods="agresti-coull")$lower
    if(safetyp>targetp){
      stopp=1
      break
    }
    
    ##Allocation algorithm
    distancec=abs(ptox.hatc-targetc)
    distancep=abs(ptox.hatp-targetp)
    
    bestc=which.is.max(-distancec)
    bestp=which.is.max(-distancep)
    best = min(bestc,bestp)
    curr=min(best,curr+1)
    
    i=i+1
  }
  if(stopc==0 & stopp==0){
    dose.select[curr]=dose.select[curr]+1;
  }
  return(list(dose.select=dose.select,tox.datac=yc,tox.datap=yp,pt.allocation=n,stopc=stopc,stopp=stopp))
}

########################### procrm.sim_gamma #######################################
#Description: Run a simulation study for modified PRO-CRM trial design with gamma priors

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_p - vector of true P-DLT rates for doses
#skeletonc - vector of C-DLT skeleton 
#skeletonp - vector of P-DLT skeleton 
#ac - numeric shape for C-DLT parameter  
#bc - numeric rate for C-DLT parameter
#ap - numeric shape for P-DLT parameter
#bp - numeric rate for P-DLT parameter
#targetc - numeric target C-DLT rate 
#targetp - numeric target P-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#n.stop - Number of patients needed on one combination to stop the trial
#start - starting dose 
#cl - confidence level for stopping rules
#ntrial - number of trial simulations to run 
#phi - correlation paramter for clayton copula between C-DLT and P-DLT observations 

#Output: list; final dose selections, matrix of C-DLTs at each dose for each trial , number of P-DLTs at each dose for each trial, patient allocation at each dose for each trial. 



procrm.sim_gamma<-function(truthc,truthp,skeletonc,skeletonp,ac,bc,ap,bp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl,ntrial, phi){
  ndose=length(truthc)
  
  dose.select<-n<-matrix(nrow=ntrial,ncol=ndose)
  n<- matrix(nrow=ntrial, ncol=ndose)
  yc<-yp<- c()
  nstopc=nstopp=0
  
  for(i in 1:ntrial){
    result<-procrm_gamma(truthc,truthp,skeletonc,skeletonp,ac,bc,ap,bp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl, phi)
    dose.select[i,]=result$dose.select
    n[i,]<- c(result$pt.allocation)
    yc[i]=sum(result$tox.datac)
    yp[i]=sum(result$tox.datap)
  }
  return(list(colMeans(dose.select), yc,yp, n))
}

