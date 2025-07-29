library(dtpcrm)
library(plyr)
library(tidyverse)
library(mvtnorm)
library(dfcrm)
library(pathviewr)

########################### efftoxcrm #######################################
#Description: Run a single (original) Eff+Tox-CRM trial with log-normal priors

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_e - vector of true efficacy rates for doses
#prior.eff.model - prior probabilities for each efficacy model
#skeletonc - vector of C-DLT skeleton 
#skeletone - vector of P-DLT skeleton 
#sc - numeric variance for C-DLT exponent 
#se - numeric variance for efficacy exponent 
#targetc - numeric target C-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#rar.prop - proportion of trial using response adaptive randomisation 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of responses at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 

efftoxcrm<-function(truthc,truthe, prior.eff.model, skeletonc,skeletone,sc,se,targetc, cohortsize,ncohort,rar.prop){
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
  yc=ye=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
  
  #set first dose 
  loc<- c(6,5,4,3,2,1,1,2,3,4,5)
  eff.m<- which(prior.eff.model==max(prior.eff.model))
  eff.m<-ifelse(length(eff.m)>1, sample(eff.m, size=1, replace=FALSE,prob=rep(1/length(eff.m), times=length(eff.m))), eff.m)
  admiss<- which(skeletonc<=targetc)
  curr<- sample(x=1:max(admiss), size=1, prob=skeletone[1,1:max(admiss)])
  
  ptox.hatc = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  i=1      
  while(((i-1)*cohortsize+1) <= ncohort)
  {
    # generate data for a new cohort of patients
    yc[curr] = yc[curr] + rbinom(1, cohortsize, truthc[curr]);
    ye[curr] = ye[curr] + rbinom(1, cohortsize, truthe[curr]);
    n[curr] = n[curr] + cohortsize;
    
    ##Create admissible doses 
    marginalc = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonc, y=yc,n=n, s2=sc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=-Inf,upper=-Inf, skeletonc, yc, n, s2=sc,abs.tol = 0)$value/marginalc
    ptox.hatc=skeletonc**exp(estc)
    
    #Assess efficacy
    marginale<- sapply(1:nrow(skeletone), function (k) integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletone[k,], y=ye,n=n, s2=se,abs.tol = 0)$value);
    w.k<- sapply(1:nrow(skeletone), function (k) prior.eff.model[k]*marginale[k])/sum(prior.eff.model*marginale)
    k_star<- ifelse(sum(w.k==max(w.k))>1,sample(which(w.k==max(w.k)), size=1, replace=FALSE,prob=rep(1/length(which(w.k==max(w.k))), times=length(which(w.k==max(w.k))))), which(w.k==max(w.k)))
    este<- integrate(bcrmht,lower=-Inf,upper=Inf, skeletone[k_star,], ye, n, se,abs.tol = 0)$value/marginale[k_star];
    peff.hat=skeletone[k_star,]**exp(este)
    if(sum(n)<=rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- ifelse(length(curr_tox)>1, sample(curr_tox, size=1,prob=peff.hat[curr_tox]), max(curr_tox))
    }
    if(sum(n)>rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- which.max(peff.hat[curr_tox])
    }
    i=i+1
  }
  return(list(dose.select=dose.select,tox.datac=yc,eff.data=ye,pt.allocation=n))
}

########################### efftoxcrm.sim #######################################
#Description: Run a simulation study for (original) Eff+Tox-CRM trial with log-normal priors

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_e - vector of true efficacy rates for doses
#prior.eff.model - prior probabilities for each efficacy model
#skeletonc - vector of C-DLT skeleton 
#skeletone - vector of P-DLT skeleton 
#sc - numeric variance for C-DLT exponent 
#se - numeric variance for efficacy exponent 
#targetc - numeric target C-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#rar.prop - proportion of trial using response adaptive randomisation 
#ntrial - number of trial simulations to run 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of responses at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 

efftoxcrm.sim<-function(truthc,truthe, prior.eff.model, skeletonc,skeletone,sc,se,targetc,cohortsize,ncohort,rar.prop, ntrial){
  ndose=length(truthc)
  dose.select<-yc<-ye<-n<-matrix(nrow=ntrial,ncol=ndose)
  for(i in 1:ntrial){
    result<-efftoxcrm(truthc, truthe,prior.eff.model, skeletonc, skeletone, sc, se, targetc, cohortsize, ncohort, rar.prop)
    dose.select[i,]=result$dose.select
    n[i,]<- c(result$pt.allocation)
    yc[i,]<-  c(result$tox.datac)
    ye[i,]<-  c(result$eff.data)
  }
  return(list(colMeans(dose.select), n, yc, ye))
}

########################### efftoxcrm #######################################
#Description: Run a single  Eff+Tox-CRM trial with half-normal prior for \zeta prior and log-normal priors for C-DLT

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_e - vector of true efficacy rates for doses
#prior.eff.model - prior probabilities for each efficacy model
#skeletonc - vector of C-DLT skeleton 
#skeletone - vector of P-DLT skeleton 
#sc - numeric variance for C-DLT exponent 
#se - numeric variance for efficacy exponent 
#targetc - numeric target C-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#rar.prop - proportion of trial using response adaptive randomisation 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of responses at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 

efftoxcrm.HN<-function(truthc,truthe, prior.eff.model, skeletonc,skeletone,sc,se,targetc, cohortsize,ncohort,rar.prop){
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
  yc=ye=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
  
  #set first dose 
  loc<- c(6,5,4,3,2,1,1,2,3,4,5)
  eff.m<- which(prior.eff.model==max(prior.eff.model))
  eff.m<-ifelse(length(eff.m)>1, sample(eff.m, size=1, replace=FALSE,prob=rep(1/length(eff.m), times=length(eff.m))), eff.m)
  admiss<- which(skeletonc<=targetc)
  curr<- sample(x=1:max(admiss), size=1, prob=skeletone[1,1:max(admiss)])
  
  ptox.hatc = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  i=1      
  while(((i-1)*cohortsize+1) <= ncohort)
  {
    # generate data for a new cohort of patients
    yc[curr] = yc[curr] + rbinom(1, cohortsize, truthc[curr]);
    ye[curr] = ye[curr] + rbinom(1, cohortsize, truthe[curr]);
    n[curr] = n[curr] + cohortsize;
    
    ##Create admissible doses 
    marginalc = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonc, y=yc,n=n, s2=sc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=-Inf,upper=-Inf, skeletonc, yc, n, s2=sc,abs.tol = 0)$value/marginalc
    ptox.hatc=skeletonc**exp(estc)
    
    #Assess efficacy
    marginale<- sapply(1:nrow(skeletone), function (k) integrate(bcrmh,lower=0,upper=Inf, p=skeletone[k,], y=ye,n=n, s2=se,abs.tol = 0)$value);
    w.k<- sapply(1:nrow(skeletone), function (k) prior.eff.model[k]*marginale[k])/sum(prior.eff.model*marginale)
    k_star<- ifelse(sum(w.k==max(w.k))>1,sample(which(w.k==max(w.k)), size=1, replace=FALSE,prob=rep(1/length(which(w.k==max(w.k))), times=length(which(w.k==max(w.k))))), which(w.k==max(w.k)))
    este<- integrate(bcrmht,lower=0,upper=Inf, skeletone[k_star,], ye, n, se,abs.tol = 0)$value/marginale[k_star];
    peff.hat=skeletone[k_star,]**exp(este)
    if(sum(n)<=rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- ifelse(length(curr_tox)>1, sample(curr_tox, size=1,prob=peff.hat[curr_tox]), max(curr_tox))
    }
    if(sum(n)>rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- which.max(peff.hat[curr_tox])
    }
    i=i+1
  }
  return(list(dose.select=dose.select,tox.datac=yc,eff.data=ye,pt.allocation=n))
}

########################### efftoxcrm.sim.HN #######################################
#Description: Run a simulation study for Eff+Tox-CRM trial with half-normal prior for \zeta prior

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_p - vector of true efficacy rates for doses
#prior.eff.model - prior probabilities for each efficacy model
#skeletonc - vector of C-DLT skeleton 
#skeletone - vector of P-DLT skeleton 
#sc - numeric variance for C-DLT exponent 
#se - numeric variance for efficacy exponent 
#targetc - numeric target C-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#rar.prop - proportion of trial using response adaptive randomisation 
#ntrial - number of trial simulations to run 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of responses at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 
efftoxcrm.sim.HN<-function(truthc,truthe, prior.eff.model, skeletonc,skeletone,sc,se,targetc, cohortsize,ncohort,rar.prop, ntrial){
  ndose=length(truthc)
  dose.select<-yc<-ye<-n<-matrix(nrow=ntrial,ncol=ndose)
  for(i in 1:ntrial){
    result<-efftoxcrm.HN(truthc, truthe,prior.eff.model, skeletonc, skeletone, sc, se, targetc, cohortsize, ncohort, rar.prop)
    dose.select[i,]=result$dose.select
    n[i,]<- c(result$pt.allocation)
    yc[i,]<-  c(result$tox.datac)
    ye[i,]<-  c(result$eff.data)
  }
  return(list(colMeans(dose.select), n, yc, ye))
}

########################### mu_kl #######################################
#Description: estimated probability a dose is the optimal dose for joint Eff+Tox-CRM dose recommendation 

#Input: 
#mtd - dose of interest 
#prob_l - probability each dose is the most effective 
#beta_shape - shape hyperparameter for C-DLT parameter 
#beta_rate - rate hyperparameter for C-DLT parameter 
#no.dose - number of doses under investigation 
#Hset_beta - matrix of the H-sets for C-DLT parameter

#Output: numeric probability specified dose is the optimal dose

mu_kl<-function(mtd, prob_l, beta_shape, beta_rate, no.doses, Hset_beta){
  if(mtd==1){
    return(1-(1-pgamma(Hset_beta[mtd+1,1], shape=beta_shape, rate=beta_rate))*sum(prob_l[(mtd+1):length(prob_l)]))
  }
  if(mtd==no.doses){
    return((1-pgamma(Hset_beta[mtd,1],shape=beta_shape, rate=beta_rate))*prob_l[mtd])
  }
  if(mtd>1 && mtd <no.doses){
    return(((1-pgamma(Hset_beta[mtd,1],shape=beta_shape, rate=beta_rate))*sum(prob_l[mtd:length(prob_l)])) - 
             ((1-pgamma(Hset_beta[mtd+1,1], shape=beta_shape, rate=beta_rate))*sum(prob_l[(mtd+1):length(prob_l)])))
  }
}

########################### KLj #######################################
#Description: KL divergence of a priori optimal dose probabilities and discrete uniform distribution for joint Eff+Tox-CRM dose recommendation 

#Input: 
#vec - vector containing C-DLT shape and rate parameters respectively. 
#no.dose - number of doses under investigation 
#Hset_b - matrix of the H-sets for C-DLT parameter
#prob_l - probability each dose is the most effective 

#Output: numeric KL divergence 


KLj<- function(vec, no.d, Hset_b, prob_l){
  b_shape<-vec[1]
  b_rate<- vec[2]
  kl<- sum(sapply(1:no.d, function (k) mu_kl(k, prob_l, b_shape, b_rate, no.d, Hset_b))*
             log(no.d*sapply(1:no.d, function (k) mu_kl(k, prob_l, b_shape, b_rate, no.d, Hset_b))))
  return(kl)
}

########################### KLja #######################################
#Description: KL divergence of a priori optimal dose probabilities and discrete uniform distribution for joint Eff+Tox-CRM dose recommendation when 
# prob_l is determined by its own shape and rate parameters  

#Input: 
#vec - vector containing C-DLT shape and rate parameters respectively. 
#l_shape - shape parameter for prob_l 
#l_rate - rate parameter for prob_l 
#no.dose - number of doses under investigation 
#Hset_b - matrix of the H-sets for C-DLT parameter

#Output: numeric KL divergence 

KLja<- function(vec, l_shape,l_rate, no.d, Hset_b){
  prob_l<- c(pbeta(1/6, l_shape, l_rate),pbeta(2/6, l_shape, l_rate)-pbeta(1/6, l_shape, l_rate),pbeta(3/6, l_shape, l_rate)-pbeta(2/6, l_shape, l_rate),
             pbeta(4/6, l_shape, l_rate)-pbeta(3/6, l_shape, l_rate),pbeta(5/6, l_shape, l_rate)-pbeta(4/6, l_shape, l_rate),1-pbeta(5/6, l_shape, l_rate))
  kl<- KLj(vec, no.d, Hset_b, prob_l)
  return(kl)
}

########################### mu_marginal #######################################
#Description: estimated probability a dose is the MTD for Eff+Tox-CRM dose recommendation

#Input: 
#mtd - dose of interest 
#shape - shape hyperparameter for C-DLT parameter 
#rate - rate hyperparameter for C-DLT parameter 
#no.dose - number of doses under investigation 
#Hset - matrix of the H-sets for C-DLT parameter

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
#Description: KL divergence of a priori MTD probabilities and discrete uniform distribution 

#Input: 
#shape - shape hyperparameter for C-DLT parameter 
#rate - rate hyperparameter for C-DLT parameter 
#no.d - number of doses under investigation 
#Hset - matrix of the H-sets for C-DLT/P-DLT parameter

#Output: numeric KL divergence 

KL_marginal<- function(shape, rate, no.d, Hset){
  kl<-sum(sapply(1:5, function(k) mu_marginal(k,shape, rate, no.d, Hset)*log(no.d*mu_marginal(k,shape, rate, no.d, Hset))))
  return(kl)
}

########################### mu_marginal_h #######################################
#Description: estimated probability a dose is the most effective dose for Eff+Tox-CRM dose recommendation

#Input: 
#mtd - dose of interest 
#shape - shape hyperparameter for prob_l 
#rate - rate hyperparameter for prob_l
#no.dose - number of doses under investigation 

#Output: numeric probability specified dose is the most effective dose

mu_marginal_h<- function(mtd, shape, rate, no.d){
  s<- seq(from=0, to =1, length.out=no.d+1)[-1]
  if(mtd==1){
    return(pbeta(s[1], shape, rate))
  }
  if(mtd==no.d){
    return(1-pbeta(s[length(s)-1], shape, rate))
  }
  if(mtd>1 && mtd < no.d){
    return(pbeta(s[mtd], shape, rate)-pbeta(s[mtd-1], shape, rate))
  }
}

########################### KL_marginal #######################################
#Description: KL divergence of a priori probabiltiy of most effective dose and discrete uniform distribution 

#Input: 
#shape - shape hyperparameter for prob_l
#rate - rate hyperparameter for prob_l
#no.d - number of doses under investigation

#Output: numeric KL divergence 


KL_marginal_h<- function(shape, rate, no.d){
  kl<-sum(sapply(1:no.d, function(k) mu_marginal_h(k,shape, rate, no.d)*log(no.d*mu_marginal_h(k,shape, rate, no.d))))
  return(kl)
}

########################### efftoxcrm.uninf #######################################
#Description: Run a single modified Eff+Tox-CRM trial with half-normal prior for \zeta prior and gamma priors for C-DLT

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_e - vector of true efficacy rates for doses
#prior.eff.model - prior probabilities for each efficacy model
#skeletonc - vector of C-DLT skeleton 
#skeletone - vector of P-DLT skeleton 
#ac - numeric shape for C-DLT parameter  
#bc - numeric rate for C-DLT parameter
#theta_s2 - variance for \zeta parameter
#targetc - numeric target C-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#rar.prop - proportion of trial using response adaptive randomisation 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of responses at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 

efftoxcrm.uninf<-function(truthc,truthe, prior.eff.model, skeletonc,skeletone,ac, bc, theta_s2, targetc, cohortsize,ncohort,rar.prop){
  bcrmh<-function(x,p,y,n,a,b){
    lik=(b^a)*(x^(a-1))*exp(-b*x)
    for(j in 1:length(p)){
      pj=p[j]**x
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  bcrmh_uninf<-function(a,p,y,n,s2){
    lik=exp(-0.5*a*a/s2)
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      #pj=p[j]**a
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
  
  bcrmht_uninf<-function(a,p,y,n,s2){
    lik=a*exp(-0.5*a*a/s2)
    for(j in 1:length(p)){
      pj=p[j]**exp(a)
      #pj=p[j]**a
      lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik);
  }
  
  ###simulate a trial       
  ndose = length(skeletonc);   #number of combos
  yc=ye=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
  
  #set first dose 
  #eff.m<-sample(x=1:11, size=1, prob=prior.eff.model)
  loc<- c(6,5,4,3,2,1,1,2,3,4,5)
  eff.m<-  which(prior.eff.model==max(prior.eff.model))
  eff.m<-ifelse(length(eff.m)>1, sample(eff.m, size=1, replace=FALSE,prob=rep(1/length(eff.m), times=length(eff.m))), eff.m)
  admiss<- which(skeletonc<=targetc)
  curr<- sample(x=1:max(admiss), size=1, prob=skeletone[1,1:max(admiss)])
  
  ptox.hatc = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  stopc=stope=0; #indicate if trial stops early
  i=1      
  while(((i-1)*cohortsize+1) <= ncohort){
    # generate data for a new cohort of patients
    yc[curr] = yc[curr] + rbinom(1, cohortsize, as.numeric(truthc[curr]));
    ye[curr] = ye[curr] + rbinom(1, cohortsize, as.numeric(truthe[curr]));
    n[curr] = n[curr] + cohortsize;
    
    ##Create admissible doses 
    marginalc = integrate(bcrmh,lower=0,upper=Inf, p=skeletonc, y=yc,n=n, a=ac, b=bc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=0,upper=Inf, skeletonc, yc, n, a=ac, b=bc,abs.tol = 0)$value/marginalc
    ptox.hatc=skeletonc**estc
    
    #Assess efficacy 
    marginale<- sapply(1:nrow(skeletone), function (k) integrate(bcrmh_uninf,lower=0,upper=Inf, p=skeletone[k,], y=ye,n=n, s2=theta_s2,abs.tol = 0)$value);
    w.k<- sapply(1:nrow(skeletone), function (k) prior.eff.model[k]*marginale[k])/sum(prior.eff.model*marginale)
    k_star<- ifelse(sum(w.k==max(w.k))>1,sample(which(w.k==max(w.k)), size=1, replace=FALSE,prob=rep(1/length(which(w.k==max(w.k))), times=length(which(w.k==max(w.k))))), which(w.k==max(w.k)))
    este<- integrate(bcrmht_uninf,lower=0,upper=Inf, skeletone[k_star,], ye, n, s2=theta_s2,abs.tol = 0)$value/marginale[k_star];
    peff.hat=skeletone[k_star,]**exp(este)
    #peff.hat=skeletone[k_star,]**este
    if(sum(n)<=rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- ifelse(length(curr_tox)>1, sample(curr_tox, size=1,prob=peff.hat[curr_tox]), max(curr_tox))
    }
    if(sum(n)>rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- which.max(peff.hat[curr_tox])
    }
    i=i+1
  }
  if(stopc==0 & stope==0){
    dose.select[curr]=dose.select[curr]+1
  }
  return(list(dose.select=dose.select,tox.datac=yc,eff.data=ye,pt.allocation=n,stopc=stopc,stope=stope))
}

########################### efftoxcrm.uninf.sim #######################################
#Description: Run a simulation study for the modified Eff+Tox-CRM trial with half-normal prior for \zeta prior and gamma priors for C-DLT

#Input: 
#truth_c - vector of true C-DLT rates for doses
#truth_e - vector of true efficacy rates for doses
#prior.eff.model - prior probabilities for each efficacy model
#skeletonc - vector of C-DLT skeleton 
#skeletone - vector of P-DLT skeleton 
#ac - numeric shape for C-DLT parameter  
#bc - numeric rate for C-DLT parameter
#theta_s2 - variance for \zeta parameter
#targetc - numeric target C-DLT rate 
#cohortsize - numeric total sample size for trial
#ncohort - number of patients in each cohort 
#rar.prop - proportion of trial using response adaptive randomisation 
#ntrial - number of trial simulations to run 

#Output: Vector; final dose selection, number of C-DLTs at each dose, number of responses at each dose, patient allocation and indicator for stopping due to C-DLT or P-DLTs. 


efftoxcrm.unif.sim<-function(truthc,truthe, prior.eff.model, skeletonc,skeletone,ac, bc, theta_s2, targetc,cohortsize,ncohort,rar.prop, ntrial){
  ndose=length(truthc)
  dose.select<-yc<-ye<-n<-matrix(nrow=ntrial,ncol=ndose)
  for(i in 1:ntrial){
    result<-efftoxcrm.uninf(truthc, truthe,prior.eff.model, skeletonc, skeletone, ac, bc, theta_s2, targetc, cohortsize, ncohort, rar.prop)
    dose.select[i,]=result$dose.select
    n[i,]<- c(result$pt.allocation)
    yc[i,]<-  c(result$tox.datac)
    ye[i,]<-  c(result$eff.data)
  }
  return(list(colMeans(dose.select), n, yc, ye))
}

