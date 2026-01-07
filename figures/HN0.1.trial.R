setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/figures")
eff1<- c(0.2, 0.3, 0.4, 0.5, 0.6)
eff2<- c(0.3, 0.4, 0.5, 0.6, 0.5)
eff3<- c(0.4, 0.5, 0.6, 0.5, 0.4)
eff4<- c(0.5, 0.6, 0.5, 0.4, 0.3)
eff5<- c(0.6, 0.5, 0.4, 0.3, 0.2)
eff6<- rep(0.6, times=5)
eff7<- c(0.5, rep(0.6, times=4))
eff8<- c(0.4, 0.5, rep(0.6, times=3))
eff9<- c(0.3, 0.4, 0.5, rep(0.6, times=2))
truthc<- c(0.05, 0.1, 0.2, 0.28, 0.5)
truthe<- c(0.6, 0.6, 0.6, 0.6, 0.6)
prior.eff.model<- rep(1/9, times=9)
skeletonc<-  c(0.02, 0.10, 0.25, 0.44, 0.62)
skeletone<- matrix(c(eff1, eff2, eff3, eff4, eff5, eff6, eff7, eff8,
                     eff9), nrow=9, byrow = TRUE)
sc<- 1.34
se<- 0.1
targetc<- 0.25
targete<- 0.2
cohortsize<- 3
ncohort<- 39
rar.prop<- 0.25

clin_mat_HN<- matrix(nrow=13, ncol=5)
eff_mat_HN<- matrix(nrow=13, ncol=5)
wk_mat_HN<-  matrix(nrow=13, ncol=9)
expect_zeta_HN<- matrix(nrow=13, ncol=9)
marginale_mat_HN<- matrix(nrow=13, ncol=9)

set.seed(20)
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
  loc<- c(5,4,3,2,1,1,2,3,4)
  eff.m<- which(prior.eff.model==max(prior.eff.model))
  eff.m<-ifelse(length(eff.m)>1, sample(eff.m, size=1, replace=FALSE,prob=rep(1/length(eff.m), times=length(eff.m))), eff.m)
  admiss<- which(skeletonc<=targetc)
  curr<- sample(x=1:max(admiss), size=1, prob=skeletone[1,1:max(admiss)])
  
  ptox.hatc = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  stopc=stope=0; #indicate if trial stops early
  i=1      
  while(((i-1)*cohortsize+1) <= ncohort){
    # generate data for a new cohort of patients
    yc[curr] = yc[curr] + rbinom(1, cohortsize, truthc[curr]);
    ye[curr] = ye[curr] + rbinom(1, cohortsize, truthe[curr]);
    n[curr] = n[curr] + cohortsize;
    
    ##Create admissible doses 
    marginalc = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonc, y=yc,n=n, s2=sc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=-Inf,upper=-Inf, skeletonc, yc, n, s2=sc,abs.tol = 0)$value/marginalc
    ptox.hatc=skeletonc**exp(estc)
    clin_mat_HN[i,]<-ptox.hatc
    #Assess efficacy
    marginale<- sapply(1:nrow(skeletone), function (k) integrate(bcrmh,lower=0,upper=Inf, p=skeletone[k,], y=ye,n=n, s2=se,abs.tol = 0)$value);
    w.k<- sapply(1:nrow(skeletone), function (k) prior.eff.model[k]*marginale[k])/sum(prior.eff.model*marginale)
    wk_mat_HN[i,]<-w.k 
    k_star<- ifelse(sum(w.k==max(w.k))>1,sample(which(w.k==max(w.k)), size=1, replace=FALSE,prob=rep(1/length(which(w.k==max(w.k))), times=length(which(w.k==max(w.k))))), which(w.k==max(w.k)))
    este<- integrate(bcrmht,lower=0,upper=Inf, skeletone[k_star,], ye, n, se,abs.tol = 0)$value/marginale[k_star];
    
    marginale_mat_HN[i,]<- marginale
    expect_zeta_HN[i,]<-sapply(1:9, function (k) integrate(bcrmht,lower=0,upper=Inf, skeletone[k,], ye, n, se,abs.tol = 0)$value/marginale[k])
    
    
    peff.hat=skeletone[k_star,]**exp(este)
    eff_mat_HN[i,]<-peff.hat
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

  # plot(1:6, clin_mat_HN[1,], ylim=c(0,1), type="l")
  # for(i in 2:12){
  #   lines(1:6, clin_mat_HN[i,])
  # }
  # 
  # plot(1:6, eff_mat_HN[1,], ylim=c(0,1), type="l")
  # for(i in 2:12){
  #   lines(1:6, eff_mat_HN[i,])
  # }
  
  data<-tidyr::pivot_longer(data.frame(wk_mat_HN[c(1,3,6,9,13),]), cols = everything(), names_to = "Column", values_to = "Value")
  
  data$Column<- rep(1:9, times=5)
  data$Cohort<- rep(1:5, each=9)
  colnames(data)<- c("Skeleton", "Value", "Cohort")
  data$Cohort<- as.factor(data$Cohort)
  
  
  pdf("HN01wk.pdf", height=4, width=7)
  ggplot(data, aes(x=Skeleton, y=Value, colour=Cohort))+geom_line(size=1.5)+ 
    theme_classic(base_size = 25)+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ labs(x = "Skeleton k", y =expression(w(k ~ "|" ~ D[n])))+
    scale_y_continuous(limits=c(0,1))+scale_x_continuous(breaks=1:11, labels=1:11)+
    scale_color_discrete(name =" ", labels = c("Cohort 1", "Cohort 3", "Cohort 6", "Cohort 9","Cohort 13"))
  dev.off()
  
  
  k<-1
  mat<-sapply(1:9, function (j) skeletone[j,]^exp(expect_zeta_HN[k,j]))
  data_mat<-tidyr::pivot_longer(data.frame(mat), cols = everything(), names_to = "Column", values_to = "Value")
  
  data_mat$Column<- rep(1:9, times=5)
  data_mat$Cohort<- rep(1:5, each=9)
  colnames(data_mat)<- c("Skeleton", "Value", "Dose")
  data_mat$Skeleton<- as.factor(data_mat$Skeleton)
  
  pdf("HN01z.pdf", height=4.5, width=7)
  ggplot(data_mat, aes(x=Dose, y=Value, colour=Skeleton))+geom_line(size=1.5)+ 
    theme_classic(base_size = 25)+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ labs(x = "Dose d", y =expression(z[d*","*k]^exp(hat(zeta))))+
    scale_y_continuous(limits=c(0,1))+scale_x_continuous(breaks=1:6, labels=1:6)+
    scale_color_brewer(palette = "Set3", name ="Skeleton", labels = c("k=1","k=2","k=3","k=4","k=5","k=6","k=7","k=8","k=9",
                                                                        "k=10","k=11"))
  dev.off()
  
  
  