setwd("/annotated")
source("efftox_functions.R")
library(clusterGeneration)
library(parallel)

skeletonc<-  c(0.01, 0.08, 0.15, 0.22, 0.29, 0.36)
eff1<- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
eff2<- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.5)
eff3<- c(0.3, 0.4, 0.5, 0.6, 0.5, 0.4)
eff4<- c(0.4, 0.5, 0.6, 0.5, 0.4, 0.3)
eff5<- c(0.5, 0.6, 0.5, 0.4, 0.3, 0.2)
eff6<- c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
eff7<- rep(0.6, times=6)
eff8<- c(0.5, rep(0.6, times=5))
eff9<- c(0.4, 0.5, rep(0.6, times=4))
eff10<- c(0.3, 0.4, 0.5, rep(0.6, times=3))
eff11<- c(0.2, 0.3, 0.4, 0.5, rep(0.6, times=2))

############################################# Marginal prior calibration #######################################################

prior.eff.model<- rep(1/11, times=11)
sc<- 1.34
se<- 1.34

all_sim<- NULL
set.seed(1001)
ntrial=5000
for(k in 1:16){
  i<- sim_scenario_k[1,k]
  j<- sim_scenario_k[2,k]
  clin_true<- as.numeric(true_tox[i,])
  eff_true<- as.numeric(true_eff[j,])
  ndose=length(clin_true)
  dose.select<-yc<-ye<-n<-matrix(nrow=ntrial,ncol=ndose)
  cl<- makeCluster(detectCores())
  clusterSetRNGStream(cl,1001)
  invisible(clusterEvalQ(cl,{
    library(clusterGeneration)
    library(parallel)
    setwd("/annotated")
    source("efftox_functions.R")  
  }))
  clusterExport(cl, c("clin_true","eff_true", "prior.eff.model", "skeletonc","mat_skeletone","sc", "se","targetc","cohortsize","rar.prop", "ntrial"))
  result<-parLapply(cl, 1:ntrial, function (k) efftoxcrm(clin_true,eff_true, prior.eff.model, skeletonc,mat_skeletone,sc, se,targetc,cohortsize,48,rar.prop))
  stopCluster(cl)
  for(m in 1:ntrial){
    dose.select[m,]<-result[[m]]$dose.select
    n[m,]<- c(result[[m]]$pt.allocation)
    yc[m,]<-  c(result[[m]]$tox.datac)
    ye[m,]<-  c(result[[m]]$eff.data)
  }
  sim<-rbind.data.frame(clin_true, eff_true, colMeans(dose.select))
  dimnames(sim)[[2]]<- c( "1","2", "3", "4", "5", "6")
  dimnames(sim)[[1]]<- c("Probability of clinician DLT", "Probability of efficacy", "Probability recommended as MTD")
  all_sim<-rbind.data.frame(all_sim, sim)
  print(all_sim)
  write.csv(n, paste0("N1.34.wages.marginal.sim", i,".",j, ".pat.allocation", ncohort, ".", cohortsize, ".csv"))
  write.csv(yc, paste0("N1.34.wages.marginal.sim", i, ".", j, ".", ncohort, ".", cohortsize, ".cdlt.csv"))
  write.csv(ye, paste0("N1.34.wages.marginal.sim", i, ".", j, ".", ncohort, ".", cohortsize, ".eff.csv"))
}
write.csv(all_sim, paste0("N1.34.marginal.joint.sim.", ncohort,".", cohortsize,".csv"))

############################################# Joint prior calibration #######################################################

alpha<- seq(from=0.5, to =4, length.out=100)
l_rate<- 1
l.df <- data.frame(b_shape=NULL,b_rate=NULL, l_alpha=NULL, l_beta=NULL, comb_div=NULL, clin_div=NULL, h_div=NULL)
no.d<- 6
clin <- exp(crmsens(skeletonc, target=0.33, model="empiric", detail=TRUE)$Hset)

for(i in 1:100){
  nm<-optim(par=c(1,1), KLja, l_shape=alpha[i], l_rate=l_rate, no.d=no.d, Hset_b=clin,lower=c(0,0.1), upper=c(10,10), method="L-BFGS-B")
  l.df<-rbind(l.df, data.frame(b_shape=nm$par[1],b_rate=nm$par[2],
                               l_alpha=alpha[i],l_beta=l_rate,comb_div=nm$value, clin_div=KL_marginal(nm$par[1],nm$par[2], 5, clin),
                               h_div=KL_marginal_h(alpha[i], l_rate,no.d)))
}

#choose priors which lead to dose agnostic a priori dose selection
df <- data.frame(x = l.df$l_alpha,
                 y = l.df$clin_div)
el<-find_curve_elbow(df, plot_curve = FALSE)

l_shape<- l.df[el,3]
l_rate<-l.df[el,4]
prob_l<- c(pbeta(1/6, l_shape, l_rate),pbeta(2/6, l_shape, l_rate)-pbeta(1/6, l_shape, l_rate),pbeta(3/6, l_shape, l_rate)-pbeta(2/6, l_shape, l_rate),
           pbeta(4/6, l_shape, l_rate)-pbeta(3/6, l_shape, l_rate),pbeta(5/6, l_shape, l_rate)-pbeta(4/6, l_shape, l_rate),1-pbeta(5/6, l_shape, l_rate))

prior.eff.model<- c(prob_l[6], prob_l[5]/2, prob_l[4]/2, prob_l[3]/2,prob_l[2]/2, prob_l[1]/2, prob_l[1]/2, prob_l[2]/2, prob_l[3]/2, prob_l[4]/2, prob_l[5]/2)

mat_skeletone<- matrix(c(eff1, eff2, eff3, eff4, eff5, eff6, eff7, eff8,
                         eff9,eff10, eff11), nrow=11, byrow = TRUE)
ac<- l.df[el,1]
bc<- l.df[el,2]
targetc<- 0.33
cohortsize<- 4
ncohort<- 48
rar.prop<- 0.25

true_tox<- rbind.data.frame(c(0.05, 0.1, 0.2, 0.28, 0.5, 0.5), c(0.05, 0.1, 0.2, 0.28, 0.4, 0.55), 
                            c(0.05, 0.1, 0.15, 0.2, 0.35, 0.4), c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05), rep(0.05, times=6), c(0.33, 0.45, 0.55, 0.6, 0.7, 0.75))
dimnames(true_tox)[[2]]<- c( "1","2", "3", "4", "5", "6")
dimnames(true_tox)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc 4", "Sc5", "Sc6")                       

true_eff<- rbind.data.frame(c(0.05, 0.13, 0.25, 0.38, 0.5, 0.63), c(0.05, 0.23, 0.47, 0.7, 0.7, 0.7), c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7),
                            c(0.7, 0.43, 0.38, 0.25, 0.23, 0.05), c(0.38, 0.50, 0.63, 0.5, 0.38, 0.25))
dimnames(true_eff)[[2]]<- c( "1","2", "3", "4", "5", "6")
dimnames(true_eff)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc4", "Sc5")   

sim_scenario_k<- matrix(nrow=2, ncol=16)
sim_scenario_k[1,]<- c(rep(1:4, each=3), rep(5:6, each=2))
sim_scenario_k[2,]<- c(rep(1:3, times=4), rep(4:5, times=2))

all_sim<- NULL
set.seed(1001)
ntrial=5000
for(k in 1:16){
    i<- sim_scenario_k[1,k]
    j<- sim_scenario_k[2,k]
    clin_true<- as.numeric(true_tox[i,])
    eff_true<- as.numeric(true_eff[j,])
    ndose=length(clin_true)
    dose.select<-yc<-ye<-n<-matrix(nrow=ntrial,ncol=ndose)
    cl<- makeCluster(detectCores())
    clusterSetRNGStream(cl,1001)
    invisible(clusterEvalQ(cl,{
      library(clusterGeneration)
      library(parallel)
      setwd("/annotated")
      source("efftox_functions.R")  
    }))
    clusterExport(cl, c("clin_true","eff_true", "prior.eff.model", "skeletonc","mat_skeletone","ac", "bc","targetc","cohortsize","rar.prop", "ntrial"))
    result<-parLapply(cl, 1:ntrial, function (k) efftoxcrm.uninf(clin_true,eff_true, prior.eff.model, skeletonc,mat_skeletone,ac, bc,0.1,targetc,cohortsize,48,rar.prop))
    stopCluster(cl)
    for(m in 1:ntrial){
      dose.select[m,]<-result[[m]]$dose.select
      n[m,]<- c(result[[m]]$pt.allocation)
      yc[m,]<-  c(result[[m]]$tox.datac)
      ye[m,]<-  c(result[[m]]$eff.data)
    }
    sim<-rbind.data.frame(clin_true, eff_true, colMeans(dose.select))
    dimnames(sim)[[2]]<- c( "1","2", "3", "4", "5", "6")
    dimnames(sim)[[1]]<- c("Probability of clinician DLT", "Probability of efficacy", "Probability recommended as MTD")
    all_sim<-rbind.data.frame(all_sim, sim)
    print(all_sim)
    write.csv(n, paste0(el,"HN0.1.wages.joint.sim", i,".",j, ".pat.allocation", ncohort, ".", cohortsize, ".csv"))
    write.csv(yc, paste0(el,"HN0.1.wages.joint.sim", i,".",j, ".", ncohort, ".", cohortsize, ".cdlt.csv"))
    write.csv(ye, paste0(el,"HN0.1.wages.joint.sim",i,".",j, ".", ncohort, ".", cohortsize, ".eff.csv"))
}
write.csv(all_sim, paste0(el,"HN0.1.wages.joint.sim.", ncohort,".", cohortsize,".csv"))





