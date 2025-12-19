setwd("/home/ealger/revision_calibrate_priors")
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

delta<- 0.05
target<- 0.33
mtd<- 5
# 
# skeletonc<- getprior(delta, target, mtd, nlevel=6)
# 
# eff1<- getprior(delta, 0.60, 6, nlevel=6)
# eff2<- c(eff1[2], eff1[3], eff1[4], eff1[5], eff1[6], eff1[5])
# eff3<- c(eff1[3], eff1[4], eff1[5], eff1[6], eff1[5], eff1[4])
# eff4<- c(eff1[4], eff1[5], eff1[6], eff1[5], eff1[4], eff1[3])
# eff5<- c(eff1[5], eff1[6], eff1[5], eff1[4], eff1[3], eff1[2])
# eff6<- c( eff1[6], eff1[5], eff1[4], eff1[3], eff1[2], eff1[1])
# eff7<- rep(eff1[6], times=6)
# eff8<- c(eff1[5], rep(eff1[6], times=5))
# eff9<- c(eff1[4], eff1[5], rep(eff1[6], times=4))
# eff10<- c(eff1[3], eff1[4], eff1[5], rep(eff1[6], times=3))
# eff11<- c(eff1[2], eff1[3], eff1[4], eff1[5], rep(eff1[6], times=2))

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

#uninformative values
df <- data.frame(x = l.df$l_alpha,
                 y = l.df$clin_div)
#el<-find_curve_elbow(df, plot_curve = FALSE)

el<- 20
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
targete<- 0.2
cohortsize<- 4
ncohort<- 48
rar.prop<- 0.25

true_tox<- rbind.data.frame(c(0.33, 0.45, 0.55, 0.60, 0.70, 0.75), c(0.25, 0.33, 0.45, 0.55, 0.60, 0.70), c(0.15, 0.25, 0.33, 0.45, 0.55, 0.60),
                            c(0.10, 0.15, 0.25, 0.33, 0.45, 0.55), c(0.05, 0.10, 0.15, 0.25, 0.33, 0.45), c(0.01, 0.05, 0.10, 0.15, 0.25, 0.33),
                            c(0.01, 0.05, 0.10, 0.15, 0.18, 0.25))
dimnames(true_tox)[[2]]<- c( "1","2", "3", "4", "5", "6")
dimnames(true_tox)[[1]]<- c("CLIN: Sc 1", "CLIN: Sc 2", "CLIN: Sc 3", "CLIN: Sc 4", "CLIN: Sc5", "CLIN: Sc6", "CLIN: Sc7")                       

true_eff<- rbind.data.frame(rep(0.63, times=6), c(0.63, 0.50, 0.38, 0.25, 0.13, 0.05), c(0.5, rep(0.63, times=5)), c(0.5, 0.63, 0.5, 0.38, 0.25, 0.13), 
                            c(0.38, 0.5, rep(0.63, times=4)), c(0.38, 0.50, 0.63, 0.50, 0.38, 0.25), c(0.25, 0.38, 0.5, rep(0.63, times=3)),
                            c(0.25, 0.38, 0.5, 0.63, 0.5, 0.38), c(0.13, 0.25, 0.38, 0.5, 0.63, 0.63), c(0.13, 0.25, 0.38, 0.50, 0.63, 0.50),
                            c(0.05, 0.13, 0.25, 0.38, 0.50, 0.63))
dimnames(true_eff)[[2]]<- c( "1","2", "3", "4", "5", "6")
dimnames(true_eff)[[1]]<- c("EFF: Sc 1", "EFF: Sc 2", "EFF: Sc 3", "EFF: Sc4", "EFF: Sc5", "EFF: Sc6", "EFF: Sc7", "EFF: Sc8", "EFF: Sc9", "EFF: Sc10", "EFF: Sc11")   

sim_scenario_k<- matrix(nrow=2, ncol=16)
sim_scenario_k[1,]<- c(rep(1:4, each=3), rep(5:6, each=2))
sim_scenario_k[2,]<- c(rep(1:3, times=4), rep(4:5, times=2))

all_sim<- NULL
set.seed(1001)
ntrial=5000
for (i in 1:7){
  for(j in 1:11){
    clin_true<- as.numeric(true_tox[i,])
    eff_true<- as.numeric(true_eff[j,])
    ndose=length(clin_true)
    dose.select<-yc<-ye<-n<-matrix(nrow=ntrial,ncol=ndose)
    cl<- makeCluster(detectCores())
    clusterSetRNGStream(cl,1001)
    invisible(clusterEvalQ(cl,{
      library(clusterGeneration)
      library(parallel)
      setwd("/home/ealger/revision_calibrate_priors")
      source("efftox_functions.R")  
    }))
    clusterExport(cl, c("clin_true","eff_true", "prior.eff.model", "skeletonc","mat_skeletone","ac", "bc","targetc","cohortsize","rar.prop", "ntrial"))
    result<-parLapply(cl, 1:ntrial, function (k) efftoxcrm.uninf(clin_true,eff_true, prior.eff.model, skeletonc,mat_skeletone,ac, bc,0.1,targetc,0.2, cohortsize,48,rar.prop))
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
    write.csv(n, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/joint/",el,"HN0.1.sim", i,".",j, ".pat.allocation", ncohort, ".", cohortsize, ".csv"))
    write.csv(yc, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/joint/",el,"HN0.1.sim", i,".",j, ".", ncohort, ".", cohortsize, ".cdlt.csv"))
    write.csv(ye, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/joint/",el,"HN0.1.sim",i,".",j, ".", ncohort, ".", cohortsize, ".eff.csv"))
  }
}
write.csv(all_sim, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/joint/",el,"HN0.1.sim.", ncohort,".", cohortsize,".csv"))

### MARGINAL APPROACH ####

prior.eff.model<- rep(1/11, times=11)
sc<- 1.34
se<- 1.34

all_sim<- NULL
set.seed(1001)
ntrial=5000
for (i in 1:7){
  for(j in 1:11){
  clin_true<- as.numeric(true_tox[i,])
  eff_true<- as.numeric(true_eff[j,])
  ndose=length(clin_true)
  dose.select<-yc<-ye<-n<-matrix(nrow=ntrial,ncol=ndose)
  cl<- makeCluster(detectCores())
  clusterSetRNGStream(cl,1001)
  invisible(clusterEvalQ(cl,{
    library(clusterGeneration)
    library(parallel)
    setwd("/home/ealger/revision_calibrate_priors")
    source("efftox_functions.R")  
  }))
  clusterExport(cl, c("clin_true","eff_true", "prior.eff.model", "skeletonc","mat_skeletone","sc", "se","targetc","cohortsize","rar.prop", "ntrial"))
  result<-parLapply(cl, 1:ntrial, function (k) efftoxcrm(clin_true,eff_true, prior.eff.model, skeletonc,mat_skeletone,sc, se,targetc,0.2, cohortsize,48,rar.prop))
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
  write.csv(n, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/N1.34/sim", i,".",j, ".pat.allocation", ncohort, ".", cohortsize, ".csv"))
  write.csv(yc, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/N1.34/sim", i, ".", j, ".", ncohort, ".", cohortsize, ".cdlt.csv"))
  write.csv(ye, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/N1.34/sim", i, ".", j, ".", ncohort, ".", cohortsize, ".eff.csv"))
  }
}
write.csv(all_sim, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/N1.34/sim.", ncohort,".", cohortsize,".csv"))



prior.eff.model<- rep(1/11, times=11)
sc<- 1.34
se<- 0.1

all_sim<- NULL
set.seed(1001)
ntrial=5000
for (i in 1:7){
  for(j in 1:11){
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
    setwd("/home/ealger/revision_calibrate_priors")
    source("efftox_functions.R")  
  }))
  clusterExport(cl, c("clin_true","eff_true", "prior.eff.model", "skeletonc","mat_skeletone","sc", "se","targetc","cohortsize","rar.prop", "ntrial"))
  result<-parLapply(cl, 1:ntrial, function (k) efftoxcrm.HN(clin_true,eff_true, prior.eff.model, skeletonc,mat_skeletone,sc, se,targetc,0.2, cohortsize,48,rar.prop))
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
  write.csv(n, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/HN0.1/sim", i,".",j, ".pat.allocation", ncohort, ".", cohortsize, ".csv"))
  write.csv(yc, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/HN0.1/sim", i, ".", j, ".", ncohort, ".", cohortsize, ".cdlt.csv"))
  write.csv(ye, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/HN0.1/sim", i, ".", j, ".", ncohort, ".", cohortsize, ".eff.csv"))
  }
}
write.csv(all_sim, paste0("/home/ealger/revision_calibrate_priors/results/efftox_results/marginal/HN0.1/sim.", ncohort,".", cohortsize,".csv"))




