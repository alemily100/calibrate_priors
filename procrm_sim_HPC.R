setwd("/home/ealger/revision_calibrate_priors")
source("procrm_functions.R")

start=1                                     ##starting dose
targetc=0.25      ##target c toxicity rate
targetp=0.35      ##target p toxicity rate

###Specify a set of skeleton values
n18.c<- c(0.02, 0.10, 0.25, 0.44, 0.62)
n18.p<- c(0.06, 0.18, 0.35, 0.53, 0.68)
n39.c<- c(0.06, 0.14, 0.25, 0.38, 0.50)
n39.p<- c(0.10, 0.21, 0.35, 0.49, 0.61)

sc.18 <- 0.938**2
sp.18 <- 0.893**2
sc.39<- 0.632**2
sp.39<- 0.696**2

cohortsize=3     ##cohort size for each inclusion
ncohort=39        ##number of cohorts
n.stop=ncohort          ##Number of patients needed on one combination to stop the trial
cl=0.9999         ##confidence level for the confidence interval
phi<- 0.9

skeletonc<- eval(parse(text=paste0("n", ncohort, ".c")))
skeletonp<- eval(parse(text=paste0("n", ncohort, ".p")))
sc<- eval(parse(text=paste0("sc.", ncohort)))
sp<- eval(parse(text=paste0("sp.", ncohort)))

true_tox_clin<- rbind.data.frame(c(0.05, 0.05, 0.25, 0.4, 0.55), c(0.05, 0.25, 0.40, 0.55, 0.7), c(0.01, 0.02, 0.05, 0.10, 0.25),
                                 c(0.02, 0.05, 0.1, 0.25, 0.4), c(0.05, 0.1, 0.16, 0.25, 0.4), c(0.05, 0.18, 0.2, 0.25, 0.4),
                                 c(0.01, 0.05, 0.1, 0.16, 0.25),c(0.25, 0.40, 0.55, 0.7, 0.8), c(0.45, 0.50, 0.55, 0.7, 0.8))
dimnames(true_tox_clin)[[2]]<- c( "1","2", "3", "4", "5")
dimnames(true_tox_clin)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc 4", "Sc 5", "Sc 6", "Sc 7", "Sc 8", "Sc 9")                       

true_tox_pro<- rbind.data.frame(c(0.17, 0.18, 0.35, 0.50, 0.65), c(0.1, 0.15, 0.35, 0.5, 0.65), c(0.04, 0.09, 0.17, 0.2, 0.35),
                                c(0.09, 0.17, 0.2, 0.35, 0.5), c(0.05, 0.2, 0.35, 0.5, 0.65), c(0.17, 0.35, 0.5, 0.65, 0.8),
                                c(0.04, 0.05, 0.2, 0.35, 0.5),c(0.35, 0.5, 0.65, 0.8, 0.85), c(0.55, 0.6, 0.65, 0.8, 0.85))
dimnames(true_tox_pro)[[2]]<- c( "1","2", "3", "4", "5")
dimnames(true_tox_pro)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc 4", "Sc 5", "Sc 6", "Sc 7","Sc 8", "Sc 9")   

all_sim<- NULL
set.seed(1001)
ntrial=5000
for (i in 1:9){
  clin_true_tox<- as.numeric(true_tox_clin[i,])
  pat_true_tox<- as.numeric(true_tox_pro[i,])
  prob_dlt<-procrm.sim(clin_true_tox,pat_true_tox,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl,ntrial, phi)
  sim<-rbind.data.frame(clin_true_tox, pat_true_tox, prob_dlt[[1]])
  dimnames(sim)[[2]]<- c( "1","2", "3", "4", "5")
  dimnames(sim)[[1]]<- c("Probability of clinician DLT", "Probability of patient DLT", "Probability recommended as MTD")
  all_sim<-rbind.data.frame(all_sim, sim)
#  print(all_sim)
  write.csv(all_sim, paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/marginal/clay", phi,".", ncohort,".", cohortsize,".csv"))
  write.csv(prob_dlt[[4]], paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/marginal/clay", phi, ".sim", i, ".", ncohort, ".", cohortsize, ".csv"))
  write.csv(prob_dlt[[2]], paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/marginal/clay", phi, ".sim", i, ".", ncohort, ".", cohortsize, ".cdlt.csv"))
  write.csv(prob_dlt[[3]], paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/marginal/clay", phi, ".sim", i, ".", ncohort, ".", cohortsize, ".pdlt.csv"))
}

############################################# New uninformative prior - add correlation between C-DLT and P-DLT #######################################################

skeletonc<- eval(parse(text=paste0("n", ncohort, ".c")))
skeletonp<- eval(parse(text=paste0("n", ncohort, ".p")))
clin<-exp(crmsens(prior=skeletonc, target=0.25, model = "empiric")$Hset)
pat<- exp(crmsens(prior=skeletonp, target=0.35, model = "empiric")$Hset)


#choose priors which lead to uninformative a priori dose selection  
rate<- seq(from=0, to =2, length.out=100)
l <- data.frame(ac=NULL,ap=NULL, bc=NULL, bp=NULL, comb_div=NULL, clin_div=NULL, pat_div=NULL)
set.seed(1001)
for(i in 1:100){
  nm<-optim(par=c(1,1,1,1), KL, no.d=5, Hset_b=clin, Hset_g=pat,lower=c(0,0.1,0,rate[i]), upper=c(10,10,10,rate[i]+0.01), method="L-BFGS-B")
  l<-rbind(l, data.frame(ac=nm$par[1],ap=nm$par[2],
                         bc=nm$par[3],bp=nm$par[4],comb_div=nm$value, clin_div=KL_marginal(nm$par[1],nm$par[3], 5, clin), pat_div=KL_marginal(nm$par[2],nm$par[4], 5, pat)))
}

eval(parse(text=paste0("write.csv(l, '/home/ealger/revision_calibrate_priors/results/procrm_results/joint/l.", ncohort,".csv')")))

best<-l[which.min(abs(l[,6]-l[,7])),]

eval(parse(text=paste0("ac.", ncohort, "<-as.numeric(best[1])")))
eval(parse(text=paste0("ap.", ncohort, "<-as.numeric(best[2])")))
eval(parse(text=paste0("bc.", ncohort, "<-as.numeric(best[3])")))
eval(parse(text=paste0("bp.", ncohort, "<-as.numeric(best[4])")))

eval(parse(text=paste0("write.csv(c(ac.",ncohort,", bc.", ncohort, ", ap.", ncohort, ", bp.", ncohort,"), '/home/ealger/revision_calibrate_priors/results/procrm_results/joint/best.sr.", ncohort,".csv')")))

all_sim<- NULL
estc<-NULL
estp<-NULL
set.seed(1001)
ntrial<-5000
cohortsize=3      ##cohort size for each inclusion
n.stop=ncohort
phi<- 0.9

skeletonc<- eval(parse(text=paste0("n", ncohort, ".c")))
skeletonp<- eval(parse(text=paste0("n", ncohort, ".p")))
eval(parse(text=paste0("ac<-ac.",ncohort)))
eval(parse(text=paste0("bc<-bc.",ncohort)))
eval(parse(text=paste0("ap<-ap.",ncohort)))
eval(parse(text=paste0("bp<-bp.",ncohort)))

for (i in 1:8){
  clin_true_tox<- as.numeric(true_tox_clin[i,])
  pat_true_tox<- as.numeric(true_tox_pro[i,])
  prob_dlt<-procrm.sim_gamma(clin_true_tox,pat_true_tox,skeletonc ,skeletonp,ac,bc, ap, bp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl,ntrial, phi)
  sim<-rbind.data.frame(clin_true_tox, pat_true_tox, prob_dlt[[1]])
  dimnames(sim)[[2]]<- c( "1","2", "3", "4", "5")
  dimnames(sim)[[1]]<- c("Probability of clinician DLT", "Probability of patient DLT", "Probability recommended as MTD")
  all_sim<-rbind.data.frame(all_sim, sim)
#  print(all_sim)
  write.csv(all_sim, paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/joint/clay",phi,".", ncohort, ".", cohortsize, ".csv"))
  write.csv(prob_dlt[[4]], paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/joint/clay", phi, ".sim", i, ".", ncohort, ".", cohortsize, ".csv"))
  write.csv(prob_dlt[[2]], paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/joint/clay", phi, ".sim", i, ".", ncohort, ".", cohortsize, ".cdlt.csv"))
  write.csv(prob_dlt[[3]], paste0("/home/ealger/revision_calibrate_priors/results/procrm_results/joint/clay", phi, ".sim", i, ".", ncohort, ".", cohortsize, ".pdlt.csv"))
}

