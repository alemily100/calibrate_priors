library(tidyverse)
library(cowplot)
library(nnet)
library(binom)
library(dfcrm)
library(tidyverse)
library(lattice)

setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")

patient_outcome<- function(dose, uniform, probability_table){
  part1<-c()
  part2<-c()
  for(i in 1:4){
    part1[i]<-ifelse(i<4, uniform[i]> probability_table[i+1, dose+1]/probability_table[i, dose+1], TRUE)
    part2[i]<-ifelse(i>1, uniform[i-1] <= probability_table[i, dose+1]/probability_table[i-1, dose+1], uniform[i] <= probability_table[i, dose+1])
  }
  return(w[min(which(part1&part2 == TRUE))])
}

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

# generate_table<- function(no.sim.table, true_tox_c, true_tox_p){
#   table<- matrix(nrow=4, ncol=6)
#   table[,1]<- 0:3
#   for(j in 1:5){
#     sim<- sapply(rep(1, times = no.sim.table), function (k) c(rbinom(1, k, true_tox_c[j]), rbinom(1,k,true_tox_p[j])))
#     wl_values<-sapply(1:no.sim.table, function(k) wl(sim[,k]))
#     for(i in 0:3){
#       table[i+1, j+1]<- sum(wl_values>=(i))/no.sim.table
#     }
#   }
#   return(table)
# }

generate_table<- function(no.sim.table, phi, true_tox_c, true_tox_p){
  table<- matrix(nrow=4, ncol=6)
  table[,1]<- 0:3
  for(j in 1:5){
    sim<- sapply(rep(phi, times = no.sim.table), function (k) time_to_dlt(j, true_tox_c, true_tox_p, k))<=1
    wl_values<-sapply(1:no.sim.table, function(k) wl(sim[,k]))
    for(i in 0:3){
      table[i+1, j+1]<- sum(wl_values>=(i))/no.sim.table
    }
  }
  return(table)
}

  
wl<- function(sim_col){
  test<-c()
  if(sum(sim_col==c(0,0))==2){
    test<-0
  }
  if(sum(sim_col==c(1,0))==2){
    test<-1
  }
  if(sum(sim_col==c(0,1))==2){
    test<-2
  }
  if(sum(sim_col==c(1,1))==2){
    test<-3
  }
  return(test)
}

true_tox_clin<- rbind.data.frame(c(0.05, 0.05, 0.25, 0.4, 0.55), c(0.05, 0.25, 0.40, 0.55, 0.7), c(0.02, 0.05, 0.10, 0.20, 0.25),
                                 c(0.02, 0.05, 0.1, 0.25, 0.4), c(0.05, 0.1, 0.16, 0.25, 0.4), c(0.05, 0.18, 0.2, 0.25, 0.4), 
                                 c(0.01, 0.05, 0.1, 0.16, 0.25),c(0.25, 0.40, 0.55, 0.7, 0.8), c(0.45, 0.50, 0.55, 0.7, 0.8))
dimnames(true_tox_clin)[[2]]<- c( "1","2", "3", "4", "5")
dimnames(true_tox_clin)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc 4", "Sc 5", "Sc 6", "Sc 7", "Sc 8", "Sc 9")                        

true_tox_pro<- rbind.data.frame(c(0.17, 0.18, 0.35, 0.50, 0.65), c(0.1, 0.15, 0.35, 0.5, 0.65), c(0.09, 0.17, 0.20, 0.25, 0.35),
                                c(0.09, 0.17, 0.2, 0.35, 0.5), c(0.05, 0.2, 0.35, 0.5, 0.65), c(0.17, 0.35, 0.5, 0.65, 0.8),
                                c(0.04, 0.05, 0.2, 0.35, 0.5),c(0.35, 0.5, 0.65, 0.8, 0.85), c(0.55, 0.6, 0.65, 0.8, 0.85))
dimnames(true_tox_pro)[[2]]<- c( "1","2", "3", "4", "5")
dimnames(true_tox_pro)[[1]]<- c("Sc 1", "Sc 2", "Sc 3", "Sc 4", "Sc 5", "Sc 6", "Sc 7","Sc 8", "Sc 9") 

procrm_complete<-function(skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl, complete_responses){
  estc_pos<- c()
  estp_pos<- c()
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
  current_dose<- c()
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
    current_dose[i]<- curr
    yc[curr] = yc[curr] + sum(sapply((sum(n)+1):(sum(n)+3), function (k) ifelse(complete_responses[curr, k]==1|complete_responses[curr, k]==3, 1, 0)));
    yp[curr] = yp[curr] + sum(sapply((sum(n)+1):(sum(n)+3), function (k) ifelse(complete_responses[curr, k]==2|complete_responses[curr, k]==3, 1, 0)));
    n[curr] = n[curr] + cohortsize;
    
    if(any(n>n.stop)){
      stop<-0
      break
    }
    ##Model-based estimation of DLT probabilities
    marginalc = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonc, y=yc,n=n, s2=sc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=-10,upper=10, skeletonc, yc, n, sc,abs.tol = 0)$value/marginalc
    
    marginalp = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonp, y=yp,n=n,s2=sp,abs.tol = 0)$value;
    estp=integrate(bcrmht,lower=-10,upper=10, skeletonp, yp, n, sp, abs.tol = 0)$value/marginalp
    
    ptox.hatc=skeletonc**exp(estc)
    ptox.hatp=skeletonp**exp(estp)
    estc_pos[i]<- estc
    estp_pos[i]<- estp
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
  return(list(current_dose, pos_est_clin=estc_pos, pos_est_pat=estp_pos))
}


procrm_gamma_complete<-function(skeletonc,skeletonp,ac, bc, ap, bp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl, alpha, complete_responses){
  estc_pos<- c()
  estp_pos<- c()
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
  current_dose<- c()
  ###simulate a trial
  ndose = length(skeletonc);   #number of combos
  yc=yp=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
  curr = start;  # current dose level
  current_dose[1]<- curr
  ptox.hatc = ptox.hatp = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  stopc=stopp=0; #indicate if trial stops early
  i=1  
  while(((i-1)*cohortsize+1) <= ncohort)
  {
    # generate data for a new cohort of patients
    current_dose[i]<- curr
    yc[curr] = yc[curr] + sum(sapply((sum(n)+1):(sum(n)+3), function (k) ifelse(complete_responses[curr, k]==1|complete_responses[curr, k]==3, 1, 0)));
    yp[curr] = yp[curr] + sum(sapply((sum(n)+1):(sum(n)+3), function (k) ifelse(complete_responses[curr, k]==2|complete_responses[curr, k]==3, 1, 0)));
    n[curr] = n[curr] + cohortsize;
    
    if(any(n>n.stop)){
      stop<-0
      break
    }
    ##Model-based estimation of DLT probabilities
    marginalc = integrate(bcrmh,lower=0,upper=Inf, p=skeletonc, y=yc,n=n, a=ac, b=bc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=0,upper=10, p=skeletonc, y=yc, n=n, a=ac, b=bc,abs.tol = 0)$value/marginalc
    
    marginalp = integrate(bcrmh,lower=0,upper=Inf, p=skeletonp, y=yp,n=n,a=ap, b=bp,abs.tol = 0)$value;
    estp=integrate(bcrmht,lower=0,upper=10, p=skeletonp, y=yp, n=n, a=ap, b=bp, abs.tol = 0)$value/marginalp
    
    ptox.hatc=skeletonc**estc
    ptox.hatp=skeletonp**estp
    estc_pos[i]<- estc
    estp_pos[i]<- estp
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
  return(list(current_dose, pos_est_clin=estc_pos, pos_est_pat=estp_pos))
}

n39.c<- c(0.06, 0.14, 0.25, 0.38, 0.50)
n39.p<- c(0.10, 0.21, 0.35, 0.49, 0.61)
skeletonc<-n39.c
skeletonp<-n39.p
sc.39<- 0.632**2
sp.39<- 0.696**2
cl=0.9999  

w<- 0:3
set.seed(17184)
u_pat<-matrix(runif(3*39), nrow=39)
table<-generate_table(10000, 0.9, as.numeric(true_tox_clin[3,]),as.numeric(true_tox_pro[3,]))
complete_responses<-sapply(1:39, function (m) sapply(1:5, function(k) patient_outcome(k, u_pat[m,], table)))
pro<-procrm_gamma_complete(n39.c, n39.p, 2.78, 2.37, 2.34, 1.96, 0.25, 0.35, 3, 39, 39, 1, cl, 2, complete_responses)
cal_pro<-procrm_gamma_complete(n39.c, n39.p, 2.42, 1.35, 2.09, 1.11, 0.25, 0.35, 3, 39, 39, 1, cl, 2,complete_responses)

procrm<- data.frame(1:39,rep(pro[[1]], each=3),  sapply(1:39, function (k) ifelse(complete_responses[rep(pro[[1]], each=3)[k], k]==1|complete_responses[rep(pro[[1]], each=3)[k], k]==3, 1, 0)), 
                    sapply(1:39, function (k) ifelse(complete_responses[rep(pro[[1]], each=3)[k], k]==2|complete_responses[rep(pro[[1]], each=3)[k], k]==3, 1, 0)))
cal_procrm<- data.frame(1:39,rep(cal_pro[[1]], each=3),  sapply(1:39, function (k) ifelse(complete_responses[rep(cal_pro[[1]], each=3)[k], k]==1|complete_responses[rep(cal_pro[[1]], each=3)[k], k]==3, 1, 0)), 
                    sapply(1:39, function (k) ifelse(complete_responses[rep(cal_pro[[1]], each=3)[k], k]==2|complete_responses[rep(cal_pro[[1]], each=3)[k], k]==3, 1, 0)))




procrm$x<- seq(from =0.5, length.out=39, by=0.5)
colnames(procrm)<- c("pat", "dose", "cdlt", "pdlt", "loc")

cal_procrm$x<- seq(from =0.5, length.out=39, by=0.5)
colnames(cal_procrm)<- c("pat", "dose", "cdlt", "pdlt", "loc")

both<- procrm%>%filter(dose==cal_procrm$dose)
procrm_dlt<- procrm%>%filter(cdlt==1 |pdlt==1)
cal_procrm_dlt<- cal_procrm%>%filter(cdlt==1 |pdlt==1)
procrm_pdlt<- procrm%>%filter(pdlt==1)
cal_procrm_pdlt<- cal_procrm%>%filter(pdlt==1)

gamma<- data.frame(c(pro[[2]], cal_pro[[2]]), c(rep("pro", times=13), rep("cal_pro", times =13)), 1:13)
colnames(gamma)<- c("val", "type", "time")
clin<-exp(crmsens(prior=skeletonc, target=0.25, model = "empiric")$Hset)
pat<- exp(crmsens(prior=skeletonp, target=0.35, model = "empiric")$Hset)

eta<- data.frame(c(pro[[3]], cal_pro[[3]]), c(rep("pro", times=13), rep("cal_pro", times =13)), 1:13)
colnames(eta)<- c("val", "type", "time")

layout(matrix(c(1,2,1,3), nrow=2))
par(mar=c(6,6,1,10)+0.1)
p<-ggplot(procrm, aes(loc, dose,  label = pat), show.legend = TRUE) +
  annotate('rect', xmin=0, xmax=19.5, ymin=4.5, ymax=5.5, alpha=.2, fill='#66C2A5')+
  geom_point(size = 10, shape = 21, fill = "white", colour="#B2182B", stroke=2) +geom_point(data=cal_procrm, aes(loc, dose), size = 10, shape = 21, fill = "white", colour="#2166AC", stroke=2)+geom_point(data=both, aes(loc, dose), size = 10, shape = 21, fill = "white", colour="black", stroke=2)+
  geom_point(data=procrm_dlt, aes(loc, dose), size = 10, shape = 21, fill = "yellow", stroke=2)+
  geom_point(data=cal_procrm_dlt, aes(loc, dose), size = 10, shape = 21, fill = "yellow", stroke=2)+
  geom_point(data=procrm_pdlt, aes(loc, dose), size = 10, shape = 21, fill = "#FC8D62", stroke=2)+
  geom_point(data=cal_procrm_pdlt, aes(loc, dose), size = 10, shape = 21, fill = "#FC8D62", stroke=2)+
  geom_point(size = 10, shape = 21, colour="#B2182B", stroke=2) +geom_point(data=cal_procrm, aes(loc, dose), size = 10, shape = 21,  colour="#2166AC", stroke=2)+geom_point(data=both, aes(loc, dose), size = 10, shape = 21,  colour="black", stroke=2)+
  geom_text(data=procrm, aes(label = pat), vjust = 0.4)+geom_text(data=cal_procrm, aes(label = pat), vjust = 0.4)+geom_text(data=both, aes(label = pat), vjust = 0.4)+
  xlab("Cohort")+ ylab("Dose")+scale_x_continuous(breaks=seq(from=1, by=1.5, length.out=13), labels=1:13)+theme_classic(base_size = 25)+ theme(panel.grid.major = element_blank(),
                                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                            panel.border = element_blank(), 
                                                                                                                            panel.background = element_blank(), plot.caption = element_text(size=15, hjust=0))
  
l_gamma<- ggplot(gamma, aes(x=time, y=val, colour=type))+geom_line(size=1.5)+ scale_colour_manual(values=c("#2166AC","#B2182B"))+
  theme_classic(base_size = 25)+theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ labs(x = "Cohort", y =expression(hat(lambda)))+annotate('rect', xmin=0, xmax=13, ymin=clin[5,1], ymax=3.5, alpha=.2, fill='#66C2A5')+
  scale_y_continuous(limits=c(0,3.5))+scale_x_continuous(breaks=1:13, labels=1:13)


l_eta<- ggplot(eta, aes(x=time, y=val, colour=type))+geom_line(size=1.5)+ scale_colour_manual(values=c("#2166AC", "#B2182B"))+
  theme_classic(base_size = 25)+theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
                    plot.caption = element_text(size=15, hjust=0))+ labs(x = "Cohort", y =expression(hat(xi)))+annotate('rect', xmin=0, xmax=13, ymin=pat[5,1], ymax=3.5, alpha=.2, fill='#66C2A5')+
  scale_y_continuous(limits=c(0.5,3.5))+scale_x_continuous(breaks=1:13, labels=1:13)

par(mar=c(0,0,0,0))
pdf("figures/legend.pdf", width=13)
plot(NULL, xlim=c(0,10), xaxt="n", ylim=c(0,10), yaxt="n", bty="n", ylab="", 
     xlab="")
legend("topleft",legend=c(expression("Dose allocation for jointly calibrated priors"),expression("Dose allocation for marginally uncalibrated priors"),
                          "Overlapping dose allocation ", "C-DLT observation", "P-DLT observation"), pch=21, cex=3, pt.cex=5, pt.lwd=2.5,col=c("#2166AC", "#B2182B", "Black", "White", "White"), pt.bg=c("white","white","white","yellow","#FC8D62"), box.col="white")
dev.off()


pdf("figures/prior_info.pdf", width=17, height=6)
par(mar=c(6,6,1,2)+0.1)
p
dev.off()

pdf("figures/gamma_est.pdf", height=4, width=5.5)
l_gamma
dev.off()

pdf("figures/eta_est.pdf", height=4, width=5.5)
l_eta
dev.off()

p
