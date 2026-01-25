library(tidyverse)
library(cowplot)
library(nnet)
library(binom)
library(dfcrm)
library(tidyverse)
library(lattice)

source("efftox_functions_5doses.R")

setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")


KLja<- function(vec, l_shape,l_rate, no.d, Hset_b){
  prob_l<- c(pbeta(1/5, l_shape, l_rate),pbeta(2/5, l_shape, l_rate)-pbeta(1/5, l_shape, l_rate),pbeta(3/5, l_shape, l_rate)-pbeta(2/5, l_shape, l_rate),
             pbeta(4/5, l_shape, l_rate)-pbeta(3/5, l_shape, l_rate),1-pbeta(4/5, l_shape, l_rate))
  kl<- KLj(vec, no.d, Hset_b, prob_l)
  return(kl)
}

patient_outcome<- function(dose, uniform, probability_table){
  part1<-c()
  part2<-c()
  for(i in 1:4){
    part1[i]<-ifelse(i<4, uniform[i]> probability_table[i+1, dose+1]/probability_table[i, dose+1], TRUE)
    part2[i]<-ifelse(i>1, uniform[i-1] <= probability_table[i, dose+1]/probability_table[i-1, dose+1], uniform[i] <= probability_table[i, dose+1])
  }
  return(w[min(which(part1&part2 == TRUE))])
}

generate_table<- function(no.sim.table, true_tox_c, true_eff_p){
  table<- matrix(nrow=4, ncol=7)
  table[,1]<- 0:3
  for(j in 1:5){
    sim<- sapply(rep(1, times = no.sim.table), function (k) c(rbinom(k, 1, true_tox_c[j]), rbinom(k,1,true_eff_p[j])))
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

true_tox_clin<- c(0.01, 0.02, 0.05, 0.10, 0.20)

true_eff<- c(0.13, 0.25, 0.38, 0.5, 0.63)
#dimnames(true_eff)[[2]]<- c( "1","2", "3", "4", "5", "6")

complete<-function(prior.eff.model, skeletonc,skeletone,sc,se,targetc,targete, cohortsize,ncohort,rar.prop, complete_responses){
  current_dose<-c()
  estc_pos<- c()
  eff_pos<-c()
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
  while(((i-1)*cohortsize+1) <= ncohort)
  {
    # generate data for a new cohort of patients
    yc[curr] = yc[curr] + sum(sapply((sum(n)+1):(sum(n)+cohortsize), function (k) ifelse(complete_responses[curr, k]==1|complete_responses[curr, k]==3, 1, 0)));
    ye[curr] = ye[curr] + sum(sapply((sum(n)+1):(sum(n)+cohortsize), function (k) ifelse(complete_responses[curr, k]==2|complete_responses[curr, k]==3, 1, 0)));
    n[curr] = n[curr] + cohortsize;
    ##Create admissible doses 
    marginalc = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonc, y=yc,n=n, s2=sc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=-10,upper=10, skeletonc, yc, n, s2=sc,abs.tol = 0)$value/marginalc
    ptox.hatc=skeletonc**exp(estc)
    estc_pos[i]<- estc
    #Assess efficacy
    marginale<- sapply(1:nrow(skeletone), function (k) integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletone[k,], y=ye,n=n, s2=se,abs.tol = 0)$value);
    w.k<- sapply(1:nrow(skeletone), function (k) prior.eff.model[k]*marginale[k])/sum(prior.eff.model*marginale)
    k_star<- ifelse(sum(w.k==max(w.k))>1,sample(which(w.k==max(w.k)), size=1, replace=FALSE,prob=rep(1/length(which(w.k==max(w.k))), times=length(which(w.k==max(w.k))))), which(w.k==max(w.k)))
    este<- integrate(bcrmht,lower=-Inf,upper=Inf, skeletone[k_star,], ye, n, se,abs.tol = 0)$value/marginale[k_star];
    peff.hat=skeletone[k_star,]**exp(este)
    eff_pos[i]<- which.max(peff.hat)
    if(sum(n)<=rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- ifelse(length(curr_tox)>1, sample(curr_tox, size=1,prob=peff.hat[curr_tox]), max(curr_tox))
    }
    if(sum(n)>rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- which.max(peff.hat[curr_tox])
    }
    current_dose[i]<- curr
    i=i+1
  }

  return(list(current_dose, pos_est_clin=estc_pos, pos_est_eff=eff_pos))
}

uninf_complete<-function(prior.eff.model, skeletonc,skeletone,ac, bc, theta_s2, targetc,targete, cohortsize,ncohort,rar.prop,complete_responses){
  current_dose<-c()
  estc_pos<- c()
  eff_pos<-c()
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
  loc<- c(5,4,3,2,1,1,2,3,4)
  eff.m<-  which(prior.eff.model==max(prior.eff.model))
  eff.m<-ifelse(length(eff.m)>1, sample(eff.m, size=1, replace=FALSE,prob=rep(1/length(eff.m), times=length(eff.m))), eff.m)
  admiss<- which(skeletonc<=targetc)
  curr<- sample(x=1:max(admiss), size=1, prob=skeletone[1,1:max(admiss)])
  ptox.hatc = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0,ndose); # a vector of indicators for dose selection
  stopc=stope=0; #indicate if trial stops early
  i=1      
  while(((i-1)*cohortsize+1) <= ncohort)
  {
    # generate data for a new cohort of patients
    yc[curr] = yc[curr] + sum(sapply((sum(n)+1):(sum(n)+cohortsize), function (k) ifelse(complete_responses[curr, k]==1|complete_responses[curr, k]==3, 1, 0)));
    ye[curr] = ye[curr] + sum(sapply((sum(n)+1):(sum(n)+cohortsize), function (k) ifelse(complete_responses[curr, k]==2|complete_responses[curr, k]==3, 1, 0)));
    n[curr] = n[curr] + cohortsize;
    ##Model-based estimation of DLT probabilities
    marginalc = integrate(bcrmh,lower=0,upper=Inf, p=skeletonc, y=yc,n=n, a=ac, b=bc,abs.tol = 0)$value;
    estc=integrate(bcrmht,lower=0,upper=Inf, skeletonc, yc, n, a=ac, b=bc,abs.tol = 0)$value/marginalc
    ptox.hatc=skeletonc**estc
    estc_pos[i]<- estc
    marginale<- sapply(1:nrow(skeletone), function (k) integrate(bcrmh_uninf,lower=0,upper=Inf, p=skeletone[k,], y=ye,n=n, s2=theta_s2,abs.tol = 0)$value);
    w.k<- sapply(1:nrow(skeletone), function (k) prior.eff.model[k]*marginale[k])/sum(prior.eff.model*marginale)
    k_star<- ifelse(sum(w.k==max(w.k))>1,sample(which(w.k==max(w.k)), size=1, replace=FALSE,prob=rep(1/length(which(w.k==max(w.k))), times=length(which(w.k==max(w.k))))), which(w.k==max(w.k)))
    este<- integrate(bcrmht_uninf,lower=0,upper=Inf, skeletone[k_star,], ye, n, s2=theta_s2,abs.tol = 0)$value/marginale[k_star];
    peff.hat=skeletone[k_star,]**exp(este)
    eff_pos[i]<- which.max(peff.hat)
    #peff.hat=skeletone[k_star,]**este
    if(sum(n)<=rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- ifelse(length(curr_tox)>1, sample(curr_tox, size=1,prob=peff.hat[curr_tox]), max(curr_tox))
    }
    if(sum(n)>rar.prop*ncohort){
      curr_tox<- which(ptox.hatc<=targetc)
      curr<- which.max(peff.hat[curr_tox])
    }
    current_dose[i]<- curr
    i=i+1
  }

  return(list(current_dose, pos_est_clin=estc_pos, pos_est_eff=eff_pos))
}

####### case study parameters
#setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")
#source("efftox_functions.R")

skeletonc<-  c(0.02, 0.10, 0.25, 0.44, 0.62)
eff1<- c(0.2, 0.3, 0.4, 0.5, 0.6)
eff2<- c(0.3, 0.4, 0.5, 0.6, 0.5)
eff3<- c(0.4, 0.5, 0.6, 0.5, 0.4)
eff4<- c(0.5, 0.6, 0.5, 0.4, 0.3)
eff5<- c(0.6, 0.5, 0.4, 0.3, 0.2)
eff6<- rep(0.6, times=5)
eff7<- c(0.5, rep(0.6, times=4))
eff8<- c(0.4, 0.5, rep(0.6, times=3))
eff9<- c(0.3, 0.4, 0.5, rep(0.6, times=2))

delta<- 0.05
target<- 0.25
mtd<- 5
alpha<- seq(from=0.5, to =4, length.out=100)
l_rate<- 1
l.df <- data.frame(b_shape=NULL,b_rate=NULL, l_alpha=NULL, l_beta=NULL, comb_div=NULL, clin_div=NULL, h_div=NULL)
no.d<- 5
clin <- exp(crmsens(skeletonc, target=0.33, model="empiric", detail=TRUE)$Hset)
for(i in 1:100){
  nm<-optim(par=c(1,1), KLja, l_shape=alpha[i], l_rate=l_rate, no.d=no.d, Hset_b=clin,lower=c(0,0.1), upper=c(10,10), method="L-BFGS-B")
  l.df<-rbind(l.df, data.frame(b_shape=nm$par[1],b_rate=nm$par[2],
                               l_alpha=alpha[i],l_beta=l_rate,comb_div=nm$value, clin_div=KL_marginal(c(nm$par[1],nm$par[2]), 5, clin),
                               h_div=KL_marginal_h(alpha[i], l_rate,no.d)))
}
df <- data.frame(x = l.df$l_alpha,
                 y = l.df$clin_div)
#el<-find_curve_elbow(df, plot_curve = FALSE)
el<-20
l_shape<- l.df[el,3]
l_rate<-l.df[el,4]
prob_l<- c(pbeta(1/5, l_shape, l_rate),pbeta(2/5, l_shape, l_rate)-pbeta(1/6, l_shape, l_rate),pbeta(3/5, l_shape, l_rate)-pbeta(2/5, l_shape, l_rate),
           pbeta(4/5, l_shape, l_rate)-pbeta(3/5, l_shape, l_rate),1-pbeta(4/5, l_shape, l_rate))
prior.eff.model<- c(prob_l[5], prob_l[4]/2, prob_l[3]/2,prob_l[2]/2, prob_l[1]/2, prob_l[1]/2, prob_l[2]/2, prob_l[3]/2, prob_l[4]/2)
mat_skeletone<- matrix(c(eff1, eff2, eff3, eff4, eff5, eff6, eff7, eff8,
                         eff9), nrow=9, byrow = TRUE)
ac<- l.df[el,1]
bc<- l.df[el,2]
targetc<- 0.25
targete<- 0.2
cohortsize<- 3
ncohort<- 39
rar.prop<- 0.25

w<- 0:3
set.seed(17013)
u_pat<-matrix(runif(3*ncohort), nrow=ncohort)
table<-generate_table(10000, true_tox_clin,true_eff)
complete_responses<-sapply(1:ncohort, function (m) sapply(1:5, function(k) patient_outcome(k, u_pat[m,], table)))
set.seed(n)
cal_wages<-uninf_complete(prior.eff.model, skeletonc,mat_skeletone,0.99,0.15,0.1,targetc,0.2, cohortsize,39,rar.prop, complete_responses)
wages<- uninf_complete(rep(1/9, times=9), skeletonc,mat_skeletone,2.78,2.37,0.1,targetc,0.2, cohortsize,39,rar.prop, complete_responses)

wages_df<- data.frame(1:39,rep(wages[[1]], each=3),  sapply(1:39, function (k) ifelse(complete_responses[rep(wages[[1]], each=4)[k], k]==1|complete_responses[rep(wages[[1]], each=4)[k], k]==3, 1, 0)), 
                    sapply(1:39, function (k) ifelse(complete_responses[rep(wages[[1]], each=4)[k], k]==2|complete_responses[rep(wages[[1]], each=4)[k], k]==3, 1, 0)))
cal_wages_df<- data.frame(1:39,rep(cal_wages[[1]], each=3),  sapply(1:39, function (k) ifelse(complete_responses[rep(cal_wages[[1]], each=4)[k], k]==1|complete_responses[rep(cal_wages[[1]], each=4)[k], k]==3, 1, 0)), 
                    sapply(1:39, function (k) ifelse(complete_responses[rep(cal_wages[[1]], each=4)[k], k]==2|complete_responses[rep(cal_wages[[1]], each=4)[k], k]==3, 1, 0)))

wages_df$x<- seq(from =0.5, length.out=39, by=0.5)
colnames(wages_df)<- c("pat", "dose", "cdlt", "eff", "loc")

cal_wages_df$x<- seq(from =0.5, length.out=39, by=0.5)
colnames(cal_wages_df)<- c("pat", "dose", "cdlt", "eff", "loc")

both<- wages_df%>%filter(dose==cal_wages_df$dose)
wages_dlt<- wages_df%>%filter(cdlt==1 |eff==1)
cal_wages_dlt<- cal_wages_df%>%filter(cdlt==1 |eff==1)
wages_cdlt<- wages_df%>%filter(cdlt==1)
cal_wages_cdlt<-cal_wages_df%>%filter(cdlt==1)

gamma<- data.frame(c(wages[[2]], cal_wages[[2]]), c(rep("wages", times=13), rep("cal_wages", times =13)), 1:13)
colnames(gamma)<- c("val", "type", "time")
clin<-exp(crmsens(prior=skeletonc, target=0.33, model = "empiric")$Hset)

cal_wages_eff<-cal_wages[[3]]+0.2

eta<- data.frame(c(wages[[3]], cal_wages_eff), c(rep("wages", times=13), rep("cal_wages", times =13)), 1:13)
colnames(eta)<- c("val", "type", "time")

layout(matrix(c(1,2,1,3), nrow=2))
par(mar=c(1,1,1,1)+0.1)
p<-ggplot(wages_df, aes(loc, dose,  label = pat), show.legend = TRUE) +
  annotate('rect', xmin=0, xmax=20, ymin=4.5, ymax=5.5, alpha=.2, fill='#66C2A5')+
  geom_point(size = 10, shape = 21, fill = "white", colour="#B2182B", stroke=2) +geom_point(data=cal_wages_df, aes(loc, dose), size = 10, shape = 21, fill = "white", colour="#2166AC", stroke=2)+geom_point(data=both, aes(loc, dose), size = 10, shape = 21, fill = "white", colour="black", stroke=2)+
  geom_point(data=wages_dlt, aes(loc, dose), size = 10, shape = 21, fill = "#A6D854", stroke=2)+
  geom_point(data=cal_wages_dlt, aes(loc, dose), size = 10, shape = 21, fill = "#A6D854", stroke=2)+
  geom_point(data=wages_cdlt, aes(loc, dose), size = 10, shape = 21, fill = "yellow", stroke=2)+
  geom_point(data=cal_wages_cdlt, aes(loc, dose), size = 10, shape = 21, fill = "yellow", stroke=2)+
  geom_point(size = 10, shape = 21, colour="#B2182B", stroke=2) +geom_point(data=cal_wages_df, aes(loc, dose), size = 10, shape = 21,  colour="#2166AC", stroke=2)+geom_point(data=both, aes(loc, dose), size = 10, shape = 21,  colour="black", stroke=2)+
  geom_text(data=wages_df, aes(label = pat), vjust = 0.4, size=5.5)+geom_text(data=cal_wages_df, aes(label = pat), vjust = 0.4, size=5.5)+geom_text(data=both, aes(label = pat), vjust = 0.4, size=5.5)+
  xlab("Cohort")+ ylab("Dose")+scale_x_continuous(breaks=seq(from=1, by=1.5, length.out=13), labels=1:13)+theme_classic(base_size = 25)+ theme(panel.grid.major = element_blank(),
                                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                            panel.border = element_blank(), 
                                                                                                                            panel.background = element_blank(), plot.caption = element_text(size=15, hjust=0))+
  scale_y_continuous(limits = c(1, 5.75), breaks=1:5, labels=1:5)+geom_vline(xintercept = 6.25, linetype = "dotted", color = "grey", size = 1)+
  annotate("text", x = 3, y = 5.7, label = "Randomisation stage", size = 7, colour="gray30")+
  annotate("text", x = 10, y = 5.7, label = "Maximisation stage", size = 7, colour="gray30")
  

p
l_gamma<- ggplot(gamma, aes(x=time, y=val, colour=type))+geom_line(size=1.5)+ scale_colour_manual(values=c("#2166AC","#B2182B"))+
  theme_classic(base_size = 25)+theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ labs(x = "Cohort", y =expression(hat(lambda)))+annotate('rect', xmin=0, xmax=13, ymin=clin[5,1], ymax=12, alpha=.2, fill='#66C2A5')+
  scale_y_continuous(limits=c(-0.5,12))+scale_x_continuous(breaks=1:13, labels=1:13)


l_eta<- ggplot(eta, aes(x=time, y=val, colour=type))+geom_point(stroke=2,shape = 4,size=3,)+ scale_colour_manual(values=c("#2166AC", "#B2182B"))+
  theme_classic(base_size = 25)+theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
                    plot.caption = element_text(size=15, hjust=0))+ labs(x = "Cohort", y ="Efficacious dose \n recommendation")+annotate('rect', xmin=0, xmax=13.5, ymin=4.5, ymax=5.5, alpha=.2, fill='#66C2A5')+
  scale_y_continuous(limits=c(0,5.5), breaks=1:5, labels=1:5)+scale_x_continuous(breaks=1:13, labels=1:13)
p

par(mar=c(0,0,0,0))
pdf("figures/wages_legend.pdf", width=13)
plot(NULL, xlim=c(0,10), xaxt="n", ylim=c(0,10), yaxt="n", bty="n", ylab="",
     xlab="")
legend("topleft",legend=c(expression("Dose allocation for jointly calibrated priors"),expression("Dose allocation for marginally calibrated priors"),
                          "Overlapping dose allocation ", "DLT observation", "Efficacy observation"), pch=21, cex=3, pt.cex=5, pt.lwd=2.5,col=c("#2166AC", "#B2182B", "Black", "White", "White"), pt.bg=c("white","white","white","yellow", "#A6D854"), box.col="white")
dev.off()

pdf("figures/wages_prior_info.pdf", width=20, height=6)

par(mar=c(1,1,1,1)+0.1)
p
dev.off()

pdf("figures/wages_gamma_est.pdf", height=4, width=5.5)
l_gamma
dev.off()

pdf("figures/wages_eta_est.pdf", height=4, width=5.5)
l_eta
dev.off()

