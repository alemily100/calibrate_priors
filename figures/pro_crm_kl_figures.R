library(RColorBrewer)
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")
source("procrm_functions.R")
###Specify a set of skeleton values
n39.c<- c(0.06, 0.14, 0.25, 0.38, 0.50)
n39.p<- c(0.10, 0.21, 0.35, 0.49, 0.61)
ncohort=39      ##number of cohorts

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

df_long <- l %>%
  mutate(rate = rate) %>%
  pivot_longer(cols = c(clin_div, pat_div, comb_div))

df_long$name <- factor(df_long$name, levels = c("clin_div", "pat_div", "comb_div"))

pdf("figures/KL_PRO_CRM_B.pdf", width=5, height=5)
ggplot(df_long, aes(x = rate, y = value, color = name)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("clin_div" = "#66C2A5", 
                                "pat_div" = "#FC8D62", 
                                "comb_div" = "#8DA0CB"),
                     labels=c("Clinician-CRM a priori dose recommendation", "Patient-CRM a priori dose recommendation",
                              "PRO-CRM-B a priori dose recommendation")) +
  geom_vline(xintercept = rate[56], linetype = "dashed", size = 1) + 
  labs(
    x = expression(beta[xi]), 
    y = expression(D[KL](Uniform*"|"*mu[theta])),
    color = "Legend"
  ) + theme_classic(base_size = 15)+theme(legend.position =c(0.35,0.95))+scale_y_continuous(limits = c(-0.2, 2.5))+
  labs(color = " ")
dev.off()

best<-l[which.min(abs(l[,6]-l[,7])),]

eval(parse(text=paste0("ac.", ncohort, "<-as.numeric(best[1])")))
eval(parse(text=paste0("ap.", ncohort, "<-as.numeric(best[2])")))
eval(parse(text=paste0("bc.", ncohort, "<-as.numeric(best[3])")))
eval(parse(text=paste0("bp.", ncohort, "<-as.numeric(best[4])")))

set.seed(100)
estimates_clin_39<-sapply(rgamma(1000, ac.39, bc.39), function(k) n39.c^k)
estimates_pat_39<-sapply(rgamma(1000, ap.39, bp.39), function(k)  n39.p^k)
prop.n39<-table(sapply(1:1000, function(j) min(which.min(abs(estimates_clin_39[,j]-0.25)), which.min(abs(estimates_pat_39[,j]-0.35)))))/1000
prop.clin.n39<- table(sapply(1:1000, function(j) which.min(abs(estimates_clin_39[,j]-0.25))))/1000
prop.pat.n39<- table(sapply(1:1000, function(j) which.min(abs(estimates_pat_39[,j]-0.35))))/1000


# creating data frame 
df <- cbind(c(prop.n39, prop.clin.n39, prop.pat.n39),rep(1:5, times=3),c(rep("uninf",times=5), rep("clin",times=5),
                                                                rep("pat", times=5)))
colnames(df)<- c("prop", "dose","type")
df<- data.frame(df)
df$prop<- as.numeric(df$prop)
df$dose<- as.numeric(df$dose)
df$type<- as.factor(df$type)

pdf("figures/a_priori_procrmb.pdf", width=5, height=5)
ggplot(df, aes(dose, prop, fill = type)) +
  geom_bar(stat="identity", position = "dodge") + 
  ylim(0,0.5)+
  geom_hline(yintercept = 0.2, color = "grey", linetype = "dashed", lwd=2)+
  theme_classic(base_size = 15)+theme(legend.position =c(0.22,0.95))+
  xlab("Dose")+ylab("A priori dose recommendation")+
  labs(fill = " ")+
  scale_fill_manual(
    labels = c("Clinician-CRM", "Patient-CRM", "PRO-CRM-B"), values=c("#66C2A5", "#FC8D62", "#8DA0CB"))
dev.off()

#######


#original 

set.seed(100)
cdlt_param<-optim(par=c(1,1), KL_marg, no.d=5, Hset=clin, lower=c(0, 0), upper=c(10,10), method="L-BFGS-B")$par
pdlt_param<-optim(par=c(1,1), KL_marg, no.d=5, Hset=pat, lower=c(0, 0), upper=c(10,10), method="L-BFGS-B")$par

estimates_clin_39<-sapply(rgamma(10000, cdlt_param[1], cdlt_param[2]), function(k) n39.c^k)
estimates_pat_39<-sapply(rgamma(10000, pdlt_param[1], pdlt_param[2]), function(k)  n39.p^k)
prop.n39<-table(sapply(1:10000, function(j) min(which.min(abs(estimates_clin_39[,j]-0.25)), which.min(abs(estimates_pat_39[,j]-0.35)))))/10000



# creating data frame 
df <- cbind(prop.n39, 1:5)
colnames(df)<- c("prop", "dose")
df<- data.frame(df)
df$prop<- as.numeric(df$prop)
df$dose<- as.numeric(df$dose)

pdf("figures/a_priori_procrmb_uncal.pdf", width=5, height=4)
ggplot(df, aes(dose, prop)) +
  geom_bar(stat="identity", position = "dodge", fill="#8DA0CB") + 
  ylim(0,0.5)+
  geom_hline(yintercept = 0.2, color = "grey", linetype = "dashed", lwd=2)+
  theme_classic(base_size = 15)+theme(legend.position =c(0.22,0.95))+
  xlab("Dose")+ylab("A priori dose recommendation")+
  labs(fill = " ")
dev.off()
