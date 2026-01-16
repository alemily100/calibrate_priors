setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")
library(RColorBrewer)
source("efftox_functions_5doses.R")

skeletonc<-  c(0.06, 0.14, 0.25, 0.38, 0.50)
eff1<- c(0.2, 0.3, 0.4, 0.5, 0.6)
eff2<- c(0.3, 0.4, 0.5, 0.6, 0.5)
eff3<- c(0.4, 0.5, 0.6, 0.5, 0.4)
eff4<- c(0.5, 0.6, 0.5, 0.4, 0.3)
eff5<- c(0.6, 0.5, 0.4, 0.3, 0.2)
eff6<- rep(0.6, times=5)
eff7<- c(0.5, rep(0.6, times=4))
eff8<- c(0.4, 0.5, rep(0.6, times=3))
eff9<- c(0.3, 0.4, 0.5, rep(0.6, times=2))



alpha<- seq(from=0.5, to =4, length.out=100)
l_rate<- 1
l.df <- data.frame(b_shape=NULL,b_rate=NULL, l_alpha=NULL, l_beta=NULL, comb_div=NULL, clin_div=NULL, h_div=NULL)
no.d<- 5
clin <- exp(crmsens(skeletonc, target=0.25, model="empiric", detail=TRUE)$Hset)

for(i in 1:100){
  nm<-optim(par=c(1,1), KLja, l_shape=alpha[i], l_rate=l_rate, no.d=no.d, Hset_b=clin,lower=c(0,0.1), upper=c(10,10), method="L-BFGS-B")
  l.df<-rbind(l.df, data.frame(b_shape=nm$par[1],b_rate=nm$par[2],
                               l_alpha=alpha[i],l_beta=l_rate,comb_div=nm$value, clin_div=KL_marginal(c(nm$par[1],nm$par[2]), 5, clin),
                               h_div=KL_marginal_h(alpha[i], l_rate,no.d)))
}

#uninformative values
df <- data.frame(x = l.df$l_alpha,
                 y = l.df$clin_div)
#el<-find_curve_elbow(df, plot_curve = TRUE)
el<-20
df_long <- l.df %>%
  mutate(l_alpha = l_alpha) %>%
  pivot_longer(cols = c(clin_div,h_div, comb_div), 
               names_to = "type", 
               values_to = "value")

palette <- brewer.pal(6, "Set2")
df_long$type <- factor(df_long$type, levels = c("clin_div", "h_div", "comb_div"))


pdf("figures/KL_wages.pdf", width=5, height=5)
ggplot(df_long, aes(x = l_alpha, y = value, color = type)) +
  geom_line(size = 1) + 
  scale_color_manual(values = c("clin_div" = "#E78AC3", 
                                "h_div" = "#A6D854", 
                                "comb_div" = "#FFD92F"),
                     labels = c("Clinician-CRM a priori dose recommendation", "Efficacy a priori dose recommendation",
                                "Eff+Tox-CRM a priori dose recommendation")) +
  geom_vline(xintercept = alpha[el], linetype = "dashed", size = 1) + 
  labs(
    x = expression(alpha[l]), 
    y = expression(D[KL](mu*"|"*Uniform)),
    color = "Legend"
  ) + theme_classic(base_size = 15)+theme(legend.position =c(0.35,0.95))+scale_y_continuous(limits = c(-0.2, 2.5))+
  labs(color = " ")
dev.off()


set.seed(100)
N<- 1000
b_shape<-l.df[el,1] 
b_rate<- l.df[el,2]
l_shape<- l.df[el,3]
l_rate<-l.df[el,4]
prob_l<- c(pbeta(1/5, l_shape, l_rate),pbeta(2/5, l_shape, l_rate)-pbeta(1/5, l_shape, l_rate),pbeta(3/5, l_shape, l_rate)-pbeta(2/5, l_shape, l_rate),
           pbeta(4/5, l_shape, l_rate)-pbeta(3/5, l_shape, l_rate),1-pbeta(4/5, l_shape, l_rate))
prior.eff.model<- c(prob_l[5], prob_l[4]/2, prob_l[3]/2, prob_l[2]/2,prob_l[1]/2, prob_l[1]/2, prob_l[2]/2, prob_l[3]/2, prob_l[4]/2)
l<-sample(1:9, size=N, prob= prior.eff.model, replace=TRUE)
targetc<- 0.25
estimates_clin<-sapply(rgamma(N,shape=b_shape, rate =b_rate), function(k) skeletonc^k)
values<- sapply(1:N, function (k) which.min(abs(estimates_clin-0.25)[,k]))
est<-rnorm(N,0, sd=sqrt(0.1))
estimates_eff<- sapply(1:N, function(k) eval(parse(text=paste0('eff',l[k])))^exp(est[k]))

prop.n<-table(sapply(1:N, function (k) which.max(estimates_eff[1:values[k],k])))/1000
prop.clin<- table(sapply(1:N, function(j) which.min(abs(estimates_clin[,j]-targetc))))/1000
prop.eff<- table(sample(1:5, size=1000, prob_l, replace=TRUE))/1000


# creating data frame 
df <- cbind(c(prop.n, prop.clin, prop.eff),rep(1:5, times=3),c(rep("uninf",times=5), rep("clin",times=5),
                                                                rep("eff", times=5)))
colnames(df)<- c("prop", "dose","type")
df<- data.frame(df)
df$prop<- as.numeric(df$prop)
df$dose<- as.numeric(df$dose)
df$type<- as.factor(df$type)

pdf("figures/a_priori_wages.pdf", width=5, height=5)
ggplot(df, aes(dose, prop, fill = type)) +
  geom_bar(stat="identity", position = "dodge") + 
  ylim(0,1.0)+
  geom_hline(yintercept = 1/5, color = "grey", linetype = "dashed", lwd=2)+
  theme_classic(base_size = 15)+theme(legend.position =c(0.22,0.95))+
  xlab("Dose")+ylab("A priori dose recommendation")+
  labs(fill = " ")+
  scale_fill_manual(
    labels = c("Clinician-CRM", "Efficacy", "Eff+Tox-CRM"), values=c("#E78AC3", "#A6D854", "#FFD92F"))
dev.off()

#######

#original a priori recommendation 

set.seed(1001)
N<- 1000
prior.eff.model<- rep(1/9, times=9)
l<-sample(1:9, size=N, prob= prior.eff.model, replace=TRUE)
targetc<- 0.25
estimates_clin<-sapply(rgamma(N,2.78, 2.37), function(k) skeletonc^k)
est<-rnorm(N,0, sd=sqrt(1.34))
estimates_eff<- sapply(1:N, function(k) eval(parse(text=paste0('eff',l[k])))^exp(est[k]))
values<- sapply(1:N, function (k) which.min(abs(estimates_clin-0.25)[,k]))
prop.n<-table(sapply(1:N, function (k) which.max(estimates_eff[1:values[k],k])))/1000

# creating data frame 
df <- cbind(prop.n,1:5)
colnames(df)<- c("prop", "dose")
df<- data.frame(df)
df$prop<- as.numeric(df$prop)
df$dose<- as.numeric(df$dose)
#df$type<- as.factor(df$type)

pdf("figures/a_priori_wages_uncal.pdf", width=5, height=4)
ggplot(df, aes(dose, prop)) +
  geom_bar(stat="identity", position = "dodge", fill="#FFD92F") + 
  ylim(0,0.5)+
  geom_hline(yintercept = 1/6, color = "grey", linetype = "dashed", lwd=2)+
  theme_classic(base_size = 15)+theme(legend.position =c(0.22,0.95))+
  xlab("Dose")+ylab("A priori dose recommendation")+
  labs(fill = " ")+ scale_x_continuous(breaks = 1:6)
dev.off()

#######
