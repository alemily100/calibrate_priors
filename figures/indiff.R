library(tidyverse)
setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")

joint<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/joint/0.07sim.39.3.csv")[,-1]
marginal<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/marginal_gamma/0.07sim.39.3.csv")[,-1]

dlt.rec<-sapply(((0:53)*3)+1, function (k) which.min(abs(joint[k,]-0.25)))
collect.eff<-sapply(1:54, function (k) joint[((k-1)*3)+2,1:dlt.rec[[k]]])
optimal.rec<-sapply(1:54, function (k) which(collect.eff[[k]] == max(collect.eff[[k]])))

joint.rec<-sapply(1:54, function(k) sum(joint[((k-1)*3)+3, optimal.rec[[k]]]))
marginal.rec<-sapply(1:54, function(k) sum(marginal[((k-1)*3)+3, optimal.rec[[k]]]))

mean(marginal.rec)
mean(joint.rec)

list.length<-sapply(optimal.rec, length)
index<-unlist(optimal.rec[list.length==1])
sd(tapply(marginal.rec[list.length==1], index, mean))

sd(tapply(joint.rec[list.length==1], index, mean))


###looking at get.prior results
joint0.03<- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/joint/0.03sim.39.3.csv")[,-1]
marginal0.03<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/marginal_gamma/0.03sim.39.3.csv")[,-1]

set.seed(1)
scenarios<-sample(x=54, size=12, replace=FALSE)
rows<-sapply(1:12, function(k) joint[(3*(k-1)+1):3*k,])
random.joint0.03<-sapply(sort(scenarios), function (k) sum(joint0.03[3*k,optimal.rec[[k]]]))
random.marginal0.03<-sapply(sort(scenarios), function (k) sum(marginal0.03[3*k,optimal.rec[[k]]]))

joint0.05<- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/joint/0.05sim.39.3.csv")[,-1]
marginal0.05<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/marginal_gamma/0.05sim.39.3.csv")[,-1]

random.joint0.05<-sapply(sort(scenarios), function (k) sum(joint0.05[3*k,optimal.rec[[k]]]))
random.marginal0.05<-sapply(sort(scenarios), function (k) sum(marginal0.05[3*k,optimal.rec[[k]]]))

joint0.07<- read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/joint/0.07sim.39.3.csv")[,-1]
marginal0.07<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/delta/marginal_gamma/0.07sim.39.3.csv")[,-1]

random.joint0.07<-sapply(sort(scenarios), function (k) sum(joint0.07[3*k,optimal.rec[[k]]]))
random.marginal0.07<-sapply(sort(scenarios), function (k) sum(marginal0.07[3*k,optimal.rec[[k]]]))

M<-cbind(c(random.joint0.03, random.joint0.05, random.joint0.07, random.marginal0.03, random.marginal0.05, random.marginal0.07), rep(1:12, times=6), rep(rep(1:3, each=12), times=2), 
         rep(c(2,1), each=36))

colnames(M)<- c("pcs", "scen", "var", "cal_type")


pdf("figures/indiff.pdf", width=8, height=6)
ggplot(M, aes(x = as.numeric(scen), y = as.numeric(pcs), color = as.factor(cal_type), 
              shape = as.factor(var))) +
  geom_point(size = 4, stroke=1.5, alpha=0.6) +  # Adjust point size
  theme_minimal() +       # Use a minimal theme
  labs(title = " ",
       x = "Scenario",
       y = "Proportion",
       color = "Prior",
       shape = "Halfwidth")+ scale_x_continuous(breaks = 1:12, 
                                                labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))+
  scale_y_continuous(limits = c(0, 1))+
  scale_color_manual(values = c("#FC8D62", "#8DA0CB"),  labels = c("Marginal priors","Dose-agnostic priors"))+
  scale_shape_discrete(labels = c(expression(paste(delta[C], "= 0.03")),
                                  expression(paste(delta[C], "= 0.05")),expression(paste(delta[C], "= 0.07"))))+theme_minimal(base_size = 15)+
  theme(legend.position = "bottom",legend.direction = "vertical")
dev.off()

