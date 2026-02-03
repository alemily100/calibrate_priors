library(ggplot2)
library(RColorBrewer)
joint<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/joint/20HN0.1.sim.39.3.csv")[,-1]

dlt.rec<-sapply(((0:53)*3)+1, function (k) which.min(abs(joint[k,]-0.25)))
collect.eff<-sapply(1:54, function (k) joint[((k-1)*3)+2,1:dlt.rec[[k]]])
optimal.rec<-sapply(1:54, function (k) which(collect.eff[[k]] == max(collect.eff[[k]])))


joint.rec<-sapply(1:54, function(k) sum(joint[((k-1)*3)+3, optimal.rec[[k]]]))

set.seed(1)
scenarios<-sample(x=54, size=12, replace=FALSE)
rows<-sapply(1:12, function(k) joint[(3*(k-1)+1):3*k,])
row_index<-sort(c(3*(scenarios-1)+1, (3*(scenarios-1)+2), (3*scenarios)))
random.joint<-joint[row_index,]

key20<- data.frame(nrow=54, ncol=2)
key16<- data.frame(nrow=54, ncol=2)
key40<- data.frame(nrow=54, ncol=2)
ind<-1
for (i in 1:6){
  for(j in 1:9){
    key20[ind,]<- c(ind,paste0("20HN0.1.sim",i,".",j,".pat.allocation39.3.csv"))
    key16[ind,]<- c(ind,paste0("16HN0.1.sim",i,".",j,".pat.allocation39.3.csv"))
    key40[ind,]<- c(ind,paste0("40HN0.1.sim",i,".",j,".pat.allocation39.3.csv"))
    ind<-ind+1}
}

joint_folder_path<-"C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/joint/"

for (i in 1:12){
  eval(parse(text=paste0('joint16.',i,'<-read.csv(paste0(joint_folder_path, key16[sort(scenarios),2][',i,']))[,-1]')))
  eval(parse(text=paste0('joint20.',i,'<-read.csv(paste0(joint_folder_path, key20[sort(scenarios),2][',i,']))[,-1]')))
  eval(parse(text=paste0('joint40.',i,'<-read.csv(paste0(joint_folder_path, key40[sort(scenarios),2][',i,']))[,-1]')))
}

prop.pat.joint16<-NULL
prop.pat.joint20<-NULL
prop.pat.joint40<-NULL
correctdose<-optimal.rec[sort(scenarios)]
underdose<- sapply(sort(scenarios), function (k) (min(1,min(optimal.rec[k][[1]])-1)):(min(optimal.rec[k][[1]])-1))
overdose<- sapply(sort(scenarios), function (k) ((max(optimal.rec[k][[1]])+1):max(5,max(optimal.rec[k][[1]]+1))))

for (i in 1:12){
  eval(parse(text=paste0("under<-sum(colMeans(joint16.",i,"/rowSums(joint16.",i,"))[underdose[",i,"][[1]]])")))
  eval(parse(text=paste0("over<-sum(colMeans(joint16.",i,"/rowSums(joint16.",i,"))[overdose[",i,"][[1]]])")))
  eval(parse(text=paste0("best<-sum(colMeans(joint16.",i,"/rowSums(joint16.",i,"))[correctdose[",i,"][[1]]])")))
  row<-cbind.data.frame(best, under, over)
  prop.pat.joint10<-rbind.data.frame(prop.pat.joint10, row)
  eval(parse(text=paste0("under<-sum(colMeans(joint20.",i,"/rowSums(joint20.",i,"))[underdose[",i,"][[1]]])")))
  eval(parse(text=paste0("over<-sum(colMeans(joint20.",i,"/rowSums(joint20.",i,"))[overdose[",i,"][[1]]])")))
  eval(parse(text=paste0("best<-sum(colMeans(joint20.",i,"/rowSums(joint20.",i,"))[correctdose[",i,"][[1]]])")))
  row<-cbind.data.frame(best, under, over)
  prop.pat.joint20<-rbind.data.frame(prop.pat.joint20, row)
  eval(parse(text=paste0("under<-sum(colMeans(joint40.",i,"/rowSums(joint40.",i,"))[underdose[",i,"][[1]]])")))
  eval(parse(text=paste0("over<-sum(colMeans(joint40.",i,"/rowSums(joint40.",i,"))[overdose[",i,"][[1]]])")))
  eval(parse(text=paste0("best<-sum(colMeans(joint40.",i,"/rowSums(joint40.",i,"))[correctdose[",i,"][[1]]])")))
  row<-cbind.data.frame(best, under, over)
  prop.pat.joint40<-rbind.data.frame(prop.pat.joint40, row)
}


M<-cbind(rbind(prop.pat.joint10,prop.pat.joint20,prop.pat.joint40), rep(1:12, times=3), rep(1:3, each=12))

colnames(M)<- c("correct", "under", "over", "scen", "category")

#eval(parse(text=paste0("write.csv(M, 'wages.over_under.",ncohort,".csv')")))

labels <-  c(expression(paste("Jointly calibrated, ", alpha[l], " = 0.82")), expression(paste("Jointly calibrated, ", alpha[l], " = 1.17")), expression(paste("Jointly calibrated, ", alpha[l], " = 1.88")))

setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")

pdf("figures/alpha_over.pdf", width=8, height=2.75)
ggplot(M, aes(x = as.numeric(scen), y = as.numeric(over), color = factor(category), shape = factor(category))) +
  geom_point(size = 4, stroke=2) + 
  theme_minimal() +   
  labs(title = " ",
       x = "Scenario",
       y = "Proportion",
       color = "Prior",
       shape = "Prior")+ scale_x_continuous(breaks = 1:12, 
                                            labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))+
  scale_y_continuous(limits = c(0, 1))+
  scale_color_manual(values=brewer.pal(3,"Set1"),labels = labels)+
  scale_shape_manual(values = c(3,5,6),labels = labels)+theme_minimal(base_size = 15)+
  theme(legend.position = "none")
dev.off()

pdf("figures/alpha_under.pdf", width=8, height=2.75)
ggplot(M, aes(x = as.numeric(scen), y = as.numeric(under), color = factor(category), shape = factor(category))) +
  geom_point(size = 4, stroke=2) +  # Adjust point size
  theme_minimal() +       # Use a minimal theme
  labs(title = " ",
       x = "Scenario",
       y = "Proportion",
       color = "Prior",
       shape = "Prior")+ scale_x_continuous(breaks = 1:12, 
                                            labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))+
  scale_y_continuous(limits = c(0, 1))+
  scale_color_manual(values=brewer.pal(3,"Set1"),labels = labels)+
  scale_shape_manual(values = c(3,5,6),labels = labels)+theme_minimal(base_size = 15)+
  theme(legend.position = "none")
dev.off()

pdf("figures/alpha_best.pdf", width=8, height=2.75)
ggplot(M, aes(x = as.numeric(scen), y = as.numeric(correct), color = factor(category), shape = factor(category))) +
  geom_point(size = 4, stroke=2) +  # Adjust point size
  theme_minimal() +       # Use a minimal theme
  labs(title = " ",
       x = "Scenario",
       y = "Proportion",
       color = "Prior",
       shape = "Prior")+ scale_x_continuous(breaks = 1:12, 
                                            labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))+
  scale_y_continuous(limits = c(0, 1))+
  scale_color_manual(values=brewer.pal(3,"Set1"),labels = labels)+
  scale_shape_manual(values = c(3,5,6), labels = labels)+theme_minimal(base_size = 15)
dev.off()
