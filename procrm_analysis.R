setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors")
joint<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/joint/clay0.9.39.3.csv")[,-1]
marginal_gamma<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/marginal_gamma/clay0.9.39.3.csv")[,-1]
marginal_lognormal<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/marginal_lognormal/clay0.9.39.3.csv")[,-1]

cdlt.rec<-sapply(((0:35)*3)+1, function (k) which.min(abs(joint[k,]-0.25)))
pdlt.rec<-sapply(((0:35)*3)+2, function (k) which.min(abs(joint[k,]-0.35)))
optimal.rec<-pmin(cdlt.rec, pdlt.rec)


joint.rec<-sapply(1:36, function(k) joint[((k-1)*3)+3, optimal.rec[k]])
marginal.rec<-sapply(1:36, function(k) marginal_gamma[((k-1)*3)+3, optimal.rec[k]])



sd(joint.rec)
sqrt(var(marginal.rec))
mean(joint.rec)
mean(marginal.rec)


list.length<-sapply(optimal.rec, length)
index<-unlist(optimal.rec[list.length==1])
sd(tapply(joint.rec[list.length==1], index, mean))

sd(tapply(marginal.rec[list.length==1], index, mean))

#random scenario generation for paper 
set.seed(1)
scenarios<-sample(x=36, size=12, replace=FALSE)
rows<-sapply(1:12, function(k) joint[(3*(k-1)+1):3*k,])
row_index<-sort(c(3*(scenarios-1)+1, (3*(scenarios-1)+2), (3*scenarios)))
random.joint<-joint[row_index,]
random.marginal<-marginal_gamma[row_index,]
write.csv(random.joint, "tables/procrm.random.joint.csv")
write.csv(random.marginal, "tables/procrm.random.marginal.csv")

key<- data.frame(nrow=36, ncol=2)
ind<-1
for (i in 1:6){
  for(j in 1:6){
    key[ind,]<- c(ind,paste0("clay0.9.sim",i,".",j,".39.3.csv"))
  ind<-ind+1}
}
key[sort(scenarios),2]


joint_folder_path<-"C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/joint/"
marginal_folder_path<-"C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/marginal_gamma/"

for (i in 1:12){
  eval(parse(text=paste0('joint',i,'<-read.csv(paste0(joint_folder_path, key[sort(scenarios),2][',i,']))[,-1]')))
  eval(parse(text=paste0('marginal',i,'<-read.csv(paste0(marginal_folder_path, key[sort(scenarios),2][',i,']))[,-1]')))
}

prop.pat.joint<-NULL
prop.pat.marg<-NULL
correctdose<-optimal.rec[sort(scenarios)]
underdose<- sapply(sort(scenarios), function (k) (min(1,optimal.rec[k]-1):(optimal.rec[k]-1)))
overdose<- sapply(sort(scenarios), function (k) ((optimal.rec[k]+1):max(5,optimal.rec[k]+1)))

for (i in 1:12){
  eval(parse(text=paste0("under<-sum(colMeans(joint",i,"/rowSums(joint",i,"))[underdose[",i,"][[1]]])")))
  eval(parse(text=paste0("over<-sum(colMeans(joint",i,"/rowSums(joint",i,"))[overdose[",i,"][[1]]])")))
  eval(parse(text=paste0("best<-sum(colMeans(joint",i,"/rowSums(joint",i,"))[correctdose[",i,"]])")))
  row<-cbind.data.frame(best, under, over)
  prop.pat.joint<-rbind.data.frame(prop.pat.joint, row)
}

for (i in 1:12){
  eval(parse(text=paste0("under<-sum(colMeans(marginal",i,"/rowSums(marginal",i,"))[underdose[",i,"][[1]]])")))
  eval(parse(text=paste0("over<-sum(colMeans(marginal",i,"/rowSums(marginal",i,"))[overdose[",i,"][[1]]])")))
  eval(parse(text=paste0("best<-sum(colMeans(marginal",i,"/rowSums(marginal",i,"))[correctdose[",i,"]])")))
  row<-cbind.data.frame(best, under, over)
  prop.pat.marg<-rbind.data.frame(prop.pat.marg, row)
}

write.csv(prop.pat.joint, "tables/procrm.prop.pat.joint.csv")
write.csv(prop.pat.marg, "tables/procrm.prop.pat.marginal.csv")


                   