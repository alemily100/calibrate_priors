
joint<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/joint/20HN0.1.sim.39.3.csv")[,-1]
marginal<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/marginal_gamma/HN0.1/sim.39.3.csv")[,-1]


dlt.rec<-sapply(((0:53)*3)+1, function (k) which.min(abs(joint[k,]-0.25)))
collect.eff<-sapply(1:54, function (k) joint[((k-1)*3)+2,1:dlt.rec[[k]]])
optimal.rec<-sapply(1:54, function (k) which(collect.eff[[k]] == max(collect.eff[[k]])))


joint.rec<-sapply(1:54, function(k) sum(joint[((k-1)*3)+3, optimal.rec[[k]]]))
marginal.rec<-sapply(1:54, function(k) sum(marginal[((k-1)*3)+3, optimal.rec[[k]]]))

sqrt(var(joint.rec))
sqrt(var(marginal.rec))
mean(joint.rec)
mean(marginal.rec)
plot(1:length(marginal.rec), marginal.rec)
points(1:length(joint.rec), joint.rec, col="blue")


list.length<-sapply(optimal.rec, length)
index<-unlist(optimal.rec[list.length==1])
sd(tapply(joint.rec[list.length==1], index, mean))

sd(tapply(marginal.rec[list.length==1], index, mean))


hist(marginal.rec)

x<-joint.rec-marginal.rec
cols <- ifelse(x >= 0, "darkgreen", "red")

barplot(
  x,
  col = cols,
  border = NA
)

brks <- pretty(range(c(joint.rec, marginal.rec)), n = 30)

hist(joint.rec,
     breaks = brks,
     col = rgb(1, 0, 0, 0.4),
     border = "red",
     freq = FALSE,
     xlim = range(c(joint.rec, marginal.rec)),
     main = "Overlayed Histograms",
     xlab = "Value")

hist(marginal.rec,
     breaks = brks,
     col = rgb(0, 0, 1, 0.4),
     border = "blue",
     freq = FALSE,
     add = TRUE)

#random scenario generation for paper 
set.seed(1)
scenarios<-sample(x=54, size=12, replace=FALSE)
rows<-sapply(1:12, function(k) joint[(3*(k-1)+1):3*k,])
row_index<-sort(c(3*(scenarios-1)+1, (3*(scenarios-1)+2), (3*scenarios)))
random.joint<-joint[row_index,]
random.marginal<-marginal[row_index,]
write.csv(random.joint, "tables/efftox.random.joint.csv")
write.csv(random.marginal, "tables/efftox.random.marginal.csv")

key<- data.frame(nrow=54, ncol=2)
ind<-1
for (i in 1:6){
  for(j in 1:9){
    key[ind,]<- c(ind,paste0("20HN0.1.sim",i,".",j,".pat.allocation39.3.csv"))
    ind<-ind+1}
}

joint_folder_path<-"C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/joint/"

for (i in 1:12){
  eval(parse(text=paste0('joint',i,'<-read.csv(paste0(joint_folder_path, key[sort(scenarios),2][',i,']))[,-1]')))
}

key<- data.frame(nrow=54, ncol=2)
ind<-1
for (i in 1:6){
  for(j in 1:9){
    key[ind,]<- c(ind,paste0("sim",i,".",j,".pat.allocation39.3.csv"))
    ind<-ind+1}
}

marginal_folder_path<-"C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/marginal_gamma/HN0.1/"

for (i in 1:12){
  eval(parse(text=paste0('marginal',i,'<-read.csv(paste0(marginal_folder_path, key[sort(scenarios),2][',i,']))[,-1]')))
}

prop.pat.joint<-NULL
prop.pat.marg<-NULL
correctdose<-optimal.rec[sort(scenarios)]
underdose<- sapply(sort(scenarios), function (k) (min(1,min(optimal.rec[k][[1]])-1)):(min(optimal.rec[k][[1]])-1))
overdose<- sapply(sort(scenarios), function (k) ((max(optimal.rec[k][[1]])+1):max(5,max(optimal.rec[k][[1]]+1))))

for (i in 1:12){
  eval(parse(text=paste0("under<-sum(colMeans(joint",i,"/rowSums(joint",i,"))[underdose[",i,"][[1]]])")))
  eval(parse(text=paste0("over<-sum(colMeans(joint",i,"/rowSums(joint",i,"))[overdose[",i,"][[1]]])")))
  eval(parse(text=paste0("best<-sum(colMeans(joint",i,"/rowSums(joint",i,"))[correctdose[",i,"][[1]]])")))
  row<-cbind.data.frame(best, under, over)
  prop.pat.joint<-rbind.data.frame(prop.pat.joint, row)
}

for (i in 1:12){
  eval(parse(text=paste0("under<-sum(colMeans(marginal",i,"/rowSums(marginal",i,"))[underdose[",i,"][[1]]])")))
  eval(parse(text=paste0("over<-sum(colMeans(marginal",i,"/rowSums(marginal",i,"))[overdose[",i,"][[1]]])")))
  eval(parse(text=paste0("best<-sum(colMeans(marginal",i,"/rowSums(marginal",i,"))[correctdose[",i,"][[1]]])")))
  row<-cbind.data.frame(best, under, over)
  prop.pat.marg<-rbind.data.frame(prop.pat.marg, row)
}

write.csv(prop.pat.joint, "tables/efftox.prop.pat.joint.csv")
write.csv(prop.pat.marg, "tables/efftox.prop.pat.marginal.csv")



