
joint<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/39.3/joint/20HN0.1.sim.39.3.csv")[,-1]
marginal<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/efftox_results/39.3/marginal_gamma/HN0.1/sim.39.3.csv")[,-1]


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

