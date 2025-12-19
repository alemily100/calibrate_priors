
joint<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/joint/clay0.9.18.3.csv")[,-1]
marginal_gamma<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/marginal_gamma/clay0.9.18.3.csv")[,-1]
marginal_lognormal<-read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/procrm_results/marginal_lognormal/clay0.9.18.3.csv")[,-1]


cdlt.rec<-sapply(((0:35)*3)+1, function (k) which.min(abs(joint[k,]-0.25)))
pdlt.rec<-sapply(((0:35)*3)+2, function (k) which.min(abs(joint[k,]-0.35)))
optimal.rec<-pmin(cdlt.rec, pdlt.rec)


joint.rec<-sapply(1:36, function(k) joint[((k-1)*3)+3, optimal.rec[k]])
marginal.rec<-sapply(1:36, function(k) marginal_gamma[((k-1)*3)+3, optimal.rec[k]])



sqrt(var(joint.rec))
sqrt(var(marginal.rec))
mean(joint.rec)
mean(marginal.rec)
