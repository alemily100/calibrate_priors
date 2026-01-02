setwd("/home/ealger/revision_calibrate_priors")

#args <- commandArgs(trailingOnly=TRUE)
#sc<-args[1]

sc <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(sc)) {
  stop("SLURM_ARRAY_TASK_ID is not set or not numeric")
}

message("Running scenario ", sc)

library(dfcrm)
sc1<- c(0.35, 0.30, 0.20, 0.10, 0.05)
sc2<- c(0.10, 0.15, 0.20, 0.25, 0.30)
sc3<- c(0.10, 0.20, 0.40, 0.20, 0.10)
sc4<- c(0.20, 0.20, 0.20, 0.20, 0.20)

prior_mat<- rbind(sc1, sc2, sc3, sc4)


skeletonc<- c(0.06, 0.14, 0.25, 0.38, 0.50)
skeletonp<- c(0.10, 0.21, 0.35, 0.49, 0.61)
clin<-exp(crmsens(prior=skeletonc, target=0.25, model = "empiric")$Hset)
pat<- exp(crmsens(prior=skeletonp, target=0.35, model = "empiric")$Hset)

KL_grid<- function(b_shape,b_rate, g_shape, g_rate, prior_vec){
  no.d<-5
  clin_param<-rgamma(100,b_shape,b_rate)
  pat_param<-rgamma(100, g_shape, g_rate)
  clin_rec<-sapply(1:100, function (k) which.min(abs(skeletonc^clin_param[k]-0.25)))
  pat_rec<-sapply(1:100, function (k) which.min(abs(skeletonp^pat_param[k]-0.35)))
  proposal<-table(pmin(clin_rec, pat_rec))/100
  kl<-sum(sapply(1:5, function(k) proposal[k]*log((1/prior_vec[k])*proposal[k])))
  return(kl)
}

KL_marginal<- function(shape, rate, no.d, Hset){
  kl<-sum(sapply(1:5, function(k) mu_marginal(k,shape, rate, no.d, Hset)*log(no.d*mu_marginal(k,shape, rate, no.d, Hset))))
  return(kl)
}

mu_k<-function(mtd, beta_shape, gamma_shape, beta_rate, gamma_rate, no.doses, Hset_beta, Hset_gamma){
  if(mtd==1){
    return(1-(1-pgamma(Hset_beta[mtd+1,1], shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd+1,1], shape=gamma_shape, rate=gamma_rate)))
  }
  if(mtd==no.doses){
    return((1-pgamma(Hset_beta[mtd,1],shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd,1], shape=gamma_shape, rate=gamma_rate)))
  }
  if(mtd>1 && mtd < no.doses){
    return((1-pgamma(Hset_beta[mtd,1],shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd,1], shape=gamma_shape, rate=gamma_rate)) -
             (1-pgamma(Hset_beta[mtd+1,1], shape=beta_shape, rate=beta_rate))*(1-pgamma(Hset_gamma[mtd+1,1], shape=gamma_shape, rate=gamma_rate)))
  }
}

mu_marginal<-function(mtd, shape, rate, no.d, Hset){
  if(mtd==1){
    return(pgamma(Hset[mtd+1,1], shape=shape, rate=rate))
  }
  if(mtd==no.d){
    return(1-pgamma(Hset[mtd,1], shape=shape, rate=rate))
  }
  if(mtd>1 && mtd < no.d){
    return(pgamma(Hset[mtd,2], shape=shape, rate=rate)-pgamma(Hset[mtd, 1], shape=shape, rate=rate))
  }
}

c.shape<- seq(from=0.01, to =3, length.out=100)
c.rate<- seq(from=0.01, to =3, length.out=100)
p.shape<- seq(from=0.01, to =3, length.out=100)
p.rate<- seq(from=0.01, to =3, length.out=100)

param_grid <- expand.grid(
  c.shape, c.rate, p.shape, p.rate
)


KL<- function(vec, no.d, Hset_b, Hset_g, prior_vec){
  b_shape<-vec[1]
  g_shape<- vec[2]
  b_rate<- vec[3]
  g_rate<- vec[4]
  kl<-sum(sapply(1:5, function(k) mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g)*log((1/prior_vec[k])*mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g))))
  return(kl)
}

collect_results<-matrix(nrow=2, ncol=8)
colnames(collect_results)<-c("type" ,"scenario", "C_shape", "C_rate", "P_shape", "P_rate", "divergence", "time_to_run")

set.seed(1001)
for(m in sc:sc){
  j<-m
  start_time_grid <- Sys.time()
  prior_vec<-prior_mat[j,]
  param_grid$score <- mapply(
    KL_grid,
    param_grid$Var1,
    param_grid$Var2,
    param_grid$Var3,
    param_grid$Var4,
    MoreArgs = list(prior_vec = prior_vec)
  )
  end_time_grid <- Sys.time()
  run_time_grid <- as.numeric(end_time_grid - start_time_grid)
  
  # Best parameters
  best_params <- param_grid[which.min(param_grid$score), ]
  input<-unlist(best_params)
  if(length(unlist(best_params))==0){
    input<-rep(NA, times=5)}
  collect_results[1,]<-c("grid", j, input, run_time_grid)
  
  ### optimisation in paper  
  rate<- seq(from=0, to =2, length.out=100)
  l <- data.frame(ac=NULL,ap=NULL, bc=NULL, bp=NULL, comb_div=NULL, clin_div=NULL, pat_div=NULL)
  start_time_opt <-Sys.time() 
  for(i in 1:100){
    nm<-optim(par=c(1,1,1,1), KL, no.d=5, Hset_b=clin, Hset_g=pat, prior_vec=prior_vec, lower=c(0,0.1,0,rate[i]), upper=c(10,10,10,rate[i]+0.01), method="L-BFGS-B")
    l<-rbind(l, data.frame(ac=nm$par[1],ap=nm$par[2],
                           bc=nm$par[3],bp=nm$par[4],comb_div=nm$value, clin_div=KL_marginal(nm$par[1],nm$par[3], 5, clin), pat_div=KL_marginal(nm$par[2],nm$par[4], 5, pat)))
  }
  best<-l[which.min(abs(l[,6]-l[,7])),]
  end_time_opt <- Sys.time()
  run_time_opt <- as.numeric(end_time_opt - start_time_opt)
  collect_results[2,]<-c("opt", j, unlist(c(best[1], best[3], best[2], best[4], best[5])), run_time_opt)
  write.csv(collect_results, paste0("/home/ealger/revision_calibrate_priors/results/gridsearch/sc",j,".csv"))
}

