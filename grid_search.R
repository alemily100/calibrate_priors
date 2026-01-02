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
  clin_param<-rgamma(1000,b_shape,b_rate)
  pat_param<-rgamma(1000, g_shape, g_rate)
  clin_rec<-sapply(1:1000, function (k) which.min(abs(skeletonc^clin_param[k]-0.25)))
  pat_rec<-sapply(1:1000, function (k) which.min(abs(skeletonp^pat_param[k]-0.35)))
  proposal<-table(pmin(clin_rec, pat_rec))/1000
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


KL<- function(vec, no.d, Hset_b, Hset_g, prior_vec){
  b_shape<-vec[1]
  g_shape<- vec[2]
  b_rate<- vec[3]
  g_rate<- vec[4]
  kl<-sum(sapply(1:5, function(k) mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g)*log((1/prior_vec[k])*mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g))))
  return(kl)
}

collect_results_grid<-matrix(nrow=5, ncol=8)
colnames(collect_results_grid)<-c("type" ,"scenario", "C_shape", "C_rate", "P_shape", "P_rate", "divergence", "time_to_run")

collect_results_optim<-matrix(nrow=5, ncol=8)
colnames(collect_results_optim)<-c("type" ,"scenario", "C_shape", "C_rate", "P_shape", "P_rate", "divergence", "time_to_run")

to_explore_grid<-c(2,3, 5, 10, 15, 20) 
to_explore_optim<-c(5, 16,20,25,50,75)
set.seed(1001)
for(m in sc:sc){
  j<-m
  for(length in 1:length(to_explore_grid)){
  prior_vec<-prior_mat[j,]
  c.shape<- seq(from=0.01, to =3, length.out=to_explore_grid[length])
  c.rate<- seq(from=0.01, to =3, length.out=to_explore_grid[length])
  p.shape<- seq(from=0.01, to =3, length.out=to_explore_grid[length])
  p.rate<- seq(from=0.01, to =3, length.out=to_explore_grid[length])
  param_grid <- expand.grid(
    c.shape, c.rate, p.shape, p.rate
  )
  t.grid<-system.time(
  param_grid$score <- mapply(
    KL_grid,
    param_grid$Var1,
    param_grid$Var2,
    param_grid$Var3,
    param_grid$Var4,
    MoreArgs = list(prior_vec = prior_vec)
  )
  )
  run_time_grid <- as.numeric(t.grid["elapsed"])
  # Best parameters
  best_params <- param_grid[which.min(param_grid$score), ]
  input<-unlist(best_params)
  if(length(unlist(best_params))==0){
    input<-rep(NA, times=5)}
  collect_results_grid[length,]<-c("grid", j, input, run_time_grid)
  write.csv(collect_results_grid, paste0("/home/ealger/revision_calibrate_priors/results/gridsearch/grid.sc",j,".csv"))
  }
  for(max.fn in to_explore_optim){
  ### optimisation in paper  
  rate<- seq(from=0, to =2, length.out=100)
  l <- data.frame(ac=NULL,ap=NULL, bc=NULL, bp=NULL, comb_div=NULL, clin_div=NULL, pat_div=NULL)
  start_par <- c(1,1,1,1)
  t.opt<-system.time(
  nm <- nlminb(
      start = start_par,
      objective = KL,
      lower = c(0, 0.1, 0, rate[i]),
      upper = c(10, 10, 10, rate[i] + 0.01),
      control = list(iter.max = max.fn,eval.max = 5000, rel.tol = 1e-8),
      no.d = 5,Hset_b = clin,Hset_g = pat,prior_vec = prior_vec)
  )
  run_time_opt <- as.numeric(t.opt["elapsed"])
  collect_results_optim[which(max.fn == to_explore_optim),]<-c("opt", j, unlist(nm$par), nm$objective, run_time_opt)
  write.csv(collect_results_optim, paste0("/home/ealger/revision_calibrate_priors/results/gridsearch/optim.sc",j,".csv"))
  }
}

