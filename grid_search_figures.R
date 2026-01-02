library(tidyverse)
set.seed(1001)

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


KL<- function(vec, no.d, Hset_b, Hset_g, prior_vec){
  b_shape<-vec[1]
  g_shape<- vec[2]
  b_rate<- vec[3]
  g_rate<- vec[4]
  kl<-sum(sapply(1:5, function(k) mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g)*log((1/prior_vec[k])*mu_k(k,b_shape, g_shape,b_rate,g_rate, no.d, Hset_b, Hset_g))))
  return(kl)
}


setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/")

skeletonc<- c(0.06, 0.14, 0.25, 0.38, 0.50)
skeletonp<- c(0.10, 0.21, 0.35, 0.49, 0.61)
clin<-exp(crmsens(prior=skeletonc, target=0.25, model = "empiric")$Hset)
pat<- exp(crmsens(prior=skeletonp, target=0.35, model = "empiric")$Hset)


scenario<-2

sc1 <- c(0.35, 0.30, 0.20, 0.10, 0.05)
sc2 <- c(0.10, 0.15, 0.20, 0.25, 0.30)
sc3 <- c(0.10, 0.20, 0.40, 0.20, 0.10)
sc4 <- c(0.20, 0.20, 0.20, 0.20, 0.20)

no.d <- 5

sc_df <- tibble(
  cat = factor(1:no.d, labels = paste0(1:no.d)),
  sc = eval(parse(text=paste0("sc",scenario)))
)
# --- ggplot: white bars with black border, #2166AC semi-transparent overlay ---
ggplot(sc_df, aes(x = cat, y = sc)) +
  geom_col(fill = "#8DA0CB", alpha = 0.45,
           color="black", fill="white", alpha = 0, position = "identity") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "",
       x = "Dose", y = "Prior probability dose is MTD") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),text = element_text(size = 16))

# plot all function evaluations requi#2166AC

set.seed(2003)
scenario<-1
sc.grid <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/grid.sc',scenario,'.csv")')))
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))
ran.1.grid.c_shape<-runif(1, min = 0, max = 10)
ran.1.grid.c_rate<-runif(1, min = 0, max = 10)
ran.1.grid.p_shape<-runif(1, min = 0, max = 10)
ran.1.grid.p_rate<-runif(1, min = 0, max = 10)
ran.div<-eval(parse(text=paste0('KL(c(ran.1.grid.c_shape, ran.1.grid.p_shape, ran.1.grid.c_rate, ran.1.grid.p_rate), 5, clin, pat, sc',scenario,')')))
sc.grid[1,4:8]<-c(ran.1.grid.c_shape,ran.1.grid.c_rate,ran.1.grid.p_shape,ran.1.grid.p_rate, ran.div)
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))
divergence<-sapply(1:nrow(sc.grid), function (k) eval(parse(text=paste0('KL(c(sc.grid[',k,', 4], sc.grid[',k,', 6], sc.grid[',k,', 5], sc.grid[',k,', 7]), no.d, clin, pat,sc',scenario,')'))))

df_grid <- data.frame(method = "grid",n.functions = sc.grid$n.functions,divergence = log(divergence))
df_optim <- data.frame(method = "optim",n.functions = sc.optim$n.functions,divergence = log(sc.optim$divergence))
df <- rbind(df_grid, df_optim)
pdf("figures/div_sc1.pdf", height=4, width=5)
ggplot(df, aes(x = n.functions, y = divergence, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(grid = "#B2182B", optim = "#2166AC")) +
  labs(
    x = "Number of function evaluations (log scale)",
    y = "KL divergence (log scale)",
    color = NULL,
    title = " ",
  ) +
  scale_x_log10(
    labels = scales::label_number()
  )+
  theme_minimal(base_size = 14) +
  theme(text = element_text(size = 14),legend.position = "none")+
  scale_color_manual(
    values = c(grid = "#B2182B", optim = "#2166AC"),
    labels = c(grid = "Grid search", optim = "Joint calibration")
  )+ scale_y_continuous(limits = c(-6.5, 1))
dev.off()

##
scenario<-2
sc.grid <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/grid.sc',scenario,'.csv")')))
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))

ran.1.grid.c_shape<-runif(1, min = 0, max = 10)
ran.1.grid.c_rate<-runif(1, min = 0, max = 10)
ran.1.grid.p_shape<-runif(1, min = 0, max = 10)
ran.1.grid.p_rate<-runif(1, min = 0, max = 10)
ran.div<-eval(parse(text=paste0('KL(c(ran.1.grid.c_shape, ran.1.grid.p_shape, ran.1.grid.c_rate, ran.1.grid.p_rate), 5, clin, pat, sc',scenario,')')))
sc.grid[1,4:8]<-c(ran.1.grid.c_shape,ran.1.grid.c_rate,ran.1.grid.p_shape,ran.1.grid.p_rate, ran.div)
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))
divergence<-sapply(1:nrow(sc.grid), function (k) eval(parse(text=paste0('KL(c(sc.grid[',k,', 4], sc.grid[',k,', 6], sc.grid[',k,', 5], sc.grid[',k,', 7]), no.d, clin, pat,sc',scenario,')'))))

df_grid <- data.frame(method = "grid",n.functions = sc.grid$n.functions,divergence = log(divergence))
df_optim <- data.frame(method = "optim",n.functions = sc.optim$n.functions,divergence = log(sc.optim$divergence))
df <- rbind(df_grid, df_optim)
pdf("figures/div_sc2.pdf", height=4, width=5)
ggplot(df, aes(x = n.functions, y = divergence, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(grid = "#B2182B", optim = "#2166AC")) +
  labs(
    x = "Number of function evaluations (log scale)",
    y = "KL divergence (log scale)",
    color = NULL,
    title = " ",
  ) +
  scale_x_log10(
    labels = scales::label_number()
  )+
  theme_minimal(base_size = 14) +
  theme(text = element_text(size = 14),legend.position = "none")+
  scale_color_manual(
    values = c(grid = "#B2182B", optim = "#2166AC"),
    labels = c(grid = "Grid search", optim = "Joint calibration")
  )+ scale_y_continuous(limits = c(-6.5, 1))
dev.off()
##
scenario<-3
sc.grid <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/grid.sc',scenario,'.csv")')))
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))

ran.1.grid.c_shape<-runif(1, min = 0, max = 10)
ran.1.grid.c_rate<-runif(1, min = 0, max = 10)
ran.1.grid.p_shape<-runif(1, min = 0, max = 10)
ran.1.grid.p_rate<-runif(1, min = 0, max = 10)
ran.div<-eval(parse(text=paste0('KL(c(ran.1.grid.c_shape, ran.1.grid.p_shape, ran.1.grid.c_rate, ran.1.grid.p_rate), 5, clin, pat, sc',scenario,')')))
sc.grid[1,4:8]<-c(ran.1.grid.c_shape,ran.1.grid.c_rate,ran.1.grid.p_shape,ran.1.grid.p_rate, ran.div)
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))
divergence<-sapply(1:nrow(sc.grid), function (k) eval(parse(text=paste0('KL(c(sc.grid[',k,', 4], sc.grid[',k,', 6], sc.grid[',k,', 5], sc.grid[',k,', 7]), no.d, clin, pat,sc',scenario,')'))))

df_grid <- data.frame(method = "grid",n.functions = sc.grid$n.functions,divergence = log(divergence))
df_optim <- data.frame(method = "optim",n.functions = sc.optim$n.functions,divergence = log(sc.optim$divergence))
df <- rbind(df_grid, df_optim)
pdf("figures/div_sc3.pdf", height=4, width=5)
ggplot(df, aes(x = n.functions, y = divergence, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(grid = "#B2182B", optim = "#2166AC")) +
  labs(
    x = "Number of function evaluations (log scale)",
    y = "KL divergence (log scale)",
    color = NULL,
    title = " ",
  ) +
  scale_x_log10(
    labels = scales::label_number()
  )+
  theme_minimal(base_size = 14) +
  theme(text = element_text(size = 14),legend.position = "none")+
  scale_color_manual(
    values = c(grid = "#B2182B", optim = "#2166AC"),
    labels = c(grid = "Grid search", optim = "Joint calibration")
  )+ scale_y_continuous(limits = c(-6.5, 1))
dev.off()
##
scenario<-4
sc.grid <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/grid.sc',scenario,'.csv")')))
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))

ran.1.grid.c_shape<-runif(1, min = 0, max = 10)
ran.1.grid.c_rate<-runif(1, min = 0, max = 10)
ran.1.grid.p_shape<-runif(1, min = 0, max = 10)
ran.1.grid.p_rate<-runif(1, min = 0, max = 10)
ran.div<-eval(parse(text=paste0('KL(c(ran.1.grid.c_shape, ran.1.grid.p_shape, ran.1.grid.c_rate, ran.1.grid.p_rate), 5, clin, pat, sc',scenario,')')))
sc.grid[1,4:8]<-c(ran.1.grid.c_shape,ran.1.grid.c_rate,ran.1.grid.p_shape,ran.1.grid.p_rate, ran.div)
sc.optim <- eval(parse(text=paste0('read.csv("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/calibrate_priors/Statistics in Medicine/revision/code/calibrate_priors/results/gridsearch/optim.sc',scenario,'.csv")')))
divergence<-sapply(1:nrow(sc.grid), function (k) eval(parse(text=paste0('KL(c(sc.grid[',k,', 4], sc.grid[',k,', 6], sc.grid[',k,', 5], sc.grid[',k,', 7]), no.d, clin, pat,sc',scenario,')'))))

df_grid <- data.frame(method = "grid",n.functions = sc.grid$n.functions,divergence = log(divergence))
df_optim <- data.frame(method = "optim",n.functions = sc.optim$n.functions,divergence = log(sc.optim$divergence))
df <- rbind(df_grid, df_optim)
pdf("figures/div_sc4.pdf", height=4, width=5)
ggplot(df, aes(x = n.functions, y = divergence, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c(grid = "#B2182B", optim = "#2166AC")) +
  labs(
    x = "Number of function evaluations (log scale)",
    y = "KL divergence (log scale)",
    color = NULL,
    title = " ",
  ) +
  scale_x_log10(
    labels = scales::label_number()
  )+
  theme_minimal(base_size = 14) +
  theme(text = element_text(size = 14),legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c("right", "top"))+
  scale_color_manual(
    values = c(grid = "#B2182B", optim = "#2166AC"),
    labels = c(grid = "Grid search", optim = "Joint calibration")
  )+ scale_y_continuous(limits = c(-6.5, 1))
dev.off()