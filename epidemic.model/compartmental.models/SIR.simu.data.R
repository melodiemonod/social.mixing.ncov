library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(gridExtra)
library("rstan")
library("rethinking")

if(1){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
  code = "~/git/social.mixing.ncov"
}

if(0){
  indir = "/rds/general/user/mm3218/home/projects/social.mixing.ncov"
  code = "/rds/general/user/mm3218/home/git/social.mixing.ncov"
}

source(file.path(code, "R", "stan_code.R"))

n_time = 14; n_age = 2; 
obs = c(1,5,10,75,225,300,260,240,190,125,60,30,20,0)
obs = c(cumsum(round(runif(25, 0, 5))), rev(cumsum(round(runif(15, 0, 5)))))
Z = matrix(nrow = n_time, ncol = n_age, rep(obs, n_age)) 
data_list = list(n_time = n_time, 
                 n_theta = 2, 
                 n_difeq = 3, 
                 n_age = n_age,
                 T = rep(764, n_age),
                 Z = Z, 
                 t0 = 0, 
                 ts = 1:n_time)

iter = 100; warmup = 50
fit = stan(model_code = model.epi.deterministic, data = data_list, chains = 1, iter = iter, warmup = warmup)
# saveRDS(fit, file=file.path(indir,'results_simu', "model.epi.deterministic_simu"))

summary(fit)
posterior = extract(fit)
S.m = matrix(ncol = n_age, nrow = n_time, 0); S.l95 = S.m; S.u95 = S.m
I.m = matrix(ncol = n_age, nrow = n_time, 0); I.l95 = S.m; I.u95 = S.m;
R.m = S.m; R.l95 = R.m; R.u95 = S.m;
for(t in 1:n_time){
  for(i in 1:n_age){
    S.m[t, i] = median(posterior$y_hat[,t, 1+3*(i-1)] * data_list$T[i])
    S.l95[t, i] = quantile(posterior$y_hat[,t, 1+3*(i-1)]* data_list$T[i], prob = 0.025)
    S.u95[t, i] = quantile(posterior$y_hat[,t, 1+3*(i-1)]* data_list$T[i], prob = 0.975)
    I.m[t, i] = median(posterior$y_hat[,t, 2+3*(i-1)]* data_list$T[i])
    I.l95[t, i] = quantile(posterior$y_hat[,t, 2+3*(i-1)]* data_list$T[i], prob = 0.025)
    I.u95[t, i] = quantile(posterior$y_hat[,t, 2+3*(i-1)]* data_list$T[i], prob = 0.975)
    R.m[t, i] = median(posterior$y_hat[,t, 3+3*(i-1)]* data_list$T[i])
    R.l95[t, i] = quantile(posterior$y_hat[,t, 3+3*(i-1)]* data_list$T[i], prob = 0.025)
    R.u95[t, i] = quantile(posterior$y_hat[,t, 3+3*(i-1)]* data_list$T[i], prob = 0.975)
  }
}


SIR.stan = data.table(t = rep(1:n_time, n_age*3), 
                      state = rep(c("S", "I", "R"), each = n_time*n_age),
                      median = c(as.vector(S.m), as.vector(I.m), as.vector(R.m)),
                      l95 = c(as.vector(S.l95), as.vector(I.l95), as.vector(R.l95)),
                      u95 = c(as.vector(S.u95), as.vector(I.u95), as.vector(R.u95)),
                      age = rep(rep(1:n_age, each = n_time), 3), 
                      X = as.vector(Z)) 

ggplot(subset(SIR.stan, state == "I"), aes(x = t, fill = state)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.5) +
  geom_line(aes(y = median, col = state)) +
  facet_wrap(~age) +
  geom_point(aes(y = Z))
