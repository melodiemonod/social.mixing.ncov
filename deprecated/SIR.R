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

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load(file.path(indir, "analyses", "outbreak_10.rda"))
load(file.path(indir, "analyses", "demographic_10.rda"))
source(file.path(code, "R", "stan_code.R"))

ggplot(outbreak, aes(x = t, y = Z)) +
  geom_point() +
  facet_wrap(~ AGE)

# From 19/01 to 27/02 
df = outbreak %>%
  select(t, Z, AGE) %>%
  dcast(t ~ AGE, value.var = "Z")

n_time = 40; n_age = 11; 
#inster row for 21/02 and 26/02 (days without observation)
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}
df = insertRow(df, c("2020-02-21", rep(NA, n_age)), 31-19+1+21)
df = insertRow(df, c("2020-02-26", rep(NA, n_age)), 31-19+1+26)
for(age in 1:n_age){
  df[1,age+1] = ifelse(is.na(df[1,age+1]), 0, df[1,age+1]) 
  for(t in 2:n_time){
    df[t,age+1] = ifelse(is.na(df[t,age+1]), df[t-1,age+1], df[t,age+1]) 
  }
}

data_list = list(n_time = n_time, 
                 n_difeq = 3, 
                 n_age = n_age,
                 T =  demographic$T,
                 Z = df[,-1], 
                 t0 = 0, 
                 ts = 1:n_time)
iter = 2000; warmup = 500
fit = stan(model_code = model.epi.deterministic, data = data_list, chains = 3, iter = iter, warmup = warmup, seed = 456789)
saveRDS(fit, file = file.path(outdir,'results', "SIR.mixing.deterministic_10.rds"))


n_time = 40; n_age = 4; 
obs = c(1,5,10,75,225,300,260,240,190,125,60,30,20,0)
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
fit = stan(model_code = model.epi.deterministic, data = data_list, chains = 1, iter = iter, warmup = warmup, seed = 45678, init = 0)
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
                      Z = as.vector(Z)) 

ggplot(subset(SIR.stan, state == "I"), aes(x = t, fill = state)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.5) +
  geom_line(aes(y = median, col = state)) +
  facet_wrap(~age) +
  geom_point(aes(y = Z))
