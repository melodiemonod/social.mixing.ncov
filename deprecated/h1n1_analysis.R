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

load(file.path(indir, "analyses", "contact_UK_H1N1.rda"))
load(file.path(indir, "analyses", "outbreak_UK_H1N1.rda"))
load(file.path(indir, "analyses", "demographic_UK_H1N1.rda"))
source(file.path(code, "R", "stan_code.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


### CONTACT MATRIX ANALYSIS WITHOUT OUTBREAK DATA ###
y = contact %>%
  group_by(AGE_i, AGE_j) %>%
  summarise(n_contact = sum(n_contact)) %>%
  dcast(AGE_i ~ AGE_j, value.var = "n_contact")
y[is.na(y)] <- 0
AGE_MIDPOINT = sort(demographic$AGE_MIDPOINT)
data_list = list(n_group = length(levels(contact$AGE_i)),
                 y = y[,-1],
                 agegroup_i = rep(AGE_MIDPOINT, each = length(AGE_MIDPOINT)), 
                 agegroup_j = rep(AGE_MIDPOINT, length(AGE_MIDPOINT)),
                 T = demographic[order(demographic$AGE_i),]$T,
                 N = c(contact_N[order(contact_N$AGE_i),]$N,0)) #careful about the 0, look if na
iter = 10000; warmup = 2000
fit = stan(model_code = model.transmission.agegroup, data = data_list, chains = 3, iter = iter, warmup = warmup, init = 0)
saveRDS(fit, file=file.path(indir,'results',"fit_contact_UK_5_h1n1.rds"))
fit = readRDS(file.path(indir,'results',"fit_contact_UK_5_h1n1.rds"))
mcmc_trace(fit, pars = c("theta", "beta", "rho", "alpha"))
posterior = extract(fit)
mu = as.vector(sapply(1:data_list$n_group, function(x){
  as.vector(posterior$mu[,x,])
}))
age_obs = levels(contact$AGE_i)
df = data.table(mu = mu, agegroup_i = rep(rep(age_obs, each = length(age_obs)), each = (iter-warmup)*3),
                agegroup_j = rep(rep(age_obs, length(age_obs)), each = (iter-warmup)*3)) %>%
  group_by(agegroup_i, agegroup_j) %>%
  summarise(mu.m = mean(mu), mu.l95 = quantile(mu, prob = .025), mu.u95 = quantile(mu, prob = .975))
df$agegroup_i = factor(df$agegroup_i, levels = c("0-4", "5-24", "25-64", "65+"))
df$agegroup_j = factor(df$agegroup_j, levels = c("0-4", "5-24", "25-64", "65+"))
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = mu.m))

c_est = as.vector(sapply(1:data_list$n_group, function(x){
  as.vector(posterior$c[,x,])
}))
df = data.table(c = c_est, agegroup_i = rep(rep(age_obs, each = length(age_obs)), each = (iter-warmup)*3),
                agegroup_j = rep(rep(age_obs, length(age_obs)), each = (iter-warmup)*3)) %>%
  group_by(agegroup_i, agegroup_j) %>%
  summarise(c.m = mean(c), c.l95 = quantile(c, prob = .025), c.u95 = quantile(c, prob = .975))
df$agegroup_i = factor(df$agegroup_i, levels = c("0-4", "5-24", "25-64", "65+"))
df$agegroup_j = factor(df$agegroup_j, levels = c("0-4", "5-24", "25-64", "65+"))
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = c.m))



## CONTACT MATRIX ANALAYSIS WITH OUTBREAK DATA AND CONTACT SURVEY
tmp = data.table(outbreak, t = 1:38) %>%
  melt(id.vars = "t")
ggplot(subset(tmp, t > 11), aes(x = t, y = value, color = variable)) + 
  geom_line()


## UNTIL BEGINNING OF NEW TERM 
n_time = 38-11; n_age = 4;
data_list = list(n_time = n_time, 
                 n_difeq = 3,
                 n_group = length(levels(contact$AGE_i)),
                 Y = y[,-1],
                 agegroup_i = rep(AGE_MIDPOINT, each = length(AGE_MIDPOINT)), 
                 agegroup_j = rep(AGE_MIDPOINT, length(AGE_MIDPOINT)),
                 T = c(demographic[demographic$AGE_i == "0-4", ]$T, demographic[demographic$AGE_i == "5-24",]$T, 
                       demographic[demographic$AGE_i == "25-64",]$T, demographic[demographic$AGE_i == "65+",]$T),
                 N = c(contact_N[contact_N$AGE_i == "0-4",]$N, contact_N[contact_N$AGE_i == "5-24",]$N, 
                       contact_N[contact_N$AGE_i == "25-64",]$N, contact_N[contact_N$AGE_i == "65+",]$N), 
                 Z = outbreak[12:(n_time+11),], 
                 t0 = 0, 
                 ts = 1:n_time)


y = polymod.tab.weekday %>%
  group_by(part.age, cont.age) %>%
  summarise(n_contact = sum(n_contact)) %>%
  dcast(AGE_i ~ AGE_j, value.var = "n_contact")
data_list = list(n_time = n_time, 
                 n_group = length(levels(contact$AGE_i)),
                 T = c(demographic[demographic$AGE_i == "0-4", ]$T, demographic[demographic$AGE_i == "5-24",]$T, 
                       demographic[demographic$AGE_i == "25-64",]$T, demographic[demographic$AGE_i == "65+",]$T),
                 Z = outbreak[12:(n_time+11),], 
                 t0 = 0, 
                 ts = 1:n_time, 
                 c)
iter = 10000; warmup = 1000
fit_overall = stan(model_code = model.overall.agegroup.deterministic, data = data_list, chains = 3, iter = iter, warmup = warmup, seed=4567, init= 0)
saveRDS(fit, file=file.path(outdir,"fit_overall_UK_H1N1.rds"))

fit_overall = readRDS(file.path(outdir1,"fit_overall_UK_H1N1_2.rds"))

mcmc_trace(fit_overall, pars = c("theta", "rho", "alpha"))

sum = summary(fit_overall)
sum$summary
posterior = extract(fit_overall)
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
                      Z = as.vector(unlist(outbreak[12:(n_time+11),]))) 

ggplot(subset(SIR.stan, state == "I"), aes(x = t, fill = state)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.5) +
  geom_line(aes(y = median, col = state)) +
  facet_wrap(~age) +
  geom_point(aes(y = Z))









y = contact %>%
  group_by(AGE_i, AGE_j) %>%
  summarise(n_contact = sum(n_contact)) %>%
  dcast(AGE_i ~ AGE_j, value.var = "n_contact")
AGE_MIDPOINT = demographic$AGE_MIDPOINT
n_time = 11; n_age = 4;
tmp = data.table(outbreak, t = 1:38) %>%
  melt(id.vars = "t")
data_list = list(n_time = n_time,
                 n_difeq = 3,
                 n_group = length(levels(contact$AGE_i)),
                 Y = y[,-1],
                 agegroup_i = rep(AGE_MIDPOINT, each = length(AGE_MIDPOINT)),
                 agegroup_j = rep(AGE_MIDPOINT, length(AGE_MIDPOINT)),
                 T = c(demographic[demographic$AGE_i == "0-4", ]$T, demographic[demographic$AGE_i == "5-24",]$T,
                       demographic[demographic$AGE_i == "25-64",]$T, demographic[demographic$AGE_i == "65+",]$T),
                 N = c(contact_N[contact_N$AGE_i == "0-4",]$N, contact_N[contact_N$AGE_i == "5-24",]$N,
                       contact_N[contact_N$AGE_i == "25-64",]$N, contact_N[contact_N$AGE_i == "65+",]$N),
                 Z = outbreak[1:11,],
                 t0 = 0,
                 ts = 1:n_time)
iter = 10000; warmup = 1000
fit_overall = stan(model_code = model.overall.agegroup.deterministic, data = data_list, chains = 3, iter = iter, warmup = warmup, seed=4567, init= 0)
saveRDS(fit_overall, file=file.path(indir,'results',"fit_overall_UK_H1N1.rds"))






