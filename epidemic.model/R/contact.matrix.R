library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(gridExtra)
library("rstan")
library("rethinking")

if(0){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
  code = "~/git/social.mixing.ncov"
}

if(1){
  indir = "/rds/general/user/mm3218/home/projects/social.mixing.ncov"
  code = "/rds/general/user/mm3218/home/git/social.mixing.ncov"
}


load(file.path(indir, "analyses", "contact.tab_UK10.rda"))
source(file.path(code, "R", "stan_code.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### age analysis ###
y = contact.tab %>%
  reshape2::dcast(part.age ~ cont.age, value.var = "y")

U = contact.tab %>%
  reshape2::dcast(part.age ~ cont.age, value.var = "U")
  
N_group = length(unique(contact.tab$cont.age))
age = unique(contact.tab$cont.age)
n_group = N_group - 2
age.midpoint = unique(contact.tab$cont.age.midpoint)

data_list = list(N_group = N_group,
                 n_group = n_group,
                 y = as.matrix(y[,-1]),
                 agegroup_i = rep(age.midpoint, each = N_group),
                 agegroup_j = rep(age.midpoint, N_group),
                 logU = log(as.matrix(U[,-1])))

iter = 10000; warmup = 2000
fit = stan(model_code = model.contact.matrix1, data = data_list, chains = 3, iter = iter, warmup = warmup)
saveRDS(fit, file=file.path(indir,'results',"fit_contact_UK_10.rds"))
fit = readRDS(file.path(indir,'results',"fit_contact_UK_10.rds"))
  
## Effective sample size and Gelman-Rubin convergence coefficient 
sum = summary(fit)
neff = sum$summary[,9]; r_1 = sum$summary[,10]
n.param = 1 + 1+ 1+1+N_group*N_group
diagnostics = matrix(nrow = 2, ncol = 2, c(min(neff[1:n.param]), max(neff[1:n.param]), min(r_1[1:n.param]), max(r_1[1:n.param])))
colnames(diagnostics) = c("neff", "R_hat")
rownames(diagnostics) = c("min", "max")

## CI of estimator
posterior = extract(fit)
# Divide c by 1e8 to go back to original scale (  deflates contact rate c)
c_est = as.vector(sapply(1:data_list$N_group, function(x){
  as.vector(posterior$c[,x,])/1e8
}))
df = data.table(c = c_est, 
                part.age = rep(rep(age, each = N_group), each = (iter-warmup)*3),
                cont.age = rep(rep(age, N_group), each = (iter-warmup)*3)) %>%
  group_by(part.age, cont.age) %>%
  summarise(c.m = mean(c), c.l95 = quantile(c, prob = .025), c.u95 = quantile(c, prob = .975), 
            c.sd = sd(c))
df <- with(df, df[order(cont.age, part.age), ])
save(df, file = file.path(indir, "results", "c.STAN"))


E_est = as.vector(sapply(1:data_list$N_group, function(x){
  as.vector(posterior$E_Y[,x,])
}))

df = data.table(E = E_est, 
                part.age = rep(rep(age, each = N_group), each = (iter-warmup)*3),
                cont.age = rep(rep(age, N_group), each = (iter-warmup)*3)) %>%
  group_by(part.age, cont.age) %>%
  summarise(E.m = mean(E), E.l95 = quantile(E, prob = .025), E.u95 = quantile(E, prob = .975), 
            E.sd = sd(E))
