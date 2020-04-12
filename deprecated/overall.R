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

load(file.path(indir, "analyses", "contact_10.rda"))
load(file.path(indir, "analyses", "outbreak_10.rda"))
load(file.path(indir, "analyses", "demographic_10.rda"))
source(file.path(code, "R", "stan_code.R"))


## prepare contact data
y = contact %>%
  select(AGE_i, c.age.under5, c.age.6to19, c.age.20to64, c.age.over65) %>%
  group_by(AGE_i) %>%
  summarise("under5" = sum(c.age.under5), 
            "6to19" = sum(c.age.6to19), 
            "20to64" = sum(c.age.20to64),
            "over65" = sum(c.age.over65))
tmp = contact %>%
  group_by(AGE) %>%
  summarise(N = n())
AGE_MIDPOINT = demographic$AGE_MIDPOINT
CORRESPONDANCE = matrix(nrow = length(levels(contact$AGE_i)), ncol = 2, c(1, 1, 2, 3, 4, 8, 9, 11), byrow=T)

## prepare outbreak data
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
                 n_group = length(levels(contact$AGE_i)),
                 n_age = length(levels(contact$AGE)),
                 Y = y[,-1],
                 age_i = rep(AGE_MIDPOINT, each = length(AGE_MIDPOINT)), 
                 age_j = rep(AGE_MIDPOINT, length(AGE_MIDPOINT)),
                 corr = CORRESPONDANCE,
                 T = demographic$T,
                 N = tmp$N, 
                 Z = df[,-1], 
                 t0 = 0, 
                 ts = 1:n_time)

