# Load packages
library(dplyr)
library(reshape2)
library(data.table)

if(1){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
  code = "~/git/social.mixing.ncov/R"
}

# HPC
if(0){
  indir = "/rds/general/user/mm3218/home/projects/social.mixing.ncov"
  code = "/rds/general/user/mm3218/home/projects/social.mixing.ncov/simu"
}


### Load data ###
# social contact estimates and population size
load(file = file.path(indir, "results", "contact.data.agg.rda"))
load(file =file.path(indir, "results", "posterior.draw.aggregate.rda"))
load(file.path(indir, "results_simu", "mu.parameters.rda"))
load(file.path(indir, "results_simu", "rho.parameters.rda"))

source(file.path(code, "functions.simu.R"))

pred.data = posterior.draw.aggregate.subset = subset(contact.data.agg, part.age.cat != "[90,95)" & part.age.cat != "[95,100]" & 
                                                       cont.age.cat != "[90,95)" & cont.age.cat != "[95,100]")
  
  
### DEMOGRAPHIC PARAMETERS ###
# age groups
age = sort(unique(pred.data$part.age.cat))
n_age = length(age)

# number of individuals by age group
T = pred.data[!duplicated(pred.data$w), ][order(pred.data[!duplicated(pred.data$w), ]$cont.age.cat),"w"]

### MODEL PARAMETERS ###
# number of infectious in pop at time t0
numb.inf.t0 = 5

#number of days after lockdown
t = 40
t.total = t + 18



### RUN MODELS
set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.9,p2=.6)
})

res.p109.p206 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = numb.inf.t0)
})

SEIIR.res.p109.p206  = do.call(rbind, res.p109.p206)
SEIIR.res.p109.p206$iteration = rep(1:100, each = nrow(res.p109.p206[[1]]))

SEIIR.table.death.p109.p206 = SEIIR.res.p109.p206 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.9, p2 = 0.6")

##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.7,p2=.6)
})

res.p107.p206 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p107.p206 = do.call(rbind, res.p107.p206)
SEIIR.res.p107.p206$iteration = rep(1:100, each = nrow(res.p107.p206[[1]]))

SEIIR.table.death.p107.p206 = SEIIR.res.p107.p206 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.7, p2 = 0.6")
##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.6,p2=.6)
})

res.p106.p206 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p106.p206 = do.call(rbind, res.p106.p206)
SEIIR.res.p106.p206$iteration = rep(1:100, each = nrow(res.p106.p206[[1]]))

SEIIR.table.death.p106.p206 = SEIIR.res.p106.p206 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.6, p2 = 0.6")

##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.9,p2=.4)
})

res.p109.p204 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p109.p204 = do.call(rbind, res.p109.p204)
SEIIR.res.p109.p204$iteration = rep(1:100, each = nrow(res.p109.p204[[1]]))

SEIIR.table.death.p109.p204 = SEIIR.res.p109.p204 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.9, p2 = 0.4")

##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.9,p2=.2)
})

res.p109.p202 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p109.p202 = do.call(rbind, res.p109.p202)
SEIIR.res.p109.p202$iteration = rep(1:100, each = nrow(res.p109.p202[[1]]))

SEIIR.table.death.p109.p202 = SEIIR.res.p109.p202 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.9, p2 = 0.2")

##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.7,p2=.4)
})

res.p107.p204 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p107.p204 = do.call(rbind, res.p107.p204)
SEIIR.res.p107.p204$iteration = rep(1:100, each = nrow(res.p107.p204[[1]]))

SEIIR.table.death.p107.p204 = SEIIR.res.p107.p204 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.7, p2 = 0.4")

##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.7,p2=.2)
})

res.p107.p202 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p107.p202 = do.call(rbind, res.p107.p202)
SEIIR.res.p107.p202$iteration = rep(1:100, each = nrow(res.p107.p202[[1]]))

SEIIR.table.death.p107.p202 = SEIIR.res.p107.p202 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.7, p2 = 0.2")

##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.6,p2=.4)
})

res.p106.p204 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p106.p204 = do.call(rbind, res.p106.p204)
SEIIR.res.p106.p204$iteration = rep(1:100, each = nrow(res.p106.p204[[1]]))

SEIIR.table.death.p106.p204 = SEIIR.res.p106.p204 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.6, p2 = 0.4")

##

set.seed(5678)
theta.draws = lapply(1:100, function(x){
  draw.theta.intervention(p1=.6,p2=.2)
})

res.p106.p202 = lapply(1:100, function(x){
  SEIIR.with.odin.intervention(t = t, theta = theta.draws[[x]], T = T, numb.inf.t0 = 5)
})

SEIIR.res.p106.p202 = do.call(rbind, res.p106.p202)
SEIIR.res.p106.p202$iteration = rep(1:100, each = nrow(res.p106.p202[[1]]))

SEIIR.table.death.p106.p202 = SEIIR.res.p106.p202 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=t.total)[t], 
         method = "p1 = 0.6, p2 = 0.2")

SEIIR.table.death = rbind(SEIIR.table.death.p109.p206, SEIIR.table.death.p107.p206, SEIIR.table.death.p106.p206, 
                          SEIIR.table.death.p109.p204, SEIIR.table.death.p107.p204, SEIIR.table.death.p106.p204, 
                          SEIIR.table.death.p109.p202, SEIIR.table.death.p107.p202, SEIIR.table.death.p106.p202)
SEIIR.table.death$method = factor(SEIIR.table.death$method, 
                                    levels = c("p1 = 0.9, p2 = 0.6", "p1 = 0.7, p2 = 0.6", "p1 = 0.6, p2 = 0.6", 
                                               "p1 = 0.9, p2 = 0.4", "p1 = 0.7, p2 = 0.4", "p1 = 0.6, p2 = 0.4",
                                               "p1 = 0.9, p2 = 0.2", "p1 = 0.7, p2 = 0.2", "p1 = 0.6, p2 = 0.2"))
save(SEIIR.table.death, file = file.path(indir, "results", "SEIIR.table.death.rda"))


###
SEIIR.table.death = SEIIR.res.noi %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) %>%
  mutate(date = seq(as.Date("2020-03-05"), by=1, len=365)[t], 
         method = "p1 = 0.9, p2 = 0.6")

  