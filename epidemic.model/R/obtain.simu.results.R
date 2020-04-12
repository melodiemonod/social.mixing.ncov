# Load packages
library(dplyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(gridExtra)

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

pred.data = subset(contact.data.agg, part.age.cat != "[90,95)" & part.age.cat != "[95,100]" & 
                                                       cont.age.cat != "[90,95)" & cont.age.cat != "[95,100]")


### DEMOGRAPHIC PARAMETERS ###
# age groups
age = sort(unique(pred.data$part.age.cat))
n_age = length(age)

# number of individuals by age group
T = pred.data[!duplicated(pred.data$w), ][order(pred.data[!duplicated(pred.data$w), ]$cont.age.cat),"w"]

### MODEL PARAMETERS ###
# number of infectious in pop at time t0
numb.inf.t0 = 100

# Length of observational period (in days)
t = 380 

# Date
min.date = "2020-03-05"; max.date = "2021-03-05"

################################ 0. Without intervention  ################################
set.seed(5678)
theta.draws.noi = lapply(1:100, function(x){
  draw.theta(intervention = 0)
})

res.noi = lapply(1:100, function(x){
  theta = theta.draws.noi[[x]]
  SEIIR.with.odin.intervention(intervention = 0, t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0)
})

SEIIR.res.noi = do.call(rbind, res.noi)
SEIIR.res.noi$iteration = rep(1:100, each = nrow(res.noi[[1]]))
SEIIR.res.noi$method = "no intervention"
max.time = max(SEIIR.res.noi$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.noi = subset(SEIIR.res.noi) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

save(theta.draws.noi, file = file.path(indir, "results", "theta.draws.noi.rda"))
save(SEIIR.res.noi, file = file.path(indir, "results", "SEIIR.res.noi.rda"))

################################ 1. Intervention 1 ################################
## prop = .9
set.seed(5678)
theta.draws.i1.p90 = lapply(1:100, function(x){
  draw.theta(intervention = 1, p = .9)
})

res.i1.p90 = lapply(1:100, function(x){
  theta = theta.draws.i1.p90[[x]]
  SEIIR.with.odin.intervention(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0, intervention = 1)
})

SEIIR.res.i1.p90 = do.call(rbind, res.i1.p90)
SEIIR.res.i1.p90$iteration = rep(1:100, each = nrow(res.i1.p90[[1]]))
SEIIR.res.i1.p90$method = "intervention 1, p = .9"
max.time = max(SEIIR.res.i1.p90$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i1.p90 = subset(SEIIR.res.i1.p90) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

# save
save(theta.draws.i1.p90, file = file.path(indir, "results", "theta.draws.i1.p90.rda"))
save(SEIIR.res.i1.p90, file = file.path(indir, "results", "SEIIR.res.i1.p90.rda"))

## prop = .6
set.seed(5678)
theta.draws.i1.p60 = lapply(1:100, function(x){
  draw.theta(intervention = 1, p = .6)
})

res.i1.p60 = lapply(1:100, function(x){
  theta = theta.draws.i1.p60[[x]]
  SEIIR.with.odin.intervention(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0, intervention = 1)
})

SEIIR.res.i1.p60 = do.call(rbind, res.i1.p60)
SEIIR.res.i1.p60$iteration = rep(1:100, each = nrow(res.i1.p60[[1]]))
SEIIR.res.i1.p60$method = "intervention 1, p = .6"
max.time = max(SEIIR.res.i1.p60$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i1.p60 = subset(SEIIR.res.i1.p60) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

# save
save(theta.draws.i1.p60, file = file.path(indir, "results", "theta.draws.i1.p60.rda"))
save(SEIIR.res.i1.p60, file = file.path(indir, "results", "SEIIR.res.i1.p60.rda"))

################################ 2. Intervention 2 ################################
set.seed(5678)
theta.draws.i2 = lapply(1:100, function(x){
  draw.theta(intervention = 2)
})

res.i2 = lapply(1:100, function(x){
  theta = theta.draws.i2[[x]]
  SEIIR.with.odin.intervention(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0, intervention = 2)
})

SEIIR.res.i2 = do.call(rbind, res.i2)
SEIIR.res.i2$iteration = rep(1:100, each = nrow(res.i2[[1]]))
SEIIR.res.i2$method = "intervention 2"
max.time = max(SEIIR.res.i2$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i2 = subset(SEIIR.res.i2) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

# save
save(theta.draws.i2, file = file.path(indir, "results", "theta.draws.i2.rda"))
save(SEIIR.res.i2, file = file.path(indir, "results", "SEIIR.res.i2.rda"))


################################ 3. Intervention 3 ################################
set.seed(5678)
theta.draws.i3 = lapply(1:100, function(x){
  draw.theta(intervention = 3)
})

res.i3 = lapply(1:100, function(x){
  theta = theta.draws.i3[[x]]
  SEIIR.with.odin.intervention(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0, intervention = 3)
})

SEIIR.res.i3 = do.call(rbind, res.i3)
SEIIR.res.i3$iteration = rep(1:100, each = nrow(res.i3[[1]]))
SEIIR.res.i3$method = "intervention 3"
max.time = max(SEIIR.res.i3$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i3 = subset(SEIIR.res.i3) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

# save
save(theta.draws.i3, file = file.path(indir, "results", "theta.draws.i3.rda"))
save(SEIIR.res.i3, file = file.path(indir, "results", "SEIIR.res.i3.rda"))


################################ 4. Intervention 4 ################################
# p1 = 0.9, p2 = 0.6
set.seed(5678)
theta.draws.i4.p190.p260 = lapply(1:100, function(x){
  draw.theta(intervention = 4, p1 = .9, p2 = .6)
})

res.i4.p190.p260 = lapply(1:100, function(x){
  theta = theta.draws.i4.p190.p260[[x]]
  SEIIR.with.odin.intervention4(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0)
})

SEIIR.res.i4.p190.p260 = do.call(rbind, res.i4.p190.p260)
SEIIR.res.i4.p190.p260$iteration = rep(1:100, each = nrow(res.i4.p190.p260[[1]]))
SEIIR.res.i4.p190.p260$method = "intervention 4, p1 = .9, p2 = .6"
max.time = max(SEIIR.res.i4.p190.p260$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i4.p190.p260 = subset(SEIIR.res.i4.p190.p260) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)
# save
save(theta.draws.i4.p190.p260, file = file.path(indir, "results", "theta.draws.i4.p190.p260.rda"))
save(SEIIR.res.i4.p190.p260, file = file.path(indir, "results", "SEIIR.res.i4.p190.p260.rda"))

# p1 = 0.9, p2 = 0.7
set.seed(5678)
theta.draws.i4.p190.p270 = lapply(1:100, function(x){
  draw.theta(intervention = 4, p1 = .9, p2 = .7)
})

res.i4.p190.p270 = lapply(1:100, function(x){
  theta = theta.draws.i4.p190.p270[[x]]
  SEIIR.with.odin.intervention4(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0)
})

SEIIR.res.i4.p190.p270 = do.call(rbind, res.i4.p190.p270)
SEIIR.res.i4.p190.p270$iteration = rep(1:100, each = nrow(res.i4.p190.p270[[1]]))
SEIIR.res.i4.p190.p270$method = "intervention 4, p1 = .9, p2 = .7"
max.time = max(SEIIR.res.i4.p190.p270$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i4.p190.p270 = subset(SEIIR.res.i4.p190.p270) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)
# save
save(theta.draws.i4.p190.p270, file = file.path(indir, "results", "theta.draws.i4.p190.p270.rda"))
save(SEIIR.res.i4.p190.p270, file = file.path(indir, "results", "SEIIR.res.i4.p190.p270.rda"))





#####
set.seed(5678)
theta.draws.i4.p190.p240 = lapply(1:100, function(x){
  draw.theta(intervention = 4, p1 = .8, p2 = .6)
})

res.i4.p190.p240 = lapply(1:100, function(x){
  theta = theta.draws.i4.p190.p240[[x]]
  SEIIR.with.odin.intervention4(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0)
})

SEIIR.res.i4.p190.p240 = do.call(rbind, res.i4.p190.p240)
SEIIR.res.i4.p190.p240$iteration = rep(1:100, each = nrow(res.i4.p190.p240[[1]]))
SEIIR.res.i4.p190.p240$method = "intervention 4, p1 = .9, p2 = .4"
max.time = max(SEIIR.res.i4.p190.p240$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i4.p190.p240 = subset(SEIIR.res.i4.p190.p240) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

# save
save(theta.draws.i4.p190.p240, file = file.path(indir, "results", "theta.draws.i4.p190.p240.rda"))
save(SEIIR.res.i4.p190.p240, file = file.path(indir, "results", "SEIIR.res.i4.p190.p240.rda"))


# p1 = 0.6, p2 = 0.2
set.seed(5678)
theta.draws.i4.p160.p220 = lapply(1:100, function(x){
  draw.theta(intervention = 4, p1 = 0.6, p2 = 0.2)
})

res.i4.p160.p220 = lapply(1:100, function(x){
  theta = theta.draws.i4.p160.p220[[x]]
  SEIIR.with.odin.intervention4(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0)
})

SEIIR.res.i4.p160.p220 = do.call(rbind, res.i4.p160.p220)
SEIIR.res.i4.p160.p220$iteration = rep(1:100, each = nrow(res.i4.p160.p220[[1]]))
SEIIR.res.i4.p160.p220$method = "intervention 4, p1 = .6, p2 = .2"
max.time = max(SEIIR.res.i4.p160.p220$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i4.p160.p220 = subset(SEIIR.res.i4.p160.p220) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

# save
save(theta.draws.i4.p160.p220, file = file.path(indir, "results", "theta.draws.i4.p160.p220.rda"))
save(SEIIR.res.i4.p160.p220, file = file.path(indir, "results", "SEIIR.res.i4.p160.p220.rda"))

# p1 = 0.9, p2 = 0.4
set.seed(5678)
theta.draws.i4.p190.p240 = lapply(1:100, function(x){
  draw.theta(intervention = 4, p1 = 0.9, p2 = 0.4)
})

res.i4.p190.p240 = lapply(1:100, function(x){
  theta = theta.draws.i4.p190.p240[[x]]
  SEIIR.with.odin.intervention4(t = t, theta = theta, T = T, numb.inf.t0 = numb.inf.t0)
})

SEIIR.res.i4.p190.p240 = do.call(rbind, res.i4.p190.p240)
SEIIR.res.i4.p190.p240$iteration = rep(1:100, each = nrow(res.i4.p190.p240[[1]]))
SEIIR.res.i4.p190.p240$method = "intervention 4, p1 = .9, p2 = .4"
max.time = max(SEIIR.res.i4.p190.p240$time.10th.death)
dates = seq((as.Date("2020-03-14")-max.time+1), by=1, len=365+max.time)
SEIIR.res.i4.p190.p240 = subset(SEIIR.res.i4.p190.p240) %>%
  group_by(age, T, variable,iteration) %>%
  mutate(date =  dates[max.time - time.10th.death + t]) %>%
  subset(date <= max.date & date >= min.date)

# save
save(theta.draws.i4.p190.p240, file = file.path(indir, "results", "theta.draws.i4.p190.p240.rda"))
save(SEIIR.res.i4.p190.p240, file = file.path(indir, "results", "SEIIR.res.i4.p190.p240.rda"))





