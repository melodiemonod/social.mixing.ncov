# Load packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
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
source(file.path(code, "functions.simu.R"))

pred.data = contact.data.agg


### DEMOGRAPHIC PARAMETERS ###
# age groups
age = sort(unique(posterior.draw.aggregate$part.age.cat))
n_age = length(age)

# number of individuals by age group
T = pred.data[!duplicated(pred.data$w), ][order(pred.data[!duplicated(pred.data$w), ]$cont.age.cat),"w"]

### MODEL PARAMETERS ###
# number of infectious in pop at time t0
numb.inf.t0 = 1


################################ 0. Without intervention  ################################
# Length of observational period (in days)
t = 300

set.seed(5678)
theta.draws.noi = lapply(1:100, function(x){
  draw.theta(intervention = 0)
})

res.noi = lapply(1:100, function(x){
  theta.b.intervention = theta.draws.noi[[x]][[1]]
  theta.a.intervention = theta.draws.noi[[x]][[2]]
  SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =0)
})

SEIIR.res.noi = do.call(rbind, res.noi)

SEIIR.table.noi = SEIIR.res.noi %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))

SEIIR.res.noi$iteration = rep(1:100, each = nrow(res.noi[[1]]))
SEIIR.table2.noi = subset(SEIIR.res.noi, variable ==  "newI" | variable ==  "newIS" |variable =="newIA" | variable == "newDeath") %>%
  group_by(iteration, age, T, variable) %>%
  summarise(total = sum(value),
            total.abs = sum(value.abs)) %>%
  group_by(age, T, variable) %>%
  summarise(total.m = mean(total),
            total.l95 = quantile(total, prob = 0.025), 
            total.u95 = quantile(total, prob = 0.975),
            total.abs.m = mean(total.abs), 
            total.abs.l95 = quantile(total.abs, prob = 0.025), 
            total.abs.u95 = quantile(total.abs, prob = 0.975))

p1 = ggplot(subset(SEIIR.table.noi, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Proportion - SEIIR - no intervention") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age)+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p2 = ggplot(subset(SEIIR.table.noi, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - SEIIR - no intervention") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age)+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p3 = ggplot(subset(SEIIR.table.noi, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Absolute - Death - no intervention") +
  geom_ribbon(aes(ymin=value.l95, ymax = value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age)+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p4 = ggplot(subset(SEIIR.table.noi, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - Death - no intervention") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax = value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age)+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p = grid.arrange(p1, p2, p3, p4, nrow = 2)

p1.2 = ggplot(subset(SEIIR.table2.noi, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total infected - no intervention") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p2.2 = ggplot(subset(SEIIR.table2.noi, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total infected - no intervention") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p3.2 = ggplot(subset(SEIIR.table2.noi, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total death - no intervention") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p4.2 = ggplot(subset(SEIIR.table2.noi, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total death - no intervention") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p2=grid.arrange(p1.2,p2.2,p3.2,p4.2,nrow = 2)

# save
save(theta.draws.noi, file = file.path(indir, "results", "theta.draws.noi.rda"))
save(SEIIR.res.noi, file = file.path(indir, "results", "SEIIR.res.noi.rda"))
ggsave(p, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results", "SEIIR.noi.pdf"))
ggsave(p2, w = 14, h = 10, device = "pdf", file  = file.path(indir, "results", "total.noi.pdf"))

################################ 1. Intervention 1 ################################
# Length of observational period (in days)
t = 300

## prop = .5
set.seed(5678)
theta.draws.i1.p50 = lapply(1:100, function(x){
  draw.theta(intervention = 1, prop = .5)
})

res.i1.p50 = lapply(1:100, function(x){
  theta.b.intervention = theta.draws.i1.p50[[x]][[1]]
  theta.a.intervention = theta.draws.i1.p50[[x]][[2]]
  SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =1)
})

SEIIR.res.i1.p50 = do.call(rbind, res.i1.p50)

time.interention.range = range(SEIIR.res.i1.p50$time.i)

SEIIR.table.i1.p50 = select(SEIIR.res.i1.p50, -time.i) %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))

SEIIR.res.i1.p50$iteration = rep(1:100, each = nrow(res.i1.p50[[1]]))
SEIIR.table2.i1.p50 = subset(SEIIR.res.i1.p50, variable ==  "newI" | variable ==  "newIS" |variable =="newIA" | variable == "newDeath") %>%
  group_by(iteration, age, T, variable) %>%
  summarise(total = sum(value),
            total.abs = sum(value.abs)) %>%
  group_by(age, T, variable) %>%
  summarise(total.m = mean(total),
            total.l95 = quantile(total, prob = 0.025), 
            total.u95 = quantile(total, prob = 0.975),
            total.abs.m = mean(total.abs), 
            total.abs.l95 = quantile(total.abs, prob = 0.025), 
            total.abs.u95 = quantile(total.abs, prob = 0.975))

p1 = ggplot(subset(SEIIR.table.i1.p50, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Proportion - SEIIR - Itervention 1, p = 0.5") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed")+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p2 = ggplot(subset(SEIIR.table.i1.p50, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - SEIIR - Itervention 1, p = 0.5") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed")+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p3 = ggplot(subset(SEIIR.table.i1.p50, variable == "death"), aes(x = t)) +
  geom_line(aes(y =value.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.5") +
  geom_ribbon(aes(ymin=value.l95, ymax = value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p4 = ggplot(subset(SEIIR.table.i1.p50, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.5") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax = value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed")+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p = grid.arrange(p1, p2, p3, p4, nrow = 2)

p1.2 = ggplot(subset(SEIIR.table2.i1.p50, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total infected - Itervention 1, p = 0.5") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p2.2 = ggplot(subset(SEIIR.table2.i1.p50, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total infected - Itervention 1, p = 0.5") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p3.2 = ggplot(subset(SEIIR.table2.i1.p50, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total death - Itervention 1, p = 0.5") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p4.2 = ggplot(subset(SEIIR.table2.i1.p50, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total death - Itervention 1, p = 0.5") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p2=grid.arrange(p1.2,p2.2,p3.2,p4.2,nrow = 2)

# save
save(theta.draws.i1.p50, file = file.path(indir, "results", "theta.draws.i1.p50.rda"))
save(SEIIR.res.i1.p50, file = file.path(indir, "results", "SEIIR.res.i1.p50.rda"))
ggsave(p, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results", "SEIIR.i1.p50.pdf"))
ggsave(p2, w = 14, h = 10, device = "pdf", file  = file.path(indir, "results", "total.i1.p50.pdf"))

## prop = .2
set.seed(5678)
theta.draws.i1.p20 = lapply(1:100, function(x){
  draw.theta(intervention = 1, prop = .2)
})

res.i1.p20 = lapply(1:100, function(x){
  theta.b.intervention = theta.draws.i1.p20[[x]][[1]]
  theta.a.intervention = theta.draws.i1.p20[[x]][[2]]
  SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =1)
})

SEIIR.res.i1.p20 = do.call(rbind, res.i1.p20)

time.interention.range = range(SEIIR.res.i1.p20$time.i)

SEIIR.table.i1.p20 = select(SEIIR.res.i1.p20, -time.i) %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))

SEIIR.res.i1.p20$iteration = rep(1:100, each = nrow(res.i1.p20[[1]]))
SEIIR.table2.i1.p20 = subset(SEIIR.res.i1.p20, variable ==  "newI" | variable ==  "newIS" |variable =="newIA" | variable == "newDeath") %>%
  group_by(iteration, age, T, variable) %>%
  summarise(total = sum(value),
            total.abs = sum(value.abs)) %>%
  group_by(age, T, variable) %>%
  summarise(total.m = mean(total),
            total.l95 = quantile(total, prob = 0.025), 
            total.u95 = quantile(total, prob = 0.975),
            total.abs.m = mean(total.abs), 
            total.abs.l95 = quantile(total.abs, prob = 0.025), 
            total.abs.u95 = quantile(total.abs, prob = 0.975))

p1 = ggplot(subset(SEIIR.table.i1.p20, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Proportion - SEIIR - Itervention 1, p = 0.2") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed")+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p2 = ggplot(subset(SEIIR.table.i1.p20, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - SEIIR - Itervention 1, p = 0.2") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed")+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p3 = ggplot(subset(SEIIR.table.i1.p20, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.2") +
  geom_ribbon(aes(ymin=value.l95, ymax = value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed")+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p4 = ggplot(subset(SEIIR.table.i1.p20, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.2") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax = value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed")+
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p = grid.arrange(p1, p2, p3, p4, nrow = 2)

p1.2 = ggplot(subset(SEIIR.table2.i1.p20, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total infected - Itervention 1, p = 0.2") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2.2 = ggplot(subset(SEIIR.table2.i1.p20, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total infected - Itervention 1, p = 0.2") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p3.2 = ggplot(subset(SEIIR.table2.i1.p20, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total death - Itervention 1, p = 0.2") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p4.2 = ggplot(subset(SEIIR.table2.i1.p20, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total death - Itervention 1, p = 0.2") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])+
  theme(text = element_text(size=15))
p2=grid.arrange(p1.2,p2.2,p3.2,p4.2,nrow = 2)

# save
save(theta.draws.i1.p20, file = file.path(indir, "results", "theta.draws.i1.p20.rda"))
save(SEIIR.res.i1.p20, file = file.path(indir, "results", "SEIIR.res.i1.p20.rda"))
ggsave(p, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results", "SEIIR.i1.p20.pdf"))
ggsave(p2, w = 14, h = 10, device = "pdf", file  = file.path(indir, "results", "total.i1.p20.pdf"))

## prop = .1
# Length of observational period (in days)
t = 450
set.seed(5678)
theta.draws.i1.p10 = lapply(1:100, function(x){
  draw.theta(intervention = 1, prop = .1)
})

res.i1.p10 = lapply(1:100, function(x){
  theta.b.intervention = theta.draws.i1.p10[[x]][[1]]
  theta.a.intervention = theta.draws.i1.p10[[x]][[2]]
  SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =1)
})

SEIIR.res.i1.p10 = do.call(rbind, res.i1.p10)

time.interention.range = range(SEIIR.res.i1.p10$time.i)

SEIIR.table.i1.p10 = select(SEIIR.res.i1.p10, -time.i) %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))

SEIIR.res.i1.p10$iteration = rep(1:100, each = nrow(res.i1.p10[[1]]))
SEIIR.table2.i1.p10 = subset(SEIIR.res.i1.p10, variable ==  "newI" | variable ==  "newIS" |variable =="newIA" | variable == "newDeath") %>%
  group_by(iteration, age, T, variable) %>%
  summarise(total = sum(value),
            total.abs = sum(value.abs)) %>%
  group_by(age, T, variable) %>%
  summarise(total.m = mean(total),
            total.l95 = quantile(total, prob = 0.025), 
            total.u95 = quantile(total, prob = 0.975),
            total.abs.m = mean(total.abs), 
            total.abs.l95 = quantile(total.abs, prob = 0.025), 
            total.abs.u95 = quantile(total.abs, prob = 0.975))

p1 = ggplot(subset(SEIIR.table.i1.p10, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Proportion - SEIIR - Itervention 1, p = 0.1") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p2 = ggplot(subset(SEIIR.table.i1.p10, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - SEIIR - Itervention 1, p = 0.1") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p3 = ggplot(subset(SEIIR.table.i1.p10, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.1") +
  geom_ribbon(aes(ymin=value.l95, ymax = value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p4 = ggplot(subset(SEIIR.table.i1.p10, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.1") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax = value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p = grid.arrange(p1, p2, p3, p4, nrow = 2)

p1.2 = ggplot(subset(SEIIR.table2.i1.p10, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total infected - Itervention 1, p = 0.1") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2.2 = ggplot(subset(SEIIR.table2.i1.p10, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total infected - Itervention 1, p = 0.1") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p3.2 = ggplot(subset(SEIIR.table2.i1.p10, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total death - Itervention 1, p = 0.1") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p4.2 = ggplot(subset(SEIIR.table2.i1.p10, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total death - Itervention 1, p = 0.1") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2=grid.arrange(p1.2,p2.2,p3.2,p4.2,nrow = 2)

# save
save(theta.draws.i1.p10, file = file.path(indir, "results", "theta.draws.i1.p10.rda"))
save(SEIIR.res.i1.p10, file = file.path(indir, "results", "SEIIR.res.i1.p10.rda"))
ggsave(p, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results", "SEIIR.i1.p10.pdf"))
ggsave(p2, w = 14, h = 10, device = "pdf", file  = file.path(indir, "results", "total.i1.p10.pdf"))

## prop = .01
set.seed(5678)
theta.draws.i1.p01 = lapply(1:100, function(x){
  draw.theta(intervention = 1, prop = .01)
})

res.i1.p01 = lapply(1:100, function(x){
  theta.b.intervention = theta.draws.i1.p01[[x]][[1]]
  theta.a.intervention = theta.draws.i1.p01[[x]][[2]]
  SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =1)
})

SEIIR.res.i1.p01 = do.call(rbind, res.i1.p01)

time.interention.range = range(SEIIR.res.i1.p01$time.i)

SEIIR.table.i1.p01 = select(SEIIR.res.i1.p01, -time.i) %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value),
            value.l95 = quantile(value, prob = 0.025),
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs),
            value.abs.l95 = quantile(value.abs, prob = 0.025),
            value.abs.u95 = quantile(value.abs, prob = 0.975))

SEIIR.res.i1.p01$iteration = rep(1:100, each = nrow(res.i1.p01[[1]]))
SEIIR.table2.i1.p01 = subset(SEIIR.res.i1.p01, variable ==  "newI" | variable ==  "newIS" |variable =="newIA" | variable == "newDeath") %>%
  group_by(iteration, age, T, variable) %>%
  summarise(total = sum(value),
            total.abs = sum(value.abs)) %>%
  group_by(age, T, variable) %>%
  summarise(total.m = mean(total),
            total.l95 = quantile(total, prob = 0.025),
            total.u95 = quantile(total, prob = 0.975),
            total.abs.m = mean(total.abs),
            total.abs.l95 = quantile(total.abs, prob = 0.025),
            total.abs.u95 = quantile(total.abs, prob = 0.975))

p1 = ggplot(subset(SEIIR.table.i1.p01, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Proportion - SEIIR - Itervention 1, p = 0.01") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) +
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p2 = ggplot(subset(SEIIR.table.i1.p01, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - SEIIR - Itervention 1, p = 0.01") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) +
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p3 = ggplot(subset(SEIIR.table.i1.p01, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.01") +
  geom_ribbon(aes(ymin=value.l95, ymax = value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) +
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p4 = ggplot(subset(SEIIR.table.i1.p01, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - Death - Itervention 1, p = 0.01") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax = value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) +
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p = grid.arrange(p1, p2, p3, p4, nrow = 2)

p1.2 = ggplot(subset(SEIIR.table2.i1.p01, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total infected - intervention 1, p = 0.01") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2.2 = ggplot(subset(SEIIR.table2.i1.p01, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total infected - intervention 1, p = 0.01") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p3.2 = ggplot(subset(SEIIR.table2.i1.p01, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total death - intervention 1, p = 0.01") +
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p4.2 = ggplot(subset(SEIIR.table2.i1.p01, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total death - intervention 1, p = 0.01") +
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2=grid.arrange(p1.2,p2.2,p3.2,p4.2,nrow = 2)

# save
save(theta.draws.i1.p01, file = file.path(indir, "results", "theta.draws.i1.p01.rda"))
save(SEIIR.res.i1.p01, file = file.path(indir, "results", "SEIIR.res.i1.p01.rda"))
ggsave(p, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results", "SEIIR.i1.p01.pdf"))
ggsave(p2, w = 14, h = 10, device = "pdf", file  = file.path(indir, "results", "total.i1.p01.pdf"))


################################ 2. Intevention 2  ################################
# Length of observational period (in days)
t = 300

set.seed(5678)
theta.draws.i2 = lapply(1:100, function(x){
  draw.theta(intervention = 2)
})

res.i2 = lapply(1:100, function(x){
  theta.b.intervention = theta.draws.i2[[x]][[1]]
  theta.a.intervention = theta.draws.i2[[x]][[2]]
  SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =2)
})

SEIIR.res.i2 = do.call(rbind, res.i2)

time.interention.range = range(SEIIR.res.i2$time.i)

SEIIR.table.i2 = select(SEIIR.res.i2, -time.i) %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))

SEIIR.res.i2$iteration = rep(1:100, each = nrow(res.i2[[1]]))
SEIIR.table2.i2 = subset(SEIIR.res.i2, variable ==  "newI" | variable ==  "newIS" |variable =="newIA" | variable == "newDeath") %>%
  group_by(iteration, age, T, variable) %>%
  summarise(total = sum(value),
            total.abs = sum(value.abs)) %>%
  group_by(age, T, variable) %>%
  summarise(total.m = mean(total),
            total.l95 = quantile(total, prob = 0.025), 
            total.u95 = quantile(total, prob = 0.975),
            total.abs.m = mean(total.abs), 
            total.abs.l95 = quantile(total.abs, prob = 0.025), 
            total.abs.u95 = quantile(total.abs, prob = 0.975))

p1 = ggplot(subset(SEIIR.table.i2, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Proportion - SEIIR - intervention 2") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15))+
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p2 = ggplot(subset(SEIIR.table.i2, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - SEIIR - intervention 2") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15))+
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p3 = ggplot(subset(SEIIR.table.i2, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Absolute - Death - intervention 2") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15))+
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p4 = ggplot(subset(SEIIR.table.i2, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - Death - intervention 2") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax = value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  theme(text = element_text(size=15))+
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))
p = grid.arrange(p1, p2, p3, p4, nrow = 2)

p1.2 = ggplot(subset(SEIIR.table2.i2, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total infected - intervention 2") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2.2 = ggplot(subset(SEIIR.table2.i2, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total infected - intervention 2") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p3.2 = ggplot(subset(SEIIR.table2.i2, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total death - intervention 2") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p4.2 = ggplot(subset(SEIIR.table2.i2, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total death - intervention 2") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2=grid.arrange(p1.2,p2.2,p3.2,p4.2,nrow = 2)

# save
save(theta.draws.i2, file = file.path(indir, "results", "theta.draws.i2.rda"))
save(SEIIR.res.i2, file = file.path(indir, "results", "SEIIR.res.i2.rda"))
ggsave(p, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results", "SEIIR.i2.pdf"))
ggsave(p2, w = 14, h = 10, device = "pdf", file  = file.path(indir, "results", "total.i2.pdf"))


################################ 3. Intevention 3  ################################
set.seed(5678)
theta.draws.i3 = lapply(1:100, function(x){
  draw.theta(intervention = 3)
})

res.i3 = lapply(1:100, function(x){
  theta.b.intervention = theta.draws.i3[[x]][[1]]
  theta.a.intervention = theta.draws.i3[[x]][[2]]
  SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =3)
})

SEIIR.res.i3 = do.call(rbind, res.i3)

time.interention.range = range(SEIIR.res.i3$time.i)

SEIIR.table.i3 = select(SEIIR.res.i3, -time.i) %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
testE2 = subset(table.2, variable == "newIA" & age == "[20,25)")$value
SEIIR.res.i3$iteration = rep(1:100, each = nrow(res.i3[[1]]))
SEIIR.table2.i3 = subset(SEIIR.res.i3, variable ==  "newI" | variable ==  "newIS" |variable =="newIA" | variable == "newDeath") %>%
  group_by(iteration, age, T, variable) %>%
  summarise(total = sum(value),
            total.abs = sum(value.abs)) %>%
  group_by(age, T, variable) %>%
  summarise(total.m = mean(total),
            total.l95 = quantile(total, prob = 0.025), 
            total.u95 = quantile(total, prob = 0.975),
            total.abs.m = mean(total.abs), 
            total.abs.l95 = quantile(total.abs, prob = 0.025), 
            total.abs.u95 = quantile(total.abs, prob = 0.975))

p1 = ggplot(subset(SEIIR.table.i3, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Proportion - SEIIR - intervention 3") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))+
  theme(text = element_text(size=15))
p2 = ggplot(subset(SEIIR.table.i3, variable == "IS" | variable == "IA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - SEIIR - intervention 3") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))+
  theme(text = element_text(size=15))
p3 = ggplot(subset(SEIIR.table.i3, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.m, col = variable)) +
  ggtitle("Absolute - Death - intervention 3") +
  geom_ribbon(aes(ymin=value.l95, ymax =value.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))+
  theme(text = element_text(size=15))
p4 = ggplot(subset(SEIIR.table.i3, variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("Absolute - Death - intervention 3") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax = value.abs.u95, fill= variable), alpha = 0.55) +
  facet_wrap(~age) + 
  geom_vline(xintercept=time.interention.range[1], linetype="dashed") +
  geom_vline(xintercept=time.interention.range[2], linetype="dashed") +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0))+
  theme(text = element_text(size=15))
p = grid.arrange(p1, p2, p3, p4, nrow = 2)

p1.2 = ggplot(subset(SEIIR.table2.i3, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total infected - intervention 3") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2.2 = ggplot(subset(SEIIR.table2.i3, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total infected - intervention 3") +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p3.2 = ggplot(subset(SEIIR.table2.i3, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.l95, ymax = total.u95)) +
  geom_point(aes(y = total.m)) +
  ggtitle("Proportion - Total death - intervention 3") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p4.2 = ggplot(subset(SEIIR.table2.i3, variable == "newDeath"),aes(x = age, col = variable)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  geom_point(aes(y = total.abs.m)) +
  ggtitle("Absolute - Total death - intervention 3") + 
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(text = element_text(size=15))
p2=grid.arrange(p1.2,p2.2,p3.2,p4.2,nrow = 2)

# save
save(theta.draws.i3, file = file.path(indir, "results", "theta.draws.i3.rda"))
save(SEIIR.res.i3, file = file.path(indir, "results", "SEIIR.res.i3.rda"))
ggsave(p, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results", "SEIIR.i3.pdf"))
ggsave(p2, w = 14, h = 10, device = "pdf", file  = file.path(indir, "results", "total.i3.pdf"))
