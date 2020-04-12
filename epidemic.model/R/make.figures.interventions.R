# Load packages
library(dplyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(gridExtra)

if(0){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
  code = "~/git/social.mixing.ncov/R"
}

# HPC
if(1){
  indir = "/rds/general/user/mm3218/home/projects/social.mixing.ncov"
  code = "/rds/general/user/mm3218/home/projects/social.mixing.ncov/simu"
}


load(file =file.path(indir, "results", "contact.data.agg.rda"))
pred.data = subset(contact.data.agg, part.age.cat != "[90,95)" & part.age.cat != "[95,100]" & 
                                                       cont.age.cat != "[90,95)" & cont.age.cat != "[95,100]")

### DEMOGRAPHIC PARAMETERS ###
# age groups
age = sort(unique(pred.data$part.age.cat))
n_age = length(age)

# number of individuals by age group
T = pred.data[!duplicated(pred.data$w), ][order(pred.data[!duplicated(pred.data$w), ]$cont.age.cat),"w"]
total.pop = sum(T)

### MODEL PARAMETERS ###
# number of infectious in pop at time t0
numb.inf.t0 = 100

# Length of observational period (in days)
t = 380 

# Date
min.date = "2020-03-05"; max.date = "2021-03-05"

## RESULTS FROM SIMULATION
load(file.path(indir, "results", "SEIIR.res.noi.rda"))
load(file.path(indir, "results", "SEIIR.res.i1.p90.rda"))
load(file.path(indir, "results", "SEIIR.res.i1.p60.rda"))
load(file.path(indir, "results", "SEIIR.res.i2.rda"))
load(file.path(indir, "results", "SEIIR.res.i3.rda"))
load(file.path(indir, "results", "SEIIR.res.i4.p190.p270.rda"))
load(file.path(indir, "results", "SEIIR.res.i4.p190.p260.rda"))

SEIIR.res = rbind(SEIIR.res.i1.p90, SEIIR.res.i1.p60, SEIIR.res.i2, SEIIR.res.i3, SEIIR.res.i4.p190.p260, SEIIR.res.i4.p190.p270)

SEIIR.res$method = ordered(SEIIR.res$method, 
                            levels = c("intervention 1, p = .9", "intervention 1, p = .6", "intervention 2", "intervention 3", 
                                       "intervention 4, p1 = .9, p2 = .7", "intervention 4, p1 = .9, p2 = .6"))


# NEW CASES, CUMULATIVE CASES & NEW DEATHS, CUULATIVE DEATHS
SEIIR.table.cum = subset(SEIIR.res, variable == "newI" | variable == "newIS" | variable == "newIA") %>%
  group_by(age,T,variable,iteration, method) %>%
  mutate(value = cumsum(value),
         value.abs = cumsum(value.abs)) %>%
  group_by(date,variable,iteration, method) %>%
  summarise(value = sum(value),
         value.abs = sum(value.abs)) %>%
  group_by(date, variable, method) %>%
  summarise(value.m = mean(value),
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975)) 
SEIIR.table.cum$variable = ordered(SEIIR.table.cum$variable, 
                                    levels = c("newIS", "newIA", "newI"))
SEIIR.table.cum$variable = factor(SEIIR.table.cum$variable, 
                                   labels=c("total IS","total IA","total I"))
SEIIR.table.new = SEIIR.res %>%
  group_by(date, variable, method) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975)) 
SEIIR.table.new$variable = ordered(SEIIR.table.new$variable, 
                                   levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
SEIIR.table.newandcumI = rbind(SEIIR.table.cum, subset(SEIIR.table.new, variable == "newI" | variable == "newIS" | variable == "newIA")) 
labels_facet = c(`newI` = "New cases", 
                 `total I` = "Cumulative cases")
p1 = ggplot(subset(SEIIR.table.newandcumI, variable == "newI"| variable == "total I"), aes(x = date, y = value.abs.m/1000000)) +
  geom_line(aes(col = method)) +
  facet_wrap(~variable, nrow =2, scales = "free", labeller = as_labeller(labels_facet)) +
  theme_bw()+
  scale_x_date(date_breaks = "2 months" , date_labels = "%b-%y") +
  scale_y_continuous("Value (in millions)", sec.axis = sec_axis(~ . *1000000/ total.pop, name = "Proportion relative to total population")) +
  theme(legend.position = "none")+ 
  scale_color_brewer(palette="Dark2")
labels_facet = c(`newdeath` = "New deaths", 
                 `death` = "Cumulative deaths")
p2 = ggplot(data = NULL, aes(x = date)) +
  geom_line(data = subset(SEIIR.table.new, variable == "newdeath"), aes(col = method, y = value.abs.m/1000)) +
  geom_line(data = subset(SEIIR.table.new, variable == "death"), aes(col = method, y = (value.abs.m)/1000)) +
  facet_wrap(~variable, nrow =2, scales = "free", labeller = as_labeller(labels_facet)) +
  theme_bw()+
  scale_x_date(date_breaks = "2 months" , date_labels = "%b-%y") +
  scale_y_continuous("Value (in thousands)", sec.axis = sec_axis(~ . *1000 / total.pop, name = "Proportion relative to total population")) + 
  scale_color_brewer(palette="Dark2")
pp= grid.arrange(p1,p2, nrow = 1, widths = c(1,1.55))
ggsave(pp, w =11, h = 6, device = "pdf", file  = file.path(indir, "results", "I.death.comp.pdf"))

## TOTAL DEATHS BY AGE AGROUPS
SEIIR.table = subset(SEIIR.res, variable == "death" & date == max.date)%>%
  group_by(variable, method,age,T) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
num.method = length(unique(SEIIR.table$method))
labels_facet = c(`intervention 1, p = .6` = "intervention 1",
                 `intervention 1, p = .9` = "intervention 1",
                 `intervention 2` = "intervention 2",
                 `intervention 3` = "intervention 3",
                 `intervention 4, p1 = .9, p2 = .7` = "intervention 4",
                 `intervention 4, p1 = .9, p2 = .6` = "intervention 4")
p =ggplot(data = SEIIR.table, aes(x = age, y = value.m, col = method)) +
  geom_errorbar(aes(ymin = value.l95, ymax = value.u95)) +
  theme_bw() +
  facet_wrap(~method, ncol = num.method, labeller = as_labeller(labels_facet))+
  scale_x_discrete(breaks=age[c(TRUE, FALSE)]) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom") +
  ylab("Proportion") +
  xlab("age group")+ 
  scale_color_brewer(palette="Dark2")
ggsave(p, w =10, h = 6, device = "pdf", file  = file.path(indir, "results", "death.age.comp.pdf"))

## REDUCTION IN TOTAL DEATHS AND TOTAL CASES COMPARED TO NO INTERVENTION
SEIIR.table.noi.cum = subset(SEIIR.res.noi, variable == "newI") %>%
  group_by(age, T, variable,iteration) %>%
  mutate(value.abs = cumsum(value.abs)) %>%
  group_by(date, variable, iteration) %>%
  summarise(total.value.abs = sum(value.abs)) %>%
  group_by(variable,date) %>%
  summarise(total.value.abs.m = mean(total.value.abs)) 
total.noi.I = subset(SEIIR.table.noi.cum, date == max.date & variable == "newI")$total.value.abs.m
over65 = c("[65,70)", "[70,75)", "[75,80)", "[80,85)", "[85,90)")
between40and64 = c("[40,45)", "[45,50)", "[50,55)", "[55,60)", "[60,65)")
SEIIR.table.noi.death.over65 = subset(SEIIR.res.noi, age %in% over65) %>%
  group_by(date, variable, iteration) %>%
  summarise(total.value.abs = sum(value.abs)) %>%
  group_by(date, variable) %>%
  summarise(total.value.abs.m = mean(total.value.abs)) 
total.noi.death.over65 = subset(SEIIR.table.noi.death.over65, date == max.date & variable == "death")$total.value.abs.m
SEIIR.table.noi.death.40to64 = subset(SEIIR.res.noi, age %in% between40and64) %>%
  group_by(date, variable, iteration) %>%
  summarise(total.value.abs = sum(value.abs)) %>%
  group_by(date, variable) %>%
  summarise(total.value.abs.m = mean(total.value.abs)) 
total.noi.death.between40and64 = subset(SEIIR.table.noi.death.40to64, date == max.date & variable == "death")$total.value.abs.m
total.I = subset(SEIIR.res, variable == "newI") %>%
  group_by(variable,iteration, method,age,T) %>%
  mutate(value.abs = cumsum(value.abs)) %>%
  group_by(date, variable, iteration,method) %>%
  summarise(value.abs = sum(value.abs)) %>%
  subset(date == max.date) %>%
  mutate(value.abs.r = value.abs/total.noi.I) %>%
  group_by(variable, method) %>%
  summarise(value.abs.r.m = mean(value.abs.r),
            value.abs.r.l95 = quantile(value.abs.r, prob = 0.025),
            value.abs.r.u95 = quantile(value.abs.r, prob = 0.975))
total.death.over65 = subset(SEIIR.res, variable == "death"& age %in% over65) %>%
  group_by(variable,iteration, method,date) %>%
  summarise(value.abs = sum(value.abs)) %>%
  subset(date == max.date) %>%
  mutate(value.abs.r = value.abs/total.noi.death.over65) %>%
  group_by(variable, method) %>%
  summarise(value.abs.r.m = mean(value.abs.r),
            value.abs.r.l95 = quantile(value.abs.r, prob = 0.025),
            value.abs.r.u95 = quantile(value.abs.r, prob = 0.975))
total.death.over65$variable = "death over65"
total.death.between40and64 = subset(SEIIR.res, variable == "death"& age %in% between40and64) %>%
  group_by(variable,iteration, method,date) %>%
  summarise(value.abs = sum(value.abs)) %>%
  subset(date == max.date ) %>%
  mutate(value.abs.r = value.abs/total.noi.death.between40and64) %>%
  group_by(variable, method) %>%
  summarise(value.abs.r.m = mean(value.abs.r),
            value.abs.r.l95 = quantile(value.abs.r, prob = 0.025),
            value.abs.r.u95 = quantile(value.abs.r, prob = 0.975))
total.death.between40and64$variable = "death between40and64"
total.I.death = rbind(total.I, total.death.over65,total.death.between40and64)
save(total.I.death, file = file.path(indir, "results", "total.I.death.rda"))

labels_facet = c(`newI` = "Total cases",
                 `death over65` = "Total death over 65",
                 `death between40and64` = "Total death between 40 and 64")
p = ggplot(total.I.death, aes(x = method, col = method)) +
  geom_point(aes(y = value.abs.r.m)) +
  geom_errorbar(aes(ymin = value.abs.r.l95, ymax = value.abs.r.u95)) +
  facet_wrap(~variable, labeller = as_labeller(labels_facet)) +
  theme(axis.text.x=element_blank())+ 
  geom_hline(yintercept=1, size=1) +
  ylab("Value") +
  theme(legend.position="bottom") + 
  scale_color_brewer(palette="Dark2")
ggsave(p, w = 9, h = 6, device = "pdf", file  = file.path(indir, "results", "I.death.r.comp.pdf"))

# subset(total.I.death, method == "intervention 1, p = .6")
# 1 - subset(total.I.death, method == "intervention 1, p = .6")$value.abs.r.m
# 1 - subset(total.I.death, method == "intervention 1, p = .6")$value.abs.r.l95
# 1 - subset(total.I.death, method == "intervention 1, p = .6")$value.abs.r.u95
# 1 - subset(total.I.death, method == "intervention 4, p1 = .9, p2 = .6")$value.abs.r.m
# 1 - subset(total.I.death, method == "intervention 4, p1 = .9, p2 = .6")$value.abs.r.l95
# 1 - subset(total.I.death, method == "intervention 4, p1 = .9, p2 = .6")$value.abs.r.u95
