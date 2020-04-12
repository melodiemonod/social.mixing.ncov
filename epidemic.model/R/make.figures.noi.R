# Load packages
library(dplyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(gridExtra)

if(1){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
}

# HPC
if(0){
  indir = "/rds/general/user/mm3218/home/projects/social.mixing.ncov"
}

# no intervention: 
load(file.path(indir, "results", "SEIIR.res.noi.rda"))
load(file = file.path(indir, "results", "contact.data.agg.rda"))

pred.data = subset(contact.data.agg, part.age.cat != "[90,95)" & part.age.cat != "[95,100]" & 
                     cont.age.cat != "[90,95)" & cont.age.cat != "[95,100]")

# age groups
age = sort(unique(pred.data$part.age.cat))
n_age = length(age)

# number of individuals by age group
T = pred.data[!duplicated(pred.data$w), ][order(pred.data[!duplicated(pred.data$w), ]$cont.age.cat),"w"]

min.date = "2020-03-05"; max.date = "2021-03-05"
total.pop = sum(T)

SEIIR.table.cum = subset(SEIIR.res.noi, variable == "newI" | variable == "newIS" | variable == "newIA") %>%
  group_by(age, T, variable,iteration) %>%
  mutate(value = cumsum(value),
         value.abs = cumsum(value.abs)) %>%
  group_by(date,age, T, variable) %>%
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
SEIIR.table.new = SEIIR.res.noi %>%
  group_by(date, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975)) 
SEIIR.table.newIandcumI = rbind(SEIIR.table.cum, subset(SEIIR.table.new, variable == "newI" | variable == "newIS" | variable == "newIA")) 
total.pop.table = SEIIR.table.newIandcumI %>%
  group_by(variable,date) %>%
  summarise(total.value.m = sum(value.m),
            total.value.abs.m = sum(value.abs.m)) 
tmp = merge(SEIIR.table.newIandcumI, total.pop.table, by = c("variable", "date")) %>%
  mutate(value.r.abs.age = ifelse(value.abs.m == 0, 0, value.abs.m/total.value.abs.m),
         value.r.age = ifelse(value.m == 0, 0, value.m/total.value.m))
total.pop = sum(T)
### NEW CASES AND CUMULATIVE CASES ###
labels_facet = c(`newI` = "New cases", 
                 `total I` = "Cumulative cases")
p1= ggplot(subset(SEIIR.table.newIandcumI, variable == "newI"| variable == "total I"), aes(x = date, y = value.abs.m/1000000)) +
  geom_area(aes(fill = age)) +
  facet_wrap(~variable, nrow =3, scales = "free", labeller = as_labeller(labels_facet)) +
  theme_bw()+
  geom_line(aes(group = age), position = "stack", size = .45, color = "azure4") +
  scale_x_date(date_breaks = "2 months" , date_labels = "%b-%y") + 
  theme(legend.position = "none")+
  scale_y_continuous("Value (in millions)", sec.axis = sec_axis(~ . *1000000/ total.pop, name = "Proportion relative to total population"))
p2 = ggplot(subset(tmp, variable == "newI"| variable == "total I"), aes(x = date, y = value.r.abs.age)) +
  geom_area(aes(fill = age)) +
  facet_wrap(~variable, nrow =3, scales = "free", labeller = as_labeller(labels_facet)) +
  theme_bw()+
  geom_line(aes(group = age), position = "stack", size = .45, color = "azure4") +
  scale_x_date(date_breaks = "2 months" , date_labels = "%b-%y") +
  ylab("Proportion")
pp=grid.arrange(p1,p2, nrow = 1, widths = c(1,1.07))
ggsave(pp, w =9, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newI.cumI.age.noi.prop.pdf"))

## number
df = subset(SEIIR.table.newIandcumI, variable == "total I" & date == max.date)
sum(df$value.abs.m)/total.pop; sum(df$value.abs.l95)/total.pop; sum(df$value.abs.u95)/total.pop 

### PROPORTION OF ASYMPTOMATIC AND SYMPTOMATIC ###
pp= ggplot(subset(SEIIR.table.newIandcumI, variable == "total IA" & date==max.date| variable == "total IS"& date == max.date), aes(x = age,y = value.m, fill = variable)) +
  geom_bar(stat = "identity") +
  theme(text = element_text(size=18)) +
  theme_bw()+
  ylab("Proportion after one year") +
  xlab("age group")+
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])
ggsave(pp, w = 8, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "propIA.IS.noi.pdf"))

## number
df = subset(SEIIR.table.newIandcumI, variable == "total IA" &  date == max.date)
df <- with(df, df[order(age), ])
sum(df$value.m[0:4])/4; sum(df$value.l95[0:4])/4; sum(df$value.u95[0:4])/4

### NEW DEATHS AND CUMULATIVE DEATHS BY AGE ###
SEIIR.table.new$variable = ordered(SEIIR.table.new$variable, 
                                       levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
total.pop.death = SEIIR.table.new %>%
  group_by(variable,date) %>%
  summarise(total.value.m = sum(value.m),
            total.value.l95 = sum(value.l95),
            total.value.u95 = sum(value.u95),
            total.value.abs.m = sum(value.abs.m),
            total.value.abs.l95 = sum(value.abs.l95),
            total.value.abs.u95 = sum(value.abs.u95)) 
tmp2 = merge(SEIIR.table.new, total.pop.death, by = c("variable", "date")) %>%
  mutate(value.r.abs.age = ifelse(value.abs.m == 0, 0, value.abs.m/total.value.abs.m),
         value.r.abs.age.l95 = ifelse(value.abs.l95 == 0, 0, value.abs.l95/total.value.abs.m),
         value.r.abs.age.u95 = ifelse(value.abs.u95 == 0, 0, value.abs.u95/total.value.abs.m),
         value.r.age = ifelse(value.m == 0, 0, value.m/total.value.m),
         value.r.age.l95 = ifelse(value.l95 == 0, 0, value.l95/total.value.m),
         value.r.age.u95 = ifelse(value.u95 == 0, 0, value.u95/total.value.m))
labels_facet = c(`newdeath` = "New deaths", 
                 `death` = "Cumulative deaths")
p1 =ggplot(data = NULL, aes(x = date)) +
  geom_area(data = subset(SEIIR.table.new, variable == "newdeath"), aes(fill = age, y = value.abs.m/1000)) +
  geom_area(data = subset(SEIIR.table.new, variable == "death"), aes(fill = age, y = (value.abs.m)/1000)) +
  facet_wrap(~variable, nrow =3, scales = "free", labeller = as_labeller(labels_facet)) +
  theme_bw()+
  scale_x_date(date_breaks = "2 months" , date_labels = "%b-%y") +
  scale_y_continuous("Value (in thousands)", sec.axis = sec_axis(~ . *1000 / total.pop, name = "Proportion relative to total population")) + 
  geom_line(data = subset(SEIIR.table.new, variable == "newdeath"), aes(group = age, y = value.abs.m/1000), position = "stack", size = .45, color = "azure4")+
  geom_line(data = subset(SEIIR.table.new, variable == "death"), aes(group = age, y = value.abs.m/1000), position = "stack", size = .45, color = "azure4") + 
  theme(legend.position = "none")
p2 = ggplot(subset(tmp2, variable == "newdeath"| variable == "death"), aes(x = date, y = value.r.abs.age)) +
  geom_area(aes(fill = age)) +
  facet_wrap(~variable, nrow =3, scales = "free", labeller = as_labeller(labels_facet)) +
  theme_bw()+
  geom_line(aes(group = age), position = "stack", size = .45, color = "azure4") +
  scale_x_date(date_breaks = "2 months" , date_labels = "%b-%y") +
  ylab("Proportion")
pp= grid.arrange(p1,p2, nrow = 1, widths = c(1,1.1))
ggsave(pp, w = 9, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.cumdeath.age.noi.prop.pdf"))

# number
df1 = subset(SEIIR.res.noi) %>%
  group_by(variable,iteration,date) %>%
  summarise(total.value.abs.m = sum(value.abs))
df2=merge(df1, SEIIR.res.noi, c("variable", "iteration", "date"))
df3 =df2 %>%
  group_by(age, T, variable,date) %>%
  mutate(value.abs.p = ifelse(value.abs ==0, value.abs, value.abs/total.value.abs.m)) %>%
  summarise(value.abs.p.m = mean(value.abs.p), 
            value.abs.p.l95 = quantile(value.abs.p, prob = 0.025), 
            value.abs.p.u95 = quantile(value.abs.p, prob = 0.975))
df = subset(df3, variable == "death" & date == max.date)

sum(df$value.abs.p.m[14:18]); sum(df$value.abs.p.l95[14:18]); sum(df$value.abs.p.u95[14:18])
sum(df$value.abs.p.m[9:13]); sum(df$value.abs.p.l95[9:13]); sum(df$value.abs.p.u95[9:13]) 

df = subset(SEIIR.table.new, variable == "death"&date==max.date)
sum(df$value.abs.m); sum(df$value.abs.l95); sum(df$value.abs.u95)



## qA and qS plots
set.seed(5678)
theta.draws.noi = lapply(1:1000, function(x){
  draw.theta(intervention = 0)
})
qA = c(); qS = c()
for(i in 1:1000){
  qA[i] = theta.draws.noi[[i]]$qA
  qS[i] = theta.draws.noi[[i]]$qA*1.5
}
df = data.table(parameter = c(rep("qA", 1000),rep("qS", 1000)), value = c(qA,qS)) 
mu <- df %>% 
  group_by(parameter) %>%
  summarise(value.m=mean(value)) 
ggplot(df, aes(x=value, col = parameter))+
  geom_density()+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  geom_point(data=mu, aes(x=value.m, y = c(0,0), color=parameter),shape=4,size=2) +
  theme_bw()
ggsave(file = file.path(indir, "results_simu", "qA.qS.pdf"), device = "pdf", w = 6, h = 4)

