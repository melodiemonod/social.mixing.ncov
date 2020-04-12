# Load packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)

if(1){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
}

## load observed data
# cumulative death 
numb.deaths = read.csv(file.path(indir, "data", "time_series_covid19_deaths_global.csv")) %>%
  subset(Country.Region == "United Kingdom") %>%
  select(-c("Province.State","Lat","Long"))
numb.deaths = apply(numb.deaths[,-1], 2, sum)

# cumulative case
numb.cases = read.csv(file.path(indir, "data", "time_series_covid19_confirmed_global.csv")) %>%
  subset(Country.Region == "United Kingdom") %>%
  select(-c("Province.State","Lat","Long"))
numb.cases = apply(numb.cases[,-1], 2, sum)

# Plot cumulative death against time 
time.deaths = data.table(n.deaths = numb.deaths, date = seq(as.Date("2020-01-22"), by=1, len=length(numb.deaths)))
time.deaths = time.deaths[-(1:which(numb.cases > 0)[1]),]
time.1th.death = max(which(time.deaths$n.deaths <1)) + 1
ggplot(time.deaths, aes(x = date, y = log10(n.deaths))) +
  geom_point()

# fit data
load(file = file.path(indir, "results", "SEIIR.table.death.rda"))
p = ggplot(subset(SEIIR.table.death, variable == "death"), aes(x = date)) +
  geom_line(aes(y = log10(value.abs.m), col = method)) +
  geom_point(data = time.deaths, aes(y = log10(n.deaths))) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-05"))), linetype = "dashed", size = .5, alpha = .5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-16"))), linetype = "dashed", size = .5, alpha = .5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-19"))), linetype = "dashed", size = .5, alpha = .5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-23"))), linetype = "dashed", size = .5, alpha = .5) +
  ylab("Cumulative number of death (log 10)") +
  theme(text = element_text(size=15)) +
  theme_bw() +
  geom_ribbon(aes(ymin = log10(value.abs.l95),ymax = log10(value.abs.u95)), alpha = .5)
ggsave(p, w = 11, h = 7, file = file.path(indir, "results_simu", "fit.data.pdf"))


ggplot(subset(SEIIR.table.death, variable == "death" & t <50), aes(x = date)) +
  geom_line(aes(y = (value.abs.m), col = method)) +
  geom_point(data = time.deaths, aes(y = (n.deaths))) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-05"))), linetype = "dashed", size = .5, alpha = .5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-16"))), linetype = "dashed", size = .5, alpha = .5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-19"))), linetype = "dashed", size = .5, alpha = .5) +
  geom_vline(aes(xintercept = as.numeric(as.Date("2020-03-23"))), linetype = "dashed", size = .5, alpha = .5) +
  ylab("Cumulative number of death (log 10)") +
  theme(text = element_text(size=15)) +
  theme_bw() 
  #geom_ribbon(aes(ymin = (value.abs.l95),ymax = (value.abs.u95)), alpha = .5)
