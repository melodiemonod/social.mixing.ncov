library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(gridExtra)

# device
if(1){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
}

## Demographic figures
pop.data <- read.csv(file = file.path(indir, "data", "demo_UK_2008.csv")) %>%
  group_by(age) %>%
  summarise(pop = sum(pop))
load(file =file.path(indir, "results", "posterior.draw.aggregate.rda")) 

pop.data.group <-contact.data.agg %>%
  select(cont.age.cat, w) %>%
  subset(!duplicated(cont.age.cat))


### DEMOGRAPHIC PARAMETERS ###
# age groups
age = sort(unique(posterior.draw.aggregate$part.age.cat))
n_age = length(age)

p1 = ggplot(pop.data, aes(x = age, y = pop/1000)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  xlab("Age") +
  ylab("Population (in thousands)")+
  theme(text = element_text(size=18))

p2 = ggplot(pop.data.group, aes(x = cont.age.cat, y = w/1000)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  xlab("Age group") +
  ylab("Population (in thousands)")+
  theme(text = element_text(size=18)) +
  scale_x_discrete(breaks=age[c(TRUE, FALSE)])

p = grid.arrange(p1,p2,nrow=2)
ggsave(p, w = 12, h = 7, device = "pdf", file = file.path(indir, "results_simu", "pop.age.pdf"))
