library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)

# device
if(1){
  indir = "~/git/social.mixing.ncov"
  outdir = "~/Box\ Sync/2020/social_mixing_ncov"
}

load(indir, "R", "stan_code.R")
load(file.path(outdir, "analyses", "contact.rda"))


## Number of participants by age
tmp = contact %>%
  group_by(age) %>%
  summarise(N = n()) 
ggplot(tmp, aes(x = age, y = N)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  labs(x = "Age", y = "Number of participants")
ggsave("age_participants.pdf", w = 8, h = 7, device = "pdf", path = file.path(outdir, "figures"))

## Contact matrix with naive estimates
tmp = contact %>%
  select(age_i, c.age.under5, c.age.6to19, c.age.20to64, c.age.over65) %>%
  group_by(age_i) %>%
  summarise(N = n(), 
            "under5" = sum(c.age.under5), 
            "6to19" = sum(c.age.6to19), 
            "20to64" = sum(c.age.20to64),
            "over65" = sum(c.age.over65)) %>%
  melt(id.vars = c("age_i", "N")) %>%
  rename(age_j = variable, contact = value) %>%
  mutate(naive_est = contact/N)
tmp$age_i = factor(tmp$age_i, levels = c("under5", "6to19", "20to64", "over65"))
tmp$age_j = factor(tmp$age_j, levels = c("under5", "6to19", "20to64", "over65"))

ggplot(tmp, aes(x = age_i, y = age_j)) +
  geom_raster(aes(fill = naive_est)) +   
  labs(x = "Age participant", y = "Age individual contacted")+ 
  guides(fill=guide_legend(title=expression(Y[ab]/N[a])))
ggsave("contact.matrix_naive.pdf", w = 8, h = 7, device = "pdf", path = file.path(outdir, "figures"))


