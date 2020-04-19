indir = "~/git/social.mixing.ncov/contact.matrix"

source(file.path(indir, "R", "obtain.contact.estimates.R"))


countries.of.interest.polymod = c("BE", "DE", "IT", "GB")
    countries.of.interest.perm = c("Austria", "Switzerland", "Denmark", "Spain", "France", "Sweden")

`%notin%` <- Negate(`%in%`)

df = NULL
for(country in countries.of.interest.polymod){
  load(file.path(indir, "results", paste0("contact.tab.agg_", country, ".rda")))
  country.name = ifelse(country == "BE", "Belgium", 
         ifelse(country == "DE", "Germany",
                ifelse(country == "IT", "Italy",
                       ifelse(country == "GB", "United Kingdom", NA))))
  contact.tab.agg = contact.tab.agg %>%
    rename(pop = T) %>%
    subset(part.age.cat %notin% c("[80,85)",  "[85,90)",  "[90,95)",  "[95,100]") &
           cont.age.cat %notin% c("[80,85)",  "[85,90)",  "[90,95)",  "[95,100]")) %>%
    mutate(country = country.name) 
  contact.tab.agg$part.age.cat = plyr::mapvalues(contact.tab.agg$part.age.cat, from = "[75,80)", to = "[75,80]") 
  contact.tab.agg$cont.age.cat = plyr::mapvalues(contact.tab.agg$cont.age.cat, from = "[75,80)", to = "[75,80]")
  df = rbind(df,contact.tab.agg)
}
for(country in countries.of.interest.perm){
  load(file.path(indir, "results_from_perm", paste0("contact.matrix.agg_",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
  df = rbind(df,tmp)
}
df <- with(df, df[order(part.age.cat, cont.age.cat), ])

df = df %>%
  mutate(age_age_combination = paste0(part.age.cat, "-", cont.age.cat))

ggplot(df, aes(x = age_age_combination, y = log(c), col = country)) +
  geom_point() +
  scale_x_discrete(breaks = unique(df$age_age_combination)[rep(c(T, rep(F, length(unique(df$part.age.cat))-1)), length(unique(df$part.age.cat)))]) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


###

countries.of.interest = c("Austria", "Belgium", "Switzerland", "Germany", "Denmark", "Spain", "Finland", "France", "Italy", "Luxembourg", "Netherlands",
                          "Poland", "Sweden", "United Kingdom")
res = NULL
for(country in countries.of.interest){
  load(file.path(indir, "results_perm", paste0("contact.matrix.agg_",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
  res = rbind(res, contact.tab.agg )
}

area = data.table(country = countries.of.interest, 
                  km2 = c(83879, 30689, 41285, 357022, 42933, 505990, 338424, 551695, 301340, 2586.4, 41865, 312696, 450295, 242495))
res2 = res %>%
  group_by(country) %>%
  summarise(total.pop = sum(pop)/16,
            median.c = median(c),
            median.m = median(m)) %>%
  merge(area, by = "country") %>%
  mutate(pop.density.km2 = total.pop/km2)

ggplot(res2, aes(x = total.pop, y = log(median.c), col = country)) +
  geom_point() +
  geom_text(aes(label=country),hjust=0.1, vjust=0) + 
  coord_cartesian(xlim =c(0, (max(res2$total.pop) + 1e7)))
ggsave(file = "~/Downloads/c.pop.pdf", device = "pdf", w = 6, h = 5)

ggplot(res2, aes(x = total.pop, y = (median.m), col = country)) +
  geom_point() +
  geom_text(aes(label=country),hjust=0.1, vjust=0) + 
  coord_cartesian(xlim =c(0, (max(res2$total.pop) + 1e7)))
ggsave(file = "~/Downloads/m.pop.pdf", device = "pdf", w = 6, h = 5)

ggplot(res2, aes(x = pop.density.km2, y = log(median.c), col = country)) +
  geom_point() +
  geom_text(aes(label=country),hjust=0.1, vjust=0) + 
  coord_cartesian(xlim =c(0, (max(res2$pop.density.km2) + 70)))
ggsave(file = "~/Downloads/c.pop.dens.pdf", device = "pdf", w = 6, h = 5)

res3 = res %>%
  group_by(country, part.age.cat) %>%
  summarise(total.pop = sum(pop)/16,
            median.c = median(c),
            median.m = median(m))

ggplot(res3, aes(x = total.pop, y = (median.m), col = country)) +
  geom_point() +
  #geom_text(aes(label=country),hjust=0.1, vjust=0) + 
  facet_wrap(~part.age.cat)
ggsave(file = "~/Downloads/c.pop.part.age.cat.pdf", device = "pdf", w = 8, h = 6)

res4 = res %>%
  group_by(country, cont.age.cat) %>%
  summarise(total.pop = sum(pop)/16,
            median.c = median(c),
            median.m = median(m))

ggplot(res4, aes(x = total.pop, y = (median.m), col = country)) +
  geom_point() +
  #geom_text(aes(label=country),hjust=0.1, vjust=0) + 
  facet_wrap(~cont.age.cat)
ggsave(file = "~/Downloads/c.pop.cont.age.cat.pdf", device = "pdf", w = 8, h = 6)

####
df1 = NULL
for(country in countries.of.interest.polymod){
  load(file.path(indir, "results", paste0("polymod.tab.bin_", country, ".rda")))
  country.name = ifelse(country == "BE", "Belgium", 
                        ifelse(country == "DE", "Germany",
                               ifelse(country == "IT", "Italy",
                                      ifelse(country == "GB", "United Kingdom", NA))))
  polymod.tab = polymod.tab %>%
    rename(pop = T) %>%
    subset(part.age < 80 & cont.age < 80) %>%
    mutate(country = country.name, 
           method = "CAR") %>%
    select(part.age, cont.age, pop, m,c, country, method)
  df1 = rbind(df1,polymod.tab)
}
countries.also.in.perm = c("Belgium", "Germany", "Italy", "United Kingdom")
df2 = NULL
for(country in countries.also.in.perm){
  load(file.path(indir, "results_from_perm", paste0("contact.matrix",  str_replace_all(string=country, pattern=" ", repl=""), ".rda")))
  tmp$method = "Perm"
  df2 = rbind(df2,tmp)
}

test = merge(df1, df2, by = c("part.age", "cont.age", "pop", "country"))
test2 = test %>%
  group_by(country, part.age) %>%
  summarise(corr.c = cor(c.x, c.y))
ggplot(test2, aes(x = part.age, y = corr.c, col = country)) +
  geom_point()
ggsave(file = "~/Downloads/corr.c.part.age.pdf", device = "pdf", w = 6, h = 4)
test2 = test %>%
  group_by(country, cont.age) %>%
  summarise(corr.c = cor(c.x, c.y))
ggplot(test2, aes(x = cont.age, y = corr.c, col = country)) +
  geom_point()
ggsave(file = "~/Downloads/corr.c.cont.age.pdf", device = "pdf", w = 6, h = 4)
