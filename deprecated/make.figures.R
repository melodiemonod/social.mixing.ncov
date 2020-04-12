
load(file.path(indir, "results_simu", "SEIIR.res.noi.rda"))
SEIIR.table.noi = SEIIR.res.noi %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) 
SEIIR.table.noi$intervention = "no intervention"
SEIIR.table2.noi = SEIIR.res.noi %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
SEIIR.table3.noi = SEIIR.res.noi %>%
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
SEIIR.table3.noi$intervention = "no intervention"
load(file.path(indir, "results_simu", "SEIIR.res.i1.p50.rda"))
SEIIR.table.i1.p50 = SEIIR.res.i1.p50 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) 
SEIIR.table.i1.p50$intervention = "intervention 1, p=.5"
SEIIR.table2.i1.p50 = SEIIR.res.i1.p50 %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
SEIIR.table3.i1.p50 = SEIIR.res.i1.p50 %>%
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
SEIIR.table3.i1.p50$intervention = "intervention 1, p=.5"
########
load(file.path(indir, "results_simu", "SEIIR.res.i1.p20.rda"))
SEIIR.table.i1.p20 = SEIIR.res.i1.p20 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) 
SEIIR.table.i1.p20$intervention = "intervention 1, p=.2"
SEIIR.table2.i1.p20 = SEIIR.res.i1.p20 %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
SEIIR.table3.i1.p20 = SEIIR.res.i1.p20 %>%
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
SEIIR.table3.i1.p20$intervention = "intervention 1, p=.2"
#####
load(file.path(indir, "results_simu", "SEIIR.res.i1.p10.rda"))
SEIIR.table.i1.p10 = SEIIR.res.i1.p10 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) 
SEIIR.table.i1.p10$intervention = "intervention 1, p=.1"
SEIIR.table2.i1.p10 = SEIIR.res.i1.p10 %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
SEIIR.table3.i1.p10 = SEIIR.res.i1.p10 %>%
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
SEIIR.table3.i1.p10$intervention = "intervention 1, p=.1"
#####
load(file.path(indir, "results_simu", "SEIIR.res.i1.p01.rda"))
SEIIR.table.i1.p01 = SEIIR.res.i1.p01 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) 
SEIIR.table.i1.p01$intervention = "intervention 1, p=.01"
SEIIR.table2.i1.p01 = SEIIR.res.i1.p01 %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
SEIIR.table3.i1.p01 = SEIIR.res.i1.p01 %>%
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
SEIIR.table3.i1.p01$intervention = "intervention 1, p=.01"
######
load(file.path(indir, "results_simu", "SEIIR.res.i2.rda"))
SEIIR.table.i2 = SEIIR.res.i2 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) 
SEIIR.table.i2$intervention = "intervention 2"
SEIIR.table2.i2 = SEIIR.res.i2 %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
SEIIR.table3.i2 = SEIIR.res.i2 %>%
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
SEIIR.table3.i2$intervention = "intervention 2"
#####
load(file.path(indir, "results_simu", "SEIIR.res.i3.rda"))
SEIIR.table.i3 = SEIIR.res.i3 %>%
  group_by(t, variable, iteration) %>%
  summarise(value.abs.2 = sum(value.abs)) %>%
  group_by(t, variable) %>%
  summarise(value.abs.m = mean(value.abs.2), 
            value.abs.l95 = quantile(value.abs.2, prob = 0.025), 
            value.abs.u95 = quantile(value.abs.2, prob = 0.975)) 
SEIIR.table.i3$intervention = "intervention 3"
SEIIR.table2.i3 = SEIIR.res.i3 %>%
  group_by(t, age, T, variable) %>%
  summarise(value.m = mean(value), 
            value.l95 = quantile(value, prob = 0.025), 
            value.u95 = quantile(value, prob = 0.975),
            value.abs.m = mean(value.abs), 
            value.abs.l95 = quantile(value.abs, prob = 0.025), 
            value.abs.u95 = quantile(value.abs, prob = 0.975))
SEIIR.table3.i3 = SEIIR.res.i3 %>%
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
SEIIR.table3.i3$intervention = "intervention 3"

SEIIR.table.all = rbind(SEIIR.table.noi, SEIIR.table.i1.p50,SEIIR.table.i1.p20,SEIIR.table.i1.p10,SEIIR.table.i1.p01,SEIIR.table.i2,SEIIR.table.i3)
SEIIR.table3.all = rbind(SEIIR.table3.noi, SEIIR.table3.i1.p50,SEIIR.table3.i1.p20,SEIIR.table3.i1.p10,SEIIR.table3.i1.p01,SEIIR.table3.i2,SEIIR.table3.i3)

# no intervention: 
#new IS, new IA, new I by age
SEIIR.table2.noi$variable = ordered(SEIIR.table2.noi$variable, 
                                   levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
p = ggplot(subset(SEIIR.table2.noi, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_wrap(~variable+age, nrow =3) +
  ggtitle("newI - no intervention") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p, w = 16, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.age.noi.pdf"))
# new death and total death by age
p2 = ggplot(subset(SEIIR.table2.noi, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_grid(rows = vars(variable), cols = vars(age), 
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - no intervention") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p2, w = 16, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.age.noi.pdf"))

ggplot(subset(SEIIR.table.noi, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95), alpha = 0.55) +
  facet_grid(rows = vars(variable),  
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - no intervention") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0))

# intervention 1 p50: 
# new IS, new IA, new I by age
SEIIR.table2.i1.p50$variable = ordered(SEIIR.table2.i1.p50$variable, 
                                    levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
p = ggplot(subset(SEIIR.table2.i1.p50, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_wrap(~variable+age, nrow =3) +
  ggtitle("newI - intervention 1, p=.5") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p, w = 16, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.age.i1.p50.pdf"))
# new death and total death by age
p2 = ggplot(subset(SEIIR.table2.i1.p50, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_grid(rows = vars(variable), cols = vars(age), 
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - intervention 1, p=.5") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p2, w = 16, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.age.i1.p50.pdf"))

# intervention 1 p20: 
# new IS, new IA, new I by age
SEIIR.table2.i1.p20$variable = ordered(SEIIR.table2.i1.p20$variable, 
                                       levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
p = ggplot(subset(SEIIR.table2.i1.p20, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_wrap(~variable+age, nrow =3) +
  ggtitle("newI - intervention 1, p=.2") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p, w = 16, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.age.i1.p20.pdf"))
# new death and total death by age
p2 = ggplot(subset(SEIIR.table2.i1.p20, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_grid(rows = vars(variable), cols = vars(age), 
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - intervention 1, p=.2") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p2, w = 16, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.age.i1.p20.pdf"))

# intervention 1 p10: 
# new IS, new IA, new I by age
SEIIR.table2.i1.p10$variable = ordered(SEIIR.table2.i1.p10$variable, 
                                       levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
p = ggplot(subset(SEIIR.table2.i1.p10, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_wrap(~variable+age, nrow =3) +
  ggtitle("newI - intervention 1, p=.1") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p, w = 16, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.age.i1.p10.pdf"))
# new death and total death by age
p2 = ggplot(subset(SEIIR.table2.i1.p10, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_grid(rows = vars(variable), cols = vars(age), 
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - intervention 1, p=.1") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p2, w = 16, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.age.i1.p10.pdf"))

# intervention 1 p01: 
# new IS, new IA, new I by age
SEIIR.table2.i1.p01$variable = ordered(SEIIR.table2.i1.p01$variable, 
                                       levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
p = ggplot(subset(SEIIR.table2.i1.p01, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_wrap(~variable+age, nrow =3) +
  ggtitle("newI - intervention 1, p=.01") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p, w = 16, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.age.i1.p01.pdf"))
# new death and total death by age
p2 = ggplot(subset(SEIIR.table2.i1.p01, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_grid(rows = vars(variable), cols = vars(age), 
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - intervention 1, p=.01") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p2, w = 16, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.age.i1.p01.pdf"))

# intervention 2: 
# new IS, new IA, new I by age
SEIIR.table2.i2$variable = ordered(SEIIR.table2.i2$variable, 
                                       levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
p = ggplot(subset(SEIIR.table2.i2, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_wrap(~variable+age, nrow =3) +
  ggtitle("newI - intervention 2") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p, w = 16, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.age.i2.pdf"))
# new death and total death by age
p2 = ggplot(subset(SEIIR.table2.i2, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_grid(rows = vars(variable), cols = vars(age), 
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - intervention 2") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p2, w = 16, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.age.i2.pdf"))

# intervention 3: 
# new IS, new IA, new I by age
SEIIR.table2.i3$variable = ordered(SEIIR.table2.i3$variable, 
                                   levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
p = ggplot(subset(SEIIR.table2.i3, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_wrap(~variable+age, nrow =3) +
  ggtitle("newI - intervention 3") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p, w = 16, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.age.i3.pdf"))
# new death and total death by age
p2 = ggplot(subset(SEIIR.table2.i3, variable == "newdeath" | variable == "death"), aes(x = t)) +
  geom_line(aes(y = value.abs.m,col = age)) +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= age), alpha = 0.55) +
  facet_grid(rows = vars(variable), cols = vars(age), 
             scales="free_y", switch = 'y') +
  ggtitle("newdeath and accumulated death - intervention 3") +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(mean(1:t), digits = 0)) 
ggsave(p2, w = 16, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "newdeath.age.i3.pdf"))


### ALL INTERVENTION
SEIIR.table.all$intervention = ordered(SEIIR.table.all$intervention, 
                                       levels = c("no intervention", "intervention 1, p=.5", "intervention 1, p=.2","intervention 1, p=.1",
                                                  "intervention 1, p=.01","intervention 2","intervention 3"))
SEIIR.table3.all$intervention = ordered(SEIIR.table3.all$intervention, 
                                       levels = c("no intervention", "intervention 1, p=.5", "intervention 1, p=.2","intervention 1, p=.1",
                                                  "intervention 1, p=.01","intervention 2","intervention 3"))
SEIIR.table.all$variable = ordered(SEIIR.table.all$variable, 
                                        levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))
SEIIR.table3.all$variable = ordered(SEIIR.table3.all$variable, 
                                    levels = c("newIS", "newIA", "newI", "newdeath", "S", "E", "IA", "IS", "R", "death"))

# S, E, IS, IA, R by intervention
p = ggplot(subset(SEIIR.table.all, variable == "S" | variable == "IS" | variable == "IA" | variable == "R" | variable == "E"), aes(x = t)) +
  geom_line(aes(y = value.abs.m, col = variable)) +
  ggtitle("S, E, IA, IS, R") +
  geom_ribbon(aes(ymin=value.abs.l95, ymax =value.abs.u95, fill= variable), alpha = 0.55) +
  theme(text = element_text(size=15)) +
  scale_x_continuous(breaks = round(seq(1, t, length.out = 3), digits = 0)) +
  facet_wrap(~intervention, nrow =1)
ggsave(p, w = 18, h = 6, device = "pdf", file  = file.path(indir, "results_simu", "SEIIR.intervention.pdf"))

# new IS, IA, I by intervention
p2 = ggplot(subset(SEIIR.table.all, variable == "newI"| variable == "newIS"| variable == "newIA"), aes(x = t)) +
  geom_ribbon(aes(ymin = value.abs.l95, ymax = value.abs.u95), alpha = 0.5) +
  geom_line(aes(y = value.abs.m)) +
  facet_wrap(~variable+intervention, nrow =3) +
  ggtitle("newIS + newIA = newI by intervention") +
  theme(text = element_text(size=15)) 
ggsave(p2, w = 12, h = 10, device = "pdf", file  = file.path(indir, "results_simu", "newI.intervention.pdf"))

p.intervention.death= ggplot(subset(SEIIR.table.all, variable == "death"|variable == "newdeath"), aes(x = t)) +
  geom_line(aes(y = value.abs.m)) +
  geom_ribbon(aes(ymin = value.abs.l95, ymax = value.abs.u95), alpha = 0.5) +
  facet_grid(rows = vars(variable), cols = vars(intervention), 
             scales="free_y", switch = 'y') +
  ggtitle("New death and accumulated death by intervention") +
  theme(text = element_text(size=15)) +
  geom_blank()
ggsave(p.intervention.death, w = 16, h = 9, device = "pdf", file  = file.path(indir, "results_simu", "death.time.intervention.pdf"))

p.intervention.death2 = ggplot(subset(SEIIR.table3.all, variable == "newdeath"), aes(x = age)) +
  geom_point(aes(y = total.abs.m)) +
  geom_errorbar(aes(ymin = total.abs.l95, ymax = total.abs.u95)) +
  facet_wrap(~intervention, nrow =1) +
  ggtitle("Total death by age and intervention") +
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_discrete(breaks=age[c(TRUE, FALSE)])
ggsave(p.intervention.death2, w = 15, h = 7, device = "pdf", file  = file.path(indir, "results_simu", "totaldeath.age.intervention.pdf"))
