#analysis contact simu 
### age analysis ###
### simulated data ###
age_obs = 1:17; age_group = 4
set.seed(789)
y = matrix(nrow = age_group, ncol = age_group, round(runif(age_group*age_group, 1, 20)))
corr = matrix(nrow = age_group, ncol = 2, c(1, 1, 2, 3, 4, 5), byrow = T)
data_list = list(n_group = age_group,
                 n_age = length(age_obs), 
                 y = y,
                 age_i = rep(age_obs, each = length(age_obs)), 
                 age_j = rep(age_obs, length(age_obs)),
                 corr = corr,
                 T = rep(100, length(age_obs)),
                 N = rep(50, length(age_obs)))

# exploratory analysis
ggplot(data = NULL, aes(x = rep(1:age_group, each = age_group), y = rep(1:age_group, age_group))) +
  geom_raster(aes(fill = as.vector(t(y))))

#analysis
iter = 1000; warmup = 500
fit = stan(model_code = model.stan.age.mixing, data = data_list, chains = 3, iter = iter, warmup = warmup)
posterior = extract(fit)
summary(fit)
mu = as.vector(sapply(1:data_list$n_group, function(x){
  as.vector(posterior$mu[,x,])
}))
df = data.table(mu = mu, agegroup_i = rep(rep(age_obs, each = length(age_obs)), each = (iter-warmup)*3), 
                agegroup_j = rep(rep(age_obs, length(age_obs)), each = (iter-warmup)*3)) %>%
  group_by(agegroup_i, agegroup_j) %>%
  summarise(mu.m = mean(mu), mu.l95 = quantile(mu, prob = .025), mu.u95 = quantile(mu, prob = .975))
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = mu.m))


## age and sex analysis##
iter = 100; warmup = 50
### simulated data ###
age_obs = 1:6;
n_age = length(age_obs)
set.seed(789)
df = data.table(cont.age = age_obs, 
                T = round(runif(n_age, 3,6)*1000000))
df2 = data.table(part.age = age_obs, 
                N = round(runif(n_age, 1,10)))
tmp1 = data.table(part.age = rep(age_obs, each = n_age), cont.age = rep(age_obs, n_age),
           y = round(runif(n_age*n_age, 0, 4)))
tmp = merge(merge(df, tmp1, by = "cont.age"), df2, by = "part.age") %>%
  mutate(U = T*N)

y = tmp %>%
  reshape2::dcast(part.age ~ cont.age, value.var = "y")

U = tmp %>%
  reshape2::dcast(part.age ~ cont.age, value.var = "U")

data_list = list(N_group = n_age,
                 n_group = n_age,
                 y = as.matrix(y[,-1]),
                 agegroup_i = rep(age_obs, each = n_age),
                 agegroup_j = rep(age_obs, n_age),
                 logU = log(as.matrix(U[,-1])))

iter = 1000; warmup = 100
fit = stan(model_code = model.contact.matrix, data = data_list, chains = 3, iter = iter, warmup = warmup)

N_age = n_age +2
age_obs = 1:8
df = data.table(cont.age = age_obs, 
                T = round(runif(N_age, 3,6)*1000000))
df2 = data.table(part.age = age_obs, 
                 N = c(round(runif(n_age, 1,10)), rep(0,2)))
tmp1 = data.table(part.age = rep(age_obs, each = N_age), cont.age = rep(age_obs, N_age),
                  y = round(runif(N_age*N_age, 0, 4)))
tmp = merge(merge(df, tmp1, by = "cont.age"), df2, by = "part.age") %>%
  mutate(U = ifelse(N == 0, 1, T*N))

y = tmp %>%
  reshape2::dcast(part.age ~ cont.age, value.var = "y")
y[c(7:8),] = 0
y[,c(8:9)] = 0

U = tmp %>%
  reshape2::dcast(part.age ~ cont.age, value.var = "U")

data_list = list(N_group = N_age,
                 n_group = n_age,
                 y = as.matrix(y[,-1]),
                 agegroup_i = rep(age_obs, each = N_age),
                 agegroup_j = rep(age_obs, N_age),
                 logU = log(as.matrix(U[,-1])))

iter = 1000; warmup = 100
fit = stan(model_code = model.contact.matrix, data = data_list, chains = 3, iter = iter, warmup = warmup)


# exploratory analysis
ggplot(tmp, aes(x = part.age, y = cont.age)) +
  geom_raster(aes(fill = y))

# analysis
posterior = extract(fit)
posterior3 = extract(fit)
# Divide c by 1e8 to go back to original scale (  deflates contact rate c)
c_est = as.vector(sapply(1:data_list$N_group, function(x){
  as.vector(posterior3$c[,x,])
}))
df2.1 = data.table(c = c_est, 
                   part.age = rep(rep(age_obs, each = N_age), each = (iter-warmup)*3),
                   cont.age = rep(rep(age_obs, N_age), each = (iter-warmup)*3)) %>%
  group_by(part.age, cont.age) %>%
  summarise(c.m = mean(c), c.l95 = quantile(c, prob = .025), c.u95 = quantile(c, prob = .975), 
            c.sd = sd(c))

ggplot(subset(df2.1, part.age != 7 & part.age != 8 & cont.age != 7 & cont.age != 8), aes(x = part.age, y = cont.age)) +
  geom_raster(aes(fill = c.m))

posterior2 = extract(fit2)
# Divide c by 1e8 to go back to original scale (  deflates contact rate c)
c_est = as.vector(sapply(1:data_list$n_group, function(x){
  as.vector(posterior2$c[,x,])
}))
df2.1 = data.table(c = c_est, 
                   part.age = rep(rep(age_obs, each = n_age), each = (iter-warmup)*3),
                   cont.age = rep(rep(age_obs, n_age), each = (iter-warmup)*3)) %>%
  group_by(part.age, cont.age) %>%
  summarise(c.m = mean(c), c.l95 = quantile(c, prob = .025), c.u95 = quantile(c, prob = .975), 
            c.sd = sd(c))

ggplot(df2.1, aes(x = part.age, y = cont.age)) +
  geom_raster(aes(fill = c.m))




posterior = extract(fit)
c = matrix(ncol = data_list$n_group, nrow = (iter-warmup)*3*data_list$n_age*data_list$n_age, 0)
for(k in 1:data_list$n_group){
  c[,k] = as.vector(sapply(1:data_list$n_age, function(x){
    as.vector(posterior$c[,x,,k])
  }))
}

df = data.table(mu = as.vector(mu), c = as.vector(c), 
                agegroup_i = rep(rep(rep(age_obs, each = length(age_obs)), each = (iter-warmup)*3), data_list$n_group), 
                agegroup_j = rep(rep(rep(age_obs, length(age_obs)), each = (iter-warmup)*3), data_list$n_group), 
                direction = rep(c("FF", "MM", "FM", "MF"), each = (iter-warmup)*3*data_list$n_age*data_list$n_age)) %>%
  group_by(agegroup_i, agegroup_j, direction) %>%
  summarise(mu.m = mean(mu), mu.l95 = quantile(mu, prob = .025), mu.u95 = quantile(mu, prob = .975), 
            c.m = mean(c), c.l95 = quantile(c, prob = .025), c.u95 = quantile(c, prob = .975))
df$direction <- factor(df$direction, levels = c("FM", "MM", "FF", "MF"))
ann_text <- data.frame(lab = c("Female", "Male", "Female", "Male"), 
                       direction = factor(c("FF", "FM", "FF", "MF"),levels = c("FM", "MM", "FF", "MF")), stringsAsFactors = F)
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = mu.m)) + 
  scale_fill_viridis_c() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title.x = element_text(vjust=-3), 
        axis.title.y = element_text(vjust=4), 
        plot.margin=unit(c(0,0,1,1),"lines")) + 
  geom_text(data = ann_text, aes(label = lab), x = c(rep(-.2,2),rep(mean(age_obs), 2)), y = c(rep(mean(age_obs), 2), rep(-.2,2)), 
            angle = c(90, 90, 360, 360)) +
  facet_wrap(~ direction, nrow = 2) +
  coord_cartesian(clip = "off") + 
  labs(fill = expression(E(Y[ij])), 
       x = "Mid-point Age-group i", 
       y = "Mid-point Age-group j") 
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = c.m)) + 
  scale_fill_viridis_c() + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.title.x = element_text(vjust=-3), 
        axis.title.y = element_text(vjust=4), 
        plot.margin=unit(c(0,0,1,1),"lines")) + 
  geom_text(data = ann_text, aes(label = lab), x = c(rep(-.2,2),rep(mean(age_obs), 2)), y = c(rep(mean(age_obs), 2), rep(-.2,2)), 
            angle = c(90, 90, 360, 360)) +
  facet_wrap(~ direction, nrow = 2) +
  coord_cartesian(clip = "off") + 
  labs(fill = expression(c[ij]), 
       x = "Mid-point Age-group i", 
       y = "Mid-point Age-group j")

