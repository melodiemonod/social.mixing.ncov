# IN THIS FILE I USE THE STAN MODEL TO OBTAIN CONTACT MATRIX WITH AGE MIXING ON SIMULATED DATA

model.stan.age.mixing	<-  "
    data{
      int<lower=1> n_group;  // number of age group
      int<lower=1> n_age;  // number of age
      int y[n_group, n_group];  // the outcome variable
      real age_i[n_age*n_age]; // age group of i
      real age_j[n_age*n_age]; // age group of j
      int corr[n_group,2]; // correspondance between age and age group
      int T[n_age]; // total number of participants by age group 
      int N[n_age]; // total number of individuals by age group 
     }
    
    transformed data {
      int n_cat = n_age*n_age; // number of categories by age
      row_vector[2] age[n_cat];
      row_vector[2] age_sym[n_cat];
      int U[n_age, n_age];
      for(i in 1:n_cat){
        age[i, 1] = age_i[i];
        age[i, 2] = age_j[i];
        age_sym[i, 1] = age_j[i];
        age_sym[i, 2] = age_i[i];
      }
      for(i in 1:n_age){
        for(j in 1:n_age){
          U[i,j] = T[i] * N[j];
        }
      }
    }
    
    parameters{
      // Global dispersion parameter 
      real ltheta[n_group];
  
      // Global intercept
      real beta;

      // Parameters for covariance matrix
      real<lower=0> rho;
      real<lower=0> alpha;
      vector[n_cat] x_delta;
    }
  
    transformed parameters{
      matrix[n_cat, n_cat] cov =  cov_exp_quad(age, age, alpha, rho) + cov_exp_quad(age_sym, age, alpha, rho)  + diag_matrix(rep_vector(1e-10, n_cat));
      matrix[n_cat, n_cat] L_cov = cholesky_decompose(cov);
      vector[n_cat] x_vec = L_cov * x_delta;
      matrix[n_age, n_age] x = to_matrix(x_vec, n_age, n_age, 0);
      real c[n_age, n_age];
      matrix[n_age, n_age] mu_age;
      real mu_group[n_group, n_group];
      for(i in 1:n_age){
        for(j in 1:n_age){
          c[i,j] = exp(beta + x[i,j]);
          mu_age[i,j] = U[i,j] + c[i,j];
        }
      }
      for(a in 1:n_group){
        for(b in 1:n_group){
          mu_group[a,b] = sum(mu_age[corr[a,1]:corr[a,2], corr[b,1]:corr[b,2]]);
        }
      }
    }

    model{
      beta ~ normal(0,0.001);
      to_vector(ltheta) ~ normal(0,0.001);

      rho ~ inv_gamma(5, 5);
      alpha ~ std_normal();
      x_delta ~ std_normal();
          
      for(a in 1:n_group){
        for(b in 1:n_group){
          y[a,b] ~ normal(mu_group[a,b], exp(ltheta[a]));
        }
      }
    }
  "

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
set.seed(789)
n_age = 5
age_obs = 1:n_age
n_group = 3
y = matrix(nrow = n_group, ncol = n_group, round(runif(n_group*n_group, 1, 20)))
corr = matrix(nrow = n_group, ncol = 2, c(1, 2, 3, 4, 5, 5), byrow = T)
data_list = list(n_group = n_group,
                 n_age = n_age, 
                 y = y,
                 age_i = rep(age_obs, each = n_age), 
                 age_j = rep(age_obs, n_age),
                 corr = corr, 
                 T = rep(1, n_age),
                 N = rep(1, n_age))
# exploratory analysis
ggplot(data = NULL, aes(x = rep(1:n_group, each = n_group), y = rep(1:n_group, n_group))) +
  geom_raster(aes(fill = as.vector(t(data_list$y))))
#analysis
fit = stan(model_code = model.stan.age.mixing, data = data_list, chains = 3, iter = iter, warmup = warmup)
posterior = extract(fit)
summary(fit)
mu = as.vector(sapply(1:data_list$n_group, function(x){
  as.vector(posterior$mu_group[,x,])
}))
df = data.table(mu = mu, agegroup_i = rep(rep(1:n_group, each = n_group), each = (iter-warmup)*3), 
                agegroup_j = rep(rep(1:n_group, n_group), each = (iter-warmup)*3)) %>%
  group_by(agegroup_i, agegroup_j) %>%
  summarise(mu.m = mean(mu), mu.l95 = quantile(mu, prob = .025), mu.u95 = quantile(mu, prob = .975))
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = mu.m)) +
  scale_fill_gradientn(colours = myPalette(100), limits=c(1, 20))
mu = as.vector(sapply(1:data_list$n_age, function(x){
  as.vector(posterior$mu_age[,x,])
}))
df = data.table(mu = mu, agegroup_i = rep(data_list$age_i, each = (iter-warmup)*3), 
                agegroup_j = rep(data_list$age_j, each = (iter-warmup)*3)) %>%
  group_by(agegroup_i, agegroup_j) %>%
  summarise(mu.m = mean(mu), mu.l95 = quantile(mu, prob = .025), mu.u95 = quantile(mu, prob = .975))
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = mu.m))+
  scale_fill_gradientn(colours = myPalette(100), limits=c(1, 20))

