model.overall.agegroup.deterministic.h1n1	<-  "
functions {
  real[,] SIR(
  int t,           // time
  real[] y,     // system state {susceptible,infected,recovered}
  matrix beta,     // parameters
  row_vector gamma,
  int n_group) {

  real total; 
  real y_hat[t+1, n_group*3];
  
  for(i in 1:n_group){
    y_hat[1,1+3*(i-1)] = y[1+3*(i-1)];
    y_hat[1,2+3*(i-1)] = y[2+3*(i-1)];
    y_hat[1,3+3*(i-1)] = y[3+3*(i-1)];
  }
  
  for(time in 2:t){
    for(i in 1:n_group){
    total = 0.0;
    for(j in 1:n_group){
      total = total + beta[i,j] * y_hat[time,2+3*(i-1)];
    }
    y_hat[time+1,1+3*(i-1)] = y_hat[time,1+3*(i-1)] - total * y_hat[time,1+3*(i-1)];
    y_hat[time+1,2+3*(i-1)] = y_hat[time,2+3*(i-1)] + total * y_hat[time,1+3*(i-1)] - gamma[i] * y_hat[time,2+3*(i-1)];
    y_hat[time+1,3+3*(i-1)] = y_hat[time,3+3*(i-1)] + gamma[i] * y_hat[time,2+3*(i-1)];
    }
  }
  return y_hat[2:(t+1), 1:(n_group*3)];
  }
}

data{
  int<lower = 1> n_time1;      // number of days observed
  int<lower = 1> n_time2;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int Y1[n_group, n_group];    // number of contact by age group
  int Y2[n_group, n_group];    // number of contact by age group
  int Z[n_time1+n_time2, n_group];       // number of infected individuals each day by age
  real agegroup_i[n_group*n_group];    // age group of i
  real agegroup_j[n_group*n_group];    // age group of j
  int T[n_group];               // total number of participants by age group 
  int N[n_group,2];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time1+n_time2];            // time points observed
}
  
transformed data {
  real x_r[0];
  int x_i[1];
  int n_cat = n_group*n_group; // number of categories by age
  row_vector[2] agegroup[n_cat];
  row_vector[2] agegroup_sym[n_cat];
  real logU1[n_group, n_group];
  real logU2[n_group, n_group];
  x_i[1] = n_group;
  for(i in 1:n_cat){
    agegroup[i, 1] = agegroup_i[i];
    agegroup[i, 2] = agegroup_j[i];
    agegroup_sym[i, 1] = agegroup_j[i];
    agegroup_sym[i, 2] = agegroup_i[i];
  }
  for(i in 1:n_group){
    for(j in 1:n_group){
      logU1[i,j] = log(T[i]) + log(N[j,1]);
      logU2[i,j] = log(T[i]) + log(N[j,2]);
    }
  }
}

parameters{
  // Transmission model parameters
  real<lower=0> theta[2];           // global dispersion parameter 
  real alpha[2];                    // global intercept
  real<lower=0> rho[2];             // length-scale
  real<lower=0> sigma[2];           // amplitude
  real x_delta[n_cat,2];         // gp function

  // Epidemic model parameters
  real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
  row_vector<lower = 0>[n_group] gamma;   // rate of recovery
  real<lower = 0, upper = 1> s0[n_group]; // initial fraction of susceptible 
  real<lower = 0, upper = 1> s_vac; // fraction of susceptible in vacation 
}
  
transformed parameters{
  // Transmission model transformed parameters first part
  matrix[n_cat, n_cat] cov1 =  cov_exp_quad(agegroup, agegroup, sigma[1], rho[1]) + cov_exp_quad(agegroup_sym, agegroup, sigma[1], rho[1])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov1 = cholesky_decompose(cov1);
  vector[n_cat] x_vec1 = L_cov1 * to_vector(x_delta[,1]);
  matrix[n_group, n_group] f1 = to_matrix(x_vec1, n_group, n_group, 0);
  matrix[n_group, n_group] c1;
  matrix[n_group, n_group] mu_age1;
  
  // Transmission model transformed parameters second part
  matrix[n_cat, n_cat] cov2 =  cov_exp_quad(agegroup, agegroup, sigma[2], rho[2]) + cov_exp_quad(agegroup_sym, agegroup, sigma[2], rho[2])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov2 = cholesky_decompose(cov2);
  vector[n_cat] x_vec2 = L_cov2 * to_vector(x_delta[,2]);
  matrix[n_group, n_group] f2 = to_matrix(x_vec2, n_group, n_group, 0);
  matrix[n_group, n_group] c2;
  matrix[n_group, n_group] mu_age2;

  // Epidemic model transformed parameters
  real y_hat[(n_time1+n_time2), n_difeq * n_group]; // solution from the ODE solver
  real y_init[n_difeq * n_group];    // initial conditions for both fractions of S, I and R
  real y_term[n_difeq * n_group];    // conditions for both fractions of S, I and R after holiday 
  matrix[n_group, n_group] beta1;
  matrix[n_group, n_group] beta2;
  real beta_gamma1[n_group + n_group*n_group]; 
  real beta_gamma2[n_group + n_group*n_group];
  
  for(i in 1:n_group){
    for(j in 1:n_group){
      c1[i,j] = exp(alpha[1] + f1[i,j]);
      mu_age1[i,j] = logU1[i,j] + c1[i,j];
      c2[i,j] = exp(alpha[2] + f2[i,j]);
      mu_age2[i,j] = logU2[i,j] + c2[i,j];
    }
    y_init[1+3*(i-1)] = s0[i]*(1-s_vac);
    y_init[2+3*(i-1)] = 1 - s0[i];
    y_init[3+3*(i-1)] = 0;
  }

  beta1 = c1 * diag_matrix(to_vector(T)) * q;
  beta2 = c2 * diag_matrix(to_vector(T)) * q;
  beta_gamma1 = to_array_1d(append_row(gamma, beta1)');
  beta_gamma2 = to_array_1d(append_row(gamma, beta2)');
  y_hat[1:n_time1,] = SIR(n_time1, y_init, beta1, gamma, n_group);
  
  for(i in 1:n_group){
    y_term[1+3*(i-1)] = y_hat[n_time1,1+3*(i-1)] + s0[i]*s_vac;
    y_term[2+3*(i-1)] = y_hat[n_time1,2+3*(i-1)];
    y_term[3+3*(i-1)] = y_hat[n_time1,3+3*(i-1)];
  }
  y_hat[(n_time1+1):(n_time1+n_time2),] = SIR(n_time2, y_term, beta2, gamma, n_group);
}
  
model{
  real lambda[(n_time1+n_time2), n_group];      //poisson parameter
  
  //priors
  // Transmission model priors
  alpha ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  rho ~ inv_gamma(5, 5);
  sigma ~ std_normal();
  to_array_1d(x_delta) ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(10, 0.1);
  s_vac ~ beta(2, 5);
  
  // Likelihood
  // Transmission model likelihood
  for(i in 1:n_group){
    for(j in 1:n_group){
      Y1[i,j] ~ normal(mu_age1[i,j], theta);
      Y2[i,j] ~ normal(mu_age2[i,j], theta);
    }
  }
  //Epidemic model likelihood
  for(t in 1:(n_time1+n_time2)){
    for(i in 1:n_group){
      lambda[t,i] = y_hat[t,2+3*(i-1)]*T[i];
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }
}

//generated quantities {
//  real R_0[n_group, n_group];  // Basic reproduction number
//  for(i in 1:n_group){
//    for(j in 1:n_group){
//      R_0[i,j] = beta[i,j]/gamma[i];   
//    }
//  }
//}
"



model.overall.agegroup.deterministic.h1n1	<-  "
functions {
  real[,] SIR(
  int t,           // time
  real[] y,     // system state {susceptible,infected,recovered}
  matrix beta,     // parameters
  row_vector gamma,
  int n_group, 
  int[] T) {

  real total; 
  real try1[2];
  real try2[2]; 
  real y_hat[t+1, n_group*3];
  
  for(i in 1:n_group){
    y_hat[1,1+3*(i-1)] = y[1+3*(i-1)];
    y_hat[1,2+3*(i-1)] = y[2+3*(i-1)];
    y_hat[1,3+3*(i-1)] = y[3+3*(i-1)];
  }
  
  for(time in 1:t){
    for(i in 1:n_group){
    total = 0.0;
    for(j in 1:n_group){
      total = total + beta[i,j] * y_hat[time,2+3*(j-1)]/T[j];
    }
    try1[1] = total * y_hat[time,1+3*(i-1)];
    try1[2] = y_hat[time,1+3*(i-1)];
    try2[1] = gamma[i] * y_hat[time,2+3*(i-1)];
    try2[2] = y_hat[time,2+3*(i-1)];
    y_hat[time+1,1+3*(i-1)] = y_hat[time,1+3*(i-1)] - min(try1);
    y_hat[time+1,2+3*(i-1)] = y_hat[time,2+3*(i-1)] + min(try1) - min(try2);
    y_hat[time+1,3+3*(i-1)] = y_hat[time,3+3*(i-1)] + min(try2);
    }
  }
  return y_hat[2:(t+1), 1:(n_group*3)];
  }
}

data{
  int<lower = 1> n_time1;      // number of days observed
  int<lower = 1> n_time2;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int Y1[n_group, n_group];    // number of contact by age group
  int Y2[n_group, n_group];    // number of contact by age group
  int Z[n_time1+n_time2, n_group];       // number of infected individuals each day by age
  real agegroup_i[n_group*n_group];    // age group of i
  real agegroup_j[n_group*n_group];    // age group of j
  int T[n_group];               // total number of participants by age group 
  int N[n_group,2];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time1+n_time2];            // time points observed
}
  
transformed data {
  int n_cat = n_group*n_group; // number of categories by age
  row_vector[2] agegroup[n_cat];
  row_vector[2] agegroup_sym[n_cat];
  real logU1[n_group, n_group];
  real logU2[n_group, n_group];
  for(i in 1:n_cat){
    agegroup[i, 1] = agegroup_i[i];
    agegroup[i, 2] = agegroup_j[i];
    agegroup_sym[i, 1] = agegroup_j[i];
    agegroup_sym[i, 2] = agegroup_i[i];
  }
  for(i in 1:n_group){
    for(j in 1:n_group){
      logU1[i,j] = log(T[i]) + log(N[j,1]);
      logU2[i,j] = log(T[i]) + log(N[j,2]);
    }
  }
}

parameters{
  // Transmission model parameters
  real<lower=0> theta;           // global dispersion parameter 
  real alpha[2];                    // global intercept
  real<lower=0> rho[2];             // length-scale
  real<lower=0> sigma[2];           // amplitude
  real x_delta[n_cat,2];         // gp function

  // Epidemic model parameters
  real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
  row_vector<lower = 0>[n_group] gamma;   // rate of recovery
  real<lower = 0, upper = 1> s0[n_group]; // initial fraction of susceptible 
  real<lower = 0, upper = 1> s_vac; // fraction of susceptible in vacation 
}
  
transformed parameters{
  // Transmission model transformed parameters first part
  matrix[n_cat, n_cat] cov1 =  cov_exp_quad(agegroup, agegroup, sigma[1], rho[1]) + cov_exp_quad(agegroup_sym, agegroup, sigma[1], rho[1])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov1 = cholesky_decompose(cov1);
  vector[n_cat] x_vec1 = L_cov1 * to_vector(x_delta[,1]);
  matrix[n_group, n_group] f1 = to_matrix(x_vec1, n_group, n_group, 0);
  matrix[n_group, n_group] c1;
  matrix[n_group, n_group] mu_age1;
  
  // Transmission model transformed parameters second part
  matrix[n_cat, n_cat] cov2 =  cov_exp_quad(agegroup, agegroup, sigma[2], rho[2]) + cov_exp_quad(agegroup_sym, agegroup, sigma[2], rho[2])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov2 = cholesky_decompose(cov2);
  vector[n_cat] x_vec2 = L_cov2 * to_vector(x_delta[,2]);
  matrix[n_group, n_group] f2 = to_matrix(x_vec2, n_group, n_group, 0);
  matrix[n_group, n_group] c2;
  matrix[n_group, n_group] mu_age2;

  // Epidemic model transformed parameters
  real<lower = 0> y_hat[(n_time1+n_time2), n_difeq * n_group]; // solution from the ODE solver
  real<lower = 0> y_init[n_difeq * n_group];    // initial conditions for both fractions of S, I and R
  real<lower = 0> y_term[n_difeq * n_group];    // conditions for both fractions of S, I and R after holiday 
  matrix[n_group, n_group] beta1;
  matrix[n_group, n_group] beta2;
  
  for(i in 1:n_group){
    for(j in 1:n_group){
      c1[i,j] = exp(alpha[1] + f1[i,j]);
      mu_age1[i,j] = exp(logU1[i,j] + alpha[1] + f1[i,j]);
      c2[i,j] = exp(alpha[2] + f2[i,j]);
      mu_age2[i,j] = exp(logU2[i,j] + alpha[2] + f2[i,j]);
    }
    y_init[1+3*(i-1)] = s0[i]*(1-s_vac)*T[i];
    y_init[2+3*(i-1)] = (1 - s0[i])*T[i];
    y_init[3+3*(i-1)] = 0;
  }

  beta1 = c1 * diag_matrix(to_vector(T)) * q;
  beta2 = c2 * diag_matrix(to_vector(T)) * q;
  y_hat[1:n_time1,] = SIR(n_time1, y_init, beta1, gamma, n_group, T);
  for(i in 1:n_group){
    y_term[1+3*(i-1)] = y_hat[n_time1,1+3*(i-1)] + (s0[i]*s_vac)*T[i];
    y_term[2+3*(i-1)] = y_hat[n_time1,2+3*(i-1)];
    y_term[3+3*(i-1)] = y_hat[n_time1,3+3*(i-1)];
  }
  y_hat[(n_time1+1):(n_time1+n_time2),] = SIR(n_time2, y_term, beta2, gamma, n_group, T);
}
  
model{
  real lambda[(n_time1+n_time2), n_group];      //poisson parameter
  
  //priors
  // Transmission model priors
  to_vector(alpha) ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  to_vector(rho) ~ inv_gamma(5, 5);
  to_vector(sigma) ~ std_normal();
  to_array_1d(x_delta) ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(.5, 0.5);
  s_vac ~ beta(2, 5);
  
  // Likelihood
  // Transmission model likelihood
  for(i in 1:n_group){
    for(j in 1:n_group){
      Y1[i,j] ~ neg_binomial_2(mu_age1[i,j], theta);
      Y2[i,j] ~ neg_binomial_2(mu_age2[i,j], theta);
    }
  }
  //Epidemic model likelihood
  for(t in 1:(n_time1+n_time2)){
    for(i in 1:n_group){
      lambda[t,i] = y_hat[t,2+3*(i-1)];
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }
}

//generated quantities {
//  real R_0[n_group, n_group];  // Basic reproduction number
//  for(i in 1:n_group){
//    for(j in 1:n_group){
//      R_0[i,j] = beta[i,j]/gamma[i];   
//    }
//  }
//}
"


model.overall.agegroup.deterministic.h1n1	<-  "
functions {
  real[,] SIR(
  int t,           // time
  real[] y,     // system state {susceptible,infected,recovered}
  matrix beta,     // parameters
  row_vector gamma,
  int n_group, 
  int[] T) {

  real total; 
  real try1[2];
  real try2[2]; 
  real y_hat[t+1, n_group*3];
  
  for(i in 1:n_group){
    y_hat[1,1+3*(i-1)] = y[1+3*(i-1)];
    y_hat[1,2+3*(i-1)] = y[2+3*(i-1)];
    y_hat[1,3+3*(i-1)] = y[3+3*(i-1)];
  }
  
  for(time in 1:t){
    for(i in 1:n_group){
    total = 0.0;
    for(j in 1:n_group){
      total = total + beta[i,j] * y_hat[time,2+3*(j-1)]/T[j];
    }
    try1[1] = total * y_hat[time,1+3*(i-1)];
    try1[2] = y_hat[time,1+3*(i-1)];
    try2[1] = gamma[i] * y_hat[time,2+3*(i-1)];
    try2[2] = y_hat[time,2+3*(i-1)];
    y_hat[time+1,1+3*(i-1)] = y_hat[time,1+3*(i-1)] - min(try1);
    y_hat[time+1,2+3*(i-1)] = y_hat[time,2+3*(i-1)] + min(try1) - min(try2);
    y_hat[time+1,3+3*(i-1)] = y_hat[time,3+3*(i-1)] + min(try2);
    }
  }
  return y_hat[2:(t+1), 1:(n_group*3)];
  }
}

data{
  int<lower = 1> n_time1;      // number of days observed
  int<lower = 1> n_time2;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int Y1[n_group, n_group];    // number of contact by age group
  int Y2[n_group, n_group];    // number of contact by age group
  int Z[n_time1+n_time2, n_group];       // number of infected individuals each day by age
  real agegroup_i[n_group*n_group];    // age group of i
  real agegroup_j[n_group*n_group];    // age group of j
  int T[n_group];               // total number of participants by age group 
  int N[n_group,2];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time1+n_time2];            // time points observed
}
  
transformed data {
  int n_cat = n_group*n_group; // number of categories by age
  row_vector[2] agegroup[n_cat];
  row_vector[2] agegroup_sym[n_cat];
  real logU1[n_group, n_group];
  real logU2[n_group, n_group];
  for(i in 1:n_cat){
    agegroup[i, 1] = agegroup_i[i];
    agegroup[i, 2] = agegroup_j[i];
    agegroup_sym[i, 1] = agegroup_j[i];
    agegroup_sym[i, 2] = agegroup_i[i];
  }
  for(i in 1:n_group){
    for(j in 1:n_group){
      logU1[i,j] = log(T[i]) + log(N[j,1]);
      logU2[i,j] = log(T[i]) + log(N[j,2]);
    }
  }
}

parameters{
  // Transmission model parameters
  real<lower=0> theta;           // global dispersion parameter 
  real alpha[2];                    // global intercept
  real<lower=0> rho[2];             // length-scale
  real<lower=0> sigma[2];           // amplitude
  real x_delta[n_cat,2];         // gp function

  // Epidemic model parameters
  real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
  row_vector<lower = 0>[n_group] gamma;   // rate of recovery
  real<lower = 0, upper = 1> s0[n_group]; // initial fraction of susceptible 
  real<lower = 0, upper = 1> s_vac; // fraction of susceptible in vacation 
}
  
transformed parameters{
  // Transmission model transformed parameters first part
  matrix[n_cat, n_cat] cov1 =  cov_exp_quad(agegroup, agegroup, sigma[1], rho[1]) + cov_exp_quad(agegroup_sym, agegroup, sigma[1], rho[1])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov1 = cholesky_decompose(cov1);
  vector[n_cat] x_vec1 = L_cov1 * to_vector(x_delta[,1]);
  matrix[n_group, n_group] f1 = to_matrix(x_vec1, n_group, n_group, 0);
  matrix[n_group, n_group] c1;
  matrix[n_group, n_group] mu_age1;
  
  // Transmission model transformed parameters second part
  matrix[n_cat, n_cat] cov2 =  cov_exp_quad(agegroup, agegroup, sigma[2], rho[2]) + cov_exp_quad(agegroup_sym, agegroup, sigma[2], rho[2])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov2 = cholesky_decompose(cov2);
  vector[n_cat] x_vec2 = L_cov2 * to_vector(x_delta[,2]);
  matrix[n_group, n_group] f2 = to_matrix(x_vec2, n_group, n_group, 0);
  matrix[n_group, n_group] c2;
  matrix[n_group, n_group] mu_age2;

  // Epidemic model transformed parameters
  real<lower = 0> y_hat[(n_time1+n_time2), n_difeq * n_group]; // solution from the ODE solver
  real<lower = 0> y_init[n_difeq * n_group];    // initial conditions for both fractions of S, I and R
  real<lower = 0> y_term[n_difeq * n_group];    // conditions for both fractions of S, I and R after holiday 
  matrix[n_group, n_group] beta1;
  matrix[n_group, n_group] beta2;
  
  for(i in 1:n_group){
    for(j in 1:n_group){
      c1[i,j] = exp(alpha[1] + f1[i,j]);
      mu_age1[i,j] = exp(logU1[i,j] + alpha[1] + f1[i,j]);
      c2[i,j] = exp(alpha[2] + f2[i,j]);
      mu_age2[i,j] = exp(logU2[i,j] + alpha[2] + f2[i,j]);
    }
    y_init[1+3*(i-1)] = s0[i]*(1-s_vac)*T[i];
    y_init[2+3*(i-1)] = (1 - s0[i])*T[i];
    y_init[3+3*(i-1)] = 0;
  }

  beta1 = c1 * diag_matrix(to_vector(T)) * q;
  beta2 = c2 * diag_matrix(to_vector(T)) * q;
  y_hat[1:n_time1,] = SIR(n_time1, y_init, beta1, gamma, n_group, T);
  for(i in 1:n_group){
    y_term[1+3*(i-1)] = y_hat[n_time1,1+3*(i-1)] + (s0[i]*s_vac)*T[i];
    y_term[2+3*(i-1)] = y_hat[n_time1,2+3*(i-1)];
    y_term[3+3*(i-1)] = y_hat[n_time1,3+3*(i-1)];
  }
  y_hat[(n_time1+1):(n_time1+n_time2),] = SIR(n_time2, y_term, beta2, gamma, n_group, T);
}
  
model{
  real lambda[(n_time1+n_time2), n_group];      //poisson parameter
  
  //priors
  // Transmission model priors
  to_vector(alpha) ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  to_vector(rho) ~ inv_gamma(5, 5);
  to_vector(sigma) ~ std_normal();
  to_array_1d(x_delta) ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(.5, 0.5);
  s_vac ~ beta(2, 5);
  
  // Likelihood
  // Transmission model likelihood
  for(i in 1:n_group){
    for(j in 1:n_group){
      Y1[i,j] ~ neg_binomial_2(mu_age1[i,j], theta);
      Y2[i,j] ~ neg_binomial_2(mu_age2[i,j], theta);
    }
  }
  //Epidemic model likelihood
  for(t in 1:(n_time1+n_time2)){
    for(i in 1:n_group){
      lambda[t,i] = y_hat[t,2+3*(i-1)];
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }

}

//generated quantities {
//  real R_0[n_group, n_group];  // Basic reproduction number
//  for(i in 1:n_group){
//    for(j in 1:n_group){
//      R_0[i,j] = beta[i,j]/gamma[i];   
//    }
//  }
//}
"


init_f <- function () list("alpha" = c(-100,-100), "s0" = rep(.99999999,4))
fit_overall = stan(model_code = model.overall.agegroup.deterministic.h1n1, data = data_list, chains = 1, iter = 10000, warmup = 1000,
                  seed=4567, init = init_f)

qstatlibrary(bayesplot)

#// Transmission model parameters
#real<lower=0> theta[2];           // global dispersion parameter 
#real alpha[2];                    // global intercept
#real<lower=0> rho[2];             // length-scale
#real<lower=0> sigma[2];           // amplitude
#real x_delta[n_cat,2];         // gp function
#
#// Epidemic model parameters
#real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
#row_vector<lower = 0>[n_group] gamma;   // rate of recovery
#real<lower = 0, upper = 1> s0[n_group]; // initial fraction of susceptible 
#real<lower = 0, upper = 1> s_vac; // fraction of susceptible in vacation 
#}

traceplot(fit_overall, pars = c("x_delta"))

mcm
sum = summary(fit_overall)
rownames(sum$summary) == "try[1]"#[1100:1300]
summary(posterior$q)
which(rownames(sum$summary) == "y_hat[1,1]")#[1100:1300]
iter = 10000; warmup = 1000
sum$summary[100:200,]
sum$summary[1200:1260,]
posterior = extract(fit_overall)

data_list$Y1
posterior$mu_age1[2000,,]
posterior$beta1
c_est = as.vector(sapply(1:data_list$n_group, function(x){
  as.vector(posterior$c1[,x,])
}))
age_obs = levels(contact$AGE_i)
df = data.table(c = c_est, agegroup_i = rep(rep(age_obs, each = length(age_obs)), each = (iter-warmup)*3),
                agegroup_j = rep(rep(age_obs, length(age_obs)), each = (iter-warmup)*3)) %>%
  group_by(agegroup_i, agegroup_j) %>%
  summarise(c.m = mean(c), c.l95 = quantile(c, prob = .025), c.u95 = quantile(c, prob = .975))
ggplot(df, aes(x = agegroup_i, y = agegroup_j)) +
  geom_raster(aes(fill = c.m))


posterior$beta1
posterior$mu_age1[,1,1]

posterior$mu_age1[700,,]
data_list$Y1
S.m = matrix(ncol = n_age, nrow = n_time, 0); S.l95 = S.m; S.u95 = S.m
I.m = matrix(ncol = n_age, nrow = n_time, 0); I.l95 = S.m; I.u95 = S.m;
R.m = S.m; R.l95 = R.m; R.u95 = S.m;
for(t in 1:n_time){
  for(i in 1:n_age){
    S.m[t, i] = median(posterior$y_hat[,t, 1+3*(i-1)] * data_list$T[i])
    S.l95[t, i] = quantile(posterior$y_hat[,t, 1+3*(i-1)]* data_list$T[i], prob = 0.025)
    S.u95[t, i] = quantile(posterior$y_hat[,t, 1+3*(i-1)]* data_list$T[i], prob = 0.975)
    I.m[t, i] = median(posterior$y_hat[,t, 2+3*(i-1)]* data_list$T[i])
    I.l95[t, i] = quantile(posterior$y_hat[,t, 2+3*(i-1)]* data_list$T[i], prob = 0.025)
    I.u95[t, i] = quantile(posterior$y_hat[,t, 2+3*(i-1)]* data_list$T[i], prob = 0.975)
    R.m[t, i] = median(posterior$y_hat[,t, 3+3*(i-1)]* data_list$T[i])
    R.l95[t, i] = quantile(posterior$y_hat[,t, 3+3*(i-1)]* data_list$T[i], prob = 0.025)
    R.u95[t, i] = quantile(posterior$y_hat[,t, 3+3*(i-1)]* data_list$T[i], prob = 0.975)
  }
}

SIR.stan = data.table(t = rep(1:n_time, n_age*3), 
                      state = rep(c("S", "I", "R"), each = n_time*n_age),
                      median = c(as.vector(S.m), as.vector(I.m), as.vector(R.m)),
                      l95 = c(as.vector(S.l95), as.vector(I.l95), as.vector(R.l95)),
                      u95 = c(as.vector(S.u95), as.vector(I.u95), as.vector(R.u95)),
                      age = rep(rep(1:n_age, each = n_time), 3), 
                      Z = unlist(c(data_list$Z))) 

ggplot(subset(SIR.stan, state == "I"), aes(x = t, fill = state)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.5) +
  geom_line(aes(y = median, col = state)) +
  facet_wrap(~age) +
  geom_line(aes(y = unlist(c(data_list$Z))))



posterior$beta1[1,1,1]
posterior$y_hat[1,]


time = 1; try1 = matrix(nrow = data_list$n_group, ncol = 2, 0); total = 0
i = 1
for(j in 1:data_list$n_group){
  try1[j,1] = posterior$y_hat[1,time,2+3*(j-1)]/data_list$T[j];
  try1[j,2] = posterior$beta1[1,i,j] * posterior$y_hat[1,time,2+3*(j-1)]/data_list$T[j];
  total = total +  min(try1);
}
posterior$y_hat[1,time,1+3*(i-1)] - total * posterior$y_hat[1,time,1+3*(i-1)];
posterior$y_hat[1,time,2+3*(i-1)] + total * posterior$y_hat[1,time,1+3*(i-1)]- posterior$gamma[i] * posterior$y_hat[1,time,2+3*(i-1)];
posterior$y_hat[1,time,3+3*(i-1)] + posterior$gamma[i] * posterior$y_hat[1,time,2+3*(i-1)];

posterior$y_hat[1,time+1,1+3*(i-1)]
posterior$y_hat[1,time+1,2+3*(i-1)]
posterior$y_hat[1,time+1,3+3*(i-1)]

model.overall.agegroup.deterministic.h1n1	<-  "
functions {
  real[,] SIR(
  int t,           // time
  real[] y,     // system state {susceptible,infected,recovered}
  matrix beta,     // parameters
  row_vector gamma,
  int n_group, 
  int[] T) {

  real total; 
  real try1[2];
  real try2[2]; 
  real y_hat[t+1, n_group*3];
  
  for(i in 1:n_group){
    y_hat[1,1+3*(i-1)] = y[1+3*(i-1)];
    y_hat[1,2+3*(i-1)] = y[2+3*(i-1)];
    y_hat[1,3+3*(i-1)] = y[3+3*(i-1)];
  }
  
  for(time in 1:t){
    for(i in 1:n_group){
    total = 0.0;
    for(j in 1:n_group){
      total = total + beta[i,j] * y_hat[time,2+3*(j-1)]/T[j];
    }
    try1[1] = total * y_hat[time,1+3*(i-1)];
    try1[2] = y_hat[time,1+3*(i-1)];
    try2[1] = gamma[i] * y_hat[time,2+3*(i-1)];
    try2[2] = y_hat[time,2+3*(i-1)];
    y_hat[time+1,1+3*(i-1)] = y_hat[time,1+3*(i-1)] - min(try1);
    y_hat[time+1,2+3*(i-1)] = y_hat[time,2+3*(i-1)] + min(try1) - min(try2);
    y_hat[time+1,3+3*(i-1)] = y_hat[time,3+3*(i-1)] + min(try2);
    }
  }
  return y_hat[2:(t+1), 1:(n_group*3)];
  }
}

data{
  int<lower = 1> n_time1;      // number of days observed
  int<lower = 1> n_time2;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int Y1[n_group, n_group];    // number of contact by age group
  int Y2[n_group, n_group];    // number of contact by age group
  int Z[n_time1+n_time2, n_group];       // number of infected individuals each day by age
  real agegroup_i[n_group*n_group];    // age group of i
  real agegroup_j[n_group*n_group];    // age group of j
  int T[n_group];               // total number of participants by age group 
  int N[n_group,2];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time1+n_time2];            // time points observed
}
  
transformed data {
  real x_r[0];
  int x_i[1];
  int n_cat = n_group*n_group; // number of categories by age
  row_vector[2] agegroup[n_cat];
  row_vector[2] agegroup_sym[n_cat];
  real logU1[n_group, n_group];
  real logU2[n_group, n_group];
  x_i[1] = n_group;
  for(i in 1:n_cat){
    agegroup[i, 1] = agegroup_i[i];
    agegroup[i, 2] = agegroup_j[i];
    agegroup_sym[i, 1] = agegroup_j[i];
    agegroup_sym[i, 2] = agegroup_i[i];
  }
  for(i in 1:n_group){
    for(j in 1:n_group){
      logU1[i,j] = log(T[i]) + log(N[j,1]);
      logU2[i,j] = log(T[i]) + log(N[j,2]);
    }
  }
}

parameters{
  // Transmission model parameters
  real<lower=0> theta[2];           // global dispersion parameter 
  real alpha[2];                    // global intercept
  real<lower=0> rho[2];             // length-scale
  real<lower=0> sigma[2];           // amplitude
  real x_delta[n_cat,2];         // gp function

  // Epidemic model parameters
  real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
  row_vector<lower = 0>[n_group] gamma;   // rate of recovery
  real<lower = 0, upper = 1> s0[n_group]; // initial fraction of susceptible 
  real<lower = 0, upper = 1> s_vac; // fraction of susceptible in vacation 
}
  
transformed parameters{
  // Transmission model transformed parameters first part
  matrix[n_cat, n_cat] cov1 =  cov_exp_quad(agegroup, agegroup, sigma[1], rho[1]) + cov_exp_quad(agegroup_sym, agegroup, sigma[1], rho[1])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov1 = cholesky_decompose(cov1);
  vector[n_cat] x_vec1 = L_cov1 * to_vector(x_delta[,1]);
  matrix[n_group, n_group] f1 = to_matrix(x_vec1, n_group, n_group, 0);
  matrix[n_group, n_group] c1;
  matrix[n_group, n_group] mu_age1;
  
  // Transmission model transformed parameters second part
  matrix[n_cat, n_cat] cov2 =  cov_exp_quad(agegroup, agegroup, sigma[2], rho[2]) + cov_exp_quad(agegroup_sym, agegroup, sigma[2], rho[2])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov2 = cholesky_decompose(cov2);
  vector[n_cat] x_vec2 = L_cov2 * to_vector(x_delta[,2]);
  matrix[n_group, n_group] f2 = to_matrix(x_vec2, n_group, n_group, 0);
  matrix[n_group, n_group] c2;
  matrix[n_group, n_group] mu_age2;

  // Epidemic model transformed parameters
  //real y_hat[(n_time1+n_time2), n_difeq * n_group]; // solution from the ODE solver
  real y_hat[n_time2+n_time1, n_difeq * n_group]; // solution from the ODE solver
  real y_init[n_difeq * n_group];    // initial conditions for both fractions of S, I and R
  //real y_term[n_difeq * n_group];    // conditions for both fractions of S, I and R after holiday 
  matrix[n_group, n_group] beta1;
  matrix[n_group, n_group] beta2;
  real beta_gamma1[n_group + n_group*n_group]; 
  real beta_gamma2[n_group + n_group*n_group];

// real y_hat[t+1, n_group*3];
  
  for(i in 1:n_group){
    for(j in 1:n_group){
      c1[i,j] = exp(alpha[1] + f1[i,j]);
      mu_age1[i,j] = logU1[i,j] + c1[i,j];
      c2[i,j] = exp(alpha[2] + f2[i,j]);
      mu_age2[i,j] = logU2[i,j] + c2[i,j];
    }
    y_init[1+3*(i-1)] = s0[i]*T[i];//*(1-s_vac);
    y_init[2+3*(i-1)] = (1 - s0[i])*T[i];
    y_init[3+3*(i-1)] = 0;
  }

  beta1 = c1 * diag_matrix(to_vector(T)) * q;
  beta2 = c2 * diag_matrix(to_vector(T)) * q;
  beta_gamma1 = to_array_1d(append_row(gamma, beta1)');
  beta_gamma2 = to_array_1d(append_row(gamma, beta2)');

   y_hat[1:(n_time1+n_time2),] = SIR(n_time1+n_time2, y_init, beta2, gamma, n_group, T);

}
  
model{
  real lambda[(n_time1+n_time2), n_group];      //poisson parameter
  
  //priors
  // Transmission model priors
  alpha ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  rho ~ inv_gamma(5, 5);
  sigma ~ std_normal();
  to_array_1d(x_delta) ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(10, 0.1);
  s_vac ~ beta(2, 5);
  
  // Likelihood
  // Transmission model likelihood
  for(i in 1:n_group){
    for(j in 1:n_group){
      Y1[i,j] ~ normal(mu_age1[i,j], theta);
      Y2[i,j] ~ normal(mu_age2[i,j], theta);
    }
  }

}

//generated quantities {
//  real R_0[n_group, n_group];  // Basic reproduction number
//  for(i in 1:n_group){
//    for(j in 1:n_group){
//      R_0[i,j] = beta[i,j]/gamma[i];   
//    }
//  }
//}
"