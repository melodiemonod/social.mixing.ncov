### 1. ESTIMATE THE CONTACT MATRIX ###
# WITH AGE CATEGORY, AND UNOBSERVED VALUE
model.contact.matrix	<-  "
    data{
      int<lower=1> n_group;  // number of age group with observed contact
      int<lower=1> N_group; // number of age group with and without observed contact
      int y[N_group, N_group];      // number of contacts
      real agegroup_i[N_group*N_group]; // age group of participant
      real agegroup_j[N_group*N_group]; // age group of contact
      real logU[N_group, N_group];
     }
    
    transformed data {
      int n_cat = N_group*N_group; // number of categories 
      row_vector[2] age[n_cat];
      row_vector[2] age_sym[n_cat]; // add symmetric age for the kernel
      for(i in 1:n_cat){
        age[i, 1] = agegroup_i[i];
        age[i, 2] = agegroup_j[i];
        age_sym[i, 1] = agegroup_j[i];
        age_sym[i, 2] = agegroup_i[i];
      }
    }
    
    parameters{
      // Global dispersion parameter 
      real<lower=0> theta;
  
      // Global intercept
      real beta;

      // Parameters for covariance matrix
      real<lower=0> rho;
      real<lower=0> alpha;
      vector[n_cat] x_delta;
    }
  
    transformed parameters{
      matrix[n_cat, n_cat] cov =  cov_exp_quad(age, age, alpha, rho) + cov_exp_quad(age_sym, age, alpha, rho) + diag_matrix(rep_vector(1e-10, n_cat));
      matrix[n_cat, n_cat] L_cov = cholesky_decompose(cov);
      vector[n_cat] x_vec = L_cov * x_delta;
      matrix[N_group, N_group] x = to_matrix(x_vec, N_group, N_group, 0);
      matrix<lower=0> [N_group, N_group] c;
      matrix<lower=0> [N_group, N_group] E_Y;
      for(i in 1:N_group){
        for(j in 1:N_group){
          c[i,j] = exp(beta + x[i,j]);
          E_Y[i,j] = exp(logU[i,j] + beta + x[i,j]);
        }
      }
    }

    model{
      beta ~ normal(0,0.001);
      theta ~ lognormal(0,0.001);

      rho ~ inv_gamma(5, 5);
      alpha ~ std_normal();
      x_delta ~ std_normal();
          
      for(i in 1:n_group){ // don't take into account unobserved age group
        for(j in 1:n_group){
          y[i,j] ~ neg_binomial_2(E_Y[i,j], theta);
        }
      }
    }
  "

# WITH AGE GROUP AND NO EXTRAPOLATION
model.transmission.agegroup	<-  "
    data{
      int<lower=1> n_group;  // number of age group
      int y[n_group, n_group];      // the outcome variable
      real agegroup_i[n_group*n_group]; // age group of i
      real agegroup_j[n_group*n_group]; // age group of j
      int T[n_group]; // total number of participants by age group 
      int N[n_group]; // total number of individuals by age group 
     }
    
    transformed data {
      int n_cat = n_group*n_group; // number of categories by age group and gender
      row_vector[2] age[n_cat];
      row_vector[2] age_sym[n_cat];
      real logU[n_group, n_group];
      for(i in 1:n_cat){
        age[i, 1] = agegroup_i[i];
        age[i, 2] = agegroup_j[i];
        age_sym[i, 1] = agegroup_j[i];
        age_sym[i, 2] = agegroup_i[i];
      }
      for(i in 1:n_group){
        for(j in 1:n_group){
          logU[i,j] = log(T[i]) + log(N[j]);
        }
      }
    }
    
    parameters{
      // Global dispersion parameter 
      real<lower=0> theta;
  
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
      matrix[n_group, n_group] x = to_matrix(x_vec, n_group, n_group, 0);
      matrix<lower=0> [n_group, n_group] c;
      matrix<lower=0> [n_group, n_group] mu;
      for(i in 1:n_group){
        for(j in 1:n_group){
          c[i,j] = exp(beta + x[i,j]);
          mu[i,j] = exp(logU[i,j] + beta + x[i,j]);
        }
      }
    }

    model{
      beta ~ normal(0,0.001);
      theta ~ lognormal(0,0.001);

      rho ~ inv_gamma(5, 5);
      alpha ~ std_normal();
      x_delta ~ std_normal();
          
      for(i in 1:n_group){
        for(j in 1:n_group){
          y[i,j] ~ normal(mu[i,j], theta);
        }
      }
    }
  "


## model with age and sex heterogenetiy 
model.transmission.agegroup.sex	<-  "
    data{
      int<lower=1> n_age;  // number of age group
      int<lower=1> n_group; // number of gender group
      int y[n_age, n_age, n_group];      // the outcome variable
      real agegroup_i[n_age*n_age]; // age group of i
      real agegroup_j[n_age*n_age]; // age group of j
      int T[n_age, 2]; // total number of participants by age group 
      int N[n_age, 2]; // total number of individuals by age group 
     }
    
    transformed data {
      int n_cat = n_age*n_age; // number of categories by age group and gender
      row_vector[2] age[n_cat];
      row_vector[2] age_sym[n_cat];
      int U[n_age, n_age, 4];
      for(i in 1:n_cat){
        age[i, 1] = agegroup_i[i];
        age[i, 2] = agegroup_j[i];
        age_sym[i, 1] = agegroup_j[i];
        age_sym[i, 2] = agegroup_i[i];
      }
      for(i in 1:n_age){
        for(j in 1:n_age){
          U[i,j,1] = N[j,1] * T[i,1]; //U_ij^FF
          U[i,j,2] = N[j,2] * T[i,2]; //U_ij^MM
          U[i,j,3] = N[j,1] * T[i,2]; //U_ij^FM
          U[i,j,4] = N[j,2] * T[i,1]; //U_ij^MF
        }
      }
    }
    
    parameters{
      // Global dispersion parameter 
      real theta[n_age];
  
      // Global intercept
      real beta;

      // Parameters for covariance matrix
      real<lower=0> rho_FF;
      real<lower=0> rho_MM;
      real<lower=0> rho_FM; // = rho_MF
      real<lower=0> alpha_FF;
      real<lower=0> alpha_MM;
      real<lower=0> alpha_FM; // = alpha_MF
      vector[n_cat] x_delta_FF;
      vector[n_cat] x_delta_MM;
      vector[n_cat] x_delta_FM; // = x_delta_MF'
    }
  
    transformed parameters{
      matrix[n_cat, n_cat] cov_FF =  cov_exp_quad(age, age, alpha_FF, rho_FF) + cov_exp_quad(age_sym, age, alpha_FF, rho_FF)+ diag_matrix(rep_vector(1e-10, n_cat));
      matrix[n_cat, n_cat] L_cov_FF = cholesky_decompose(cov_FF);
      vector[n_cat] x_vec_FF = L_cov_FF * x_delta_FF;
      matrix[n_age, n_age] x_FF = to_matrix(x_vec_FF, n_age, n_age, 0);
      
      matrix[n_cat, n_cat] cov_MM =  cov_exp_quad(age, age, alpha_MM, rho_MM) + cov_exp_quad(age_sym, age, alpha_MM, rho_MM)  + diag_matrix(rep_vector(1e-10, n_cat));
      matrix[n_cat, n_cat] L_cov_MM = cholesky_decompose(cov_MM);
      vector[n_cat] x_vec_MM = L_cov_MM * x_delta_MM;
      matrix[n_age, n_age] x_MM = to_matrix(x_vec_MM, n_age, n_age, 0);
      
      
      matrix[n_cat, n_cat] cov_FM =  cov_exp_quad(age, age, alpha_FM, rho_FM) + diag_matrix(rep_vector(1e-10, n_cat));
      matrix[n_cat, n_cat] L_cov_FM = cholesky_decompose(cov_FM);
      vector[n_cat] x_vec_FM = L_cov_FM * x_delta_FM;
      matrix[n_age, n_age] x_FM = to_matrix(x_vec_FM, n_age, n_age, 0);
      
      matrix[n_age, n_age] x_MF = x_FM';
      
      real c[n_age, n_age, n_group];
      real mu[n_age, n_age, n_group];
      for(i in 1:n_age){
        for(j in 1:n_age){
          c[i,j,1] = exp(beta + x_FF[i,j]);
          mu[i,j,1] = exp(log(U[i,j,1]) + beta + x_FF[i,j]);
          c[i,j,2] = exp(beta + x_MM[i,j]);
          mu[i,j,2] = exp(log(U[i,j,2]) + beta + x_MM[i,j]);
          c[i,j,3] = exp(beta + x_FM[i,j]);
          mu[i,j,3] = exp(log(U[i,j,3]) + beta + x_FM[i,j]);
          c[i,j,4] = exp(beta + x_MF[i,j]);
          mu[i,j,4] = exp(log(U[i,j,4]) + beta + x_MF[i,j]);
        }
      }
    }

    model{
      beta ~ normal(0,0.001);
      to_vector(theta) ~ lognormal(0,0.001);

      rho_FF ~ inv_gamma(5, 5);
      rho_MM ~ inv_gamma(5, 5);
      rho_FM ~ inv_gamma(5, 5);
      alpha_FF ~ std_normal();
      alpha_MM ~ std_normal();
      alpha_FM ~ std_normal();
      x_delta_FF ~ std_normal();
      x_delta_MM ~ std_normal();
      x_delta_FM ~ std_normal();
      
      for(i in 1:n_age){
        for(j in 1:n_age){
          for(k in 1:n_group){
            y[i,j,k] ~ neg_binomial_2(mu[i,j,k], theta[i]);
          }
        }
      }
    }
  "


### 2. EPIDEMIC MODEL ##
## DETERMINISTIC WITH AGE 
model.epi.deterministic <- "
functions {
  real[] SIR(
  real t,  // time
  real[] y_raw,           // system state {susceptible,infected,recovered}
  real[] theta,       // parameters
  real[] x_r,
  int[] x_i) {
  
  int n_age = x_i[1];
  real dy_dt[3*n_age];
  matrix[3, n_age] y = to_matrix(y_raw, 3, n_age);
  matrix[n_age, n_age] beta = to_matrix(theta[(n_age+1):(n_age+n_age*n_age)], n_age, n_age); 
  real total;

  for(i in 1:n_age){
    dy_dt[1+3*(i-1)] = - beta[i,] * to_vector(y[2,]) * y[1,i];
    dy_dt[2+3*(i-1)] = beta[i,] * to_vector(y[2,]) * y[1,i] - theta[i] * y[2,i];
    dy_dt[3+3*(i-1)] = theta[i] * y[2,i];
  }

  return dy_dt;
  }
}
data {
  int<lower = 1> n_time;       // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower = 1> n_age;      // number of age groups
  int<lower = 1> T[n_age];       // population 
  int Z[n_time, n_age];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_time];         // time points observed
  }
  
  transformed data {
  real x_r[0];
  int x_i[1];
  x_i[1] = n_age;
}
  
  parameters {
  row_vector<lower = 0>[n_age] gamma; // model parameters {gamma}
  real<lower = 0, upper = 1> s0[n_age];  // initial fraction of susceptible individuals
  //real<lower = 0, upper = 1> q; // probability of transmission given a contact
  vector<lower=0>[n_age] beta_raw;
  }
  
  transformed parameters{
  real<lower = 0> y_hat[n_time, n_difeq * n_age]; // solution from the ODE solver
  real<lower = 0> y_init[n_difeq * n_age];    // initial conditions for both fractions of S, I and R
  matrix[n_age, n_age] beta = diag_matrix(beta_raw);
  real theta[n_age + n_age*n_age] = to_array_1d(append_row(gamma, to_matrix(beta))'); 
  
  for(i in 1:n_age){
      y_init[1+3*(i-1)] = s0[i];
      y_init[2+3*(i-1)] = 1 - s0[i];
      y_init[3+3*(i-1)] = 0;
  }
    
  y_hat = integrate_ode_bdf(SIR, y_init, t0, ts, theta, x_r, x_i);
}
  
model {
  real lambda[n_time, n_age];      //poisson parameter
  
  //priors
  //q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(0.5, 0.5);
  
  //likelihood
  for(i in 1:n_age){
  beta_raw[i] ~ lognormal(0,1);
    for(t in 1:n_time){
      lambda[t,i] = y_hat[t,2+3*(i-1)]*T[i];
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }
}
  
generated quantities {
  real R_0[n_age, n_age];  // Basic reproduction number
  for(i in 1:n_age){
    for(j in 1:n_age){
      R_0[i,j] = beta[i,j]/gamma[i];   
    }
  }
}
"

model.epi.deterministic.euler <- "
data {
  int<lower = 1> n_time;       // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower = 1> n_age;      // number of age groups
  int<lower = 1> T[n_age];       // population 
  int Z[n_time, n_age];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_time];         // time points observed
}

parameters {
  real<lower = 0> gamma[n_age] ; // model parameters {gamma}
  real<lower = 0, upper = 1> s0[n_age];  // initial fraction of susceptible individuals
  //real<lower = 0, upper = 1> q; // probability of transmission given a contact
  real<lower=0> beta[n_age];
  }
  
transformed parameters{
  real S[n_age,n_time]; 
  real I[n_age,n_time]; 
  real R[n_age,n_time]; 
  real<lower=0> newI[n_age,n_time-1];
  real<lower=0> newR[n_age,n_time-1];
  
  // Initial 
  for(age in 1:n_age){
      S[age,1] = s0[age]*T[age];
      I[age,1] = (1 - s0[age])*T[age];
      R[age,1] = 0;
  }
    
  // SIR
   for (t in 2:n_time) {
    for (age in 1:n_age) {
      newI[age,t-1] = fmin(S[age,t-1], S[age,t-1]*beta[age]*I[age,t-1]/T[age]);
      
      S[age,t] = S[age,t-1] - newI[age,t-1];
      I[age,t] = I[age,t-1] + newI[age,t-1];
      
      newR[age,t-1] = fmin(I[age,t], I[age,t]*gamma[age]);
      I[age,t] = I[age,t] - newR[age,t-1];
      R[age,t] = R[age,t-1] + newR[age,t-1];
    }
  }
}
  
model {
  real lambda[n_time, n_age];      //poisson parameter
  
  //priors
  //q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(0.5, 0.5);
  to_vector(beta) ~ lognormal(0,1);
  //likelihood
  for(i in 1:n_age){
    for(t in 1:n_time){
      lambda[t,i] = I[i,t];
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }
}
  
//generated quantities {
//  real R_0[n_age, n_age];  // Basic reproduction number
//  for(i in 1:n_age){
//    for(j in 1:n_age){
//      R_0[i,j] = beta[i,j]/gamma[i];   
//    }
//  }
//}
"


model.epi.stochastic <- "
functions {
  real[] SIR(
  real t,  // time
  real[] y_raw,           // system state {susceptible,infected,recovered}
  real[] theta,       // parameters
  real[] x_r,
  int[] x_i) {
  
  int n_age = size(y_raw)/3;
  real dy_dt[3*n_age];
  matrix[3, n_age] y = to_matrix(y_raw, 3, n_age);
  matrix[n_age, n_age] beta = to_matrix(theta[(n_age+1):(n_age+n_age*n_age)], n_age, n_age); 
  real total;

  for(i in 1:n_age){
    total = 0.0;
    for(j in 1:n_age){
      total = total + beta[i,j] * y[2,j];
    }
    dy_dt[1+3*(i-1)] = - total * y[1,i];
    dy_dt[2+3*(i-1)] = total * y[1,i] - theta[i] * y[2,i];
    dy_dt[3+3*(i-1)] = theta[i] * y[2,i];
  }

  return dy_dt;
  }
}
  
data {
  int<lower = 1> n_time;       // number of days observed
  int<lower = 1> n_theta;     // number of model parameters
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower = 1> n_age;      // number of age groups
  int<lower = 1> T[n_age];       // population 
  int Z[n_time, n_age];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_time];         // time points observed
  }
  
transformed data {
  real x_r[0];
  int x_i[0];
  }
  
parameters {
  row_vector<lower = 0>[n_age] gamma; // model parameters {gamma}
  real<lower = 0, upper = 1> s0[n_age];  // initial fraction of susceptible individuals
  real kappa[n_time, n_age];              // logarithm of poisson parameter
  real<lower=0> phi[n_age];              // speed of reversion
  real<lower=0> s_sq[n_age];             // square of instantaneous diffusion term
  vector<lower=0>[n_age] beta_raw;
}
  
transformed parameters{
  real mu[n_time, n_age];    
  real y_hat[n_time, n_difeq * n_age]; // solution from the ODE solver
  real y_init[n_difeq * n_age];    // initial conditions for both fractions of S, I and R
  matrix[n_age, n_age] beta = diag_matrix(beta_raw);
  real theta[n_age + n_age*n_age] = to_array_1d(append_row(gamma, to_matrix(beta))'); 
  real sigma[n_age];                  // variance of OU process

  for(i in 1:n_age){
    y_init[1+3*(i-1)] = s0[i];
    y_init[2+3*(i-1)] = 1 - s0[i];
    y_init[3+3*(i-1)] = 0;
    sigma[i]=(1-exp(-2*phi[i]))*(s_sq[i]/(2*phi[i]));
  }
    
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
  
  for(i in 1:n_age){
    for(t in 1:n_time){
      mu[t,i] = log(y_hat[t,2+3*(i-1)]*T[i]);
    }
  }
}
  
model {
  real lambda[n_time, n_age];      //poisson parameter
  
  //priors
  //q ~ beta(20, 0.1); 
  gamma ~ gamma(0.004,0.02);
  to_vector(s0) ~ beta(.5,.5);
  to_vector(phi) ~ normal(0,10);
  to_vector(s_sq) ~ inv_gamma(0.1,0.1);

  //likelihood
  for(i in 1:n_age){
    beta_raw[i] ~ lognormal(0,1);
    kappa[1,i]~normal(mu[1,i],sigma[i]); // kappa=log(lambda) is an Ornstein-Uhlenbeck process
    for(t in 2:n_time){
      kappa[t,i]~normal(mu[t,i]+(kappa[t-1,i]-mu[t,i])*exp(-phi[i]),sigma[i]);   
      lambda[t,i]=exp(kappa[t,i]);
      //lambda[t,i]=(exp(kappa[t,i]) + exp(kappa[t-1,i]))/2;
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }
}
  
generated quantities {
  real R_0[n_age, n_age];  // Basic reproduction number
  for(i in 1:n_age){
    for(j in 1:n_age){
      R_0[i,j] = beta[i,j]/gamma[i];   
    }
  }
}
"


## TRANSMISSION AND EPIDEMIC MODEL WITH AGE MIXING
## DETERMINISTIC and age group
model.overall.agegroup.deterministic	<-  "
data{
  int<lower = 1> n_time;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int Y[n_group, n_group];    // number of contact by age group
  int Z[n_time, n_group];       // number of infected individuals each day by age
  real agegroup_i[n_group*n_group];    // age group of i
  real agegroup_j[n_group*n_group];    // age group of j
  int T[n_group];               // total number of participants by age group 
  int N[n_group];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time];            // time points observed
}
  
transformed data {
  real x_r[0];
  int x_i[1];
  int n_cat = n_group*n_group; // number of categories by age
  row_vector[2] agegroup[n_cat];
  row_vector[2] agegroup_sym[n_cat];
  real logU[n_group, n_group];
  x_i[1] = n_group;
  for(i in 1:n_cat){
    agegroup[i, 1] = agegroup_i[i];
    agegroup[i, 2] = agegroup_j[i];
    agegroup_sym[i, 1] = agegroup_j[i];
    agegroup_sym[i, 2] = agegroup_i[i];
  }
  for(i in 1:n_group){
    for(j in 1:n_group){
      logU[i,j] = log(T[i]) + log(N[j]);
    }
  }
}

parameters{
  // Transmission model parameters
  real<lower=0> theta;           // global dispersion parameter 
  real alpha;                    // global intercept
  real<lower=0> rho;             // length-scale
  real<lower=0> sigma;           // amplitude
  vector[n_cat] x_delta;         // gp function

  // Epidemic model parameters
  real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
  row_vector<lower = 0>[n_group] gamma;   // rate of recovery
  real<lower = 0, upper = 1> s0[n_group]; // initial fraction of susceptible 
}
  
transformed parameters{
  // Transmission model transformed parameters
  matrix[n_cat, n_cat] cov =  cov_exp_quad(agegroup, agegroup, sigma, rho) + cov_exp_quad(agegroup_sym, agegroup, sigma, rho)  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov = cholesky_decompose(cov);
  vector[n_cat] x_vec = L_cov * x_delta;
  matrix[n_group, n_group] f = to_matrix(x_vec, n_group, n_group, 0);
  matrix[n_group, n_group] c;
  matrix[n_group, n_group] mu_age;

  // Epidemic model transformed parameters
  real y_hat[n_time, n_difeq * n_group]; // solution from the ODE solver
  real y_init[n_difeq * n_group];    // initial conditions for both fractions of S, I and R
  matrix[n_group, n_group] beta;
  real beta_gamma[n_group + n_group*n_group]; 
  
  for(i in 1:n_group){
    for(j in 1:n_group){
      c[i,j] = exp(alpha + f[i,j]);
      mu_age[i,j] = logU[i,j] + c[i,j];
    }
    y_init[1+3*(i-1)] = s0[i];
    y_init[2+3*(i-1)] = 1 - s0[i];
    y_init[3+3*(i-1)] = 0;
  }

  beta = c * diag_matrix(to_vector(T)) * q;
  beta_gamma = to_array_1d(append_row(gamma, beta)');
  y_hat = integrate_ode_bdf(SIR, y_init, t0, ts, beta_gamma, x_r, x_i);
}
  
model{
  real lambda[n_time, n_group];      //poisson parameter
  
  //priors
  // Transmission model priors
  alpha ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  rho ~ inv_gamma(5, 5);
  sigma ~ std_normal();
  x_delta ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(10, 0.1);
  
  // Likelihood
  // Transmission model likelihood
  for(i in 1:n_group){
    for(j in 1:n_group){
      Y[i,j] ~ normal(mu_age[i,j], theta);
    }
  }
  //Epidemic model likelihood
  for(t in 1:n_time){
    for(i in 1:n_group){
      lambda[t,i] = y_hat[t,2+3*(i-1)]*T[i];
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }
}

generated quantities {
  real R_0[n_group, n_group];  // Basic reproduction number
  for(i in 1:n_group){
    for(j in 1:n_group){
      R_0[i,j] = beta[i,j]/gamma[i];   
    }
  }
}
"

## DETERMINISTIC
model.overall.age.deterministic	<-  "
functions {
  real[] SIR(
  real t,           // time
  real[] y_raw,     // system state {susceptible,infected,recovered}
  real[] theta,     // parameters
  real[] x_r,
  int[] x_i) {
  
  int n_age = size(y_raw)/3;
  real dy_dt[3*n_age];
  matrix[3, n_age] y = to_matrix(y_raw, 3, n_age);
  matrix[n_age, n_age] beta = to_matrix(theta[(n_age+1):(n_age+n_age*n_age)], n_age, n_age); 
  real total;

  for(i in 1:n_age){
    total = 0.0;
    for(j in 1:n_age){
      total = total + beta[i,j] * y[2,j];
    }
    dy_dt[1+3*(i-1)] = - total * y[1,i];
    dy_dt[2+3*(i-1)] = total * y[1,i] - theta[i] * y[2,i];
    dy_dt[3+3*(i-1)] = theta[i] * y[2,i];
  }

  return dy_dt;
  }
}

data{
  int<lower = 1> n_time;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int<lower=1> n_age;         // number of age
  int Y[n_group, n_group];    // number of contact by age group
  int Z[n_time, n_age];       // number of infected individuals each day by age
  real age_i[n_age*n_age];    // age group of i
  real age_j[n_age*n_age];    // age group of j
  int corr[n_group,2];        // correspondance between age and age group
  int T[n_age];               // total number of participants by age group 
  int N[n_age];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time];            // time points observed
}
  
transformed data {
  real x_r[0];
  int x_i[0];
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
  // Transmission model parameters
  real<lower=0> theta;           // global dispersion parameter 
  real alpha;                    // global intercept
  real<lower=0> rho;             // length-scale
  real<lower=0> sigma;           // amplitude
  vector[n_cat] x_delta;         // gp function

  // Epidemic model parameters
  real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
  row_vector<lower = 0>[n_age] gamma;   // rate of recovery
  real<lower = 0, upper = 1> s0[n_age]; // initial fraction of susceptible 
}
  
transformed parameters{
  // Transmission model transformed parameters
  matrix[n_cat, n_cat] cov =  cov_exp_quad(age, age, sigma, rho) + cov_exp_quad(age_sym, age, sigma, rho)  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov = cholesky_decompose(cov);
  vector[n_cat] x_vec = L_cov * x_delta;
  matrix[n_age, n_age] f = to_matrix(x_vec, n_age, n_age, 0);
  matrix[n_age, n_age] c;
  matrix[n_age, n_age] mu_age;
  real mu_group[n_group, n_group];
  
  // Epidemic model transformed parameters
  real y_hat[n_time, n_difeq * n_age]; // solution from the ODE solver
  real y_init[n_difeq * n_age];    // initial conditions for both fractions of S, I and R
  matrix[n_age, n_age] beta = c * diag_matrix(to_vector(T)) * q;
  real beta_gamma[n_age + n_age*n_age] = to_array_1d(append_row(gamma, beta)'); 
  
  
  for(i in 1:n_age){
    for(j in 1:n_age){
      c[i,j] = exp(alpha + f[i,j]);
      mu_age[i,j] = U[i,j] + c[i,j];
    }
    y_init[1+3*(i-1)] = s0[i];
    y_init[2+3*(i-1)] = 1 - s0[i];
    y_init[3+3*(i-1)] = 0;
  }
  for(a in 1:n_group){
    for(b in 1:n_group){
      mu_group[a,b] = sum(mu_age[corr[a,1]:corr[a,2], corr[b,1]:corr[b,2]]);
    }
  }
  
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, beta_gamma, x_r, x_i);
}
  
model{
  real lambda[n_time, n_age];      //poisson parameter
  
  //priors
  // Transmission model priors
  alpha ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  rho ~ inv_gamma(5, 5);
  sigma ~ std_normal();
  x_delta ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(10, 0.1);
  
  // Likelihood
  // Transmission model likelihood
  for(a in 1:n_group){
    for(b in 1:n_group){
      Y[a,b] ~ normal(mu_group[a,b], theta);
    }
  }
  //Epidemic model likelihood
  for(t in 1:n_time){
    for(i in 1:n_age){
      lambda[t,i] = y_hat[t,2+3*(i-1)]*T[i];
      Z[t,i] ~ poisson(lambda[t,i]);
    }
  }
}

generated quantities {
  real R_0[n_age, n_age];  // Basic reproduction number
  for(i in 1:n_age){
    for(j in 1:n_age){
      R_0[i,j] = beta[i,j]/gamma[i];   
    }
  }
}
"


model.overall.agegroup.deterministic.h1n1	<-  "
functions {
  real[] SIR(
  real t,           // time
  real[] y_raw,     // system state {susceptible,infected,recovered}
  real[] theta,     // parameters
  real[] x_r,
  int[] x_i) {
  
  int n_group = x_i[1];
  real dy_dt[3*n_group];
  matrix[3, n_group] y = to_matrix(y_raw, 3, n_group);
  matrix[n_group, n_group] beta = to_matrix(theta[(n_group+1):(n_group+n_group*n_group)], n_group, n_group); 
  real total; 

    for(i in 1:n_group){
    total = 0.0;
    for(j in 1:n_group){
      total = total + beta[i,j] * y[2,j];
    }
    dy_dt[1+3*(i-1)] = - total * y[1,i];
    dy_dt[2+3*(i-1)] = total * y[1,i] - theta[i] * y[2,i];
    dy_dt[3+3*(i-1)] = theta[i] * y[2,i];
  }

  return dy_dt;
  }
}

data{
  int<lower = 1> n_time;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int Y[n_group, n_group];    // number of contact by age group
  int Z[n_time, n_group];       // number of infected individuals each day by age
  real agegroup_i[n_group*n_group];    // age group of i
  real agegroup_j[n_group*n_group];    // age group of j
  int T[n_group];               // total number of participants by age group 
  int N[n_group];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time];            // time points observed
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
      logU1[i,j] = log(T[i,1]) + log(N[j,1]);
      logU2[i,j] = log(T[i,2]) + log(N[j,2]);
    }
  }
}

parameters{
  // Transmission model parameters
  real<lower=0> theta[2];           // global dispersion parameter 
  real alpha[2];                    // global intercept
  real<lower=0> rho[2];             // length-scale
  real<lower=0> sigma[2];           // amplitude
  vector[n_cat] x_delta[2];         // gp function

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
  vector[n_cat] x_vec1 = L_cov1 * x_delta[,1];
  matrix[n_group, n_group] f1 = to_matrix(x_vec1, n_group, n_group, 0);
  matrix[n_group, n_group] c1;
  matrix[n_group, n_group] mu_age1;
  
  // Transmission model transformed parameters second part
  matrix[n_cat, n_cat] cov2 =  cov_exp_quad(agegroup, agegroup, sigma[2], rho[2]) + cov_exp_quad(agegroup_sym, agegroup, sigma[2], rho[2])  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov2 = cholesky_decompose(cov2);
  vector[n_cat] x_vec2 = L_cov2 * x_delta[,2];
  matrix[n_group, n_group] f2 = to_matrix(x_vec2, n_group, n_group, 0);
  matrix[n_group, n_group] c2;
  matrix[n_group, n_group] mu_age2;

  // Epidemic model transformed parameters
  real y_hat[n_time, n_difeq * n_group]; // solution from the ODE solver
  real y_init[n_difeq * n_group];    // initial conditions for both fractions of S, I and R
  matrix[n_group, n_group] beta1;
  matrix[n_group, n_group] beta2;
  real beta_gamma[n_group + n_group*n_group]; 
  real beta_gamma2[n_group + n_group*n_group];
  
  for(i in 1:n_group){
    for(j in 1:n_group){
      c1[i,j] = exp(alpha[1] + f1[i,j]);
      mu_age1[i,j] = logU1[i,j] + c1[i,j];
      c2[i,j] = exp(alpha[2] + f2[i,j]);
      mu_age2[i,j] = logU2[i,j] + c2[i,j];
    }
    y_init[1+3*(i-1)] = s0[i]*(1-s_vac[i]);
    y_init[2+3*(i-1)] = 1 - s0[i];
    y_init[3+3*(i-1)] = 0;
  }

  beta1 = c1 * diag_matrix(to_vector(T)) * q;
  beta2 = c2 * diag_matrix(to_vector(T)) * q;
  beta_gamma1 = to_array_1d(append_row(gamma, beta1)');
  beta_gamma2 = to_array_1d(append_row(gamma, beta2)');
  y_hat[1:n_time1,] = integrate_ode_bdf(SIR, y_init, t0, ts[1:n_time1], beta_gamma1, x_r, x_i);
  for(i in 1:n_group){
    y_term[1+3*(i-1)] = y_hat[n_time1,1+3*(i-1)] + s0[i]*s_vac[i];
    y_term[2+3*(i-1)] = y_hat[n_time1,2+3*(i-1)];
    y_term[3+3*(i-1)] = y_hat[n_time1,3+3*(i-1)];
  }
  y_hat[(n_time1+1):(n_time1+n_time2),] = integrate_ode_bdf(SIR, y_term, t0, ts[1:n_time2], beta_gamma2, x_r, x_i);
}
  
model{
  real lambda[n_time, n_group];      //poisson parameter
  
  //priors
  // Transmission model priors
  alpha ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  rho ~ inv_gamma(5, 5);
  sigma ~ std_normal();
  x_delta ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(10, 0.1);
  
  // Likelihood
  // Transmission model likelihood
  for(i in 1:n_group){
    for(j in 1:n_group){
      Y1[i,j] ~ normal(mu_age1[i,j], theta);
      Y2[i,j] ~ normal(mu_age2[i,j], theta);
    }
  }
  //Epidemic model likelihood
  for(t in 1:n_time){
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

model.overall.agegroup.deterministic.euler	<-  "
data{
  int<lower = 1> n_time;      // number of days observed
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower=1> n_group;       // number of age group
  int Y[n_group, n_group];    // number of contact by age group
  int Z[n_time, n_group];       // number of infected individuals each day by age
  real agegroup_i[n_group*n_group];    // age group of i
  real agegroup_j[n_group*n_group];    // age group of j
  int T[n_group];               // total number of participants by age group 
  int N[n_group];               // total number of individuals by age group 
  real t0;                    // initial time point (zero)
  real ts[n_time];            // time points observed
}
  
transformed data {
  real<lower = 0> Tinv[n_group]; 
  int n_cat = n_group*n_group; // number of categories by age
  row_vector[2] agegroup[n_cat];
  row_vector[2] agegroup_sym[n_cat];
  real logU[n_group, n_group];
  for(i in 1:n_cat){
    agegroup[i, 1] = agegroup_i[i];
    agegroup[i, 2] = agegroup_j[i];
    agegroup_sym[i, 1] = agegroup_j[i];
    agegroup_sym[i, 2] = agegroup_i[i];
  }
  for(i in 1:n_group){
  Tinv[i] = T[i]^(-1);
    for(j in 1:n_group){
      logU[i,j] = log(T[i]) + log(N[j]);
    }
  }
}

parameters{
  // Transmission model parameters
  real<lower=0> theta;           // global dispersion parameter 
  real alpha;                    // global intercept
  real<lower=0> rho;             // length-scale
  real<lower=0> sigma;           // amplitude
  vector[n_cat] x_delta;         // gp function

  // Epidemic model parameters
  real<lower = 0, upper = 1> q;         // probability of transmission given a contact 
  real<lower = 0, upper = 1> report[n_group]; // reporting rate
  row_vector<lower = 0>[n_group] gamma;   // rate of recovery
  real<lower = 0, upper = 1> s0[n_group]; // initial fraction of susceptible 
}
  
transformed parameters{
  // Transmission model transformed parameters
  matrix[n_cat, n_cat] cov =  cov_exp_quad(agegroup, agegroup, sigma, rho) + cov_exp_quad(agegroup_sym, agegroup, sigma, rho)  + diag_matrix(rep_vector(1e-10, n_cat));
  matrix[n_cat, n_cat] L_cov = cholesky_decompose(cov);
  vector[n_cat] x_vec = L_cov * x_delta;
  matrix[n_group, n_group] f = to_matrix(x_vec, n_group, n_group, 0);
  matrix[n_group, n_group] c;
  matrix[n_group, n_group] mu_age;

  // Epidemic model transformed parameters
  real S[n_group,n_time]; 
  real I[n_group,n_time]; 
  real R[n_group,n_time]; 
  real<lower=0> newI[n_group,n_time-1];
  real<lower=0> newR[n_group,n_time-1];
  matrix[n_group, n_group] beta;
  
  // Initial 
  for(i in 1:n_group){
    for(j in 1:n_group){
      c[i,j] = exp(alpha + f[i,j]);
      mu_age[i,j] = logU[i,j] + c[i,j];
    }
      S[i,1] = (s0[i])*T[i];
      I[i,1] = (1 - s0[i])*T[i];
      R[i,1] = 0;
  }
    
  // SIR
  beta = c * diag_matrix(to_vector(T)) * q;
   for (t in 2:n_time) {
    for (age in 1:n_group) {
      newI[age,t-1] = fmin(S[age,t-1], S[age,t-1]*sum(to_vector(beta[age,]).*to_vector(I[,t-1]).*to_vector(Tinv)));
      
      S[age,t] = S[age,t-1] - newI[age,t-1];
      I[age,t] = I[age,t-1] + newI[age,t-1];
      
      newR[age,t-1] = fmin(I[age,t], I[age,t]*gamma[age]);
      I[age,t] = I[age,t] - newR[age,t-1];
      R[age,t] = R[age,t-1] + newR[age,t-1];
    }
  }
}
  
model{
  real lambda[n_time, n_group];      //poisson parameter
  
  //priors
  // Transmission model priors
  alpha ~ normal(0,0.001);
  theta ~ lognormal(0,0.001);
  rho ~ inv_gamma(5, 5);
  sigma ~ std_normal();
  x_delta ~ std_normal();
  // Epidemic model priors
  q ~ beta(0.5, 0.5); 
  to_vector(report) ~ beta(0.5, 0.5); 
  to_vector(gamma) ~ lognormal(0,1);
  to_vector(s0) ~ beta(10, 0.1);
  
  // Likelihood
  // Transmission model likelihood
  for(i in 1:n_group){
    for(j in 1:n_group){
      Y[i,j] ~ normal(mu_age[i,j], theta);
    }
  }
  //Epidemic model likelihood
  for(t in 1:n_time){
    for(i in 1:n_group){
      lambda[t,i] = I[i,t]*report[i];
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
