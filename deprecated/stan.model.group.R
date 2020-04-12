mod1_stat <- '
functions {
  real[] SIR(real t,
  real[] y0,  // time
  real[] theta,
  real[] x_r,
  int[] x_i) {
  
  real S = y0[1];
  real I = y0[2];
  real R = y0[3];
  
  real dS_dt = - theta[1] * S * I;
  real dI_dt = theta[1] * S * I - theta[2] * I;
  real dR_dt = theta[2] * I;
  
  return {dS_dt, dI_dt, dR_dt};
  }
  
  }
  data {
  int<lower = 1> n_obs;       // number of days observed
  int<lower = 1> n_theta;     // number of model parameters
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower = 1> n_pop;       // population 
  int<lower=1> n_age;
  int y[n_obs];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_obs];         // time points observed
  }
  
  transformed data {
  real x_r[0,0];
  int x_i[0,0];
  int len[n_age] = rep_array(n_obs, n_age);
  real ts2[n_age*n_obs];
  for(i in 1:n_age){
    ts2[((i-1)*n_obs +1):(i*n_obs)] = ts;
  }
  }
  
  parameters {
  real<lower = 0> theta[n_age,n_theta]; // model parameters {beta,gamma}
  real<lower = 0, upper = 1> S0;  // initial fraction of susceptible individuals
  }
  
  transformed parameters{
  matrix[n_difeq,n_obs*n_age] y_hat; // solution from the ODE solver
  real y_init[n_age,n_difeq];     // initial conditions for both fractions of S and I
  
  for(i in 1:n_age){
  y_init[i,1] = S0;
  y_init[i,2] = 1 - S0;
  y_init[i,3] = 0;
  }
  
  y_hat = pmx_integrate_ode_group_rk45(SIR, y_init, 0, len, ts2, theta, rep_array(rep_array(0.0,0),n_obs), rep_array(rep_array(0,0),n_obs));
  
  }
  
  model {
  real lambda[n_obs];      //poisson parameter
  
  //priors
  to_vector(theta[1,]) ~ lognormal(0,1);
  to_vector(theta[2,]) ~ gamma(0.004,0.02);  //Assume mean infectious period = 5 days 
  S0 ~ beta(0.5, 0.5);
  
  //likelihood
  for (i in 1:n_obs){
  lambda[i] = y_hat[i,2]*n_pop;
  }
  y ~ poisson(lambda);
  }

'