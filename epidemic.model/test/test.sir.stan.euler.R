### COMPARE STAN CODE AND EULER METHOD ###

n_time = 20; n_age = 3; 
X = matrix(nrow = n_time, ncol = n_age, rep(cumsum(round(runif(n_time, 0,5))), n_age)) 
beta = matrix(nrow = n_age, ncol = n_age, c(0.1,0.1,0.1,4,5,10,5,10,15), byrow = T)

S = matrix(ncol = n_age, nrow = n_time+1, 0); I = S; R = S
s0 = rep(.9, n_age)
S[1,] = s0; I[1,] = 1 - s0; R[1,] = 0
gamma = rep(0.01,3)
q = 0.01
T = c(100, 500,1000)
S = S*T; I = I*T; R = R*T
c = beta/T
for(time in 1:n_time){
  for(i in 1:n_age){
    S[time + 1,i] = S[time,i] - q * sum(c[i,]*I[time,])*S[time,i]
    I[time + 1,i] = I[time,i] + q*sum(c[i,]*I[time,])*S[time,i] - gamma[i]*I[time,i]
    R[time + 1,i] = R[time,i] + gamma[i]*I[time,i] 
  }
}
SIR.DISCRETE = data.table(t = 1:n_time, S = as.vector(S[-1,]), I = as.vector(I[-1,]), R = as.vector(R[-1,]), age = rep(1:n_age, each = n_time)) %>%
  melt(id.vars = c("t", "age"))

ggplot(SIR.DISCRETE, aes(x = t, y = value, col = variable)) +
  geom_line() +
  facet_wrap(~age)

mod1_stat <- "
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
  matrix[n_age, n_age] beta = to_matrix(theta[(n_age+1):(n_age+n_age*n_age)], n_age, n_age,0); // no 0 at the end
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
  int X[n_time, n_age];           // data, total number of infected individuals each day
  matrix[n_age, n_age] c; // contact matrix
  real t0;                // initial time point (zero)
  real ts[n_time];         // time points observed
  row_vector[n_age] loggamma; // model parameters {beta,gamma}
  real<lower = 0, upper = 1> s0[n_age];  // initial fraction of susceptible individuals
  }
  
transformed data {
  real x_r[0];
  int x_i[0];
  real<lower = 0, upper = 1> q; // probability of transmission given a contact
  q = 0.01;
}
  
  
transformed parameters{
  real y_hat[n_time, n_difeq * n_age]; // solution from the ODE solver
  real y_init[n_difeq * n_age];    // initial conditions for both fractions of S and I
  matrix[n_age, n_age] beta = c * diag_matrix(to_vector(T)) * q;
  real theta[n_age + n_age*n_age] = to_array_1d(append_row(exp(loggamma), beta)');
  
  for(i in 1:n_age){
      y_init[1+3*(i-1)] = s0[i];
      y_init[2+3*(i-1)] = 1 - s0[i];
      y_init[3+3*(i-1)] = 0;
  }
    
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
  }
"

data_list = list(n_time = n_time, 
                 n_theta = 2, 
                 n_difeq = 3, 
                 n_age = n_age,
                 T = rep(1000, 3),
                 c = c,
                 X = X, 
                 t0 = 0, 
                 ts = 1:n_time, 
                 loggamma = log(gamma), 
                 s0 = s0)
iter = 1; warmup = 0
fit = stan(model_code = mod1_stat, data = data_list, chains = 1, iter = iter, warmup = warmup, algorithm = "Fixed_param")
posterior = extract(fit)
S.stan = matrix(ncol = n_age, nrow = n_time, 0); I.stan = S.stan; R.stan = S.stan
for(t in 1:n_time){
  for(i in 1:n_age){
    S.stan[t, i] = posterior$y_hat[,t, 1+3*(i-1)]
    I.stan[t, i] = posterior$y_hat[,t,2+3*(i-1)]
    R.stan[t, i] = posterior$y_hat[,t,3+3*(i-1)]
  }
}

SIR.stan = data.table(t = rep(1:n_time, n_age), S = as.vector(S.stan), I = as.vector(I.stan), 
                          R = as.vector(R.stan), age = rep(1:n_age, each = n_time)) %>%
  melt(id.vars = c("t", "age"))

ggplot(SIR.stan, aes(x = t, y = value, col = variable)) +
  geom_line() +
  facet_wrap(~age)


