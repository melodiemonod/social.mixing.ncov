#### CHAIN BINOMIAL STOCHASTIC SIR IMPLEMENTED IN STAN AND ODIN ####

## SUITABLE FOR SINGLE ARRAY SIR
## COMPARISON USING INTERACTIVE PLOT WITH SHINY

library(dplyr)
library(data.table)
library(ggplot2)

odin.stochastic = function(theta, s0, t, N)
{
  library(dde)
  library(odin)
  sir_generator <- odin::odin({
    ## Core equations for transitions between compartments:
    update(S) <- S - n_SI
    update(I) <- I + n_SI - n_IR
    update(R) <- R + n_IR
    
    ## Individual probabilities of transition:
    p_SI <- 1 - exp(-beta * I / N) # S to I
    p_IR <- 1 - exp(-gamma) # I to R
    
    ## Draws from binomial distributions for numbers changing between
    ## compartments:
    n_SI <- rbinom(S, p_SI)
    n_IR <- rbinom(I, p_IR)
    
    ## Initial states:
    initial(S) <- N*s0
    initial(I) <- N*(1 - s0)
    initial(R) <- 0
    
    ## User defined parameters - default in parentheses:
    N <- user(100)
    s0 <- user(.99)
    beta <- user(0.2)
    gamma <- user(0.1)
  }, verbose = FALSE)
  cat("start ODE solver with ODIN \n")
  time.start_odin <- Sys.time()
  tmp = NULL
  res = lapply(1:1000, function(x){
    model.ODIN = sir_generator(s0 = s0, beta = theta[1], gamma = theta[2], N = N)
    fit.ODIN = model.ODIN$run(0:t)
    as.data.table(fit.ODIN[-1,]) %>%
      melt(id.vars = "step") %>%
      rename(t = step)%>%
      mutate(method = "ODIN", 
             value = value/N)
  })
  time.end_odin <- Sys.time()
  cat("end ODE solver with ODIN \n")
  duration_odin <- time.end_odin - time.start_odin
  cat("duration ODE solver with ODIN:", duration_odin, "\n")
  SIR.ODIN = do.call(rbind, res) 
  SIR.ODIN = SIR.ODIN %>%
    group_by(t, variable) %>%
    summarise(mean = mean(value), l95 = quantile(value, prob = 0.025), u95 = quantile(value, prob = 0.975)) %>%
    rename(state = variable)
  return(list(duration_odin, SIR.ODIN))
}

stan.stochastic = function(theta, s0, t, N)
{
  library("rstan")
  library("rethinking")
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  mod2_stat <- "
  functions {
    real[] SIR(real t,  // time
    real[] y,           // system state {susceptible,infected,recovered}
    real[] theta,       // theta
    real[] x_r,
    int[] x_i) {
    
    real dy_dt[3];
    
    dy_dt[1] = - theta[1] * y[1] * y[2] / sum(y);
    dy_dt[2] = theta[1] * y[1] * y[2] / sum(y) - theta[2] * y[2];
    dy_dt[3] = theta[2] * y[2];
    
    return dy_dt;
    }
  }
  
  data {
    int<lower = 1> n_obs;       // number of days observed
    int<lower = 1> n_theta;     // number of model parameters
    int<lower = 1> n_pop;     // number of individual
    int<lower = 1> n_difeq;     // number of differential equations
    real t0;                // initial time point (zero)
    real ts[n_obs];         // time points observed
    real<lower = 0> theta[n_theta];
    real<lower = 0, upper = 1> s0;  // initial fraction of susceptible individuals
  }
    
  transformed data {
    real x_r[0];
    int x_i[0];
  }
  
  transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];     // initial conditions for both fractions of S and I
    
    y_init[1] = s0;
    y_init[2] = (1 - s0);
    y_init[3] = 0;
    y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
  }
  
  generated quantities {
    real y[n_obs+1, n_difeq];
    y[1,1] = s0*n_pop;
    y[1,2] = (1 - s0)*n_pop;
    y[1,3] = 0;
    
    y[2,2] = normal_rng(y_init[1]*n_pop*(1- exp(- theta[1] * y_init[2])), sqrt(y_init[1]*n_pop*(1- exp(- theta[1] * y_init[2]))*exp(- theta[1] * y_init[2])));
    y[2,3] = normal_rng(y_init[2]*n_pop*(1-exp(- theta[2])), sqrt(y_init[2]*n_pop*(1-exp(- theta[2]))*exp(- theta[2])));
    y[2,1] = n_pop - y[2,2] - y[2,3];
      
    for (i in 3:(n_obs+1)){
      y[i,1] = normal_rng(y_hat[i-1,1]*n_pop*exp(- theta[1] * y_hat[i-1,2]), sqrt(y_hat[i-1,1]*n_pop*(1- exp(- theta[1] * y_hat[i-1,2]))*exp(- theta[1] * y_hat[i-1,2])));
      y[i,3] = y_hat[i-1,3]*n_pop + normal_rng(y_hat[i-1,2]*n_pop*(1-exp(- theta[2])), sqrt(y_hat[i-1,2]*n_pop*(1-exp(- theta[2]))*exp(- theta[2])));
      y[i,2] = n_pop - y[i,1] - y[i,3];
    }
  }
  "
  df = list(n_obs = t, 
            n_theta = length(theta), 
            n_pop = N,
            n_difeq = 3, 
            t0 = 0, 
            ts = 1:t, 
            theta = theta, 
            s0 = s0)
  
  cat("start ODE solver with STAN \n")
  time.start_stan <- Sys.time()
  fit.STAN = stan(model_code = mod2_stat, data = df, algorithm='Fixed_param', iter = 1000, warmup = 0, chain = 1)
  time.end_stan <- Sys.time()
  cat("end ODE solver with STAN \n")
  duration_stan <- time.end_stan - time.start_stan
  cat("duration ODE solver with STAN:", duration_stan, "\n")
  sample = extract.samples(fit.STAN)
  
  S.STAN = data.table(t(sample$y[,-1,1])/N) %>%
    mutate(t = 1:t) %>%
    reshape2::melt(id.vars = c("t")) %>%
    group_by(t) %>%
    summarise(mean = mean(as.numeric(value)), l95 = quantile(as.numeric(value), prob = 0.025), 
              u95 = quantile(as.numeric(value), prob = .975), state = "S")
  I.STAN = data.table(t(sample$y[,-1,2])/N) %>%
    mutate(t = 1:t) %>%
    reshape2::melt(id.vars = c("t")) %>%
    group_by(t) %>%
    summarise(mean = mean(as.numeric(value)), l95 = quantile(as.numeric(value), prob = 0.025), 
              u95 = quantile(as.numeric(value), prob = .975), state = "I")
  R.STAN = data.table(t(sample$y[,-1,3])/N) %>%
    mutate(t = 1:t) %>%
    reshape2::melt(id.vars = c("t")) %>%
    group_by(t) %>%
    summarise(mean = mean(as.numeric(value)), l95 = quantile(as.numeric(value), prob = 0.025), 
              u95 = quantile(as.numeric(value), prob = .975), state = "R")
  SIR.STAN = rbind(S.STAN, I.STAN, R.STAN)
  SIR.STAN$state = factor(SIR.STAN$state, levels = c("S", "I", "R"))
return(list(duration_stan, SIR.STAN))
}

library(shiny)
ui <- fluidPage(
  inputPanel(
    sliderInput('beta', label = 'beta:',
                min = 0, max = 2, value = .2, step = .01),
    
    sliderInput('gamma', label = 'gamma:',
                min = 0, max = 2, value = .1, step = .01),
    
    sliderInput('s0', label = 's0:',
                min = 0, max = 1, value = .99, step = .01),
    
    sliderInput('N', label = 'N:',
                min = 1, max = 10000, value = 1000, step = 10),
    
    sliderInput('t', label = 't:',
                min = 20, max = 200, value = 120, step = 10)
  ),
  
  mainPanel(plotOutput("outplot"))
)


server = function(input, output) 
{
  output$outplot <- renderPlot({
    t = input$t
    res.odin = odin.stochastic(c(input$beta, input$gamma), input$s0, t, input$N)
    p.odin = ggplot(res.odin[[2]], aes(x = t, ymin = l95, ymax = u95, fill = state)) +
      geom_ribbon(alpha = .4) +
      geom_line(aes(y = mean, col = state))+
      labs(x = "time", y = "% of the population") +
      theme_minimal() +
      theme(text = element_text(size=18)) +
      annotate(geom="text", x=95, y=.93, label=paste("Computing time ODIN:", round(res.odin[[1]], 4), "s")) 
    
    res.stan = stan.stochastic(c(input$beta, input$gamma), input$s0, t, input$N)
    p.stan = ggplot(res.stan[[2]], aes(x = t, ymin = l95, ymax = u95, fill = state)) +
      geom_ribbon(alpha = .4) +
      geom_line(aes(y = mean, col = state))+
      labs(x = "time", y = "% of the population") +
      theme_minimal() +
      theme(text = element_text(size=18)) +
      annotate(geom="text", x=95, y=.93, label=paste("Computing time STAN:", round(res.stan[[1]], 4), "s")) 
    par(mar = c(4, 4, .1, .5))
    gridExtra::grid.arrange(p.odin, p.stan, ncol = 1)
  })
}

shinyApp(ui = ui, server = server)




  
