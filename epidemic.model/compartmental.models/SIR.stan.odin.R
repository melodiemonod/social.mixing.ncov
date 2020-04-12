### DETERMINISTIC SIR WITH ODIN AND STAN, SUITABLE FOR SINGLE ARRAY.
library(data.table)
library(dplyr)
library(reshape2)

sir.comp = function(theta, S0, t)
{
  ## STAN ###
  library("rstan")
  library("rethinking")
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  mod1_stat <- "
  functions {
    real[] SIR(real t,  // time
    real[] y,           // system state {susceptible,infected,recovered}
    real[] theta,       // parameters
    real[] x_r,
    int[] x_i) {
    
    real dy_dt[3];
    
    dy_dt[1] = - theta[1] * y[1] * y[2];
    dy_dt[2] = theta[1] * y[1] * y[2] - theta[2] * y[2];
    dy_dt[3] = theta[2] * y[2];
    
    return dy_dt;
    }
  }
  
  data {
    int<lower = 1> n_obs;       // number of days observed
    int<lower = 1> n_theta;     // number of model parameters
    int<lower = 1> n_difeq;     // number of differential equations
    real t0;                // initial time point (zero)
    real ts[n_obs];         // time points observed
    real<lower = 0> theta[n_theta];
    real<lower = 0, upper = 1> S0;  // initial fraction of susceptible individuals
  }
    
  transformed data {
    real x_r[0];
    int x_i[0];
  }

  transformed parameters{
    real y_hat[n_obs, n_difeq]; // solution from the ODE solver
    real y_init[n_difeq];     // initial conditions for both fractions of S and I
    
    y_init[1] = S0;
    y_init[2] = 1 - S0;
    y_init[3] = 0;
    y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
  }
  "
  df = list(n_obs = length(t), 
              n_theta = length(theta), 
              n_difeq = 3, 
              t0 = 0, 
              ts = t, 
              theta = theta, 
              S0 = S0)
  
  cat("start ODE solver with STAN \n")
  time.start_stan <- Sys.time()
  fit.STAN = stan(model_code = mod1_stat, data = df, algorithm='Fixed_param', iter = 1, warmup = 0, chain = 1)
  time.end_stan <- Sys.time()
  cat("end ODE solver with STAN \n")
  duration_stan <- time.end_stan - time.start_stan
  cat("duration ODE solver with STAN:", duration_stan, "\n")
  sample = extract.samples(fit.STAN)
  SIR.STAN = data.table(t = t, 
                        S = sample$y_hat[,,1], 
                        I = sample$y_hat[,,2], 
                        R = sample$y_hat[,,3]) %>%
    melt(id.vars = "t")%>%
    mutate(method = "STAN")


  ### ODIN ###
  library(dde)
  library(odin)
  
  sir_generator <- odin::odin({
    ## Core equations for transitions between compartments:
    update(S) <- S - beta * I * S
    update(I) <- I + beta * I * S - gamma * I
    update(R) <- R + gamma * I
    
    ## Initial states:
    initial(S) <- S0
    initial(I) <- 1 - S0
    initial(R) <- 0
    
    ## User defined parameters - default in parentheses:
    S0 <- user(.99)
    beta <- user(.2)
    gamma <- user(.1)
  }, verbose = FALSE)
  cat("start ODE solver with ODIN \n")
  time.start_odin <- Sys.time()
  model.ODIN = sir_generator(S0 = S0, beta = theta[1], gamma = theta[2])
  fit.ODIN = model.ODIN$run(c(0,t))
  time.end_odin <- Sys.time()
  cat("end ODE solver with ODIN \n")
  duration_odin <- time.end_odin - time.start_odin
  cat("duration ODE solver with ODIN:", duration_odin, "\n")
  SIR.ODIN <- as.data.table(fit.ODIN[-1,]) %>%
    melt(id.vars = "step") %>%
    rename(t = step)%>%
    mutate(method = "ODIN")
  
  ### DISCRETE ###
  S = c(); I = c(); R = c()
  S[1] = S0; I[1] = 1- S0; R[1] = 0
  for(time in t){
    S[time + 1] = S[time] - theta[1]*S[time]*I[time]
    I[time + 1] = I[time] + theta[1]*S[time]*I[time] - theta[2]*I[time]
    R[time + 1] = R[time] + theta[2]*I[time]
  }
  SIR.DISCRETE = data.table(t = t, S = S[-1], I = I[-1], R = R[-1]) %>%
    melt(id.vars = "t") %>%
    mutate(method = "Discrete")
  
  ### comparison ###
  SIR = rbind(SIR.ODIN, SIR.STAN, SIR.DISCRETE) %>%
    rename(state = variable)
  
return(list(c(duration_stan, duration_odin), SIR))
}

# define parameters
t = 1:120 # observation period
S0 = .99 # proportion of susceptible at t0, i.e., S(0)/N
theta = c(.2,.1) # beta and gamma

# run models
res.cont = sir.comp(theta, S0, t)

# Plot
ggplot(data = res.cont[[2]], aes(x = t, y = value, col = state, linetype = method)) +
  geom_line() + 
  labs(x = "time", y = "Proportion of the population") +
  theme_minimal() +
  theme(text = element_text(size=18))+ 
  annotate(geom="text", x=95, y=.98, label=paste("Computing time STAN:", round(res.cont[[1]][1], 2), "s")) +
  annotate(geom="text", x=95, y=.93, label=paste("Computing time ODIN:", round(res.cont[[1]][2], 4), "s")) + 
  annotate(geom="text", x=95, y=.88, label = paste("R0 =", round(theta[1]/theta[2], 2)))
  
# Interactive plot
library(shiny)

ui <- fluidPage(
  inputPanel(
    sliderInput('beta', label = 'beta:',
                min = 0, max = 2, value = .2, step = .01),
    
    sliderInput('gamma', label = 'gamma:',
                min = 0, max = 2, value = .1, step = .01),
    
    sliderInput('s0', label = 's0:',
                min = 0, max = 1, value = .99, step = .01)
  ),
  
  mainPanel(plotOutput("outplot"))
)

server = function(input, output) 
  {
  output$outplot <- renderPlot({
    res.cont = sir.comp(c(input$beta, input$gamma), input$s0, t)
    par(mar = c(4, 4, .1, .5))
    ggplot(data = subset(res.cont[[2]], method != "Discrete"), aes(x = t, y = value, col = state, linetype = method)) +
      geom_line(data = subset(res.cont[[2]], method == "Discrete"), aes(x = t, y = value, group = state, linetype = method), col = "black") +
      geom_line() +
      labs(x = "time", y = "Proportion of the population") +
      theme_minimal() +
      theme(text = element_text(size=18))+
      annotate(geom="text", x=95, y=.98, label=paste("Computing time STAN:", round(res.cont[[1]][1], 2), "s")) +
      annotate(geom="text", x=95, y=.93, label=paste("Computing time ODIN:", round(res.cont[[1]][2], 4), "s")) +
      annotate(geom="text", x=95, y=.88, label = paste("R0 =", round(theta[1]/theta[2], 2)))
  })
}

shinyApp(ui = ui, server = server)



## discrete with a shifting parameter
### DISCRETE ###
ui <- fluidPage(
  inputPanel(
    sliderInput('beta', label = 'beta:',
                min = 0, max = 100, value = .2, step = .01),
    
    sliderInput('gamma', label = 'gamma:',
                min = 0, max = 2, value = .1, step = .01),
    
    sliderInput('s0', label = 's0:',
                min = 0, max = 1, value = .99, step = .01),
    
    sliderInput('shift', label = 'shift:',
                min = 0, max = 3, value = .1, step = .01),
    
    sliderInput('newsusceptible', label = 'newsusceptible:',
                min = 0, max = 1, value = .3, step = .01)
  ),
  
  mainPanel(plotOutput("outplot"))
)

server = function(input, output) 
{
  output$outplot <- renderPlot({
    SIR.DISCRETE = SIR(input$beta, input$gamma, input$shift, input$s0, input$newsusceptible)
    
    ggplot(data = SIR.DISCRETE, aes(x = t, y = value, col = state)) +
      geom_line() + 
      labs(x = "time", y = "Proportion of the population") +
      theme_minimal() +
      theme(text = element_text(size=18))
  })
}


shinyApp(ui = ui, server = server)


SIR.DISCRETE = SIR(input$beta, input$gamma, input$shift, input$s0)

ggplot(data = SIR.DISCRETE, aes(x = t, y = value, col = state)) +
  geom_line() + 
  labs(x = "time", y = "Proportion of the population") +
  theme_minimal() +
  theme(text = element_text(size=18))

SIR = function(beta, gamma, shift, s0, newsusceptible){
  t = 1:120 
  S = c(); I = c(); R = c()
  S[1] = s0; I[1] = 1- s0; R[1] = 0
  for(time in t){
    #if(time == 50){
    #  S[time] = S[time] + newsusceptible
    #}
    #if(time >= 50){
    #  S[time + 1] = S[time] - (beta+shift)*S[time]*I[time]
    #  I[time + 1] = I[time] + (beta+shift)*S[time]*I[time] - gamma*I[time]
    #  R[time + 1] = R[time] + gamma*I[time]
    #} else{
      S[time + 1] = S[time] - min(beta*I[time], I[time])*S[time]
      I[time + 1] = I[time] + min(beta*I[time], I[time])*S[time] - gamma*I[time]
      R[time + 1] = R[time] + gamma*I[time] 
   # }
  }
  SIR.DISCRETE = data.table(t = t, S = S[-1], I = I[-1], R = R[-1]) %>%
    melt(id.vars = "t") %>%
    rename(state = variable)
  return(SIR.DISCRETE)
}


ggplot(data = SIR.DISCRETE, aes(x = t, y = value, col = state)) +
  geom_line() + 
  labs(x = "time", y = "Proportion of the population") +
  theme_minimal() +
  theme(text = element_text(size=18))
