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

model.ODIN = sir_generator(S0 = 0.99, beta = 0.2, gamma = 0.1)
fit.ODIN = model.ODIN$run(0:300)



sir_s_a_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  
  ## Individual probabilities of transition:
  p_SI[] <- 1 - exp(- sum(force[i,]) / N[i])
  p_IR <- 1 - exp(-gamma)
  force[,] <- beta[i,j] * I[i]
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SI[] <- rbinom(S[i], p_SI[i])
  n_IR[] <- rbinom(I[i], p_IR)
  
  ## Total population size
  N[] <- S[i] + I[i] + R[i]
  
  ## Initial states:
  initial(S[]) <- S_ini
  initial(I[]) <- I_ini
  initial(R[]) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  I_ini <- user(1)
  beta[,] <- user()
  gamma <- user(0.1)
  
  ## Number of replicates
  nsim <- user(2)
  dim(beta) <- c(nsim,nsim)
  dim(N) <- nsim
  dim(S) <- nsim
  dim(I) <- nsim
  dim(R) <- nsim
  dim(p_SI) <- nsim
  dim(n_SI) <- nsim
  dim(n_IR) <- nsim
  dim(force) <- c(nsim,nsim)
  }, verbose = FALSE)
model.ODIN = sir_s_a_generator()
fit.ODIN = model.ODIN$run(0:300)



sir_s_a_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  
  ## Force of transmission:
  betaI[,] <- beta[i,j] * I[i]
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SI[] <- S[i]*sum(betaI[i,]) 
  n_IR[] <- I[i]*gamma
  
  ## Initial states:
  initial(S[]) <- s0[i]
  initial(I[]) <- 1 - s0[i]
  initial(R[]) <- 0
  
  ## User defined parameters - default in parentheses:
  s0[] <- user()
  beta[,] <- user()
  gamma <- user()
  n_age <- user()
  
  ## Number of replicates
  dim(S) <- n_age
  dim(I) <- n_age
  dim(R) <- n_age
  dim(n_SI) <- n_age
  dim(n_IR) <- n_age
  dim(beta) <- c(n_age,n_age)
  dim(betaI) <- c(n_age,n_age)
  dim(s0) <- n_age
}, verbose = FALSE)
model.ODIN = sir_s_a_generator(beta = beta, gamma = gamma, s0 = s0, n_age = n_age)
fit.ODIN = data.table(model.ODIN$run(0:t))[-1,]
SIR.odin = data.table(t = rep(1:t, n_age), S = unlist(c(select(fit.ODIN,"S[1]","S[2]"))), I = unlist(c(select(fit.ODIN,"I[1]","I[2]"))), 
                       R = unlist(c(select(fit.ODIN,"R[1]","R[2]"))), age = rep(1:n_age, each = t)) %>%
  melt(id.vars = c("t", "age"))
SIR.odin$method = "odin"

seiir_s_a_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  update(S[]) <- S[i] - n_SE[i]
  update(E[]) <- E[i] + n_SE[i] - n_EIS[i] - n_EIA[i]
  update(IS[]) <- IS[i] + n_EIS[i] - n_ISR[i] - n_ISD[i]
  update(IA[]) <- IA[i] + n_EIA[i] - n_IAR[i]
  update(R[]) <- R[i] + n_ISR[i] + n_IAR[i]
  update(death[]) <- death[i] + n_ISD[i]
  
  ## Force of transmission:
  betaS_I[,] <- betaS[i,j] * IS[i]
  betaA_I[,] <- betaA[i,j] * IA[i]
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE[] <- S[i]*(sum(betaS_I[i,]) + sum(betaA_I[i,])) 
  n_EIS[] <- E[i]*a*rho[i] 
  n_EIA[] <- E[i]*a*(1-rho[i])
  n_ISR[] <- IS[i]*gammaS
  n_ISD[] <- IS[i]*mu[i]
  n_IAR[] <- IA[i]*gammaA
  
  ## Initial states:
  initial(S[]) <- s0[i]
  initial(E[]) <- 0
  initial(IS[]) <- (1 - s0[i])*rho[i]
  initial(IA[]) <- (1 - s0[i])*(1-rho[i])
  initial(R[]) <- 0
  initial(death[]) <- 0
  
  ## User defined parameters - default in parentheses:
  s0[] <- user()
  rho[] <- user()
  mu[] <- user()
  a <- user()
  betaS[,] <- user()
  betaA[,] <- user()
  gammaA <- user()
  gammaS <- user()
  n_age <- user()
  
  ## Number of replicates
  dim(S) <- n_age
  dim(E) <- n_age
  dim(IS) <- n_age
  dim(IA) <- n_age
  dim(R) <- n_age
  dim(death) <- n_age
  dim(n_SE) <- n_age
  dim(n_EIS) <- n_age
  dim(n_EIA) <- n_age
  dim(n_ISR) <- n_age
  dim(n_ISD) <- n_age
  dim(n_IAR) <- n_age
  dim(s0) <- n_age
  dim(rho) <- n_age
  dim(mu) <- n_age
  dim(betaS) <- c(n_age,n_age)
  dim(betaA) <- c(n_age,n_age)
  dim(betaS_I) <- c(n_age,n_age)
  dim(betaA_I) <- c(n_age,n_age)
}, verbose = FALSE)

betaA = select(theta.b.intervention$beta, part.age.cat, cont.age.cat, betaA) %>%
  dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
  select(-part.age.cat)
betaS = select(theta.b.intervention$beta, part.age.cat, cont.age.cat, betaS) %>%
  dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
  select(-part.age.cat)
model.ODIN = seiir_s_a_generator(s0 = s0(T), rho = theta.b.intervention$rho, betaA = as.matrix(betaA), betaS = as.matrix(betaS),
                                 gammaS = theta.b.intervention$gammaS, gammaA = theta.b.intervention$gammaA, 
                                 mu = theta.b.intervention$mu, a = theta.b.intervention$a, n_age = n_age)
fit.ODIN = data.table(model.ODIN$run(0:t))[-1,]
SIR.odin = data.table(t = rep(1:t, n_age), S = unlist(c(select(fit.ODIN, contains("S["),-contains("IS[")))), E = unlist(c(select(fit.ODIN,contains("E[")))), 
                      IA = unlist(c(select(fit.ODIN,contains("IA[")))), IS = unlist(c(select(fit.ODIN,contains("IS[")))), R = unlist(c(select(fit.ODIN,contains("R[")))), 
                      age = rep(age, each = t), T = rep(T, each = t)) %>%
  melt(id.vars = c("t", "age", "T")) %>%
  mutate(value.abs = value*T)
SIR.odin$method = "odin"
table = SEIIR(t = t, theta.b.intervention = theta.b.intervention, theta.a.intervention = theta.a.intervention, T = T, numb.inf.t0 = 1, intervention =0)
table$method = "discrete"
df = subset(rbind(table, SIR.odin), variable == "S" | variable == "IS" | variable == "IA" | variable == "R" | variable == "E")
ggplot(data = df, aes(x = t, y = value, col = variable, linetype = method)) +
  geom_line() + 
  labs(x = "time", y = "Proportion of the population") +
  theme_minimal() +
  theme(text = element_text(size=18)) +
  facet_wrap(~age)


