### EPI PARAMETERS ###
draw.theta=function(intervention,  p = NULL, p1 = NULL, p2 = NULL){
  # contact matrix 
  u = round(runif(1, 1, 1000))
  m = paste0("m.V",u)
  posterior.draw.aggregate.subset = subset(posterior.draw.aggregate, part.age.cat != "[90,95)" & part.age.cat != "[95,100]" & 
                                             cont.age.cat != "[90,95)" & cont.age.cat != "[95,100]")
  posterior.m = select(posterior.draw.aggregate.subset, m, part.age.cat, cont.age.cat)
  names(posterior.m)[1] = c("m")
  
  # reproduction number - Li et al. (2020)
  R0 = 2.2
  
  #Transition rate of exposed individuals to the infected class - Riou and Althaus (2020)
  ainv = rlnorm(1,1.612,0.418)
  a = 1/ainv
  
  #probability of having symptoms among infected individuals - Tang et al. (2020)
  rho.0to20 = rbeta(1, rho.parameters$alpha[1], rho.parameters$beta[1])
  rho.20to60 = rbeta(1, rho.parameters$alpha[2], rho.parameters$beta[2])
  rho.60to90 = rbeta(1, rho.parameters$alpha[3], rho.parameters$beta[3])
  rho = c(rep(rho.0to20,4), rep(rho.20to60,8), rep(rho.60to90,6))
  
  # recovery rate of symptomatic infected individuals - Tang et al. (2020)
  gammaS = rgamma(1, 0.33029^2/(0.052135^2), 0.33029/(0.052135^2))
  
  # recovery rate of asymptomatic infected individuals - Tang et al. (2020)
  gammaA = rgamma(1,0.13978^2/(0.034821^2), 0.13978/(0.034821^2))

  # Fatality rate fo symptotic individuals by 5y age group - Verity et al. (2020) & Neil M Ferguson (2020)
  mu.0to10 = rbeta(1, mu.parameters$alpha[1], mu.parameters$beta[1])
  mu.10to20 = rbeta(1, mu.parameters$alpha[2], mu.parameters$beta[2])
  mu.20to30 = rbeta(1, mu.parameters$alpha[3], mu.parameters$beta[3])
  mu.30to40 = rbeta(1, mu.parameters$alpha[4], mu.parameters$beta[4])
  mu.40to50 = rbeta(1, mu.parameters$alpha[5], mu.parameters$beta[5])
  mu.50to60 = rbeta(1, mu.parameters$alpha[6], mu.parameters$beta[6])
  mu.60to70 = rbeta(1, mu.parameters$alpha[7], mu.parameters$beta[7])
  mu.70to80 = rbeta(1, mu.parameters$alpha[8], mu.parameters$beta[8])
  mu.80to90 = rbeta(1, mu.parameters$alpha[9], mu.parameters$beta[9])
  
  mu = c(rep(mu.0to10, 2), rep(mu.10to20, 2),rep(mu.20to30, 2),rep(mu.30to40, 2),rep(mu.40to50, 2),
        rep(mu.50to60, 2), rep(mu.60to70, 2),rep(mu.70to80, 2),rep(mu.80to90, 2))
  
  # probability of transmission per contact for asymptotic
  M = pred.data %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var = "m")
  M = M[,-1]
  p_age = T/sum(T); n = length(T)
  w=T/sum(T); Mw = matrix(nrow=n,ncol=n,rep(w,each=n), byrow=T)
  F = matrix(nrow = n*3, ncol = n*3, 0); V = matrix(nrow = n*3, ncol = n*3, 0)
  F[1:n,(n+1):(n+n)] = 1.5*as.matrix(M)*Mw
  F[1:n,(2*n+1):(2*n+n)] = as.matrix(M)*Mw
  V[1:n, 1:n] = diag(rep(a, n)*w)
  V[(n+1):(n+n), 1:n] = -diag(rho*a*w)
  V[(2*n+1):(2*n+n), 1:n] = -diag((1-rho)*a*w)
  V[(n+1):(n+n), (n+1):(n+n)] = diag((gammaS+mu)*w)
  V[(2*n+1):(2*n+n), (2*n+1):(2*n+n)] = diag((gammaA)*w)
  C = F%*%solve(V)
  eig = eigen(C)
  qA = R0/max(Re(eig$values)) 
  qS = 1.5*qA
  posterior.m$betaA = posterior.m$m*qA
  posterior.m$betaS = posterior.m$m*qS
  theta = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a, qA = qA)
  
  if(intervention == 0){
    theta = theta
  } else if(intervention == 1){
    posterior.m$betaA = posterior.m$m*p*qA
    posterior.m$betaS = posterior.m$m*p*qS
    theta.i = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
    theta = list(theta, theta.i)
  } else if(intervention == 2){
    young = c("[0,5)", "[5,10)", "[10,15)","[15,20)")
    posterior.m[which(posterior.m$part.age.cat %in% young & posterior.m$cont.age.cat %in% young),]$m = 0
    theta.i = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
    theta = list(theta, theta.i)
  } else if(intervention == 3){
    elderly = c("[65,70)", "[70,75)", "[75,80)","[80,85)","[85,90]","[85,90)")
    posterior.m[which(posterior.m$part.age.cat %in% elderly | posterior.m$cont.age.cat %in% elderly),]$m = 0
    posterior.m$betaA = posterior.m$m*qA
    posterior.m$betaS = posterior.m$m*qS
    theta.i = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
    theta = list(theta, theta.i)
  } else if(intervention == 4){
    # FROM 14/03 TO 15/03
    posterior.m$betaA = posterior.m$m*qA
    posterior.m$betaS = posterior.m$m*qS
    theta.1 = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
    # FROM 16/03 TO 18/03: contacted reduced by p1 and ederly isolates
    elderly = c("[65,70)", "[70,75)", "[75,80)","[80,85)","[85,90]", "[85,90)")
    posterior.m[which(posterior.m$part.age.cat %in% elderly | posterior.m$cont.age.cat %in% elderly),]$m = 0
    posterior.m$betaA = posterior.m$m*p1*qA
    posterior.m$betaS = posterior.m$m*p1*qS
    theta.2 = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
    # FROM 19/03 TO 22/03: contacted reduced by p1, ederly isolates and school are closed
    young = c("[0,5)", "[5,10)", "[10,15)","[15,20)")
    posterior.m[which(posterior.m$part.age.cat %in% young & posterior.m$cont.age.cat %in% young),]$m = 0
    posterior.m$betaA = posterior.m$m*p1*qA
    posterior.m$betaS = posterior.m$m*p1*qS
    theta.3 = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
    # FROM 23/03 onward: contacted reduced by p2 and school are closed
    posterior.m$betaA = posterior.m$m*p2*qA
    posterior.m$betaS = posterior.m$m*p2*qS
    theta.4 = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
    theta = list(theta.1, theta.2, theta.3, theta.4)
  }
  return(theta)
}

### MODEL ###
## SEIIR without intervention
SEIIR.with.odin = function(t, theta, T, numb.inf.t0){
  s0 = function(T) 1-numb.inf.t0/T
  
  seiir_s_a_generator <- odin::odin({
    ## Core equations for transitions between compartments:
    update(S[]) <- S[i] - n_SE[i]
    update(E[]) <- E[i] + n_SE[i] - n_EIS[i] - n_EIA[i]
    update(IS[]) <- IS[i] + n_EIS[i] - n_ISR[i] - n_ISD[i]
    update(IA[]) <- IA[i] + n_EIA[i] - n_IAR[i]
    update(RS[]) <- RS[i] + n_ISR[i]
    update(RA[]) <- RA[i] + n_IAR[i]
    update(death[]) <- death[i] + n_ISD[i]
    
    ## Force of transmission:
    betaS_I[,] <- betaS[i,j] * IS[j]
    betaA_I[,] <- betaA[i,j] * IA[j]
    
    ## Compartments:
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
    initial(RS[]) <- 0
    initial(RA[]) <- 0
    initial(death[]) <- 0
    
    ## User defined parameters:
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
    dim(RS) <- n_age
    dim(RA) <- n_age
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
  
  betaA = select(theta$beta, part.age.cat, cont.age.cat, betaA) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
    select(-part.age.cat)
  betaS = select(theta$beta, part.age.cat, cont.age.cat, betaS) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
    select(-part.age.cat)
  
  model.ODIN = seiir_s_a_generator(s0 = c(rep(1,5), s0(T[6:13]), rep(1, 5)), 
                                   rho = theta$rho, betaA = as.matrix(betaA), betaS = as.matrix(betaS),
                                   gammaS = theta$gammaS, gammaA = theta$gammaA, 
                                       mu = theta$mu, a = theta$a, n_age = n_age)
  fit.ODIN = data.table(model.ODIN$run(0:t))[-1,]
  
  newIS = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("IS[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("IS[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("death[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("death[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("RS[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("RS[")))[1:(t-1),x]
  })
  newIS = rbind(rep(0,18),newIS)
  newIA = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("IA[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("IA[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("RA[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("RA[")))[1:(t-1),x]
  })
  newIA = rbind(rep(0,18),newIA)
  newI = newIA + newIS
  newdeath = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("death[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("death[")))[1:(t-1),x]
  })
  newdeath = rbind(rep(0,18),newdeath)
  time.10th.death = min(which(apply(select(fit.ODIN,contains("death["))*matrix(nrow = t, ncol = n_age, T, byrow = T), 1, sum) >= 1021))
  SEIIR.table = data.table(t = rep(1:t, n_age), S = unlist(c(select(fit.ODIN, contains("S["),-contains("IS["),-contains("RS[")))), E = unlist(c(select(fit.ODIN,contains("E[")))), 
                           IA = unlist(c(select(fit.ODIN,contains("IA[")))), IS = unlist(c(select(fit.ODIN,contains("IS[")))), 
                           R = unlist(c(select(fit.ODIN,contains("RS["))))+unlist(c(select(fit.ODIN,contains("RA[")))), 
                           death = unlist(c(select(fit.ODIN,contains("death[")))), 
                           newIS = c(newIS), newIA = c(newIA), newI = c(newI), newdeath = c(newdeath),
                           age = rep(age, each = t), T = rep(T, each = t), time.10th.death = time.10th.death) %>%
    melt(id.vars = c("t", "age", "T", "time.10th.death")) %>%
    mutate(value.abs = value*T)
  
  return(SEIIR.table)
}

## SEIIR for intervention 1, 2 or 3
SEIIR.with.odin.intervention = function(intervention, t, theta, T, numb.inf.t0){
  
  s0 = function(T) 1-numb.inf.t0/T
  
  if(intervention == 0){
    theta.1 = theta
    theta.2 = theta
  } else{
    theta.1 = theta[[1]]
    theta.2 = theta[[2]]
  }
  
  seiir_s_a_generator <- odin::odin({
    ## Core equations for transitions between compartments:
    update(S[]) <- S[i] - n_SE[i]
    update(E[]) <- E[i] + n_SE[i] - n_EIS[i] - n_EIA[i]
    update(IS[]) <- IS[i] + n_EIS[i] - n_ISR[i] - n_ISD[i]
    update(IA[]) <- IA[i] + n_EIA[i] - n_IAR[i]
    update(RS[]) <- RS[i] + n_ISR[i]
    update(RA[]) <- RA[i] + n_IAR[i]
    update(death[]) <- death[i] + n_ISD[i]
    
    ## Force of transmission:
    betaS_I[,] <- betaS[i,j] * IS[j]
    betaA_I[,] <- betaA[i,j] * IA[j]
    
    ## Draws from binomial distributions for numbers changing between
    ## compartments:
    n_SE[] <- S[i]*(sum(betaS_I[i,]) + sum(betaA_I[i,])) 
    n_EIS[] <- E[i]*a*rho[i] 
    n_EIA[] <- E[i]*a*(1-rho[i])
    n_ISR[] <- IS[i]*gammaS
    n_ISD[] <- IS[i]*mu[i]
    n_IAR[] <- IA[i]*gammaA
    
    ## Initial states:
    initial(S[]) <- S0[i]
    initial(E[]) <- E0[i]
    initial(IS[]) <- IS0[i]
    initial(IA[]) <- IA0[i]
    initial(RS[]) <- RS0[i]
    initial(RA[]) <- RA0[i]
    initial(death[]) <- death0[i]
    
    ## User defined parameters - default in parentheses:
    S0[] <- user()
    E0[] <- user()
    IS0[] <- user()
    IA0[] <- user()
    RS0[] <- user()
    RA0[] <- user()
    death0[] <- user()
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
    dim(RS) <- n_age
    dim(RA) <- n_age
    dim(death) <- n_age
    dim(n_SE) <- n_age
    dim(n_EIS) <- n_age
    dim(n_EIA) <- n_age
    dim(n_ISR) <- n_age
    dim(n_ISD) <- n_age
    dim(n_IAR) <- n_age
    dim(S0) <- n_age
    dim(E0) <- n_age
    dim(IS0) <- n_age
    dim(IA0) <- n_age
    dim(RS0) <- n_age
    dim(RA0) <- n_age
    dim(death0) <- n_age
    dim(rho) <- n_age
    dim(mu) <- n_age
    dim(betaS) <- c(n_age,n_age)
    dim(betaA) <- c(n_age,n_age)
    dim(betaS_I) <- c(n_age,n_age)
    dim(betaA_I) <- c(n_age,n_age)
  }, verbose = FALSE)
  
   # beta without intervention
  betaA.1 = select(theta.1$beta, part.age.cat, cont.age.cat, betaA) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
    select(-part.age.cat)
  betaS.1 = select(theta.1$beta, part.age.cat, cont.age.cat, betaS) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
    select(-part.age.cat)
  # beta with intervention
  betaA.2 = select(theta.2$beta, part.age.cat, cont.age.cat, betaA) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
    select(-part.age.cat)
  betaS.2 = select(theta.2$beta, part.age.cat, cont.age.cat, betaS) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
    select(-part.age.cat)
  
  # Start the model and identify where the 10th death happenned
  t0 = 20
  model.ODIN.1 = seiir_s_a_generator(S0 = c(rep(1,5), s0(T[6:13]), rep(1, 5)), E0 = rep(0, n_age), IS0 = c(rep(0,5), (1-s0(T[6:13]))*theta.1$rho[6:13], rep(0, 5)),
                                     IA0 = c(rep(0,5), (1-s0(T[6:13]))*(1-theta.1$rho[6:13]), rep(0, 5)),
                                     RS0 = rep(0, n_age), RA0 = rep(0, n_age), death0 = rep(0, n_age),
                                     rho = theta.1$rho, betaA = as.matrix(betaA.1), betaS = as.matrix(betaS.1),
                                     gammaS = theta.1$gammaS, gammaA = theta.1$gammaA, mu = theta.1$mu, 
                                     a = theta.1$a, n_age = n_age)
  fit.ODIN.1 = data.table(model.ODIN.1$run(0:t0))[-1,]
  # which is march 14?
  death = apply(select(fit.ODIN.1, contains("death"))*matrix(nrow = t0, ncol = n_age, T, byrow=T), 1, sum)
  time.10th.death = min(which(death > 10))
  if(intervention == 1 | intervention == 3){
    # From 14/03 to 15/03
    t1 = time.10th.death + 1
  } else if(intervention == 2){
    # From 14/03 to 18/03
    t1 = time.10th.death + 4
  } else{
    t1 = time.10th.death
  }
  fit.ODIN.1 = fit.ODIN.1[1:t1,]
  S0.1 = select(fit.ODIN.1, contains("S["),-contains("IS["),-contains("RS["))[t1,]
  E.1 = select(fit.ODIN.1, contains("E["))[t1,]
  IS.1 = select(fit.ODIN.1, contains("IS["))[t1,]
  IA.1 = select(fit.ODIN.1, contains("IA["))[t1,]
  RS.1 = select(fit.ODIN.1, contains("RS["))[t1,]
  RA.1 = select(fit.ODIN.1, contains("RA["))[t1,]
  death.1 = select(fit.ODIN.1, contains("death"))[t1,]
  
  t2 = t - t1
  model.ODIN.2 = seiir_s_a_generator(S0 = as.numeric(S0.1), E0 = as.numeric(E.1), IS = as.numeric(IS.1), IA = as.numeric(IA.1), 
                                     RS = as.numeric(RS.1), RA = as.numeric(RA.1), death = as.numeric(death.1),
                                     rho = theta.2$rho, betaA = as.matrix(betaA.2), betaS = as.matrix(betaS.2),
                                     gammaS = theta.2$gammaS, gammaA = theta.2$gammaA, 
                                     mu = theta.2$mu, a = theta.2$a, n_age = n_age)
  fit.ODIN.2 = data.table(model.ODIN.2$run(0:t2))[-1,]
  
  fit.ODIN = rbind(fit.ODIN.1,fit.ODIN.2)
  newIS = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("IS[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("IS[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("death[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("death[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("RS[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("RS[")))[1:(t-1),x]
  })
  # movement begins at t1
  newIS = rbind(rep(0,n_age),newIS) 
  newIA = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("IA[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("IA[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("RA[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("RA[")))[1:(t-1),x]
  })
  newIA = rbind(rep(0,n_age),newIA)
  newI = newIA + newIS
  newdeath = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("death[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("death[")))[1:(t-1),x]
  })
  newdeath = rbind(rep(0,n_age),newdeath)
  SEIIR.table = data.table(t = rep(1:t, n_age), S = unlist(c(select(fit.ODIN, contains("S["),-contains("IS["),-contains("RS[")))), E = unlist(c(select(fit.ODIN,contains("E[")))), 
                           IA = unlist(c(select(fit.ODIN,contains("IA[")))), IS = unlist(c(select(fit.ODIN,contains("IS[")))), 
                           R = unlist(c(select(fit.ODIN,contains("RS["))))+unlist(c(select(fit.ODIN,contains("RA[")))), 
                           death = unlist(c(select(fit.ODIN,contains("death[")))), 
                           newIS = c(newIS), newIA = c(newIA), newI = c(newI), newdeath = c(newdeath),
                           age = rep(age, each = t), T = rep(T, each = t), time.10th.death = time.10th.death) %>%
    melt(id.vars = c("t", "age", "T", "time.10th.death")) %>%
    mutate(value.abs = value*T)
  return(SEIIR.table)
}

## SEIIR for intervention 4
SEIIR.with.odin.intervention4 = function(t,theta,T,numb.inf.t0){
  s0 = function(T) 1-numb.inf.t0/T
  
  theta.1 = theta[[1]]
  theta.2 = theta[[2]]
  theta.3 = theta[[3]]
  theta.4 = theta[[4]]
  
  seiir_s_a_generator <- odin::odin({
    ## Core equations for transitions between compartments:
    update(S[]) <- S[i] - n_SE[i]
    update(E[]) <- E[i] + n_SE[i] - n_EIS[i] - n_EIA[i]
    update(IS[]) <- IS[i] + n_EIS[i] - n_ISR[i] - n_ISD[i]
    update(IA[]) <- IA[i] + n_EIA[i] - n_IAR[i]
    update(RS[]) <- RS[i] + n_ISR[i]
    update(RA[]) <- RA[i] + n_IAR[i]
    update(death[]) <- death[i] + n_ISD[i]
    
    ## Force of transmission:
    betaS_I[,] <- betaS[i,j] * IS[j]
    betaA_I[,] <- betaA[i,j] * IA[j]
    
    ## Draws from binomial distributions for numbers changing between
    ## compartments:
    n_SE[] <- S[i]*(sum(betaS_I[i,]) + sum(betaA_I[i,])) 
    n_EIS[] <- E[i]*a*rho[i] 
    n_EIA[] <- E[i]*a*(1-rho[i])
    n_ISR[] <- IS[i]*gammaS
    n_ISD[] <- IS[i]*mu[i]
    n_IAR[] <- IA[i]*gammaA
    
    ## Initial states:
    initial(S[]) <- S0[i]
    initial(E[]) <- E0[i]
    initial(IS[]) <- IS0[i]
    initial(IA[]) <- IA0[i]
    initial(RS[]) <- RS0[i]
    initial(RA[]) <- RA0[i]
    initial(death[]) <- death0[i]
    
    ## User defined parameters - default in parentheses:
    S0[] <- user()
    E0[] <- user()
    IS0[] <- user()
    IA0[] <- user()
    RS0[] <- user()
    RA0[] <- user()
    death0[] <- user()
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
    dim(RS) <- n_age
    dim(RA) <- n_age
    dim(death) <- n_age
    dim(n_SE) <- n_age
    dim(n_EIS) <- n_age
    dim(n_EIA) <- n_age
    dim(n_ISR) <- n_age
    dim(n_ISD) <- n_age
    dim(n_IAR) <- n_age
    dim(S0) <- n_age
    dim(E0) <- n_age
    dim(IS0) <- n_age
    dim(IA0) <- n_age
    dim(RS0) <- n_age
    dim(RA0) <- n_age
    dim(death0) <- n_age
    dim(rho) <- n_age
    dim(mu) <- n_age
    dim(betaS) <- c(n_age,n_age)
    dim(betaA) <- c(n_age,n_age)
    dim(betaS_I) <- c(n_age,n_age)
    dim(betaA_I) <- c(n_age,n_age)
  }, verbose = FALSE)
  
  betaA.1 = select(theta.1$beta, part.age.cat, cont.age.cat, betaA) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
    select(-part.age.cat)
  betaS.1 = select(theta.1$beta, part.age.cat, cont.age.cat, betaS) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
    select(-part.age.cat)
  betaA.2 = select(theta.2$beta, part.age.cat, cont.age.cat, betaA) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
    select(-part.age.cat)
  betaS.2 = select(theta.2$beta, part.age.cat, cont.age.cat, betaS) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
    select(-part.age.cat)
  betaA.3 = select(theta.3$beta, part.age.cat, cont.age.cat, betaA) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
    select(-part.age.cat)
  betaS.3 = select(theta.3$beta, part.age.cat, cont.age.cat, betaS) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
    select(-part.age.cat)
  betaA.4 = select(theta.4$beta, part.age.cat, cont.age.cat, betaA) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaA") %>%
    select(-part.age.cat)
  betaS.4 = select(theta.4$beta, part.age.cat, cont.age.cat, betaS) %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var="betaS") %>%
    select(-part.age.cat)
  
  # From the start of the model to the 10th death
  t0 = 20
  model.ODIN.1 = seiir_s_a_generator(S0 = c(rep(1,5), s0(T[6:13]), rep(1, 5)), E0 = rep(0, n_age), IS0 = c(rep(0,5), (1-s0(T[6:13]))*theta.1$rho[6:13], rep(0, 5)),
                                     IA0 = c(rep(0,5), (1-s0(T[6:13]))*(1-theta.1$rho[6:13]), rep(0, 5)),
                                     RS0 = rep(0, n_age), RA0 = rep(0, n_age), death0 = rep(0, n_age),
                                     rho = theta.1$rho, betaA = as.matrix(betaA.1), betaS = as.matrix(betaS.1),
                                     gammaS = theta.1$gammaS, gammaA = theta.1$gammaA, mu = theta.1$mu, 
                                     a = theta.1$a, n_age = n_age)
  fit.ODIN.1 = data.table(model.ODIN.1$run(0:t0))[-1,]
  # which is march 14?
  death = apply(select(fit.ODIN.1, contains("death"))*matrix(nrow = t0, ncol = n_age, T, byrow=T), 1, sum)
  time.10th.death = min(which(death > 10))
  # From start to to 15/03
  t1 = time.10th.death + 1
  fit.ODIN.1 = fit.ODIN.1[1:t1,]
  S0.1 = select(fit.ODIN.1, contains("S["),-contains("IS["),-contains("RS["))[t1,]
  E.1 = select(fit.ODIN.1, contains("E["))[t1,]
  IS.1 = select(fit.ODIN.1, contains("IS["))[t1,]
  IA.1 = select(fit.ODIN.1, contains("IA["))[t1,]
  RS.1 = select(fit.ODIN.1, contains("RS["))[t1,]
  RA.1 = select(fit.ODIN.1, contains("RA["))[t1,]
  death.1 = select(fit.ODIN.1, contains("death"))[t1,]
  
  # From 16/03 to 18/03
  t2 = 3
  model.ODIN.2 = seiir_s_a_generator(S0 = as.numeric(S0.1), E0 = as.numeric(E.1), IS = as.numeric(IS.1), IA = as.numeric(IA.1), 
                                     RS = as.numeric(RS.1), RA = as.numeric(RA.1), death = as.numeric(death.1),
                                     rho = theta.2$rho, betaA = as.matrix(betaA.2), betaS = as.matrix(betaS.2),
                                     gammaS = theta.2$gammaS, gammaA = theta.2$gammaA, 
                                     mu = theta.2$mu, a = theta.2$a, n_age = n_age)
  fit.ODIN.2 = data.table(model.ODIN.2$run(0:t2))[-1,]
  S0.2 = select(fit.ODIN.2, contains("S["),-contains("IS["),-contains("RS["))[t2,]
  E.2 = select(fit.ODIN.2, contains("E["))[t2,]
  IS.2 = select(fit.ODIN.2, contains("IS["))[t2,]
  IA.2 = select(fit.ODIN.2, contains("IA["))[t2,]
  RS.2 = select(fit.ODIN.2, contains("RS["))[t2,]
  RA.2 = select(fit.ODIN.2, contains("RA["))[t2,]
  death.2 = select(fit.ODIN.2, contains("death"))[t2,]
  
  # From 19/03 to 22/03
  t3 = 4
  model.ODIN.3 = seiir_s_a_generator(S0 = as.numeric(S0.2), E0 = as.numeric(E.2), IS = as.numeric(IS.2), IA = as.numeric(IA.2), 
                                     RS = as.numeric(RS.2), RA = as.numeric(RA.2), death = as.numeric(death.2),
                                     rho = theta.3$rho, betaA = as.matrix(betaA.3), betaS = as.matrix(betaS.3),
                                     gammaS = theta.3$gammaS, gammaA = theta.3$gammaA, 
                                     mu = theta.3$mu, a = theta.3$a, n_age = n_age)
  fit.ODIN.3 = data.table(model.ODIN.3$run(0:t3))[-1,]
  S0.3 = select(fit.ODIN.3, contains("S["),-contains("IS["),-contains("RS["))[t3,]
  E.3 = select(fit.ODIN.3, contains("E["))[t3,]
  IS.3 = select(fit.ODIN.3, contains("IS["))[t3,]
  IA.3 = select(fit.ODIN.3, contains("IA["))[t3,]
  RS.3 = select(fit.ODIN.3, contains("RS["))[t3,]
  RA.3 = select(fit.ODIN.3, contains("RA["))[t3,]
  death.3 = select(fit.ODIN.3, contains("death"))[t3,]
  
  # From 23/03 onwards
  t4 = t - t1 - t2 - t3
  model.ODIN.4 = seiir_s_a_generator(S0 = as.numeric(S0.3), E0 = as.numeric(E.3), IS = as.numeric(IS.3), IA = as.numeric(IA.3), 
                                     RS = as.numeric(RS.3), RA = as.numeric(RA.3), death = as.numeric(death.3),
                                     rho = theta.4$rho, betaA = as.matrix(betaA.4), betaS = as.matrix(betaS.4),
                                     gammaS = theta.4$gammaS, gammaA = theta.4$gammaA, 
                                     mu = theta.4$mu, a = theta.4$a, n_age = n_age)
  fit.ODIN.4 = data.table(model.ODIN.4$run(0:t4))[-1,]
  fit.ODIN = rbind(fit.ODIN.1,fit.ODIN.2,fit.ODIN.3,fit.ODIN.4)
  newIS = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("IS[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("IS[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("death[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("death[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("RS[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("RS[")))[1:(t-1),x]
  })
  newIS = rbind(rep(0,n_age),newIS)
  newIA = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("IA[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("IA[")))[1:(t-1),x] +
      as.matrix(select(fit.ODIN,contains("RA[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("RA[")))[1:(t-1),x]
  })
  newIA = rbind(rep(0,n_age),newIA)
  newI = newIA + newIS
  newdeath = sapply(1:n_age, function(x){
    as.matrix(select(fit.ODIN,contains("death[")))[2:t,x] - as.matrix(select(fit.ODIN,contains("death[")))[1:(t-1),x]
  })
  newdeath = rbind(rep(0,n_age),newdeath)
  SEIIR.table = data.table(t = rep(1:t, n_age), S = unlist(c(select(fit.ODIN, contains("S["),-contains("IS["),-contains("RS[")))), E = unlist(c(select(fit.ODIN,contains("E[")))), 
                           IA = unlist(c(select(fit.ODIN,contains("IA[")))), IS = unlist(c(select(fit.ODIN,contains("IS[")))), 
                           R = unlist(c(select(fit.ODIN,contains("RS["))))+unlist(c(select(fit.ODIN,contains("RA[")))), 
                           death = unlist(c(select(fit.ODIN,contains("death[")))), 
                           newIS = c(newIS), newIA = c(newIA), newI = c(newI), newdeath = c(newdeath),
                           age = rep(age, each = t), T = rep(T, each = t), time.10th.death = time.10th.death) %>%
    melt(id.vars = c("t", "age", "T", "time.10th.death")) %>%
    mutate(value.abs = value*T)
  return(SEIIR.table)
}

# draw.theta=function(){
#   # contact matrix 
#   u = round(runif(1, 1, 1000))
#   m = paste0("m.V",u)
#   posterior.draw.aggregate.subset = subset(posterior.draw.aggregate, part.age.cat != "[90,95)" & part.age.cat != "[95,100]" & 
#                                              cont.age.cat != "[90,95)" & cont.age.cat != "[95,100]")
#   posterior.m = select(posterior.draw.aggregate.subset, m, part.age.cat, cont.age.cat)
#   names(posterior.m)[1] = c("m")
#   
#   # reproduction number - Li et al. (2020)
#   R0 = 2.4
#   
#   #Transition rate of exposed individuals to the infected class - Lauer et al. (2020)
#   dL = rlnorm(1,1.612,0.418)
#   a = 1-exp(-1/dL)
#   a = 1/dL
#   #probability of having symptoms among infected individuals - Tang et al. (2020)
#   rho.0to20 = rbeta(1, rho.parameters$alpha[1], rho.parameters$beta[1])
#   rho.20to60 = rbeta(1, rho.parameters$alpha[2], rho.parameters$beta[2])
#   rho.60to90 = rbeta(1, rho.parameters$alpha[3], rho.parameters$beta[3])
#   rho = c(rep(rho.0to20,4), rep(rho.20to60,8), rep(rho.60to90,6))
#   
#   # recovery rate of symptomatic infected individuals - Tang et al. (2020)
#   gammaS = rgamma(1, 0.33029^2/(0.052135^2), 0.33029/(0.052135^2))
#   
#   # recovery rate of asymptomatic infected individuals - Tang et al. (2020)
#   gammaA = rgamma(1,0.13978^2/(0.034821^2), 0.13978/(0.034821^2))
#   #gammaA = rlnorm(1, log(0.33029), 0.052135)
#   #dI = runif(1, 3, 6)
#   dI = rnorm(1, 4.5, 1)
#   #gammaA = gammaS= rgamma(1,0.13978^2/(0.034821^2), 0.13978/(0.034821^2))
#   #gammaA = gammaS = 1-exp(-1/dI)
#   gammaA = gammaS = 1/dI
#   # duration between symptoms to death
#   dD = rgamma(1,15^2/(6.9^2), 15/(6.9^2))
#   
#   # Fatality rate fo symptotic individuals by 5y age group - Verity et al. (2020) & Neil M Ferguson (2020)
#   mu.0to10 = rbeta(1, mu.parameters$alpha[1], mu.parameters$beta[1])
#   mu.10to20 = rbeta(1, mu.parameters$alpha[2], mu.parameters$beta[2])
#   mu.20to30 = rbeta(1, mu.parameters$alpha[3], mu.parameters$beta[3])
#   mu.30to40 = rbeta(1, mu.parameters$alpha[4], mu.parameters$beta[4])
#   mu.40to50 = rbeta(1, mu.parameters$alpha[5], mu.parameters$beta[5])
#   mu.50to60 = rbeta(1, mu.parameters$alpha[6], mu.parameters$beta[6])
#   mu.60to70 = rbeta(1, mu.parameters$alpha[7], mu.parameters$beta[7])
#   mu.70to80 = rbeta(1, mu.parameters$alpha[8], mu.parameters$beta[8])
#   mu.80to90 = rbeta(1, mu.parameters$alpha[9], mu.parameters$beta[9])
#   
#   #mu = c(rep(mu.0to10, 2), rep(mu.10to20, 2),rep(mu.20to30, 2),rep(mu.30to40, 2),rep(mu.40to50, 2),
#    #      rep(mu.50to60, 2), rep(mu.60to70, 2),rep(mu.70to80, 2),rep(mu.80to90, 2))
#   #mu = c(rep(0.0026/100, 2), rep(0.0148/100, 2),rep(0.06/100, 2),rep(0.146/100, 2),rep(0.3/100, 2),
#   #            rep(1.3/100, 2), rep(4/100, 2),rep(8.6/100, 2),rep(13.4/100, 2))
#   mu = c(rep(mu.0to10, 2), rep(mu.10to20, 2),rep(mu.20to30, 2),rep(mu.30to40, 2),rep(mu.40to50, 2),
#          rep(mu.50to60, 2), rep(mu.60to70, 2),rep(mu.70to80, 2),rep(mu.80to90, 2))
#   #mu = c(rep(0.0026/100, 2), rep(0.0148/100, 2),rep(0.06/100, 2),rep(0.146/100, 2),rep(0.3/100, 2),
#    #           rep(1.3/100, 2), rep(4/100, 2),rep(8.6/100, 2),rep(13.4/100, 2))
#   #mu = c(rep(0.002/100, 2), rep(0.006/100, 2),rep(0.03/100, 2),rep(0.08/100, 2),rep(0.15/100, 2),
#   #      rep(0.6/100, 2), rep(2.2/100, 2),rep(5.1/100, 2),rep(9.3/100, 2))
#   
#   # probability of transmission per contact for asymptotic
#   m.subset = subset(pred.data, part.age.cat == cont.age.cat)
#   qA = R0 * (sum(T/sum(T)*m.subset$m*((1-rho)*(gammaS+mu)+1.5*rho*gammaA)/((gammaS+mu)*gammaA)))^(-1)
#   #qA = .038
#   qS = 1.5*qA
#   
#   # C = pred.data %>%
#   #   reshape2::dcast(part.age.cat ~ cont.age.cat, value.var = "m")
#   # C = C[,-1]
#   # p_age = T/sum(T); n = length(T)
#   # M = C
#   # for(i in 1:n)
#   # {
#   #   for(j in 1:n){
#   #     M[i,j] = C[i,j]*p_age[i]/p_age[j]*((1-rho[i])*(gammaS+mu[i])+1.5*rho[i]*gammaA)/((gammaS+mu[i])*gammaA)
#   #   }
#   # }
#   # 
#   # eig = eigen(M)
#   # qA = R0/max(Re(eig$values))  # reverse engineer beta from the R0 and gamma
#   # qS = 1.5*qA
# 
#   posterior.m$betaA = posterior.m$m*qA
#   posterior.m$betaS = posterior.m$m*qS
#   theta = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a, qA = qA)
#   
#   return(list(theta))
# }

