SEIIR = function(t, theta.b.intervention, theta.a.intervention, T, numb.inf.t0, intervention){
  S = matrix(nrow = t+1, ncol = n_age,0); E = S; IS = S; IA = S; R = S; Death = S
  newE = S; newI = S; newIS = S; newIA = S; newI = S; newDeath = S
  beta = theta.b.intervention$beta
  rho = theta.b.intervention$rho
  gammaS = theta.b.intervention$gammaS
  gammaA = theta.b.intervention$gammaA
  mu = theta.b.intervention$mu
  a = theta.b.intervention$a
  s0 = function(T) 1-numb.inf.t0/T
  
  for(i in 1:n_age){
    S[1,i] = s0(T[i]); IS[1,i] = (1 - s0(T[i]))*rho[i]; IA[1,i] = (1 - s0(T[i]))*(1-rho[i]); R[1,i] = 0; E[1,i] = 0; Death[1,i] = 0
    newI[1,i] = (1 - s0(T[i])); newIS[1,i] = (1 - s0(T[i]))*rho[i]; newIA[1,i] = (1 - s0(T[i]))*(1-rho[i]); newE[1,i] = 0; newDeath[1,i] = 0
  }
  for(time in 1:t){
    for(i in 1:n_age){
      betaA = subset(beta, part.age.cat == age[i])$betaA
      betaS = subset(beta, part.age.cat == age[i])$betaS
      
      newE[time+1,i] = min(S[time,i], sum(betaA*IA[time,] + betaS*IS[time,])*S[time,i])
      S[time+1,i] = S[time,i] - newE[time+1,i]
      E[time+1,i] = E[time,i] + newE[time+1,i] 
      
      newI[time+1,i] = min(a*E[time,i], E[time,i])
      newIS[time+1,i] = newI[time+1,i]*rho[i] 
      newIA[time+1,i] = newI[time+1,i]*(1-rho[i])
      E[time+1,i] = E[time+1,i] - newIS[time+1,i] - newIA[time+1,i]
      IS[time+1,i] = IS[time,i] + newIS[time+1,i]
      IA[time+1,i] = IA[time,i] + newIA[time+1,i]
      
      moveIS = c()
      moveIS[1] = ifelse((gammaS + mu[i])*IS[time,i] > IS[time,i], gammaS/(gammaS + mu[i])*IS[time,i], gammaS*IS[time,i])
      moveIS[2] = ifelse((gammaS + mu[i])*IS[time,i] > IS[time,i], mu[i]/(gammaS + mu[i])*IS[time,i], mu[i]*IS[time,i])
      
      IS[time+1,i] = IS[time+1,i] - sum(moveIS)
      
      newDeath[time+1,i] = moveIS[2]
      Death[time+1,i] = Death[time,i] + newDeath[time+1,i]
      
      moveIA = min(gammaA*IA[time,i], IA[time,i])
      IA[time+1,i] = IA[time+1,i] - moveIA 
      
      R[time+1,i] = R[time,i] + moveIS[1] + moveIA 
    }
    if(sum(Death[time+1,]*T)>=10 & sum(Death[time,]*T)<10 & intervention != 0){
      beta = theta.a.intervention$beta
      rho = theta.a.intervention$rho
      gammaS = theta.a.intervention$gammaS
      gammaA = theta.a.intervention$gammaA
      mu = theta.a.intervention$mu
      a = theta.a.intervention$a
      time.intervention = time+1
    }
  }
  
  SEIIR.table = data.table(t = rep(1:t, n_age), S = c(S[-1,]), E = c(E[-1,]), IA = c(IA[-1,]), 
                           IS = c(IS[-1,]), R = c(R[-1,]), age = rep(age, each = t), T = rep(T, each = t),
                           death = c(Death[-1,]), 
                           newI = c(newI[-1,]), newE = c(newE[-1,]), newIS = c(newIS[-1,]), newIA = c(newIA[-1,]), 
                           newDeath = c(newDeath[-1,])) %>%
    melt(id.vars = c("t", "age", "T")) %>%
    mutate(value.abs = value*T)
  if(intervention != 0) SEIIR.table$time.i = time.intervention
  return(SEIIR.table)
}


draw.theta.intervention=function(p1,p2){
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
  #gammaA = rlnorm(1, log(0.33029), 0.052135)
  
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
  #mu = c(rep(0.0026/100, 2), rep(0.0148/100, 2),rep(0.06/100, 2),rep(0.146/100, 2),rep(0.3/100, 2),
  #            rep(1.3/100, 2), rep(4/100, 2),rep(8.6/100, 2),rep(13.4/100, 2))
  #mu = c(rep(0.002/100, 2), rep(0.006/100, 2),rep(0.03/100, 2),rep(0.08/100, 2),rep(0.15/100, 2),
  #      rep(0.6/100, 2), rep(2.2/100, 2),rep(5.1/100, 2),rep(9.3/100, 2))
  
  # probability of transmission per contact for asymptotic
  #m.subset = subset(pred.data, part.age.cat == cont.age.cat)
  #qA = R0 * (sum(T/sum(T)*m.subset$m*((1-rho)*(gammaS+mu)+1.5*rho*gammaA)/((gammaS+mu)*gammaA)))^(-1)
  #qS = 1.5*qA
  
  C = pred.data %>%
    reshape2::dcast(part.age.cat ~ cont.age.cat, value.var = "m")
  C = C[,-1]
  p_age = T/sum(T); n = length(T)
  M = C
  for(i in 1:n)
  {
    for(j in 1:n){
      M[i,j] = C[i,j]*p_age[i]/p_age[j]*((1-rho[i])*(gammaS+mu[i])+1.5*rho[i]*gammaA)/((gammaS+mu[i])*gammaA)
    }
  }
  eig = eigen(M)
  qA = R0/max(Re(eig$values))  # reverse engineer beta from the R0 and gamma 
  qS = 1.5*qA
  
  # FROM 14/03 TO 15/03
  posterior.m$betaA = posterior.m$m*qA
  posterior.m$betaS = posterior.m$m*qS
  theta.1 = list(R0=R0,beta = posterior.m, rho = rho, gammaS = gammaS, gammaA = gammaA, mu= mu, a=a)
  
  # FROM 16/03 TO 18/03: contacted reduced by p1 and ederly isolates
  elderly = c("[65,70)", "[70,75)", "[75,80)","[80,85)","[85,90]")
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
  
  return(list(theta.1, theta.2, theta.3, theta.4))
}