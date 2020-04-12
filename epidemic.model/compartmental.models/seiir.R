### SEIIR WITH EULER METHOD AND SUITABLE FOR MULTIPLE ARRAY ###

SEIIR = function(t, theta, T, numb.inf.t0){
  S = matrix(nrow = t+1, ncol = n_age,0); E = S; IS = S; IA = S; R = S; Death = S
  beta = theta$beta
  rho = theta$rho
  gammaS = theta$gammaS
  gammaA = theta$gammaA
  mu = theta$mu
  a = theta$a
  s0 = function(T) 1-numb.inf.t0/T
  
  for(i in 1:n_age){
    S[1,i] = s0(T[i]); IS[1,i] = (1 - s0(T[i]))*rho; IA[1,i] = (1 - s0(T[i]))*(1-rho); R[1,i] = 0; E[1,i] = 0; Death[1,i] = 0
  }
  for(time in 1:t){
    for(i in 1:n_age){
      betaA = subset(beta, part.age.cat == age[i])$betaA
      betaS = subset(beta, part.age.cat == age[i])$betaS
      
      newE = min(S[time,i], sum(betaA*IA[time,] + betaS*IS[time,])*S[time,i])
      S[time+1,i] = S[time,i] - newE
      E[time+1,i] = E[time,i] + newE 
      
      newI = min(a*E[time,i], E[time,i])
      E[time+1,i] = E[time+1,i] - newI
      IS[time+1,i] = IS[time,i] + newI*rho
      IA[time+1,i] = IA[time,i] + newI*(1-rho)
      
      moveIS = c()
      moveIS[1] = ifelse((gammaS + mu[i])*IS[time,i] > IS[time,i], gammaS/(gammaS + mu[i])*IS[time,i], gammaS*IS[time,i])
      moveIS[2] = ifelse((gammaS + mu[i])*IS[time,i] > IS[time,i], mu[i]/(gammaS + mu[i])*IS[time,i], mu[i]*IS[time,i])
      
      IS[time+1,i] = IS[time+1,i] - sum(moveIS)
      Death[time+1,i] = Death[time,i] + moveIS[2]
      
      moveIA = min(gammaA*IA[time,i], IA[time,i])
      IA[time+1,i] = IA[time+1,i] - moveIA 
      
      R[time+1,i] = R[time,i] + moveIS[1] + moveIA 
    }
  }
  SEIIR.table = data.table(t = rep(1:t, n_age), S = c(S[-1,]), E = c(E[-1,]), IA = c(IA[-1,]), 
                           IS = c(IS[-1,]), R = c(R[-1,]), age = rep(age, each = t), T = rep(T, each = t)) %>%
    melt(id.vars = c("t", "age", "T")) %>%
    mutate(value.abs = value*T)
  Death.table = data.table(t = rep(1:t, n_age), death = c(Death[-1,]), age = rep(age, each = t), T = rep(T, each = t)) %>%
    melt(id.vars = c("t", "age", "T")) %>%
    mutate(value.abs = value*T)
  return(list(SEIIR.table, Death.table))
}
