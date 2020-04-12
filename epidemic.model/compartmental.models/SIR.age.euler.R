### SIR WITH EULER METHOD ###
# SUITABLE FOR MULTIPLE ARRAY 

SIR = function(t, n_age, beta, gamma, s0){
  S = matrix(nrow = t+1, ncol = n_age,0); I = S; R = S
  for(i in 1:n_age){
    S[1,i] = s0; I[1,i] = 1 - s0; R[1,i] = 0
  }
  for(time in 1:t){
    for(i in 1:n_age){
      newI = min(S[time,i], sum(beta[i,]*I[time,])*S[time,i])
      S[time+1,i] = S[time,i] - newI
      I[time+1,i] = I[time,i] + newI
      
      newR = min(gamma*I[time,i], I[time,i])
      I[time+1,i] = I[time+1,i] - newR
      R[time+1,i] = R[time,i] + newR 
    }
  }
  SIR.table = data.table(t = rep(1:t, n_age), S = c(S[-1,]), I = c(I[-1,]), R = c(R[-1,]), age = rep(1:n_age, each = t)) %>%
    melt(id.vars = c("t", "age"))
  return(SIR.table)
}
