library(prevalence)
library(data.table)

if(1){
  indir = "~/Box\ Sync/2020/social_mixing_ncov"
}

## The methodology behind the betaExpert function is presented by Branscum et al. (2005) and implemented in 
## the BetaBuster software, written by Chun-Lung Su.
## Mean value is assumed to be best
## Expert states with p certainty that true value lies in between lower and upper

### MU ###
mu.parameters = data.table(age = 0:8, alpha = 0, beta = 0)

# 0-10
mu.parameters[1,2] = betaExpert(best = 0.0026*0.002/0.0026/100, lower = 0.0003*0.002/0.0026/100, upper = 0.038*0.002/0.0026/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[1,3] = betaExpert(best = 0.0026*0.002/0.0026/100, lower = 0.0003*0.002/0.0026/100, upper = 0.038*0.002/0.0026/100, p = 0.95, method = "mean")[["beta"]]

# 10-20
mu.parameters[2,2] = betaExpert(best = 0.0148*0.006/0.0148/100, lower = 0.003*0.006/0.0148/100, upper = 0.076*0.006/0.0148/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[2,3] = betaExpert(best = 0.0148*0.006/0.0148/100, lower = 0.003*0.006/0.0148/100, upper = 0.076*0.006/0.0148/100, p = 0.95, method = "mean")[["beta"]]

# 20-30
mu.parameters[3,2] = betaExpert(best = 0.06*0.03/0.06/100, lower = 0.032*0.03/0.06/100, upper = 0.132*0.03/0.06/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[3,3] = betaExpert(best = 0.06*0.03/0.06/100, lower = 0.032*0.03/0.06/100, upper = 0.132*0.03/0.06/100, p = 0.95, method = "mean")[["beta"]]

# 30-40
mu.parameters[4,2] = betaExpert(best = 0.146*0.08/0.146/100, lower = 0.103*0.08/0.146/100, upper = 0.255*0.08/0.146/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[4,3] = betaExpert(best = 0.146*0.08/0.146/100, lower = 0.103*0.08/0.146/100, upper = 0.255*0.08/0.146/100, p = 0.95, method = "mean")[["beta"]]

# 40-50
mu.parameters[5,2] = betaExpert(best = 0.3*0.15/0.3/100, lower = 0.22*0.15/0.3/100, upper = 0.42*0.15/0.3/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[5,3] = betaExpert(best = 0.3*0.15/0.3/100, lower = 0.22*0.15/0.3/100, upper = 0.42*0.15/0.3/100, p = 0.95, method = "mean")[["beta"]]

#50-60
mu.parameters[6,2] = betaExpert(best = 1.3*0.6/1.3/100, lower = 1*0.6/1.3/100, upper = 1.6*0.6/1.3/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[6,3] = betaExpert(best = 1.3*0.6/1.3/100, lower = 1*0.6/1.3/100, upper = 1.6*0.6/1.3/100, p = 0.95, method = "mean")[["beta"]]

#60-70
mu.parameters[7,2] = betaExpert(best = 4*2.2/4/100, lower = 3.4*2.2/4/100, upper = 4.6*2.2/4/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[7,3] = betaExpert(best = 4*2.2/4/100, lower = 3.4*2.2/4/100, upper = 4.6*2.2/4/100, p = 0.95, method = "mean")[["beta"]]

#70-80
mu.parameters[8,2] = betaExpert(best = 8.6*5.1/8.6/100, lower = 7.5*5.1/8.6/100, upper = 10*5.1/8.6/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[8,3] = betaExpert(best = 8.6*5.1/8.6/100, lower = 7.5*5.1/8.6/100, upper = 10*5.1/8.6/100, p = 0.95, method = "mean")[["beta"]]
         
#80-90
mu.parameters[9,2] = betaExpert(best = 13.4*9.3/13.4/100, lower = 11.2*9.3/13.4/100, upper = 15.9*9.3/13.4/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[9,3] = betaExpert(best = 13.4*9.3/13.4/100, lower = 11.2*9.3/13.4/100, upper = 15.9*9.3/13.4/100, p = 0.95, method = "mean")[["beta"]]


### RHO ###
rho.parameters = data.table(age = 0:2, alpha = 0, beta = 0)

# 0-20
rho.parameters[1,2] = betaExpert(best = 0.333, lower = 0.09, upper = 0.906, p = 0.95, method = "mean")[["alpha"]]
rho.parameters[1,3] = betaExpert(best = 0.333, lower = 0.09, upper = 0.906, method = "mean")[["beta"]]

# 20-60
rho.parameters[2,2] = betaExpert(best = 0.702, lower = 0.6, upper = 0.774, p = 0.95, method = "mean")[["alpha"]]
rho.parameters[2,3] = betaExpert(best = 0.702, lower = 0.6, upper = 0.774, p = 0.95, method = "mean")[["beta"]]

# 60-90
rho.parameters[3,2] = betaExpert(best = 0.528, lower = 0.425, upper = 0.581, p = 0.95, method = "mean")[["alpha"]]
rho.parameters[3,3] = betaExpert(best = 0.528, lower = 0.425, upper = 0.581, p = 0.95, method = "mean")[["beta"]]

## save
save(mu.parameters, file = file.path(indir, "results_simu", "mu.parameters.rda"))
save(rho.parameters, file = file.path(indir, "results_simu", "rho.parameters.rda"))


0-10
mu.parameters[1,2] = betaExpert(best = 0.0026/100, lower = 0.0003/100, upper = 0.038/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[1,3] = betaExpert(best = 0.0026/100, lower = 0.0003/100, upper = 0.038/100, p = 0.95, method = "mean")[["beta"]]

# 10-20
mu.parameters[2,2] = betaExpert(best = 0.0148/100, lower = 0.003/100, upper = 0.076/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[2,3] = betaExpert(best = 0.0148/100, lower = 0.003/100, upper = 0.076/100, p = 0.95, method = "mean")[["beta"]]

# 20-30
mu.parameters[3,2] = betaExpert(best = 0.06/100, lower = 0.032/100, upper = 0.132/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[3,3] = betaExpert(best = 0.06/100, lower = 0.032/100, upper = 0.132/100, p = 0.95, method = "mean")[["beta"]]

# 30-40
mu.parameters[4,2] = betaExpert(best = 0.146/100, lower = 0.103/100, upper = 0.255/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[4,3] = betaExpert(best = 0.146/100, lower = 0.103/100, upper = 0.255/100, p = 0.95, method = "mean")[["beta"]]

# 40-50
mu.parameters[5,2] = betaExpert(best = 0.3/100, lower = 0.22/100, upper = 0.42/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[5,3] = betaExpert(best = 0.3/100, lower = 0.22/100, upper = 0.42/100, p = 0.95, method = "mean")[["beta"]]

#50-60
mu.parameters[6,2] = betaExpert(best = 1.3/100, lower = 1/100, upper = 1.6/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[6,3] = betaExpert(best = 1.3/100, lower = 1/100, upper = 1.6/100, p = 0.95, method = "mean")[["beta"]]

#60-70
mu.parameters[7,2] = betaExpert(best = 4/100, lower = 3.4/100, upper = 4.6/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[7,3] = betaExpert(best = 4/100, lower = 3.4/100, upper = 4.6/100, p = 0.95, method = "mean")[["beta"]]

#70-80
mu.parameters[8,2] = betaExpert(best = 8.6/100, lower = 7.5/100, upper = 10/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[8,3] = betaExpert(best = 8.6/100, lower = 7.5/100, upper = 10/100, p = 0.95, method = "mean")[["beta"]]

#80-90
mu.parameters[9,2] = betaExpert(best = 13.4/100, lower = 11.2/100, upper = 15.9/100, p = 0.95, method = "mean")[["alpha"]]
mu.parameters[9,3] = betaExpert(best = 13.4/100, lower = 11.2/100, upper = 15.9/100, p = 0.95, method = "mean")[["beta"]]
