library(tidyverse)
read_csv("https://kingaa.github.io/sbied/mif/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) %>%
  filter(week<=42) -> dat

dat %>%
  ggplot(aes(x=week,y=reports))+
  geom_line()

library(pomp)

sir_step <- function (S, I, R, H, N, Beta, mu_IR, delta.t, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR;
  c(S = S, I = I, R = R, H = H)
}

sir_init <- function (N, eta, ...) {
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)), H = 0)
}

dmeas <- function (reports, H, rho, log, ...) {
  dbinom(x=reports, size=H, prob=rho, log=log)
}

rmeas <- function (H, rho, ...) {
  c(reports=rbinom(n=1, size=H, prob=rho))
}

dat %>%
  pomp(
    times="week",
    t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    accumvars="H",
    rmeasure=rmeas,
    dmeasure=dmeas,
    statenames = c("S", "I", "R", "H"),
    paramnames = c("Beta", "mu_IR", "eta", "rho", "N"),
    partrans=parameter_trans(log=c("Beta","mu_IR"),logit=c("rho","eta"))
  ) -> mearSIR

fixed_params = c(N = 38000)
params = c(Beta = 20, mu_IR =2, rho = .5, eta=.1, N = 38000)

pmcmc(pomp(mf, dprior=Csnippet(priorDens),
           paramnames=c("sigma","phi","r")),
      Nmcmc = 500, Np = 1000,
      proposal = mvn.diag.rw(
        rw.sd=c(N_0=0.1, sigma=0.02, r=0.02, phi=0.02)
      )) -> pmh


params = c(N = 38000)
pf = mearSIR %>%
  pfilter(params = params, Np = 10000)
pf %>%
  coef()

pf %>%
  logLik() %>%
  logmeanexp(se = TRUE)

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  H = 0;
")

dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
")

rmeas <- Csnippet("
  reports = rbinom(H,rho);
")


dat %>%
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    partrans=parameter_trans(log=c("Beta","mu_IR"),logit=c("rho","eta")),
    statenames=c("S","I","H"),
    paramnames=c("Beta","mu_IR","eta","rho","N")
  ) -> measSIR


params <- c(Beta=20,mu_IR=2,rho=0.5,eta=0.1,N=38000)
measSIR %>%
  simulate(params=params,nsim=10,format="data.frame") -> y
y %>%
  ggplot(aes(x=week,y=reports,group=.id))+
  geom_line()
measSIR %>%
  pfilter(Np=1000,params=params) -> pf
plot(pf)


####################################

obs = c(1,5,10,75,225,300,260,240,190,125,60,30,20,0)

rmeas <- "
  cases = rnbinom_mu(theta, rho * H);
"

dmeas <- "
  lik = dnbinom_mu(cases, theta, rho * H, give_log);
"

sir.step <- "
  double rate[6];
  double dN[6];
  double P;
  P = S + I + R;
  rate[0] = mu * P;       // birth
  rate[1] = Beta * I / P; // transmission
  rate[2] = mu;           // death from S
  rate[3] = gamma;        // recovery
  rate[4] = mu;           // death from I
  rate[5] = mu;           // death from R
  dN[0] = rpois(rate[0] * dt);
  reulermultinom(2, S, &rate[1], dt, &dN[1]);
  reulermultinom(2, I, &rate[3], dt, &dN[3]);
  reulermultinom(1, R, &rate[5], dt, &dN[5]);
  S += dN[0] - dN[1] - dN[2];
  I += dN[1] - dN[3] - dN[4];
  R += dN[3] - dN[5];
  H += dN[1];
"
tmp = data.table(Z = obs, t = 1:14)
sir1 <- pomp(
  data = tmp,
  times = "t",
  t0 = -1/7,
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
  rprocess = euler(step.fun = Csnippet(sir.step), delta.t = 1/7),
  obsnames="cases",
  statenames = c("S", "I", "R", "H"),
  paramnames = c("gamma", "mu", "theta", "Beta", "popsize",
                 "rho", "S.0", "I.0", "R.0"),
  accumvars = "H",
  rinit = Csnippet("
    double sum = S_0 + I_0 + R_0;
    S = nearbyint(popsize * S_0 / sum);
    I = nearbyint(popsize * I_0 / sum);
    R = nearbyint(popsize * R_0 / sum);
    H = 0;
    "),
  partrans=parameter_trans(log=c("Beta","mu", "gamma", "theta"),logit=c("rho","eta")),
  params = c(popsize = 763, Beta = 1, gamma = 1,
             mu = 1/50, rho = 0.1, theta = 10, S.0 = 761/763,
             I.0 = 1/763, R.0 = 1/763))

mif2(sir1, Nmif=30, Np=1000,
     cooling.fraction.50=0.8,cooling.type="geometric",
     rw.sd=rw.sd(r=0.02,sigma=0.02,phi=0.02,N_0=ivp(0.1))
) -> mf


plot(sir1)
