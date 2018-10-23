rm(list = ls(all = TRUE))
library(mvtnorm)
library(Formula)
library(rstan)
library(loo)

source('../R/ordinal_splogit.R')
rng_seed <- 20180306
set.seed(rng_seed)

n <- 400
## parameters shared by all models
beta <- c(2, 2)      # regression parameters
eta <- c(0, 3)       # cut points
x1 <- rnorm(n, 0, 1) # covariate
r <- 0.5

## time series ##
t <- seq(-1, 1, length = n)
## GP hyperparameters
sig_sq <- 2           # magnitude parameter
rho_sq <- 0.2         # length-scale parameter
jitter <- 0.0001      # jitter

subjects <- data.frame(id = 1:n, x1 = x1, t = t)

dd <- subjects
dd$Y <- NA
setup <- setup.ordinal.timeseries(Y ~ 1 + x1 | t, data = dd)
#                                     X*beta | location
dd2 = simulate.ordinal.timeseries(setup,
                                  beta = beta,
                                  eta = eta,
                                  a = sig_sq,
                                  b = rho_sq,
                                  r = r,
                                  jitter = jitter,
                                  seed = rng_seed)

dd$Y <- dd2$y
write.table(dd, file = "data_example.txt", row.names = F)

head(dd, 20)
table(dd$Y)

#### Fit splogit model ####
PAR_splogit <- c('beta', 'cuts', 'eta_sq', 'rho_sq', 'r', 'f', 'log_lik')
fit_splogit <- stan(file = paste0('../inst/stan/ordinal_fit_splogit_ts.stan'),
                    seed = rng_seed,
                    data = setup.ordinal.timeseries(Y ~ 1 + x1 | t, data = dd)$data,
                    chains = 2, 
                    iter = 1000, 
                    refresh = 50,
                    cores = 2, 
                    open_progress = FALSE,
                    pars = PAR_splogit)

fit.sum_splogit <- summary(fit_splogit)$summary
save(fit.sum_splogit, file = 'results/splogit_results.rdata')
save(fit_splogit, file = 'results/splogit_results_full.rdata')

print(fit.sum_splogit[1:7, c('mean', '2.5%', "97.5%", "Rhat")])

#### Fit logit model ####
PAR_logit <- c('beta', 'cuts', 'eta_sq', 'rho_sq', 'f', 'log_lik')
fit_logit <- stan(file = paste0('../inst/stan/ordinal_fit_logit_ts.stan'),
                  seed = rng_seed,
                  data = setup.ordinal.timeseries(Y ~ 1 + x1 | t, data = dd)$data,
                  chains = 2, 
                  iter = 1000, 
                  refresh = 50,
                  cores = 2, 
                  open_progress = FALSE,
                  pars = PAR_logit)

fit.sum_logit <- summary(fit_logit)$summary
save(fit.sum_logit, file = 'results/logit_results.rdata')
save(fit_logit, file = 'results/logit_results_full.rdata')

print(fit.sum_logit[1:6, c('mean', '2.5%', "97.5%", "Rhat")])

#### trace plot ####
traceplot(fit_splogit, 
          pars = c('beta[1]', 'beta[2]', 'cuts[2]', 'eta_sq', 'rho_sq', 'r'), 
          inc_warmup = TRUE)

traceplot(fit_logit, 
          pars = c('beta[1]', 'beta[2]', 'cuts[2]', 'eta_sq', 'rho_sq'), 
          inc_warmup = TRUE)


#### Model comparison criteria: DIC and LPML ####
criteria.splogit <- comparison.splogit(fit_splogit, dd, dependent = TRUE)
criteria.logit <- comparison.logit(fit_logit, dd, dependent = TRUE)
tab <- data.frame(DIC = rep(NA, 2), LPML = rep(NA, 2))
rownames(tab) <- c("splogit", "logit")
tab[1,1] <- criteria.splogit$DIC
tab[1,2] <- criteria.splogit$LPML
tab[2,1] <- criteria.logit$DIC
tab[2,2] <- criteria.logit$LPML

print(tab)
