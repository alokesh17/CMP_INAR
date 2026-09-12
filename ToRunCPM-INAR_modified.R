################################################################################
### Generating samples from ZI-INAR(1) with Generalized Poisson Inovations
################################################################################

library(HMMpa)

# Thinning operator: Auxiliary function

thin_operator <- function(x, alpha) {
  sum(rbinom(x, size = 1, prob = alpha))
}

simul_zinarp<-function(n, alpha, lambda, theta, zi_prob){

# Initialize series
zinp_inar <- integer(n)
zinp_inar[1] <- rgenpois(1, lambda, theta)

# Simulation loop
for (t in 2:n) {
  thinned <- thin_operator(zinp_inar[t - 1], alpha)
  
  # Generate zero-inflated innovation
  if (runif(1) < zi_prob) {
    innovation <- 0
  } else {
    innovation <- rgenpois(1, lambda, theta)
  }
  
  zinp_inar[t] <- thinned + innovation
}
zinp_inar
}

################################################################################
### Generating samples from ZI-INAR(1) with CMP distribution
################################################################################
library(COMPoissonReg)

simul_zinarCMP<-function(n, alpha, lambda, nu, zi_prob){
  
  # Initialize series
  zinp_inar <- integer(n)
  zinp_inar[1] <- rcmp(1, lambda, nu)
  
  # Simulation loop
  for (t in 2:n) {
    thinned <- thin_operator(zinp_inar[t - 1], alpha)
    
    # Generate zero-inflated innovation
    if (runif(1) < zi_prob) {
      innovation <- 0
    } else {
      innovation <- rcmp(1, lambda, nu)
    }
    
    zinp_inar[t] <- thinned + innovation
  }
  zinp_inar
}
######################################################################


# Parameters
alpha <- 0.6              # INAR(1) thinning parameter
lambda <- 2.5#1.5             # CPM lambda
nu <- 1.2#0.4              # CMP nu
n <- 500                 # series length
rho <- 0.8            # probability of zero inflation

#y1<-simul_zinarp(n, alpha, lambda1, lambda2, rho)
y1<-simul_zinarCMP(n, alpha, lambda, nu, rho)

# Plot the result
plot(y1, type = "o", col = "darkgreen", main = "Zero-Inflated INAR(1) with GP Innovations", xlab = "Time", ylab = "Value")
hist(y1)

library(rstan, quietly = T)
library(shinystan)

ff<-10 

y<-y1[1:(n-ff)]
T<-length(y)

#setwd("C:/Users/vid09002/Dropbox/ResearchLASSO/COM-Poisson/")

M = 300             # increase if tail is heavy
hybrid_tol = 1e-6    # COMPoissonReg-style switching tolerance

fitCMP_stan <- stan(file='ZIINAR1-CMP.stan', 
                 data = list(y=c(y), T=T,ff=ff, M=M, hybrid_tol=hybrid_tol),
                 thin = 2, chains = 1, iter = 1000, warmup = 100,
                 seed = 9955)

qoi <- c("lambda", "nu", "alpha","rho", "aic","bic")
print(fitCMP_stan, pars=qoi)




fitPred <-summary(fitCMP_stan, pars = "y_pred",  probs = c(0.1, 0.9))$summary
fitPredM<-fitPred[,1]
fitPred05<-fitPred[,4]
fitPred95<-fitPred[,5]
plot(y1[(n-ff+1):n],ylim=c(0,10))
lines(fitPredM)
lines(fitPred05)
lines(fitPred95)

########################################################################
# faster version (Alokesh)
########################################################################
fitCMP_stan_fast <- stan(file='ZIINAR1-CMP-fast.stan', 
                    data = list(y=c(y), T=T,ff=ff, M=M, hybrid_tol=hybrid_tol),
                    thin = 5, chains = 2, iter = 10000, warmup = 1000,
                    seed = 9955)

qoi <- c("lambda", "nu", "alpha","rho", "bic", "aic")
print(fitCMP_stan_fast, pars=qoi)



# Basic summary — shows Rhat (should be <1.01) and n_eff for each parameter
print(fitCMP_stan_fast)

# Rhat specifically — values close to 1 mean chains converged
rhat_vals <- summary(fitCMP_stan_fast)$summary[, "Rhat"]
#print(rhat_vals)

# Trace plots — visually check both chains mix well
traceplot(fitCMP_stan_fast, pars = c("lambda", "nu", "alpha","rho"))


library(bayesplot)

# Extract posterior draws
posterior <- as.array(fitCMP_stan_fast)

# Density overlay — both chains overlaid, great for convergence check
mcmc_dens_overlay(posterior, pars = c("lambda", "nu", "alpha", "rho"))

# Plain density (chains combined)
mcmc_dens(posterior, pars = c("lambda", "nu", "alpha", "rho"))

# Histogram
mcmc_hist(posterior, pars = c("lambda", "nu", "alpha", "rho"))

# Pairs plot — joint distributions, useful for spotting correlations
mcmc_pairs(posterior, pars = c("lambda", "nu", "alpha", "rho"))


# Extract all draws (both chains combined)
draws <- extract(fitCMP_stan_fast)

# e.g. posterior mean and 95% CI for lambda
mean(draws$lambda)
quantile(draws$lambda, c(0.025, 0.975))

# n_eff should be reasonably large (say >100 per parameter)
# With 2 chains, iter=10000, warmup=1000, thin=5:
# usable draws per chain = (10000-1000)/5 = 1800
# total = 3600 draws across both chains

summary(fitCMP_stan_fast)$summary[, c("mean", "sd", "2.5%", "97.5%", "Rhat", "n_eff")]

fitPred_fast <-summary(fitCMP_stan_fast, pars = "y_pred",  probs = c(0.1, 0.9))$summary
fitPredM_fast<-fitPred_fast[,1]
fitPred05_fast<-fitPred_fast[,4]
fitPred95_fast<-fitPred_fast[,5]
plot(y1[(n-ff+1):n],ylim=c(0,10))
lines(fitPredM_fast)
lines(fitPred05_fast)
lines(fitPred95_fast)


#setwd("C:/Users/vid09002/Dropbox/ResearchLASSO/ZINAR-GeneralizedPoisson/")

fit_stanGP <- stan(file='ZIINAR1-GP.stan', 
                   data = list(y=c(y), T=T,ff=ff),
                   thin = 10, chains = 1, iter = 20000, warmup = 1000,
                   seed = 9955)
                   
qoi <- c("lambda1", "lambda2", "alpha","rho","bic", "aic")
print(fit_stanGP, pars=qoi)

fitPred <-summary(fit_stanGP, pars = "y_pred",  probs = c(0.1, 0.9))$summary
fitPredM<-fitPred[,1]
fitPred05<-fitPred[,4]
fitPred95<-fitPred[,5]
plot(y1[(n-ff+1):n],ylim=c(0,20))
lines(fitPredM)
lines(fitPred05)
lines(fitPred95)
 
 
