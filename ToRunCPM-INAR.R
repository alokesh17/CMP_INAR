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

# ── Overdispersed (original)
# Parameters 
# alpha <- 0.6              # INAR(1) thinning parameter
# lambda1 <- 2.5             # GP parameter ?1
# lambda2 <- 0.2              # GP parameter ?2
# n <- 510                 # series length
# rho <- 0.8            # probability of zero inflation



# ── Equidispersed ─────────────────────────────────────────────
#alpha <- 0.6; lambda1 <- 2.5; lambda2 <- 0.0; rho <- 0.05; n <- 510

# ── Underdispersed ────────────────────────────────────────────
alpha <- 0.6; lambda1 <- 2.5; lambda2 <- -0.2; rho <- 0.05; n <- 510

y1<-simul_zinarp(n, alpha, lambda1, lambda2, rho)

# Plot the result
plot(y1, type = "o", col = "darkgreen", main = "Zero-Inflated INAR(1) with GP Innovations", xlab = "Time", ylab = "Value")
hist(y1)

library(rstan, quietly = T)
library(shinystan)

ff<-10 

y<-y1[1:(n-ff)]
T<-length(y)

#setwd("C:/Users/vid09002/Dropbox/ResearchLASSO/COM-Poisson/")

M = 200             # increase if tail is heavy
hybrid_tol = 1e-6    # COMPoissonReg-style switching tolerance

#older stan version

# fitCMP_stan <- stan(file='ZIINAR1-CMP.stan', 
#                  data = list(y=c(y), T=T,ff=ff, M=M, hybrid_tol=hybrid_tol),
#                  thin = 2, chains = 1, iter = 1000, warmup = 100,
#                  seed = 9955)

#optimized stan version
fitCMP_stan <- stan(file='ZIINAR1-CMP-fast.stan', 
                    data = list(y=c(y), T=T,ff=ff, M=M, hybrid_tol=hybrid_tol),
                    thin = 2, chains = 1, iter = 5000, warmup = 1000,
                    seed = 9955)

qoi <- c("lambda", "nu", "alpha","rho")
print(fitCMP_stan, pars=qoi)

fitPred <-summary(fitCMP_stan, pars = "y_pred",  probs = c(0.1, 0.9))$summary
fitPredM<-fitPred[,1]
fitPred05<-fitPred[,4]
fitPred95<-fitPred[,5]
plot(y1[(n-ff+1):n],ylim=c(0,10))
lines(fitPredM)
lines(fitPred05)
lines(fitPred95)


#setwd("C:/Users/vid09002/Dropbox/ResearchLASSO/ZINAR-GeneralizedPoisson/")

fit_stan <- stan(file='ZIINAR1-GP.stan', 
                   data = list(y=c(y), T=T,ff=ff),
                   thin = 10, chains = 1, iter = 20000, warmup = 1000,
                   seed = 9955)
                   
qoi <- c("lambda1", "lambda2", "alpha","rho")
print(fit_stan, pars=qoi)

fitPred <-summary(fit_stan, pars = "y_pred",  probs = c(0.1, 0.9))$summary
fitPredM<-fitPred[,1]
fitPred05<-fitPred[,4]
fitPred95<-fitPred[,5]
plot(y1[(n-ff+1):n],ylim=c(0,20))
lines(fitPredM)
lines(fitPred05)
lines(fitPred95)


# CMP
aic_cmp <- mean(extract(fitCMP_stan)$aic)
bic_cmp <- mean(extract(fitCMP_stan)$bic)

# GP
aic_gp  <- mean(extract(fit_stan)$aic)
bic_gp  <- mean(extract(fit_stan)$bic)

cat(sprintf("CMP: AIC = %.2f, BIC = %.2f\n", aic_cmp, bic_cmp))
cat(sprintf("GP:  AIC = %.2f, BIC = %.2f\n", aic_gp,  bic_gp))

 
