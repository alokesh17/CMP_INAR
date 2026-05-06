library(COMPoissonReg)

set.seed(123)

n <- 200
X <- cbind(1, runif(n, -1, 1))
beta <- c(1, -3)
lambda <- as.vector(exp(X %*% beta))
nu <- 2

y <- rcmp(n, lambda, nu)
hist(y)

fit <- glm.cmp(
  formula.lambda = y ~ X - 1,
  formula.nu = ~ 1
)

summary(fit)

### Stan

#setwd("C:/Users/vid09002/Dropbox/ResearchLASSO/COM-Poisson")
M = 200             # increase if tail is heavy
hybrid_tol = 1e-6    # COMPoissonReg-style switching tolerance

library(rstan, quietly = T)
library(shinystan)

data = list(N = n, y = y, x = X[,2], M=M, hybrid_tol=hybrid_tol)

## Normal

fit.n_stan <- stan(file='CMP-regression.stan',
                   data =data ,
                   thin = 5, chains = 1, iter = 10000, warmup = 1000,
                   seed = 9900)

qoi=c("beta0","beta1","nu")
print(fit.n_stan,par=qoi)



