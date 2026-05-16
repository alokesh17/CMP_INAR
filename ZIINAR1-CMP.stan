/**
 * ZINAR(1) with Conway-Maxwell Poisson (CMP) innovations
 * I use the pmf given in the library COMPoissonReg 
 * an also the paper by Raim and  Sellers (2022) 
 * https://www.census.gov/content/dam/Census/library/working-papers/2022/adrm/RRC2022-01.pdf
 */

functions {

  real cmp_log_Z_asymp(real lambda, real nu) {
    real log_lambda = log(lambda);

    return nu * exp(log_lambda / nu)
           - ((nu - 1.0) / (2.0 * nu)) * log_lambda
           - ((nu - 1.0) / 2.0) * log(2.0 * pi())
           - 0.5 * log(nu);
  }

  real cmp_log_Z_trunc(real lambda, real nu, int M) {
    vector[M + 1] log_terms;
    real log_lambda = log(lambda);

    for (r in 0:M) {
      log_terms[r + 1] = r * log_lambda - nu * lgamma(r + 1);
    }

    return log_sum_exp(log_terms);
  }

  real cmp_log_Z(real lambda, real nu, int M, real hybrid_tol) {
    /*
      COMPoissonReg idea:
      if lambda^(-1/nu) < hybrid_tol, use asymptotic approximation.
      Otherwise use truncated series.
    */

    real test = exp(-log(lambda) / nu);

    if (test < hybrid_tol) {
      return cmp_log_Z_asymp(lambda, nu);
    } else {
      return cmp_log_Z_trunc(lambda, nu, M);
    }
  }

  real cmp_lpmf(int y, real lambda, real nu, int M, real hybrid_tol) {
    if (y < 0) {
      return negative_infinity();
    }

    return y * log(lambda)
           - nu * lgamma(y + 1)
           - cmp_log_Z(lambda, nu, M, hybrid_tol);
  }

  real cmp_lcdf(int q, real lambda, real nu, int M, real hybrid_tol) {
    vector[M + 1] log_terms;
    real log_Z;
    int upper;
    real log_lambda = log(lambda);

    if (q < 0) {
      return negative_infinity();
    }

    upper = min(q, M);

    for (r in 0:M) {
      log_terms[r + 1] = r * log_lambda - nu * lgamma(r + 1);
    }

    log_Z = cmp_log_Z(lambda, nu, M, hybrid_tol);

    return log_sum_exp(log_terms[1:(upper + 1)]) - log_Z;
  }
  
int cmp_rng(real lambda, real nu, int M) {
  vector[M + 1] log_terms;
  real log_Z;
  real u;
  real cdf;
  real log_lambda;

  if (lambda <= 0) reject("lambda must be positive.");
  if (nu <= 0) reject("nu must be positive.");
  if (M < 1) reject("M must be at least 1.");

  log_lambda = log(lambda);

  for (r in 0:M) {
    log_terms[r + 1] = r * log_lambda - nu * lgamma(r + 1);
  }

  log_Z = log_sum_exp(log_terms);
  u = uniform_rng(0, 1);
  cdf = 0;

  for (y in 0:M) {
    cdf += exp(log_terms[y + 1] - log_Z);

    if (u <= cdf) {
      return y;
    }
  }

  reject("CMP RNG failed: truncated CDF did not reach 1. Increase M.");
  return -1;
}
}

data{
  int<lower=0> T;
  int y[T];
  int<lower=1> M;
  real<lower=0> hybrid_tol;
    int<lower=0> ff; // forecast
}
parameters{
  real<lower=0, upper=1> alpha;           // thinning parameter
  real<lower=0> lambda;    // CMP parameter λ
  real<lower=0> nu;    // CPM parameter nu

  real<lower=0, upper=1> rho; // probability of zero inflation
}

transformed parameters{
    vector[T] mu;
    mu[1] = y[1];
    for (t in 2:T){
       int pp=min(y[t-1:t]);
          if(y[t]==0){
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      +  cmp_lpmf(y[t]| lambda, nu, M, hybrid_tol)));
           }
           else{
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+cmp_lpmf(y[t]| lambda, nu, M, hybrid_tol));
            }
        for (j in 1:pp){
             if(y[t]==j){
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + cmp_lpmf(y[t]-j| lambda, nu, M, hybrid_tol)));
           }
           else{
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+cmp_lpmf(y[t]-j| lambda, nu, M, hybrid_tol));
            }
        }
    }
}
model{
    alpha  ~ uniform(0,1);
    rho ~   uniform(0,1);
  //lambda  ~ lognormal(0, 2);
  //nu  ~ lognormal(0, 2);
  lambda ~ student_t(4, 0, 5);
  nu ~ student_t(4, 0, 5);

for (t in 1:T) {
               target += log(mu[t]);
               }
}

generated quantities { // FOR PREDICTION
    // Generate posterior predictives
    int y_pred[ff+1];
    int aa;

    // First P points are known
    y_pred[1] = y[T];

    // Posterior predictive
   for (t in 2:ff+1){
        y_pred[t] = binomial_rng(y_pred[t-1], alpha);
          aa = bernoulli_rng(rho);
          if(aa==1){
            y_pred[t] += 0;
            }
            else{
            y_pred[t] +=cmp_rng(lambda, nu, M);
            }
    }

    vector[T] log_lik;

    for (t in 1:T) {
      log_lik[t] = log(mu[t]);
    }

    real ll = sum(log_lik[2:T]);

    real aic = -2*ll+2*3;
    real bic = -2*ll+3*log(T-1);

}


