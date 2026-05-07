/**
 * ZINAR(1) with generalized Poission (GP) innovations
 * I use the pmf given in the library HMMpa (dgenpois)
 * an also the paper by Harry Joe, Rong Zhu (2005) Bioemtrica Journal
 * https://onlinelibrary.wiley.com/doi/epdf/10.1002/bimj.200410102 ()
 */

functions {
  real genpoisson_lpmf(int y, real lambda1, real lambda2) {
    if (y < 0) return negative_infinity();
    real lambda = lambda1 + lambda2 * y;
    if (lambda <= 0) return negative_infinity();
    return log(lambda1) + (y - 1) * log(lambda) - lambda - lgamma(y + 1);
  }

  int generalized_poisson_rng(real lambda1, real lambda2) {
    real u = uniform_rng(0, 1);
    real sum = 0;
    //int y;
    for (y in 0:999) {
      real p = exp(genpoisson_lpmf(y | lambda1, lambda2));
      sum += p;
      if (u < sum) {
        return y;
      }
    }
    reject("RNG failed to converge within 1000 iterations.");
    return -1;  // this line is never reached, but required for return type
  }
}

data{
  int<lower=0> T;
  int y[T];
    int<lower=0> ff; // forecast
}
parameters{
  real<lower=0, upper=1> alpha;           // thinning parameter
  real<lower=0> lambda1;    // GP parameter λ1
  real<lower=0,upper=1> lambda2;    //GP parameter λ2
  //real<lower=-1,upper=1> lambda2;    //GP parameter λ2
  real<lower=0, upper=1> rho; // probability of zero inflation
}
transformed parameters{
    vector[T] mu;
    mu[1] = y[1];
    for (t in 2:T){
       int pp=min(y[t-1:t]);
          if(y[t]==0){
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + genpoisson_lpmf(y[t]| lambda1, lambda2)));
           }
           else{
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+genpoisson_lpmf(y[t]| lambda1, lambda2));
            }
        for (j in 1:pp){
             if(y[t]==j){
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + genpoisson_lpmf(y[t]-j| lambda1, lambda2)));
           }
           else{
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+genpoisson_lpmf(y[t]-j| lambda1, lambda2));
            }
        }
    }
}
model{
    alpha  ~ uniform(0,1);
    rho ~   uniform(0,1);
  lambda1  ~ lognormal(0, 2);
  lambda2  ~ uniform(0,1);

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
            y_pred[t] +=generalized_poisson_rng(lambda1, lambda2);
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


