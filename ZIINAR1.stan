/**
 * Time Series of Count Data: ZINAR(1) 
 * Log-Likelihood derived from:
 * Aldo Paper JSCS
 */

data{
  int<lower=0> T;
  int y[T];
    int<lower=0> ff; // forecast
}
parameters{
  real<lower=0, upper=1> alpha;           // intercept
  real<lower=0> lambda;    // poisson parameter innovation
  real<lower=0, upper=1> rho;
}
transformed parameters{
    vector[T] mu;
    mu[1] = y[1];
    for (t in 2:T){
       int pp=min(y[t-1:t]);
          if(y[t]==0){
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + poisson_lpmf(y[t] | lambda)));                  
           }
           else{
            mu[t] = exp(binomial_lpmf(0| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+poisson_lpmf(y[t]|lambda));
            }           
        for (j in 1:pp){
             if(y[t]==j){
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(log_sum_exp(bernoulli_lpmf(1 | rho), bernoulli_lpmf(0 | rho)
                      + poisson_lpmf(y[t]-j | lambda)));                  
           }
           else{
            mu[t] += exp(binomial_lpmf(j| y[t-1], alpha))*exp(bernoulli_lpmf(0 | rho)+poisson_lpmf(y[t]-j|lambda));
            }
        }
    } 
}
model{
    alpha  ~ uniform(0,1);
    rho ~   uniform(0,1);
  lambda  ~ lognormal(0, 2); 
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
            y_pred[t] +=poisson_rng(lambda);
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
