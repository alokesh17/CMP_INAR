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
}


data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  vector[N] x;
  int<lower=1> M;
  real<lower=0> hybrid_tol;
}

parameters {
  real beta0;
  real beta1;
  real<lower=0> nu;
}

model {
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  nu ~ lognormal(0, 1);

  for (i in 1:N) {
    real lambda_i = exp(beta0 + beta1 * x[i]);
    y[i] ~ cmp(lambda_i, nu, M, hybrid_tol);
  }
}
