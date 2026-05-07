functions {
  real cmp_log_Z_asymp(real lambda, real nu) {
    real log_lambda = log(lambda);
    return nu * exp(log_lambda / nu)
           - ((nu - 1.0) / (2.0 * nu)) * log_lambda
           - ((nu - 1.0) / 2.0) * log(2.0 * pi())
           - 0.5 * log(nu);
  }

  real cmp_log_Z_trunc(real log_lambda, real nu, int M,
                       vector lgam) {
    // ← accept log_lambda and precomputed lgam directly
    vector[M + 1] log_terms;
    for (r in 0:M) {
      log_terms[r + 1] = r * log_lambda - nu * lgam[r + 1];
    }
    return log_sum_exp(log_terms);
  }

  real cmp_log_Z(real lambda, real nu, int M,
                 real hybrid_tol, vector lgam) {
    real log_lambda = log(lambda);
    real test       = exp(-log_lambda / nu);
    if (test < hybrid_tol) {
      return cmp_log_Z_asymp(lambda, nu);
    } else {
      return cmp_log_Z_trunc(log_lambda, nu, M, lgam);
    }
  }

  real cmp_lpmf(int y, real lambda, real nu, int M,
                real hybrid_tol, vector lgam) {
    if (y < 0) return negative_infinity();
    return y * log(lambda)
           - nu * lgam[y + 1]           // ← precomputed lgamma
           - cmp_log_Z(lambda, nu, M, hybrid_tol, lgam);
  }
}

data {
  int<lower=1> N;
  array[N] int<lower=0> y;
  vector[N] x;
  int<lower=1> M;
  real<lower=0> hybrid_tol;
}

transformed data {
  // ── Precompute log(lambda) per observation ────────────────
  // (x is fixed data — log_x could be precomputed if x > 0)

  // ── Precompute lgamma(0:M+1) once before sampling ─────────
  vector[M + 2] lgam;
  for (k in 0:(M + 1)) {
    lgam[k + 1] = lgamma(k + 1);
  }
}

parameters {
  real beta0;
  real beta1;
  real<lower=0> nu;
}

transformed parameters {
  // ── Precompute lambda_i and log_lambda_i once per gradient ─
  vector[N] log_lambda;
  vector[N] log_Z;

  for (i in 1:N) {
    log_lambda[i] = beta0 + beta1 * x[i];   // = log(lambda_i)
  }

  // ── Precompute log_Z for each unique lambda ────────────────
  // (if many x values repeat, could deduplicate — skip for now)
  for (i in 1:N) {
    real lambda_i = exp(log_lambda[i]);
    real test     = exp(-log_lambda[i] / nu);
    if (test < hybrid_tol) {
      log_Z[i] = cmp_log_Z_asymp(lambda_i, nu);
    } else {
      // pass log_lambda directly — avoids recomputing log()
      vector[M + 1] log_terms;
      for (r in 0:M) {
        log_terms[r + 1] = r * log_lambda[i] - nu * lgam[r + 1];
      }
      log_Z[i] = log_sum_exp(log_terms);
    }
  }
}

model {
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  nu    ~ lognormal(0, 1);

  // ── Vectorized log likelihood ──────────────────────────────
  for (i in 1:N) {
    target += y[i] * log_lambda[i]
              - nu * lgam[y[i] + 1]    // ← table lookup O(1)
              - log_Z[i];              // ← precomputed
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = y[i] * log_lambda[i]
                 - nu * lgam[y[i] + 1]
                 - log_Z[i];
  }
}