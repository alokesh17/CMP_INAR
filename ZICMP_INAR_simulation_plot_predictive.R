############################################################################
# ZICMP-INAR(1): ONE script for everything --
#   PART 1 -- simulation: parameter-recovery grid, 5 dispersion regimes
#             (original nu in {0.5,1.0,1.5} + 2 new extremes nu in {0.2,2.5})
#   PART 2 -- plots: recovery, bias, model-comparison win rate (base R)
#   PART 3 -- predictive analysis: qoi + y_pred forecast plot, this
#             model's real parameters (alpha, lambda, nu, rho), not the
#             "lambda1, lambda2" naming from a different model
#   PART 4 -- model comparison vs. ZIP/ZINB/ZIGP (via ZIHINAR1), same 5
#             regimes CROSSED WITH the same n-grid as Part 1 (100,200,
#             400,600) -- EAIC/EBIC/DIC/WAIC1/WAIC2 + win rates per cell
#
# Built around ZIINAR1-CMP-fast-reparam.stan (parameters: alpha, mu_cmp,
# nu, rho; lambda := mu_cmp^nu recovered in transformed parameters;
# generated quantities: y_pred, log_lik, ll, aic, bic).
#
# Nothing fits or plots by itself just from sourcing this file: Part 1's
# grid and Part 4's comparison loop both run once you set RUN_MODE and
# execute; Part 2 only plots once Part 1 has populated `summary_tab`;
# Part 3's worked example is gated with `if (FALSE)`.
#
# INSTALL (once, for Part 4 only): install.packages("ZIHINAR1")
############################################################################

library(rstan)
library(COMPoissonReg)   # for rcmp()
library(ZIHINAR1)        # only needed for Part 4 (ZIP/ZINB/ZIGP comparison)
library(parallel)        # base package, no install needed
rstan_options(auto_write = TRUE)

# ---------------------------------------------------------------------
# PARALLELISM: reps within a cell are independent, so they're run
# concurrently across cores via parallel::mclapply (fork-based) --
# Mac/Linux only, since Windows has no fork(). Each individual
# sampling() call is forced to run its own chains SEQUENTIALLY
# (mc.cores=1 below) so cores aren't double-booked: parallelism happens
# ACROSS reps instead, which uses the machine far more fully than the
# old setup (which only ever kept 2 cores busy, one fit at a time).
# ---------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")   # gives mclapply's forked workers statistically
                            # independent RNG streams (R's own recommended
                            # setup for parallel simulation -- see ?mclapply)
set.seed(2026)

HAS_FORK  <- .Platform$OS.type == "unix"   # mclapply needs Mac/Linux (fork())
N_WORKERS <- if (HAS_FORK) max(1, parallel::detectCores() - 1) else 1
if (!HAS_FORK) {
  message("Note: fork-based parallel reps (parallel::mclapply) need Mac/Linux; ",
          "Windows detected here, so reps will run one at a time below. For ",
          "real speedup on Windows, run this script inside WSL2.")
} else {
  message("Parallel reps enabled: ", N_WORKERS, " worker(s) (", parallel::detectCores(),
          " cores detected, 1 held back for the OS).")
}
options(mc.cores = 1)   # forces each sampling() call's own chains to run
                         # sequentially -- see comment block above

STAN_FILE <- "ZIINAR1-CMP-fast-reparam.stan"   # samples (alpha, mu_cmp, nu, rho),
                                                # lambda := mu_cmp^nu in
                                                # transformed parameters

RUN_MODE <- "full_r10" # <- "quick" | "pilot" | "full_r10" | "full"  (drives Parts 1 and 4)
                        #    "full_r10": every (nu,n) cell, R=10 reps (requested tradeoff)
                        #    "full":     every (nu,n) cell, R=30 reps (original, most expensive)


############################################################################
# PART 1 -- SIMULATION GRID (parameter recovery)
############################################################################

# ---------------------------------------------------------------------
# 1.0 Config
# ---------------------------------------------------------------------
M           <- 300
HYBRID_TOL  <- 1e-6
FF          <- 0        # no forecast needed for the parameter-recovery grid;
                         # PART 3 below uses its own ff > 0 fit
ALPHA_TRUE  <- 0.30
RHO_TRUE    <- 0.30
TARGET_MU   <- 3.0      # common CMP innovation mean held across regimes

if (RUN_MODE == "quick") {
  NUS <- c(2.5); NS <- c(100); R <- 3
  ITER <- 800; WARMUP <- 400; CHAINS <- 2
} else if (RUN_MODE == "pilot") {
  NUS <- c(0.2); NS <- c(600); R <- 15
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else if (RUN_MODE == "full_r10") {
  # Full (nu, n) grid -- all 5 regimes x all 4 sample sizes -- but with
  # R=10 reps/cell instead of 30, to cut cost roughly 3x while still
  # covering every combination. R=10 also matches the pilot rep count
  # already quoted in the draft's Table 2/3 captions.
  NUS <- c(0.2, 0.5, 1.0, 1.5, 2.5)
  NS  <- c(100, 200, 400, 600)
  R   <- 10
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
} else {
  NUS <- c(0.2, 0.5, 1.0, 1.5, 2.5)   # original 3 regimes + 2 new extremes:
                                        # nu=0.2 (strong overdispersion),
                                        # nu=2.5 (strong underdispersion)
  NS  <- c(100, 200, 400, 600)
  R   <- 30
  ITER <- 2000; WARMUP <- 1000; CHAINS <- 2
}
QUICK_TEST <- !(RUN_MODE %in% c("full", "full_r10"))

# ---------------------------------------------------------------------
# 1.1 Simulator (reused by Part 4 too)
# ---------------------------------------------------------------------
thin_operator <- function(x, alpha) sum(rbinom(x, size = 1, prob = alpha))

simul_zinarCMP <- function(n, alpha, lambda, nu, zi_prob) {
  zinp_inar <- integer(n)
  zinp_inar[1] <- rcmp(1, lambda, nu)
  for (t in 2:n) {
    thinned <- thin_operator(zinp_inar[t - 1], alpha)
    if (runif(1) < zi_prob) {
      innovation <- 0
    } else {
      innovation <- rcmp(1, lambda, nu)
    }
    zinp_inar[t] <- thinned + innovation
  }
  zinp_inar
}

# ---------------------------------------------------------------------
# 1.2 Exact CMP moments (reused by Part 4; kmax raised to 600: nu=0.2's
#     heavier tail needs more terms to converge than the original
#     nu in {0.5,1,1.5} did)
# ---------------------------------------------------------------------
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

cmp_moments <- function(lambda, nu, kmax = 600) {
  k <- 0:kmax
  logw <- k * log(lambda) - nu * lgamma(k + 1)
  logZ <- logsumexp(logw)
  p <- exp(logw - logZ)
  mu <- sum(k * p)
  c(mean = mu, var = sum(k^2 * p) - mu^2)
}

solve_lambda_for_mean <- function(nu, target_mu, kmax = 600) {
  f <- function(lam) cmp_moments(lam, nu, kmax)["mean"] - target_mu
  uniroot(f, c(1e-6, 1e6), tol = 1e-8)$root
}

# ---------------------------------------------------------------------
# 1.2b R port of the CMP transition log-likelihood (exact match to
#      ZIINAR1-CMP-fast-reparam.stan's transformed-parameters block) and
#      the EAIC/EBIC/DIC/WAIC1/WAIC2 model-selection criteria for a CMP
#      fit -- exact same formulas as ZIHINAR1::get_mod_sel(), so all four
#      models (CMP + Part 4's ZIP/ZINB/ZIGP) are judged identically.
#      Moved up here (out of Part 4) so THIS PART's own grid loop below
#      can compute and stash these numbers at fit time -- Part 4 then
#      reuses them instead of re-simulating + re-fitting CMP from
#      scratch at every (nu,n) cell it shares with this grid.
# ---------------------------------------------------------------------
cmp_log_Z_R <- function(lambda, nu, M, hybrid_tol, lgam) {
  log_lambda <- log(lambda)
  test <- exp(-log_lambda / nu)
  if (test < hybrid_tol) {
    nu * exp(log_lambda / nu) - ((nu - 1) / (2 * nu)) * log_lambda -
      ((nu - 1) / 2) * log(2 * pi) - 0.5 * log(nu)
  } else {
    r <- 0:M
    logsumexp(r * log_lambda - nu * lgam[r + 1])
  }
}

cmp_transition_loglik <- function(y, alpha, lambda, nu, rho, M = 300, hybrid_tol = 1e-6) {
  Tt <- length(y)
  lgam <- lgamma(0:(M + 1) + 1)
  log_Z <- cmp_log_Z_R(lambda, nu, M, hybrid_tol, lgam)
  log_lam <- log(lambda); log_rho <- log(rho); log1mrho <- log1p(-rho)
  out <- numeric(Tt - 1)
  for (t in 2:Tt) {
    yt <- y[t]; yt1 <- y[t - 1]; p <- min(yt1, yt)
    lbin0 <- dbinom(0, yt1, alpha, log = TRUE)
    lcmp0 <- yt * log_lam - nu * lgam[yt + 1] - log_Z
    lterm0 <- if (yt == 0) lbin0 + logsumexp(c(log_rho, log1mrho + lcmp0)) else
                            lbin0 + log1mrho + lcmp0
    if (p == 0) {
      out[t - 1] <- lterm0
    } else {
      lterms <- numeric(p + 1); lterms[1] <- lterm0
      for (j in 1:p) {
        lbinj <- dbinom(j, yt1, alpha, log = TRUE)
        diff  <- yt - j
        lcmpj <- diff * log_lam - nu * lgam[diff + 1] - log_Z
        lterms[j + 1] <- if (yt == j) lbinj + logsumexp(c(log_rho, log1mrho + lcmpj)) else
                                       lbinj + log1mrho + lcmpj
      }
      out[t - 1] <- logsumexp(lterms)
    }
  }
  out
}

get_mod_sel_cmp <- function(y, stan_fit, M = 300, hybrid_tol = 1e-6) {
  aic <- rstan::extract(stan_fit, pars = "aic")[[1]]; eaic <- mean(aic)
  bic <- rstan::extract(stan_fit, pars = "bic")[[1]]; ebic <- mean(bic)

  ph <- summary(stan_fit, pars = c("alpha", "rho", "lambda", "nu"))$summary
  logphat <- sum(cmp_transition_loglik(y, ph["alpha", "mean"], ph["lambda", "mean"],
                                        ph["nu", "mean"], ph["rho", "mean"],
                                        M, hybrid_tol))

  ll <- rstan::extract(stan_fit, pars = "ll")[[1]]
  pdic <- 2 * (logphat - mean(ll))
  dic  <- -2 * logphat + 2 * pdic

  loglik_mat <- rstan::extract(stan_fit, pars = "log_lik")[[1]][, 2:length(y), drop = FALSE]
  lik_mat    <- exp(loglik_mat)
  lppd   <- sum(log(colMeans(lik_mat)))
  pwaic1 <- 2 * sum(log(colMeans(lik_mat)) - colMeans(loglik_mat))
  pwaic2 <- sum(matrixStats::colVars(loglik_mat))
  waic1  <- -2 * (lppd - pwaic1)
  waic2  <- -2 * (lppd - pwaic2)

  data.frame(EAIC = eaic, EBIC = ebic, DIC = dic, WAIC1 = waic1, WAIC2 = waic2)
}

# ---------------------------------------------------------------------
# 1.3 Comparators independent of the Stan fit (Part 1 only)
# ---------------------------------------------------------------------
cls_alpha <- function(y) {
  y1 <- y[-length(y)]; y2 <- y[-1]
  m1 <- mean(y1); m2 <- mean(y2)
  sum((y1 - m1) * (y2 - m2)) / sum((y1 - m1)^2)
}

mom_oracle_lambda_nu <- function(y, alpha_true, rho_true, kmax = 600) {
  mu_hat <- mean(y); s2_hat <- var(y)
  muU_target <- mu_hat * (1 - alpha_true) / (1 - rho_true)
  rhs <- (1 - alpha_true^2) * s2_hat - alpha_true * (1 - alpha_true) * mu_hat
  sigU2_target <- (rhs - rho_true * (1 - rho_true) * muU_target^2) / (1 - rho_true)
  if (!is.finite(sigU2_target) || sigU2_target <= 0 || muU_target <= 0) {
    return(c(lambda = NA_real_, nu = NA_real_))
  }
  resid <- function(nu) {
    lam <- tryCatch(solve_lambda_for_mean(nu, muU_target, kmax), error = function(e) NA)
    if (is.na(lam)) return(NA_real_)
    cmp_moments(lam, nu, kmax)["var"] - sigU2_target
  }
  nu_hat <- tryCatch({
    r_lo <- resid(0.02); r_hi <- resid(12)
    if (!is.finite(r_lo) || !is.finite(r_hi) || sign(r_lo) == sign(r_hi)) return(NA_real_)
    uniroot(resid, c(0.02, 12), tol = 1e-6)$root
  }, error = function(e) NA_real_)
  if (is.na(nu_hat)) return(c(lambda = NA_real_, nu = NA_real_))
  lam_hat <- tryCatch(solve_lambda_for_mean(nu_hat, muU_target, kmax), error = function(e) NA_real_)
  c(lambda = lam_hat, nu = nu_hat)
}

# ---------------------------------------------------------------------
# 1.4 Calibrate lambda per regime
# ---------------------------------------------------------------------
lambdas <- setNames(sapply(NUS, solve_lambda_for_mean, target_mu = TARGET_MU), as.character(NUS))
cat("Calibrated lambdas (target CMP mean =", TARGET_MU, "):\n"); print(lambdas)
for (nu in NUS) {
  mv <- cmp_moments(lambdas[[as.character(nu)]], nu)
  cat(sprintf("  nu=%.2f  lambda=%.4f  mean=%.3f  var=%.3f  d_U=%.3f\n",
              nu, lambdas[[as.character(nu)]], mv["mean"], mv["var"], mv["var"] / mv["mean"]))
}

# ---------------------------------------------------------------------
# 1.5 Compile once (mod is reused by Part 4's CMP fits too)
# ---------------------------------------------------------------------
cat("\nCompiling", STAN_FILE, "...\n")
mod <- stan_model(STAN_FILE)

# ---------------------------------------------------------------------
# 1.6 Main grid -- fills `results` (raw per-replicate draws) used by
#     both the LaTeX table emitter and the plots in PART 2 below.
#     NOTE: the reparam Stan file's `parameters` block has `mu_cmp`, not
#     `lambda` -- but `lambda` is still available to pull posterior
#     summaries from because it's declared in `transformed parameters`
#     (lambda := mu_cmp^nu), so `pars = c("alpha","lambda","nu","rho")`
#     below works unchanged against ZIINAR1-CMP-fast-reparam.stan.
# ---------------------------------------------------------------------
results <- list()
# cmp_cache stores, per (nu,n) cell and rep, the exact y that was
# simulated and this fit's EAIC/EBIC/DIC/WAIC1/WAIC2 -- Part 4 reuses
# these for cells it shares with this grid, so the CMP model is never
# re-simulated + re-fit from scratch there (only ZIP/ZINB/ZIGP are).
cmp_cache <- list()
t_start <- Sys.time()

# One rep's full unit of work -- called in parallel across reps via
# mclapply below. Returns NULL on any failure (sampling error or the
# fit simply didn't happen); returns list(row=..., cache=...) otherwise.
# Must be self-contained: everything it uses (mod, ALPHA_TRUE, M,
# HYBRID_TOL, FF, CHAINS, ITER, WARMUP, cls_alpha, mom_oracle_lambda_nu,
# get_mod_sel_cmp, RHO_TRUE) is inherited by the forked child processes
# automatically (fork copies the parent's whole memory), so nothing
# needs to be passed in explicitly beyond what varies per call.
fit_one_rep_part1 <- function(r, n, nu_true, lam_true) {
  y <- simul_zinarCMP(n, ALPHA_TRUE, lam_true, nu_true, RHO_TRUE)

  fit <- tryCatch(
    sampling(mod,
             data = list(T = n, y = y, M = M, hybrid_tol = HYBRID_TOL, ff = FF),
             chains = CHAINS, iter = ITER, warmup = WARMUP,
             seed = 1000 + r, refresh = 0,
             control = list(adapt_delta = 0.95)),
    error = function(e) { message("  [rep ", r, "] sampling failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(NULL)

  s <- summary(fit, pars = c("alpha", "lambda", "nu", "rho"))$summary
  rhat_ok <- all(s[, "Rhat"] < 1.05, na.rm = TRUE)

  post <- rstan::extract(fit, pars = c("lambda", "nu"))
  lamnu_cor <- suppressWarnings(cor(post$lambda, post$nu))

  a_cls <- cls_alpha(y)
  mom <- mom_oracle_lambda_nu(y, ALPHA_TRUE, RHO_TRUE)

  # Cache this fit's model-selection criteria for Part 4 to reuse -- if
  # this fails for some numerical reason, just skip caching (Part 4
  # falls back to fitting CMP fresh for this rep).
  crit_cmp <- tryCatch(get_mod_sel_cmp(y, fit, M, HYBRID_TOL), error = function(e) NULL)

  row <- data.frame(
    rep = r,
    alpha_mean = s["alpha", "mean"],   alpha_sd = s["alpha", "sd"],
    lambda_mean = s["lambda", "mean"], lambda_sd = s["lambda", "sd"],
    nu_mean = s["nu", "mean"],         nu_sd = s["nu", "sd"],
    rho_mean = s["rho", "mean"],       rho_sd = s["rho", "sd"],
    rhat_ok = rhat_ok,
    lambda_nu_cor = lamnu_cor,
    alpha_cls = a_cls,
    lambda_mom = unname(mom["lambda"]), nu_mom = unname(mom["nu"])
  )
  cat(sprintf("  rep %d/%d  Rhat_ok=%s  alpha=%.3f  lambda=%.3f  nu=%.3f  rho=%.3f  cor(lam,nu)=%.3f\n",
              r, R, rhat_ok, s["alpha", "mean"], s["lambda", "mean"],
              s["nu", "mean"], s["rho", "mean"], lamnu_cor))

  list(row = row, cache = list(y = y, crit_cmp = crit_cmp))
}

ok_result <- function(x) !is.null(x) && !inherits(x, "try-error")

for (nu_true in NUS) {
  lam_true <- lambdas[[as.character(nu_true)]]
  for (n in NS) {
    key <- paste0("nu=", nu_true, "_n=", n)
    cat("\n===", key, "=== (", N_WORKERS, "worker(s) in parallel; output below may interleave)\n")

    rep_results <- mclapply(seq_len(R), fit_one_rep_part1, n = n, nu_true = nu_true,
                             lam_true = lam_true, mc.cores = N_WORKERS, mc.preschedule = FALSE)

    rows <- lapply(rep_results, function(x) if (ok_result(x)) x$row else NULL)
    cmp_cache[[key]] <- lapply(rep_results, function(x) if (ok_result(x)) x$cache else NULL)

    results[[key]] <- do.call(rbind, rows[!sapply(rows, is.null)])
    saveRDS(results, "sim_results_partial.rds")
    saveRDS(cmp_cache, "cmp_cache_partial.rds")   # checkpointed every cell now, not just at the end
    cat(sprintf("  elapsed so far: %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}
saveRDS(results, "sim_results_final.rds")
saveRDS(cmp_cache, "cmp_cache.rds")
cat("\nTOTAL TIME (min):", as.numeric(difftime(Sys.time(), t_start, units = "mins")), "\n")

# ---------------------------------------------------------------------
# 1.7 Summarize + write CSV/LaTeX (generic over however many regimes
#     are in NUS)
# ---------------------------------------------------------------------
summarize_cell <- function(df) {
  data.frame(
    n_reps = nrow(df),
    alpha_MCmean = mean(df$alpha_mean),  alpha_MCsd = sd(df$alpha_mean),
    lambda_MCmean = mean(df$lambda_mean), lambda_MCsd = sd(df$lambda_mean),
    nu_MCmean = mean(df$nu_mean),        nu_MCsd = sd(df$nu_mean),
    rho_MCmean = mean(df$rho_mean),      rho_MCsd = sd(df$rho_mean),
    alpha_cls_mean = mean(df$alpha_cls, na.rm = TRUE),
    nu_mom_mean = mean(df$nu_mom, na.rm = TRUE),
    rhat_ok_frac = mean(df$rhat_ok),
    lambda_nu_cor_mean = mean(df$lambda_nu_cor, na.rm = TRUE)
  )
}

summary_tab <- do.call(rbind, lapply(names(results), function(k) {
  cbind(cell = k, summarize_cell(results[[k]]))
}))
write.csv(summary_tab, "sim_summary_table.csv", row.names = FALSE)
cat("\n--- Summary table (also written to sim_summary_table.csv) ---\n")
print(summary_tab)

fmt <- function(m, s) sprintf("%.3f (%.3f)", m, s)
nu_regime_label <- function(nu) {
  if (nu < 0.5) return(sprintf("Strong overdispersion regime ($\\nu=%.1f$)", nu))
  if (nu < 1.0) return(sprintf("Overdispersion regime ($\\nu=%.1f$)", nu))
  if (nu == 1.0) return("Equidispersion regime ($\\nu=1.0$)")
  if (nu <= 1.5) return(sprintf("Underdispersion regime ($\\nu=%.1f$)", nu))
  return(sprintf("Strong underdispersion regime ($\\nu=%.1f$)", nu))
}
emit_latex <- function(summary_tab, lambdas, ns, nus) {
  lines <- c()
  for (nu in sort(nus)) {
    lines <- c(lines, paste0("\\multicolumn{6}{l}{\\textbf{", nu_regime_label(nu), "}} \\\\"))
    for (par in c("alpha", "lambda", "nu", "rho")) {
      true_val <- switch(par, alpha = 0.30, rho = 0.30, nu = nu, lambda = lambdas[[as.character(nu)]])
      cells <- sapply(ns, function(n) {
        key <- paste0("nu=", nu, "_n=", n)
        row <- summary_tab[summary_tab$cell == key, ]
        if (nrow(row) == 0) return("--")
        fmt(row[[paste0(par, "_MCmean")]], row[[paste0(par, "_MCsd")]])
      })
      symb <- switch(par, alpha = "\\alpha", lambda = "\\lambda", nu = "\\nu", rho = "\\rho")
      lines <- c(lines, sprintf("$%s$ & %.2f & %s & %s & %s & %s \\\\",
                                 symb, true_val, cells[1], cells[2], cells[3], cells[4]))
    }
    lines <- c(lines, "\\midrule")
  }
  writeLines(lines, "sim_table_body.tex")
  cat("\nWrote sim_table_body.tex\n")
}
if (!QUICK_TEST) emit_latex(summary_tab, lambdas, NS, NUS)


############################################################################
# PART 2 -- PLOTS (parameter recovery, bias, model-comparison win rate)
#           Uses `summary_tab` and `NUS`/`NS` from Part 1 directly --
#           no CSV re-read needed since it's the same script/session.
#           Nothing is plotted until you actually run Part 1 above (or
#           point `summary_tab` at a saved sim_summary_table.csv via
#           read.csv() instead).
############################################################################

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"
YELLOW <- "#eda100"; PURPLE <- "#a259d9"
regime_colors <- c(BLUE, ORANGE, AQUA, YELLOW, PURPLE)

build_regime_list <- function(summary_tab, n_vals) {
  tab <- summary_tab
  tab$nu <- as.numeric(sub("nu=([^_]+)_n=.*", "\\1", tab$cell))
  tab$n  <- as.numeric(sub(".*_n=", "", tab$cell))
  nus <- sort(unique(tab$nu))
  out <- list()
  for (i in seq_along(nus)) {
    nu <- nus[i]
    sub_tab <- tab[tab$nu == nu, ]
    sub_tab <- sub_tab[match(n_vals, sub_tab$n), ]
    out[[as.character(nu)]] <- list(
      color = regime_colors[((i - 1) %% length(regime_colors)) + 1],
      true_nu = nu, true_lambda = lambdas[[as.character(nu)]],
      alpha  = list(mean = sub_tab$alpha_MCmean,  sd = sub_tab$alpha_MCsd),
      lambda = list(mean = sub_tab$lambda_MCmean, sd = sub_tab$lambda_MCsd),
      nu     = list(mean = sub_tab$nu_MCmean,     sd = sub_tab$nu_MCsd),
      rho    = list(mean = sub_tab$rho_MCmean,    sd = sub_tab$rho_MCsd)
    )
  }
  out
}

param_labels <- c(alpha = "hat(alpha)", lambda = "hat(lambda)", nu = "hat(nu)", rho = "hat(rho)")
panel_order  <- c("nu", "alpha", "lambda", "rho")

plot_recovery <- function(data_list, n_vals, file_stem, title_suffix = "") {
  regimes <- names(data_list); n_regimes <- length(regimes)
  jitters <- (seq_len(n_regimes) - (n_regimes + 1) / 2) * 0.03 * diff(range(n_vals))
  make_plot <- function() {
    op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), oma = c(4.5, 0, 2, 0)); on.exit(par(op))
    for (pname in panel_order) {
      all_means <- unlist(lapply(data_list, function(d) d[[pname]]$mean))
      all_sds   <- unlist(lapply(data_list, function(d) d[[pname]]$sd))
      ylim <- range(c(all_means - all_sds, all_means + all_sds), na.rm = TRUE)
      plot(NA, xlim = range(n_vals), ylim = ylim, xaxt = "n",
           xlab = "Sample size n", ylab = paste0(param_labels[pname], " (posterior mean +/- SD)"),
           main = param_labels[pname])
      axis(1, at = n_vals)
      for (i in seq_along(regimes)) {
        reg <- regimes[i]; d <- data_list[[reg]]; x <- n_vals + jitters[i]
        abline(h = if (pname == "nu") d$true_nu else if (!is.na(d$true_lambda) && pname == "lambda") d$true_lambda else NA,
               col = d$color, lty = 2, lwd = 1)
        arrows(x, d[[pname]]$mean - d[[pname]]$sd, x, d[[pname]]$mean + d[[pname]]$sd,
               angle = 90, code = 3, length = 0.03, col = d$color, lwd = 1.2)
        lines(x, d[[pname]]$mean, col = d$color, lwd = 1.4)
        points(x, d[[pname]]$mean, col = d$color, pch = 16, cex = 1)
      }
    }
    legend_labels <- sapply(regimes, function(r) sprintf("nu = %s", r))
    legend_cols <- sapply(data_list, function(d) d$color)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = legend_labels, col = legend_cols, lwd = 1.4, pch = 16,
           horiz = TRUE, bty = "n", inset = c(0, -0.01), xpd = TRUE)
    mtext(paste0("Parameter recovery across dispersion regimes and sample sizes", title_suffix),
          outer = TRUE, cex = 1.1, line = 0.5)
  }
  pdf(paste0(file_stem, ".pdf"), width = 8.5, height = 7); make_plot(); dev.off()
  png(paste0(file_stem, ".png"), width = 2000, height = 1650, res = 220); make_plot(); dev.off()
}

plot_bias <- function(data_list, n_vals, file_stem) {
  regimes <- names(data_list); n_regimes <- length(regimes)
  jitters <- (seq_len(n_regimes) - (n_regimes + 1) / 2) * 0.03 * diff(range(n_vals))
  make_plot <- function() {
    op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), oma = c(4.5, 0, 2, 0)); on.exit(par(op))
    for (pname in panel_order) {
      biases <- lapply(data_list, function(d) {
        true_val <- if (pname == "nu") d$true_nu else if (pname %in% c("alpha", "rho")) 0.30 else d$true_lambda
        d[[pname]]$mean - true_val
      })
      ylim <- range(unlist(biases), na.rm = TRUE)
      plot(NA, xlim = range(n_vals), ylim = ylim, xaxt = "n",
           xlab = "Sample size n", ylab = paste0("Bias: ", param_labels[pname], " - true value"),
           main = param_labels[pname])
      axis(1, at = n_vals); abline(h = 0, col = "grey40", lwd = 1)
      for (i in seq_along(regimes)) {
        reg <- regimes[i]; d <- data_list[[reg]]; x <- n_vals + jitters[i]
        lines(x, biases[[reg]], col = d$color, lwd = 1.4)
        points(x, biases[[reg]], col = d$color, pch = 16, cex = 1)
      }
    }
    legend_labels <- sapply(regimes, function(r) sprintf("nu = %s", r))
    legend_cols <- sapply(data_list, function(d) d$color)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = legend_labels, col = legend_cols, lwd = 1.4, pch = 16,
           horiz = TRUE, bty = "n", inset = c(0, -0.01), xpd = TRUE)
    mtext("Estimation bias (posterior mean - true value) by dispersion regime and n",
          outer = TRUE, cex = 1.1, line = 0.5)
  }
  pdf(paste0(file_stem, ".pdf"), width = 8.5, height = 7); make_plot(); dev.off()
  png(paste0(file_stem, ".png"), width = 2000, height = 1650, res = 220); make_plot(); dev.off()
}

# Win-rate plot: reads directly from Part 4's `win_rate` data.frame once
# that has run; falls back to the original 3-regime numbers already in
# the draft (Table~\ref{tab:model_comparison_winrate}) if Part 4 hasn't
# been run yet, so this still produces a figure either way.
plot_winrate <- function(win_rate_df = NULL) {
  models_wr  <- c("ZICMP", "ZIP", "ZINB", "ZIGP")
  model_colors <- c(BLUE, ORANGE, AQUA, YELLOW)
  if (is.null(win_rate_df)) {
    regimes_wr <- c("Overdispersion\n(nu=0.5)", "Equidispersion\n(nu=1.0)", "Underdispersion\n(nu=1.5)")
    winrate <- rbind(c(0.70, 0.00, 0.30, 0.00), c(0.00, 0.90, 0.10, 0.00), c(0.90, 0.10, 0.00, 0.00))
  } else {
    # win_rate_df may be old-style (just a "regime" column, e.g. "nu=0.5",
    # one row per regime) or Part 4's new-style (separate "nu"/"n" columns,
    # possibly several n's per nu) -- prefer the "nu" column when present
    # so ordering/labels don't depend on parsing "nu=0.5_n=600" strings.
    nu_num <- if (!is.null(win_rate_df$nu)) as.numeric(win_rate_df$nu) else as.numeric(sub("nu=", "", win_rate_df$regime))
    ord <- order(nu_num)
    win_rate_df <- win_rate_df[ord, ]
    regimes_wr <- if (!is.null(win_rate_df$nu)) sprintf("nu=%s", win_rate_df$nu) else win_rate_df$regime
    winrate <- as.matrix(win_rate_df[, models_wr])
  }
  colnames(winrate) <- models_wr; rownames(winrate) <- regimes_wr
  bp <- barplot(t(winrate), beside = TRUE, col = model_colors, ylim = c(0, 1.05),
                ylab = "Fraction of replicates with lowest WAIC2",
                names.arg = regimes_wr, legend.text = models_wr,
                args.legend = list(x = "top", horiz = TRUE, bty = "n", inset = c(0, -0.08)),
                main = "Which model wins, by dispersion regime")
  for (i in seq_len(nrow(winrate))) for (j in seq_len(ncol(winrate))) {
    v <- winrate[i, j]; if (!is.na(v) && v > 0) text(bp[j, i], v + 0.03, sprintf("%.2f", v), cex = 0.75)
  }
}

# Run Part 2 (only meaningful once Part 1 has populated `summary_tab`):
if (exists("summary_tab") && nrow(summary_tab) > 0) {
  n_vals <- NS
  data_all <- build_regime_list(summary_tab, n_vals)
  plot_recovery(data_all, n_vals, "fig_recovery")
  plot_bias(data_all, n_vals, "fig_bias")
  pdf("fig_modelcomparison.pdf", width = 7.5, height = 4.2)
  plot_winrate(if (exists("win_rate")) win_rate else NULL)
  dev.off()
  png("fig_modelcomparison.png", width = 1900, height = 1050, res = 220)
  plot_winrate(if (exists("win_rate")) win_rate else NULL)
  dev.off()
  cat("\nWrote fig_recovery.{pdf,png}, fig_bias.{pdf,png}, fig_modelcomparison.{pdf,png}\n")
}


############################################################################
# PART 3 -- PREDICTIVE ANALYSIS
#           qoi + print(fit, pars=qoi) + y_pred-based forecast plot,
#           adapted to this model's real parameters/generated quantities.
############################################################################

qoi <- c("alpha", "lambda", "nu", "rho", "aic", "bic")
# print(fit, pars = qoi)   # once you have a fitted `fit` in scope

plot_predictive <- function(fit, y, ff, probs = c(0.1, 0.9),
                             series_name = "", model_name = "ZICMP", ylim = NULL) {
  n <- length(y)
  stopifnot(ff > 0, ff < n)

  # y_pred has ff+1 entries: y_pred[1] = y[T] (last TRAINING point, an
  # anchor, not a forecast); y_pred[2:(ff+1)] are the ff chained
  # one-step-ahead forecasts, matching y[(n-ff+1):n].
  fitPred   <- summary(fit, pars = "y_pred", probs = probs)$summary
  fitPredM  <- fitPred[-1, "mean"]
  loName    <- grep(paste0(probs[1] * 100, "%"), colnames(fitPred), value = TRUE)[1]
  hiName    <- grep(paste0(probs[2] * 100, "%"), colnames(fitPred), value = TRUE)[1]
  fitPredLo <- fitPred[-1, loName]
  fitPredHi <- fitPred[-1, hiName]

  y_obs <- y[(n - ff + 1):n]
  if (is.null(ylim)) {
    ylim <- range(c(y_obs, fitPredM, fitPredLo, fitPredHi), na.rm = TRUE)
    ylim[2] <- ylim[2] * 1.05
  }

  plot(y_obs, type = "o", pch = 16, ylim = ylim,
       xlab = "Held-out time index", ylab = "Count",
       main = sprintf("%s%s one-step-ahead forecasts (%d%% band)",
                       model_name, if (nzchar(series_name)) paste0(" -- ", series_name) else "",
                       round(100 * (probs[2] - probs[1]))))
  lines(fitPredM, col = "#2a78d6", lwd = 1.8)
  lines(fitPredLo, col = "#2a78d6", lty = 2)
  lines(fitPredHi, col = "#2a78d6", lty = 2)
  legend("topleft", legend = c("Observed", "Posterior predictive mean",
                                sprintf("%d%%-%d%% interval", 100 * probs[1], 100 * probs[2])),
         col = c("black", "#2a78d6", "#2a78d6"), lty = c(1, 1, 2), pch = c(16, NA, NA),
         bty = "n", cex = 0.85)

  invisible(data.frame(t = seq_len(ff), observed = y_obs,
                        pred_mean = fitPredM, pred_lo = fitPredLo, pred_hi = fitPredHi))
}

# Worked example: fit one real (or simulated) series with a genuine
# holdout (ff > 0) and run the predictive plot. Gated with `if (FALSE)`
# so sourcing this script doesn't try to fit anything on its own --
# fill in your own series and flip to TRUE, or just copy the body out.
if (FALSE) {
  y_full <- scan("sexoffences.txt")      # replace with your real loader
  n      <- length(y_full)
  ff     <- round(0.20 * n)              # same 20% holdout as forecast_evaluation.R
  n_tr   <- n - ff
  y_train <- y_full[1:n_tr]

  fit_pred <- sampling(mod,              # reuses the `mod` compiled in Part 1
                        data = list(T = n_tr, y = y_train, M = M,
                                    hybrid_tol = HYBRID_TOL, ff = ff),
                        chains = 4, iter = 2000, warmup = 1000, seed = 1,
                        control = list(adapt_delta = 0.95))

  print(fit_pred, pars = qoi)

  pred_tab    <- plot_predictive(fit_pred, y_full, ff, probs = c(0.1, 0.9),
                                  series_name = "sex offenses")
  pred_tab_95 <- plot_predictive(fit_pred, y_full, ff, probs = c(0.025, 0.975),
                                  series_name = "sex offenses")
}


############################################################################
# PART 4 -- MODEL COMPARISON vs. ZIP / ZINB / ZIGP (via ZIHINAR1)
#           Same 5 dispersion regimes AND the same n-grid as Part 1 --
#           EAIC/EBIC/DIC/WAIC1/WAIC2 for all four models + win rates,
#           crossed over every (nu, n) cell. This is the "does the model
#           actually get SELECTED, not just recovered" question, and is
#           what feeds fig_modelcomparison above and Table 2/3 in the
#           draft. Reuses thin_operator/simul_zinarCMP/cmp_moments/
#           solve_lambda_for_mean/mod from Part 1 -- run Part 1 first
#           (at least through section 1.5, the Stan compile) so those
#           are in scope, or move this block above Part 1 if you'd
#           rather run comparison-only.
#
#           NOTE ON COST: this now fits 4 models (CMP + ZIP + ZINB + ZIGP)
#           at EVERY (nu, n) cell, vs. Part 1's single CMP fit per cell --
#           so "full" mode here is roughly 4x (models) the cost of Part 1
#           at the same R, summed over the same 5x4=20 cells. Use
#           CMP_NS/CMP_NUS below to shrink the grid (e.g. drop small n,
#           or only the 2 new extreme regimes) if that's too slow.
# ---------------------------------------------------------------------
# 4.0 Config (own regime/n/rep settings; defaults to the SAME grid as
#     Part 1 -- set CMP_NS to a subset, e.g. c(600), to go back to the
#     single-n comparison this used to do)
# ---------------------------------------------------------------------
COMPARE_EXTRA_ONLY <- FALSE  # TRUE = only the 2 new regimes (nu=0.2, 2.5),
                              # e.g. if you already have results for the
                              # original 3 from an earlier compare_innovations.R run

if (RUN_MODE == "quick") {
  CMP_NUS <- c(2.5); CMP_NS <- c(100); CMP_R <- 2
  CMP_ITER <- 800; CMP_WARMUP <- 400; CMP_CHAINS <- 2
} else if (RUN_MODE == "pilot") {
  CMP_NUS <- if (COMPARE_EXTRA_ONLY) c(0.2, 2.5) else c(0.2, 0.5, 1.0, 1.5, 2.5)
  CMP_NS <- c(600)
  CMP_R <- 10
  CMP_ITER <- 2000; CMP_WARMUP <- 1000; CMP_CHAINS <- 2
} else if (RUN_MODE == "full_r10") {
  # Full (nu, n) grid, R=10 reps/cell -- same tradeoff as Part 1 above.
  CMP_NUS <- if (COMPARE_EXTRA_ONLY) c(0.2, 2.5) else c(0.2, 0.5, 1.0, 1.5, 2.5)
  CMP_NS <- c(100, 200, 400, 600)   # same n-grid as Part 1
  CMP_R <- 10
  CMP_ITER <- 2000; CMP_WARMUP <- 1000; CMP_CHAINS <- 2
} else {
  CMP_NUS <- if (COMPARE_EXTRA_ONLY) c(0.2, 2.5) else c(0.2, 0.5, 1.0, 1.5, 2.5)
  CMP_NS <- c(100, 200, 400, 600)   # same n-grid as Part 1
  CMP_R <- 30
  CMP_ITER <- 2000; CMP_WARMUP <- 1000; CMP_CHAINS <- 2
}

# ---------------------------------------------------------------------
# 4.1 cmp_log_Z_R / cmp_transition_loglik / get_mod_sel_cmp -- MOVED to
#     Part 1, section 1.2b (so Part 1's own grid loop can compute and
#     cache these once, at fit time, and Part 4 below can reuse them
#     instead of re-simulating + re-fitting CMP from scratch). Nothing
#     else here changes.
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# 4.2 Calibrate lambda for the comparison regimes (own lookup, since
#     CMP_NUS need not match Part 1's NUS if you set COMPARE_EXTRA_ONLY)
# ---------------------------------------------------------------------
cmp_lambdas <- setNames(sapply(CMP_NUS, solve_lambda_for_mean, target_mu = TARGET_MU),
                         as.character(CMP_NUS))
cat("\n[Part 4] Calibrated lambdas (target CMP mean =", TARGET_MU, "):\n"); print(cmp_lambdas)

# ---------------------------------------------------------------------
# 4.3 Main comparison loop -- reuses `mod` (compiled in Part 1, section
#     1.5) for the CMP fits; ZIP/ZINB/ZIGP compile lazily inside
#     ZIHINAR1::get_stanfit() the first time each is called.
#
#     AVOIDING REPEATED WORK: when this (nu, n, r) cell was already
#     fit by Part 1's own grid (same y, same CMP model, cached in
#     `cmp_cache`), reuse that y and its already-computed EAIC/EBIC/
#     DIC/WAIC criteria instead of re-simulating and re-fitting CMP
#     from scratch -- only ZIP/ZINB/ZIGP are fit fresh in that case.
#     Falls back to simulating + fitting CMP here (the original
#     behavior) for any cell Part 1 didn't cover -- e.g. CMP_NUS/
#     CMP_NS wider than NUS/NS, or r > R, which happens whenever
#     RUN_MODE gives Part 1 and Part 4 different grids (quick/pilot).
# ---------------------------------------------------------------------
compare_results <- list()
t0 <- Sys.time()
if (!exists("cmp_cache")) cmp_cache <- list()   # empty if Part 1 hasn't run this session

# One rep's full unit of work for Part 4 -- called in parallel across
# reps via mclapply below, same pattern as Part 1's fit_one_rep_part1.
# Returns NULL on failure; otherwise list(row=..., reused=TRUE/FALSE).
fit_one_rep_part4 <- function(r, nu_true, n_val, lam_true, key) {
  cached_cell <- cmp_cache[[key]]   # NULL if Part 1 never ran this cell
  cached_rep <- if (!is.null(cached_cell) && r <= length(cached_cell)) cached_cell[[r]] else NULL
  reuse_cmp <- !is.null(cached_rep) && !is.null(cached_rep$crit_cmp)

  if (reuse_cmp) {
    y <- cached_rep$y
    crit_cmp <- cached_rep$crit_cmp
  } else {
    y <- simul_zinarCMP(n_val, ALPHA_TRUE, lam_true, nu_true, RHO_TRUE)
    fit_cmp <- tryCatch(
      sampling(mod, data = list(T = n_val, y = y, M = M,
                                 hybrid_tol = HYBRID_TOL, ff = 0),
               chains = CMP_CHAINS, iter = CMP_ITER, warmup = CMP_WARMUP,
               seed = 2000 + r, refresh = 0, control = list(adapt_delta = 0.95)),
      error = function(e) { message("  [rep ", r, "] ZICMP fit failed: ", conditionMessage(e)); NULL })
    if (is.null(fit_cmp)) return(NULL)
    crit_cmp <- get_mod_sel_cmp(y, fit_cmp, M, HYBRID_TOL)
  }

  fit_poi <- tryCatch(
    ZIHINAR1::get_stanfit(mod_type = "zi", distri = "poi", y = y,
                           n_pred = 0, chains = CMP_CHAINS, iter = CMP_ITER,
                           warmup = CMP_WARMUP, seed = 2000 + r),
    error = function(e) { message("  [rep ", r, "] ZIP fit failed: ", conditionMessage(e)); NULL })

  fit_nb <- tryCatch(
    ZIHINAR1::get_stanfit(mod_type = "zi", distri = "nb", y = y,
                           n_pred = 0, chains = CMP_CHAINS, iter = CMP_ITER,
                           warmup = CMP_WARMUP, seed = 2000 + r),
    error = function(e) { message("  [rep ", r, "] ZINB fit failed: ", conditionMessage(e)); NULL })

  fit_gp <- tryCatch(
    ZIHINAR1::get_stanfit(mod_type = "zi", distri = "gp", y = y,
                           n_pred = 0, chains = CMP_CHAINS, iter = CMP_ITER,
                           warmup = CMP_WARMUP, seed = 2000 + r),
    error = function(e) { message("  [rep ", r, "] ZIGP fit failed: ", conditionMessage(e)); NULL })

  if (is.null(fit_poi) || is.null(fit_nb) || is.null(fit_gp)) return(NULL)

  crit_poi <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "poi", stan_fit = fit_poi)
  crit_nb  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "nb",  stan_fit = fit_nb)
  crit_gp  <- ZIHINAR1::get_mod_sel(y = y, mod_type = "zi", distri = "gp",  stan_fit = fit_gp)

  row <- rbind(
    cbind(model = "ZICMP", rep = r, crit_cmp),
    cbind(model = "ZIP",   rep = r, crit_poi),
    cbind(model = "ZINB",  rep = r, crit_nb),
    cbind(model = "ZIGP",  rep = r, crit_gp)
  )
  cat(sprintf("  rep %d/%d  [%s]  EAIC[cmp,poi,nb,gp] = %.1f, %.1f, %.1f, %.1f  |  WAIC2[cmp,poi,nb,gp] = %.1f, %.1f, %.1f, %.1f\n",
              r, CMP_R, if (reuse_cmp) "CMP reused from Part 1" else "CMP fit fresh",
              crit_cmp$EAIC, crit_poi$EAIC, crit_nb$EAIC, crit_gp$EAIC,
              crit_cmp$WAIC2, crit_poi$WAIC2, crit_nb$WAIC2, crit_gp$WAIC2))

  list(row = row, reused = reuse_cmp)
}

total_reused <- 0L; total_fresh <- 0L

for (nu_true in CMP_NUS) {
  lam_true <- cmp_lambdas[[as.character(nu_true)]]

  for (n_val in CMP_NS) {
    key <- paste0("nu=", nu_true, "_n=", n_val)
    cat("\n=== [Part 4]", key, "=== (", N_WORKERS, "worker(s) in parallel; output below may interleave)\n")

    rep_results <- mclapply(seq_len(CMP_R), fit_one_rep_part4, nu_true = nu_true, n_val = n_val,
                             lam_true = lam_true, key = key, mc.cores = N_WORKERS, mc.preschedule = FALSE)

    rows <- lapply(rep_results, function(x) if (ok_result(x)) x$row else NULL)
    cell_reused <- sum(sapply(rep_results, function(x) ok_result(x) && isTRUE(x$reused)))
    cell_fresh  <- sum(sapply(rep_results, function(x) ok_result(x) && !isTRUE(x$reused)))
    total_reused <- total_reused + cell_reused
    total_fresh  <- total_fresh + cell_fresh

    compare_results[[key]] <- do.call(rbind, rows[!sapply(rows, is.null)])
    saveRDS(compare_results, "compare_results_partial.rds")
    cat(sprintf("  elapsed so far: %.1f min  (this cell: %d reused, %d fresh CMP)\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), cell_reused, cell_fresh))
  }
}
saveRDS(compare_results, "compare_results_final.rds")
cat(sprintf("\n[Part 4] CMP reused from Part 1's cache: %d reps total  |  CMP fit fresh here: %d reps total\n",
            total_reused, total_fresh))

# ---------------------------------------------------------------------
# 4.4 Summarize: mean criteria per model per (regime, n) cell, and win
#     rate -- `win_rate` here feeds PART 2's plot_winrate() above if you
#     rerun that section afterward (or just call
#     plot_winrate(win_rate) / pdf(...); plot_winrate(win_rate); dev.off()
#     directly once this part has run).
# ---------------------------------------------------------------------
compare_summary_tab <- do.call(rbind, lapply(names(compare_results), function(k) {
  df <- compare_results[[k]]
  agg <- aggregate(cbind(EAIC, EBIC, DIC, WAIC1, WAIC2) ~ model, df, mean)
  parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
  cbind(regime = k, nu = parts[2], n = parts[3], agg)
}))
write.csv(compare_summary_tab, "compare_summary_table.csv", row.names = FALSE)

win_rate <- do.call(rbind, lapply(names(compare_results), function(k) {
  df <- compare_results[[k]]
  reps <- unique(df$rep)
  winners <- sapply(reps, function(rr) {
    sub <- df[df$rep == rr, ]
    sub$model[which.min(sub$WAIC2)]
  })
  parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
  # NOTE: table()/prop.table() results keep R's "table" class even after
  # t(); data.frame() silently unpacks a "table"-classed argument into
  # long format (Var1/Var2/Freq) instead of wide columns -- unclass()
  # first so ZICMP/ZIP/ZINB/ZIGP land as actual named columns.
  wr_tab <- prop.table(table(factor(winners, levels = c("ZICMP", "ZIP", "ZINB", "ZIGP"))))
  data.frame(regime = k, nu = parts[2], n = parts[3],
             as.list(unclass(wr_tab)), check.names = FALSE)
}))
write.csv(win_rate, "compare_win_rate.csv", row.names = FALSE)

cat("\n--- [Part 4] Mean model-selection criteria by regime (lower = better) ---\n")
print(compare_summary_tab)
cat("\n--- [Part 4] Fraction of replicates where each model has the LOWEST WAIC2 ---\n")
print(win_rate)

# ---------------------------------------------------------------------
# 4.5 Emit LaTeX. The draft's current Table~\ref{tab:model_comparison} /
#     Table~\ref{tab:model_comparison_winrate} are captioned "n=600" only
#     (3 regimes, no n dimension) -- with n now crossed in, two sets of
#     files are written:
#       compare_table_n600_body.tex / compare_winrate_n600_body.tex
#         -- the n=600 slice only, SAME shape as the tables already in
#            the draft (drop-in paste, just update the caption's R= and
#            regime count -- now 5 regimes, not 3 -- and rerun bibtex/
#            pdflatex; no table structure edits needed).
#       compare_table_full_body.tex / compare_winrate_full_body.tex
#         -- all 5x4=20 (regime,n) cells, one \multirow block per cell.
#            This is long (20 blocks); decide whether it replaces the
#            n=600 table in the main text or becomes a supplementary
#            table in the appendix -- both are just written out below,
#            nothing is auto-inserted into the .tex draft itself.
# ---------------------------------------------------------------------
model_order <- c("ZICMP", "ZIGP", "ZINB", "ZIP")
crit_cols   <- c("EAIC", "EBIC", "DIC", "WAIC1", "WAIC2")

emit_compare_block <- function(sub, label) {
  sub <- sub[match(model_order, sub$model), ]
  best <- sapply(crit_cols, function(cc) which.min(sub[[cc]]))
  lines <- paste0("\\multirow{4}{*}{", label, "}")
  rows <- sapply(seq_len(nrow(sub)), function(i) {
    vals <- sapply(seq_along(crit_cols), function(j) {
      v <- sprintf("%.2f", sub[[crit_cols[j]]][i])
      if (best[j] == i) paste0("\\textbf{", v, "}") else v
    })
    paste0(" & ", sub$model[i], " & ", paste(vals, collapse = " & "), " \\\\")
  })
  c(paste0(lines, rows[1]), rows[-1], "\\midrule")
}

emit_winrate_row <- function(wr_row, label) {
  vals <- sapply(model_order, function(m) sprintf("%.2f", wr_row[[m]]))
  paste0(label, " & ", paste(vals, collapse = " & "), " \\\\")
}

# -- n=600 slice: same shape as the existing draft tables --------------
tab600 <- compare_summary_tab[compare_summary_tab$n == "600", ]
wr600  <- win_rate[win_rate$n == "600", ]
lines_tab <- unlist(lapply(sort(unique(as.numeric(tab600$nu))), function(nu) {
  emit_compare_block(tab600[tab600$nu == as.character(nu), ], nu_regime_label(nu))
}))
writeLines(lines_tab, "compare_table_n600_body.tex")
lines_wr <- sapply(sort(unique(as.numeric(wr600$nu))), function(nu) {
  emit_winrate_row(wr600[wr600$nu == as.character(nu), ], nu_regime_label(nu))
})
writeLines(lines_wr, "compare_winrate_n600_body.tex")

# -- full (regime, n) grid: every cell, label includes n ----------------
lines_tab_full <- unlist(lapply(names(compare_results), function(k) {
  parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
  nu_val <- as.numeric(parts[2]); n_val <- parts[3]
  emit_compare_block(compare_summary_tab[compare_summary_tab$regime == k, ],
                      paste0(nu_regime_label(nu_val), ", $n=", n_val, "$"))
}))
writeLines(lines_tab_full, "compare_table_full_body.tex")
lines_wr_full <- sapply(names(compare_results), function(k) {
  parts <- regmatches(k, regexec("^nu=([^_]+)_n=(.+)$", k))[[1]]
  nu_val <- as.numeric(parts[2]); n_val <- parts[3]
  emit_winrate_row(win_rate[win_rate$regime == k, ],
                    paste0(nu_regime_label(nu_val), ", $n=", n_val, "$"))
})
writeLines(lines_wr_full, "compare_winrate_full_body.tex")

cat("\nWrote compare_table_n600_body.tex + compare_winrate_n600_body.tex",
    "(n=600 slice, drop-in for the existing Table 2/3 shape) and",
    "compare_table_full_body.tex + compare_winrate_full_body.tex",
    "(all 20 (regime,n) cells, for a supplementary/appendix table).\n")

# ---------------------------------------------------------------------
# 4.6 Re-plot fig_modelcomparison with the REAL win rates. Part 2 (above,
#     earlier in the script) plots this figure right after Part 1, before
#     Part 4 has run -- on a single top-to-bottom "full" run that call
#     only sees the hardcoded 3-regime fallback (`win_rate` doesn't exist
#     yet at that point), so it would otherwise silently ship the OLD
#     placeholder numbers instead of this run's actual results. Redo it
#     here, now that `win_rate` is populated, using the n=600 slice
#     (same regimes/meaning as the original figure).
# ---------------------------------------------------------------------
pdf("fig_modelcomparison.pdf", width = 7.5, height = 4.2); plot_winrate(wr600); dev.off()
png("fig_modelcomparison.png", width = 1900, height = 1050, res = 220); plot_winrate(wr600); dev.off()
cat("\nRe-wrote fig_modelcomparison.{pdf,png} using this run's actual win rates",
    "(n=600 slice) -- this supersedes the copy Part 2 wrote earlier from",
    "the fallback numbers.\n")

cat("\nDone. Part 1 fits the recovery grid; Part 2 plots off `summary_tab`",
    "(and `win_rate` once Part 4 has run); Part 3's plot_predictive()/qoi",
    "are ready for any fit with ff > 0; Part 4 compares against ZIP/ZINB/ZIGP.\n")
