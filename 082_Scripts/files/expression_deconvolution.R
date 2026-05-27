# =============================================================================
# Bayesian Gene Expression Deconvolution
# Van Loo model implemented via Stan / rstan
# =============================================================================
# Install dependencies (run once):
#   install.packages("rstan")
#   install.packages("posterior")   # tidy posterior summaries
#   install.packages("bayesplot")   # posterior / PPC plots
# =============================================================================

library(rstan)
library(posterior)
library(bayesplot)

# Use multiple cores and cache compiled models
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# -----------------------------------------------------------------------------
# 1.  Point-estimate helper (Van Loo closed-form solution)
#     Useful for sanity-checking the MCMC output.
# -----------------------------------------------------------------------------

#' Closed-form point estimates of eH and eT from a single variant
#'
#' @param eA   Integer. Reference-allele read count.
#' @param eB   Integer. Alternate-allele read count.
#' @param nAT  Numeric. Tumour copies of the A allele.
#' @param nBT  Numeric. Tumour copies of the B allele.
#' @param nAH  Numeric. Host copies of the A allele.
#' @param nBH  Numeric. Host copies of the B allele.
#' @param rho  Numeric in (0,1). Tumour purity (aberrant cell fraction).
#' @return Named numeric vector with elements eH and eT.
point_estimate <- function(eA, eB, nAT, nBT, nAH, nBH, rho) {
  denom_common <- nAT * nBH - nAH * nBT   # shared denominator term

  eH_est <- (nAT * eB - nBT * eA) /
              ((1 - rho) * (nAT * nBH - nAH * nBT))

  eT_est <- (nAH * eB - nBH * eA) /
              (rho * (nAH * nBT - nAT * nBH))

  c(eH = eH_est, eT = eT_est)
}


# -----------------------------------------------------------------------------
# 2.  Fit the Bayesian model
# -----------------------------------------------------------------------------

#' Fit the expression deconvolution Stan model
#'
#' @param eA    Integer vector length V. Reference-allele counts per variant.
#' @param eB    Integer vector length V. Alternate-allele counts per variant.
#' @param nAT   Numeric vector length V. Tumour A-allele copy numbers.
#' @param nBT   Numeric vector length V. Tumour B-allele copy numbers.
#' @param nAH   Numeric vector length V. Host A-allele copy numbers.
#' @param nBH   Numeric vector length V. Host B-allele copy numbers.
#' @param rho   Numeric vector length V. Per-variant tumour purity.
#' @param stan_file  Path to the .stan model file.
#' @param iter   Total MCMC iterations per chain (including warmup).
#' @param warmup Warmup iterations per chain.
#' @param chains Number of MCMC chains.
#' @param seed   Random seed for reproducibility.
#' @return A stanfit object.
fit_expression_model <- function(eA, eB,
                                 nAT, nBT, nAH, nBH,
                                 rho,
                                 stan_file = "expression_deconvolution.stan",
                                 iter   = 2000,
                                 warmup = 1000,
                                 chains = 4,
                                 seed   = 42) {

  stopifnot(
    length(eA) == length(eB),
    length(eA) == length(nAT),
    all(c(nAT, nBT, nAH, nBH) >= 0),
    all(rho > 0 & rho < 1)
  )

  stan_data <- list(
    V   = length(eA),
    eA  = as.integer(eA),
    eB  = as.integer(eB),
    nAT = nAT,
    nBT = nBT,
    nAH = nAH,
    nBH = nBH,
    rho = rho
  )

  fit <- stan(
    file    = stan_file,
    data    = stan_data,
    iter    = iter,
    warmup  = warmup,
    chains  = chains,
    seed    = seed,
    control = list(adapt_delta = 0.95)   # reduce divergences
  )

  fit
}


# -----------------------------------------------------------------------------
# 3.  Summarise and plot results
# -----------------------------------------------------------------------------

#' Print a tidy posterior summary for eH, eT, and phi
summarise_fit <- function(fit) {
  draws <- as_draws_df(fit)

  summary_tbl <- summarise_draws(
    draws,
    mean, median, sd,
    ~quantile(.x, probs = c(0.025, 0.975)),
    rhat, ess_bulk
  )

  params_of_interest <- summary_tbl[
    summary_tbl$variable %in% c("eH", "eT", "phi"), ]

  print(params_of_interest)
  invisible(params_of_interest)
}

#' Trace plots for the three main parameters
plot_traces <- function(fit) {
  posterior_array <- as.array(fit, pars = c("eH", "eT", "phi"))
  mcmc_trace(posterior_array)
}

#' Posterior density plots for eH and eT
plot_posteriors <- function(fit) {
  posterior_array <- as.array(fit, pars = c("eH", "eT"))
  mcmc_areas(posterior_array, prob = 0.95)
}

#' Posterior predictive check — observed vs replicated counts
plot_ppc <- function(fit, eA_obs, eB_obs) {
  draws_df <- as.data.frame(fit)

  # Extract replicated A counts (eA_rep[1] ... eA_rep[V])
  eA_rep_cols <- grep("^eA_rep\\[", names(draws_df), value = TRUE)
  eA_rep_mat  <- as.matrix(draws_df[, eA_rep_cols])

  ppc_dens_overlay(y = eA_obs, yrep = eA_rep_mat[1:50, ]) +
    ggtitle("PPC: reference-allele counts (eA)")
}


# =============================================================================
# 4.  Example usage with simulated data
# =============================================================================

set.seed(123)

# True parameter values used to simulate data
eH_true  <- 150   # host expression per copy
eT_true  <- 300   # tumour expression per copy
phi_true <- 20    # overdispersion

V   <- 10                          # number of variants
rho <- rep(0.7, V)                 # constant purity across variants
nAT <- sample(1:3, V, replace = TRUE)
nBT <- sample(0:2, V, replace = TRUE)
nAH <- rep(1, V)                   # diploid host, A allele
nBH <- rep(1, V)                   # diploid host, B allele

# Simulate expected means
muA <- rho * nAT * eT_true + (1 - rho) * nAH * eH_true
muB <- rho * nBT * eT_true + (1 - rho) * nBH * eH_true

# Simulate observed counts from negative binomial
eA_sim <- rnbinom(V, mu = muA, size = phi_true)
eB_sim <- rnbinom(V, mu = muB, size = phi_true)

cat("Simulated data summary:\n")
cat("  eA:", eA_sim, "\n")
cat("  eB:", eB_sim, "\n\n")

# --- Per-variant point estimates (closed form) ---
cat("Closed-form point estimates per variant:\n")
for (v in seq_len(V)) {
  pe <- point_estimate(eA_sim[v], eB_sim[v],
                       nAT[v], nBT[v], nAH[v], nBH[v], rho[v])
  cat(sprintf("  Variant %2d:  eH = %7.2f,  eT = %7.2f\n",
              v, pe["eH"], pe["eT"]))
}
cat(sprintf("\n  True values: eH = %.2f,  eT = %.2f\n\n",
            eH_true, eT_true))

# --- Bayesian fit ---
# Uncomment the block below to run MCMC (requires rstan to be installed
# and the .stan file to be present in the working directory):

fit <- fit_expression_model(
  eA  = eA_sim, eB  = eB_sim,
  nAT = nAT,    nBT = nBT,
  nAH = nAH,    nBH = nBH,
  rho = rho,
  stan_file = "~/Downloads/files/expression_deconvolution.stan"
)

summarise_fit(fit)
plot_traces(fit)
plot_posteriors(fit)
plot_ppc(fit, eA_sim, eB_sim)

fitAdrian <- fit_expression_model(
  eA  = adrianDf$eA_nrm[1:1000], 
  eB  = adrianDf$eB_nrm[1:1000], 
  nAT = adrianDf$nAT[1:1000], 
  nBT = adrianDf$nBT[1:1000], 
  nAH = adrianDf$nAH[1:1000], 
  nBH = adrianDf$nBH[1:1000], 
  rho = adrianDf$p[1:1000], 
  stan_file = "~/Downloads/files/expression_deconvolution.stan"
)

summarise_fit(fitAdrian)
plot_traces(fitAdrian)
plot_posteriors(fitAdrian)
plot_ppc(fitAdrian, eA_sim, eB_sim)

# > test <- point_estimate(41.622483, 6.091095, 3, 0, 1, 1, 0.5958787)
# > test
# eH       eT 
# 15.07244 19.87619 --> Adrian's table contains point estimates?

