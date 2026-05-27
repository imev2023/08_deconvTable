// Bayesian gene expression deconvolution model
// Based on Van Loo equations, as described in the thesis
//
// For each variant v = 1..V:
//   eAv ~ NegBinomial(muAv, phi)
//   eBv ~ NegBinomial(muBv, phi)
//
//   muAv = rho_v * nATv * eT + (1 - rho_v) * nAHv * eH
//   muBv = rho_v * nBTv * eT + (1 - rho_v) * nBHv * eH
//
// Priors (half-normal on [0, inf)):
//   eH ~ Normal(0, 1e3)  T+
//   eT ~ Normal(0, 1e3)  T+
//   phi ~ Normal(0, 1e3) T+

data {
  int<lower=1> V;             // Number of variants

  // Observed allele read counts (integers >= 0)
  array[V] int<lower=0> eA;   // Reference allele counts
  array[V] int<lower=0> eB;   // Alternate allele counts

  // Copy numbers (can be non-integer due to subclonality estimates)
  vector<lower=0>[V] nAT;     // Tumour copies of A allele
  vector<lower=0>[V] nBT;     // Tumour copies of B allele
  vector<lower=0>[V] nAH;     // Host copies of A allele
  vector<lower=0>[V] nBH;     // Host copies of B allele

  // Per-variant tumour purity
  vector<lower=0, upper=1>[V] rho;  // Aberrant cell fraction
}

parameters {
  real<lower=0> eH;   // Host expression per gene copy
  real<lower=0> eT;   // Tumour expression per gene copy
  real<lower=0> phi;  // Overdispersion parameter (shared across variants)
}

model {
  // Weak half-normal priors
  eH  ~ normal(0, 1e3);
  eT  ~ normal(0, 1e3);
  phi ~ normal(0, 1e3);

  // Likelihood

for (v in 1:V) {
    real muAv = fmax(rho[v] * nAT[v] * eT + (1 - rho[v]) * nAH[v] * eH, 1e-9);
    real muBv = fmax(rho[v] * nBT[v] * eT + (1 - rho[v]) * nBH[v] * eH, 1e-9);

    eA[v] ~ neg_binomial_2(muAv, phi);
    eB[v] ~ neg_binomial_2(muBv, phi);
}
}

generated quantities {
  // Posterior predictive samples (useful for model checking)
  array[V] int eA_rep;
  array[V] int eB_rep;

  for (v in 1:V) {
    real muAv = rho[v] * nAT[v] * eT + (1 - rho[v]) * nAH[v] * eH;
    real muBv = rho[v] * nBT[v] * eT + (1 - rho[v]) * nBH[v] * eH;

    eA_rep[v] = neg_binomial_2_rng(muAv, phi);
    eB_rep[v] = neg_binomial_2_rng(muBv, phi);
  }
}
