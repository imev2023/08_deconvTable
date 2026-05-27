pclamp <- function(p, eps=1e-3) {
  pmin(1-eps, pmax(eps, p))
}

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

genotyper <- function(host_ploidy,
                      tumour_ploidy,
                      host_observed_alt_reads,
                      host_observed_total_reads,
                      tumour_observed_alt_reads,
                      tumour_observed_total_reads,
                      tumour_purity,
                      eps = 1e-3) {
  
  # Local purity adjusted by relative copy numbers
  wa <- (1 - tumour_purity) * host_ploidy
  wb <- tumour_purity * tumour_ploidy
  w <- wb / (wa + wb) # w := tumour read fraction at this locus
  
  # likelihood function
  # probability of each combination of tumour and host alt allele count
  # g_host := host alt allele count
  # g_tumour := tumour alt allele count
  ll <- function(g_host, g_tumour) {
    
    stopifnot(g_host >= 0 & g_host <= host_ploidy)
    stopifnot(g_tumour >= 0 & g_tumour <= tumour_ploidy)
    
    p_host_raw <- g_host / host_ploidy
    p_tumour_raw <- g_tumour / tumour_ploidy
    
    p_host <- pclamp(p_host_raw, eps)
    
    p_mixture_raw <- (1 - w) * p_host_raw +
      w * p_tumour_raw
    
    p_mixture <- pclamp(p_mixture_raw, eps)
    
    logLikA <- dbinom(host_observed_alt_reads,
                      host_observed_total_reads,
                      prob = p_host,
                      log = TRUE)
    
    logLikB <- dbinom(tumour_observed_alt_reads,
                      tumour_observed_total_reads,
                      prob = p_mixture,
                      log = TRUE)
    
    logLikA + logLikB
  }
  
  # put all genotype combinations into a table
  grid <- expand.grid(
    g_host = seq(0, host_ploidy),
    g_tumour = seq(0, tumour_ploidy)
  )
  
  # compute likelihood of each combination
  grid$loglik <- mapply(ll, grid$g_host, grid$g_tumour)
  
  # normalise the likelihoods to get "posterior prob" (no explicit prior)
  logZ <- logsumexp(grid$loglik) # normalizing factor
  grid$posterior <- exp(grid$loglik - logZ) # prob of each grid parameter given log likelihoods 
  
  # MAP := outcome with maximum posterior prob
  map_idx <- which.max(grid$loglik) # maximum probability = most likely scenario. --> will max posterior always == max loglikelihood?
  map <- grid[map_idx, ]
# 
#   if(which.max(grid$posterior)!=which.max(grid$loglik)){
#     print(paste0("Check this combination of variables, "))
#   }
#   
  # expected values := probability-weighted average of genotypes
  expected_host_alt_alleles <- sum(grid$g_host * grid$posterior)
  expected_tumour_alt_alleles <- sum(grid$g_tumour * grid$posterior)
  
  list(
    ll = ll,
    grid = grid,
    logZ = logZ,
    w = w,
    params = list(
      host_ploidy = host_ploidy,
      tumour_ploidy = tumour_ploidy,
      tumour_purity = tumour_purity,
      eps = eps
    ),
    map = map,
    E = list(g_host = expected_host_alt_alleles, 
             g_tumour = expected_tumour_alt_alleles)
  )
}

# Example (host alt = 1/2, tumour alt = 2/2)
# est <- genotyper(2, 2, 10, 30, 51, 60, 0.4985)
# est

genotyper_somatic <- function(host_ploidy,
                              tumour_ploidy,
                              tumour_observed_alt_reads,
                              tumour_observed_total_reads,
                              tumour_purity,
                              eps = 1e-3) {
  
  # --- sanity ---
  if (!is.finite(tumour_ploidy) || tumour_ploidy < 0) {
    stop("Invalid tumour ploidy")
  }
  
  tumour_ploidy <- as.integer(round(tumour_ploidy))
  
  # --- mixture weight ---
  wa <- (1 - tumour_purity) * host_ploidy
  wb <- tumour_purity * tumour_ploidy
  w <- wb / (wa + wb)
  
  # host is FIXED: no alt alleles
  p_host_raw <- 0
  
  # --- likelihood over tumour genotypes only ---
  g_tumour_vals <- seq(0, tumour_ploidy)
  
  loglik <- sapply(g_tumour_vals, function(g_tumour) {
    
    p_tumour_raw <- g_tumour / tumour_ploidy
    
    p_mixture_raw <- (1 - w) * p_host_raw + w * p_tumour_raw
    p_mixture <- pclamp(p_mixture_raw, eps)
    
    dbinom(tumour_observed_alt_reads,
           tumour_observed_total_reads,
           prob = p_mixture,
           log = TRUE)
  })
  
  # --- normalize ---
  logZ <- logsumexp(loglik)
  posterior <- exp(loglik - logZ)
  
  # --- MAP ---
  map_idx <- which.max(loglik)
  
  # --- expectation ---
  expected_tumour_alt <- sum(g_tumour_vals * posterior)
  
  list(
    grid = data.frame(
      g_host = 0,
      g_tumour = g_tumour_vals,
      loglik = loglik,
      posterior = posterior
    ),
    logZ = logZ,
    w = w,
    map = data.frame(
      g_host = 0,
      g_tumour = g_tumour_vals[map_idx],
      loglik = loglik[map_idx],
      posterior = posterior[map_idx]
    ),
    E = list(
      g_host = 0,
      g_tumour = expected_tumour_alt
    ),
    params = list(
      host_ploidy = host_ploidy,
      tumour_ploidy = tumour_ploidy,
      tumour_purity = tumour_purity,
      eps = eps
    )
  )
}
