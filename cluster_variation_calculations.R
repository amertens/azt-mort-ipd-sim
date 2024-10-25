# Function to estimate CV from reported estimates
estimate_cv <- function(
    point_est,      # Point estimate of mortality rate
    ci_lower,       # Lower 95% CI
    ci_upper,       # Upper 95% CI
    n_clusters,     # Number of clusters
    method = "simple"  # Calculation method
) {
  
  if(method == "simple") {
    # Simple method using CI width
    cv <- (ci_upper - ci_lower)/(4 * point_est)
    
  } else if(method == "nb") {
    # Method based on negative binomial model
    # For negative binomial: Var(Y) = μ + μ^2/r where r is size parameter
    # CV^2 = Var(Y)/μ^2 = 1/μ + 1/r
    
    # Estimate standard error from CI
    se <- (ci_upper - ci_lower)/(2 * 1.96)
    
    # Variance of log rate
    var_log <- (se/point_est)^2
    
    # Account for cluster sampling
    var_between <- var_log - 1/(n_clusters * point_est)
    
    # Convert to CV
    cv <- sqrt(exp(var_between) - 1)
    
  } else if(method == "bootstrap") {
    # Simulation-based approach
    # Generate many datasets with different CVs
    # Find CV that produces similar CI width
    
    cvs <- seq(0.1, 1, by=0.1)
    matches <- sapply(cvs, function(cv) {
      ci_widths <- replicate(1000, {
        # Generate cluster effects
        cluster_effects <- rnorm(n_clusters, 0, sqrt(log(cv^2 + 1)))
        rates <- point_est * exp(cluster_effects)
        # Calculate CI width
        quantile(rates, 0.975) - quantile(rates, 0.025)
      })
      mean(ci_widths)
    })
    
    # Find CV that produces closest match to observed CI width
    cv <- cvs[which.min(abs(matches - (ci_upper - ci_lower)))]
  }
  
  return(cv)
}

# Calculate CVs for each study
studies <- list(
  AVENIR = list(
    rate = 13.9,
    ci_lower = 13.0,
    ci_upper = 14.8,
    n_clusters = 1455
  ),
  TANA = list(
    rate = 8.3,
    ci_lower = 5.3, 
    ci_upper = 13.1,
    n_clusters = 24
  ),
  MORDOR_II = list(
    rate = 26.3,
    ci_lower = 24.2,
    ci_upper = 28.8,
    n_clusters = 594
  )
)

# Calculate CV estimates using different methods
cv_estimates <- lapply(names(studies), function(study) {
  params <- studies[[study]]
  
  cvs <- sapply(c("simple", "nb", "bootstrap"), function(method) {
    estimate_cv(
      point_est = params$rate,
      ci_lower = params$ci_lower,
      ci_upper = params$ci_upper,
      n_clusters = params$n_clusters,
      method = method
    )
  })
  
  names(cvs) <- c("simple", "nb", "bootstrap")
  return(cvs)
})
names(cv_estimates) <- names(studies)

# Print results
print(cv_estimates)

# Compare to MORDOR reported CV
cat("\nMORDOR reported CV: 0.34\n")