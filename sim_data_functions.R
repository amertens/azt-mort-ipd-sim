
# params = study_params$MORDOR_Niger
# n_clusters = params$n_clusters
# cluster_size = params$cluster_size
# baseline_rate = params$baseline_rate
# effect_size = params$effect_size
# cv = params$cv
# study_duration = 2  
# loss_rate = 0.10 
# study_name = "MORDOR_Niger"
# region = "rural"

#-------------------------------------------------------------------------------
# Data sim function
#-------------------------------------------------------------------------------

simulate_trial_data <- function(
    n_clusters,           
    cluster_size,         
    baseline_rate,        
    effect_size,          
    cv,                   
    study_duration = 2,   
    loss_rate = 0.10,    
    study_name,
    region = "rural",
    # # New interaction effect parameters -manually programmed-need to update
    # age_interaction = 0.2,     # Stronger effect in younger children
    # distance_interaction = 0.1, # Stronger effect in remote areas
    # wealth_interaction = -0.1,  # Weaker effect in wealthy households
    # season_interaction = 0.15,   # Stronger effect in wet season
    A=NULL
){
  
  # Calculate between-cluster variance
  sigma2 <- log(cv^2 + 1)
  
  # Generate cluster-level random effects
  cluster_effects <- rnorm(n_clusters * 2, mean = -sigma2/2, sd = sqrt(sigma2))
  
  # Generate individual-level data
  n_total <- n_clusters * 2 * cluster_size
  
  # Generate covariates
  covariates <- generate_covariates(n_total, region)
  
  # Create base data frame
  data <- data.frame(
    study = study_name,
    cluster_id = rep(1:(n_clusters * 2), each = cluster_size),
    treatment = rep(rep(c(0,1), each = n_clusters), each = cluster_size),
    cluster_effect = rep(cluster_effects, each = cluster_size)
  )
  if(!is.null(A)){
    data$treatment <- A
  }
  
  # Add covariates
  data <- cbind(data, covariates)
  
  # Define relative effect modifications (will be scaled)
  effect_mods <- with(data, {
    # Age effects (relative to 24-59mo)
    age_mod <- case_when(
      age_group == "1-5mo"  ~ 1.5,    # 50% stronger effect
      age_group == "6-11mo" ~ 1.3,    # 30% stronger effect
      age_group == "12-23mo" ~ 1.1,   # 10% stronger effect
      TRUE ~ 1.0                      # Reference group
    )
    
    # Sex effect
    sex_mod <- ifelse(sex == 1, 1.1, 1.0)  # 10% stronger in males
    
    # Distance effect (standardize to mean 1)
    dist_mod <- log(distance)/mean(log(distance))
    dist_mod <- 1 + 0.2 * (dist_mod - 1)  # 20% variation by distance
    
    # Multiply modifications together
    total_mod <- age_mod * sex_mod * dist_mod
    
    # Scale to maintain overall effect
    scale_factor <- mean(total_mod)
    total_mod/scale_factor
  })
  

  
  # Calculate covariate effects on log hazard
  log_hazard  <- with(data, {
    # Base rate
    log(baseline_rate) +
      
      # Main effects
      # Age effects (reference: 24-59mo)
      1.2 * (age_group == "1-5mo") +    
      0.8 * (age_group == "6-11mo") +   
      0.4 * (age_group == "12-23mo") +
      
      # Sex effect (reference: female)
      0.1 * (sex == 1) +                
      
      # Stunting effect
      0.4 * stunted +                   
      
      # Distance effect (log-linear)
      0.1 * log(distance) +             
      
      # Wealth effect (reference: quintile 3)
      0.3 * (wealth_quintile == 1) +    
      0.2 * (wealth_quintile == 2) +
      0 * (wealth_quintile == 3) +
      -0.2 * (wealth_quintile == 4) +
      -0.3 * (wealth_quintile == 5) +
      
      # Season effect
      0.2 * season +                    
      
      # Treatment  effect with interactions
      log(effect_size) * treatment * effect_mods +
      
      
      # Cluster effect
      cluster_effect
  })
  
  # Generate survival times from exponential distribution
  hazard <- exp(log_hazard)
  data$survival_time <- rexp(nrow(data), rate = hazard)
  
  # Generate censoring times
  data$censoring_time <- rexp(nrow(data), rate = loss_rate)
  
  # Calculate observed time and status
  data$time <- pmin(data$survival_time, data$censoring_time, study_duration)
  data$status <- as.numeric(data$survival_time <= pmin(data$censoring_time, study_duration))
  
  # Add time periods for time-varying analysis if needed
  data$period <- cut(data$time, 
                     breaks = seq(0, study_duration, by = 0.5),  # 6-month intervals
                     labels = FALSE)
  
  # Calculate person-time for rate calculations
  data$person_time <- data$time
  data$log_offset <- log(data$person_time)
  
  # For compatibility with existing code, add deaths column
  data$dead <- data$status
  #mean(data$deaths)
    
  return(data)
}


# Function to generate base cluster and treatment data
generate_base_data <- function(n_clusters, cluster_size, study_name) {
  
  # Calculate between-cluster variance on log scale
  sigma2 <- log(cv^2 + 1)
  
  # Generate cluster-level random effects
  cluster_effects <- rnorm(n_clusters * 2, mean = -sigma2/2, sd = sqrt(sigma2))
  
  # Calculate total sample size
  n_total <- n_clusters * 2 * cluster_size
  
  # Create base data frame
  data <- data.frame(
    study = study_name,
    cluster_id = rep(1:(n_clusters * 2), each = cluster_size),
    treatment = rep(rep(c(0,1), each = n_clusters), each = cluster_size),
    cluster_effect = rep(cluster_effects, each = cluster_size)
  )
  
  return(data)
}

# And here's the complementary generate_outcomes() function that was implied:
generate_outcomes <- function(data, log_rate, cv, study_duration, loss_rate) {
  
  # Generate followup time
  n_total <- nrow(data)
  followup_time <- rexp(n_total, rate = loss_rate)
  followup_time[followup_time > study_duration] <- study_duration
  data$log_offset <- log(followup_time)
  
  # Generate expected counts
  lambda <- exp(log_rate + data$log_offset)
  
  # Generate observed deaths using negative binomial
  phi <- cv^2  # Approximate relationship between CV and NB dispersion
  data$deaths <- rnbinom(n_total, mu = lambda, size = 1/phi)
  data$deaths <- ifelse(data$deaths==0,0,1)
  return(data)
}






# Function to generate baseline covariates
generate_covariates <- function(n, region="rural") {
  
  # Age distribution based on MORDOR/AVENIR
  age_group <- sample(c("1-5mo", "6-11mo", "12-23mo", "24-59mo"), 
                      size=n, replace=TRUE,
                      prob=c(0.06, 0.10, 0.23, 0.61))
  
  # Convert to continuous months for modeling
  age_months <- case_when(
    age_group == "1-5mo" ~ runif(n, 1, 5),
    age_group == "6-11mo" ~ runif(n, 6, 11),
    age_group == "12-23mo" ~ runif(n, 12, 23),
    age_group == "24-59mo" ~ runif(n, 24, 59)
  )
  
  # Sex distribution (~51% male across studies)
  sex <- rbinom(n, 1, 0.51)
  
  # Stunting status - using typical rates in rural sub-Saharan Africa
  # WHO estimates ~30% stunting prevalence
  stunted <- rbinom(n, 1, 0.30)
  
  # Distance to healthcare facility (km)
  # Log-normal distribution with median ~5km for rural areas
  if(region == "rural") {
    distance <- rlnorm(n, log(5), 0.5) 
  } else {
    distance <- rlnorm(n, log(2), 0.5)
  }
  
  # Wealth index - generate quintiles
  wealth_score <- rnorm(n)
  wealth_quintile <- cut(wealth_score, 
                         breaks=quantile(wealth_score, probs=seq(0,1,0.2)),
                         labels=1:5,
                         include.lowest=TRUE)
  
  # Season - assign based on enrollment timing
  # Simplified to wet vs dry season
  season <- rbinom(n, 1, 0.5)
  
  # Return as data frame
  data.frame(
    age_months = age_months,
    age_group = age_group,
    sex = sex,
    stunted = stunted,
    distance = distance,
    wealth_quintile = wealth_quintile,
    season = season
  )
}



check_simulated_variation <- function(sim_data) {
  # Calculate empirical CVs from simulated data
  empirical_cvs <- by(sim_data, sim_data$study, function(x) {
    # Calculate CV of mortality rates
    rates <- x$deaths/exp(x$log_offset)
    sd(rates)/mean(rates)
  })
  
  # Compare to input CVs
  cv_comparison <- data.frame(
    study = names(study_params),
    input_cv = sapply(study_params, `[[`, "cv"),
    simulated_cv = as.numeric(empirical_cvs)
  )
  
  return(cv_comparison)
}





