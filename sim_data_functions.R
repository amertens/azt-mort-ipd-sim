



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
    # New interaction effect parameters
    age_interaction = 0.2,     # Stronger effect in younger children
    distance_interaction = 0.1, # Stronger effect in remote areas
    wealth_interaction = -0.1,  # Weaker effect in wealthy households
    season_interaction = 0.15,   # Stronger effect in wet season
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
  
  # Calculate covariate effects on log mortality rate
  log_rate <- with(data, {
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
      
      # Treatment main effect
      log(effect_size) * treatment +
      
      # Treatment interactions
      treatment * (
        # Age interaction (stronger effect in younger)
        age_interaction * ((age_group == "1-5mo") + 
                             0.75 * (age_group == "6-11mo") +
                             0.5 * (age_group == "12-23mo")) #+
          
          # # Distance interaction (stronger effect with distance)
          # distance_interaction * (log(distance) - mean(log(distance))) +
          # 
          # # Wealth interaction (weaker effect in wealthy)
          # wealth_interaction * (as.numeric(wealth_quintile) - 3) +
          # 
          # # Season interaction (stronger in wet season)
          # season_interaction * season
      ) +
      
      # Cluster effect
      cluster_effect
  })
  
  # Generate observed person-time accounting for loss to follow-up
  followup_time <- rexp(n_total, rate = loss_rate)
  followup_time[followup_time > study_duration] <- study_duration
  data$log_offset <- log(followup_time)
  
  # Generate death counts
  lambda <- exp(log_rate + data$log_offset)
  phi <- cv^2  
  data$deaths <- rnbinom(n_total, mu = lambda, size = 1/phi)
  
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






