
rm(list=ls())
source("sim_data_functions.R")
source("study_params.R")
params <- study_params$MORDOR_Niger


n_clusters = params$n_clusters
cluster_size = params$cluster_size
baseline_rate = params$baseline_rate
#baseline_rate = 0.0275/2
effect_size = params$effect_size
cv = params$cv
study_name = "MORDOR_Niger"

       
    study_duration = 2
    loss_rate = 0.10   
    study_name
    region = "rural"
    # New interaction effect parameters
    age_interaction = 0.2     # Stronger effect in younger children
    distance_interaction = 0.1 # Stronger effect in remote areas
    wealth_interaction = -0.1  # Weaker effect in wealthy households
    season_interaction = 0.15   # Stronger effect in wet season

  
  # Calculate between-cluster variance
  sigma2 <- log(cv^2 + 1)
  
  # Generate cluster-level random effects
  cluster_effects <- rnorm(n_clusters * 2, mean = -sigma2/2, sd = sqrt(sigma2))
  
  # Generate individual-level data
  n_total <- n_clusters * 2 * cluster_size
  
  # Generate covariates
  set.seed(2334)
  covariates <- generate_covariates(n_total, region)
  
  # Create base data1 frame
  data1 <- data.frame(
    study = study_name,
    cluster_id = rep(1:(n_clusters * 2), each = cluster_size),
    treatment = rep(rep(c(0,1), each = n_clusters), each = cluster_size),
    cluster_effect = rep(cluster_effects, each = cluster_size)
  )

    data1$treatment <- 1
  
  
  # Add covariates
  data1 <- cbind(data1, covariates)
  
  set.seed(2334)
  # Calculate covariate effects on log mortality rate
  log_rate1 <- with(data1, {
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
  data1$log_offset <- log(followup_time)
  
  # Generate death counts
  lambda <- exp(log_rate1 + data1$log_offset)
  phi <- cv^2  
  data1$deaths <- rnbinom(n_total, mu = lambda, size = 1/phi)
  data1$deaths <- ifelse(data1$deaths>0,1,0)
  prop.table(table(data1$deaths))
  
  set.seed(2334)
  # Create base data frame
  data0 <- data.frame(
    study = study_name,
    cluster_id = rep(1:(n_clusters * 2), each = cluster_size),
    treatment = rep(rep(c(0,1), each = n_clusters), each = cluster_size),
    cluster_effect = rep(cluster_effects, each = cluster_size)
  )

    data0$treatment <- 0
  
  
  # Add covariates
  data0 <- cbind(data0, covariates)
  
  # Calculate covariate effects on log mortality rate
  log_rate0 <- with(data0, {
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
  data0$log_offset <- log(followup_time)
  
  # Generate death counts
  lambda <- exp(log_rate0 + data0$log_offset)
  phi <- cv^2  
  data0$deaths <- rnbinom(n_total, mu = lambda, size = 1/phi)
  data0$deaths <- ifelse(data0$deaths>0,1,0)
  
  data0$A=0
  data1$A=1
  
  summary(exp(log_rate1-log_rate0))
  
  prop.table(table(data1$deaths))
  prop.table(table(data0$deaths))
  

  res  <- bind_rows(data0,data1)
  tab <- table(res$deaths,res$A)/100
  tab[2,2]*tab[1,1]/tab[1,2]/tab[2,1]
  