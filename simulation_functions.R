




# Function to analyze meta-analysis
analyze_meta_analysis <- function(studies_data,  method = "REML", adjusted=TRUE) {
  
  # Get results for each study
  study_results <- lapply(split(studies_data, studies_data$study), function(x) analyze_trial_data(data=x, adjusted=adjusted))
  
  # Combine into matrix
  results_matrix <- do.call(rbind, study_results)
  
  # Random effects meta-analysis
  meta_fit <- rma(yi = results_matrix[,"coef.treatment"], 
                  sei = results_matrix[,"se.treatment"],
                  method = method)
  
  
  return(list(
    individual_results = results_matrix,
    meta_results = meta_fit
  ))
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

# Update trial simulation to include prognostic effects
simulate_trial_data <- function(
    n_clusters,           
    cluster_size,         
    baseline_rate,        
    effect_size,          
    cv,                   
    study_duration = 2,   
    loss_rate = 0.10,    
    study_name,
    region = "rural"){
  
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
  
  # Add covariates
  data <- cbind(data, covariates)
  
  # Calculate covariate effects on log mortality rate
  # Based on literature and trial subgroup analyses
  log_rate <- with(data, {
    # Base rate
    log(baseline_rate) +
      
      # Age effects (reference: 24-59mo)
      1.2 * (age_group == "1-5mo") +    # Highest risk in youngest
      0.8 * (age_group == "6-11mo") +   
      0.4 * (age_group == "12-23mo") +
      
      # Sex effect (reference: female)
      0.1 * (sex == 1) +                # Slightly higher male mortality
      
      # Stunting effect
      0.4 * stunted +                   # Higher mortality if stunted
      
      # Distance effect (log-linear)
      0.1 * log(distance) +             # Higher mortality with distance
      
      # Wealth effect (reference: quintile 3)
      0.3 * (wealth_quintile == 1) +    # Higher mortality in poorest
      0.2 * (wealth_quintile == 2) +
      0 * (wealth_quintile == 3) +
      -0.2 * (wealth_quintile == 4) +
      -0.3 * (wealth_quintile == 5) +
      
      # Season effect
      0.2 * season +                    # Higher in wet season
      
      # Treatment effect
      log(effect_size) * treatment +
      
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

# Update analysis function to allow for covariate adjustment
analyze_trial_data <- function(data, adjusted=TRUE) {
  if(adjusted) {
    # Adjusted analysis
    fit <- MASS::glm.nb(deaths ~ treatment + age_group + sex + stunted +
                          log(distance) + wealth_quintile + season + 
                          offset(log_offset), data = data)
  } else {
    # Unadjusted analysis
    fit <- MASS::glm.nb(deaths ~ treatment + offset(log_offset), 
                        data = data)
  }
  
  # Extract treatment effect
  coef <- coef(fit)["treatment"]
  se <- sqrt(diag(vcov(fit)))["treatment"]
  
  return(c(coef = coef, se = se))
}

# Rest of simulation code remains same but can now compare
# adjusted vs unadjusted analyses



# Function to analyze subgroups for a single study
analyze_subgroups <- function(data) {
  
  # Main effect
  main <- analyze_trial_data(data, adjusted=TRUE)
  
  # Age subgroups - analyze each separately
  age_groups <- c("1-5mo", "6-11mo", "12-23mo", "24-59mo")
  age_results <- lapply(age_groups, function(age) {
    sub_data <- subset(data, age_group == age)
    fit <- try({
      MASS::glm.nb(deaths ~ treatment + sex + stunted +
                     log(distance) + wealth_quintile + season + 
                     offset(log_offset), data = sub_data)
    }, silent=TRUE)
    
    if(!inherits(fit, "try-error")) {
      c(coef = coef(fit)["treatment"],
        se = sqrt(diag(vcov(fit)))["treatment"])
    } else {
      c(coef = NA, se = NA)
    }
  })
  names(age_results) <- paste0("age_", age_groups)
  
  # Sex subgroups
  sex_results <- lapply(c(0,1), function(sex_val) {
    sub_data <- subset(data, sex == sex_val)
    fit <- MASS::glm.nb(deaths ~ treatment + age_group + stunted +
                          log(distance) + wealth_quintile + season + 
                          offset(log_offset), data = sub_data)
    c(coef = coef(fit)["treatment"],
      se = sqrt(diag(vcov(fit)))["treatment"])
  })
  names(sex_results) <- c("sex_female", "sex_male")
  
  # Return all results
  c(
    main = main,
    do.call(c, age_results),
    do.call(c, sex_results)
  )
}

# Function to meta-analyze subgroup results
meta_analyze_subgroups <- function(studies_data, method="REML") {
  
  # Split by study
  by_study <- split(studies_data, studies_data$study)
  
  # Get subgroup results for each study
  study_results <- lapply(by_study, analyze_subgroups)
  
  # Organize results into matrices for meta-analysis
  results <- list()
  
  # Names of all analyses (main + subgroups)
  analysis_names <- c("main", 
                      paste0("age_", c("1-5mo", "6-11mo", "12-23mo", "24-59mo")),
                      c("sex_female", "sex_male"))
  
  # Do meta-analysis for each type of analysis
  for(i in seq_along(analysis_names)){
    name <- analysis_names[i]
    
    # Extract coefficients and SEs for this analysis
    coefs <- sapply(study_results, function(x) x[paste0(name, ".coef.treatment")])
    ses <- sapply(study_results, function(x) x[paste0(name, ".se.treatment")])
    
    # Remove any NA results
    valid <- !is.na(coefs) & !is.na(ses)
    
    if(sum(valid) >= 2) {  # Need at least 2 studies for meta-analysis
      
      meta_fit=NULL
      try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method=method, measure="GEN"))
      if(is.null(meta_fit)){try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method="ML", measure="GEN"))}
      if(is.null(meta_fit)){try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method="DL", measure="GEN"))}
      if(is.null(meta_fit)){try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method="HE", measure="GEN"))}
      results[[name]] <- meta_fit
    } else {
      results[[name]] <- NULL
    }
  }
  
  #bind list result:
  study_res_df <- data.frame(t(do.call(rbind, study_results)))
  
  study_res_df$temp <- rownames(study_res_df)
  study_res_df$temp <- gsub(".treatment", "", study_res_df$temp)
  study_res_df$group <- str_split_i(study_res_df$temp, "\\.", 1)
  study_res_df$est <- str_split_i(study_res_df$temp, "\\.", 2)
  study_res_df_coef <- study_res_df %>% filter(est=="coef") %>% select(-c(temp,est))  
  study_res_df_se <- study_res_df %>% filter(est=="se") %>% select(-c(temp,est))    
  study_res_df_coef <- study_res_df_coef %>%
    pivot_longer(cols = -c(group), names_to = "study",values_to = "effect")
  study_res_df_se <- study_res_df_se %>%
    pivot_longer(cols = -c(group), names_to = "study",values_to = "se")
  
  study_res_df <- left_join(study_res_df_coef,study_res_df_se,by=c("group","study"))
  
  return(list(results=results, study_results=study_res_df))
}