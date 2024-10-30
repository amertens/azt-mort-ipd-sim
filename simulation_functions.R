




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



# Function to analyze subgroups for a single study
analyze_subgroups <- function(data) {
  
  # Main effect
  main <- analyze_trial_data(data, adjusted=TRUE)
  names(main) <- c("coef", "se" )
  
  # Age subgroups - analyze each separately
  age_groups <- c("1-5mo", "6-11mo", "12-23mo", "24-59mo")
  age_results <- lapply(age_groups, function(age) {
    sub_data <- subset(data, age_group == age)
    fit <- try({analyze_trial_data(sub_data, adjusted=TRUE, age_subgroup=TRUE)}, silent=TRUE)
    
    if(!inherits(fit, "try-error")) {
      names(fit) <- c("coef", "se" )
      fit
    } else {
      c(coef = NA, se = NA)
    }
  })
  names(age_results) <- paste0("age_", age_groups)
  
  # Sex subgroups
  sex_results <- lapply(c(0,1), function(sex_val) {
    sub_data <- subset(data, sex == sex_val)
    fit <- try({analyze_trial_data(sub_data, adjusted=TRUE, sex_subgroup=TRUE)}, silent=TRUE)
    names(fit) <- c("coef", "se" )
    fit
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
    coefs <- sapply(study_results, function(x) x[paste0(name, ".coef")])
    ses <- sapply(study_results, function(x) x[paste0(name, ".se")])
    
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














# Function to test for treatment interactions
test_interactions <- function(data) {
  
  # Fit models with interactions
  fits <- list(
    # Age interaction
    age = MASS::glm.nb(deaths ~ treatment * age_group + sex + stunted +
                         log(distance) + wealth_quintile + season + 
                         offset(log_offset), data = data),
    
    # Sex interaction  
    sex = MASS::glm.nb(deaths ~ treatment * sex + age_group + stunted +
                         log(distance) + wealth_quintile + season + 
                         offset(log_offset), data = data),
    
    # Distance interaction
    distance = MASS::glm.nb(deaths ~ treatment * log(distance) + age_group + 
                              sex + stunted + wealth_quintile + season + 
                              offset(log_offset), data = data)
  )
  
  # Extract interaction tests
  interaction_tests <- lapply(fits, function(fit) {
    # Get ANOVA comparison
    null_fit <- update(fit, . ~ . - treatment:age_group)
    test <- anova(null_fit, fit)
    
    # Return test results
    c(
      chi_sq = test$Deviance[2],
      df = test$Df[2],
      p_value = test$`Pr(>Chi)`[2]
    )
  })
  
  return(interaction_tests)
}

# Function to calculate required sample size for subgroup
calculate_subgroup_power <- function(
    baseline_rate,
    effect_size,
    cv,
    cluster_size,
    alpha = 0.05,
    power = 0.80,
    subgroup_prop = 1.0  # Proportion of sample in subgroup
) {
  
  # Function to calculate power for given n_clusters
  power_for_n <- function(n) {
    # Simplified calculation based on cluster-randomized design
    effective_n <- n * cluster_size * subgroup_prop
    sigma2 <- log(cv^2 + 1)
    var_effect <- 2 * (1/effective_n + sigma2)
    
    # Calculate power
    crit_value <- qnorm(1 - alpha/2)
    z_power <- abs(log(effect_size))/sqrt(var_effect) - crit_value
    pnorm(z_power)
  }
  
  # Binary search for required n_clusters
  n_low <- 2
  n_high <- 10000
  
  while(n_high - n_low > 1) {
    n_mid <- floor((n_low + n_high)/2)
    achieved_power <- power_for_n(n_mid)
    
    if(achieved_power < power) {
      n_low <- n_mid
    } else {
      n_high <- n_mid
    }
  }
  
  return(n_high)  # Return conservative estimate
}


# Function to analyze trial data with clustered SEs
analyze_trial_data <- function(data, adjusted=TRUE, age_subgroup=FALSE, sex_subgroup=FALSE){
  
  # Fit model
  if(adjusted){
    # Adjusted analysis
    if(age_subgroup){
      fit <- MASS::glm.nb(deaths ~ treatment + sex + stunted +
                            log(distance) + wealth_quintile + season + 
                            offset(log_offset), data = data)
    }
    if(sex_subgroup){
      fit <- MASS::glm.nb(deaths ~ treatment + age_group  + stunted +
                            log(distance) + wealth_quintile + season + 
                            offset(log_offset), data = data)
    }
    
    if(!age_subgroup & !sex_subgroup){
      fit <- MASS::glm.nb(deaths ~ treatment + age_group + sex + stunted +
                            log(distance) + wealth_quintile + season + 
                            offset(log_offset), data = data)
    }
    
  } else {
    # Unadjusted analysis
    fit <- MASS::glm.nb(deaths ~ treatment + offset(log_offset), 
                        data = data)
  }
  
  # Get clustered variance-covariance matrix
  cluster_vcov <- sandwich::vcovCL(fit, cluster = data$cluster_id, type = "HC1")
  
  # Extract treatment effect and clustered SE
  coef <- coef(fit)["treatment"]
  se <- sqrt(cluster_vcov["treatment", "treatment"])
  
  return(c(coef = coef, se = se))
}