




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







# Function to analyze subgroups for a single study
analyze_subgroups <- function(data, adjusted=TRUE, run_tmle=TRUE) {
  
  # Main effect
  main <- analyze_trial_data(data, adjusted=adjusted, run_tmle=run_tmle)

  # Age subgroup analysis
  age_results <- lapply(unique(data$age_group), function(age) {
    analyze_trial_data(data, subgroup = "age_group", subgroup_level = age, adjusted=adjusted, run_tmle=run_tmle)
  })
  names(age_results) <- unique(data$age_group)
  age_results=data.table::rbindlist(age_results, use.names = TRUE, idcol = "level") %>% mutate(group="age_group")
  
  # Sex subgroup analysis
  sex_results <- lapply(c(0,1), function(sex) {
    analyze_trial_data(data, subgroup = "sex", subgroup_level = sex, adjusted=adjusted, run_tmle=run_tmle)
  })
  names(sex_results) <- c("female", "male")
  sex_results=data.table::rbindlist(sex_results, use.names = TRUE, idcol = "level") %>% mutate(group="sex")
  
  results = bind_rows(main %>% mutate(group="main", level="main"),
                      age_results,
                      sex_results)

  # Return all results
  return(results)
}



# Function to get all estimator and subgroup results
run_analysis <- function(studies_data, adjusted=TRUE, run_tmle=TRUE) {

  full_res=NULL
  # Get subgroup results for each study
  for(i in unique(studies_data$study)){
    df <- studies_data %>% filter(study==i)
    res=NULL
    try(res <- analyze_subgroups(data=df, adjusted=adjusted, run_tmle=run_tmle))
    full_res=bind_rows(full_res,res)
  }
  
  #NEED TO DEBUG UNADJUSTED COX PH
  
  res= studies_data %>% group_by(study) %>%
    do(res=analyze_subgroups(data=., adjusted=adjusted, run_tmle=run_tmle))

  full_res=NULL
  for(i in 1:length(res$res)){
    temp=res$res[[i]]
    temp$study=res$study[i]
    full_res=bind_rows(full_res,temp)
  }
  

  return(full_res)
}

# 
# # Function to meta-analyze subgroup results
# meta_analyze_subgroups <- function(studies_data, method="REML", adjusted=TRUE, run_tmle=TRUE) {
#   
#   # Split by study
#   by_study <- split(studies_data, studies_data$study)
#   
#   # Get subgroup results for each study
#   study_results <- lapply(by_study, analyze_subgroups(adjusted=adjusted, run_tmle=run_tmle))
#   
#   # Organize results into matrices for meta-analysis
#   results <- list()
#   
#   # Names of all analyses (main + subgroups)
#   analysis_names <- c("main", 
#                       paste0("age_", c("1-5mo", "6-11mo", "12-23mo", "24-59mo")),
#                       c("sex_female", "sex_male"))
#   
#   # Do meta-analysis for each type of analysis
#   for(i in seq_along(analysis_names)){
#     name <- analysis_names[i]
#     
#     # Extract coefficients and SEs for this analysis
#     coefs <- sapply(study_results, function(x) x[paste0(name, ".coef")])
#     ses <- sapply(study_results, function(x) x[paste0(name, ".se")])
#     
#     # Remove any NA results
#     valid <- !is.na(coefs) & !is.na(ses)
#     
#     if(sum(valid) >= 2) {  # Need at least 2 studies for meta-analysis
#       
#       meta_fit=NULL
#       try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method=method, measure="GEN"))
#       if(is.null(meta_fit)){try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method="ML", measure="GEN"))}
#       if(is.null(meta_fit)){try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method="DL", measure="GEN"))}
#       if(is.null(meta_fit)){try(meta_fit<-rma(yi=coefs[valid], sei= ses[valid],  method="HE", measure="GEN"))}
#       results[[name]] <- meta_fit
#     } else {
#       results[[name]] <- NULL
#     }
#   }
#   
#   #bind list result:
#   study_res_df <- data.frame(t(do.call(rbind, study_results)))
#   
#   study_res_df$temp <- rownames(study_res_df)
#   study_res_df$temp <- gsub(".treatment", "", study_res_df$temp)
#   study_res_df$group <- str_split_i(study_res_df$temp, "\\.", 1)
#   study_res_df$est <- str_split_i(study_res_df$temp, "\\.", 2)
#   study_res_df_coef <- study_res_df %>% filter(est=="coef") %>% select(-c(temp,est))  
#   study_res_df_se <- study_res_df %>% filter(est=="se") %>% select(-c(temp,est))    
#   study_res_df_coef <- study_res_df_coef %>%
#     pivot_longer(cols = -c(group), names_to = "study",values_to = "effect")
#   study_res_df_se <- study_res_df_se %>%
#     pivot_longer(cols = -c(group), names_to = "study",values_to = "se")
#   
#   study_res_df <- left_join(study_res_df_coef,study_res_df_se,by=c("group","study"))
#   
#   return(list(results=results, study_results=study_res_df))
# }














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


# # Function to analyze trial data with clustered SEs
# analyze_trial_data <- function(data, adjusted=TRUE, age_subgroup=FALSE, sex_subgroup=FALSE){
#   
#   browser()
#   
#   # Fit model
#   if(adjusted){
#     # Adjusted analysis
#     if(age_subgroup){
#       fit <- MASS::glm.nb(deaths ~ treatment + sex + stunted +
#                             log(distance) + wealth_quintile + season + 
#                             offset(log_offset), data = data)
#     }
#     if(sex_subgroup){
#       fit <- MASS::glm.nb(deaths ~ treatment + age_group  + stunted +
#                             log(distance) + wealth_quintile + season + 
#                             offset(log_offset), data = data)
#     }
#     
#     if(!age_subgroup & !sex_subgroup){
#       fit <- MASS::glm.nb(deaths ~ treatment + age_group + sex + stunted +
#                             log(distance) + wealth_quintile + season + 
#                             offset(log_offset), data = data)
#     }
#     
#   } else {
#     # Unadjusted analysis
#     fit <- MASS::glm.nb(deaths ~ treatment + offset(log_offset), 
#                         data = data)
#   }
#   
#   # Get clustered variance-covariance matrix
#   cluster_vcov <- sandwich::vcovCL(fit, cluster = data$cluster_id, type = "HC1")
#   
#   # Extract treatment effect and clustered SE
#   coef <- coef(fit)["treatment"]
#   se <- sqrt(cluster_vcov["treatment", "treatment"])
#   
#   
#   #get risk difference
#   # Get predicted rates under treatment and control
#   newdata0 <- newdata1 <- data
#   newdata0$treatment <- 0
#   newdata1$treatment <- 1
#   
#   # Predict rates and transform to natural scale
#   pred0 <- predict(fit, newdata = newdata0, type = "response")
#   pred1 <- predict(fit, newdata = newdata1, type = "response")
#   
#   # Calculate rate difference (per 1000 person-years)
#   rd <- mean((pred1 - pred0)/exp(data$log_offset)) * 1000
#   
#   # Get SE using delta method
#   X0 <- model.matrix(~ treatment + age_group + sex + stunted +
#                        log(distance) + wealth_quintile + season, 
#                      data = newdata0)
#   X1 <- model.matrix(~ treatment + age_group + sex + stunted +
#                        log(distance) + wealth_quintile + season, 
#                      data = newdata1)
#   
#   # Get clustered variance matrix
#   vcov_cluster <- sandwich::vcovCL(fit, cluster = data$cluster_id)
#   
#   # Calculate SE using delta method
#   delta <- colMeans((X1 - X0)/exp(data$log_offset)) * 1000
#   se_rd <- sqrt(t(delta) %*% vcov_cluster %*% delta)
#   
#   return(c(coef = coef, se = se, rd=rd, se_rd=se_rd))
# }



analyze_trial_data <- function(data, 
                              subgroup = NULL,     # "age_group", "sex", or NULL
                              subgroup_level = NULL, # specific level for subgroup
                              adjusted = TRUE,
                              run_tmle=FALSE) {
  
  # Filter for subgroup if specified
  if(!is.null(subgroup) && !is.null(subgroup_level)) {
    data <- data[data[[subgroup]] == subgroup_level,]
  }
  
  # Define adjustment variables based on subgroup
  if(adjusted){
    
    data$log_distance <- log(data$distance)
    
    if(is.null(subgroup)){
      adjust_vars <- c("age_group", "sex", "stunted", "log_distance", 
                       "wealth_quintile", "season")
    }else{
      
      if(subgroup == "age_group"){
        adjust_vars <- c("sex", "stunted", "log_distance", 
                         "wealth_quintile", "season")} 
      if(subgroup == "sex"){
        adjust_vars <- c("age_group", "stunted", "log_distance", 
                         "wealth_quintile", "season")} 
    }
  }else{
    adjust_vars <- NULL
  }
  
  # 1. Cox Model for Hazard Ratio
  hr_formula <- if(adjusted) {
    as.formula(paste("Surv(time, status) ~ treatment +", 
                     paste(adjust_vars, collapse = " + "),
                     "+ cluster(cluster_id)"))
  }else{
    Surv(time, status) ~ treatment# + cluster(cluster_id)
  }
  
  cox_fit <- survival::coxph(hr_formula, data = data)
  
  # Get robust standard errors
  try(vcov_robust_hr <- vcovCL(cox_fit, cluster = data$cluster_id, type = "HC1"))

  # Extract results
  coef_hr <- coef(cox_fit)["treatment"]
  se_hr =NA
  try(se_hr <- sqrt(vcov_robust_hr["treatment", "treatment"]))
  hr_pval=NA
  try(hr_pval <- summary(cox_fit)$coefficients["treatment", 6])
  
  # # 4. Incidence Rate Ratio using Poisson GEE with offset
  # irr_formula <- if(adjusted) {
  #   as.formula(paste("dead ~ treatment +", 
  #                    paste(adjust_vars, collapse = " + "),
  #                    "+ offset(log_offset)"))
  # } else {
  #   dead ~ treatment + offset(log_offset)
  # }
  
  #5 CIR
  cir_formula <- if(adjusted) {
    as.formula(paste("dead ~ treatment +", 
                     paste(adjust_vars, collapse = " + ")))
  } else {
    dead ~ treatment 
  }
  
  
  fit_cir <- glm(as.formula(cir_formula),
                 family = poisson(link = "log"),
                 data = data)
  
  # Get robust standard errors
  vcov_robust_cir <- vcovCL(fit_cir, cluster = data$cluster_id, type = "HC1")
  
  # Extract results
  coef_cir <- coef(fit_cir)["treatment"]
  se_cir <- sqrt(vcov_robust_cir["treatment", "treatment"])
  
  
  #6 CID
  fit_cid <- glm(as.formula(cir_formula), family = "gaussian", data = data)
  
  # Get robust standard errors
  vcov_robust_cid <- vcovCL(fit_cid, cluster = data$cluster_id, type = "HC1")
  
  # Extract results
  coef_cid <- coef(fit_cid)["treatment"]
  se_cid <- sqrt(vcov_robust_cid["treatment", "treatment"])
  
  res = data.frame(
    hr=coef_hr,
    hr_se=se_hr,
    hr_ci_lb=coef_hr - 1.96*se_hr,
    hr_ci_ub=coef_hr + 1.96*se_hr,
    hr_pval=hr_pval,
    cir=coef_cir,
    cir_se=se_cir,
    cir_ci_lb=coef_cir - 1.96*se_cir,
    cir_ci_ub=coef_cir + 1.96*se_cir,
    cir_pval=summary(fit_cir)$coefficients["treatment", 4],
    cid=coef_cid,
    cid_se=se_cid,
    cid_ci_lb=coef_cid - 1.96*se_cid,
    cid_ci_ub=coef_cid + 1.96*se_cid,
    cid_pval=summary(fit_cid)$coefficients["treatment", 4]
  )
  
  #tmle estimates
  if(run_tmle){
    
    if(is.null(adjust_vars)){
      Wdf=NULL
    }else{
      Wdf=data.frame(W1=rep(1, nrow(data)), W2=rep(1, nrow(data)))
    }
    
    #SL.lib = c("SL.glm","SL.glmnet","SL.ranger")
    #SL.lib = c("SL.mean","SL.glm","SL.biglasso","SL.ranger")
    SL.lib = c("SL.mean","SL.glm","SL.step.interaction","SL.ranger")
    
    tmle_fit <- NULL
    
    try(tmle_fit <- tmle::tmle(
      A = data$treatment,
      Y = data$dead,
      W = Wdf,
      id= data$cluster_id,
      family = "binomial",
      g.SL.library = SL.lib,
      Q.SL.library = SL.lib,
      g.Delta.SL.library = SL.lib,  
      V.Q = 5, V.g = 5, V.Delta = 5, V.Z = 5,
      alpha = 0.05,
      verbose = FALSE
    ))
    if(!is.null(tmle_fit)){
      tmle_results <- summary(tmle_fit)$estimates
      tmle_res = data.frame(tmle_rr=tmle_results$RR$psi, 
                            tmle_log_rr=tmle_results$RR$log.psi, 
                            tmle_rr_log_se=sqrt(tmle_results$RR$var.log.psi),
                            tmle_rr_ci_lb=tmle_results$RR$CI[1],
                            tmle_rr_ci_ub=tmle_results$RR$CI[2],
                            tmle_rr_pval=tmle_results$RR$pvalue,
                            tmle_ate=tmle_results$ATE$psi,
                            tmle_ate_se=sqrt(tmle_results$ATE$var.psi),
                            tmle_ate_ci_lb=tmle_results$ATE$CI[1],
                            tmle_ate_ci_ub=tmle_results$ATE$CI[2],
                            tmle_ate_pval=tmle_results$ATE$pvalue
      )
    }else{
      tmle_results <- summary(tmle_fit)$estimates
      tmle_res = data.frame(tmle_rr=NA, 
                            tmle_log_rr=NA, 
                            tmle_rr_log_se=NA,
                            tmle_rr_ci_lb=NA,
                            tmle_rr_ci_ub=NA,
                            tmle_rr_pval=NA,
                            tmle_ate=NA,
                            tmle_ate_se=NA,
                            tmle_ate_ci_lb=NA,
                            tmle_ate_ci_ub=NA,
                            tmle_ate_pval=NA)
    }

    res = bind_cols(res, tmle_res)
  }
  
  
  # Add crude rates for context
  crude_rates <- with(data, {
    # Treatment group
    rate1 <- sum(status[treatment == 1]) / sum(person_time[treatment == 1])
    rate0 <- sum(status[treatment == 0]) / sum(person_time[treatment == 0])
    
    data.frame(
      treated_rate = rate1 * 1000,    # per 1000 person-years
      control_rate = rate0 * 1000,
      crude_rr = rate1/rate0
    )
  })
  
  
  res = bind_cols(crude_rates, res)
  
  return(res)
}
