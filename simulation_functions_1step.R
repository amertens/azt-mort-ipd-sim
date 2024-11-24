

# Load required packages
rm(list=ls())
library(tidyverse)
library(survival)
library(lme4)
library(here)
library(metafor)
library(foreach)
library(doParallel)
library(sandwich)
source("simulation_functions.R")
source("sim_data_functions.R")
source("study_params.R")



# simulation parameters
full_res = NULL
iter_df =NULL
iter_range=c(1,1000)

sim=1
n_cores=20



full_res = NULL
sim =NULL
adjusted=FALSE
run_tmle=FALSE

  set.seed(sim)
  # Generate data for all studies
  sim_study_data <- do.call(rbind, lapply(names(study_params), function(study_name) {
    params <- study_params[[study_name]]
    simulate_trial_data(
      n_clusters = params$n_clusters,
      cluster_size = params$cluster_size,
      baseline_rate = params$baseline_rate,
      effect_size = params$effect_size,
      cv = params$cv,
      study_name = study_name
    )
  }))

  sim_study_data$study <- factor(sim_study_data$study)

# Run analyses
studies_data=sim_study_data
adjusted=TRUE
run_tmle=FALSE
one_step=TRUE


results <- run_analysis_one_step(studies_data=sim_study_data, adjusted=TRUE, run_tmle=FALSE){

 #  full_res=NULL
 #  data=studies_data
 #  
 # adjusted=TRUE
 # run_tmle=TRUE
 #  
 # subgroup = NULL # "age_group", "sex", or NULL
 # subgroup_level = NULL # specific level for subgroup
 # adjusted = TRUE
 # run_tmle=FALSE
  
  # Filter for subgroup if specified
  if(!is.null(subgroup) && !is.null(subgroup_level)) {
    data <- data[data[[subgroup]] == subgroup_level,]
  }
  
  
  # Define adjustment variables based on subgroup
  if(adjusted){
    
    data$log_distance <- log(data$distance)
    
    if(is.null(subgroup)){
      adjust_vars <- c("study", "age_group", "sex", "stunted", "log_distance", 
                       "wealth_quintile", "season")
    }else{
      
      if(subgroup == "age_group"){
        adjust_vars <- c("study", "sex", "stunted", "log_distance", 
                         "wealth_quintile", "season")} 
      if(subgroup == "sex"){
        adjust_vars <- c("study", "age_group", "stunted", "log_distance", 
                         "wealth_quintile", "season")} 
    }
  }else{
    adjust_vars <- NULL
  }
  

  
  #5 CIR
  cir_formula <- if(adjusted){
    as.formula(paste("dead ~ treatment +", 
                     paste(adjust_vars, collapse = " + "), " + (cluster_id|study)"))
  }else{
    dead ~ treatment + (cluster_id|study)
  }
  
  
  #need to check how to get the CID and CIR
  fit_cir <- glmer(as.formula(cir_formula), family = binomial(link = "log"), data = data)
  res= summary(fit_cir)$coefficients

  # Extract results
  coef_cir <- res[2,1]
  se_cir <- res[2,2]
  
  #6 CID
  fit_cid <- lmer(as.formula(cir_formula),data = data)
  res_cid = summary(fit_cid)$coefficients
  
  # Extract results
  coef_cid <- res_cid[2,1]
  se_cid <- res_cid[2,2]
  
  res = data.frame(
    # hr=coef_hr,
    # hr_se=se_hr,
    # hr_ci_lb=coef_hr - 1.96*se_hr,
    # hr_ci_ub=coef_hr + 1.96*se_hr,
    # hr_pval=hr_pval,
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