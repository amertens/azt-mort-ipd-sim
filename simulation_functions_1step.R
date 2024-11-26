

# Load required packages
rm(list=ls())
library(tidyverse)
library(survival)
library(lme4)
library(lmerTest)
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


run_analysis_one_step<-function(studies_data=sim_study_data, adjusted=TRUE, run_tmle=FALSE){
  
  require(lme4)
  require(lmerTest)

  full_res=NULL
  data=studies_data

 adjusted=TRUE
 run_tmle=TRUE

 subgroup = NULL # "age_group", "sex", or NULL
 subgroup_level = NULL # specific level for subgroup
 adjusted = TRUE
 run_tmle=FALSE
  
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
  pval_cir <- res[2,4]
  
  #6 CID
  fit_cid <- lmer(as.formula(cir_formula),data = data)
  res_cid = summary(fit_cid)$coefficients
  
  # Extract results
  coef_cid <- res_cid[2,1]
  se_cid <- res_cid[2,2]
  pval_cid <- 2*pt(-abs(res_cid[2,3]), df=nrow(data)-1) 
  
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
    cir_pval=pval_cir,
    cid=coef_cid,
    cid_se=se_cid,
    cid_ci_lb=coef_cid - 1.96*se_cid,
    cid_ci_ub=coef_cid + 1.96*se_cid,
    cid_pval=pval_cid
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









# Modified TMLE with 1-step random effects
tmle_random_effects <- function(data, 
                                g.SL.library, 
                                Q.SL.library,
                                family = "binomial",
                                weighting = "variance") {  # Options: "size", "variance", "none"
  
  # Calculate study-specific weights
  get_study_weights <- function(data, type = weighting) {
    studies <- unique(data$study)
    
    weights <- switch(type,
                      "size" = {
                        # Weights proportional to study size
                        ns <- table(data$study)
                        setNames(ns/sum(ns), names(ns))
                      },
                      "variance" = {
                        # Inverse variance weights based on initial outcome model
                        study_vars <- sapply(studies, function(s) {
                          study_data <- data[data$study == s,]
                          # Fit initial model to get variance estimate
                          fit <- glm(dead ~ treatment, family = binomial, data = study_data)
                          var(predict(fit, type = "response"))
                        })
                        w <- 1/study_vars
                        setNames(w/sum(w), studies)
                      },
                      "none" = {
                        # Equal weights
                        setNames(rep(1/length(studies), length(studies)), studies)
                      }
    )
    
    # Map weights to individual observations
    weights[as.character(data$study)]
  }
  
  
  # Add mixed effects model to libraries
  # Q.SL.library <- c(Q.SL.library, "SL.glmer")
  # g.SL.library <- c(g.SL.library, "SL.glmer")
  
  # Custom targeting step that incorporates random effects
  targeting_step <- function(Q, g, H, id) {
    # Create epsilon model with random effects
    eps_model <- glmer(dead ~ offset(qlogis(Q)) + H + (1|study),
                       family = binomial(),
                       data = data)
    
    # Update Q incorporating random effects
    Q_star <- plogis(qlogis(Q) + ranef(eps_model)$study[id,1] + 
                       fixef(eps_model)[2] * H)
    
    return(Q_star)
  }
  
  
  # Get study weights
  study_weights <- get_study_weights(data, weighting)
  
  # Initial SuperLearner fits
  # Need to pass study ID to SuperLearner
  data$log_distance <- log(data$distance)
  Q <- SuperLearner(Y = data$dead,
                    X = data[, c("age_months", "sex", "season", "wealth_quintile", "log_distance", "treatment")],
                    family = family,
                    cvControl = SuperLearner.CV.control(V = 4, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                    SL.library = Q.SL.library,
                    id = as.numeric(data$study))
  
  # Treatment mechanism
  g <- SuperLearner(Y = data$treatment,
                    X = data[, c("age_months", "sex", "season", "wealth_quintile", "log_distance")],
                    family = binomial(),
                    cvControl = SuperLearner.CV.control(V = 4, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                    SL.library = g.SL.library,
                    id = as.numeric(data$study))
  
  # Get predictions
  Q_pred <- predict(Q)$pred
  g_pred <- predict(g)$pred
  
  # Get clever covariate
  # And incorporate study-specific weights
  H1W <- study_weights * data$treatment/g_pred
  H0W <- -study_weights * (1-data$treatment)/(1-g_pred)
  H <- H1W + H0W
  
  
  # Update step incorporating random effects
  Q_star <- targeting_step(Q_pred, g_pred, H, data$study)
  
  
  # Get treatment-specific means
  mu1 <- mean(Q_star[data$treatment == 1])
  mu0 <- mean(Q_star[data$treatment == 0])
  
  # Risk Difference
  rd <- mu1 - mu0
  
  # Risk Ratio
  rr <- mu1/mu0
  log_rr <- log(rr)  # Log RR for inference
  
  # Influence curves
  # For Risk Difference
  IC_rd <- H1W * (data$dead - Q_star) - H0W * (data$dead - Q_star) + 
    Q_star[data$treatment == 1] - Q_star[data$treatment == 0] - rd
  
  # For Log Risk Ratio (using delta method)
  IC_log_rr <- (1/mu1) * (H1W * (data$dead - Q_star)) - 
    (1/mu0) * (H0W * (data$dead - Q_star)) + 
    (Q_star[data$treatment == 1]/mu1 - Q_star[data$treatment == 0]/mu0 - log_rr)
  
  # Variance components from ICs
  study <- data$study
  var_components_rd <- lmer(IC_rd ~ 1 + (1|study))
  var_components_log_rr <- lmer(IC_log_rr ~ 1 + (1|study))
  
  # Total variance incorporating between-study heterogeneity
  variance_rd <- VarCorr(var_components_rd)$study[1] + 
    attr(VarCorr(var_components_rd), "sc")^2/nrow(data)
  
  variance_log_rr <- VarCorr(var_components_log_rr)$study[1] + 
    attr(VarCorr(var_components_log_rr), "sc")^2/nrow(data)
  
  
  # Test statistics
  z_rd <- rd/sqrt(variance_rd)
  z_log_rr <- log_rr/sqrt(variance_log_rr)
  
  # P-values (two-sided)
  p_rd <- 2 * (1 - pnorm(abs(z_rd)))
  p_rr <- 2 * (1 - pnorm(abs(z_log_rr)))
  
  # Confidence intervals
  ci_rd <- rd + c(-1.96, 1.96) * sqrt(variance_rd)
  ci_rr <- exp(log_rr + c(-1.96, 1.96) * sqrt(variance_log_rr))
  
  return(data.frame(
    tmle_rr=rr, 
    tmle_log_rr= log(rr), 
    tmle_rr_log_se=variance_log_rr,
    tmle_rr_ci_lb=ci_rr[1],
    tmle_rr_ci_ub=ci_rr[2],
    tmle_rr_pval=p_rr,
    tmle_ate=rd,
    tmle_ate_se=sqrt(variance_rd),
    tmle_ate_ci_lb=ci_rd[1],
    tmle_ate_ci_ub=ci_rd[2],
    tmle_ate_pval=p_rd))
  
}

# Simulation to test
set.seed(123)

# Define SuperLearner libraries
Q.SL.library <- c("SL.glm", "SL.mean")
g.SL.library <- c("SL.glm", "SL.mean")

# Fit models
results_1step <- tmle_random_effects(
  data = sim_study_data,
  Q.SL.library = Q.SL.library,
  g.SL.library = g.SL.library
)