


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
library(tmle)
library(SuperLearner)
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

data=sim_study_data




# First create a custom SuperLearner wrapper for mixed effects
SL.glmer <- function(Y, X, newX, family, obsWeights, id, ...) {
  # Ensure id is in the data
  X$study <- id
  if(!is.null(newX)) newX$study <- id
  
  # Fit mixed model
  fit <- glmer(Y ~ . + (1|study), 
               data = X, 
               family = family, 
               weights = obsWeights)
  
  # Predictions
  pred <- predict(fit, newdata = if(!is.null(newX)) newX else X, 
                  type = "response", 
                  allow.new.levels = TRUE)
  
  fit <- list(object = fit)
  class(fit) <- "SL.glmer"
  out <- list(pred = pred, fit = fit)
  return(out)
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

results_1step

