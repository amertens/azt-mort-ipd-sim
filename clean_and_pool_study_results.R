
# Load required packages
rm(list=ls())
library(tidyverse)
library(survival)
library(lme4)
library(here)
library(metafor)
library(foreach)
library(doParallel)
source("sim_data_functions.R")
source("study_params.R")

study_results <- readRDS(file=here("results/sim_results_interim_par.rds")) %>% mutate(adjusted=TRUE)
study_results_unadj <- readRDS(file=here("results/sim_results_interim_par_unadj.rds")) %>% mutate(adjusted=FALSE)
study_results <- bind_rows(study_results, study_results_unadj)
study_results <- clean_sim_results(study_results)

saveRDS(study_results, file=here("results/sim_results_clean.rds"))

#-------------------------------------------------------------------------------
# pooled results wrapper function around rma() that allows specification of method and measure
#-------------------------------------------------------------------------------


pool_results <- function(data, metric, method="REML"){
  
  fit<-NULL
  if(metric %in% c("cid", "ird","tmle_ate")){
    try(fit<-rma(yi=effect , sei=se, data=data, method=method, measure="GEN"))
  }else if(metric %in% c("cir", "irr","tmle_log_rr","hr" )){
    try(fit<-rma(yi=effect , sei=se, data=data, method=method, measure="RR"))
    if(method=="REML"){
      if(is.null(fit)){try(fit<-rma(yi=effect, sei=se, data=data, method="ML", measure="RR"))}
      if(is.null(fit)){try(fit<-rma(yi=effect, sei=se, data=data, method="DL", measure="RR"))}
      if(is.null(fit)){try(fit<-rma(yi=effect, sei=se, data=data, method="HE", measure="RR"))}
    }
  }
  
  out <- data %>%
    ungroup() %>%
    summarise(
      nstudies = length(study)#,
      #nmeas = sum(data[[ni]])
    ) %>%
    mutate(
      est = fit$beta[,1],
      se = fit$se,
      lb = fit$ci.lb,
      ub = fit$ci.ub,
      method.used=method,
      Qstat=fit$QE,
      tau2=fit$tau2,
      I2=fit$I2
    ) %>% as.data.frame()
  
  return(out)
  
}


#pool by iteration, metric, and subgroup
pooled_results_RE <- study_results %>%
  group_by(iteration, metric, level, adjusted) %>%
  do(pool_results(data=., metric=.$metric[1], method="REML")) 


pooled_results_FE <- study_results %>%
  group_by(iteration, metric, level, adjusted) %>%
  do(pool_results(data=., metric=.$metric[1], method="FE")) 

pooled_results_RE <- pooled_results_RE %>% rename(
  effect=est,
  lower=lb,
  upper=ub)
pooled_results_FE <- pooled_results_FE %>% rename(
  effect=est,
  lower=lb,
  upper=ub)

saveRDS(pooled_results_RE, file=here("results/pooled_results_RE.rds"))
saveRDS(pooled_results_FE, file=here("results/pooled_results_FE.rds"))


