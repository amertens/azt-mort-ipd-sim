
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

#transform data to long format with columns for effect and se, and the group
study_result_est <- study_results %>% 
  mutate(crude_rr = log(crude_rr)) %>% #make sure all ratios are on the log scale
  select(iteration, study, group, level, adjusted, crude_rr, hr, cid, cir, tmle_log_rr, tmle_ate) %>%
  gather(key="metric", value="effect", -iteration, -study, -group, -level, -adjusted) %>% distinct()

study_result_se <- study_results %>% 
  select(iteration, study, group, level, adjusted, hr_se, cid_se, cir_se, tmle_rr_log_se, tmle_ate_se) %>%
  gather(key="metric", value="se", -iteration, -study, -group, -level, -adjusted) %>%
  mutate(metric = gsub("_se", "", metric)) %>% 
  mutate(metric = gsub("tmle_rr_log", "tmle_log_rr", metric)) %>% 
  distinct()

study_result_pval <- study_results %>% 
  select(iteration, study, group, level, adjusted, hr_pval,  cir_pval, cid_pval, tmle_rr_pval, tmle_ate_pval) %>%
  gather(key="metric", value="pval", -iteration, -study, -group, -level, -adjusted) %>%
  mutate(metric = gsub("_pval", "", metric)) %>% 
  mutate(metric = gsub("tmle_rr", "tmle_log_rr", metric)) %>% 
  distinct()

study_results <- left_join(study_result_est, study_result_se, by=c("iteration", "study", "group", "level", "metric", "adjusted"))
study_results <- left_join(study_results, study_result_pval, by=c("iteration", "study", "group", "level", "metric", "adjusted"))
study_results <- study_results %>% mutate(lower=effect-1.96*se, upper=effect+1.96*se)

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
pooled_results_REML <- study_results %>% filter(level=="main") %>%
  group_by(iteration, metric, adjusted) %>%
  do(pool_results(data=., metric=.$metric[1], method="REML")) %>%
  mutate(pool_method="REML")

pooled_results_ML <- study_results %>% filter(level=="main") %>%
  group_by(iteration, metric, adjusted) %>%
  do(pool_results(data=., metric=.$metric[1], method="ML")) %>%
  mutate(pool_method="ML") 

pooled_results_DL <- study_results %>% filter(level=="main") %>%
  group_by(iteration, metric, adjusted) %>%
  do(pool_results(data=., metric=.$metric[1], method="DL")) %>%
  mutate(pool_method="DL") 

pooled_results_HE <- study_results %>% filter(level=="main") %>%
  group_by(iteration, metric, adjusted) %>%
  do(pool_results(data=., metric=.$metric[1], method="HE")) %>%
  mutate(pool_method="HE") 


pooled_results <- bind_rows(pooled_results_REML, pooled_results_ML, 
                            pooled_results_DL, pooled_results_HE)
saveRDS(pooled_results, file=here("results/pooling_comp_results.rds"))


