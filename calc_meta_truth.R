
#calculate truth

# Load required packages
rm(list=ls())
library(tidyverse)
library(here)
library(metafor)
source("sim_data_functions.R")
source("study_params.R")

#---------------------------------------------------------------------------
# Notes on truth calculation

#For the fixed effect meta-analysis, need to increase the sample size of all studies, 
#keeping the relative sizes the same. (how do I change clusters number versus cluster size?)

#For the random effects meta-analysis, need to increase the number of studies so they are all equally larger



#---------------------------------------------------------------------------
# make large sample size
#---------------------------------------------------------------------------

study_params_RE <- study_params_FE <- study_params

study_params_RE$AVENIR$n_clusters <- 1000
study_params_RE$AVENIR$cluster_size <- 1000
study_params_RE$MORDOR_Niger$n_clusters <- 1000
study_params_RE$MORDOR_Niger$cluster_size <- 1000
study_params_RE$MORDOR_Malawi$n_clusters <- 1000
study_params_RE$MORDOR_Malawi$cluster_size <- 1000
study_params_RE$MORDOR_Tanzania$n_clusters <- 1000
study_params_RE$MORDOR_Tanzania$cluster_size <- 1000
study_params_RE$MORDOR_II_Niger$n_clusters <- 1000
study_params_RE$MORDOR_II_Niger$cluster_size <- 1000
study_params_RE$TANA_I$n_clusters <- 1000
study_params_RE$TANA_I$cluster_size <- 1000

#increase size of fixed effect studies
mult_factor_n_clusters = 1
mult_factor_cluster_size = 10
study_params_FE$AVENIR$n_clusters <- study_params_FE$AVENIR$n_clusters *mult_factor_n_clusters
study_params_FE$AVENIR$cluster_size <- study_params_FE$AVENIR$cluster_size *mult_factor_cluster_size
study_params_FE$MORDOR_Niger$n_clusters <- study_params_FE$MORDOR_Niger$n_clusters *mult_factor_n_clusters
study_params_FE$MORDOR_Niger$cluster_size <- study_params_FE$MORDOR_Niger$cluster_size *mult_factor_cluster_size
study_params_FE$MORDOR_Malawi$n_clusters <- study_params_FE$MORDOR_Malawi$n_clusters *mult_factor_n_clusters
study_params_FE$MORDOR_Malawi$cluster_size <- study_params_FE$MORDOR_Malawi$cluster_size *mult_factor_cluster_size
study_params_FE$MORDOR_Tanzania$n_clusters <- study_params_FE$MORDOR_Tanzania$n_clusters *mult_factor_n_clusters
study_params_FE$MORDOR_Tanzania$cluster_size <- study_params_FE$MORDOR_Tanzania$cluster_size *mult_factor_cluster_size
study_params_FE$MORDOR_II_Niger$n_clusters <- study_params_FE$MORDOR_II_Niger$n_clusters *mult_factor_n_clusters
study_params_FE$MORDOR_II_Niger$cluster_size <- study_params_FE$MORDOR_II_Niger$cluster_size *mult_factor_cluster_size
study_params_FE$TANA_I$n_clusters <- study_params_FE$TANA_I$n_clusters *mult_factor_n_clusters
study_params_FE$TANA_I$cluster_size <- study_params_FE$TANA_I$cluster_size *mult_factor_cluster_size

#---------------------------------------------------------------------------
# Generate counterfactual data
#---------------------------------------------------------------------------

  set.seed(123456)
  # Generate data for all studies
  sim_study_data_A1_FE <- do.call(rbind, lapply(names(study_params_FE), function(study_name) {
    params <- study_params_FE[[study_name]]
    simulate_trial_data(
      n_clusters = params$n_clusters,
      cluster_size = params$cluster_size, 
      baseline_rate = params$baseline_rate,
      effect_size = params$effect_size,
      cv = params$cv,
      study_name = study_name,
      A=1
    )
  }))
  
  
  set.seed(123456)
  # Generate data for all studies
  sim_study_data_A0_FE <- do.call(rbind, lapply(names(study_params_FE), function(study_name) {
    params <- study_params_FE[[study_name]]
    simulate_trial_data(
      n_clusters = params$n_clusters,
      cluster_size = params$cluster_size, 
      baseline_rate = params$baseline_rate,
      effect_size = params$effect_size,
      cv = params$cv,
      study_name = study_name,
      A=0
    )
  }))
  
  
  set.seed(123456)
  # Generate data for all studies
  sim_study_data_A1_RE <- do.call(rbind, lapply(names(study_params_RE), function(study_name) {
    params <- study_params_RE[[study_name]]
    simulate_trial_data(
      n_clusters = params$n_clusters,
      cluster_size = params$cluster_size, 
      baseline_rate = params$baseline_rate,
      effect_size = params$effect_size,
      cv = params$cv,
      study_name = study_name,
      A=1
    )
  }))
  
  
  set.seed(123456)
  # Generate data for all studies
  sim_study_data_A0_RE <- do.call(rbind, lapply(names(study_params_RE), function(study_name) {
    params <- study_params_RE[[study_name]]
    simulate_trial_data(
      n_clusters = params$n_clusters,
      cluster_size = params$cluster_size, 
      baseline_rate = params$baseline_rate,
      effect_size = params$effect_size,
      cv = params$cv,
      study_name = study_name,
      A=0
    )
  }))
  
#---------------------------------------------------------------------------
# Calculate truth by study and subgroup
#---------------------------------------------------------------------------
  
  # A1 = sim_study_data_A1 %>% mutate(level = "main")
  # A0 = sim_study_data_A0 %>% mutate(level = "main")
  # level_var = "level"

calc_truth <- function(A1, A0, level_var) {
  
  head(A1)
    
    truth_A1 <- A1 %>% group_by(!!sym(level_var)) %>%
      summarize(
        incidence = mean(dead),
        person_time = sum(exp(log_offset)),
        IR=incidence/person_time * 1000000) %>%
       rename(level = !!sym(level_var))
    
    head(truth_A1)
    
    truth_A0 <- A0 %>% group_by(!!sym(level_var)) %>%
      summarize(
        incidence = mean(dead),
        person_time = sum(exp(log_offset)),
        IR=incidence/person_time * 1000000)%>%
      rename(level = !!sym(level_var))
    
    truth <- left_join(truth_A1, truth_A0,  by=c("level")) %>%
      mutate(
        cid= incidence.x - incidence.y,
        cir= incidence.x / incidence.y,
        ird = IR.x - IR.y,
        irr = IR.x / IR.y
      ) %>%
      select(cid,  cir, ird, irr, level)
    
    return(truth)
}
  
  truth_main_FE = calc_truth(sim_study_data_A1_FE %>% mutate(level = "main"), 
                                sim_study_data_A0_FE %>% mutate(level = "main"), "level") 
  truth_age_FE = calc_truth(sim_study_data_A1_FE, sim_study_data_A0_FE, "age_group") %>% mutate(level=factor(level, levels=c("1-5mo", "6-11mo", "12-23mo", "24-59mo"))) %>% arrange(level) 
  truth_sex_FE = calc_truth(sim_study_data_A1_FE, sim_study_data_A0_FE, "sex") %>% mutate(level=case_when(level==1 ~ "male",level==0 ~ "female"))

  
  #combine
  truth_FE <- bind_rows(truth_main_FE, truth_age_FE, truth_sex_FE)
  
  
  truth_main_RE = calc_truth(sim_study_data_A1_RE %>% mutate(level = "main"), 
                                   sim_study_data_A0_RE %>% mutate(level = "main"), "level") 
  truth_age_RE = calc_truth(sim_study_data_A1_RE, sim_study_data_A0_RE, "age_group") %>% mutate(level=factor(level, levels=c("1-5mo", "6-11mo", "12-23mo", "24-59mo"))) %>% arrange(level) 
  truth_sex_RE = calc_truth(sim_study_data_A1_RE, sim_study_data_A0_RE, "sex") %>% mutate(level=case_when(level==1 ~ "male",level==0 ~ "female"))
  
  
  #combine
  truth_RE <- bind_rows(truth_main_RE, truth_age_RE, truth_sex_RE)
  
  
  #transform for merge with results
  truth_FE <- truth_FE %>% pivot_longer(cols=c(cid, cir, ird, irr), names_to="metric", values_to="true_value")
  truth_RE <- truth_RE %>% pivot_longer(cols=c(cid, cir, ird, irr), names_to="metric", values_to="true_value")
  table(truth_RE$metric)
  
  saveRDS(truth_FE, file=here("results/meta_truth_FE.rds"))
  saveRDS(truth_RE, file=here("results/meta_truth_RE.rds"))
  
  

  summary(truth_FE$true_value)
  summary(truth_RE$true_value)
  