
#calculate truth

# Load required packages
rm(list=ls())
library(tidyverse)
library(here)
library(metafor)
source("sim_data_functions.R")
source("study_params.R")

#---------------------------------------------------------------------------
# make large sample size
#---------------------------------------------------------------------------

study_params$AVENIR$n_clusters <- 1000
study_params$AVENIR$cluster_size <- 1000
study_params$MORDOR_Niger$n_clusters <- 1000
study_params$MORDOR_Niger$cluster_size <- 1000
study_params$MORDOR_Malawi$n_clusters <- 1000
study_params$MORDOR_Malawi$cluster_size <- 1000
study_params$MORDOR_Tanzania$n_clusters <- 1000
study_params$MORDOR_Tanzania$cluster_size <- 1000
study_params$MORDOR_II_Niger$n_clusters <- 1000
study_params$MORDOR_II_Niger$cluster_size <- 1000
study_params$TANA_I$n_clusters <- 1000
study_params$TANA_I$cluster_size <- 1000


#---------------------------------------------------------------------------
# Generate counterfactual data
#---------------------------------------------------------------------------

  set.seed(123456)
  # Generate data for all studies
  sim_study_data_A1 <- do.call(rbind, lapply(names(study_params), function(study_name) {
    params <- study_params[[study_name]]
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
  sim_study_data_A0 <- do.call(rbind, lapply(names(study_params), function(study_name) {
    params <- study_params[[study_name]]
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
  
  table(sim_study_data_A1$deaths)
  table(sim_study_data_A0$deaths)
  
  #---------------------------------------------------------------------------
  # Calculate truth by study and subgroup
  #---------------------------------------------------------------------------

calc_group_truth <- function(A1, A0, group_var) {
    
    truth_A1 <- A1 %>% group_by(study, !!sym(group_var)) %>%
      summarize(
        incidence = mean(deaths),
        person_time = sum(exp(log_offset)),
        IR=incidence/person_time * 1000000) %>%
       rename(group = !!sym(group_var))
    truth_A0 <- A0 %>% group_by(study, !!sym(group_var)) %>%
      summarize(
        incidence = mean(deaths),
        person_time = sum(exp(log_offset)),
        IR=incidence/person_time * 1000000)%>%
      rename(group = !!sym(group_var))
    
    truth <- left_join(truth_A1, truth_A0,  by=c("study","group")) %>%
      mutate(
        IR_diff = IR.x - IR.y,
        IR_ratio = IR.x / IR.y
      ) %>%
      select(study, IR_diff, IR_ratio, group) 
    
    return(truth)
}
  
  truth_main = calc_group_truth(sim_study_data_A1 %>% mutate(group = "main"), 
                                sim_study_data_A0 %>% mutate(group = "main"), "group") 
  truth_age = calc_group_truth(sim_study_data_A1, sim_study_data_A0, "age_group") %>% mutate(group=paste0("age_", group), group=factor(group, levels=c("age_1-5mo", "age_6-11mo", "age_12-23mo", "age_24-59mo"))) %>% arrange(group) 
  truth_sex = calc_group_truth(sim_study_data_A1, sim_study_data_A0, "sex") %>% mutate(group=case_when(group==1 ~ "sex_male",group==0 ~ "sex_female"))

  
  #combine
  truth <- bind_rows(truth_main, truth_age, truth_sex)
  saveRDS(truth, file=here("results/truth.rds"))
  
  
  truth_main %>% filter(study=="MORDOR_Niger") 
  truth_age %>% filter(study=="MORDOR_Niger") 
  