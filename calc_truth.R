
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

names(study_params)
study_params$AVENIR$n_clusters <- 1000
study_params$AVENIR$cluster_size <- 1000
study_params$MORDOR_Niger$n_clusters <- 1000
study_params$MORDOR_Niger$cluster_size <- 1000
study_params$MORDOR_Malawi$n_clusters <- 1000
study_params$MORDOR_Malawi$cluster_size <- 1000
study_params$MORDOR_Tanzania$n_clusters <- 1000
study_params$MORDOR_Tanzania$cluster_size <- 1000
study_params$CHAT$n_clusters <- 1000
study_params$CHAT$cluster_size <- 1000
study_params$TANA_I$n_clusters <- 1000
study_params$TANA_I$cluster_size <- 1000


#---------------------------------------------------------------------------
# Generate counterfactual data
#---------------------------------------------------------------------------

study_name=""
  
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
  
  # table(sim_study_data_A1$dead)
  # table(sim_study_data_A0$dead)
  
#---------------------------------------------------------------------------
# Calculate truth by study and subgroup
#---------------------------------------------------------------------------
  
  A1 = sim_study_data_A1 %>% mutate(level = "main")
  A0 = sim_study_data_A0 %>% mutate(level = "main")
  level_var = "level"

calc_group_truth <- function(A1, A0, level_var) {
  
  head(A1)
    
    truth_A1 <- A1 %>% group_by(study, !!sym(level_var)) %>%
      summarize(
        incidence = mean(dead),
        person_time = sum(exp(log_offset)),
        IR=incidence/person_time * 1000000) %>%
       rename(level = !!sym(level_var))
    
    head(truth_A1)
    
    truth_A0 <- A0 %>% group_by(study, !!sym(level_var)) %>%
      summarize(
        incidence = mean(dead),
        person_time = sum(exp(log_offset)),
        IR=incidence/person_time * 1000000)%>%
      rename(level = !!sym(level_var))
    
    truth <- left_join(truth_A1, truth_A0,  by=c("study","level")) %>%
      mutate(
        cid= incidence.x - incidence.y,
        cir= incidence.x / incidence.y,
        ird = IR.x - IR.y,
        irr = IR.x / IR.y
      ) %>%
      select(study, cid,  cir, ird, irr, level)
    


    # # HR calculation: Ratio of event rates
    # delta1 <- sum(A1$dead) / sum(A1$time)
    # delta0 <- sum(A0$dead) / sum(A0$time)
    # true_hr<- delta1/delta0
    # 
    return(truth)
}
  
  truth_main = calc_group_truth(sim_study_data_A1 %>% mutate(level = "main"), 
                                sim_study_data_A0 %>% mutate(level = "main"), "level") 
  truth_age = calc_group_truth(sim_study_data_A1, sim_study_data_A0, "age_group") %>% mutate(level=factor(level, levels=c("1-5mo", "6-11mo", "12-23mo", "24-59mo"))) %>% arrange(level) 
  truth_sex = calc_group_truth(sim_study_data_A1, sim_study_data_A0, "sex") %>% mutate(level=case_when(level==1 ~ "male",level==0 ~ "female"))

  
  #combine
  truth <- bind_rows(truth_main, truth_age, truth_sex)
  
  #transform for merge with results
  truth <- truth %>% pivot_longer(cols=c(cid, cir, ird, irr), names_to="metric", values_to="true_value")
  table(truth$metric)
  
  saveRDS(truth, file=here("results/truth.rds"))
  
  
  truth_main %>% filter(study=="MORDOR_Niger") 
  truth_age %>% filter(study=="MORDOR_Niger") 
  
  