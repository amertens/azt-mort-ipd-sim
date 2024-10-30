
#calculate truth

# Load required packages
rm(list=ls())
library(tidyverse)
library(here)
library(metafor)
source("sim_data_functions.R")
source("study_params.R")

names(study_params)

#make large sample size
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


# Run simulations
  gc()
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
  
  
  # Calculate truth by study
  truth_A1 <- sim_study_data_A1 %>%
    group_by(study) %>%
    summarize(
      incidence = mean(deaths),
      person_time = sum(exp(log_offset)),
      IR=incidence/person_time * 1000000
    )
  
  truth_A0 <- sim_study_data_A0 %>%
    group_by(study) %>%
    summarize(
      incidence = mean(deaths),
      person_time = sum(exp(log_offset)),
      IR=incidence/person_time * 1000000
    )
  
  truth <- left_join(truth_A1, truth_A0, by="study") %>%
    mutate(
      IR_diff = IR.x - IR.y,
      IR_ratio = IR.x / IR.y
    ) %>%
    select(study, IR_diff, IR_ratio) 

    truth%>%knitr::kable()
  
  saveRDS(truth, file=here("study"))
  
  #need to figure out how to include person-time
  
  
  