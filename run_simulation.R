# Load required packages
rm(list=ls())
library(tidyverse)
library(survival)
library(lme4)
library(here)
library(metafor)
library(foreach)
library(doParallel)
source("simulation_functions.R")
source("sim_data_functions.R")
source("study_params.R")


# simulation parameters
full_res = NULL
iter_df =NULL
iter_range=c(1,1000)

sim=1

    
# Run simulations
for(sim in iter_range[1]:iter_range[2]){

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
      prop.table(table(sim_study_data$dead))*100
      
      # Test for interactions
      #interaction_results[[sim]] <- lapply(split(sim_study_data, sim_data$study), test_interactions)
      # interaction_results <- lapply(split(sim_study_data, sim_study_data$study), test_interactions)
      # head(interaction_results)

      # Run analyses
      results <- run_analysis(studies_data=sim_study_data, adjusted=TRUE, run_tmle=TRUE)
      
      results$iteration = sim
      
      full_res=bind_rows(full_res, results)
      
      # Save and print progress
      saveRDS(full_res, file=here("results/sim_results_interim.rds"))
      if(sim %% 10 == 0) cat(sprintf("Completed simulation %d of %d\n", sim, iter_range[2]))
}
    


saveRDS(full_res, file=here("results/sim_results.rds"))



