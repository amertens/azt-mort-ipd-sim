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
source("study_params.R")


# simulation parameters
full_res = NULL
n_analyses = 7  # Main + 4 age groups + 2 sex groups
iter_df =NULL
iter_range=c(7,1000)
meta_method="REML"

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
      
      # Run meta-analyses
      results_list  <- meta_analyze_subgroups(sim_study_data, method=meta_method)
      meta_results <- results_list[[1]]
      study_results <- results_list[[2]]
      
      pooled_res <- NULL
      for(j in seq_along(meta_results)) {
        if(!is.null(meta_results[[j]])) {
          res <- data.frame(
            effect = exp(meta_results[[j]]$b[1]),
            se = sqrt(meta_results[[j]]$vb[1,1]),
            tau2 = meta_results[[j]]$tau2,
            pvalue = meta_results[[j]]$pval
          )
          res$group <- names(meta_results)[j]
        }
        pooled_res <- bind_rows(pooled_res, res)
      }
      pooled_res$study="pooled"
      
      res=bind_rows(pooled_res, study_results)
      res$iteration = sim
      
      full_res=bind_rows(full_res, res)
      
      # Save and print progress
      saveRDS(full_res, file=here("results/sim_results_interim.rds"))
      if(sim %% 10 == 0) cat(sprintf("Completed simulation %d of %d\n", sim, iter_range[2]))
}
    


saveRDS(full_res, file=here("results/sim_results.rds"))



