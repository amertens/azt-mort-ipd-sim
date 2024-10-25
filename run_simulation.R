# Load required packages
library(tidyverse)
library(survival)
library(lme4)
library(here)
library(metafor)
library(foreach)
library(doParallel)
source("simulation_functions.R")
source("study_params.R")

n_sims = 100
parallel=FALSE

# Run simulation incorporating real trial parameters
run_power_simulation <- function(n_sims = 1000, 
                                 parallel=FALSE,     
                                 n_analyses = 7  # Main + 4 age groups + 2 sex groups
                                 ) {
  
  if(parallel){
    
  # Setup parallel processing
  n_cores <- parallel::detectCores()/2
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  results <- foreach(
    sim = 1:n_sims,
    .combine = rbind,
    .packages = c( "lme4", "metafor")
  ) %dopar% {
    
    # Simulate each study
    study_data <- do.call(rbind, lapply(names(study_params), function(study_name) {
      params <- study_params[[study_name]]
      simulate_trial_data(
        n_clusters = params$n_clusters,
        cluster_size = params$cluster_size,
        baseline_rate = params$baseline_rate,
        effect_size = params$effect_size,
        study_name = study_name
      )
    }))
    
    # Analyze meta-analysis
    analysis <- analyze_meta_analysis(study_data)
    
    # Return key results
    c(
      pooled_effect = exp(analysis$meta_results$b[1]),
      se = sqrt(analysis$meta_results$vb[1,1]),
      tau2 = analysis$meta_results$tau2,
      pvalue = analysis$meta_results$pval,
      # Add study-specific results
      analysis$individual_results[,"coef.treatment"]
    )
  }
  
  stopCluster(cl)
  }else{

    # Initialize results matrices
    results <- array(NA, 
                     dim=c(n_sims, 4, n_analyses),
                     dimnames=list(
                       NULL,
                       c("effect", "se", "tau2", "pvalue"),
                       c("main", paste0("age_", c("1-5mo", "6-11mo", "12-23mo", "24-59mo")),
                         "sex_female", "sex_male")
                     ))
    # Run simulations
    sim=1
    # Run simulations
    for(sim in 1:n_sims){
      
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
      meta_results <- meta_analyze_subgroups(sim_study_data)
      
      # Store results
      for(i in seq_along(meta_results)) {
        if(!is.null(meta_results[[i]])) {
          results[sim, , names(meta_results)[i]] <- c(
            effect = exp(meta_results[[i]]$b[1]),
            se = sqrt(meta_results[[i]]$vb[1,1]),
            tau2 = meta_results[[i]]$tau2,
            pvalue = meta_results[[i]]$pval
          )
        }
      }
      
      # Print progress
      if(sim %% 10 == 0) cat(sprintf("Completed simulation %d of %d\n", sim, n_sims))
    }
    
  }
  
  # Calculate power for each analysis
  power_results <- apply(results, 3, function(x) {
    data.frame(
      power = mean(x[,"pvalue"] < 0.05, na.rm=TRUE),
      mean_effect = mean(x[,"effect"], na.rm=TRUE),
      mean_se = mean(x[,"se"], na.rm=TRUE),
      mean_tau2 = mean(x[,"tau2"], na.rm=TRUE)
    )
  })
  
  # Convert to nice data frame
  power_summary <- do.call(rbind, power_results)
  power_summary$analysis <- rownames(power_summary)
  print(power_summary)
  
  return(power_summary)
}

# Run simulation
set.seed(123)
sim_results <- run_power_simulation(n_sims = 1, parallel=FALSE)  # Increase for final analysis
print(sim_results)

saveRDS(sim_results, file=here("results/sim_results.rds"))

# Plot results



# Visualization of power by subgroup
library(ggplot2)

ggplot(sim_results, aes(x=analysis, y=power)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Statistical Power by Analysis",
       x="Analysis", 
       y="Power") +
  theme_minimal()

# Effect size comparison
ggplot(sim_results, aes(x=analysis, y=mean_effect)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_effect-1.96*mean_se, 
                    ymax=mean_effect+1.96*mean_se)) +
  coord_flip() +
  labs(title="Effect Sizes by Analysis",
       x="Analysis",
       y="Risk Ratio") +
  theme_minimal()


#TO DO:
#Need to compare doing a traditional meta-analysis with an IPD meta-analysis
#Need to add in subgroup analyses
#Need to come up with a list of secondary analyses (AMR)
#Need to come up with a list of sensitivity analyses
  #A