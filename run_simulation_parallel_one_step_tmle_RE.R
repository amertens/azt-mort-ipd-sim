
# Load required packages
rm(list=ls())
library(tidyverse)
library(survival)
library(lme4)
library(here)
library(metafor)
library(foreach)
library(doParallel)
library(sandwich)
library(SuperLearner)
source("simulation_functions_1step.R")
source("sim_data_functions.R")
source("study_params.R")



# simulation parameters
full_res = NULL
iter_df =NULL
iter_range=c(1,1000)

sim=1
n_cores=50

run_sim_par <- function(full_res = NULL,sim =NULL){
  
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
  
  # Define SuperLearner libraries
  #Q.SL.library <- g.SL.library <- c("SL.glm", "SL.mean")
  Q.SL.library <- g.SL.library <- c("SL.mean","SL.glm","SL.step.interaction","SL.ranger")

  # Run analyses
  results <- tmle_random_effects(
    data = sim_study_data,
    Q.SL.library = Q.SL.library,
    g.SL.library = g.SL.library)

  results$iteration = sim
  
  
  full_res=bind_rows(full_res, results)
  
  return(full_res)
}


cl <- makeCluster(n_cores)
clusterExport(cl, c("run_sim_par","study_params","simulate_trial_data","generate_covariates",
                    "tmle_random_effects","iter_range"), envir=environment())


clusterEvalQ(cl,lapply(c("tmle","SuperLearner","data.table","tidyverse","survival","lme4","sandwich","here","ranger"), FUN = function(X) {
  do.call("require", list(X))
}))

# temp = run_sim_par(full_res = NULL,
#             sim =1)
# temp

full_res=NULL
full_res <- readRDS(here("results/sim_results_par_1step_tmle_RE.rds"))
 max(full_res$iteration)/50

#for(i in iter_range[1]:(iter_range[2]/n_cores)){
for(i in 19:20){
  cat(i,"\n")
  res=parLapply(cl=cl, c(((i-1)*50+1):(50*i)), function(z) run_sim_par(full_res = NULL,
                                                        sim =z))
  resdf=data.table::rbindlist(res)
  full_res=bind_rows(full_res, resdf)
  saveRDS(full_res, file=here("results/sim_results_interim_par_1step_tmle_RE.rds"))
}

length(unique(full_res$iteration))

saveRDS(full_res, file=here("results/sim_results_par_1step_tmle_RE.rds"))
stopCluster(cl)


