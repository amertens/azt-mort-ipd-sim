
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

truth <- readRDS(here("results/truth.rds"))

sim_results_full <- readRDS(file=here("results/sim_results_interim.rds"))
sim_results_pooled <- sim_results_full %>% filter(study=="pooled")

study_results <- sim_results_full %>% filter(study!="pooled")
#calculate Pvalue from effect and se
study_results$pvalue <- 2 * (1 - pnorm(abs(study_results$effect/ study_results$se)))

d <- left_join(study_results, truth, by=c("study", "group"))

unique(study_results$group)
unique(truth$group)

d <- d %>% rename(estimate=effect, std.err=se, trueRD=IR_diff, trueRR=IR_ratio) %>%
  mutate(lower=estimate-1.96*std.err, upper=estimate+1.96*std.err)

head(d)

# Calculate performance metrics

#Need to add RD to estimation
# resRD <- d %>% 
#     summarise(abs_bias_RD=mean(abs(estimate-trueRD)),
#               estimator_variance_RD=mean(((estimate)-mean((estimate)))^2),
#               mean_variance_RD=mean((std.err)^2),
#               bias_se_ratio_RD=abs_bias_RD/sqrt(mean_variance_RD),
#               coverage_RD=mean(lower<=trueRD & trueRD<=upper)*100,
#               O_coverage_RD=mean(estimate-1.96*sd(estimate)< trueRD & trueRD < estimate+1.96*sd(estimate))*100
#     )


  
resRR <- d %>% group_by(group) %>%
    summarise(est_RR=exp(mean(estimate)),
              trueRR=exp(mean(log(trueRR))),
              abs_log_bias_RR=mean(abs(estimate-log(trueRR))),
              estimator_variance_RR=mean(((estimate)-mean((estimate)))^2),
              mean_variance_RR=mean((std.err)^2),
              bias_se_ratio_RR=abs_log_bias_RR/sqrt(mean_variance_RR),
              power = mean(pvalue<0.05) * 100,
              coverage_RR=mean(lower<=log(trueRR) & log(trueRR)<=(upper))*100,
              O_coverage_RR=mean((estimate)-1.96*sd((estimate))< log(trueRR) & 
                                   log(trueRR) < (estimate)+1.96*sd((estimate)))*100)
resRR

resRR_study <- d %>% group_by(group, study) %>%
  summarise(est_RR=exp(mean(estimate)),
            trueRR=exp(mean(log(trueRR))),
            abs_log_bias_RR=mean(abs(estimate-log(trueRR))),
            estimator_variance_RR=mean(((estimate)-mean((estimate)))^2),
            mean_variance_RR=mean((std.err)^2),
            power = mean(pvalue<0.05) * 100,
            bias_se_ratio_RR=abs_log_bias_RR/sqrt(mean_variance_RR),
            coverage_RR=mean(lower<=log(trueRR) & log(trueRR)<=(upper))*100,
            O_coverage_RR=mean((estimate)-1.96*sd((estimate))< log(trueRR) & 
                                 log(trueRR) < (estimate)+1.96*sd((estimate)))*100)
resRR_study

#-------------------------------------------------------------------------------
# plot_performance
#-------------------------------------------------------------------------------

#transform data to long format
resRR_study_long <- resRR_study %>% 
  select(group, study,  abs_log_bias_RR, estimator_variance_RR, 
         power, bias_se_ratio_RR, coverage_RR, O_coverage_RR ) %>%
  gather(key="metric", value="value", -group, -study)

ggplot(resRR_study_long, aes(x=study, y=value, color=group, shape=group)) +
  coord_flip() +
  geom_point() +
  facet_wrap(~metric, scales="free") +
  theme_minimal()
  
#tab = merge(resRD, resRR, by=c("Estimator","estimator"))
  

