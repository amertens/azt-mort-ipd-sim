
# Plot results
library(tidyverse)
library(metafor)
library(here)
source("meta_functions.R")

sim_results_full<-readRDS(file=here("results/sim_results_interim.rds"))
sim_results_pooled <- sim_results_full %>% filter(study=="pooled")
study_results <- sim_results_full %>% filter(study!="pooled")

study_results <- study_results %>% group_by(study) %>% mutate(med_effect=median(exp(effect))) 
table(study_results$med_effect)

ggplot(study_results, 
       aes(x=exp(effect), y=se, color=group)) +
  geom_point() +
  facet_wrap(~study) +
  scale_x_continuous(trans="log10") +
  coord_cartesian(xlim=c(0.1, 2)) +
  geom_vline(xintercept=1, linetype="dashed") +
  geom_vline(aes(xintercept=med_effect), linetype="dashed", color="red") +
  labs(title="Study Results",
       x="Effect Size",
       y="Standard Error") +
  theme_minimal()


power_results <- sim_results_pooled %>% group_by(group) %>% 
  summarize(power=mean(pvalue<0.05, na.rm=TRUE), 
            mean_effect=mean(effect, na.rm=TRUE), 
            mean_se=mean(se, na.rm=TRUE), 
            mean_tau2=mean(tau2, na.rm=TRUE))

#need to add oracle coverage.

main <- sim_results_pooled[,,1]
age <- sim_results_pooled[,,2]

head(main)
head(age)

# Convert to nice data frame
power_summary <- do.call(rbind, power_results)
power_summary$analysis <- rownames(power_summary)


plotdf <- power_summary

# ggplot(plotdf, aes(x=analysis, y=power)) +
#   geom_bar(stat="identity") +
#   coord_flip() +
#   labs(title="Statistical Power by Analysis",
#        x="Analysis", 
#        y="Power") +
#   theme_minimal()
# 
# # Effect size comparison
# ggplot(plotdf, aes(x=analysis, y=mean_effect)) +
#   geom_point() +
#   geom_errorbar(aes(ymin=mean_effect-1.96*mean_se, 
#                     ymax=mean_effect+1.96*mean_se)) +
#   coord_flip() +
#   labs(title="Effect Sizes by Analysis",
#        x="Analysis",
#        y="Risk Ratio") +
#   theme_minimal()


#TO DO:
#Need to compare doing a traditional meta-analysis with an IPD meta-analysis
#Need to add in subgroup analyses
#Need to come up with a list of secondary analyses (AMR)
#Need to come up with a list of sensitivity analyses



res = study_results %>% group_by(iteration, group) %>%
  do(nb_meta(log_irr=.$effect, se_log_irr=.$se, study_labels=.$study, method="REML"))

res

power_results <- res %>% group_by(group) %>% 
  summarize(power=mean(pval<0.05, na.rm=TRUE), 
            mean_effect=exp(mean(log(irr), na.rm=TRUE))#, 
            # mean_se=mean(se, na.rm=TRUE), 
            # mean_tau2=mean(tau2, na.rm=TRUE)  $add to function
            )
power_results

res_FE = study_results %>% group_by(iteration, group) %>%
  do(nb_meta(log_irr=.$effect, se_log_irr=.$se, study_labels=.$study, method="FE"))

res_FE


power_results_FE <- res_FE %>% group_by(group) %>% 
  summarize(power=mean(pval<0.01, na.rm=TRUE), 
            mean_effect=exp(mean(log(irr), na.rm=TRUE)))
power_results_FE

saveRDS(res, file=here("results/pooled_res.rds"))
saveRDS(res_FE, file=here("results/pooled_res_FE.rds"))



