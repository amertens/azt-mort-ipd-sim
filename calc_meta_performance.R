
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

truth_FE <- readRDS(here("results/meta_truth_FE.rds"))
truth_RE <- readRDS(here("results/meta_truth_RE.rds"))

results_RE <- readRDS(file=here("results/pooled_results_RE.rds")) %>% filter(!(metric %in% c("crude_rr", "hr")))
results_FE <- readRDS(file=here("results/pooled_results_FE.rds")) %>% filter(!(metric %in% c("crude_rr", "hr")))

truth_FE <- truth_FE %>% mutate(true_value = case_when(
  metric %in% c("cid", "ird") ~ true_value,
  metric %in% c("cir", "irr") ~ log(true_value)
)) %>% filter(!is.na(metric))
truth_RE <- truth_RE %>% mutate(true_value = case_when(
  metric %in% c("cid", "ird") ~ true_value,
  metric %in% c("cir", "irr") ~ log(true_value)
)) %>% filter(!is.na(metric))


truth_tmle_FE <- truth_FE %>% mutate(metric=case_when(
  metric=="cid" ~ "tmle_ate",
  metric=="cir" ~ "tmle_log_rr"
)) %>% filter(!is.na(metric))
truth_FE <- bind_rows(truth_FE, truth_tmle_FE)

truth_tmle_RE <- truth_RE %>% mutate(metric=case_when(
  metric=="cid" ~ "tmle_ate",
  metric=="cir" ~ "tmle_log_rr"
)) %>% filter(!is.na(metric))
truth_RE <- bind_rows(truth_RE, truth_tmle_RE)

table(results_FE$level)
table(results_FE$metric)


d_FE <- left_join(results_FE, truth_FE, by=c("level", "metric"))
d_RE <- left_join(results_RE, truth_RE, by=c("level", "metric"))


res_FE <- d_FE %>% group_by(level, metric, adjusted) %>%
  summarise(n=n(),
            true_value=true_value[1],
            abs_bias=mean(abs(effect - true_value)),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100)
res_RE <- d_RE %>% group_by(level, metric, adjusted) %>%
  summarise(n=n(),
            true_value=true_value[1],
            abs_bias=mean(abs(effect - true_value)),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100)


#need to get 100 iterations for all
d_RE %>% group_by(level,metric, adjusted) %>% summarize(mean(iteration))
  
res_RE_summarized <- d_RE %>% group_by(level,metric, adjusted) %>% 
  #TEMP
  filter(iteration <901) %>%
  summarise(n=n(),
            true_value=true_value[1],
            abs_bias=mean(abs(effect - true_value)),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100) %>%
  group_by(metric, adjusted) %>%
  summarise(n=sum(n),
            abs_bias=mean(abs_bias),
            estimator_variance=mean(estimator_variance),
            mean_variance=mean(mean_variance),
            bias_se_ratio=mean(bias_se_ratio),
            power=mean(power),
            coverage=mean(coverage),
            O_coverage=mean(O_coverage)) %>%
  filter(!is.na(abs_bias))


#transform data to long format
resRR_study_long_FE <- res_FE %>% 
  select(level, metric, adjusted, abs_bias, estimator_variance, mean_variance, bias_se_ratio, power, coverage, O_coverage ) %>%
  gather(key="performance_metric", value="value", -level, -metric, -adjusted) %>% filter(!is.na(value))
resRR_study_long_RE <- res_RE %>% 
  select(level, metric, adjusted, abs_bias, estimator_variance, mean_variance, bias_se_ratio, power, coverage, O_coverage ) %>%
  gather(key="performance_metric", value="value", -level, -metric, -adjusted) %>% filter(!is.na(value))

#-------------------------------------------------------------------------------
# tabulate performance
#-------------------------------------------------------------------------------
tabRD_FE <- resRR_study_long_FE %>% filter( metric %in% c("cid","tmle_ate")) %>% 
  spread(key="performance_metric", value="value") %>% arrange(level, abs(O_coverage-95), bias_se_ratio )
tabRD_FE

tabRR_FE <- resRR_study_long_FE %>% filter( metric %in% c("cir","tmle_log_rr")) %>% 
  spread(key="performance_metric", value="value") %>% arrange(level, abs(O_coverage-95), bias_se_ratio )
tabRR_FE

tabRD_RE <- resRR_study_long_RE %>% filter( metric %in% c("cid","tmle_ate")) %>% 
  spread(key="performance_metric", value="value") %>% arrange(level, abs(O_coverage-95), bias_se_ratio )
tabRD_RE

tabRR_RE <- resRR_study_long_RE %>% filter( metric %in% c("cir","tmle_log_rr")) %>% 
  spread(key="performance_metric", value="value") %>% arrange(level, abs(O_coverage-95), bias_se_ratio )
tabRR_RE

tabRD_FE %>% filter(level=="main")
tabRR_FE %>% filter(level=="main")
tabRD_RE %>% filter(level=="main")
tabRR_RE %>% filter(level=="main")

#-------------------------------------------------------------------------------
# plot_performance
#-------------------------------------------------------------------------------

#transform data to long format
resRR_long <- tabRD_RE %>% 
  select(level, metric,  abs_bias, estimator_variance, mean_variance, bias_se_ratio, coverage, O_coverage, adjusted) %>%
  gather(key="performance_metric", value="value", -level, -metric,, -adjusted)
resRR_long$metric[resRR_long$metric=="cid" & resRR_long$adjusted ] <- "GLM"
resRR_long$metric[resRR_long$metric=="cid" & !resRR_long$adjusted ] <- "Unadjusted GLM"
resRR_long$metric[resRR_long$metric=="tmle_ate" & resRR_long$adjusted ] <- "TMLE"

resRR_long <- resRR_long %>% filter(performance_metric!="estimator_variance") %>%
  mutate(level=factor(level, levels=c("24-59mo",  "12-23mo", "6-11mo" ,"1-5mo", "female","male","main"))) 

ggplot(resRR_long ,
       aes(x=level, y=value, color=metric, shape=metric  )) +
  coord_flip() +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~performance_metric, scales="free") +
  theme_minimal() + theme(legend.position = c(0.8,0.25)) + 
  guides(color=guide_legend(title="Estimator"), shape=guide_legend(title="Estimator"))


save(tabRD_FE, tabRR_FE, tabRD_RE, tabRR_RE, res_RE_summarized, resRR_long, file=here("results/meta_performance.rdata"))
