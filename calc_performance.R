
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
#study_results <- readRDS(file=here("results/sim_results_interim_par.rds"))
study_results <- readRDS(here("results/sim_results_clean.rds"))


#join with truth
table(study_results$metric)
table(truth$metric)

truth <- truth %>% mutate(true_value = case_when(
  metric %in% c("cid", "ird") ~ true_value,
  metric %in% c("cir", "irr") ~ log(true_value)
)) %>% filter(!is.na(metric))
truth

truth_tmle <- truth %>% mutate(metric=case_when(
  metric=="cid" ~ "tmle_ate",
  metric=="cir" ~ "tmle_log_rr"
)) %>% filter(!is.na(metric))
truth <- bind_rows(truth, truth_tmle)

table(study_results$study)
table(study_results$group)
table(study_results$level)
table(study_results$metric)
table(truth$study)
table(truth$group)
table(truth$metric)

d <- left_join(study_results, truth, by=c("study", "level", "metric"))
table(d$metric,!is.na(d$true_value))
table(d$level,!is.na(d$true_value))
# Calculate performance metrics
tail(d)

summary( d$effect[d$metric=="tmle_log_rr"] )
truth

res <- d %>% group_by(study,level, metric, adjusted) %>%
  summarise(n=n(),
            true_value=true_value[1],
            abs_bias=mean(abs(effect - true_value)),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100)
head(res)

#temp drop HR
res <- res %>% filter(metric!="hr", metric!="crude_rr")


#transform data to long format
resRR_study_long <- res %>% 
  select(level, metric, study, adjusted,  abs_bias, estimator_variance, mean_variance, bias_se_ratio, power, coverage, O_coverage ) %>%
  gather(key="performance_metric", value="value", -level, -metric, -study)

#-------------------------------------------------------------------------------
# tabulate performance
#-------------------------------------------------------------------------------
tabRD <- resRR_study_long %>% filter(level=="main", metric %in% c("cid","tmle_ate")) %>% 
  spread(key="performance_metric", value="value") %>% arrange(abs(O_coverage-95), bias_se_ratio )
tabRD

tabRR <- resRR_study_long %>% filter(level=="main", metric %in% c("cir","tmle_log_rr")) %>% 
  spread(key="performance_metric", value="value") %>% arrange(abs(O_coverage-95), bias_se_ratio )
tabRR


#-------------------------------------------------------------------------------
# plot performance
#-------------------------------------------------------------------------------



ggplot(resRR_study_long %>% filter(level=="main"),
         aes(x=study, y=value, color=metric  , shape=metric  )) +
  coord_flip() +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~performance_metric, scales="free") +
  theme_minimal()

temp <- resRR_study_long %>% filter(level=="main", metric %in% c("cid","tmle_ate"))

ggplot(resRR_study_long %>% filter(level=="main", metric %in% c("cid","tmle_ate")),
       aes(x=study, y=value, color=metric  , shape=metric  )) +
  coord_flip() +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~performance_metric, scales="free") +
  theme_minimal()

ggplot(resRR_study_long %>% filter(level=="main", !(metric %in% c("cid","tmle_ate"))),
       aes(x=study, y=value, color=metric  , shape=metric  )) +
  coord_flip() +
  geom_point(position = position_dodge(0.5)) +
  facet_wrap(~performance_metric, scales="free") +
  theme_minimal()

unique(resRR_study_long$level)
resRR_study_long<-resRR_study_long %>% mutate(level=factor(level, levels=c("24-59mo",  "12-23mo", "6-11mo" ,"1-5mo", "female","male","main" )))

temp=resRR_study_long %>% filter(performance_metric=="O_coverage")
ggplot(resRR_study_long %>% filter(performance_metric=="O_coverage"), aes(x=level, y=value, color=metric  , shape=metric  )) +
  coord_flip() +
  geom_hline(yintercept = 95, linetype="dashed") +
  geom_point(position=position_dodge(0.5)) +
  facet_wrap(~study, scales="fixed") +
  theme_minimal()

ggplot(resRR_study_long %>% filter(performance_metric=="coverage"), aes(x=level, y=value, color=metric  , shape=metric  )) +
  coord_flip() +
  geom_hline(yintercept = 95, linetype="dashed") +
  geom_point(position=position_dodge(0.5)) +
  facet_wrap(~study, scales="fixed") +
  theme_minimal()


ggplot(resRR_study_long %>% filter(performance_metric=="abs_bias", metric %in% c("cid","tmle_ate")), 
       aes(x=level, y=value, color=metric  , shape=metric  )) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_point(position=position_dodge(0.5)) +
  facet_wrap(~study, scales="fixed") +
  theme_minimal()

ggplot(resRR_study_long %>% filter(performance_metric=="bias_se_ratio", metric %in% c("cir","tmle_log_rr")),
       aes(x=level, y=value, color=metric  , shape=metric  )) +
  coord_flip() +
  geom_point(position=position_dodge(0.5)) +
  facet_wrap(~study, scales="fixed") +
  theme_minimal()
  
#tab = merge(resRD, resRR, by=c("Estimator","estimator"))
  

