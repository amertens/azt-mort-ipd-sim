
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
source("simulation_functions.R")
source("study_params.R")

truth_FE <- readRDS(here("results/meta_truth_FE.rds"))
truth_RE <- readRDS(here("results/meta_truth_RE.rds"))

results_RE_2step <- readRDS(file=here("results/pooled_results_RE.rds")) %>% filter(!(metric %in% c("crude_rr", "hr")))
results_FE_2step <- readRDS(file=here("results/pooled_results_FE.rds")) %>% filter(!(metric %in% c("crude_rr", "hr")))

results_FE_1step <- readRDS(file=here("results/sim_results_interim_par_1step_FE.rds")) %>%  filter(level=="main") %>%
  mutate(study="pooled", group="main", level="main", adjusted=TRUE, crude_rr=NA, hr=NA)

results_RE_1step <- readRDS(file=here("results/sim_results_interim_par_1step_RE.rds")) %>% 
  mutate(study="pooled", group="main", level="main", adjusted=TRUE, hr=NA, tmle_log_rr=NULL, tmle_ate=NULL, hr_se=NA, tmle_rr_log_se=NULL, tmle_ate_se=NULL,hr_pval=NA, tmle_rr_pval=NULL, tmle_ate_pval=NULL)
results_RE_tmle_1step <- readRDS(file=here("results/sim_results_par_1step_tmle_RE.rds")) %>% rename(iteration_tmle=iteration) %>%
  mutate(tmle_rr_log_se =(tmle_rr_log_se))
results_RE_1step <- bind_cols(results_RE_1step, results_RE_tmle_1step)



results_RE_1step <- clean_sim_results(results_RE_1step)
results_FE_1step <- clean_sim_results(results_FE_1step)



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



d_FE_2step <- left_join(results_FE_2step, truth_FE, by=c("level", "metric"))
d_RE_2step <- left_join(results_RE_2step, truth_RE, by=c("level", "metric"))
d_FE_1step <- left_join(results_FE_1step, truth_FE, by=c("level", "metric"))
d_RE_1step <- left_join(results_RE_1step, truth_RE, by=c("level", "metric"))

d_RE_1step %>% filter(metric=="cid") %>% summarise(mean(lower<true_value & true_value<upper))
d_RE_1step %>% filter(metric=="tmle_ate") %>% summarise(mean(lower<true_value & true_value<upper))
d_RE_1step %>% filter(metric=="tmle_log_rr") %>% summarise(mean(lower<true_value & true_value<upper))

res_FE_2step <- d_FE_2step %>% 
  filter(level=="main") %>%
  group_by(level, metric, adjusted) %>% 
  summarise(n=n(),
            true_value=true_value[1],
            mean_value=mean(effect),
            abs_bias=mean(abs(effect - true_value)),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100) %>%
  filter(!is.na(abs_bias), metric %in% c("cid", "cir","tmle_ate","tmle_log_rr")) %>%
  mutate(pooling="2-step FE")

res_RE_2step <- d_RE_2step %>% group_by(level, metric, adjusted) %>%
  filter(level=="main") %>%
  summarise(n=n(),
            true_value=true_value[1],
            mean_value=mean(effect),
            abs_bias=mean(abs(effect - true_value)),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100) %>%
  filter(!is.na(abs_bias), metric %in% c("cid", "cir","tmle_ate","tmle_log_rr")) %>%
  mutate(pooling="2-step RE")

res_FE_1step <- d_FE_1step %>% group_by(level, metric, adjusted) %>% 
  filter(level=="main") %>%
  summarise(n=n(),
            true_value=true_value[1],
            mean_value=mean(effect),
            abs_bias=mean(abs(effect - true_value)),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100) %>%
  filter(!is.na(abs_bias), metric %in% c("cid", "cir","tmle_ate","tmle_log_rr")) %>%
  mutate(pooling="1-step FE")

res_RE_1step <- d_RE_1step %>% group_by(level, metric, adjusted) %>%
  filter(level=="main") %>%
  summarise(n=n(),
            true_value=true_value[1],
            mean_value=mean(effect),
            abs_bias=mean(abs(effect - true_value[1])),
            estimator_variance=mean(((effect )-mean((effect )))^2),
            mean_variance=mean((se)^2),
            bias_se_ratio=abs_bias/sqrt(mean_variance),
            power=mean((lower<0 & upper<0)|(lower>0 & upper>0))*100,
            coverage=mean(lower<=true_value & true_value<=upper)*100,
            O_coverage=mean(effect -1.96*sd(effect )< true_value & true_value < effect +1.96*sd(effect ))*100) %>%
  filter(!is.na(abs_bias), metric %in% c("cid", "cir","tmle_ate","tmle_log_rr")) %>%
  mutate(pooling="1-step RE")

res <- bind_rows(res_FE_2step, res_RE_2step, res_FE_1step, res_RE_1step)

res %>% filter(metric=="cid"|metric=="tmle_ate") %>% select(pooling, everything())
res %>% filter(metric=="cir"|metric=="tmle_log_rr") %>% select(pooling, everything())

res %>% filter(metric=="cid"|metric=="tmle_ate", adjusted, grepl("RE",pooling)) %>% select(pooling, everything())
res %>% filter(metric=="cir"|metric=="tmle_log_rr", adjusted, grepl("RE",pooling)) %>% select(pooling, everything())

tab <- res %>% ungroup() %>%
  filter(level =="main", adjusted==TRUE) %>% 
  mutate(conversion=ifelse(metric %in% c("cid", "tmle_ate"),1000,100)) %>%
  mutate( abs_bias=abs_bias*conversion, mean_variance=mean_variance*conversion) %>%
  select(pooling, metric, abs_bias, mean_variance, bias_se_ratio, power, coverage, O_coverage) %>%
  mutate(mean_variance =round(mean_variance ,3), bias_se_ratio  =round(bias_se_ratio  ,3)) %>%
  mutate(metric=case_when(
    metric=="cid" ~ "GLM cumulative incidence difference",
    metric=="cir" ~ "GLM cumulative incidence ratio",
    metric=="tmle_ate" ~ "TMLE cumulative incidence difference",
    metric=="tmle_log_rr" ~ "TMLE cumulative incidence ratio"
  )) %>%
  rename(Pooling=pooling, Parameter=metric, `Absolute bias`=abs_bias , `Mean variance`=mean_variance , `Bias/SE ratio`=bias_se_ratio,
         `Power`=power , `CI coverage`=coverage , `Oracle coverage`=O_coverage)
tab 

tab %>% filter(grepl("GLM",Parameter),grepl("FE",Pooling))
tab %>% filter(grepl("TMLE",Parameter),grepl("FE",Pooling))
tab %>% filter(grepl("GLM",Parameter),grepl("RE",Pooling))
tab %>% filter(grepl("TMLE",Parameter),grepl("RE",Pooling))

saveRDS(tab, file=here::here("results/one_step_comparison.RDS"))
