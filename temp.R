
data = sim_study_data %>% filter(study == "TANA_I")
subgroup = "age_group"     # "age_group", "sex", or NULL
subgroup_level = "1-5mo" # specific level for subgroup
adjusted = TRUE




# Example usage:
# Overall analysis
overall_results <- analyze_mortality(data, adjusted = TRUE, run_tmle=TRUE)

# Age subgroup analysis
age_results <- lapply(unique(data$age_group), function(age) {
  analyze_mortality(data, subgroup = "age_group", subgroup_level = age, adjusted = TRUE)
})
names(age_results) <- unique(data$age_group)
age_results=data.table::rbindlist(age_results, use.names = TRUE, idcol = "level") %>% mutate(group="age_group")

# Sex subgroup analysis
sex_results <- lapply(c(0,1), function(sex) {
  analyze_mortality(data, subgroup = "sex", subgroup_level = sex, adjusted = TRUE)
})
names(sex_results) <- c("female", "male")
sex_results=data.table::rbindlist(sex_results, use.names = TRUE, idcol = "level") %>% mutate(group="sex")

