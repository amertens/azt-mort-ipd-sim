
#note: multiply all cluster_size by 1.1 because these numbers are based on 

loss_rate_correction = 1.1

# Define study parameters with multi-site trials split
study_params <- list(
  MORDOR_Niger = list(
    n_clusters = 297,    # About 1/3 of total 1533 communities
    cluster_size = floor(128*loss_rate_correction),  
    baseline_rate = 0.0275,  # 27.5 deaths per 1000 person-years
    effect_size = 0.819,     # 18.1% reduction (95% CI: 10.0-25.5)
    cv = 0.34                # Reported CV for Niger
  ),
  
  MORDOR_Malawi = list(
    n_clusters = 152,        # About 1/3 of communities
    cluster_size = floor(205*loss_rate_correction),      
    baseline_rate = 0.0096,  # 9.6 deaths per 1000 person-years
    effect_size = 0.943,     # 5.7% reduction (95% CI: -9.7 to 18.9)
    cv = 0.34                # Using Niger CV as not reported
  ),
  
  MORDOR_Tanzania = list(
    n_clusters = 307,        # About 1/3 of communities
    cluster_size = floor(58*loss_rate_correction),      
    baseline_rate = 0.0055,  # 5.5 deaths per 1000 person-years 
    effect_size = 0.966,     # 3.4% reduction (95% CI: -21.2 to 23.0)
    cv = 0.34                # Using Niger CV as not reported
  ),
  
  MORDOR_II_Niger = list(
    n_clusters = 594,
    cluster_size = floor(119*loss_rate_correction),
    baseline_rate = 0.0263,  # 26.3 deaths per 1000 person-years
    effect_size = 0.84,      # 16% reduction in year 1
    cv = 0.044               # Back-calculated CV
  ),
  
  TANA_I = list(
    n_clusters = 50,
    cluster_size = floor(79*loss_rate_correction),
    baseline_rate = 0.0083,  # 8.3 deaths per 1000 person-years
    effect_size = 0.49,      # 51% reduction
    cv = 0.232               # Back-calculated CV
  ),
  
  AVENIR = list(
    n_clusters = 463,
    cluster_size = floor(95*loss_rate_correction), 
    baseline_rate = 0.0139,  # 13.9 deaths per 1000 person-years
    effect_size = 0.86,      # 14% reduction
    cv = 0.032               # Back-calculated CV
  ),
  
  #Note this is bi-annual vs annual
  PRET_biannual = list(
    n_clusters = 24,
    cluster_size = floor(107*loss_rate_correction),
    baseline_rate = 0.0290,  # 29.0 deaths per 1000 person-years
    effect_size = 0.81,      # From reported IRR
    cv = 0.15                # Using same as annual arm
  )
)