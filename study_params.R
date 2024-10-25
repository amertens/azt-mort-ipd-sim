
# Define study parameters with CVs
study_params <- list(
  MORDOR_I = list(
    n_clusters = 767,
    cluster_size = 124,
    baseline_rate = 0.0165,
    effect_size = 0.865,
    cv = 0.34            # Directly reported CV
  ),
  MORDOR_II = list(
    n_clusters = 594,
    cluster_size = 119,
    baseline_rate = 0.0263,
    effect_size = 0.84,
    cv = 0.044           # Back-calculated CV
  ),
  TANA_I = list(
    n_clusters = 24,
    cluster_size = 89,
    baseline_rate = 0.0083,
    effect_size = 0.49,
    cv = 0.232           # Back-calculated CV
  ),
  AVENIR = list(
    n_clusters = 1455,
    cluster_size = 95, 
    baseline_rate = 0.0139,
    effect_size = 0.86,
    cv = 0.032           # Back-calculated CV
  )
)