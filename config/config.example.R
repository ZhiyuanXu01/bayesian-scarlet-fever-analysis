# Example project configuration
#
# Copy this file to config/config.R and edit the local paths.
# config/config.R is excluded from version control.

analysis_config <- list(
  data = list(
    analysis_file = file.path(
      "data",
      "restricted",
      "scarlet_fever_analysis.csv"
    ),
    hubei_boundary = file.path(
      "data",
      "restricted",
      "hubei_counties.geojson"
    ),
    hunan_boundary = file.path(
      "data",
      "restricted",
      "hunan_counties.geojson"
    )
  ),

  output = list(
    results_dir = "results",
    figures_dir = "figures"
  ),

  computation = list(
    inla_num_threads = "4:1",
    posterior_samples = 500L,
    random_seed = 2026L
  )
)