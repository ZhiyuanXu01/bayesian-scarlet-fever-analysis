# Results

This directory contains selected aggregate outputs that do not disclose
individual-level study records.

- `analysis_overview.csv`: analytical sample dimensions and date coverage;
- `annual_case_counts.csv`: annual aggregate case counts;
- `monthly_case_counts.csv`: monthly aggregate case counts and trailing
  12-month means;
- `case_characteristics.csv`: aggregate demographic, geographic, and seasonal
  summaries;
- `exposure_summary.csv`: descriptive statistics for current-day environmental
  exposures.
- `conditional_poisson_models.csv`: posterior relative-risk estimates from 189
  single-exposure Bayesian conditional Poisson models;
- `pooled_bdlnm_metrics.csv`: fit diagnostics and computational information for
  18 pooled Bayesian distributed lag nonlinear models;
- `pooled_bdlnm_slices.csv`: posterior nonlinear exposure-response estimates at
  single lag 3 and moving average 4, based on 500 posterior samples;
- `pooled_bdlnm_surfaces.csv`: posterior-mean exposure-lag-response surfaces
  across all nine environmental exposures.
- `interaction_spearman_screen.csv`: Spearman-correlation screening of 72
  prespecified exposure pairs at single lag 3 and moving average 4;
- `interaction_models.csv`: posterior estimates from 64 IQR-scaled Bayesian
  interaction models retained after correlation screening.
- `subgroup_sample_sizes.csv`: aggregate analytical-row and matched-stratum
  counts for 12 prespecified population, age, sex, period, and seasonal groups;
- `subgroup_models.csv`: posterior estimates from 216 Bayesian conditional
  Poisson subgroup models at single lag 3 and moving average 4;
- `subgroup_iqr_scaled_estimates.csv`: the same subgroup estimates rescaled to
  each subgroup's exposure IQR for comparable public presentation.
- `sensitivity_linear_priors.csv`: estimates from 54 fixed-window models under
  three PC-prior specifications;
- `sensitivity_bdlnm_metrics.csv`: model diagnostics for 54 pooled B-DLNM
  structural sensitivity models;
- `sensitivity_bdlnm_slices.csv`: posterior-mean exposure-response slices under
  alternative maximum-lag, lag-spline, and exposure-knot specifications;
- `sensitivity_bdlnm_quantile_contrasts.csv`: 144 P10-versus-P50 and
  P90-versus-P50 contrasts comparing the primary and alternative pooled
  B-DLNM specifications.

All files in this directory are aggregate outputs selected for public release.
Restricted row-level data and fitted model objects are not included.
