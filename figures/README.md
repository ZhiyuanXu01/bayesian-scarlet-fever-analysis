# Figures

This directory contains selected aggregate figures suitable for public
release.

- `monthly_case_counts.png`: monthly reported case counts with a trailing
  12-month mean.
- `conditional_poisson_single_lag_rr.png`: relative-risk estimates for
  single-day lags 0-10;
- `conditional_poisson_moving_average_rr.png`: relative-risk estimates for
  cumulative moving-average windows 1-10.
- `pooled_bdlnm_single_lag3_slices.png`: nonlinear exposure-response curves at
  single lag 3 with 95% credible intervals;
- `pooled_bdlnm_mavg4_slices.png`: nonlinear exposure-response curves at moving
  average 4 with 95% credible intervals.
- `spatial_postprob_pm10_lag0.png`: county-level posterior probability of
  increased risk for PM10 at lag 0, comparing P90 with P50;
- `spatial_postprob_surface_pressure_lag3.png`: county-level posterior
  probability of increased risk for surface pressure at lag 3, comparing P90
  with P50.
- `interaction_ratio_heatmap.png`: IQR-scaled multiplicative interaction
  ratios for exposure pairs at single lag 3 and moving average 4; grey cells
  identify pairs excluded by the prespecified correlation screen.
- `subgroup_iqr_scaled_rr_heatmap.png`: subgroup-specific relative risks per
  exposure IQR across age, sex, period, and seasonal groups; matched strata
  were assigned from the case record and retained with all referent days.
- `sensitivity_bdlnm_quantile_contrasts.png`: posterior central relative-risk
  estimates at exposure P10 and P90 under the primary pooled B-DLNM and three
  alternative structural specifications.

The spatial probability maps are aggregate outputs from the original
high-memory SB-DLNM run. Full spatial models were not refitted on the local
desktop during repository restructuring.

All figures in this directory contain aggregate results selected for public
release.
