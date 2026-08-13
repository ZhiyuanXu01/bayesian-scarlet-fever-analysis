# Computational requirements

## Locally validated components

The following components were validated on a standard desktop environment:

- restricted analytical-data loading and validation;
- reconstruction of moving-average exposure variables;
- descriptive analysis;
- 189 Bayesian conditional Poisson models;
- 18 pooled Bayesian distributed lag nonlinear models;
- county-boundary processing, point-to-county assignment, and adjacency graph
  construction;
- spatial B-DLNM design-matrix construction.

The spatial preparation retained 128,508 of 128,557 analytical rows. The 49
records that did not intersect a study county boundary are excluded only from
the spatial models.

## Full spatial B-DLNM

A spatial B-DLNM contains 12 BYM2 spatially varying cross-basis components in
addition to the matched-stratum random effect. The resulting models are
substantially more memory intensive than the pooled B-DLNMs.

The full spatial analysis used a high-memory cloud environment with 128 GB of
RAM. It is intentionally not executed by the default local workflow. Users
with an appropriate environment can call `run_spatial_bdlnm_model()` for one
exposure specification at a time.

Posterior configuration, CPO computation, and fitted-value storage increase
memory requirements and are disabled by default. They should be enabled only
when required for posterior maps or formal model comparison.

No model object or individual-level spatial output is included in the public
repository.
