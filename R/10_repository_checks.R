# ==============================================================================
# Public repository integrity checks
#
# These checks do not read restricted analytical data. They validate public
# scripts, aggregate outputs, dependency metadata, and tracked-file privacy.
# ==============================================================================

check_results <- data.frame(
  check = character(),
  passed = logical(),
  details = character(),
  stringsAsFactors = FALSE
)

record_check <- function(
    check,
    passed,
    details = ""
) {
  check_results <<- rbind(
    check_results,
    data.frame(
      check = check,
      passed = isTRUE(passed),
      details = as.character(details),
      stringsAsFactors = FALSE
    )
  )

  invisible(passed)
}

read_public_csv <- function(
    path
) {
  if (!file.exists(path)) {
    stop(
      "Missing public result file: ",
      path,
      call. = FALSE
    )
  }

  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

expected_files <- c(
  "README.md",
  "DATA_AVAILABILITY.md",
  "COMPUTE_REQUIREMENTS.md",
  "LICENSE",
  "renv.lock",
  "config/config.example.R",
  "data/data_dictionary.csv",
  file.path(
    "R",
    sprintf(
      "%02d_%s.R",
      1:10,
      c(
        "data_contract",
        "prepare_analysis_data",
        "descriptive_analysis",
        "bayesian_conditional_poisson",
        "pooled_bdlnm",
        "spatial_bdlnm",
        "interaction_analysis",
        "subgroup_analysis",
        "sensitivity_analysis",
        "repository_checks"
      )
    )
  ),
  "results/conditional_poisson_models.csv",
  "results/pooled_bdlnm_metrics.csv",
  "results/pooled_bdlnm_slices.csv",
  "results/pooled_bdlnm_surfaces.csv",
  "results/interaction_spearman_screen.csv",
  "results/interaction_models.csv",
  "results/subgroup_sample_sizes.csv",
  "results/subgroup_models.csv",
  "results/subgroup_iqr_scaled_estimates.csv",
  "results/sensitivity_linear_priors.csv",
  "results/sensitivity_bdlnm_metrics.csv",
  "results/sensitivity_bdlnm_slices.csv",
  "results/sensitivity_bdlnm_quantile_contrasts.csv",
  "figures/monthly_case_counts.png",
  "figures/conditional_poisson_single_lag_rr.png",
  "figures/conditional_poisson_moving_average_rr.png",
  "figures/pooled_bdlnm_single_lag3_slices.png",
  "figures/pooled_bdlnm_mavg4_slices.png",
  "figures/spatial_postprob_pm10_lag0.png",
  "figures/spatial_postprob_surface_pressure_lag3.png",
  "figures/interaction_ratio_heatmap.png",
  "figures/subgroup_iqr_scaled_rr_heatmap.png",
  "figures/sensitivity_bdlnm_quantile_contrasts.png"
)

missing_files <- expected_files[
  !file.exists(
    expected_files
  )
]

record_check(
  check = "Expected public files",
  passed = length(missing_files) == 0L,
  details = if (length(missing_files) == 0L) {
    paste(
      length(expected_files),
      "files present"
    )
  } else {
    paste(
      missing_files,
      collapse = "; "
    )
  }
)

r_files <- list.files(
  "R",
  pattern = "\\.R$",
  full.names = TRUE
)

parse_errors <- vapply(
  r_files,
  function(path) {
    tryCatch(
      {
        invisible(
          parse(
            file = path
          )
        )
        ""
      },
      error = function(error) {
        conditionMessage(error)
      }
    )
  },
  FUN.VALUE = character(1)
)

record_check(
  check = "R script syntax",
  passed = all(
    parse_errors == ""
  ),
  details = if (all(parse_errors == "")) {
    paste(
      length(r_files),
      "scripts parsed"
    )
  } else {
    paste(
      names(parse_errors)[
        parse_errors != ""
      ],
      parse_errors[
        parse_errors != ""
      ],
      sep = ": ",
      collapse = "; "
    )
  }
)

dictionary <- read_public_csv(
  "data/data_dictionary.csv"
)

record_check(
  check = "Data dictionary",
  passed =
    nrow(dictionary) == 203L &&
    sum(dictionary$storage == "input") == 109L &&
    sum(dictionary$storage == "derived") == 94L,
  details = paste(
    nrow(dictionary),
    "variables"
  )
)

conditional_models <- read_public_csv(
  "results/conditional_poisson_models.csv"
)

record_check(
  check = "Conditional Poisson models",
  passed =
    nrow(conditional_models) == 189L &&
    length(
      unique(
        conditional_models$model_key
      )
    ) == 189L &&
    all(
      conditional_models$status == "ok"
    ),
  details = paste(
    nrow(conditional_models),
    "successful models"
  )
)

pooled_metrics <- read_public_csv(
  "results/pooled_bdlnm_metrics.csv"
)

record_check(
  check = "Pooled B-DLNM models",
  passed =
    nrow(pooled_metrics) == 18L &&
    length(
      unique(
        pooled_metrics$model_key
      )
    ) == 18L &&
    all(
      pooled_metrics$status == "ok"
    ),
  details = paste(
    nrow(pooled_metrics),
    "successful models"
  )
)

interaction_screen <- read_public_csv(
  "results/interaction_spearman_screen.csv"
)

interaction_models <- read_public_csv(
  "results/interaction_models.csv"
)

record_check(
  check = "Interaction analysis",
  passed =
    nrow(interaction_screen) == 72L &&
    sum(
      interaction_screen$decision == "included"
    ) == 64L &&
    nrow(interaction_models) == 64L &&
    length(
      unique(
        interaction_models$model_key
      )
    ) == 64L &&
    all(
      interaction_models$status == "ok"
    ),
  details = paste(
    nrow(interaction_models),
    "successful models after screening"
  )
)

subgroup_sizes <- read_public_csv(
  "results/subgroup_sample_sizes.csv"
)

subgroup_models <- read_public_csv(
  "results/subgroup_models.csv"
)

subgroup_iqr <- read_public_csv(
  "results/subgroup_iqr_scaled_estimates.csv"
)

record_check(
  check = "Subgroup analysis",
  passed =
    nrow(subgroup_sizes) == 12L &&
    all(
      subgroup_sizes$
        strata_without_one_case_day == 0L
    ) &&
    nrow(subgroup_models) == 216L &&
    nrow(subgroup_iqr) == 216L &&
    length(
      unique(
        subgroup_models$model_key
      )
    ) == 216L &&
    all(
      subgroup_models$status == "ok"
    ),
  details = paste(
    nrow(subgroup_models),
    "successful models across",
    nrow(subgroup_sizes),
    "groups"
  )
)

linear_sensitivity <- read_public_csv(
  "results/sensitivity_linear_priors.csv"
)

bdlnm_sensitivity_metrics <- read_public_csv(
  "results/sensitivity_bdlnm_metrics.csv"
)

bdlnm_sensitivity_slices <- read_public_csv(
  "results/sensitivity_bdlnm_slices.csv"
)

bdlnm_contrasts <- read_public_csv(
  "results/sensitivity_bdlnm_quantile_contrasts.csv"
)

record_check(
  check = "Sensitivity analyses",
  passed =
    nrow(linear_sensitivity) == 54L &&
    all(
      linear_sensitivity$status == "ok"
    ) &&
    nrow(bdlnm_sensitivity_metrics) == 54L &&
    all(
      bdlnm_sensitivity_metrics$status == "ok"
    ) &&
    length(
      unique(
        bdlnm_sensitivity_slices$model_key
      )
    ) == 54L &&
    nrow(bdlnm_contrasts) == 144L,
  details = paste(
    "54 prior models, 54 structural models,",
    nrow(bdlnm_contrasts),
    "contrasts"
  )
)

figure_files <- expected_files[
  grepl(
    "\\.png$",
    expected_files
  )
]

figure_sizes <- file.info(
  figure_files
)$size

record_check(
  check = "Aggregate figures",
  passed =
    all(
      is.finite(
        figure_sizes
      )
    ) &&
    all(
      figure_sizes > 0
    ),
  details = paste(
    length(figure_files),
    "non-empty PNG files"
  )
)

lockfile_text <- readLines(
  "renv.lock",
  warn = FALSE
)

record_check(
  check = "Dependency lockfile",
  passed =
    any(
      grepl(
        '"Version": "4.4.1"',
        lockfile_text,
        fixed = TRUE
      )
    ) &&
    any(
      grepl(
        '"INLA"',
        lockfile_text,
        fixed = TRUE
      )
    ) &&
    any(
      grepl(
        "https://inla.r-inla-download.org/R/stable",
        lockfile_text,
        fixed = TRUE
      )
    ),
  details = "R 4.4.1 with named INLA repository"
)

tracked_files <- tryCatch(
  {
    output <- system2(
      "git",
      "ls-files",
      stdout = TRUE,
      stderr = TRUE
    )

    status <- attr(
      output,
      "status"
    )

    if (!is.null(status) &&
        status != 0L) {
      character()
    } else {
      output
    }
  },
  error = function(error) {
    character()
  }
)

record_check(
  check = "Git tracked-file inventory",
  passed = length(tracked_files) > 0L,
  details = paste(
    length(tracked_files),
    "tracked files"
  )
)

forbidden_tracked_patterns <- c(
  "^config/config\\.R$",
  "^data/(raw|restricted|private|interim|processed)/",
  "^results/private/",
  "\\.(rds|RDS|rda|RData|graph)$"
)

forbidden_tracked <- tracked_files[
  vapply(
    tracked_files,
    function(path) {
      any(
        vapply(
          forbidden_tracked_patterns,
          function(pattern) {
            grepl(
              pattern,
              path
            )
          },
          FUN.VALUE = logical(1)
        )
      )
    },
    FUN.VALUE = logical(1)
  )
]

forbidden_basenames <- c(
  "hb.csv",
  "hn.csv",
  "hb_sheet3.csv",
  "hb.tsv",
  "hn.tsv",
  "scarlet_fever_hb_hn_combined.csv"
)

forbidden_tracked <- unique(
  c(
    forbidden_tracked,
    tracked_files[
      tolower(
        basename(
          tracked_files
        )
      ) %in%
        forbidden_basenames
    ]
  )
)

record_check(
  check = "Restricted-file exclusion",
  passed = length(forbidden_tracked) == 0L,
  details = if (length(forbidden_tracked) == 0L) {
    "No restricted or fitted-model files are tracked"
  } else {
    paste(
      forbidden_tracked,
      collapse = "; "
    )
  }
)

tracked_sizes <- file.info(
  tracked_files
)$size

record_check(
  check = "Tracked file sizes",
  passed =
    length(tracked_sizes) > 0L &&
    all(
      tracked_sizes <
        10 * 1024^2,
      na.rm = TRUE
    ),
  details = paste(
    "largest tracked file:",
    round(
      max(
        tracked_sizes,
        na.rm = TRUE
      ) /
        1024^2,
      2
    ),
    "MiB"
  )
)

cat(
  "\nPublic repository check summary\n\n"
)

print(
  check_results,
  row.names = FALSE
)

failed_checks <- check_results[
  !check_results$passed,
  ,
  drop = FALSE
]

if (nrow(failed_checks) > 0L) {
  stop(
    nrow(failed_checks),
    " public repository check(s) failed.",
    call. = FALSE
  )
}

cat(
  "\nAll public repository checks passed.\n"
)
