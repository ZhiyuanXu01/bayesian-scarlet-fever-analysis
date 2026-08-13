# ==============================================================================
# Sensitivity analyses
#
# 1. PC-prior sensitivity for Bayesian conditional Poisson models at the
#    prespecified single-lag 3 and moving-average 4 windows.
# 2. Structural sensitivity specifications for pooled B-DLNMs.
# ==============================================================================

source(
  file.path(
    "R",
    "04_bayesian_conditional_poisson.R"
  )
)

source(
  file.path(
    "R",
    "05_pooled_bdlnm.R"
  )
)

SENSITIVITY_WINDOWS <- data.table::data.table(
  exposure_type = c(
    "lag",
    "mavg"
  ),
  window = c(
    3L,
    4L
  ),
  window_label = c(
    "Single lag 3",
    "Moving average 4"
  ),
  window_order = c(
    1L,
    2L
  )
)

LINEAR_PRIOR_SPECS <- data.table::data.table(
  prior_key = c(
    "pc_u2",
    "pc_u3_baseline",
    "pc_u4"
  ),
  prior_U = c(
    2,
    3,
    4
  ),
  prior_alpha = c(
    0.01,
    0.01,
    0.01
  ),
  prior_label = c(
    "PC prior U = 2",
    "PC prior U = 3 (baseline)",
    "PC prior U = 4"
  ),
  prior_order = c(
    1L,
    2L,
    3L
  )
)

BDLNM_SENSITIVITY_SPECS <- data.table::data.table(
  sensitivity_key = c(
    "lag5",
    "lag_df4",
    "knots_p33_p66"
  ),
  sensitivity_label = c(
    "Maximum lag = 5 days",
    "Lag spline df = 4",
    "Exposure knots = P33/P66"
  ),
  max_lag = c(
    5L,
    10L,
    10L
  ),
  lag_df = c(
    3L,
    4L,
    3L
  ),
  exposure_knot_probabilities = list(
    c(
      0.25,
      0.50,
      0.75
    ),
    c(
      0.25,
      0.50,
      0.75
    ),
    c(
      0.33,
      0.66
    )
  ),
  sensitivity_order = c(
    1L,
    2L,
    3L
  )
)

run_linear_prior_sensitivity_model <- function(
    variable,
    exposure_type = c(
      "lag",
      "mavg"
    ),
    window,
    prior_U,
    prior_alpha = 0.01,
    config_path = file.path(
      "config",
      "config.R"
    ),
    return_model = FALSE
) {
  requested_exposure_type <- match.arg(
    exposure_type
  )
  
  requested_window <- as.integer(
    window
  )
  
  requested_prior_U <- as.numeric(
    prior_U
  )
  
  requested_prior_alpha <- as.numeric(
    prior_alpha
  )
  
  if (length(requested_prior_U) != 1L ||
      !is.finite(requested_prior_U) ||
      requested_prior_U <= 0) {
    stop(
      "`prior_U` must be one positive finite number.",
      call. = FALSE
    )
  }
  
  if (length(requested_prior_alpha) != 1L ||
      !is.finite(requested_prior_alpha) ||
      requested_prior_alpha <= 0 ||
      requested_prior_alpha >= 1) {
    stop(
      "`prior_alpha` must be between zero and one.",
      call. = FALSE
    )
  }
  
  window_information <- SENSITIVITY_WINDOWS[
    SENSITIVITY_WINDOWS[["exposure_type"]] ==
      requested_exposure_type &
      SENSITIVITY_WINDOWS[["window"]] ==
      requested_window
  ]
  
  if (nrow(window_information) != 1L) {
    stop(
      "The requested sensitivity-analysis window is not prespecified.",
      call. = FALSE
    )
  }
  
  prior_information <- LINEAR_PRIOR_SPECS[
    LINEAR_PRIOR_SPECS[["prior_U"]] ==
      requested_prior_U &
      LINEAR_PRIOR_SPECS[["prior_alpha"]] ==
      requested_prior_alpha
  ]
  
  if (nrow(prior_information) != 1L) {
    stop(
      "The requested PC prior is not prespecified.",
      call. = FALSE
    )
  }
  
  config <- load_analysis_config(
    config_path
  )
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "linear"
  )
  
  fit <- fit_conditional_poisson(
    data = prepared$data,
    variable = variable,
    exposure_type =
      requested_exposure_type,
    window =
      requested_window,
    pc_prior = c(
      requested_prior_U,
      requested_prior_alpha
    ),
    inla_num_threads =
      config$computation$inla_num_threads,
    return_model =
      return_model
  )
  
  summary <- data.table::copy(
    fit$summary
  )
  
  summary[
    ,
    `:=`(
      prior_key =
        prior_information$prior_key,
      prior_label =
        prior_information$prior_label,
      prior_order =
        prior_information$prior_order,
      window_label =
        window_information$window_label,
      window_order =
        window_information$window_order
    )
  ]
  
  output <- list(
    summary = summary
  )
  
  if (isTRUE(return_model)) {
    output$model <- fit$model
  }
  
  output
}

build_linear_prior_sensitivity_specifications <- function() {
  specifications <- data.table::rbindlist(
    lapply(
      seq_len(nrow(LINEAR_PRIOR_SPECS)),
      function(prior_index) {
        current_prior <- LINEAR_PRIOR_SPECS[
          prior_index
        ]
        
        data.table::rbindlist(
          lapply(
            seq_len(nrow(SENSITIVITY_WINDOWS)),
            function(window_index) {
              current_window <- SENSITIVITY_WINDOWS[
                window_index
              ]
              
              data.table::data.table(
                prior_key =
                  current_prior$prior_key,
                prior_U =
                  current_prior$prior_U,
                prior_alpha =
                  current_prior$prior_alpha,
                prior_label =
                  current_prior$prior_label,
                prior_order =
                  current_prior$prior_order,
                exposure_type =
                  current_window$exposure_type,
                window =
                  current_window$window,
                window_label =
                  current_window$window_label,
                window_order =
                  current_window$window_order,
                variable =
                  ANALYSIS_EXPOSURES$variable,
                exposure =
                  ANALYSIS_EXPOSURES$label,
                exposure_order =
                  seq_len(
                    nrow(ANALYSIS_EXPOSURES)
                  )
              )
            }
          ),
          use.names = TRUE
        )
      }
    ),
    use.names = TRUE
  )
  
  data.table::setorder(
    specifications,
    prior_order,
    window_order,
    exposure_order
  )
  
  specifications[
    ,
    specification_order := .I
  ]
  
  specifications[
    ,
    model_key := paste(
      prior_key,
      variable,
      exposure_type,
      sprintf(
        "%02d",
        window
      ),
      sep = "__"
    )
  ]
  
  if (data.table::uniqueN(
    specifications$model_key
  ) != nrow(specifications)) {
    stop(
      "Linear sensitivity model keys are not unique.",
      call. = FALSE
    )
  }
  
  specifications[]
}

run_all_linear_prior_sensitivity_models <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    output_file = NULL,
    resume = TRUE
) {
  check_conditional_poisson_packages()
  
  config <- load_analysis_config(
    config_path
  )
  
  if (is.null(output_file)) {
    output_file <- file.path(
      config$output$results_dir,
      "sensitivity_linear_priors.csv"
    )
  }
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "linear"
  )
  
  specifications <-
    build_linear_prior_sensitivity_specifications()
  
  checkpoint <- if (isTRUE(resume)) {
    read_model_checkpoint(
      output_file
    )
  } else {
    data.table::data.table()
  }
  
  completed_keys <- if (
    nrow(checkpoint) > 0L &&
    all(
      c(
        "model_key",
        "status"
      ) %in%
      names(checkpoint)
    )
  ) {
    checkpoint[
      checkpoint[["status"]] == "ok",
      "model_key",
      with = FALSE
    ][["model_key"]]
  } else {
    character()
  }
  
  total_models <- nrow(
    specifications
  )
  
  set.seed(
    config$computation$random_seed
  )
  
  for (index in seq_len(total_models)) {
    specification <- specifications[
      index
    ]
    
    current_model_key <- as.character(
      specification$model_key
    )
    
    current_specification_order <- as.integer(
      specification$specification_order
    )
    
    if (current_model_key %in%
        completed_keys) {
      message(
        sprintf(
          "[%d/%d] Skipping completed model: %s",
          index,
          total_models,
          current_model_key
        )
      )
      
      next
    }
    
    message(
      sprintf(
        "[%d/%d] Fitting prior sensitivity: %s",
        index,
        total_models,
        current_model_key
      )
    )
    
    result_row <- tryCatch(
      {
        fit <- fit_conditional_poisson(
          data = prepared$data,
          variable =
            specification$variable,
          exposure_type =
            specification$exposure_type,
          window =
            specification$window,
          pc_prior = c(
            specification$prior_U,
            specification$prior_alpha
          ),
          inla_num_threads =
            config$computation$inla_num_threads,
          return_model = FALSE
        )
        
        result <- data.table::copy(
          fit$summary
        )
        
        result[
          ,
          `:=`(
            prior_key =
              as.character(
                specification$prior_key
              ),
            prior_label =
              as.character(
                specification$prior_label
              ),
            prior_order =
              as.integer(
                specification$prior_order
              ),
            window_label =
              as.character(
                specification$window_label
              ),
            window_order =
              as.integer(
                specification$window_order
              ),
            model_key =
              current_model_key,
            specification_order =
              current_specification_order,
            status = "ok",
            error_message =
              NA_character_
          )
        ]
        
        result
      },
      error = function(error) {
        data.table::data.table(
          prior_key =
            as.character(
              specification$prior_key
            ),
          prior_label =
            as.character(
              specification$prior_label
            ),
          prior_order =
            as.integer(
              specification$prior_order
            ),
          pc_prior_U =
            as.numeric(
              specification$prior_U
            ),
          pc_prior_alpha =
            as.numeric(
              specification$prior_alpha
            ),
          window_label =
            as.character(
              specification$window_label
            ),
          window_order =
            as.integer(
              specification$window_order
            ),
          variable =
            as.character(
              specification$variable
            ),
          exposure =
            as.character(
              specification$exposure
            ),
          exposure_type =
            as.character(
              specification$exposure_type
            ),
          window =
            as.integer(
              specification$window
            ),
          model_key =
            current_model_key,
          specification_order =
            current_specification_order,
          status = "failed",
          error_message =
            conditionMessage(error)
        )
      }
    )
    
    if (
      nrow(checkpoint) > 0L &&
      "model_key" %in%
      names(checkpoint)
    ) {
      checkpoint <- checkpoint[
        checkpoint[["model_key"]] !=
          current_model_key
      ]
    }
    
    checkpoint <- data.table::rbindlist(
      list(
        checkpoint,
        result_row
      ),
      use.names = TRUE,
      fill = TRUE
    )
    
    write_model_checkpoint(
      checkpoint,
      output_file
    )
    
    if (result_row$status == "ok") {
      message(
        sprintf(
          "Completed %s in %.1f seconds.",
          current_model_key,
          result_row$elapsed_seconds
        )
      )
    } else {
      message(
        "Failed ",
        current_model_key,
        ": ",
        result_row$error_message
      )
    }
    
    rm(
      result_row
    )
    
    invisible(
      gc()
    )
  }
  
  successful_models <- checkpoint[
    checkpoint[["status"]] == "ok"
  ]
  
  failed_models <- checkpoint[
    checkpoint[["status"]] != "ok"
  ]
  
  diagnostics <- checkpoint[
    ,
    .(
      models = .N,
      total_minutes = round(
        sum(
          elapsed_seconds,
          na.rm = TRUE
        ) /
          60,
        1
      )
    ),
    by = status
  ]
  
  message(
    sprintf(
      paste0(
        "Linear prior-sensitivity batch complete: ",
        "%d successful, %d failed."
      ),
      nrow(successful_models),
      nrow(failed_models)
    )
  )
  
  list(
    models = checkpoint,
    diagnostics = diagnostics
  )
}

compute_bdlnm_mean_slice <- function(
    fit_light,
    lag_window,
    probabilities = seq(
      0.01,
      0.99,
      by = 0.01
    ),
    reference_probability = 0.50
) {
  grid <- make_bdlnm_exposure_grid(
    fit_light$exposure_values,
    probabilities
  )
  
  reference_value <- as.numeric(
    stats::quantile(
      fit_light$exposure_values,
      probs = reference_probability,
      names = FALSE,
      na.rm = TRUE
    )
  )
  
  rows <- lapply(
    seq_len(nrow(grid)),
    function(grid_index) {
      contrast <- make_bdlnm_basis_contrast(
        fit_light = fit_light,
        exposure_value =
          grid$exposure_value[grid_index],
        reference_value =
          reference_value,
        lag_window =
          lag_window
      )
      
      data.table::data.table(
        variable =
          fit_light$variable,
        exposure_type =
          fit_light$exposure_type,
        lag_window =
          as.integer(lag_window),
        exposure_value =
          grid$exposure_value[grid_index],
        exposure_probability =
          grid$exposure_probability[grid_index],
        reference_value =
          reference_value,
        RR = exp(
          sum(
            fit_light$coefficient_mean *
              contrast
          )
        )
      )
    }
  )
  
  data.table::rbindlist(
    rows,
    use.names = TRUE
  )
}

build_bdlnm_sensitivity_specifications <- function() {
  specifications <- data.table::rbindlist(
    lapply(
      seq_len(nrow(BDLNM_SENSITIVITY_SPECS)),
      function(sensitivity_index) {
        current_sensitivity <-
          BDLNM_SENSITIVITY_SPECS[
            sensitivity_index
          ]
        
        data.table::rbindlist(
          lapply(
            seq_len(nrow(ANALYSIS_EXPOSURES)),
            function(exposure_index) {
              data.table::data.table(
                sensitivity_key =
                  current_sensitivity$sensitivity_key,
                sensitivity_label =
                  current_sensitivity$sensitivity_label,
                sensitivity_order =
                  current_sensitivity$sensitivity_order,
                max_lag =
                  current_sensitivity$max_lag,
                lag_df =
                  current_sensitivity$lag_df,
                exposure_knot_probabilities =
                  rep(
                    list(
                      current_sensitivity$
                        exposure_knot_probabilities[[
                          1L
                        ]]
                    ),
                    2L
                  ),
                variable =
                  ANALYSIS_EXPOSURES$variable[
                    exposure_index
                  ],
                exposure =
                  ANALYSIS_EXPOSURES$label[
                    exposure_index
                  ],
                exposure_order =
                  exposure_index,
                exposure_type = c(
                  "lag",
                  "mavg"
                ),
                exposure_type_order = c(
                  1L,
                  2L
                ),
                slice_window = c(
                  3L,
                  4L
                )
              )
            }
          ),
          use.names = TRUE
        )
      }
    ),
    use.names = TRUE
  )
  
  data.table::setorder(
    specifications,
    sensitivity_order,
    exposure_type_order,
    exposure_order
  )
  
  specifications[
    ,
    specification_order := .I
  ]
  
  specifications[
    ,
    model_key := paste(
      sensitivity_key,
      variable,
      exposure_type,
      sep = "__"
    )
  ]
  
  if (data.table::uniqueN(
    specifications$model_key
  ) != nrow(specifications)) {
    stop(
      "B-DLNM sensitivity model keys are not unique.",
      call. = FALSE
    )
  }
  
  specifications[]
}

run_all_bdlnm_sensitivity_models <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    metrics_file = file.path(
      "results",
      "sensitivity_bdlnm_metrics.csv"
    ),
    slices_file = file.path(
      "results",
      "sensitivity_bdlnm_slices.csv"
    ),
    resume = TRUE
) {
  check_bdlnm_packages()
  
  config <- load_analysis_config(
    config_path
  )
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "bdlnm"
  )
  
  specifications <-
    build_bdlnm_sensitivity_specifications()
  
  metrics <- if (isTRUE(resume)) {
    read_bdlnm_checkpoint(
      metrics_file
    )
  } else {
    data.table::data.table()
  }
  
  slices <- if (isTRUE(resume)) {
    read_bdlnm_checkpoint(
      slices_file
    )
  } else {
    data.table::data.table()
  }
  
  successful_metric_keys <- if (
    nrow(metrics) > 0L &&
    all(
      c(
        "model_key",
        "status"
      ) %in%
      names(metrics)
    )
  ) {
    metrics[
      metrics[["status"]] == "ok",
      "model_key",
      with = FALSE
    ][["model_key"]]
  } else {
    character()
  }
  
  completed_slice_keys <- if (
    nrow(slices) > 0L &&
    "model_key" %in%
    names(slices)
  ) {
    unique(
      slices$model_key
    )
  } else {
    character()
  }
  
  completed_keys <- intersect(
    successful_metric_keys,
    completed_slice_keys
  )
  
  total_models <- nrow(
    specifications
  )
  
  set.seed(
    config$computation$random_seed
  )
  
  for (index in seq_len(total_models)) {
    specification <- specifications[
      index
    ]
    
    current_model_key <- as.character(
      specification$model_key
    )
    
    current_specification_order <- as.integer(
      specification$specification_order
    )
    
    knot_probabilities <-
      specification$
      exposure_knot_probabilities[[
        1L
      ]]
    
    if (current_model_key %in%
        completed_keys) {
      message(
        sprintf(
          "[%d/%d] Skipping completed model: %s",
          index,
          total_models,
          current_model_key
        )
      )
      
      next
    }
    
    message(
      sprintf(
        "[%d/%d] Fitting B-DLNM sensitivity: %s",
        index,
        total_models,
        current_model_key
      )
    )
    
    outcome <- tryCatch(
      {
        fit <- fit_pooled_bdlnm(
          data = prepared$data,
          variable =
            specification$variable,
          exposure_type =
            specification$exposure_type,
          max_lag =
            specification$max_lag,
          lag_df =
            specification$lag_df,
          exposure_knot_probabilities =
            knot_probabilities,
          inla_num_threads =
            config$computation$inla_num_threads,
          compute_config = FALSE,
          compute_cpo = FALSE,
          return_model = FALSE
        )
        
        model_slice <- compute_bdlnm_mean_slice(
          fit_light =
            fit$fit_light,
          lag_window =
            specification$slice_window
        )
        
        model_slice[
          ,
          `:=`(
            sensitivity_key =
              as.character(
                specification$sensitivity_key
              ),
            sensitivity_label =
              as.character(
                specification$sensitivity_label
              ),
            sensitivity_order =
              as.integer(
                specification$sensitivity_order
              ),
            model_key =
              current_model_key,
            specification_order =
              current_specification_order
          )
        ]
        
        model_metrics <- data.table::copy(
          fit$summary
        )
        
        model_metrics[
          ,
          `:=`(
            sensitivity_key =
              as.character(
                specification$sensitivity_key
              ),
            sensitivity_label =
              as.character(
                specification$sensitivity_label
              ),
            sensitivity_order =
              as.integer(
                specification$sensitivity_order
              ),
            exposure_knot_probabilities =
              paste(
                knot_probabilities,
                collapse = ";"
              ),
            slice_window =
              as.integer(
                specification$slice_window
              ),
            model_key =
              current_model_key,
            specification_order =
              current_specification_order,
            status = "ok",
            error_message =
              NA_character_
          )
        ]
        
        list(
          metrics = model_metrics,
          slice = model_slice
        )
      },
      error = function(error) {
        list(
          metrics = data.table::data.table(
            sensitivity_key =
              as.character(
                specification$sensitivity_key
              ),
            sensitivity_label =
              as.character(
                specification$sensitivity_label
              ),
            sensitivity_order =
              as.integer(
                specification$sensitivity_order
              ),
            variable =
              as.character(
                specification$variable
              ),
            exposure =
              as.character(
                specification$exposure
              ),
            exposure_type =
              as.character(
                specification$exposure_type
              ),
            max_lag =
              as.integer(
                specification$max_lag
              ),
            lag_df =
              as.integer(
                specification$lag_df
              ),
            exposure_knot_probabilities =
              paste(
                knot_probabilities,
                collapse = ";"
              ),
            slice_window =
              as.integer(
                specification$slice_window
              ),
            model_key =
              current_model_key,
            specification_order =
              current_specification_order,
            status = "failed",
            error_message =
              conditionMessage(error)
          ),
          slice = data.table::data.table()
        )
      }
    )
    
    metrics <- replace_bdlnm_checkpoint_rows(
      existing = metrics,
      replacement =
        outcome$metrics,
      key_to_replace =
        current_model_key
    )
    
    if (
      nrow(slices) > 0L &&
      "model_key" %in%
      names(slices)
    ) {
      slices <- slices[
        slices[["model_key"]] !=
          current_model_key
      ]
    }
    
    if (nrow(outcome$slice) > 0L) {
      slices <- data.table::rbindlist(
        list(
          slices,
          outcome$slice
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }
    
    write_bdlnm_checkpoint(
      data = metrics,
      path = metrics_file,
      order_columns =
        "specification_order"
    )
    
    write_bdlnm_checkpoint(
      data = slices,
      path = slices_file,
      order_columns = c(
        "specification_order",
        "exposure_probability"
      )
    )
    
    if (outcome$metrics$status == "ok") {
      message(
        sprintf(
          "Completed %s in %.1f seconds.",
          current_model_key,
          outcome$metrics$elapsed_seconds
        )
      )
    } else {
      message(
        "Failed ",
        current_model_key,
        ": ",
        outcome$metrics$error_message
      )
    }
    
    rm(
      outcome
    )
    
    invisible(
      gc()
    )
  }
  
  successful_models <- metrics[
    metrics[["status"]] == "ok"
  ]
  
  failed_models <- metrics[
    metrics[["status"]] != "ok"
  ]
  
  diagnostics <- metrics[
    ,
    .(
      models = .N,
      total_minutes = round(
        sum(
          elapsed_seconds,
          na.rm = TRUE
        ) /
          60,
        1
      )
    ),
    by = status
  ]
  
  message(
    sprintf(
      paste0(
        "B-DLNM sensitivity batch complete: ",
        "%d successful, %d failed."
      ),
      nrow(successful_models),
      nrow(failed_models)
    )
  )
  
  list(
    metrics = metrics,
    slices = slices,
    diagnostics = diagnostics
  )
}

build_bdlnm_sensitivity_contrasts <- function(
    primary_surfaces,
    sensitivity_slices,
    target_probabilities = c(
      0.10,
      0.90
    )
) {
  primary <- data.table::copy(
    primary_surfaces[
      (
        primary_surfaces[["exposure_type"]] == "lag" &
          primary_surfaces[["lag_window"]] == 3L
      ) |
        (
          primary_surfaces[["exposure_type"]] == "mavg" &
            primary_surfaces[["lag_window"]] == 4L
        )
    ]
  )
  
  primary[
    ,
    `:=`(
      sensitivity_key = "baseline",
      sensitivity_label = "Primary specification",
      sensitivity_order = 0L
    )
  ]
  
  alternative <- data.table::copy(
    sensitivity_slices
  )
  
  combined <- data.table::rbindlist(
    list(
      primary[
        ,
        .(
          variable,
          exposure_type,
          lag_window,
          exposure_probability,
          RR,
          sensitivity_key,
          sensitivity_label,
          sensitivity_order
        )
      ],
      alternative[
        ,
        .(
          variable,
          exposure_type,
          lag_window,
          exposure_probability,
          RR,
          sensitivity_key,
          sensitivity_label,
          sensitivity_order
        )
      ]
    ),
    use.names = TRUE
  )
  
  contrasts <- combined[
    ,
    {
      ordered_rows <- order(
        exposure_probability
      )
      
      ordered_probability <-
        exposure_probability[
          ordered_rows
        ]
      
      ordered_rr <- RR[
        ordered_rows
      ]
      
      interpolated_rr <- stats::approx(
        x = ordered_probability,
        y = ordered_rr,
        xout = target_probabilities,
        method = "linear",
        ties = mean,
        rule = 2
      )$y
      
      .(
        exposure_probability =
          target_probabilities,
        RR = interpolated_rr
      )
    },
    by = .(
      variable,
      exposure_type,
      lag_window,
      sensitivity_key,
      sensitivity_label,
      sensitivity_order
    )
  ]
  
  contrasts[
    ,
    exposure := ANALYSIS_EXPOSURES$label[
      match(
        variable,
        ANALYSIS_EXPOSURES$variable
      )
    ]
  ]
  
  contrasts[
    ,
    exposure_order := match(
      variable,
      ANALYSIS_EXPOSURES$variable
    )
  ]
  
  contrasts[
    ,
    window_label := data.table::fcase(
      exposure_type == "lag" &
        lag_window == 3L,
      "Single lag 3",
      exposure_type == "mavg" &
        lag_window == 4L,
      "Moving average 4",
      default = NA_character_
    )
  ]
  
  contrasts[
    ,
    window_order := match(
      window_label,
      SENSITIVITY_WINDOWS$window_label
    )
  ]
  
  contrasts[
    ,
    contrast_label := paste0(
      "Exposure P",
      sprintf(
        "%02d",
        round(
          exposure_probability *
            100
        )
      ),
      " vs P50"
    )
  ]
  
  data.table::setorder(
    contrasts,
    window_order,
    exposure_probability,
    exposure_order,
    sensitivity_order
  )
  
  contrasts[]
}

run_bdlnm_sensitivity_contrasts <- function(
    primary_surface_file = file.path(
      "results",
      "pooled_bdlnm_surfaces.csv"
    ),
    sensitivity_slice_file = file.path(
      "results",
      "sensitivity_bdlnm_slices.csv"
    ),
    output_file = file.path(
      "results",
      "sensitivity_bdlnm_quantile_contrasts.csv"
    )
) {
  required_files <- c(
    primary_surface_file,
    sensitivity_slice_file
  )
  
  missing_files <- required_files[
    !file.exists(
      required_files
    )
  ]
  
  if (length(missing_files) > 0L) {
    stop(
      "Cannot find required result files: ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  
  contrasts <- build_bdlnm_sensitivity_contrasts(
    primary_surfaces =
      data.table::fread(
        primary_surface_file,
        showProgress = FALSE
      ),
    sensitivity_slices =
      data.table::fread(
        sensitivity_slice_file,
        showProgress = FALSE
      )
  )
  
  dir.create(
    dirname(
      output_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  data.table::fwrite(
    contrasts,
    output_file
  )
  
  message(
    "Wrote ",
    nrow(contrasts),
    " B-DLNM sensitivity contrasts."
  )
  
  contrasts
}

plot_bdlnm_sensitivity_contrasts <- function(
    contrasts,
    output_file = file.path(
      "figures",
      "sensitivity_bdlnm_quantile_contrasts.png"
    )
) {
  if (!requireNamespace(
    "ggplot2",
    quietly = TRUE
  )) {
    stop(
      "Package `ggplot2` is required.",
      call. = FALSE
    )
  }
  
  plot_data <- data.table::copy(
    contrasts
  )
  
  exposure_levels <-
    ANALYSIS_EXPOSURES$label
  
  sensitivity_levels <- c(
    "Primary specification",
    BDLNM_SENSITIVITY_SPECS$
      sensitivity_label
  )
  
  contrast_levels <- c(
    "Exposure P10 vs P50",
    "Exposure P90 vs P50"
  )
  
  plot_data[
    ,
    exposure := factor(
      exposure,
      levels = exposure_levels
    )
  ]
  
  plot_data[
    ,
    sensitivity_label := factor(
      sensitivity_label,
      levels = sensitivity_levels
    )
  ]
  
  plot_data[
    ,
    window_label := factor(
      window_label,
      levels =
        SENSITIVITY_WINDOWS$window_label
    )
  ]
  
  plot_data[
    ,
    contrast_label := factor(
      contrast_label,
      levels = contrast_levels
    )
  ]
  
  sensitivity_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = exposure,
      y = RR,
      colour = sensitivity_label,
      shape = sensitivity_label
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(
        width = 0.65
      ),
      size = 2.4
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(
        contrast_label
      ),
      cols = ggplot2::vars(
        window_label
      ),
      scales = "free_y"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Primary specification" = "#222222",
        "Maximum lag = 5 days" = "#3B6FB6",
        "Lag spline df = 4" = "#D17A22",
        "Exposure knots = P33/P66" = "#338A5A"
      ),
      name = "Model specification"
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Primary specification" = 16,
        "Maximum lag = 5 days" = 17,
        "Lag spline df = 4" = 15,
        "Exposure knots = P33/P66" = 18
      ),
      name = "Model specification"
    ) +
    ggplot2::labs(
      title = "Pooled B-DLNM structural sensitivity analysis",
      subtitle = paste0(
        "Posterior central relative-risk estimates under ",
        "alternative lag and spline specifications"
      ),
      x = "Environmental exposure",
      y = "Relative risk",
      caption = paste0(
        "Primary curves use posterior-mean surfaces. ",
        "Alternative specifications vary maximum lag, lag-spline df, ",
        "or exposure-knot placement."
      )
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid.minor =
        ggplot2::element_blank(),
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
      strip.text =
        ggplot2::element_text(
          face = "bold"
        ),
      plot.title =
        ggplot2::element_text(
          face = "bold"
        ),
      plot.caption =
        ggplot2::element_text(
          hjust = 0
        ),
      legend.position = "bottom"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      ),
      shape = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
    )
  
  dir.create(
    dirname(
      output_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggplot2::ggsave(
    filename = output_file,
    plot = sensitivity_plot,
    width = 14,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  message(
    "Wrote B-DLNM sensitivity figure to ",
    output_file,
    "."
  )
  
  invisible(
    sensitivity_plot
  )
}