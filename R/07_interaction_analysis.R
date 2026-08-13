# ==============================================================================
# Environmental-exposure interaction analysis
#
# Exposure pairs are screened using Spearman correlations before interaction
# models are fitted. The two prespecified windows are single lag 3 and moving
# average 4.
# ==============================================================================

source(
  file.path(
    "R",
    "04_bayesian_conditional_poisson.R"
  )
)

INTERACTION_WINDOWS <- data.table::data.table(
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
  )
)

build_interaction_term_table <- function(
    exposure_type,
    window
) {
  term_table <- data.table::data.table(
    variable = ANALYSIS_EXPOSURES$variable,
    exposure = ANALYSIS_EXPOSURES$label
  )
  
  term_table[
    ,
    term := vapply(
      variable,
      function(current_variable) {
        get_exposure_term(
          variable = current_variable,
          exposure_type = exposure_type,
          window = window
        )
      },
      FUN.VALUE = character(1)
    )
  ]
  
  term_table
}

compute_spearman_screen <- function(
    data,
    exposure_type,
    window,
    window_label,
    correlation_cutoff = 0.70
) {
  if (length(correlation_cutoff) != 1L ||
      !is.finite(correlation_cutoff) ||
      correlation_cutoff <= 0 ||
      correlation_cutoff >= 1) {
    stop(
      "`correlation_cutoff` must be between zero and one.",
      call. = FALSE
    )
  }
  
  term_table <- build_interaction_term_table(
    exposure_type = exposure_type,
    window = window
  )
  
  missing_columns <- setdiff(
    term_table$term,
    names(data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing interaction-screening columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  pair_indices <- utils::combn(
    seq_len(nrow(term_table)),
    2L
  )
  
  screening_rows <- lapply(
    seq_len(ncol(pair_indices)),
    function(pair_number) {
      first_index <- pair_indices[
        1L,
        pair_number
      ]
      
      second_index <- pair_indices[
        2L,
        pair_number
      ]
      
      first_exposure <- term_table[
        first_index
      ]
      
      second_exposure <- term_table[
        second_index
      ]
      
      first_values <- data[[
        first_exposure$term
      ]]
      
      second_values <- data[[
        second_exposure$term
      ]]
      
      complete_rows <- is.finite(
        first_values
      ) &
        is.finite(
          second_values
        )
      
      n_complete <- sum(
        complete_rows
      )
      
      spearman_rho <- if (n_complete < 10L) {
        NA_real_
      } else {
        stats::cor(
          first_values[
            complete_rows
          ],
          second_values[
            complete_rows
          ],
          method = "spearman"
        )
      }
      
      included <- is.finite(
        spearman_rho
      ) &&
        abs(
          spearman_rho
        ) <= correlation_cutoff
      
      data.table::data.table(
        exposure_type = exposure_type,
        window = as.integer(
          window
        ),
        window_label = window_label,
        variable_1 = first_exposure$variable,
        exposure_1 = first_exposure$exposure,
        term_1 = first_exposure$term,
        variable_2 = second_exposure$variable,
        exposure_2 = second_exposure$exposure,
        term_2 = second_exposure$term,
        spearman_rho = spearman_rho,
        absolute_rho = abs(
          spearman_rho
        ),
        n_complete = n_complete,
        correlation_cutoff = correlation_cutoff,
        decision = if (included) {
          "included"
        } else if (is.finite(spearman_rho)) {
          "excluded_high_correlation"
        } else {
          "excluded_unavailable"
        }
      )
    }
  )
  
  data.table::rbindlist(
    screening_rows,
    use.names = TRUE
  )
}

run_interaction_screen <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    correlation_cutoff = 0.70,
    output_path = file.path(
      "results",
      "interaction_spearman_screen.csv"
    )
) {
  if (!requireNamespace(
    "data.table",
    quietly = TRUE
  )) {
    stop(
      "Package `data.table` is required.",
      call. = FALSE
    )
  }
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "linear"
  )
  
  screening_results <- lapply(
    seq_len(nrow(INTERACTION_WINDOWS)),
    function(window_index) {
      current_window <- INTERACTION_WINDOWS[
        window_index
      ]
      
      compute_spearman_screen(
        data = prepared$data,
        exposure_type =
          current_window$exposure_type,
        window = current_window$window,
        window_label =
          current_window$window_label,
        correlation_cutoff =
          correlation_cutoff
      )
    }
  )
  
  screening_results <- data.table::rbindlist(
    screening_results,
    use.names = TRUE
  )
  
  diagnostics <- screening_results[
    ,
    .(
      exposure_pairs = .N,
      included_pairs = sum(
        decision == "included"
      ),
      excluded_pairs = sum(
        decision != "included"
      ),
      minimum_rho = min(
        spearman_rho,
        na.rm = TRUE
      ),
      maximum_rho = max(
        spearman_rho,
        na.rm = TRUE
      ),
      maximum_absolute_rho = max(
        absolute_rho,
        na.rm = TRUE
      )
    ),
    by = .(
      exposure_type,
      window,
      window_label
    )
  ]
  
  dir.create(
    dirname(
      output_path
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  data.table::fwrite(
    screening_results,
    output_path
  )
  
  message(
    "Interaction screen completed: ",
    nrow(screening_results),
    " exposure pairs."
  )
  
  list(
    screen = screening_results,
    diagnostics = diagnostics
  )
}

fit_iqr_scaled_interaction <- function(
    data,
    screening_row,
    pc_prior = c(
      3,
      0.01
    ),
    inla_num_threads = "1:1",
    return_model = FALSE
) {
  check_conditional_poisson_packages()
  
  if (nrow(screening_row) != 1L) {
    stop(
      "`screening_row` must contain exactly one exposure pair.",
      call. = FALSE
    )
  }
  
  if (screening_row$decision != "included") {
    stop(
      "The selected pair was excluded during correlation screening.",
      call. = FALSE
    )
  }
  
  first_term <- screening_row$term_1
  second_term <- screening_row$term_2
  
  required_columns <- c(
    "case",
    "global_id",
    first_term,
    second_term,
    "holiday",
    "outbreak_effect"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing interaction-model columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  model_data <- data.table::copy(
    data[
      stats::complete.cases(
        data[
          ,
          ..required_columns
        ]
      ),
      ..required_columns
    ]
  )
  
  first_centre <- stats::median(
    model_data[[
      first_term
    ]]
  )
  
  second_centre <- stats::median(
    model_data[[
      second_term
    ]]
  )
  
  first_iqr <- stats::IQR(
    model_data[[
      first_term
    ]]
  )
  
  second_iqr <- stats::IQR(
    model_data[[
      second_term
    ]]
  )
  
  if (!is.finite(first_iqr) ||
      first_iqr <= 0 ||
      !is.finite(second_iqr) ||
      second_iqr <= 0) {
    stop(
      "Both exposures must have a positive finite IQR.",
      call. = FALSE
    )
  }
  
  model_data[
    ,
    exposure_1_iqr :=
      (
        get(first_term) -
          first_centre
      ) /
      first_iqr
  ]
  
  model_data[
    ,
    exposure_2_iqr :=
      (
        get(second_term) -
          second_centre
      ) /
      second_iqr
  ]
  
  model_data[
    ,
    gid_inla := as.integer(
      factor(global_id)
    )
  ]
  
  hyper_stratum <- list(
    prec = list(
      prior = "pc.prec",
      param = pc_prior
    )
  )
  
  model_formula <- stats::as.formula(
    paste0(
      "case ~ 1 + ",
      "exposure_1_iqr * exposure_2_iqr",
      " + holiday + outbreak_effect",
      " + f(",
      "gid_inla, ",
      "model = 'iid', ",
      "hyper = hyper_stratum",
      ")"
    )
  )
  
  environment(
    model_formula
  ) <- environment()
  
  started_at <- proc.time()[[
    "elapsed"
  ]]
  
  model <- INLA::inla(
    model_formula,
    family = "poisson",
    data = model_data,
    control.fixed = list(
      mean.intercept = 0,
      prec.intercept = 1,
      mean = 0,
      prec = 1
    ),
    control.compute = list(
      config = FALSE,
      dic = TRUE,
      waic = TRUE,
      cpo = FALSE
    ),
    control.predictor = list(
      compute = FALSE
    ),
    control.inla = list(
      strategy = "simplified.laplace"
    ),
    num.threads = inla_num_threads,
    verbose = FALSE
  )
  
  elapsed_seconds <- unname(
    proc.time()[["elapsed"]] -
      started_at
  )
  
  fixed_summary <- model$summary.fixed
  
  interaction_term <-
    "exposure_1_iqr:exposure_2_iqr"
  
  expected_terms <- c(
    "exposure_1_iqr",
    "exposure_2_iqr",
    interaction_term
  )
  
  missing_terms <- setdiff(
    expected_terms,
    rownames(fixed_summary)
  )
  
  if (length(missing_terms) > 0L) {
    stop(
      "INLA did not return expected terms: ",
      paste(missing_terms, collapse = ", "),
      call. = FALSE
    )
  }
  
  centre_column <- if (
    "0.5quant" %in%
    colnames(fixed_summary)
  ) {
    "0.5quant"
  } else {
    "mean"
  }
  
  first_coefficient <- as.numeric(
    fixed_summary[
      "exposure_1_iqr",
      centre_column
    ]
  )
  
  second_coefficient <- as.numeric(
    fixed_summary[
      "exposure_2_iqr",
      centre_column
    ]
  )
  
  interaction_coefficient <- as.numeric(
    fixed_summary[
      interaction_term,
      centre_column
    ]
  )
  
  result <- data.table::data.table(
    exposure_type =
      screening_row$exposure_type,
    window =
      screening_row$window,
    window_label =
      screening_row$window_label,
    variable_1 =
      screening_row$variable_1,
    exposure_1 =
      screening_row$exposure_1,
    term_1 =
      screening_row$term_1,
    variable_2 =
      screening_row$variable_2,
    exposure_2 =
      screening_row$exposure_2,
    term_2 =
      screening_row$term_2,
    spearman_rho =
      screening_row$spearman_rho,
    exposure_1_median =
      first_centre,
    exposure_1_IQR =
      first_iqr,
    exposure_2_median =
      second_centre,
    exposure_2_IQR =
      second_iqr,
    exposure_1_coefficient =
      first_coefficient,
    exposure_1_RR =
      exp(first_coefficient),
    exposure_1_lower_95_CrI =
      exp(
        fixed_summary[
          "exposure_1_iqr",
          "0.025quant"
        ]
      ),
    exposure_1_upper_95_CrI =
      exp(
        fixed_summary[
          "exposure_1_iqr",
          "0.975quant"
        ]
      ),
    exposure_2_coefficient =
      second_coefficient,
    exposure_2_RR =
      exp(second_coefficient),
    exposure_2_lower_95_CrI =
      exp(
        fixed_summary[
          "exposure_2_iqr",
          "0.025quant"
        ]
      ),
    exposure_2_upper_95_CrI =
      exp(
        fixed_summary[
          "exposure_2_iqr",
          "0.975quant"
        ]
      ),
    interaction_coefficient =
      interaction_coefficient,
    interaction_coefficient_sd =
      fixed_summary[
        interaction_term,
        "sd"
      ],
    interaction_ratio =
      exp(interaction_coefficient),
    interaction_lower_95_CrI =
      exp(
        fixed_summary[
          interaction_term,
          "0.025quant"
        ]
      ),
    interaction_upper_95_CrI =
      exp(
        fixed_summary[
          interaction_term,
          "0.975quant"
        ]
      ),
    n_rows = nrow(
      model_data
    ),
    n_strata = data.table::uniqueN(
      model_data$gid_inla
    ),
    DIC = extract_model_metric(
      model,
      "dic",
      "dic"
    ),
    WAIC = extract_model_metric(
      model,
      "waic",
      "waic"
    ),
    pc_prior_U = pc_prior[1],
    pc_prior_alpha = pc_prior[2],
    elapsed_seconds =
      elapsed_seconds
  )
  
  output <- list(
    summary = result
  )
  
  if (isTRUE(return_model)) {
    output$model <- model
  }
  
  output
}

run_interaction_model <- function(
    variable_1,
    variable_2,
    exposure_type = c(
      "lag",
      "mavg"
    ),
    window,
    config_path = file.path(
      "config",
      "config.R"
    ),
    correlation_cutoff = 0.70,
    return_model = FALSE
) {
  requested_exposure_type <- match.arg(
    exposure_type
  )
  
  requested_window <- as.integer(
    window
  )
  
  requested_variable_1 <- as.character(
    variable_1
  )
  
  requested_variable_2 <- as.character(
    variable_2
  )
  
  config <- load_analysis_config(
    config_path
  )
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "linear"
  )
  
  window_row <- INTERACTION_WINDOWS[
    INTERACTION_WINDOWS$exposure_type ==
      requested_exposure_type &
      INTERACTION_WINDOWS$window ==
      requested_window
  ]
  
  if (nrow(window_row) != 1L) {
    stop(
      "The requested interaction window is not prespecified.",
      call. = FALSE
    )
  }
  
  screen <- compute_spearman_screen(
    data = prepared$data,
    exposure_type =
      requested_exposure_type,
    window =
      requested_window,
    window_label =
      window_row$window_label,
    correlation_cutoff =
      correlation_cutoff
  )
  
  selected_pair_index <- (
    screen[["variable_1"]] ==
      requested_variable_1 &
      screen[["variable_2"]] ==
      requested_variable_2
  ) |
    (
      screen[["variable_1"]] ==
        requested_variable_2 &
        screen[["variable_2"]] ==
        requested_variable_1
    )
  
  selected_pair <- screen[
    selected_pair_index
  ]
  
  if (nrow(selected_pair) != 1L) {
    stop(
      "Could not identify exactly one requested exposure pair.",
      call. = FALSE
    )
  }
  
  fit_iqr_scaled_interaction(
    data = prepared$data,
    screening_row = selected_pair,
    inla_num_threads =
      config$computation$inla_num_threads,
    return_model = return_model
  )
}

build_interaction_specifications <- function(
    screening_results
) {
  specifications <- data.table::copy(
    screening_results[
      screening_results$decision ==
        "included"
    ]
  )
  
  data.table::setorder(
    specifications,
    exposure_type,
    window,
    variable_1,
    variable_2
  )
  
  specifications[
    ,
    specification_order := .I
  ]
  
  specifications[
    ,
    model_key := paste(
      exposure_type,
      sprintf(
        "%02d",
        window
      ),
      variable_1,
      variable_2,
      sep = "__"
    )
  ]
  
  specifications[]
}

run_all_interaction_models <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    screening_file = NULL,
    output_file = NULL,
    correlation_cutoff = 0.70,
    resume = TRUE
) {
  check_conditional_poisson_packages()
  
  config <- load_analysis_config(
    config_path
  )
  
  if (is.null(screening_file)) {
    screening_file <- file.path(
      config$output$results_dir,
      "interaction_spearman_screen.csv"
    )
  }
  
  if (is.null(output_file)) {
    output_file <- file.path(
      config$output$results_dir,
      "interaction_models.csv"
    )
  }
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "linear"
  )
  
  screening_results <- lapply(
    seq_len(nrow(INTERACTION_WINDOWS)),
    function(window_index) {
      current_window <- INTERACTION_WINDOWS[
        window_index
      ]
      
      compute_spearman_screen(
        data = prepared$data,
        exposure_type =
          current_window$exposure_type,
        window =
          current_window$window,
        window_label =
          current_window$window_label,
        correlation_cutoff =
          correlation_cutoff
      )
    }
  )
  
  screening_results <- data.table::rbindlist(
    screening_results,
    use.names = TRUE
  )
  
  dir.create(
    dirname(
      screening_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  data.table::fwrite(
    screening_results,
    screening_file
  )
  
  specifications <-
    build_interaction_specifications(
      screening_results
    )
  
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
    
    current_model_key <-
      specification$model_key
    
    current_specification_order <-
      specification$specification_order
    
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
        "[%d/%d] Fitting interaction: %s",
        index,
        total_models,
        current_model_key
      )
    )
    
    result_row <- tryCatch(
      {
        fit <- fit_iqr_scaled_interaction(
          data = prepared$data,
          screening_row = specification,
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
          exposure_type =
            specification$exposure_type,
          window =
            specification$window,
          window_label =
            specification$window_label,
          variable_1 =
            specification$variable_1,
          exposure_1 =
            specification$exposure_1,
          term_1 =
            specification$term_1,
          variable_2 =
            specification$variable_2,
          exposure_2 =
            specification$exposure_2,
          term_2 =
            specification$term_2,
          spearman_rho =
            specification$spearman_rho,
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
        "Interaction batch complete: ",
        "%d successful, %d failed."
      ),
      nrow(successful_models),
      nrow(failed_models)
    )
  )
  
  list(
    screen = screening_results,
    models = checkpoint,
    diagnostics = diagnostics
  )
}

plot_interaction_heatmap <- function(
    screening_results,
    model_results,
    output_file = file.path(
      "figures",
      "interaction_ratio_heatmap.png"
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
  
  screening_plot <- data.table::copy(
    screening_results
  )
  
  model_plot <- data.table::copy(
    model_results[
      model_results[["status"]] == "ok"
    ]
  )
  
  join_columns <- c(
    "exposure_type",
    "window",
    "variable_1",
    "variable_2"
  )
  
  plot_data <- merge(
    screening_plot,
    model_plot[
      ,
      c(
        join_columns,
        "interaction_ratio",
        "interaction_lower_95_CrI",
        "interaction_upper_95_CrI"
      ),
      with = FALSE
    ],
    by = join_columns,
    all.x = TRUE,
    sort = FALSE
  )
  
  short_labels <- c(
    t2m = "Temperature",
    rh = "Humidity",
    ssr = "Solar radiation",
    uv = "Wind speed",
    sp = "Pressure",
    tp = "Precipitation",
    PM25 = "PM2.5",
    PM10 = "PM10",
    O3 = "O3"
  )
  
  exposure_levels <- unname(
    short_labels[
      ANALYSIS_EXPOSURES$variable
    ]
  )
  
  plot_data[
    ,
    exposure_1_short :=
      unname(
        short_labels[
          variable_1
        ]
      )
  ]
  
  plot_data[
    ,
    exposure_2_short :=
      unname(
        short_labels[
          variable_2
        ]
      )
  ]
  
  plot_data[
    ,
    interval_excludes_one :=
      !is.na(
        interaction_lower_95_CrI
      ) &
      (
        interaction_lower_95_CrI > 1 |
          interaction_upper_95_CrI < 1
      )
  ]
  
  plot_data[
    ,
    cell_label := data.table::fcase(
      decision != "included",
      "\u00d7",
      interval_excludes_one,
      sprintf(
        "%.3f*",
        interaction_ratio
      ),
      default = sprintf(
        "%.3f",
        interaction_ratio
      )
    )
  ]
  
  plot_data[
    ,
    exposure_1_short := factor(
      exposure_1_short,
      levels = exposure_levels
    )
  ]
  
  plot_data[
    ,
    exposure_2_short := factor(
      exposure_2_short,
      levels = rev(
        exposure_levels
      )
    )
  ]
  
  plot_data[
    ,
    window_label := factor(
      window_label,
      levels =
        INTERACTION_WINDOWS$window_label
    )
  ]
  
  maximum_deviation <- max(
    abs(
      plot_data$interaction_ratio -
        1
    ),
    na.rm = TRUE
  )
  
  fill_limits <- c(
    1 - maximum_deviation,
    1 + maximum_deviation
  )
  
  interaction_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = exposure_1_short,
      y = exposure_2_short,
      fill = interaction_ratio
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.6
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = cell_label
      ),
      size = 3.2
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(
        window_label
      ),
      ncol = 2
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#3B6FB6",
      mid = "white",
      high = "#C7473A",
      midpoint = 1,
      limits = fill_limits,
      na.value = "grey85",
      name = "Interaction\nratio"
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "IQR-scaled environmental exposure interactions",
      subtitle = paste0(
        "Posterior interaction ratios; ",
        "\u00d7 denotes pairs excluded when ",
        "|Spearman rho| > 0.70"
      ),
      x = "First exposure",
      y = "Second exposure",
      caption = paste0(
        "Ratios represent multiplicative interaction terms ",
        "for joint one-IQR increases. ",
        "* indicates a 95% credible interval excluding 1."
      )
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0
      ),
      legend.position = "right"
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
    plot = interaction_plot,
    width = 13,
    height = 7.5,
    dpi = 300,
    bg = "white"
  )
  
  message(
    "Wrote interaction heatmap to ",
    output_file,
    "."
  )
  
  invisible(
    interaction_plot
  )
}