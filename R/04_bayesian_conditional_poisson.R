# ==============================================================================
# Bayesian conditional Poisson models
#
# Each model includes one environmental exposure, holiday and outbreak-effect
# covariates, and an iid effect for the matched case-control stratum.
# ==============================================================================

source(
  file.path(
    "R",
    "02_prepare_analysis_data.R"
  )
)

check_conditional_poisson_packages <- function() {
  required_packages <- c(
    "data.table",
    "INLA"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

get_exposure_term <- function(
  variable,
  exposure_type = c(
    "lag",
    "mavg"
  ),
  window
) {
  exposure_type <- match.arg(
    exposure_type
  )

  if (!variable %in%
      ANALYSIS_EXPOSURES$variable) {
    stop(
      "Unknown exposure variable: ",
      variable,
      call. = FALSE
    )
  }

  if (length(window) != 1L ||
      is.na(window) ||
      window != as.integer(window)) {
    stop(
      "`window` must be one integer.",
      call. = FALSE
    )
  }

  window <- as.integer(window)

  if (exposure_type == "lag") {
    if (window < 0L ||
        window > ANALYSIS_MAX_LAG) {
      stop(
        "Single-day lag must be between 0 and ",
        ANALYSIS_MAX_LAG,
        ".",
        call. = FALSE
      )
    }

    if (window == 0L) {
      variable
    } else {
      paste0(
        variable,
        "_lag",
        window
      )
    }
  } else {
    if (window < 1L ||
        window > ANALYSIS_MAX_LAG) {
      stop(
        "Moving-average window must be between 1 and ",
        ANALYSIS_MAX_LAG,
        ".",
        call. = FALSE
      )
    }

    paste0(
      variable,
      "_mavg",
      window
    )
  }
}

extract_model_metric <- function(
  model,
  component,
  metric
) {
  value <- tryCatch(
    model[[component]][[metric]],
    error = function(error) {
      NULL
    }
  )

  if (is.null(value) ||
      length(value) != 1L) {
    return(NA_real_)
  }

  as.numeric(value)
}

fit_conditional_poisson <- function(
  data,
  variable,
  exposure_type = c(
    "lag",
    "mavg"
  ),
  window,
  pc_prior = c(
    3,
    0.01
  ),
  inla_num_threads = "1:1",
  return_model = FALSE
) {
  check_conditional_poisson_packages()

  exposure_type <- match.arg(
    exposure_type
  )

  term <- get_exposure_term(
    variable = variable,
    exposure_type = exposure_type,
    window = window
  )

  required_columns <- c(
    "case",
    "global_id",
    term,
    "holiday",
    "outbreak_effect"
  )

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing model columns: ",
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
      term,
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

  if (!term %in%
      rownames(fixed_summary)) {
    stop(
      "Exposure term was not returned by INLA: ",
      term,
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

  coefficient_centre <- fixed_summary[
    term,
    centre_column
  ]

  result <- data.table::data.table(
    variable = variable,
    exposure = ANALYSIS_EXPOSURES$label[
      match(
        variable,
        ANALYSIS_EXPOSURES$variable
      )
    ],
    exposure_type = exposure_type,
    window = as.integer(window),
    term = term,
    coefficient = as.numeric(
      coefficient_centre
    ),
    coefficient_sd = as.numeric(
      fixed_summary[
        term,
        "sd"
      ]
    ),
    RR = exp(
      coefficient_centre
    ),
    lower_95_CrI = exp(
      fixed_summary[
        term,
        "0.025quant"
      ]
    ),
    upper_95_CrI = exp(
      fixed_summary[
        term,
        "0.975quant"
      ]
    ),
    n_rows = nrow(model_data),
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
    elapsed_seconds = elapsed_seconds
  )

  output <- list(
    summary = result
  )

  if (isTRUE(return_model)) {
    output$model <- model
  }

  output
}

run_conditional_poisson_model <- function(
  variable,
  exposure_type = c(
    "lag",
    "mavg"
  ),
  window,
  config_path = file.path(
    "config",
    "config.R"
  ),
  return_model = FALSE
) {
  exposure_type <- match.arg(
    exposure_type
  )

  config <- load_analysis_config(
    config_path
  )

  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "linear"
  )

  set.seed(
    config$computation$random_seed
  )

  fit_conditional_poisson(
    data = prepared$data,
    variable = variable,
    exposure_type = exposure_type,
    window = window,
    inla_num_threads =
      config$computation$inla_num_threads,
    return_model = return_model
  )
}

build_conditional_poisson_specifications <- function(
    max_lag = ANALYSIS_MAX_LAG
) {
  specifications <- lapply(
    ANALYSIS_EXPOSURES$variable,
    function(variable) {
      data.table::rbindlist(
        list(
          data.table::data.table(
            variable = variable,
            exposure_type = "lag",
            window = 0L:max_lag
          ),
          data.table::data.table(
            variable = variable,
            exposure_type = "mavg",
            window = seq_len(max_lag)
          )
        )
      )
    }
  )
  
  specifications <- data.table::rbindlist(
    specifications
  )
  
  specifications[
    ,
    specification_order := .I
  ]
  
  specifications[
    ,
    model_key := sprintf(
      "%s__%s__%02d",
      variable,
      exposure_type,
      window
    )
  ]
  
  specifications[]
}

read_model_checkpoint <- function(
    output_file
) {
  if (!file.exists(output_file)) {
    return(
      data.table::data.table()
    )
  }
  
  checkpoint <- data.table::fread(
    output_file,
    na.strings = c(
      "",
      "NA"
    ),
    showProgress = FALSE
  )
  
  if (!"model_key" %in%
      names(checkpoint)) {
    stop(
      "Existing checkpoint lacks `model_key`: ",
      output_file,
      call. = FALSE
    )
  }
  
  checkpoint
}

write_model_checkpoint <- function(
    checkpoint,
    output_file
) {
  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if ("specification_order" %in%
      names(checkpoint)) {
    data.table::setorder(
      checkpoint,
      specification_order
    )
  }
  
  data.table::fwrite(
    checkpoint,
    output_file,
    na = ""
  )
  
  invisible(output_file)
}

run_all_conditional_poisson_models <- function(
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
      "conditional_poisson_models.csv"
    )
  }
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "linear"
  )
  
  specifications <-
    build_conditional_poisson_specifications()
  
  checkpoint <- if (isTRUE(resume)) {
    read_model_checkpoint(
      output_file
    )
  } else {
    data.table::data.table()
  }
  
  completed_keys <- if (
    nrow(checkpoint) > 0L &&
    "status" %in% names(checkpoint)
  ) {
    checkpoint[
      status == "ok",
      model_key
    ]
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
    
    model_key <- specification$model_key
    
    if (model_key %in%
        completed_keys) {
      message(
        sprintf(
          "[%d/%d] Skipping completed model: %s",
          index,
          total_models,
          model_key
        )
      )
      
      next
    }
    
    message(
      sprintf(
        "[%d/%d] Fitting: %s",
        index,
        total_models,
        model_key
      )
    )
    
    result_row <- tryCatch(
      {
        fit <- fit_conditional_poisson(
          data = prepared$data,
          variable = specification$variable,
          exposure_type =
            specification$exposure_type,
          window = specification$window,
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
            model_key = model_key,
            specification_order =
              specification$specification_order,
            status = "ok",
            error_message = NA_character_
          )
        ]
        
        result
      },
      error = function(error) {
        data.table::data.table(
          model_key = model_key,
          specification_order =
            specification$specification_order,
          variable = specification$variable,
          exposure = ANALYSIS_EXPOSURES$label[
            match(
              specification$variable,
              ANALYSIS_EXPOSURES$variable
            )
          ],
          exposure_type =
            specification$exposure_type,
          window = specification$window,
          status = "failed",
          error_message = conditionMessage(
            error
          )
        )
      }
    )
    
    if (nrow(checkpoint) > 0L) {
      checkpoint <- checkpoint[
        model_key !=
          result_row$model_key[1L]
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
          model_key,
          result_row$elapsed_seconds
        )
      )
    } else {
      message(
        "Failed ",
        model_key,
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
    status == "ok"
  ]
  
  failed_models <- checkpoint[
    status != "ok"
  ]
  
  message(
    "Batch complete: ",
    nrow(successful_models),
    " successful, ",
    nrow(failed_models),
    " failed."
  )
  
  checkpoint[]
}

check_conditional_poisson_plot_packages <- function() {
  required_packages <- c(
    "ggplot2",
    "scales"
  )
  
  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  
  if (length(missing_packages) > 0L) {
    stop(
      "Missing plotting packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

plot_conditional_poisson_results <- function(
    results,
    exposure_type = c(
      "lag",
      "mavg"
    )
) {
  check_conditional_poisson_plot_packages()
  
  exposure_type <- match.arg(
    exposure_type
  )
  
  plot_data <- data.table::copy(
    results[
      status == "ok" &
        get("exposure_type") ==
        exposure_type
    ]
  )
  
  if (nrow(plot_data) == 0L) {
    stop(
      "No successful models found for: ",
      exposure_type,
      call. = FALSE
    )
  }
  
  plot_data[
    ,
    exposure := factor(
      exposure,
      levels =
        ANALYSIS_EXPOSURES$label
    )
  ]
  
  data.table::setorder(
    plot_data,
    exposure,
    window
  )
  
  plot_settings <- if (
    exposure_type == "lag"
  ) {
    list(
      title =
        "Single-day environmental associations",
      subtitle = paste(
        "Bayesian conditional Poisson models,",
        "lag 0 through lag 10"
      ),
      x_label = "Single-day lag",
      colour = "#397D8C",
      breaks = 0L:ANALYSIS_MAX_LAG
    )
  } else {
    list(
      title =
        "Moving-average environmental associations",
      subtitle = paste(
        "Bayesian conditional Poisson models,",
        "cumulative windows from lag 0 through lag k"
      ),
      x_label = "Moving-average window (k)",
      colour = "#C8553D",
      breaks = seq_len(
        ANALYSIS_MAX_LAG
      )
    )
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = window,
      y = RR
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      colour = "grey45",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower_95_CrI,
        ymax = upper_95_CrI
      ),
      colour = plot_settings$colour,
      linewidth = 0.45,
      width = 0.14,
      alpha = 0.75
    ) +
    ggplot2::geom_line(
      colour = plot_settings$colour,
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      colour = plot_settings$colour,
      fill = "white",
      shape = 21,
      stroke = 0.8,
      size = 2
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(
        exposure
      ),
      ncol = 3,
      scales = "free_y"
    ) +
    ggplot2::scale_x_continuous(
      breaks = plot_settings$breaks
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.0001
      )
    ) +
    ggplot2::labs(
      title = plot_settings$title,
      subtitle = plot_settings$subtitle,
      x = plot_settings$x_label,
      y = "Relative risk",
      caption = paste(
        "Points show posterior median RRs per one-unit",
        "increase; error bars show 95% credible intervals.",
        "Y-axis scales vary across exposures."
      )
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid.minor =
        ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.caption = ggplot2::element_text(
        colour = "grey35"
      )
    )
}

save_conditional_poisson_figures <- function(
    results_file = file.path(
      "results",
      "conditional_poisson_models.csv"
    ),
    figures_dir = "figures"
) {
  if (!file.exists(results_file)) {
    stop(
      "Cannot find model results: ",
      results_file,
      call. = FALSE
    )
  }
  
  results <- data.table::fread(
    results_file,
    showProgress = FALSE
  )
  
  dir.create(
    figures_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  lag_plot <- plot_conditional_poisson_results(
    results,
    exposure_type = "lag"
  )
  
  moving_average_plot <-
    plot_conditional_poisson_results(
      results,
      exposure_type = "mavg"
    )
  
  lag_file <- file.path(
    figures_dir,
    "conditional_poisson_single_lag_rr.png"
  )
  
  moving_average_file <- file.path(
    figures_dir,
    "conditional_poisson_moving_average_rr.png"
  )
  
  ggplot2::ggsave(
    filename = lag_file,
    plot = lag_plot,
    width = 10,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = moving_average_file,
    plot = moving_average_plot,
    width = 10,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  message(
    "Saved conditional Poisson result figures."
  )
  
  invisible(
    c(
      single_lag = lag_file,
      moving_average =
        moving_average_file
    )
  )
}