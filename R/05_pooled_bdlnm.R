# ==============================================================================
# Pooled Bayesian distributed lag nonlinear models
# ==============================================================================

source(
  file.path(
    "R",
    "02_prepare_analysis_data.R"
  )
)

check_bdlnm_packages <- function() {
  required_packages <- c(
    "data.table",
    "dlnm",
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

make_bdlnm_exposure_columns <- function(
  variable,
  exposure_type = c(
    "lag",
    "mavg"
  ),
  max_lag = ANALYSIS_MAX_LAG
) {
  exposure_type <- match.arg(
    exposure_type
  )

  if (exposure_type == "lag") {
    make_lag_columns(
      variable,
      max_lag
    )
  } else {
    c(
      variable,
      paste0(
        variable,
        "_mavg",
        seq_len(max_lag)
      )
    )
  }
}

build_crossbasis <- function(
  exposure_history,
  max_lag = ANALYSIS_MAX_LAG,
  lag_df = 3L,
  exposure_knot_probabilities = c(
    0.25,
    0.50,
    0.75
  )
) {
  exposure_values <- as.numeric(
    exposure_history
  )

  exposure_values <- exposure_values[
    is.finite(exposure_values)
  ]

  if (length(unique(exposure_values)) < 10L) {
    stop(
      "Too few unique exposure values.",
      call. = FALSE
    )
  }

  exposure_knots <- as.numeric(
    stats::quantile(
      exposure_values,
      probs = exposure_knot_probabilities,
      names = FALSE,
      na.rm = TRUE
    )
  )

  exposure_boundary <- range(
    exposure_values,
    na.rm = TRUE
  )

  exposure_basis <- list(
    fun = "ns",
    knots = exposure_knots,
    Boundary.knots = exposure_boundary
  )

  lag_basis <- list(
    fun = "ns",
    df = lag_df,
    intercept = TRUE
  )

  crossbasis <- dlnm::crossbasis(
    x = exposure_history,
    lag = c(
      0,
      max_lag
    ),
    argvar = exposure_basis,
    arglag = lag_basis
  )

  crossbasis_matrix <- as.matrix(
    crossbasis
  )

  crossbasis_names <- paste0(
    "cb_",
    seq_len(
      ncol(crossbasis_matrix)
    )
  )

  colnames(
    crossbasis_matrix
  ) <- crossbasis_names

  list(
    matrix = crossbasis_matrix,
    names = crossbasis_names,
    exposure_basis = exposure_basis,
    lag_basis = lag_basis,
    exposure_knots = exposure_knots,
    exposure_boundary = exposure_boundary,
    exposure_values = exposure_values
  )
}

extract_bdlnm_metrics <- function(
  model
) {
  cpo_values <- if (
    !is.null(model$cpo) &&
      !is.null(model$cpo$cpo)
  ) {
    as.numeric(
      model$cpo$cpo
    )
  } else {
    numeric()
  }

  valid_cpo <- is.finite(cpo_values) &
    cpo_values > 0

  data.table::data.table(
    DIC = if (
      !is.null(model$dic$dic)
    ) {
      as.numeric(model$dic$dic)
    } else {
      NA_real_
    },

    effective_parameters_DIC = if (
      !is.null(model$dic$p.eff)
    ) {
      as.numeric(model$dic$p.eff)
    } else {
      NA_real_
    },

    WAIC = if (
      !is.null(model$waic$waic)
    ) {
      as.numeric(model$waic$waic)
    } else {
      NA_real_
    },

    effective_parameters_WAIC = if (
      !is.null(model$waic$p.eff)
    ) {
      as.numeric(model$waic$p.eff)
    } else {
      NA_real_
    },

    CPO_valid = sum(
      valid_cpo
    ),

    LPML = if (
      any(valid_cpo)
    ) {
      sum(
        log(
          cpo_values[
            valid_cpo
          ]
        )
      )
    } else {
      NA_real_
    }
  )
}

fit_pooled_bdlnm <- function(
  data,
  variable,
  exposure_type = c(
    "lag",
    "mavg"
  ),
  max_lag = ANALYSIS_MAX_LAG,
  lag_df = 3L,
  exposure_knot_probabilities = c(
    0.25,
    0.50,
    0.75
  ),
  pc_prior = c(
    3,
    0.01
  ),
  inla_num_threads = "1:1",
  compute_config = TRUE,
  compute_cpo = FALSE,
  return_model = FALSE
) {
  check_bdlnm_packages()

  exposure_type <- match.arg(
    exposure_type
  )

  exposure_columns <-
    make_bdlnm_exposure_columns(
      variable = variable,
      exposure_type = exposure_type,
      max_lag = max_lag
    )

  required_columns <- c(
    "case",
    "global_id",
    "holiday",
    "outbreak_effect",
    exposure_columns
  )

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing B-DLNM columns: ",
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
    case := as.integer(case)
  ]

  model_data[
    ,
    gid_inla := as.integer(
      factor(global_id)
    )
  ]

  exposure_history <- as.matrix(
    model_data[
      ,
      ..exposure_columns
    ]
  )

  storage.mode(
    exposure_history
  ) <- "numeric"

  basis <- build_crossbasis(
    exposure_history = exposure_history,
    max_lag = max_lag,
    lag_df = lag_df,
    exposure_knot_probabilities =
      exposure_knot_probabilities
  )

  basis_data <- data.table::as.data.table(
    basis$matrix
  )

  model_data <- cbind(
    model_data,
    basis_data
  )

  hyper_stratum <- list(
    prec = list(
      prior = "pc.prec",
      param = pc_prior
    )
  )

  model_formula <- stats::as.formula(
    paste0(
      "case ~ 1 + ",
      paste(
        c(
          basis$names,
          "holiday",
          "outbreak_effect"
        ),
        collapse = " + "
      ),
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
      config = compute_config,
      dic = TRUE,
      waic = TRUE,
      cpo = compute_cpo
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

  metrics <- extract_bdlnm_metrics(
    model
  )

  summary <- data.table::data.table(
    variable = variable,
    exposure = ANALYSIS_EXPOSURES$label[
      match(
        variable,
        ANALYSIS_EXPOSURES$variable
      )
    ],
    exposure_type = exposure_type,
    max_lag = max_lag,
    lag_df = lag_df,
    exposure_knots = paste(
      signif(
        basis$exposure_knots,
        7
      ),
      collapse = ";"
    ),
    exposure_minimum =
      basis$exposure_boundary[1],
    exposure_maximum =
      basis$exposure_boundary[2],
    basis_columns = length(
      basis$names
    ),
    n_rows = nrow(
      model_data
    ),
    n_strata = data.table::uniqueN(
      model_data$gid_inla
    ),
    compute_config = compute_config,
    compute_cpo = compute_cpo,
    model_size_mb = as.numeric(
      utils::object.size(model)
    ) / 1024^2,
    elapsed_seconds = elapsed_seconds
  )

  summary <- cbind(
    summary,
    metrics
  )

  fit_light <- list(
    variable = variable,
    exposure_type = exposure_type,
    max_lag = max_lag,
    lag_df = lag_df,
    exposure_knot_probabilities =
      exposure_knot_probabilities,
    crossbasis_names = basis$names,
    exposure_basis = basis$exposure_basis,
    lag_basis = basis$lag_basis,
    exposure_knots = basis$exposure_knots,
    exposure_boundary =
      basis$exposure_boundary,
    exposure_values =
      basis$exposure_values,
    coefficient_mean =
      model$summary.fixed[
        basis$names,
        "mean"
      ]
  )

  output <- list(
    summary = summary,
    fit_light = fit_light
  )

  if (isTRUE(return_model)) {
    output$model <- model
  }

  output
}

run_pooled_bdlnm_model <- function(
  variable,
  exposure_type = c(
    "lag",
    "mavg"
  ),
  config_path = file.path(
    "config",
    "config.R"
  ),
  compute_config = TRUE,
  compute_cpo = FALSE,
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
    modules = "bdlnm"
  )

  set.seed(
    config$computation$random_seed
  )

  fit_pooled_bdlnm(
    data = prepared$data,
    variable = variable,
    exposure_type = exposure_type,
    inla_num_threads =
      config$computation$inla_num_threads,
    compute_config = compute_config,
    compute_cpo = compute_cpo,
    return_model = return_model
  )
}

make_bdlnm_basis_contrast <- function(
    fit_light,
    exposure_value,
    reference_value,
    lag_window
) {
  max_lag <- fit_light$max_lag
  
  if (lag_window < 0L ||
      lag_window > max_lag) {
    stop(
      "`lag_window` is outside the fitted range.",
      call. = FALSE
    )
  }
  
  reference_history <- matrix(
    reference_value,
    nrow = 1L,
    ncol = max_lag + 1L
  )
  
  comparison_history <-
    reference_history
  
  comparison_history[
    1L,
    lag_window + 1L
  ] <- exposure_value
  
  reference_basis <- dlnm::crossbasis(
    x = reference_history,
    lag = c(0, max_lag),
    argvar = fit_light$exposure_basis,
    arglag = fit_light$lag_basis
  )
  
  comparison_basis <- dlnm::crossbasis(
    x = comparison_history,
    lag = c(0, max_lag),
    argvar = fit_light$exposure_basis,
    arglag = fit_light$lag_basis
  )
  
  as.numeric(
    as.matrix(comparison_basis) -
      as.matrix(reference_basis)
  )
}

get_inla_config_contents <- function(model) {
  contents <- as.data.frame(
    model$misc$configs$contents
  )
  
  required_columns <- c(
    "tag",
    "start",
    "length"
  )
  
  if (!all(required_columns %in%
           names(contents))) {
    stop(
      paste(
        "Posterior configuration is unavailable.",
        "Fit with `compute_config = TRUE`."
      ),
      call. = FALSE
    )
  }
  
  contents$clean_tag <- gsub(
    "`",
    "",
    as.character(contents$tag)
  )
  
  contents
}

get_posterior_latent_vector <- function(
    posterior_sample
) {
  latent <- posterior_sample$latent
  
  if (is.matrix(latent)) {
    values <- as.numeric(
      latent[, 1L]
    )
    
    names(values) <- rownames(
      latent
    )
  } else {
    values <- as.numeric(
      latent
    )
    
    names(values) <- names(
      latent
    )
  }
  
  values
}

find_inla_config_row <- function(
    contents,
    target
) {
  clean_target <- gsub(
    "`",
    "",
    target
  )
  
  matching_rows <- which(
    contents$clean_tag ==
      clean_target
  )
  
  if (length(matching_rows) == 1L) {
    return(
      matching_rows
    )
  }
  
  matching_rows <- grep(
    paste0(
      "(^|[:.])",
      clean_target,
      "($|[:.])"
    ),
    contents$clean_tag
  )
  
  if (length(matching_rows) != 1L) {
    stop(
      "Cannot locate posterior component: ",
      target,
      call. = FALSE
    )
  }
  
  matching_rows
}

extract_posterior_component <- function(
    posterior_sample,
    contents,
    target,
    expected_length = 1L
) {
  component_row <- find_inla_config_row(
    contents,
    target
  )
  
  latent <- get_posterior_latent_vector(
    posterior_sample
  )
  
  raw_start <- as.integer(
    contents$start[
      component_row
    ]
  )
  
  component_length <- as.integer(
    contents$length[
      component_row
    ]
  )
  
  minimum_start <- min(
    as.integer(contents$start),
    na.rm = TRUE
  )
  
  start_adjustment <- if (
    minimum_start == 0L
  ) {
    1L
  } else {
    0L
  }
  
  indices <- raw_start +
    start_adjustment +
    seq_len(component_length) -
    1L
  
  values <- as.numeric(
    latent[
      indices
    ]
  )
  
  if (!is.null(expected_length) &&
      length(values) != expected_length) {
    stop(
      target,
      " returned ",
      length(values),
      " values; expected ",
      expected_length,
      ".",
      call. = FALSE
    )
  }
  
  values
}

sample_crossbasis_coefficients <- function(
    model,
    crossbasis_names,
    n_samples = 500L,
    inla_num_threads = "1:1",
    seed = 2026L
) {
  if (n_samples < 1L) {
    stop(
      "`n_samples` must be positive.",
      call. = FALSE
    )
  }
  
  set.seed(seed)
  
  contents <- get_inla_config_contents(
    model
  )
  
  started_at <- proc.time()[[
    "elapsed"
  ]]
  
  posterior_samples <-
    INLA::inla.posterior.sample(
      n = n_samples,
      result = model,
      num.threads = inla_num_threads
    )
  
  coefficient_samples <- matrix(
    NA_real_,
    nrow = n_samples,
    ncol = length(crossbasis_names),
    dimnames = list(
      NULL,
      crossbasis_names
    )
  )
  
  for (sample_id in seq_len(n_samples)) {
    for (coefficient_id in
         seq_along(crossbasis_names)) {
      coefficient_samples[
        sample_id,
        coefficient_id
      ] <- extract_posterior_component(
        posterior_sample =
          posterior_samples[[sample_id]],
        contents = contents,
        target =
          crossbasis_names[coefficient_id],
        expected_length = 1L
      )
    }
  }
  
  sampling_seconds <- unname(
    proc.time()[["elapsed"]] -
      started_at
  )
  
  rm(posterior_samples)
  
  invisible(
    gc()
  )
  
  attr(
    coefficient_samples,
    "sampling_seconds"
  ) <- sampling_seconds
  
  coefficient_samples
}

make_bdlnm_exposure_grid <- function(
    exposure_values,
    probabilities = seq(
      0.01,
      0.99,
      by = 0.01
    )
) {
  values <- as.numeric(
    stats::quantile(
      exposure_values,
      probs = probabilities,
      names = FALSE,
      na.rm = TRUE
    )
  )
  
  keep <- !duplicated(values)
  
  data.table::data.table(
    exposure_value = values[keep],
    exposure_probability =
      probabilities[keep]
  )
}

compute_bdlnm_slice <- function(
    fit_light,
    coefficient_samples,
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
    function(index) {
      contrast <- make_bdlnm_basis_contrast(
        fit_light = fit_light,
        exposure_value =
          grid$exposure_value[index],
        reference_value = reference_value,
        lag_window = lag_window
      )
      
      log_relative_risk <- as.numeric(
        coefficient_samples %*%
          contrast
      )
      
      relative_risk <- exp(
        log_relative_risk
      )
      
      data.table::data.table(
        variable = fit_light$variable,
        exposure_type =
          fit_light$exposure_type,
        lag_window =
          as.integer(lag_window),
        exposure_value =
          grid$exposure_value[index],
        exposure_probability =
          grid$exposure_probability[index],
        reference_value =
          reference_value,
        RR_median = stats::median(
          relative_risk
        ),
        RR_lower = as.numeric(
          stats::quantile(
            relative_risk,
            0.025,
            names = FALSE
          )
        ),
        RR_upper = as.numeric(
          stats::quantile(
            relative_risk,
            0.975,
            names = FALSE
          )
        )
      )
    }
  )
  
  data.table::rbindlist(
    rows
  )
}

build_pooled_bdlnm_specifications <- function() {
  specifications <- data.table::CJ(
    variable = ANALYSIS_EXPOSURES$variable,
    exposure_type = c(
      "lag",
      "mavg"
    ),
    sorted = FALSE
  )
  
  specifications[
    ,
    specification_order := .I
  ]
  
  specifications[
    ,
    model_key := paste(
      variable,
      exposure_type,
      sep = "__"
    )
  ]
  
  specifications[]
}

compute_bdlnm_surface_mean <- function(
    fit_light,
    probabilities = seq(
      0.01,
      0.99,
      by = 0.03
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
  
  rows <- vector(
    "list",
    nrow(grid) *
      (fit_light$max_lag + 1L)
  )
  
  row_id <- 1L
  
  for (exposure_id in
       seq_len(nrow(grid))) {
    for (lag_window in
         0L:fit_light$max_lag) {
      contrast <- make_bdlnm_basis_contrast(
        fit_light = fit_light,
        exposure_value =
          grid$exposure_value[exposure_id],
        reference_value = reference_value,
        lag_window = lag_window
      )
      
      rows[[row_id]] <-
        data.table::data.table(
          variable = fit_light$variable,
          exposure_type =
            fit_light$exposure_type,
          exposure_value =
            grid$exposure_value[exposure_id],
          exposure_probability =
            grid$exposure_probability[
              exposure_id
            ],
          lag_window = lag_window,
          reference_value =
            reference_value,
          RR = exp(
            sum(
              fit_light$coefficient_mean *
                contrast
            )
          )
        )
      
      row_id <- row_id + 1L
    }
  }
  
  data.table::rbindlist(
    rows
  )
}

read_bdlnm_checkpoint <- function(
    path
) {
  if (!file.exists(path)) {
    return(
      data.table::data.table()
    )
  }
  
  data.table::fread(
    path,
    na.strings = c(
      "",
      "NA"
    ),
    showProgress = FALSE
  )
}

replace_bdlnm_checkpoint_rows <- function(
    existing,
    replacement,
    key_to_replace
) {
  if (nrow(existing) > 0L &&
      "model_key" %in%
      names(existing)) {
    rows_to_keep <- existing[[
      "model_key"
    ]] != key_to_replace
    
    existing <- existing[
      rows_to_keep
    ]
  }
  
  data.table::rbindlist(
    list(
      existing,
      replacement
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

write_bdlnm_checkpoint <- function(
    data,
    path,
    order_columns
) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  available_order_columns <- intersect(
    order_columns,
    names(data)
  )
  
  if (length(available_order_columns) > 0L) {
    data.table::setorderv(
      data,
      available_order_columns
    )
  }
  
  data.table::fwrite(
    data,
    path,
    na = ""
  )
  
  invisible(path)
}

run_all_pooled_bdlnm_models <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    metrics_file = file.path(
      "results",
      "pooled_bdlnm_metrics.csv"
    ),
    slices_file = file.path(
      "results",
      "pooled_bdlnm_slices.csv"
    ),
    surfaces_file = file.path(
      "results",
      "pooled_bdlnm_surfaces.csv"
    ),
    resume = TRUE
) {
  config <- load_analysis_config(
    config_path
  )
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "bdlnm"
  )
  
  specifications <-
    build_pooled_bdlnm_specifications()
  
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
  
  surfaces <- if (isTRUE(resume)) {
    read_bdlnm_checkpoint(
      surfaces_file
    )
  } else {
    data.table::data.table()
  }
  
  completed_keys <- if (
    nrow(metrics) > 0L &&
    "status" %in% names(metrics)
  ) {
    metrics[
      status == "ok",
      model_key
    ]
  } else {
    character()
  }
  
  total_models <- nrow(
    specifications
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
        "[%d/%d] Fitting pooled B-DLNM: %s",
        index,
        total_models,
        model_key
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
          inla_num_threads =
            config$computation$inla_num_threads,
          compute_config = TRUE,
          compute_cpo = FALSE,
          return_model = TRUE
        )
        
        coefficient_samples <-
          sample_crossbasis_coefficients(
            model = fit$model,
            crossbasis_names =
              fit$fit_light$
              crossbasis_names,
            n_samples =
              config$computation$
              posterior_samples,
            inla_num_threads =
              config$computation$
              inla_num_threads,
            seed =
              config$computation$
              random_seed +
              specification$
              specification_order
          )
        
        slice_window <- if (
          specification$exposure_type ==
          "lag"
        ) {
          3L
        } else {
          4L
        }
        
        model_slice <- compute_bdlnm_slice(
          fit_light = fit$fit_light,
          coefficient_samples =
            coefficient_samples,
          lag_window = slice_window
        )
        
        model_surface <-
          compute_bdlnm_surface_mean(
            fit$fit_light
          )
        
        model_slice[
          ,
          `:=`(
            model_key = model_key,
            specification_order =
              specification$
              specification_order
          )
        ]
        
        model_surface[
          ,
          `:=`(
            model_key = model_key,
            specification_order =
              specification$
              specification_order
          )
        ]
        
        model_metrics <- data.table::copy(
          fit$summary
        )
        
        model_metrics[
          ,
          `:=`(
            model_key = model_key,
            specification_order =
              specification$
              specification_order,
            posterior_samples =
              nrow(coefficient_samples),
            sampling_seconds = attr(
              coefficient_samples,
              "sampling_seconds"
            ),
            slice_window =
              slice_window,
            status = "ok",
            error_message =
              NA_character_
          )
        ]
        
        rm(
          fit,
          coefficient_samples
        )
        
        invisible(
          gc()
        )
        
        list(
          success = TRUE,
          metrics = model_metrics,
          slice = model_slice,
          surface = model_surface
        )
      },
      error = function(error) {
        list(
          success = FALSE,
          error_message =
            conditionMessage(error)
        )
      }
    )
    
    if (isTRUE(outcome$success)) {
      slices <-
        replace_bdlnm_checkpoint_rows(
          existing = slices,
          replacement = outcome$slice,
          key_to_replace = model_key
        )
      
      write_bdlnm_checkpoint(
        slices,
        slices_file,
        c(
          "specification_order",
          "exposure_probability"
        )
      )
      
      surfaces <-
        replace_bdlnm_checkpoint_rows(
          existing = surfaces,
          replacement =
            outcome$surface,
          key_to_replace = model_key
        )
      
      write_bdlnm_checkpoint(
        surfaces,
        surfaces_file,
        c(
          "specification_order",
          "exposure_probability",
          "lag_window"
        )
      )
      
      # Metrics are written last. A model is considered
      # complete only after its aggregate outputs exist.
      metrics <-
        replace_bdlnm_checkpoint_rows(
          existing = metrics,
          replacement =
            outcome$metrics,
          key_to_replace = model_key
        )
      
      write_bdlnm_checkpoint(
        metrics,
        metrics_file,
        "specification_order"
      )
      
      message(
        sprintf(
          "Completed %s: fit %.1f s, sampling %.1f s.",
          model_key,
          outcome$metrics$
            elapsed_seconds,
          outcome$metrics$
            sampling_seconds
        )
      )
    } else {
      failed_row <- data.table::data.table(
        variable =
          specification$variable,
        exposure_type =
          specification$exposure_type,
        model_key = model_key,
        specification_order =
          specification$
          specification_order,
        status = "failed",
        error_message =
          outcome$error_message
      )
      
      metrics <-
        replace_bdlnm_checkpoint_rows(
          existing = metrics,
          replacement = failed_row,
          key_to_replace = model_key
        )
      
      write_bdlnm_checkpoint(
        metrics,
        metrics_file,
        "specification_order"
      )
      
      message(
        "Failed ",
        model_key,
        ": ",
        outcome$error_message
      )
    }
    
    rm(outcome)
    
    invisible(
      gc()
    )
  }
  
  message(
    "Pooled B-DLNM batch complete: ",
    metrics[status == "ok", .N],
    " successful, ",
    metrics[status != "ok", .N],
    " failed."
  )
  
  invisible(
    list(
      metrics = metrics,
      slices = slices,
      surfaces = surfaces
    )
  )
}

plot_pooled_bdlnm_slices <- function(
    slices,
    exposure_type = c(
      "lag",
      "mavg"
    )
) {
  if (!requireNamespace(
    "ggplot2",
    quietly = TRUE
  ) ||
  !requireNamespace(
    "scales",
    quietly = TRUE
  )) {
    stop(
      "Packages `ggplot2` and `scales` are required.",
      call. = FALSE
    )
  }
  
  selected_exposure_type <- match.arg(
    exposure_type
  )
  
  rows_to_plot <- slices[[
    "exposure_type"
  ]] == selected_exposure_type
  
  plot_data <- data.table::copy(
    slices[
      rows_to_plot
    ]
  )
  
  exposure_labels <- paste0(
    ANALYSIS_EXPOSURES$label,
    " (",
    ANALYSIS_EXPOSURES$unit,
    ")"
  )
  
  names(exposure_labels) <-
    ANALYSIS_EXPOSURES$variable
  
  plot_data[
    ,
    facet_label := factor(
      exposure_labels[variable],
      levels = exposure_labels
    )
  ]
  
  settings <- if (
    selected_exposure_type == "lag"
  ) {
    list(
      title =
        "Nonlinear exposure-response associations at lag 3",
      subtitle = paste(
        "Pooled Bayesian distributed lag nonlinear models;",
        "reference exposure = P50"
      ),
      colour = "#397D8C"
    )
  } else {
    list(
      title =
        "Nonlinear exposure-response associations at moving average 4",
      subtitle = paste(
        "Pooled Bayesian distributed lag nonlinear models;",
        "reference exposure = P50"
      ),
      colour = "#C8553D"
    )
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = exposure_value,
      y = RR_median
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      colour = "grey45",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(
        xintercept = reference_value
      ),
      linetype = "dotted",
      colour = "grey55",
      linewidth = 0.4
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = RR_lower,
        ymax = RR_upper
      ),
      fill = settings$colour,
      alpha = 0.18
    ) +
    ggplot2::geom_line(
      colour = settings$colour,
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(
        facet_label
      ),
      ncol = 3,
      scales = "free"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.01
      )
    ) +
    ggplot2::labs(
      title = settings$title,
      subtitle = settings$subtitle,
      x = "Exposure value",
      y = "Relative risk",
      caption = paste(
        "Lines show posterior median RRs;",
        "bands show 95% credible intervals.",
        "Vertical dotted lines indicate exposure-specific P50 values.",
        "Axis scales vary across exposures."
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

save_pooled_bdlnm_slice_figures <- function(
    slices_file = file.path(
      "results",
      "pooled_bdlnm_slices.csv"
    ),
    figures_dir = "figures"
) {
  if (!file.exists(slices_file)) {
    stop(
      "Cannot find B-DLNM slices: ",
      slices_file,
      call. = FALSE
    )
  }
  
  slices <- data.table::fread(
    slices_file,
    showProgress = FALSE
  )
  
  lag_plot <- plot_pooled_bdlnm_slices(
    slices,
    exposure_type = "lag"
  )
  
  moving_average_plot <-
    plot_pooled_bdlnm_slices(
      slices,
      exposure_type = "mavg"
    )
  
  dir.create(
    figures_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  lag_file <- file.path(
    figures_dir,
    "pooled_bdlnm_single_lag3_slices.png"
  )
  
  moving_average_file <- file.path(
    figures_dir,
    "pooled_bdlnm_mavg4_slices.png"
  )
  
  ggplot2::ggsave(
    lag_file,
    lag_plot,
    width = 10,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  ggplot2::ggsave(
    moving_average_file,
    moving_average_plot,
    width = 10,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  message(
    "Saved pooled B-DLNM slice figures."
  )
  
  invisible(
    c(
      single_lag = lag_file,
      moving_average =
        moving_average_file
    )
  )
}