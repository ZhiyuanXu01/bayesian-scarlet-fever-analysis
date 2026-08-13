# ==============================================================================
# Subgroup analysis
#
# Subgroups are assigned using the case record in each matched stratum. All
# case and referent days belonging to an eligible stratum are then retained.
# ==============================================================================

source(
  file.path(
    "R",
    "04_bayesian_conditional_poisson.R"
  )
)

SUBGROUP_WINDOWS <- data.table::data.table(
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

build_subgroup_membership <- function(
    data
) {
  required_columns <- c(
    "case",
    "global_id",
    "age",
    "gender",
    "year",
    "season"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing subgroup columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  subgroup_data <- data.table::copy(
    data
  )
  
  stratum_check <- subgroup_data[
    ,
    .(
      case_days = sum(
        case == 1L,
        na.rm = TRUE
      )
    ),
    by = global_id
  ]
  
  invalid_strata <- stratum_check[
    case_days != 1L
  ]
  
  if (nrow(invalid_strata) > 0L) {
    stop(
      nrow(invalid_strata),
      " matched strata do not contain exactly one case day.",
      call. = FALSE
    )
  }
  
  case_records <- subgroup_data[
    case == 1L,
    .(
      global_id,
      age = suppressWarnings(
        as.numeric(age)
      ),
      gender = toupper(
        trimws(
          as.character(gender)
        )
      ),
      year = as.integer(year),
      season = as.character(season)
    )
  ]
  
  case_records[
    ,
    age_group := data.table::fcase(
      is.finite(age) & age < 6,
      "Age <6",
      is.finite(age) & age >= 6,
      "Age >=6",
      default = NA_character_
    )
  ]
  
  case_records[
    ,
    sex_group := data.table::fcase(
      gender == "M",
      "Male",
      gender == "F",
      "Female",
      default = NA_character_
    )
  ]
  
  case_records[
    ,
    period_group := data.table::fcase(
      year >= 2005L & year <= 2010L,
      "2005-2010",
      year >= 2011L & year <= 2019L,
      "2011-2019",
      year >= 2020L & year <= 2024L,
      "2020-2024",
      default = NA_character_
    )
  ]
  
  valid_seasons <- c(
    "Spring",
    "Summer",
    "Autumn",
    "Winter"
  )
  
  case_records[
    !season %in% valid_seasons,
    season := NA_character_
  ]
  
  overall_membership <- case_records[
    ,
    .(
      global_id,
      group_type = "Overall",
      subgroup = "Overall",
      group_order = 1L,
      subgroup_order = 1L
    )
  ]
  
  age_membership <- case_records[
    !is.na(age_group),
    .(
      global_id,
      group_type = "Age",
      subgroup = age_group,
      group_order = 2L,
      subgroup_order = match(
        age_group,
        c(
          "Age <6",
          "Age >=6"
        )
      )
    )
  ]
  
  sex_membership <- case_records[
    !is.na(sex_group),
    .(
      global_id,
      group_type = "Sex",
      subgroup = sex_group,
      group_order = 3L,
      subgroup_order = match(
        sex_group,
        c(
          "Male",
          "Female"
        )
      )
    )
  ]
  
  period_membership <- case_records[
    !is.na(period_group),
    .(
      global_id,
      group_type = "Period",
      subgroup = period_group,
      group_order = 4L,
      subgroup_order = match(
        period_group,
        c(
          "2005-2010",
          "2011-2019",
          "2020-2024"
        )
      )
    )
  ]
  
  season_membership <- case_records[
    !is.na(season),
    .(
      global_id,
      group_type = "Season",
      subgroup = season,
      group_order = 5L,
      subgroup_order = match(
        season,
        valid_seasons
      )
    )
  ]
  
  membership <- data.table::rbindlist(
    list(
      overall_membership,
      age_membership,
      sex_membership,
      period_membership,
      season_membership
    ),
    use.names = TRUE
  )
  
  membership[
    ,
    group_key := paste(
      group_type,
      subgroup,
      sep = "__"
    )
  ]
  
  data.table::setorder(
    membership,
    group_order,
    subgroup_order,
    global_id
  )
  
  duplicated_memberships <- membership[
    duplicated(
      membership[
        ,
        .(
          group_key,
          global_id
        )
      ]
    )
  ]
  
  if (nrow(duplicated_memberships) > 0L) {
    stop(
      "At least one stratum was assigned twice to the same subgroup.",
      call. = FALSE
    )
  }
  
  missingness <- case_records[
    ,
    .(
      case_records = .N,
      missing_age = sum(
        is.na(age_group)
      ),
      unclassified_gender = sum(
        is.na(sex_group)
      ),
      unclassified_period = sum(
        is.na(period_group)
      ),
      missing_season = sum(
        is.na(season)
      )
    )
  ]
  
  list(
    case_records = case_records,
    membership = membership,
    missingness = missingness
  )
}

summarise_subgroup_samples <- function(
    data,
    membership
) {
  stratum_summary <- data[
    ,
    .(
      analytical_rows = .N,
      case_days = sum(
        case == 1L,
        na.rm = TRUE
      )
    ),
    by = global_id
  ]
  
  joined <- merge(
    membership,
    stratum_summary,
    by = "global_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  diagnostics <- joined[
    ,
    .(
      matched_strata =
        data.table::uniqueN(global_id),
      analytical_rows = sum(
        analytical_rows
      ),
      case_days = sum(
        case_days
      ),
      minimum_rows_per_stratum = min(
        analytical_rows
      ),
      median_rows_per_stratum = as.numeric(
        stats::median(
          analytical_rows
        )
      ),
      maximum_rows_per_stratum = max(
        analytical_rows
      ),
      strata_without_one_case_day = sum(
        case_days != 1L
      )
    ),
    by = .(
      group_key,
      group_type,
      subgroup,
      group_order,
      subgroup_order
    )
  ]
  
  data.table::setorder(
    diagnostics,
    group_order,
    subgroup_order
  )
  
  diagnostics[]
}

subset_subgroup_strata <- function(
    data,
    membership,
    group_key
) {
  requested_group_key <- as.character(
    group_key
  )
  
  selected_strata <- membership[
    membership[["group_key"]] ==
      requested_group_key,
    "global_id",
    with = FALSE
  ][["global_id"]]
  
  if (length(selected_strata) == 0L) {
    stop(
      "No matched strata were assigned to subgroup: ",
      requested_group_key,
      call. = FALSE
    )
  }
  
  subgroup_data <- data[
    data[["global_id"]] %in%
      selected_strata
  ]
  
  subgroup_check <- subgroup_data[
    ,
    .(
      case_days = sum(
        case == 1L,
        na.rm = TRUE
      )
    ),
    by = global_id
  ]
  
  if (nrow(subgroup_check) !=
      length(selected_strata)) {
    stop(
      "Not all selected matched strata were retained.",
      call. = FALSE
    )
  }
  
  if (any(
    subgroup_check$case_days != 1L
  )) {
    stop(
      "Subgroup filtering broke the matched-stratum structure.",
      call. = FALSE
    )
  }
  
  subgroup_data
}

run_subgroup_diagnostics <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    output_file = NULL
) {
  config <- load_analysis_config(
    config_path
  )
  
  if (is.null(output_file)) {
    output_file <- file.path(
      config$output$results_dir,
      "subgroup_sample_sizes.csv"
    )
  }
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "subgroup"
  )
  
  subgroup_objects <- build_subgroup_membership(
    prepared$data
  )
  
  diagnostics <- summarise_subgroup_samples(
    data = prepared$data,
    membership =
      subgroup_objects$membership
  )
  
  dir.create(
    dirname(
      output_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Only aggregate subgroup counts are written. Stratum membership remains local.
  data.table::fwrite(
    diagnostics,
    output_file
  )
  
  message(
    "Subgroup diagnostics completed for ",
    nrow(diagnostics),
    " groups."
  )
  
  list(
    diagnostics = diagnostics,
    missingness =
      subgroup_objects$missingness
  )
}

run_subgroup_model <- function(
    group_key,
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
  requested_group_key <- as.character(
    group_key
  )
  
  requested_exposure_type <- match.arg(
    exposure_type
  )
  
  requested_window <- as.integer(
    window
  )
  
  window_information <- SUBGROUP_WINDOWS[
    SUBGROUP_WINDOWS[["exposure_type"]] ==
      requested_exposure_type &
      SUBGROUP_WINDOWS[["window"]] ==
      requested_window
  ]
  
  if (nrow(window_information) != 1L) {
    stop(
      "The requested subgroup-analysis window is not prespecified.",
      call. = FALSE
    )
  }
  
  config <- load_analysis_config(
    config_path
  )
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "subgroup"
  )
  
  subgroup_objects <- build_subgroup_membership(
    prepared$data
  )
  
  selected_membership <-
    subgroup_objects$membership[
      subgroup_objects$membership[["group_key"]] ==
        requested_group_key
    ]
  
  if (nrow(selected_membership) == 0L) {
    stop(
      "Unknown subgroup key: ",
      requested_group_key,
      call. = FALSE
    )
  }
  
  group_information <- unique(
    selected_membership[
      ,
      .(
        group_key,
        group_type,
        subgroup,
        group_order,
        subgroup_order
      )
    ]
  )
  
  if (nrow(group_information) != 1L) {
    stop(
      "The subgroup key does not identify exactly one group.",
      call. = FALSE
    )
  }
  
  subgroup_data <- subset_subgroup_strata(
    data = prepared$data,
    membership =
      subgroup_objects$membership,
    group_key =
      requested_group_key
  )
  
  fit <- fit_conditional_poisson(
    data = subgroup_data,
    variable = variable,
    exposure_type =
      requested_exposure_type,
    window =
      requested_window,
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
      group_key =
        group_information$group_key,
      group_type =
        group_information$group_type,
      subgroup =
        group_information$subgroup,
      group_order =
        group_information$group_order,
      subgroup_order =
        group_information$subgroup_order,
      window_label =
        window_information$window_label
    )
  ]
  
  data.table::setcolorder(
    summary,
    c(
      "group_key",
      "group_type",
      "subgroup",
      "group_order",
      "subgroup_order",
      "window_label",
      setdiff(
        names(summary),
        c(
          "group_key",
          "group_type",
          "subgroup",
          "group_order",
          "subgroup_order",
          "window_label"
        )
      )
    )
  )
  
  output <- list(
    summary = summary
  )
  
  if (isTRUE(return_model)) {
    output$model <- fit$model
  }
  
  output
}

build_subgroup_model_specifications <- function(
    subgroup_diagnostics
) {
  group_specifications <- subgroup_diagnostics[
    ,
    .(
      group_key,
      group_type,
      subgroup,
      group_order,
      subgroup_order
    )
  ]
  
  specifications <- data.table::rbindlist(
    lapply(
      seq_len(nrow(group_specifications)),
      function(group_index) {
        current_group <- group_specifications[
          group_index
        ]
        
        data.table::rbindlist(
          lapply(
            seq_len(nrow(SUBGROUP_WINDOWS)),
            function(window_index) {
              current_window <- SUBGROUP_WINDOWS[
                window_index
              ]
              
              data.table::data.table(
                group_key =
                  current_group$group_key,
                group_type =
                  current_group$group_type,
                subgroup =
                  current_group$subgroup,
                group_order =
                  current_group$group_order,
                subgroup_order =
                  current_group$subgroup_order,
                exposure_type =
                  current_window$exposure_type,
                window =
                  current_window$window,
                window_label =
                  current_window$window_label,
                window_order =
                  window_index,
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
    group_order,
    subgroup_order,
    window_order,
    exposure_order
  )
  
  specifications[
    ,
    specification_order := .I
  ]
  
  specifications[
    ,
    model_key := sprintf(
      "g%02d_s%02d__%s__%s__%02d",
      group_order,
      subgroup_order,
      variable,
      exposure_type,
      window
    )
  ]
  
  if (data.table::uniqueN(
    specifications$model_key
  ) != nrow(specifications)) {
    stop(
      "Subgroup model keys are not unique.",
      call. = FALSE
    )
  }
  
  specifications[]
}

run_all_subgroup_models <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    output_file = NULL,
    sample_size_file = NULL,
    resume = TRUE
) {
  check_conditional_poisson_packages()
  
  config <- load_analysis_config(
    config_path
  )
  
  if (is.null(output_file)) {
    output_file <- file.path(
      config$output$results_dir,
      "subgroup_models.csv"
    )
  }
  
  if (is.null(sample_size_file)) {
    sample_size_file <- file.path(
      config$output$results_dir,
      "subgroup_sample_sizes.csv"
    )
  }
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "subgroup"
  )
  
  subgroup_objects <- build_subgroup_membership(
    prepared$data
  )
  
  subgroup_diagnostics <-
    summarise_subgroup_samples(
      data = prepared$data,
      membership =
        subgroup_objects$membership
    )
  
  dir.create(
    dirname(
      sample_size_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  data.table::fwrite(
    subgroup_diagnostics,
    sample_size_file
  )
  
  specifications <-
    build_subgroup_model_specifications(
      subgroup_diagnostics
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
  
  cached_group_key <- NA_character_
  cached_group_data <- NULL
  
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
    
    current_group_key <- as.character(
      specification$group_key
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
    
    if (!identical(
      cached_group_key,
      current_group_key
    )) {
      cached_group_data <- subset_subgroup_strata(
        data = prepared$data,
        membership =
          subgroup_objects$membership,
        group_key =
          current_group_key
      )
      
      cached_group_key <- current_group_key
      
      message(
        sprintf(
          "Prepared subgroup %s: %d rows, %d strata.",
          current_group_key,
          nrow(cached_group_data),
          data.table::uniqueN(
            cached_group_data$global_id
          )
        )
      )
    }
    
    message(
      sprintf(
        "[%d/%d] Fitting subgroup model: %s",
        index,
        total_models,
        current_model_key
      )
    )
    
    result_row <- tryCatch(
      {
        fit <- fit_conditional_poisson(
          data = cached_group_data,
          variable =
            specification$variable,
          exposure_type =
            specification$exposure_type,
          window =
            specification$window,
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
            group_key =
              current_group_key,
            group_type =
              as.character(
                specification$group_type
              ),
            subgroup =
              as.character(
                specification$subgroup
              ),
            group_order =
              as.integer(
                specification$group_order
              ),
            subgroup_order =
              as.integer(
                specification$subgroup_order
              ),
            window_label =
              as.character(
                specification$window_label
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
          group_key =
            current_group_key,
          group_type =
            as.character(
              specification$group_type
            ),
          subgroup =
            as.character(
              specification$subgroup
            ),
          group_order =
            as.integer(
              specification$group_order
            ),
          subgroup_order =
            as.integer(
              specification$subgroup_order
            ),
          window_label =
            as.character(
              specification$window_label
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
  
  batch_diagnostics <- checkpoint[
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
        "Subgroup batch complete: ",
        "%d successful, %d failed."
      ),
      nrow(successful_models),
      nrow(failed_models)
    )
  )
  
  list(
    sample_sizes =
      subgroup_diagnostics,
    models = checkpoint,
    diagnostics =
      batch_diagnostics
  )
}

compute_subgroup_iqr_estimates <- function(
    data,
    membership,
    model_results
) {
  successful_models <- data.table::copy(
    model_results[
      model_results[["status"]] == "ok"
    ]
  )
  
  if (nrow(successful_models) == 0L) {
    stop(
      "No successful subgroup models were provided.",
      call. = FALSE
    )
  }
  
  required_columns <- c(
    "group_key",
    "term",
    "coefficient",
    "lower_95_CrI",
    "upper_95_CrI"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(successful_models)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing subgroup-result columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  group_keys <- unique(
    successful_models$group_key
  )
  
  scaled_results <- lapply(
    group_keys,
    function(current_group_key) {
      group_models <- data.table::copy(
        successful_models[
          successful_models[["group_key"]] ==
            current_group_key
        ]
      )
      
      group_data <- subset_subgroup_strata(
        data = data,
        membership = membership,
        group_key = current_group_key
      )
      
      missing_terms <- setdiff(
        unique(
          group_models$term
        ),
        names(group_data)
      )
      
      if (length(missing_terms) > 0L) {
        stop(
          "Missing exposure terms for ",
          current_group_key,
          ": ",
          paste(missing_terms, collapse = ", "),
          call. = FALSE
        )
      }
      
      group_models[
        ,
        exposure_IQR := vapply(
          term,
          function(current_term) {
            stats::IQR(
              group_data[[
                current_term
              ]],
              na.rm = TRUE
            )
          },
          FUN.VALUE = numeric(1)
        )
      ]
      
      if (any(
        !is.finite(
          group_models$exposure_IQR
        ) |
        group_models$exposure_IQR <= 0
      )) {
        stop(
          "A non-positive exposure IQR was found for ",
          current_group_key,
          ".",
          call. = FALSE
        )
      }
      
      group_models[
        ,
        `:=`(
          RR_per_IQR = exp(
            coefficient *
              exposure_IQR
          ),
          lower_95_CrI_per_IQR = exp(
            log(
              lower_95_CrI
            ) *
              exposure_IQR
          ),
          upper_95_CrI_per_IQR = exp(
            log(
              upper_95_CrI
            ) *
              exposure_IQR
          )
        )
      ]
      
      group_models[
        ,
        interval_excludes_one :=
          lower_95_CrI_per_IQR > 1 |
          upper_95_CrI_per_IQR < 1
      ]
      
      group_models
    }
  )
  
  scaled_results <- data.table::rbindlist(
    scaled_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  data.table::setorder(
    scaled_results,
    group_order,
    subgroup_order,
    exposure_type,
    window,
    variable
  )
  
  scaled_results[]
}

run_subgroup_iqr_scaling <- function(
    config_path = file.path(
      "config",
      "config.R"
    ),
    model_file = NULL,
    output_file = NULL
) {
  config <- load_analysis_config(
    config_path
  )
  
  if (is.null(model_file)) {
    model_file <- file.path(
      config$output$results_dir,
      "subgroup_models.csv"
    )
  }
  
  if (is.null(output_file)) {
    output_file <- file.path(
      config$output$results_dir,
      "subgroup_iqr_scaled_estimates.csv"
    )
  }
  
  if (!file.exists(model_file)) {
    stop(
      "Cannot find subgroup model results: ",
      model_file,
      call. = FALSE
    )
  }
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = "subgroup"
  )
  
  subgroup_objects <- build_subgroup_membership(
    prepared$data
  )
  
  model_results <- data.table::fread(
    model_file,
    showProgress = FALSE
  )
  
  scaled_results <- compute_subgroup_iqr_estimates(
    data = prepared$data,
    membership =
      subgroup_objects$membership,
    model_results =
      model_results
  )
  
  dir.create(
    dirname(
      output_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  data.table::fwrite(
    scaled_results,
    output_file
  )
  
  diagnostics <- scaled_results[
    ,
    .(
      models = .N,
      minimum_RR_per_IQR = min(
        RR_per_IQR
      ),
      maximum_RR_per_IQR = max(
        RR_per_IQR
      ),
      intervals_excluding_one = sum(
        interval_excludes_one
      )
    ),
    by = window_label
  ]
  
  message(
    "Wrote ",
    nrow(scaled_results),
    " IQR-scaled subgroup estimates."
  )
  
  list(
    estimates = scaled_results,
    diagnostics = diagnostics
  )
}

plot_subgroup_iqr_heatmap <- function(
    scaled_results,
    output_file = file.path(
      "figures",
      "subgroup_iqr_scaled_rr_heatmap.png"
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
    scaled_results[
      scaled_results[["status"]] == "ok"
    ]
  )
  
  if (nrow(plot_data) == 0L) {
    stop(
      "No successful subgroup estimates were provided.",
      call. = FALSE
    )
  }
  
  short_exposure_labels <- c(
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
    short_exposure_labels[
      ANALYSIS_EXPOSURES$variable
    ]
  )
  
  plot_data[
    ,
    exposure_short := unname(
      short_exposure_labels[
        variable
      ]
    )
  ]
  
  plot_data[
    ,
    subgroup_display := data.table::fcase(
      group_type == "Overall",
      "Overall",
      group_type == "Age",
      paste0(
        "Age: ",
        sub(
          "^Age ",
          "",
          subgroup
        )
      ),
      group_type == "Sex",
      paste0(
        "Sex: ",
        subgroup
      ),
      group_type == "Period",
      paste0(
        "Period: ",
        subgroup
      ),
      group_type == "Season",
      paste0(
        "Season: ",
        subgroup
      ),
      default = paste(
        group_type,
        subgroup,
        sep = ": "
      )
    )
  ]
  
  subgroup_levels <- unique(
    plot_data[
      order(
        group_order,
        subgroup_order
      ),
      subgroup_display
    ]
  )
  
  plot_data[
    ,
    exposure_short := factor(
      exposure_short,
      levels = exposure_levels
    )
  ]
  
  plot_data[
    ,
    subgroup_display := factor(
      subgroup_display,
      levels = rev(
        subgroup_levels
      )
    )
  ]
  
  plot_data[
    ,
    window_label := factor(
      window_label,
      levels =
        SUBGROUP_WINDOWS$window_label
    )
  ]
  
  plot_data[
    ,
    cell_label := ifelse(
      interval_excludes_one,
      sprintf(
        "%.3f*",
        RR_per_IQR
      ),
      sprintf(
        "%.3f",
        RR_per_IQR
      )
    )
  ]
  
  maximum_deviation <- max(
    abs(
      plot_data$RR_per_IQR -
        1
    ),
    na.rm = TRUE
  )
  
  fill_limits <- c(
    1 - maximum_deviation,
    1 + maximum_deviation
  )
  
  subgroup_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = exposure_short,
      y = subgroup_display,
      fill = RR_per_IQR
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
      size = 2.7
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
      name = "RR per\nsubgroup IQR"
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Environmental associations across prespecified subgroups",
      subtitle = paste0(
        "Bayesian conditional Poisson estimates ",
        "scaled to each subgroup's exposure IQR"
      ),
      x = "Environmental exposure",
      y = NULL,
      caption = paste0(
        "Matched strata were assigned using the case record and retained ",
        "with all referent days. ",
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
    plot = subgroup_plot,
    width = 15,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  message(
    "Wrote subgroup heatmap to ",
    output_file,
    "."
  )
  
  invisible(
    subgroup_plot
  )
}