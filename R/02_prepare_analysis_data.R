# ==============================================================================
# Prepare the restricted analytical dataset
#
# The workflow reads only variables required by the requested modules.
# Moving-average exposures are reconstructed from lag 0 through lag k.
# No individual-level data are written to disk.
# ==============================================================================

source(
  file.path(
    "R",
    "01_data_contract.R"
  )
)

load_analysis_config <- function(
  path = file.path(
    "config",
    "config.R"
  )
) {
  if (!file.exists(path)) {
    stop(
      "Cannot find ",
      path,
      ". Copy config/config.example.R to config/config.R first.",
      call. = FALSE
    )
  }

  config_environment <- new.env(
    parent = baseenv()
  )

  sys.source(
    path,
    envir = config_environment
  )

  if (!exists(
    "analysis_config",
    envir = config_environment,
    inherits = FALSE
  )) {
    stop(
      "`analysis_config` was not defined in ",
      path,
      ".",
      call. = FALSE
    )
  }

  config_environment$analysis_config
}

resolve_input_path <- function(path) {
  if (!is.character(path) ||
      length(path) != 1L ||
      !nzchar(path)) {
    stop(
      "The configured analysis-data path is invalid.",
      call. = FALSE
    )
  }

  is_absolute <- grepl(
    "^([A-Za-z]:[/\\\\]|/|\\\\\\\\)",
    path
  )

  resolved_path <- if (is_absolute) {
    path
  } else {
    file.path(
      getwd(),
      path
    )
  }

  if (!file.exists(resolved_path)) {
    stop(
      "Cannot find analysis data: ",
      resolved_path,
      call. = FALSE
    )
  }

  normalizePath(
    resolved_path,
    winslash = "/",
    mustWork = TRUE
  )
}

read_restricted_analysis_data <- function(
  data_file,
  modules = c(
    "descriptive",
    "linear",
    "bdlnm",
    "spatial",
    "subgroup"
  ),
  max_lag = ANALYSIS_MAX_LAG
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

  input_path <- resolve_input_path(
    data_file
  )

  available_columns <- names(
    data.table::fread(
      input_path,
      nrows = 0L,
      showProgress = FALSE
    )
  )

  required_columns <- analysis_required_columns(
    modules = modules,
    max_lag = max_lag
  )

  missing_columns <- setdiff(
    required_columns,
    available_columns
  )

  if (length(missing_columns) > 0L) {
    stop(
      "The analytical dataset is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  message(
    "Reading ",
    length(required_columns),
    " required columns from the restricted dataset."
  )

  data.table::fread(
    input_path,
    select = required_columns,
    na.strings = c("", "NA"),
    showProgress = interactive()
  )
}

standardise_analysis_variables <- function(data) {
  data.table::setDT(data)

  data[
    ,
    case := as.integer(case)
  ]

  data[
    ,
    global_id := as.character(global_id)
  ]

  if ("time" %in% names(data)) {
    data[
      ,
      time := as.Date(time)
    ]
  }

  for (column in intersect(
    c("province", "gender"),
    names(data)
  )) {
    data.table::set(
      data,
      j = column,
      value = as.character(
        data[[column]]
      )
    )
  }

  data
}

derive_time_variables <- function(data) {
  if (!"time" %in% names(data)) {
    return(data)
  }

  data[
    ,
    year := as.integer(
      format(time, "%Y")
    )
  ]

  data[
    ,
    month := as.integer(
      format(time, "%m")
    )
  ]

  data[
    ,
    season := data.table::fcase(
      month %in% c(3L, 4L, 5L),
      "Spring",

      month %in% c(6L, 7L, 8L),
      "Summer",

      month %in% c(9L, 10L, 11L),
      "Autumn",

      month %in% c(12L, 1L, 2L),
      "Winter",

      default = NA_character_
    )
  ]

  data
}

derive_moving_average_exposures <- function(
  data,
  max_lag = ANALYSIS_MAX_LAG
) {
  data.table::setDT(data)

  for (variable in ANALYSIS_EXPOSURES$variable) {
    lag_columns <- make_lag_columns(
      variable,
      max_lag
    )

    missing_columns <- setdiff(
      lag_columns,
      names(data)
    )

    if (length(missing_columns) > 0L) {
      stop(
        "Cannot derive moving averages for ",
        variable,
        "; missing: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }

    lag_matrix <- as.matrix(
      data[
        ,
        ..lag_columns
      ]
    )

    storage.mode(lag_matrix) <- "numeric"

    for (lag_id in seq_len(max_lag)) {
      window_columns <- seq_len(
        lag_id + 1L
      )

      exposure_window <- lag_matrix[
        ,
        window_columns,
        drop = FALSE
      ]

      available_days <- rowSums(
        !is.na(exposure_window)
      )

      moving_average <- rowMeans(
        exposure_window,
        na.rm = TRUE
      )

      moving_average[
        available_days == 0L
      ] <- NA_real_

      output_column <- paste0(
        variable,
        "_mavg",
        lag_id
      )

      data.table::set(
        data,
        j = output_column,
        value = moving_average
      )
    }
  }

  data
}

summarise_analysis_data <- function(data) {
  stratum_summary <- data[
    ,
    .(
      case_days = sum(
        case == 1L,
        na.rm = TRUE
      )
    ),
    by = global_id
  ]

  exposure_columns <- analysis_exposure_input_columns()

  diagnostics <- list(
    rows = nrow(data),
    matched_strata = data.table::uniqueN(
      data$global_id
    ),
    case_days = sum(
      data$case == 1L,
      na.rm = TRUE
    ),
    strata_without_one_case_day = sum(
      stratum_summary$case_days != 1L
    ),
    complete_exposure_rows = sum(
      stats::complete.cases(
        data[
          ,
          ..exposure_columns
        ]
      )
    )
  )

  if ("time" %in% names(data)) {
    diagnostics$first_matched_date <- min(
      data$time,
      na.rm = TRUE
    )

    diagnostics$last_matched_date <- max(
      data$time,
      na.rm = TRUE
    )

    case_dates <- data[
      case == 1L & !is.na(time),
      time
    ]

    diagnostics$first_case_date <- min(
      case_dates
    )

    diagnostics$last_case_date <- max(
      case_dates
    )
  }

  diagnostics
}

prepare_analysis_data <- function(
  config_path = file.path(
    "config",
    "config.R"
  ),
  modules = c(
    "descriptive",
    "linear",
    "bdlnm",
    "spatial",
    "subgroup"
  ),
  max_lag = ANALYSIS_MAX_LAG
) {
  config <- load_analysis_config(
    config_path
  )

  data <- read_restricted_analysis_data(
    data_file = config$data$analysis_file,
    modules = modules,
    max_lag = max_lag
  )

  data <- standardise_analysis_variables(
    data
  )

  validate_analysis_schema(
    data,
    modules = modules,
    max_lag = max_lag
  )

  data <- derive_time_variables(
    data
  )

  data <- derive_moving_average_exposures(
    data,
    max_lag = max_lag
  )

  diagnostics <- summarise_analysis_data(
    data
  )

  list(
    data = data,
    diagnostics = diagnostics
  )
}