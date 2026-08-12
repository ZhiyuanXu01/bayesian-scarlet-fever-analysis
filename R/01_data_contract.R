# ==============================================================================
# Analysis data contract
#
# This repository starts from a restricted analytical dataset in which
# environmental exposures have already been linked to case and referent days.
# Individual-level records must never be committed to this repository.
# ==============================================================================

ANALYSIS_MAX_LAG <- 10L

ANALYSIS_EXPOSURES <- data.frame(
  variable = c(
    "t2m", "rh", "ssr", "uv", "sp", "tp",
    "PM25", "PM10", "O3"
  ),
  label = c(
    "Temperature",
    "Relative humidity",
    "Solar radiation",
    "Wind speed",
    "Surface pressure",
    "Total precipitation",
    "PM2.5",
    "PM10",
    "O3"
  ),
  unit = c(
    "deg C",
    "%",
    "MJ/m2",
    "m/s",
    "kPa",
    "mm",
    "ug/m3",
    "ug/m3",
    "ug/m3"
  ),
  source_group = c(
    rep("ERA5-Land-derived meteorology", 6),
    rep("CHAP air-pollution estimates", 3)
  ),
  stringsAsFactors = FALSE
)

make_lag_columns <- function(
  variable,
  max_lag = ANALYSIS_MAX_LAG
) {
  stopifnot(
    length(variable) == 1L,
    length(max_lag) == 1L,
    max_lag >= 1L
  )

  c(
    variable,
    paste0(variable, "_lag", seq_len(max_lag))
  )
}

make_moving_average_columns <- function(
  variable,
  max_lag = ANALYSIS_MAX_LAG
) {
  stopifnot(
    length(variable) == 1L,
    length(max_lag) == 1L,
    max_lag >= 1L
  )

  paste0(
    variable,
    "_mavg",
    seq_len(max_lag)
  )
}

analysis_exposure_input_columns <- function(
  max_lag = ANALYSIS_MAX_LAG
) {
  unlist(
    lapply(
      ANALYSIS_EXPOSURES$variable,
      make_lag_columns,
      max_lag = max_lag
    ),
    use.names = FALSE
  )
}

analysis_moving_average_columns <- function(
  max_lag = ANALYSIS_MAX_LAG
) {
  unlist(
    lapply(
      ANALYSIS_EXPOSURES$variable,
      make_moving_average_columns,
      max_lag = max_lag
    ),
    use.names = FALSE
  )
}

analysis_required_columns <- function(
  modules = c("linear", "bdlnm"),
  max_lag = ANALYSIS_MAX_LAG
) {
  allowed_modules <- c(
    "descriptive",
    "linear",
    "bdlnm",
    "spatial",
    "subgroup"
  )

  modules <- match.arg(
    modules,
    allowed_modules,
    several.ok = TRUE
  )

  required <- character()

  model_modules <- c(
    "linear",
    "bdlnm",
    "spatial",
    "subgroup"
  )

  if (any(modules %in% model_modules)) {
    required <- c(
      required,
      "case",
      "global_id",
      "holiday",
      "outbreak_effect",
      analysis_exposure_input_columns(max_lag)
    )
  }

  if ("descriptive" %in% modules) {
    required <- c(
      required,
      "case",
      "time",
      "province",
      "age",
      "gender"
    )
  }

  if ("spatial" %in% modules) {
    required <- c(
      required,
      "longitude",
      "latitude"
    )
  }

  if ("subgroup" %in% modules) {
    required <- c(
      required,
      "time",
      "age",
      "gender"
    )
  }

  unique(required)
}

make_core_dictionary <- function() {
  data.frame(
    variable = c(
      "case",
      "global_id",
      "time",
      "holiday",
      "outbreak_effect",
      "province",
      "age",
      "gender",
      "longitude",
      "latitude",
      "county_id",
      "year",
      "month",
      "season"
    ),
    storage = c(
      rep("input", 10),
      rep("derived", 4)
    ),
    data_type = c(
      "integer",
      "character",
      "Date",
      "integer",
      "integer",
      "character",
      "numeric",
      "character",
      "numeric",
      "numeric",
      "integer",
      "integer",
      "integer",
      "character"
    ),
    unit = c(
      "0/1",
      NA,
      "YYYY-MM-DD",
      "0/1",
      "count",
      NA,
      "years",
      NA,
      "decimal degrees",
      "decimal degrees",
      NA,
      "year",
      "1-12",
      NA
    ),
    source_group = c(
      rep("Matched health records", 10),
      "Spatial join",
      rep("Derived from time", 3)
    ),
    role = c(
      "outcome",
      "matched-set identifier",
      "calendar time",
      "covariate",
      "covariate",
      "geographic grouping",
      "subgroup variable",
      "subgroup variable",
      "spatial coordinate",
      "spatial coordinate",
      "spatial index",
      "temporal grouping",
      "temporal grouping",
      "temporal grouping"
    ),
    required_for = c(
      "All modules",
      "All model modules",
      "Descriptive and subgroup modules",
      "All model modules",
      "All model modules",
      "Descriptive module",
      "Descriptive and subgroup modules",
      "Descriptive and subgroup modules",
      "Spatial module",
      "Spatial module",
      "Spatial model after geographic assignment",
      "Derived temporal analyses",
      "Derived temporal analyses",
      "Derived temporal analyses"
    ),
    description = c(
      "Indicator equal to 1 on a case day and 0 on a matched referent day.",
      paste(
        "Pseudonymous identifier for a matched case-control stratum;",
        "repeated across its case and referent days."
      ),
      "Calendar date of the case or matched referent day.",
      "Public-holiday indicator used as an adjustment covariate.",
      "Precomputed local outbreak adjustment covariate used in the analysis.",
      "Province containing the record.",
      "Age at disease onset.",
      "Recorded sex category, coded M or F in the analytical dataset.",
      "Restricted record longitude in WGS84.",
      "Restricted record latitude in WGS84.",
      "Sequential county index created after geographic assignment.",
      "Calendar year derived from time.",
      "Calendar month derived from time.",
      "Season derived from calendar month."
    ),
    privacy = "Restricted at row level",
    stringsAsFactors = FALSE
  )
}

make_exposure_dictionary <- function(
  max_lag = ANALYSIS_MAX_LAG
) {
  dictionaries <- lapply(
    seq_len(nrow(ANALYSIS_EXPOSURES)),
    function(i) {
      definition <- ANALYSIS_EXPOSURES[i, ]

      current <- data.frame(
        variable = definition$variable,
        storage = "input",
        data_type = "numeric",
        unit = definition$unit,
        source_group = definition$source_group,
        role = "current-day exposure",
        required_for = "Environmental model modules",
        description = paste(
          definition$label,
          "on the case or referent day (lag 0)."
        ),
        privacy = "Restricted at row level",
        stringsAsFactors = FALSE
      )

      lag_ids <- seq_len(max_lag)

      lags <- data.frame(
        variable = paste0(
          definition$variable,
          "_lag",
          lag_ids
        ),
        storage = "input",
        data_type = "numeric",
        unit = definition$unit,
        source_group = definition$source_group,
        role = "single-day lagged exposure",
        required_for = "Environmental model modules",
        description = paste0(
          definition$label,
          " ",
          lag_ids,
          " day(s) before the case or referent day."
        ),
        privacy = "Restricted at row level",
        stringsAsFactors = FALSE
      )

      moving_averages <- data.frame(
        variable = paste0(
          definition$variable,
          "_mavg",
          lag_ids
        ),
        storage = "derived",
        data_type = "numeric",
        unit = definition$unit,
        source_group = "Derived from matched lag columns",
        role = "cumulative moving-average exposure",
        required_for = paste(
          "Linear and moving-average",
          "B-DLNM modules"
        ),
        description = paste0(
          definition$label,
          " mean from lag 0 through lag ",
          lag_ids,
          ", recomputed by the public workflow."
        ),
        privacy = "Restricted at row level",
        stringsAsFactors = FALSE
      )

      rbind(
        current,
        lags,
        moving_averages
      )
    }
  )

  do.call(rbind, dictionaries)
}

build_data_dictionary <- function(
  max_lag = ANALYSIS_MAX_LAG
) {
  dictionary <- rbind(
    make_core_dictionary(),
    make_exposure_dictionary(max_lag)
  )

  rownames(dictionary) <- NULL
  dictionary
}

validate_analysis_schema <- function(
  data,
  modules = c("linear", "bdlnm"),
  max_lag = ANALYSIS_MAX_LAG
) {
  if (!is.data.frame(data)) {
    stop(
      "`data` must be a data frame or data.table.",
      call. = FALSE
    )
  }

  duplicate_columns <- unique(
    names(data)[duplicated(names(data))]
  )

  if (length(duplicate_columns) > 0L) {
    stop(
      "Duplicated columns: ",
      paste(duplicate_columns, collapse = ", "),
      call. = FALSE
    )
  }

  required <- analysis_required_columns(
    modules,
    max_lag
  )

  missing_columns <- setdiff(
    required,
    names(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  exposure_columns <- intersect(
    analysis_exposure_input_columns(max_lag),
    names(data)
  )

  non_numeric_exposures <- exposure_columns[
    !vapply(
      exposure_columns,
      function(column) {
        is.numeric(data[[column]])
      },
      logical(1)
    )
  ]

  if (length(non_numeric_exposures) > 0L) {
    stop(
      "Exposure columns must be numeric: ",
      paste(non_numeric_exposures, collapse = ", "),
      call. = FALSE
    )
  }

  if ("case" %in% names(data)) {
    case_values <- unique(
      as.character(
        stats::na.omit(data$case)
      )
    )

    if (!all(case_values %in% c("0", "1"))) {
      stop(
        "`case` must contain only 0, 1, or missing values.",
        call. = FALSE
      )
    }
  }

  if ("time" %in% names(data) &&
      !inherits(data$time, "Date")) {
    parsed_time <- suppressWarnings(
      as.Date(
        as.character(data$time)
      )
    )

    invalid_time <- is.na(parsed_time) &
      !is.na(data$time)

    if (any(invalid_time)) {
      stop(
        "`time` contains values that cannot be parsed as dates.",
        call. = FALSE
      )
    }
  }

  if ("global_id" %in% names(data)) {
    missing_id <- is.na(data$global_id) |
      !nzchar(as.character(data$global_id))

    if (any(missing_id)) {
      warning(
        "`global_id` contains missing or empty values.",
        call. = FALSE
      )
    }
  }

  message(
    "Schema validation passed for modules: ",
    paste(modules, collapse = ", "),
    "."
  )

  invisible(
    list(
      modules = modules,
      required_columns = required,
      exposure_input_columns = exposure_columns
    )
  )
}

write_data_dictionary <- function(
  path = file.path(
    "data",
    "data_dictionary.csv"
  )
) {
  dictionary <- build_data_dictionary()

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::write.csv(
    dictionary,
    path,
    row.names = FALSE,
    na = ""
  )

  message(
    "Wrote ",
    nrow(dictionary),
    " variables to ",
    path,
    "."
  )

  invisible(dictionary)
}