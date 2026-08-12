# ==============================================================================
# Privacy-safe descriptive analysis
#
# Only aggregate tables and figures are written. Case characteristics are
# calculated from case-day rows rather than repeated case-control records.
# ==============================================================================

source(
  file.path(
    "R",
    "02_prepare_analysis_data.R"
  )
)

check_descriptive_packages <- function() {
  required_packages <- c(
    "data.table",
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
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

summarise_case_characteristic <- function(
  data,
  column,
  characteristic
) {
  values <- as.character(
    data[[column]]
  )

  values[
    is.na(values) |
      !nzchar(values)
  ] <- "Missing"

  summary <- data.table::data.table(
    characteristic = characteristic,
    category = values
  )[
    ,
    .(
      n = .N
    ),
    by = .(
      characteristic,
      category
    )
  ]

  summary[
    ,
    percent := round(
      100 * n / sum(n),
      1
    ),
    by = characteristic
  ]

  summary[]
}

build_case_characteristics <- function(data) {
  case_data <- data.table::copy(
    data[
      case == 1L
    ]
  )

  case_data[
    ,
    age_group := as.character(
      cut(
        age,
        breaks = c(
          -Inf,
          2,
          6,
          9,
          14,
          Inf
        ),
        labels = c(
          "0-2",
          "3-6",
          "7-9",
          "10-14",
          "15+"
        ),
        right = TRUE
      )
    )
  ]

  case_data[
    ,
    gender_label := data.table::fcase(
      gender == "M",
      "Male",

      gender == "F",
      "Female",

      default = "Missing or other"
    )
  ]

  summaries <- list(
    summarise_case_characteristic(
      case_data,
      "province",
      "Province"
    ),
    summarise_case_characteristic(
      case_data,
      "gender_label",
      "Sex"
    ),
    summarise_case_characteristic(
      case_data,
      "age_group",
      "Age group (years)"
    ),
    summarise_case_characteristic(
      case_data,
      "season",
      "Season"
    )
  )

  result <- data.table::rbindlist(
    summaries,
    use.names = TRUE
  )

  result[
    ,
    display := sprintf(
      "%s (%.1f%%)",
      format(
        n,
        big.mark = ",",
        scientific = FALSE
      ),
      percent
    )
  ]

  result[]
}

build_exposure_summary <- function(data) {
  summaries <- lapply(
    seq_len(nrow(ANALYSIS_EXPOSURES)),
    function(i) {
      definition <- ANALYSIS_EXPOSURES[i, ]
      values <- data[[
          definition$variable
      ]]

      observed <- values[
        !is.na(values)
      ]

      data.table::data.table(
        variable = definition$variable,
        exposure = definition$label,
        unit = definition$unit,
        n_rows = length(values),
        n_observed = length(observed),
        n_missing = sum(is.na(values)),
        mean = mean(observed),
        sd = stats::sd(observed),
        median = stats::median(observed),
        p25 = as.numeric(
          stats::quantile(
            observed,
            0.25,
            names = FALSE
          )
        ),
        p75 = as.numeric(
          stats::quantile(
            observed,
            0.75,
            names = FALSE
          )
        ),
        minimum = min(observed),
        maximum = max(observed)
      )
    }
  )

  data.table::rbindlist(
    summaries,
    use.names = TRUE
  )
}

build_temporal_counts <- function(data) {
  case_data <- data[
    case == 1L &
      !is.na(time)
  ]

  annual_counts <- case_data[
    ,
    .(
      case_count = .N
    ),
    by = .(
      year
    )
  ][
    order(year)
  ]

  monthly_counts <- case_data[
    ,
    .(
      case_count = .N
    ),
    by = .(
      month_date = as.Date(
        format(
          time,
          "%Y-%m-01"
        )
      )
    )
  ][
    order(month_date)
  ]

  calendar <- data.table::data.table(
    month_date = seq.Date(
      min(monthly_counts$month_date),
      max(monthly_counts$month_date),
      by = "month"
    )
  )

  monthly_counts <- merge(
    calendar,
    monthly_counts,
    by = "month_date",
    all.x = TRUE,
    sort = TRUE
  )

  monthly_counts[
    is.na(case_count),
    case_count := 0L
  ]

  monthly_counts[
    ,
    rolling_12_month_mean :=
      data.table::frollmean(
        case_count,
        n = 12L,
        align = "right",
        fill = NA_real_
      )
  ]

  list(
    annual = annual_counts,
    monthly = monthly_counts
  )
}

build_analysis_overview <- function(
  data,
  diagnostics
) {
  data.table::data.table(
    metric = c(
      "Analytical rows",
      "Matched strata",
      "Case days",
      "First case date",
      "Last case date"
    ),
    value = c(
      as.character(
        diagnostics$rows
      ),
      as.character(
        diagnostics$matched_strata
      ),
      as.character(
        diagnostics$case_days
      ),
      as.character(
        diagnostics$first_case_date
      ),
      as.character(
        diagnostics$last_case_date
      )
    )
  )
}

plot_monthly_case_counts <- function(
  monthly_counts
) {
  ggplot2::ggplot(
    monthly_counts,
    ggplot2::aes(
      x = month_date,
      y = case_count
    )
  ) +
    ggplot2::geom_line(
      colour = "#7A9E9F",
      linewidth = 0.45,
      alpha = 0.8
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = rolling_12_month_mean
      ),
      colour = "#C8553D",
      linewidth = 1,
      na.rm = TRUE
    ) +
    ggplot2::scale_x_date(
      date_breaks = "2 years",
      date_labels = "%Y",
      expand = ggplot2::expansion(
        mult = c(0.01, 0.02)
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(
        mult = c(0, 0.05)
      )
    ) +
    ggplot2::labs(
      title = "Monthly scarlet fever case counts",
      subtitle = paste(
        "Monthly observations and trailing",
        "12-month mean"
      ),
      x = NULL,
      y = "Reported cases",
      caption = paste(
        "Blue: monthly cases;",
        "orange: trailing 12-month mean."
      )
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.caption = ggplot2::element_text(
        colour = "grey40"
      )
    )
}

run_descriptive_analysis <- function(
  config_path = file.path(
    "config",
    "config.R"
  )
) {
  check_descriptive_packages()

  config <- load_analysis_config(
    config_path
  )

  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = c(
      "descriptive",
      "linear"
    )
  )

  data <- prepared$data

  overview <- build_analysis_overview(
    data,
    prepared$diagnostics
  )

  case_characteristics <- build_case_characteristics(
    data
  )

  exposure_summary <- build_exposure_summary(
    data
  )

  temporal_counts <- build_temporal_counts(
    data
  )

  results_dir <- config$output$results_dir
  figures_dir <- config$output$figures_dir

  dir.create(
    results_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    figures_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  data.table::fwrite(
    overview,
    file.path(
      results_dir,
      "analysis_overview.csv"
    )
  )

  data.table::fwrite(
    case_characteristics,
    file.path(
      results_dir,
      "case_characteristics.csv"
    )
  )

  data.table::fwrite(
    exposure_summary,
    file.path(
      results_dir,
      "exposure_summary.csv"
    )
  )

  data.table::fwrite(
    temporal_counts$annual,
    file.path(
      results_dir,
      "annual_case_counts.csv"
    )
  )

  data.table::fwrite(
    temporal_counts$monthly,
    file.path(
      results_dir,
      "monthly_case_counts.csv"
    )
  )

  monthly_plot <- plot_monthly_case_counts(
    temporal_counts$monthly
  )

  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "monthly_case_counts.png"
    ),
    plot = monthly_plot,
    width = 9,
    height = 5,
    dpi = 300,
    bg = "white"
  )

  message(
    "Descriptive analysis completed."
  )

  invisible(
    list(
      overview = overview,
      case_characteristics = case_characteristics,
      exposure_summary = exposure_summary,
      annual_counts = temporal_counts$annual,
      monthly_counts = temporal_counts$monthly,
      monthly_plot = monthly_plot
    )
  )
}