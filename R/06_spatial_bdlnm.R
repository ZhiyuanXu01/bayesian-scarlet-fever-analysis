# ==============================================================================
# Exploratory spatial Bayesian distributed lag nonlinear models
#
# Spatial preparation can be validated locally. Full SB-DLNM fitting requires
# a high-memory computing environment and is not run automatically.
# ==============================================================================

source(
  file.path(
    "R",
    "05_pooled_bdlnm.R"
  )
)

CHINA_ALBERS_EQUAL_AREA <- paste(
  "+proj=aea",
  "+lat_1=25",
  "+lat_2=47",
  "+lat_0=0",
  "+lon_0=105",
  "+datum=WGS84",
  "+units=m",
  "+no_defs"
)

check_spatial_packages <- function() {
  required_packages <- c(
    "data.table",
    "dplyr",
    "sf",
    "spdep",
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
      "Missing spatial packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

read_province_boundaries <- function(
    path,
    province
) {
  if (!file.exists(path)) {
    stop(
      "Cannot find county boundaries: ",
      path,
      call. = FALSE
    )
  }
  
  boundary <- sf::st_read(
    path,
    quiet = TRUE,
    stringsAsFactors = FALSE
  )
  
  if (!"name" %in%
      names(boundary)) {
    stop(
      basename(path),
      " must contain a `name` field.",
      call. = FALSE
    )
  }
  
  if (is.na(sf::st_crs(boundary))) {
    sf::st_crs(boundary) <- 4326
  }
  
  boundary$province <- province
  boundary$county_name <-
    as.character(boundary$name)
  boundary$source_county_id <-
    seq_len(nrow(boundary))
  
  boundary[
    ,
    c(
      "province",
      "county_name",
      "source_county_id",
      "geometry"
    )
  ]
}

build_county_layer <- function(
    hubei_path,
    hunan_path
) {
  check_spatial_packages()
  
  hubei <- read_province_boundaries(
    hubei_path,
    "Hubei"
  )
  
  hunan <- read_province_boundaries(
    hunan_path,
    "Hunan"
  )
  
  counties <- dplyr::bind_rows(
    hubei,
    hunan
  )
  
  counties <- sf::st_transform(
    counties,
    CHINA_ALBERS_EQUAL_AREA
  )
  
  counties <- sf::st_make_valid(
    counties
  )
  
  counties <- sf::st_collection_extract(
    counties,
    "POLYGON",
    warn = FALSE
  )
  
  counties <- counties |>
    dplyr::group_by(
      province,
      county_name,
      source_county_id
    ) |>
    dplyr::summarise(
      geometry = sf::st_union(geometry),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      province,
      source_county_id
    ) |>
    dplyr::mutate(
      county_id = dplyr::row_number()
    ) |>
    dplyr::select(
      county_id,
      province,
      county_name,
      geometry
    )
  
  counties <- sf::st_make_valid(
    sf::st_buffer(
      counties,
      0
    )
  )
  
  counties
}

assign_counties_to_records <- function(
    data,
    counties
) {
  required_columns <- c(
    "longitude",
    "latitude"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing coordinates: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  spatial_data <- data.table::copy(
    data
  )
  
  existing_spatial_columns <- intersect(
    c(
      "county_id",
      "county_name",
      "spatial_province"
    ),
    names(spatial_data)
  )
  
  if (length(existing_spatial_columns) > 0L) {
    spatial_data[
      ,
      (existing_spatial_columns) := NULL
    ]
  }
  
  spatial_data[
    ,
    spatial_row_id := .I
  ]
  
  points <- sf::st_as_sf(
    spatial_data,
    coords = c(
      "longitude",
      "latitude"
    ),
    crs = 4326,
    remove = FALSE
  )
  
  points <- sf::st_transform(
    points,
    sf::st_crs(counties)
  )
  
  county_join_layer <- counties |>
    dplyr::rename(
      spatial_province = province
    )
  
  old_s2 <- sf::sf_use_s2()
  
  on.exit(
    sf::sf_use_s2(old_s2),
    add = TRUE
  )
  
  sf::sf_use_s2(FALSE)
  
  joined <- sf::st_join(
    points,
    county_join_layer,
    join = sf::st_intersects,
    left = TRUE
  )
  
  joined_data <- data.table::as.data.table(
    sf::st_drop_geometry(joined)
  )
  
  if (data.table::uniqueN(
    joined_data$spatial_row_id
  ) != nrow(joined_data)) {
    stop(
      paste(
        "At least one record matched multiple counties.",
        "Inspect boundary overlaps before modelling."
      ),
      call. = FALSE
    )
  }
  
  data.table::setorder(
    joined_data,
    spatial_row_id
  )
  
  diagnostics <- list(
    input_rows = nrow(spatial_data),
    matched_rows = sum(
      !is.na(joined_data$county_id)
    ),
    unmatched_rows = sum(
      is.na(joined_data$county_id)
    ),
    counties_with_records =
      data.table::uniqueN(
        joined_data[
          !is.na(county_id),
          county_id
        ]
      )
  )
  
  joined_data[
    ,
    spatial_row_id := NULL
  ]
  
  list(
    data = joined_data[
      !is.na(county_id)
    ],
    diagnostics = diagnostics
  )
}

repair_islands <- function(
    neighbours,
    counties
) {
  island_ids <- which(
    spdep::card(neighbours) == 0L
  )
  
  if (length(island_ids) == 0L) {
    return(neighbours)
  }
  
  representative_points <-
    sf::st_point_on_surface(
      sf::st_geometry(counties)
    )
  
  all_ids <- seq_along(
    neighbours
  )
  
  for (island_id in island_ids) {
    candidates <- setdiff(
      all_ids,
      island_id
    )
    
    nearest_position <-
      sf::st_nearest_feature(
        representative_points[island_id],
        representative_points[candidates]
      )
    
    neighbour_id <-
      candidates[nearest_position]
    
    neighbours[[island_id]] <-
      as.integer(neighbour_id)
    
    neighbours[[neighbour_id]] <- sort(
      unique(
        as.integer(
          c(
            neighbours[[neighbour_id]],
            island_id
          )
        )
      )
    )
  }
  
  attr(
    neighbours,
    "region.id"
  ) <- as.character(
    seq_along(neighbours)
  )
  
  class(neighbours) <- "nb"
  
  neighbours
}

build_county_adjacency <- function(
    counties,
    graph_file = file.path(
      "results",
      "private",
      "county_adjacency.graph"
    )
) {
  dir.create(
    dirname(graph_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ordered_counties <- counties |>
    dplyr::arrange(
      county_id
    ) |>
    sf::st_make_valid()
  
  old_s2 <- sf::sf_use_s2()
  
  on.exit(
    sf::sf_use_s2(old_s2),
    add = TRUE
  )
  
  sf::sf_use_s2(FALSE)
  
  neighbours <- spdep::poly2nb(
    ordered_counties,
    queen = TRUE,
    snap = 100
  )
  
  islands_before <- which(
    spdep::card(neighbours) == 0L
  )
  
  neighbours <- repair_islands(
    neighbours,
    ordered_counties
  )
  
  islands_after <- which(
    spdep::card(neighbours) == 0L
  )
  
  spdep::nb2INLA(
    file = graph_file,
    nb = neighbours
  )
  
  graph <- INLA::inla.read.graph(
    graph_file
  )
  
  diagnostics <- list(
    counties = nrow(ordered_counties),
    islands_before_repair =
      length(islands_before),
    islands_after_repair =
      length(islands_after),
    minimum_neighbours = min(
      spdep::card(neighbours)
    ),
    median_neighbours = stats::median(
      spdep::card(neighbours)
    ),
    maximum_neighbours = max(
      spdep::card(neighbours)
    )
  )
  
  list(
    neighbours = neighbours,
    graph = graph,
    diagnostics = diagnostics,
    graph_file = graph_file
  )
}

prepare_spatial_objects <- function(
    config_path = file.path(
      "config",
      "config.R"
    )
) {
  check_spatial_packages()
  
  config <- load_analysis_config(
    config_path
  )
  
  prepared <- prepare_analysis_data(
    config_path = config_path,
    modules = c(
      "bdlnm",
      "spatial"
    )
  )
  
  counties <- build_county_layer(
    hubei_path =
      config$data$hubei_boundary,
    hunan_path =
      config$data$hunan_boundary
  )
  
  assigned <- assign_counties_to_records(
    prepared$data,
    counties
  )
  
  adjacency <- build_county_adjacency(
    counties
  )
  
  diagnostics <- list(
    boundary_counties =
      nrow(counties),
    input_rows =
      assigned$diagnostics$input_rows,
    matched_rows =
      assigned$diagnostics$matched_rows,
    unmatched_rows =
      assigned$diagnostics$unmatched_rows,
    counties_with_records =
      assigned$diagnostics$
      counties_with_records,
    islands_before_repair =
      adjacency$diagnostics$
      islands_before_repair,
    islands_after_repair =
      adjacency$diagnostics$
      islands_after_repair,
    minimum_neighbours =
      adjacency$diagnostics$
      minimum_neighbours,
    median_neighbours =
      adjacency$diagnostics$
      median_neighbours,
    maximum_neighbours =
      adjacency$diagnostics$
      maximum_neighbours
  )
  
  list(
    data = assigned$data,
    counties = counties,
    neighbours = adjacency$neighbours,
    graph = adjacency$graph,
    diagnostics = diagnostics
  )
}

build_spatial_bdlnm_formula <- function(
    crossbasis_names
) {
  fixed_terms <- paste(
    c(
      crossbasis_names,
      "holiday",
      "outbreak_effect"
    ),
    collapse = " + "
  )
  
  spatial_terms <- vapply(
    seq_along(crossbasis_names),
    function(index) {
      paste0(
        "f(",
        "area_cb_",
        index,
        ", ",
        crossbasis_names[index],
        ", model = 'bym2'",
        ", graph = county_graph",
        ", scale.model = TRUE",
        ", constr = TRUE",
        ", hyper = hyper_spatial",
        ")"
      )
    },
    character(1)
  )
  
  model_formula <- stats::as.formula(
    paste0(
      "case ~ 1 + ",
      fixed_terms,
      " + f(",
      "gid_inla",
      ", model = 'iid'",
      ", hyper = hyper_stratum",
      ") + ",
      paste(
        spatial_terms,
        collapse = " + "
      )
    )
  )
  
  environment(model_formula) <-
    parent.frame()
  
  model_formula
}

inspect_spatial_bdlnm_design <- function(
    spatial_objects,
    variable = "t2m",
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
    )
) {
  exposure_type <- match.arg(
    exposure_type
  )
  
  exposure_columns <-
    make_bdlnm_exposure_columns(
      variable,
      exposure_type,
      max_lag
    )
  
  required_columns <- c(
    "case",
    "global_id",
    "county_id",
    "holiday",
    "outbreak_effect",
    exposure_columns
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(spatial_objects$data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing spatial model columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  model_data <- spatial_objects$data[
    stats::complete.cases(
      spatial_objects$data[
        ,
        ..required_columns
      ]
    ),
    ..required_columns
  ]
  
  exposure_history <- as.matrix(
    model_data[
      ,
      ..exposure_columns
    ]
  )
  
  storage.mode(exposure_history) <-
    "numeric"
  
  basis <- build_crossbasis(
    exposure_history,
    max_lag = max_lag,
    lag_df = lag_df,
    exposure_knot_probabilities =
      exposure_knot_probabilities
  )
  
  data.table::data.table(
    variable = variable,
    exposure_type = exposure_type,
    rows = nrow(model_data),
    matched_strata =
      data.table::uniqueN(
        model_data$global_id
      ),
    counties_with_records =
      data.table::uniqueN(
        model_data$county_id
      ),
    graph_counties =
      spatial_objects$graph$n,
    crossbasis_columns =
      length(basis$names),
    BYM2_components =
      length(basis$names),
    maximum_lag = max_lag,
    lag_df = lag_df
  )
}

fit_spatial_bdlnm <- function(
    spatial_objects,
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
    pc_prior_stratum = c(
      3,
      0.01
    ),
    inla_num_threads = "4:1",
    compute_config = FALSE,
    compute_cpo = FALSE,
    compute_fitted = FALSE,
    return_model = FALSE
) {
  check_spatial_packages()
  
  exposure_type <- match.arg(
    exposure_type
  )
  
  spatial_data <- spatial_objects$data
  county_graph <- spatial_objects$graph
  
  exposure_columns <-
    make_bdlnm_exposure_columns(
      variable,
      exposure_type,
      max_lag
    )
  
  required_columns <- c(
    "case",
    "global_id",
    "county_id",
    "holiday",
    "outbreak_effect",
    exposure_columns
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(spatial_data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing SB-DLNM columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  model_data <- data.table::copy(
    spatial_data[
      stats::complete.cases(
        spatial_data[
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
  
  model_data[
    ,
    county_id := as.integer(
      county_id
    )
  ]
  
  if (max(model_data$county_id) >
      county_graph$n) {
    stop(
      "County identifiers exceed the INLA graph size.",
      call. = FALSE
    )
  }
  
  exposure_history <- as.matrix(
    model_data[
      ,
      ..exposure_columns
    ]
  )
  
  storage.mode(exposure_history) <-
    "numeric"
  
  basis <- build_crossbasis(
    exposure_history,
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
  
  for (index in seq_along(
    basis$names
  )) {
    model_data[[
        paste0(
          "area_cb_",
          index
        )
      ]] <- model_data$county_id
  }
  
  hyper_stratum <- list(
    prec = list(
      prior = "pc.prec",
      param = pc_prior_stratum
    )
  )
  
  hyper_spatial <- list(
    prec = list(
      prior = "pc.prec",
      param = c(
        1,
        0.01
      )
    ),
    phi = list(
      prior = "pc",
      param = c(
        0.5,
        2 / 3
      )
    )
  )
  
  model_formula <-
    build_spatial_bdlnm_formula(
      basis$names
    )
  
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
      compute = compute_fitted
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
    n_rows = nrow(model_data),
    n_strata =
      data.table::uniqueN(
        model_data$gid_inla
      ),
    n_counties =
      data.table::uniqueN(
        model_data$county_id
      ),
    graph_counties =
      county_graph$n,
    crossbasis_columns =
      length(basis$names),
    BYM2_components =
      length(basis$names),
    compute_config = compute_config,
    compute_cpo = compute_cpo,
    compute_fitted = compute_fitted,
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
    crossbasis_names = basis$names,
    exposure_basis =
      basis$exposure_basis,
    lag_basis = basis$lag_basis,
    exposure_knots =
      basis$exposure_knots,
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
    
    if (isTRUE(compute_fitted)) {
      output$county_id <-
        model_data$county_id
    }
  } else {
    rm(model)
    
    invisible(
      gc()
    )
  }
  
  output
}

run_spatial_bdlnm_model <- function(
    variable,
    exposure_type = c(
      "lag",
      "mavg"
    ),
    config_path = file.path(
      "config",
      "config.R"
    ),
    compute_config = FALSE,
    compute_cpo = FALSE,
    compute_fitted = FALSE,
    return_model = FALSE
) {
  exposure_type <- match.arg(
    exposure_type
  )
  
  config <- load_analysis_config(
    config_path
  )
  
  spatial_objects <-
    prepare_spatial_objects(
      config_path
    )
  
  fit_spatial_bdlnm(
    spatial_objects =
      spatial_objects,
    variable = variable,
    exposure_type = exposure_type,
    inla_num_threads =
      config$computation$
      inla_num_threads,
    compute_config = compute_config,
    compute_cpo = compute_cpo,
    compute_fitted = compute_fitted,
    return_model = return_model
  )
}