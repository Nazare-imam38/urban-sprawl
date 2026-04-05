#' Phase 6 — Spatial drivers: growth raster, OSM roads, distances, optional WorldPop/DEM,
#' NDVI loss, predictor stack, stratified samples, glm + lagsarlm + GWR, correlations, plots.
#'
#' **Speed:** Driver rasters use a **coarsened** LULC grid (`phase6.coarsen_fact`, default 8) so distance
#' layers are tractable; road distance uses rasterize + grid distance (not per-cell geometry). OSM roads
#' cache to `data/processed/phase6_osm_roads.gpkg`. Use `coarsen_fact: 1` only if you need native resolution.
#'
#' **Memory-friendly workflow:** run in steps (writes GeoTIFFs, then frees RAM):
#' 1. `run_phase6_build()` — build stack, write `spatial_predictors.tif` + `phase6_growth_*.tif`
#' 2. `run_phase6_sample()` — sample + `regression_samples.csv` (reads rasters from disk)
#' 3. `run_phase6_models()` — GLM / Moran / optional GWR from CSV (`phase6.run_gwr` in yml)
#'
#' Or `run_phase6(step = "all")` runs 1→2→3 with `gc()` between steps.
#'
#' Optional packages: **spdep**, **spatialreg**, **sp**; **GWmodel** only if `run_gwr = TRUE`.
#' Working directory = project root.

source("R/00_setup.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

get_phase6_cfg <- function(cfg) {
  p <- cfg$phase6
  if (is.null(p)) stop("config: add `phase6` block (growth years, sample sizes, optional paths).", call. = FALSE)
  p
}

phase6_dirs <- function(cfg) {
  list(
    proc = file.path(project_root(), cfg$paths$processed_dir),
    outd = file.path(project_root(), cfg$paths$outputs_dir)
  )
}

phase6_predictors_path <- function(cfg) {
  file.path(phase6_dirs(cfg)$proc, "spatial_predictors.tif")
}

phase6_growth_path <- function(cfg) {
  p6 <- get_phase6_cfg(cfg)
  y0 <- as.integer(p6$growth_year_t0 %||% 2000)
  y1 <- as.integer(p6$growth_year_t1 %||% 2024)
  file.path(phase6_dirs(cfg)$proc, sprintf("phase6_growth_%d_%d.tif", y0, y1))
}

phase6_samples_path <- function(cfg) {
  file.path(phase6_dirs(cfg)$proc, "regression_samples.csv")
}

path_lulc <- function(year, cfg) {
  file.path(cfg$paths$processed_dir, sprintf("lulc_%d.tif", year))
}

path_features <- function(year, cfg) {
  file.path(cfg$paths$processed_dir, sprintf("features_%d.tif", year))
}

#' Coarse grid for Phase 6 (distances scale ~ with n pixels; 10 m LULC → hours otherwise).
phase6_working_template <- function(cfg, lg = function(...) NULL) {
  p6 <- get_phase6_cfg(cfg)
  y1 <- as.integer(p6$growth_year_t1 %||% 2024)
  r <- terra::rast(path_lulc(y1, cfg))
  res_m <- suppressWarnings(as.numeric(p6$working_resolution_m))
  fact <- as.integer(p6$coarsen_fact %||% 8L)
  fact <- max(1L, fact)
  if (length(res_m) && is.finite(res_m[1]) && res_m[1] > 0) {
    cur <- terra::res(r)[1]
    fact <- max(1L, as.integer(round(res_m[1] / cur)))
    lg("Phase 6: working_resolution_m=", res_m[1], " → aggregate fact=", fact, ".")
  }
  if (fact > 1L) {
    lg("Phase 6: aggregating LULC (coarsen_fact=", fact, " → ~", fact * fact, "× fewer cells) …")
    r <- terra::aggregate(r, fact = fact, fun = "modal", na.rm = TRUE)
  }
  nc <- ncol(r)
  nr <- nrow(r)
  lg("Phase 6: working grid ", nc, " × ", nr, " cells; res ≈ ", format(mean(terra::res(r)), digits = 5), " m.")
  r
}

study_aoi_polygon <- function(cfg) {
  sa <- cfg$study_area
  bb <- sf::st_bbox(
    c(xmin = sa$xmin, ymin = sa$ymin, xmax = sa$xmax, ymax = sa$ymax),
    crs = sa$crs_wgs84
  )
  sf::st_as_sf(sf::st_as_sfc(bb))
}

city_centre_sf <- function(cfg) {
  s <- cfg$sprawl
  if (is.null(s$centre_lon)) stop("config$sprawl centre_lon/lat required.", call. = FALSE)
  p <- sf::st_point(c(s$centre_lon, s$centre_lat))
  sf::st_sf(geometry = sf::st_sfc(p, crs = cfg$study_area$crs_wgs84))
}

# --------------------------------------------------------------------------- #
# Growth Y and distances
# --------------------------------------------------------------------------- #

growth_binary_rast <- function(year_t0, year_t1, cfg, template) {
  r0 <- terra::round(terra::rast(path_lulc(year_t0, cfg)))
  r1 <- terra::round(terra::rast(path_lulc(year_t1, cfg)))
  r0 <- terra::resample(r0, template, method = "near")
  r1 <- terra::resample(r1, template, method = "near")
  y <- terra::ifel((r0 != 1L) & (r1 == 1L), 1L, 0L)
  names(y) <- "growth"
  y
}

dist_to_urban_class <- function(lulc_on_template, template) {
  lu <- terra::round(lulc_on_template)
  x <- terra::ifel(lu == 1L, 1L, NA_integer_)
  d <- terra::distance(x)
  terra::ifel(lu == 1L, 0, d)
}

fetch_osm_roads_lines <- function(cfg, progress_log = NULL) {
  lg <- if (is.function(progress_log)) progress_log else function(...) NULL
  crs_p <- cfg$study_area$crs_projected
  cache_gpkg <- file.path(project_root(), cfg$paths$processed_dir, "phase6_osm_roads.gpkg")
  p6 <- cfg$phase6 %||% list()
  refresh <- isTRUE(p6$refresh_osm_cache %||% FALSE)
  oc <- p6$osm_roads_cache
  if (!is.null(oc) && nzchar(as.character(oc)[1])) {
    cache_gpkg <- file.path(project_root(), as.character(oc)[1])
  }
  if (!isTRUE(refresh) && file.exists(cache_gpkg)) {
    lg("Phase 6: OSM roads from cache ", cache_gpkg, " (delete file or set phase6.refresh_osm_cache: true to refresh).")
    ln <- sf::st_read(cache_gpkg, quiet = TRUE)
    return(sf::st_transform(ln, terra::crs(crs_p)))
  }

  aoi <- study_aoi_polygon(cfg)
  bb <- sf::st_bbox(sf::st_transform(aoi, 4326))
  highways <- c("motorway", "trunk", "primary", "secondary", "tertiary", "residential")
  overpass_urls <- c(
    "https://overpass-api.de/api/interpreter",
    "https://api.openstreetmap.fr/oapi/interpreter",
    "https://overpass.osm.vi-di.fr/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter"
  )
  last_err <- NULL
  for (u in overpass_urls) {
    try(osmdata::set_overpass_url(u), silent = TRUE)
    q <- osmdata::opq(bbox = bb, timeout = 180L) |>
      osmdata::add_osm_feature(key = "highway", value = highways)
    od <- tryCatch(osmdata::osmdata_sf(q), error = function(e) e)
    if (inherits(od, "error")) {
      last_err <- od
      next
    }
    ln <- od$osm_lines
    if (!is.null(ln) && nrow(ln) > 0) {
      ln_p <- sf::st_transform(ln, terra::crs(crs_p))
      dir.create(dirname(cache_gpkg), showWarnings = FALSE, recursive = TRUE)
      tryCatch(
        sf::st_write(ln_p, cache_gpkg, delete_dsn = TRUE, quiet = TRUE),
        error = function(e) warning("Phase 6: could not cache OSM roads: ", conditionMessage(e), call. = FALSE)
      )
      return(ln_p)
    }
    last_err <- simpleError(paste0("No OSM road lines from server: ", u))
  }
  if (is.null(last_err)) {
    last_err <- simpleError("Overpass request failed (no details)")
  }
  warning(
    "OSM roads download failed (Overpass). Continuing Phase 6 without roads. Last error: ",
    conditionMessage(last_err)
  )
  NULL
}

dist_to_lines_rast <- function(template, lines_sf) {
  if (is.null(lines_sf) || nrow(lines_sf) == 0) {
    out <- template
    terra::values(out) <- 5e4
    names(out) <- "dist_roads"
    return(out)
  }
  v <- terra::vect(lines_sf)
  # Rasterize lines then grid distance — avoids terra::distance(raster, lines), which is O(n_cells × geometry).
  rw <- tryCatch(
    suppressWarnings(terra::rasterize(v, template, touches = TRUE)),
    error = function(e) NULL
  )
  vals <- if (!is.null(rw)) terra::values(rw, mat = FALSE) else NULL
  if (is.null(rw) || !length(vals) || !any(!is.na(vals))) {
    warning("Phase 6: roads did not rasterize; using flat distance fill.", call. = FALSE)
    out <- template
    terra::values(out) <- 5e4
    names(out) <- "dist_roads"
    return(out)
  }
  hit <- !is.na(rw)
  x <- terra::ifel(hit, 1L, NA_integer_)
  d <- tryCatch(terra::distance(x), error = function(e) NULL)
  if (is.null(d)) {
    out <- template
    terra::values(out) <- 5e4
    names(out) <- "dist_roads"
    return(out)
  }
  d <- terra::ifel(hit, 0, d)
  d <- terra::ifel(is.infinite(d), NA, d)
  mx <- suppressWarnings(max(terra::values(d, mat = FALSE), na.rm = TRUE))
  if (!is.finite(mx) || mx <= 0) mx <- 5e4
  d <- terra::ifel(is.na(d), mx, d)
  names(d) <- "dist_roads"
  d
}

dist_to_point_rast <- function(template, pt_sf_proj) {
  xy_pt <- terra::crds(terra::vect(pt_sf_proj))[1L, , drop = TRUE]
  xr <- terra::init(template, "x")
  yr <- terra::init(template, "y")
  d <- sqrt((xr - xy_pt[1])^2 + (yr - xy_pt[2])^2)
  names(d) <- "dist_centre"
  d
}

# --------------------------------------------------------------------------- #
# Optional WorldPop + DEM + NDVI loss
# --------------------------------------------------------------------------- #

load_resample_raster <- function(path, template, method = "bilinear") {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  r <- terra::rast(path)
  terra::resample(r, template, method = method)
}

pop_density_per_km2 <- function(pop_rast, template) {
  if (is.null(pop_rast)) return(NULL)
  pr <- terra::resample(pop_rast, template, method = "bilinear")
  px_km2 <- prod(terra::res(template)) / 1e6
  pr / px_km2
}

dem_and_slope <- function(dem_path, template) {
  if (is.null(dem_path) || !nzchar(dem_path) || !file.exists(dem_path)) {
    return(list(elev = NULL, slope = NULL))
  }
  dem <- terra::rast(dem_path)
  dem <- terra::resample(dem, template, method = "bilinear")
  sl <- terra::terrain(dem, v = "slope", unit = "degrees")
  list(elev = dem, slope = sl)
}

ndvi_loss_rast <- function(cfg, template, y0, y1) {
  f0 <- path_features(y0, cfg)
  f1 <- path_features(y1, cfg)
  if (!file.exists(f0) || !file.exists(f1)) {
    warning("Missing features rasters for NDVI loss.")
    return(NULL)
  }
  s0 <- terra::rast(f0)[["NDVI"]]
  s1 <- terra::rast(f1)[["NDVI"]]
  s0 <- terra::resample(s0, template, method = "bilinear")
  s1 <- terra::resample(s1, template, method = "bilinear")
  loss <- s0 - s1
  terra::ifel(loss < 0, 0, loss)
}

# --------------------------------------------------------------------------- #
# Full stack + stratified samples
# --------------------------------------------------------------------------- #

build_phase6_stack <- function(cfg, progress_log = NULL) {
  lg <- if (is.function(progress_log)) progress_log else function(...) NULL
  p6 <- get_phase6_cfg(cfg)
  y0 <- as.integer(p6$growth_year_t0 %||% 2000)
  y1 <- as.integer(p6$growth_year_t1 %||% 2024)
  lg("Phase 6 stack: build working template (coarse grid for speed) …")
  template <- phase6_working_template(cfg, lg)
  crs_t <- cfg$study_area$crs_projected

  lg("Phase 6 stack: growth binary (", y0, "→", y1, ") …")
  growth <- growth_binary_rast(y0, y1, cfg, template)

  lg("Phase 6 stack: OSM roads (network or cache) …")
  roads <- fetch_osm_roads_lines(cfg, progress_log = lg)
  lg("Phase 6 stack: distance to roads …")
  d_road <- dist_to_lines_rast(template, roads)

  lg("Phase 6 stack: distance to city centre …")
  ctr <- sf::st_transform(city_centre_sf(cfg), crs_t)
  d_ctr <- dist_to_point_rast(template, ctr)

  lg("Phase 6 stack: distance to urban (t0) …")
  l0 <- terra::resample(terra::round(terra::rast(path_lulc(y0, cfg))), template, method = "near")
  d_urb <- dist_to_urban_class(l0, template)

  lg("Phase 6 stack: optional WorldPop / DEM / NDVI loss …")
  pop <- load_resample_raster(p6$worldpop_2000_path %||% "", template)
  pop_d <- pop_density_per_km2(pop, template)

  ds <- dem_and_slope(p6$dem_path %||% "", template)
  elev <- ds$elev
  slope <- ds$slope

  ndloss <- ndvi_loss_rast(cfg, template, y0, y1)

  lg("Phase 6 stack: combining predictor layers …")
  layers <- list(
    dist_roads = d_road,
    dist_centre = d_ctr,
    dist_urban_2000 = d_urb
  )
  if (!is.null(pop_d)) layers$pop_density_2000 <- pop_d
  if (!is.null(elev)) {
    layers$elevation <- elev
    layers$slope <- slope
  }
  if (!is.null(ndloss)) layers$ndvi_loss <- ndloss

  stk <- terra::rast(layers)
  names(stk) <- names(layers)
  lg("Phase 6 stack: predictor stack ready (", terra::nlyr(stk), " layers).")
  list(growth = growth, predictors = stk, template = template)
}

stratified_sample_points <- function(growth_r, n_per = 1000L, seed = 42L) {
  set.seed(seed)
  n_per <- as.integer(n_per)
  g1 <- growth_r == 1L
  g0 <- growth_r == 0L
  n1 <- sum(terra::values(g1, mat = FALSE) == 1, na.rm = TRUE)
  n0 <- sum(terra::values(g0, mat = FALSE) == 1, na.rm = TRUE)
  s1 <- min(n_per, n1)
  s0 <- min(n_per, n0)
  if (s1 < 10L || s0 < 10L) {
    stop("Too few growth or non-growth cells for stratified sampling.", call. = FALSE)
  }
  p1 <- terra::spatSample(g1, s1, method = "random", na.rm = TRUE, as.points = TRUE)
  p0 <- terra::spatSample(g0, s0, method = "random", na.rm = TRUE, as.points = TRUE)
  list(vect = rbind(p1, p0), n_pos = s1, n_neg = s0)
}

extract_regression_frame <- function(growth_r, pred_stk, sample_vect) {
  eg <- terra::extract(growth_r, sample_vect, ID = FALSE, na.rm = FALSE)
  g <- eg[[1]]
  x <- terra::extract(pred_stk, sample_vect, ID = FALSE, na.rm = FALSE)
  xy <- terra::crds(sample_vect)
  df <- cbind.data.frame(growth = g, x = xy[, 1], y = xy[, 2], x)
  # Keep only fully observed rows for modelling, but report how many are dropped.
  ok <- stats::complete.cases(df)
  n_drop <- sum(!ok)
  if (n_drop > 0L) {
    pipeline_log("Phase 6: dropped ", n_drop, " / ", nrow(df), " sampled points due to NA predictors.")
  }
  df <- df[ok, , drop = FALSE]
  df
}

# --------------------------------------------------------------------------- #
# Models
# --------------------------------------------------------------------------- #

glm_logistic_fit <- function(df) {
  rhs <- setdiff(names(df), c("growth", "x", "y"))
  f <- stats::as.formula(paste("growth ~", paste(rhs, collapse = " + ")))
  stats::glm(f, data = df, family = stats::binomial())
}

glm_results_table <- function(m) {
  s <- summary(m)$coefficients
  or <- exp(s[, 1])
  data.frame(
    term = rownames(s),
    estimate = s[, 1],
    std_error = s[, 2],
    z_value = s[, 3],
    p_value = s[, 4],
    odds_ratio = or,
    pseudo_r2 = 1 - m$deviance / m$null.deviance,
    stringsAsFactors = FALSE
  )
}

lagsarlm_fit <- function(df, listw) {
  rhs <- setdiff(names(df), c("growth", "x", "y"))
  f <- stats::as.formula(paste("growth ~", paste(rhs, collapse = " + ")))
  spatialreg::lagsarlm(f, data = df, listw = listw, zero.policy = TRUE, NA.action = stats::na.exclude)
}

moran_glm_residuals <- function(m, df, listw) {
  r <- stats::residuals(m, type = "response")
  spdep::moran.test(r, listw, zero.policy = TRUE)
}

knn_listw <- function(df, k = 8L) {
  xy <- as.matrix(df[, c("x", "y")])
  # Identical sample coords (common on coarse grids) break kd-tree kNN; tiny jitter fixes Moran/lagsarlm.
  dup <- duplicated(xy) | duplicated(xy, fromLast = TRUE)
  if (any(dup)) {
    ext <- apply(xy, 2L, function(z) diff(range(z)))
    ext <- ext[is.finite(ext) & ext > 0]
    eps <- if (length(ext)) max(sqrt(.Machine$double.eps), 1e-6 * min(ext)) else sqrt(.Machine$double.eps)
    u <- stats::runif(nrow(xy), 0, eps)
    v <- stats::runif(nrow(xy), 0, eps)
    xy[dup, 1L] <- xy[dup, 1L] + u[dup]
    xy[dup, 2L] <- xy[dup, 2L] + v[dup]
  }
  k <- as.integer(k)
  k_eff <- min(k, max(1L, nrow(xy) - 1L))
  kn <- spdep::knearneigh(xy, k = k_eff)
  nb <- spdep::knn2nb(kn)
  spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
}

run_gwr_safe <- function(df_sp, formula_str) {
  f <- stats::as.formula(formula_str)
  bw <- tryCatch(
    GWmodel::bw.gwr(
      f,
      data = df_sp,
      approach = "AIC",
      kernel = "bisquare",
      adaptive = TRUE,
      longlat = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(bw)) return(NULL)
  tryCatch(
    GWmodel::gwr.basic(
      f,
      data = df_sp,
      bw = bw,
      kernel = "bisquare",
      adaptive = TRUE,
      longlat = FALSE
    ),
    error = function(e) NULL
  )
}

# --------------------------------------------------------------------------- #
# Main — stepped runs (low RAM): build → sample → models
# --------------------------------------------------------------------------- #

phase6_verbose_tools <- function(verbose) {
  log_step <- function(...) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    pipeline_log(paste0(..., collapse = ""))
  }
  step_time <- function(label, expr) {
    t0 <- Sys.time()
    log_step(label, " …")
    on.exit(
      {
        dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        log_step(label, " done (", sprintf("%.1f", dt), "s).")
      },
      add = TRUE
    )
    force(expr)
  }
  list(log_step = log_step, step_time = step_time)
}

# terra::terraOptions() prints a long table to the console; capture.output keeps logs readable.
phase6_push_terra_progress <- function(verbose) {
  old <- NULL
  utils::capture.output({
    old <- tryCatch(terra::terraOptions()$progress, error = function(e) NULL)
  })
  utils::capture.output(
    terra::terraOptions(progress = if (isTRUE(verbose)) 1L else 0L),
    type = "output"
  )
  old
}

phase6_pop_terra_progress <- function(old) {
  if (is.null(old)) return(invisible(NULL))
  utils::capture.output(
    try(terra::terraOptions(progress = old), silent = TRUE),
    type = "output"
  )
  invisible(NULL)
}

phase6_require_models <- function(include_gwr = FALSE) {
  req <- c("spdep", "spatialreg", "sp")
  if (isTRUE(include_gwr)) req <- c(req, "GWmodel")
  miss <- req[!vapply(req, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss)) {
    inst <- if (isTRUE(include_gwr)) {
      "install.packages(c(\"spdep\",\"spatialreg\",\"sp\",\"GWmodel\"))"
    } else {
      "install.packages(c(\"spdep\",\"spatialreg\",\"sp\"))"
    }
    stop("Install Phase 6 modelling packages: ", paste(miss, collapse = ", "), "\n", inst, call. = FALSE)
  }
  invisible(NULL)
}

#' Build predictor stack and write GeoTIFFs (then frees RAM). Run this first, or use `run_phase6()`.
run_phase6_build <- function(cfg = load_config(), verbose = TRUE) {
  with_pipeline_phase("Phase 6 — BUILD stack (LULC, roads, predictors → GeoTIFF)", {
  d <- phase6_dirs(cfg)
  dir.create(d$proc, showWarnings = FALSE, recursive = TRUE)
  dir.create(d$outd, showWarnings = FALSE, recursive = TRUE)
  old_prog <- phase6_push_terra_progress(verbose)
  on.exit(phase6_pop_terra_progress(old_prog), add = TRUE)
  vt <- phase6_verbose_tools(verbose)
  log_step <- vt$log_step
  step_time <- vt$step_time
  p6 <- get_phase6_cfg(cfg)
  plog <- if (isTRUE(verbose)) pipeline_log else function(...) NULL
  log_step(
    "Phase 6 build: growth ",
    p6$growth_year_t0,
    "→",
    p6$growth_year_t1,
    " | writing predictors + growth to disk (then gc)."
  )

  built <- step_time(
    "Build predictor stack (LULC/features; may fetch OSM)",
    build_phase6_stack(cfg, progress_log = plog)
  )
  pred_path <- phase6_predictors_path(cfg)
  growth_path <- phase6_growth_path(cfg)
  log_step("Predictor layers: ", paste(names(built$predictors), collapse = ", "))

  step_time(
    paste0("Write predictors -> ", pred_path),
    terra::writeRaster(
      built$predictors,
      pred_path,
      overwrite = TRUE,
      gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
    )
  )
  step_time(
    paste0("Write growth raster -> ", growth_path),
    terra::writeRaster(
      built$growth,
      growth_path,
      overwrite = TRUE,
      gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
    )
  )
  rm(built)
  gc(verbose = FALSE)
  log_step("Phase 6 build done. Next: run_phase6_sample() (or run_phase6(step = \"sample\")).")
  invisible(list(predictors_path = pred_path, growth_path = growth_path))
  })
}

#' Stratified sample from on-disk rasters; writes `regression_samples.csv`. Middle step.
run_phase6_sample <- function(cfg = load_config(), seed = 42L, verbose = TRUE) {
  with_pipeline_phase("Phase 6 — SAMPLE (points + extract → regression_samples.csv)", {
  d <- phase6_dirs(cfg)
  dir.create(d$proc, showWarnings = FALSE, recursive = TRUE)
  old_prog <- phase6_push_terra_progress(verbose)
  on.exit(phase6_pop_terra_progress(old_prog), add = TRUE)
  p6 <- get_phase6_cfg(cfg)
  n_samp <- as.integer(p6$n_sample_per_stratum %||% 1000L)
  vt <- phase6_verbose_tools(verbose)
  log_step <- vt$log_step
  step_time <- vt$step_time

  pred_path <- phase6_predictors_path(cfg)
  growth_path <- phase6_growth_path(cfg)
  if (!file.exists(pred_path) || !file.exists(growth_path)) {
    stop(
      "Missing GeoTIFFs. Run run_phase6_build() first.\n  ",
      pred_path,
      "\n  ",
      growth_path,
      call. = FALSE
    )
  }
  log_step("Phase 6 sample: n_sample_per_stratum=", n_samp, " | reading rasters from disk.")

  growth <- terra::rast(growth_path)
  pred <- terra::rast(pred_path)
  sp_list <- step_time("Stratified sampling points", stratified_sample_points(growth, n_samp, seed))
  df <- step_time(
    "Extract regression predictors at sample points",
    extract_regression_frame(growth, pred, sp_list$vect)
  )
  rm(growth, pred, sp_list)
  gc(verbose = FALSE)

  samp_path <- phase6_samples_path(cfg)
  step_time(paste0("Write sample table -> ", samp_path), utils::write.csv(df, samp_path, row.names = FALSE))
  log_step(
    "Samples: n=",
    nrow(df),
    " | positives=",
    sum(df$growth == 1, na.rm = TRUE),
    " | negatives=",
    sum(df$growth == 0, na.rm = TRUE)
  )
  invisible(df)
  })
}

#' GLM, Moran, lagsarlm, correlations, plot; optional GWR (heavy — off by default in yml).
run_phase6_models <- function(cfg = load_config(), verbose = TRUE, run_gwr = NULL) {
  with_pipeline_phase("Phase 6 — MODELS (GLM, Moran, lagsarlm, plots; optional GWR)", {
  d <- phase6_dirs(cfg)
  dir.create(d$outd, showWarnings = FALSE, recursive = TRUE)
  old_prog <- phase6_push_terra_progress(verbose)
  on.exit(phase6_pop_terra_progress(old_prog), add = TRUE)
  p6 <- get_phase6_cfg(cfg)
  if (is.null(run_gwr)) run_gwr <- isTRUE(p6$run_gwr %||% FALSE)
  phase6_require_models(include_gwr = isTRUE(run_gwr))

  k_knn <- as.integer(p6$knn_k %||% 8L)
  vt <- phase6_verbose_tools(verbose)
  log_step <- vt$log_step
  step_time <- vt$step_time

  samp_path <- phase6_samples_path(cfg)
  if (!file.exists(samp_path)) {
    stop("Missing ", samp_path, ". Run run_phase6_sample() first.", call. = FALSE)
  }
  df <- utils::read.csv(samp_path, check.names = FALSE)
  if (nrow(df) < 10L) {
    stop("Too few regression rows (", nrow(df), "). Check sampling and NA predictors.", call. = FALSE)
  }

  y1 <- as.integer(p6$growth_year_t1 %||% 2024)
  template <- terra::rast(file.path(project_root(), path_lulc(y1, cfg)))

  log_step(
    "Phase 6 models: knn_k=",
    k_knn,
    " | run_gwr=",
    isTRUE(run_gwr),
    " (GWR is memory/CPU heavy)."
  )

  rhs <- setdiff(names(df), c("growth", "x", "y"))
  formula_str <- paste("growth ~", paste(rhs, collapse = " + "))

  m_glm <- step_time("Fit logistic regression (glm)", glm_logistic_fit(df))
  glm_tab <- glm_results_table(m_glm)
  glm_path <- file.path(d$outd, "logistic_regression_results.csv")
  step_time(paste0("Write GLM results -> ", glm_path), utils::write.csv(glm_tab, glm_path, row.names = FALSE))

  lw <- step_time("Build kNN spatial weights (listw)", knn_listw(df, k_knn))
  w_path <- file.path(d$proc, "spatial_weights.rds")
  step_time(paste0("Write spatial weights -> ", w_path), saveRDS(list(k = k_knn, listw = lw), w_path))

  m_lag <- step_time("Fit spatial lag model (lagsarlm)", tryCatch(lagsarlm_fit(df, lw), error = function(e) NULL))
  if (!is.null(m_lag)) {
    lag_path <- file.path(d$outd, "lagsarlm_summary.txt")
    step_time(
      paste0("Write lagsarlm summary -> ", lag_path),
      {
        sink(lag_path)
        print(summary(m_lag))
        sink()
      }
    )
  }

  mi <- step_time("Moran's I on GLM residuals", tryCatch(moran_glm_residuals(m_glm, df, lw), error = function(e) NULL))
  if (!is.null(mi)) {
    mor_path <- file.path(d$outd, "moran_residuals_glm.csv")
    utils::write.csv(
      data.frame(
        moran_statistic = unname(mi$estimate[1]),
        p_value = mi$p.value,
        stringsAsFactors = FALSE
      ),
      mor_path,
      row.names = FALSE
    )
    log_step("Wrote ", mor_path)
  }

  cor_cols <- c("growth", rhs)
  cor_df <- df[, intersect(cor_cols, names(df)), drop = FALSE]
  C <- stats::cor(cor_df, use = "pairwise.complete.obs", method = "pearson")
  cor_path <- file.path(d$outd, "correlation_matrix.csv")
  utils::write.csv(
    cbind(variable = rownames(C), as.data.frame(C)),
    cor_path,
    row.names = FALSE
  )
  log_step("Wrote ", cor_path)

  gwr <- NULL
  if (isTRUE(run_gwr)) {
    sf_pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = terra::crs(template))
    sp_pts <- as(sf_pts, "Spatial")
    gwr <- step_time("GWR (bandwidth search + fit) [heavy]", run_gwr_safe(sp_pts, formula_str))
    if (!is.null(gwr) && !is.null(gwr$SDF)) {
      sf_gwr <- sf::st_as_sf(gwr$SDF)
      gwr_path <- file.path(d$outd, "gwr_results.shp")
      sf::st_write(sf_gwr, gwr_path, delete_dsn = TRUE, quiet = TRUE)
      log_step("Wrote ", gwr_path)
      cn <- names(gwr$SDF@data)
      lr2col <- cn[grepl("localR2|^LocalR2$", cn, ignore.case = TRUE)][1]
      if (!is.na(lr2col)) {
        coords <- sp::coordinates(gwr$SDF)
        lr2 <- gwr$SDF@data[[lr2col]]
        gv <- terra::vect(
          data.frame(x = coords[, 1], y = coords[, 2], lr2 = lr2),
          geom = c("x", "y"),
          crs = terra::crs(template)
        )
        r2map <- terra::rasterize(gv, template, field = "lr2", fun = "mean")
        names(r2map) <- "local_R2"
        r2_path <- file.path(d$outd, "gwr_local_r2.tif")
        terra::writeRaster(
          r2map,
          r2_path,
          overwrite = TRUE,
          gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
        )
        log_step("Wrote ", r2_path)
      }
    }
  } else {
    log_step("Skipping GWR (set phase6.run_gwr: true or run_phase6_models(run_gwr = TRUE)).")
  }

  cf <- stats::coef(m_glm)
  cf <- cf[names(cf) != "(Intercept)"]
  imp <- data.frame(variable = names(cf), abs_coef = abs(cf), stringsAsFactors = FALSE)
  imp <- imp[order(-imp$abs_coef), ]

  p <- ggplot2::ggplot(imp, ggplot2::aes(x = reorder(variable, abs_coef), y = abs_coef)) +
    ggplot2::geom_col(fill = "darkred") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Driver importance (|logistic coefficient|)",
      x = NULL,
      y = "|Coefficient|"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(
    file.path(d$outd, "driver_importance_chart.png"),
    p,
    width = 7,
    height = 5,
    dpi = 150
  )

  log_step("Phase 6 models complete. See ", d$proc, " and ", d$outd)
  invisible(list(glm = m_glm, samples = df, gwr = gwr))
  })
}

#' @param step `"all"` runs build → sample → models with `gc()` between. Or run one of `"build"`, `"sample"`, `"models"`.
#' @param run_gwr Passed to models step when `step` is `"all"` or `"models"`; `NULL` uses `phase6.run_gwr` in config (default FALSE).
run_phase6 <- function(
    cfg = load_config(),
    seed = 42L,
    verbose = TRUE,
    step = c("all", "build", "sample", "models"),
    run_gwr = NULL) {
  step <- match.arg(step)
  if (identical(step, "all")) {
    run_phase6_build(cfg, verbose)
    gc(verbose = FALSE)
    run_phase6_sample(cfg, seed, verbose)
    gc(verbose = FALSE)
    return(invisible(run_phase6_models(cfg, verbose, run_gwr = run_gwr)))
  }
  if (identical(step, "build")) return(invisible(run_phase6_build(cfg, verbose)))
  if (identical(step, "sample")) return(invisible(run_phase6_sample(cfg, seed, verbose)))
  invisible(run_phase6_models(cfg, verbose, run_gwr = run_gwr))
}

if (interactive()) {
  message(
    "Phase 6: run_phase6_build(); run_phase6_sample(); run_phase6_models() — or run_phase6() for all. ",
    "GWR off by default (phase6.run_gwr in yml)."
  )
}
