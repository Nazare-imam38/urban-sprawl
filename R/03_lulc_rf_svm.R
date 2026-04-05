#' Phase 3 — LULC: Random Forest + SVM (4 classes), dual sensor, validation, maps, stats
#'
#' Option A: separate models for Sentinel-2 stacks vs Landsat stacks (feature names differ).
#' Training: `sf` with integer `class` (1–4). **Default:** `classification.auto_training` builds
#' points from **ESA WorldCover** (Planetary Computer STAC) or optional **GlobeLand30** GeoTIFF —
#' no manual digitizing. Manual: `collect_training_mapedit()`.
#' Outputs: `lulc_{year}.tif`, saved models, `accuracy_report.csv`, `variable_importance.csv`,
#'          `urban_area_stats.csv`.
#'
#' Post-classification: 3×3 focal majority; spectral rules on **denormalized** NDVI/MNDWI
#' (from `normalization_params.csv` + normalized feature stack).
#' Large AOIs: set `classification.prediction` in `study_area.yml` (memory_safe, tiles) or
#' `predict_phase3_year()` / `run_phase3(predict=FALSE)` to split training vs prediction.
#'
#' Working directory = project root.

source("R/00_setup.R")

suppressPackageStartupMessages(library(caret))

CLASS_LABELS <- c("Urban", "Vegetation", "BareSoil", "Water")
CLASS_CODES <- 1L:4L
names(CLASS_LABELS) <- as.character(CLASS_CODES)

`%??%` <- function(x, y) {
  if (length(x) == 1L && !is.na(x)) x else y
}

# --------------------------------------------------------------------------- #
# Paths & sensor detection
# --------------------------------------------------------------------------- #

path_features <- function(year, cfg) {
  file.path(cfg$paths$processed_dir, sprintf("features_%d.tif", year))
}

path_lulc <- function(year, cfg) {
  file.path(cfg$paths$processed_dir, sprintf("lulc_%d.tif", year))
}

path_training_default <- function(cfg) {
  file.path(cfg$paths$processed_dir, "training_samples.gpkg")
}

sensor_from_stack <- function(r) {
  nm <- names(r)
  if (any(grepl("^SR_B", nm))) return("landsat")
  if (any(nm == "B02")) return("sentinel")
  stop("Cannot detect sensor from layer names: ", paste(nm, collapse = ", "), call. = FALSE)
}

feature_columns <- function(r) {
  names(r)
}

landsat_training_year <- function(cfg) {
  acq <- cfg$planetary_computer$landsat$acquisitions
  ys <- vapply(acq, function(z) as.integer(z$year), integer(1))
  max(ys)
}

sentinel_training_year <- function(cfg) {
  max(unlist(cfg$planetary_computer$sentinel2$years, use.names = FALSE))
}

phase3_years <- function(cfg) {
  pc <- cfg$planetary_computer
  c(
    vapply(pc$landsat$acquisitions, function(z) as.integer(z$year), integer(1)),
    unlist(pc$sentinel2$years, use.names = FALSE)
  )
}

# --------------------------------------------------------------------------- #
# Interactive training collection (mapedit)
# --------------------------------------------------------------------------- #

#' Draw points on an RGB from a Phase 2 stack; add `class` in QGIS or R after draw, or
#' use separate draws per class and `sf::st_as_sf(rbind(...))` with class column.
#' Saves GeoPackage (recommended). For Shapefile use path ending in .shp.
#'
#' @param class_value Integer 1–4 for all features drawn in this session (draw multiple times per class).
collect_training_mapedit <- function(
    cfg = load_config(),
    ref_year = NULL,
    class_value = 1L,
    out_path = NULL) {
  if (!requireNamespace("mapedit", quietly = TRUE)) stop("Install mapedit", call. = FALSE)
  if (!requireNamespace("mapview", quietly = TRUE)) {
    stop("Install mapview for interactive training: install.packages(\"mapview\")", call. = FALSE)
  }
  if (is.null(ref_year)) ref_year <- sentinel_training_year(cfg)
  fp <- path_features(ref_year, cfg)
  if (!file.exists(fp)) stop("Missing features raster: ", fp, call. = FALSE)
  r <- terra::rast(fp)
  rgb <- r[[c("B04", "B03", "B02")]]
  if (any(!c("B04", "B03", "B02") %in% names(r))) {
    rgb <- r[[c("SR_B4", "SR_B3", "SR_B2")]]
    names(rgb) <- c("R", "G", "B")
  } else {
    names(rgb) <- c("R", "G", "B")
  }
  message("Draw training geometry; assign class ", class_value, " (repeat per class with different class_value).")
  drawn <- mapedit::drawFeatures(mapview::mapview(rgb))
  if (inherits(drawn, "try-error") || nrow(drawn) == 0) stop("No features drawn.", call. = FALSE)
  drawn$class <- as.integer(class_value)
  if (is.null(out_path)) out_path <- path_training_default(cfg)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(drawn, out_path, delete_dsn = file.exists(out_path), quiet = TRUE)
  message("Appended/saved: ", out_path, " — merge multiple class draws with rbind or in GIS.")
  invisible(drawn)
}

# --------------------------------------------------------------------------- #
# Automated training samples (WorldCover STAC or GlobeLand30 raster)
# --------------------------------------------------------------------------- #

`%||%` <- function(x, y) if (is.null(x)) y else x

get_auto_training_cfg <- function(cfg) {
  cl <- cfg$classification %||% list()
  at <- cl$auto_training %||% list()
  dt <- at$datetime
  if (is.null(dt) || (length(dt) == 1L && !nzchar(trimws(as.character(dt))))) {
    dt <- NULL
  } else {
    dt <- as.character(dt)[1]
  }
  list(
    enable = isTRUE(at$enable %||% TRUE),
    source = tolower(as.character(at$source %||% "esa_worldcover")),
    # Planetary Computer collection id is "esa-worldcover" (not esa-worldcover-v100).
    stac_collection = as.character(at$stac_collection %||% "esa-worldcover"),
    datetime = dt %||% "2020-01-01/2020-12-31",
    asset_name = as.character(at$asset_name %||% "map"),
    n_per_class = as.integer(at$n_per_class %||% 400L),
    globeland30_path = at$globeland30_path,
    ref_align_year = at$ref_align_year
  )
}

sign_pc_href <- function(href, sign_endpoint) {
  if (!nzchar(href)) stop("Empty asset href", call. = FALSE)
  q <- paste0(sign_endpoint, "?href=", utils::URLencode(href, reserved = TRUE))
  r <- httr::GET(q)
  if (httr::status_code(r) >= 400) {
    stop("PC sign failed (", httr::status_code(r), ")", call. = FALSE)
  }
  txt <- httr::content(r, "text", encoding = "UTF-8")
  j <- jsonlite::fromJSON(txt)
  out <- j$href
  if (is.null(out) || !nzchar(out)) stop("Sign response missing href", call. = FALSE)
  out
}

phase3_stac_search_feats <- function(cfg, collection, datetime, limit) {
  sa <- cfg$study_area
  bbox <- c(sa$xmin, sa$ymin, sa$xmax, sa$ymax)
  pc <- cfg$planetary_computer
  s <- rstac::stac(pc$stac_url)
  do_req <- function(dt_str) {
    if (is.null(dt_str) || !nzchar(dt_str)) {
      rstac::stac_search(
        s,
        collections = collection,
        bbox = bbox,
        limit = as.integer(limit)
      )
    } else {
      rstac::stac_search(
        s,
        collections = collection,
        bbox = bbox,
        datetime = dt_str,
        limit = as.integer(limit)
      )
    }
  }
  doc <- tryCatch(rstac::get_request(do_req(datetime)), error = function(e) NULL)
  feats <- if (!is.null(doc)) doc$features else NULL
  if (!is.null(feats) && length(feats)) return(feats)
  if (!is.null(datetime) && nzchar(datetime)) {
    message("STAC returned no items with datetime ", datetime, "; retrying without datetime …")
    doc2 <- tryCatch(rstac::get_request(do_req(NULL)), error = function(e) NULL)
    if (!is.null(doc2) && !is.null(doc2$features) && length(doc2$features)) return(doc2$features)
  }
  NULL
}

phase3_pick_raster_asset <- function(f, asset_names) {
  ast <- f$assets
  if (is.null(ast)) return(NULL)
  nms_all <- names(ast)
  if (is.null(nms_all)) return(NULL)
  want <- unique(c(asset_names, "map", "classification", "data", "landcover"))
  want <- want[!vapply(want, function(w) is.null(w) || !nzchar(w), logical(1))]
  for (nm in want) {
    if (!nm %in% nms_all) next
    a <- ast[[nm]]
    if (is.null(a)) next
    h <- a$href
    if (!is.null(h) && nzchar(as.character(h)[1])) {
      return(list(name = nm, href = as.character(h)[1]))
    }
  }
  for (nm in nms_all) {
    a <- ast[[nm]]
    if (is.null(a)) next
    h <- a$href
    if (is.null(h)) next
    hs <- as.character(h)[1]
    if (grepl("\\.(tif|tiff)(\\?|$)", hs, ignore.case = TRUE)) {
      return(list(name = nm, href = hs))
    }
  }
  NULL
}

phase3_stac_feature_with_assets <- function(cfg, collection, datetime, asset_names, limit = 200L) {
  feats <- phase3_stac_search_feats(cfg, collection, datetime, limit)
  if (is.null(feats) || !length(feats)) return(NULL)
  an <- if (length(asset_names)) asset_names else "map"
  for (f in feats) {
    hit <- phase3_pick_raster_asset(f, an)
    if (!is.null(hit)) {
      return(list(feature = f, href = hit$href, asset = hit$name))
    }
  }
  NULL
}

rast_from_pc_href <- function(href, sign_endpoint) {
  signed <- sign_pc_href(href, sign_endpoint)
  terra::rast(paste0("/vsicurl/", signed))
}

remap_esa_worldcover_to_lulc4 <- function(r) {
  f <- function(z, ...) {
    u <- suppressWarnings(as.integer(round(z)))
    y <- rep(NA_integer_, length(u))
    y[u == 50L] <- 1L
    y[u %in% c(10L, 20L, 30L, 40L, 90L, 95L, 96L)] <- 2L
    y[u %in% c(60L, 70L)] <- 3L
    y[u == 80L] <- 4L
    y
  }
  terra::app(r, fun = f)
}

remap_globeland30_to_lulc4 <- function(r) {
  f <- function(z, ...) {
    u <- suppressWarnings(as.integer(round(z)))
    y <- rep(NA_integer_, length(u))
    y[u == 80L] <- 1L
    y[u %in% c(10L, 20L, 30L, 40L, 50L, 70L, 100L)] <- 2L
    y[u == 90L] <- 3L
    y[u == 60L] <- 4L
    y
  }
  terra::app(r, fun = f)
}

#' Build `training_samples.gpkg` from a reference land-cover grid aligned to Phase 2 features.
#'
#' @param cfg Config from `load_config()`.
#' @param out_path Output GeoPackage path (default `data/processed/training_samples.gpkg`).
#' @param ref_year Year whose `features_{year}.tif` defines grid and CRS (default: latest Sentinel year).
#' @details Sources: **esa_worldcover** — STAC collection + datetime in `classification.auto_training`;
#'   **globeland30** — set `auto_training.globeland30_path` to a merged GL30 GeoTIFF.
auto_build_training_samples <- function(
    cfg = load_config(),
    out_path = NULL,
    ref_year = NULL) {
  if (!requireNamespace("rstac", quietly = TRUE)) stop("Install rstac", call. = FALSE)
  at <- get_auto_training_cfg(cfg)
  if (!isTRUE(at$enable)) {
    stop("classification.auto_training.enable is FALSE; provide training_samples.gpkg or enable auto_training.", call. = FALSE)
  }
  ys <- sentinel_training_year(cfg)
  if (is.null(ref_year)) ref_year <- at$ref_align_year %||% ys
  ref_year <- as.integer(ref_year)
  fp <- path_features(ref_year, cfg)
  if (!file.exists(fp)) {
    stop("Need Phase 2 features for ref year ", ref_year, ": ", fp, call. = FALSE)
  }
  feat <- terra::rast(fp)
  tmpl <- feat[[1]]
  sign_ep <- cfg$planetary_computer$sign_endpoint

  src <- at$source
  if (src == "globeland30") {
    gp <- at$globeland30_path
    if (is.null(gp) || !nzchar(gp) || !file.exists(gp)) {
      stop("auto_training.source=globeland30 requires classification.auto_training.globeland30_path to an existing .tif", call. = FALSE)
    }
    message("Auto-training: reading GlobeLand30 raster ", gp)
    raw <- terra::rast(gp)
    r4 <- remap_globeland30_to_lulc4(raw)
  } else {
    if (!src %in% c("esa_worldcover", "worldcover")) {
      stop("Unknown auto_training.source: ", src, " (use esa_worldcover or globeland30).", call. = FALSE)
    }
    anames <- trimws(strsplit(at$asset_name, ",", fixed = TRUE)[[1]])
    message(
      "Auto-training: STAC ", at$stac_collection, " (", at$datetime, ") assets: ",
      paste(anames, collapse = ", "), " …"
    )
    hit <- phase3_stac_feature_with_assets(
      cfg,
      at$stac_collection,
      at$datetime,
      anames,
      limit = 250L
    )
    if (is.null(hit)) {
      stop(
        "No STAC items or raster assets for collection ", at$stac_collection,
        " over the study bbox. On Planetary Computer use stac_collection: esa-worldcover (not esa-worldcover-v100). ",
        "Or set classification.auto_training.source: globeland30 with globeland30_path.",
        call. = FALSE
      )
    }
    message("Using STAC item ", hit$feature$id %||% "<id>", " asset: ", hit$asset)
    raw <- tryCatch(
      rast_from_pc_href(hit$href, sign_ep),
      error = function(e) stop("Could not open WorldCover COG: ", conditionMessage(e), call. = FALSE)
    )
    if (terra::nlyr(raw) > 1L) raw <- raw[[1L]]
    r4 <- remap_esa_worldcover_to_lulc4(raw)
  }

  if (terra::nlyr(r4) > 1L) r4 <- r4[[1L]]
  sa <- cfg$study_area
  bb <- sf::st_bbox(
    c(xmin = sa$xmin, ymin = sa$ymin, xmax = sa$xmax, ymax = sa$ymax),
    crs = sa$crs_wgs84
  )
  aoi_w <- sf::st_as_sf(sf::st_as_sfc(bb))
  r4 <- terra::crop(r4, terra::vect(aoi_w), snap = "out")

  r4 <- terra::project(r4, tmpl, method = "near", gdal = TRUE)
  r4 <- terra::resample(r4, tmpl, method = "near")
  names(r4) <- "lulc4"

  n <- max(50L, at$n_per_class)
  tab <- terra::freq(r4, digits = 6)
  if (is.null(tab) || !nrow(tab)) stop("Reference map has no classified pixels in AOI.", call. = FALSE)
  if ("layer" %in% names(tab)) tab <- tab[tab$layer == 1L, , drop = FALSE]
  tab <- tab[tab$value %in% 1:4 & tab$count > 0L, , drop = FALSE]
  if (!nrow(tab)) stop("No pixels mapped to classes 1–4 after reclass; check AOI and source product.", call. = FALSE)

  message(
    "Class counts on grid: ",
    paste(sprintf("%s=%s", tab$value, tab$count), collapse = ", ")
  )

  pts <- tryCatch(
    terra::spatSample(r4, n, method = "stratified", na.rm = TRUE, as.points = TRUE),
    error = function(e) {
      message("Stratified spatSample failed (", conditionMessage(e), "); trying class-wise random samples …")
      NULL
    }
  )
  if (is.null(pts)) {
    lst <- vector("list", nrow(tab))
    for (i in seq_len(nrow(tab))) {
      cls <- as.integer(tab$value[i])
      ni <- min(n, as.integer(tab$count[i]))
      if (ni < 1L) next
      m <- terra::ifel(r4 == cls, 1L, NA_integer_)
      lst[[i]] <- terra::spatSample(m, ni, method = "random", na.rm = TRUE, as.points = TRUE)
    }
    lst <- lst[!vapply(lst, is.null, logical(1))]
    if (!length(lst)) stop("Could not draw any training points.", call. = FALSE)
    pts <- do.call(rbind, lst)
  }

  sf_pts <- sf::st_as_sf(pts)
  ex <- terra::extract(r4, terra::vect(sf_pts), ID = FALSE, na.rm = TRUE)
  cls_col <- names(ex)[1]
  sf_pts$class <- as.integer(ex[[cls_col]])
  sf_pts <- sf_pts[!is.na(sf_pts$class) & sf_pts$class >= 1L & sf_pts$class <= 4L, ]
  if (nrow(sf_pts) < 40L) {
    stop("Too few auto training points (", nrow(sf_pts), "); widen AOI or increase n_per_class.", call. = FALSE)
  }

  if (is.null(out_path)) out_path <- path_training_default(cfg)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  gcol <- attr(sf_pts, "sf_column")
  sf_out <- sf_pts[, c("class", gcol)]
  sf::st_write(sf_out, out_path, delete_dsn = TRUE, quiet = TRUE)
  message("Wrote ", nrow(sf_pts), " auto-generated training points: ", out_path)
  invisible(sf_pts)
}

# --------------------------------------------------------------------------- #
# Extract training matrix from features + sf samples
# --------------------------------------------------------------------------- #

prepare_samples_sf <- function(samples_sf, class_col = "class") {
  if (!class_col %in% names(samples_sf)) {
    stop("Samples must have column `", class_col, "` with values 1–4.", call. = FALSE)
  }
  g <- sf::st_geometry(samples_sf)
  if (inherits(g, "sfc_POINT")) {
    return(samples_sf)
  }
  if (inherits(g, "sfc_POLYGON") || inherits(g, "sfc_MULTIPOLYGON")) {
    p <- sf::st_point_on_surface(samples_sf)
    return(p)
  }
  stop("Geometry must be POINT, POLYGON, or MULTIPOLYGON.", call. = FALSE)
}

extract_training_table <- function(samples_sf, features_r, class_col = "class") {
  v <- terra::vect(prepare_samples_sf(samples_sf, class_col))
  if (!terra::same.crs(v, features_r)) {
    v <- terra::project(v, terra::crs(features_r))
  }
  ex <- terra::extract(features_r, v, fun = mean, na.rm = TRUE, ID = FALSE)
  y <- as.integer(sf::st_drop_geometry(samples_sf)[[class_col]])
  if (length(y) != nrow(ex)) stop("Sample / extract row mismatch.", call. = FALSE)
  bad <- y < 1L | y > 4L | is.na(y)
  if (any(bad)) stop("All class values must be integers 1–4.", call. = FALSE)
  df <- as.data.frame(ex)
  df$class <- factor(y, levels = CLASS_CODES, labels = CLASS_LABELS)
  df
}

stratified_train_test <- function(df, class_col = "class", p_train = 0.7, seed = 42) {
  set.seed(seed)
  idx <- caret::createDataPartition(df[[class_col]], p = p_train, list = FALSE)[, 1]
  list(train = df[idx, , drop = FALSE], test = df[-idx, , drop = FALSE])
}

#' Same row indices for Sentinel + Landsat tables (paired geographic samples).
stratified_train_test_paired <- function(df_a, df_b, class_col = "class", p_train = 0.7, seed = 42) {
  if (nrow(df_a) != nrow(df_b)) {
    stop(
      "Paired split requires same n: Sentinel ", nrow(df_a), " vs Landsat ", nrow(df_b),
      call. = FALSE
    )
  }
  set.seed(seed)
  idx <- caret::createDataPartition(df_a[[class_col]], p = p_train, list = FALSE)[, 1]
  list(
    train_a = df_a[idx, , drop = FALSE],
    test_a = df_a[-idx, , drop = FALSE],
    train_b = df_b[idx, , drop = FALSE],
    test_b = df_b[-idx, , drop = FALSE]
  )
}

check_min_samples <- function(df, class_col = "class", min_per_class = 100L) {
  tab <- table(df[[class_col]])
  low <- names(tab)[tab < min_per_class]
  if (length(low)) {
    warning("Classes with < ", min_per_class, " samples: ", paste(low, collapse = ", "),
            " — aim for 200–300 per class.", call. = FALSE)
  }
  invisible(tab)
}

# --------------------------------------------------------------------------- #
# Train RF + SVM
# --------------------------------------------------------------------------- #

train_rf_lulc <- function(df_train, class_col = "class") {
  feats <- setdiff(names(df_train), class_col)
  f <- stats::as.formula(paste(class_col, "~", paste(feats, collapse = " + ")))
  randomForest::randomForest(
    f,
    data = df_train,
    ntree = 500L,
    mtry = 4L,
    importance = TRUE,
    na.action = stats::na.omit
  )
}

train_svm_lulc <- function(df_train, class_col = "class", cost = 1, gamma = NULL, tune = FALSE) {
  feats <- setdiff(names(df_train), class_col)
  x <- as.matrix(df_train[, feats, drop = FALSE])
  y <- df_train[[class_col]]
  if (is.null(gamma)) gamma <- 1 / ncol(x)
  if (tune) {
    tun <- e1071::tune(
      e1071::svm,
      train.x = x,
      train.y = y,
      kernel = "radial",
      ranges = list(cost = c(0.1, 1, 10, 100), gamma = c(0.001, 0.01, 0.1, 1)),
      tunecontrol = e1071::tune.control(cross = 5)
    )
    return(tun$best.model)
  }
  e1071::svm(
    x = x,
    y = y,
    kernel = "radial",
    cost = cost,
    gamma = gamma,
    scale = TRUE,
    probability = TRUE
  )
}

predict_test_factor <- function(model, df_test, class_col = "class", type = c("rf", "svm")) {
  type <- match.arg(type)
  feats <- setdiff(names(df_test), class_col)
  lv <- levels(df_test[[class_col]])
  nd <- as.data.frame(df_test[, feats, drop = FALSE])
  if (type == "rf") {
    p <- stats::predict(model, newdata = nd)
  } else {
    p <- stats::predict(model, as.matrix(nd), type = "class")
  }
  if (is.matrix(p)) {
    p <- apply(p, 1L, function(z) names(which.max(z)))
  }
  if (is.factor(p)) {
    pc <- as.character(p)
  } else if (is.numeric(p)) {
    pc <- lv[as.integer(p)]
  } else {
    pc <- as.character(p)
  }
  if (length(pc) != nrow(nd)) {
    stop(
      "predict length ", length(pc), " != nrow(test) ", nrow(nd), " (", type, ")",
      call. = FALSE
    )
  }
  factor(pc, levels = lv)
}

#' Drop test rows with NA predictors; return aligned pred + ref for confusionMatrix.
predict_test_aligned <- function(model, df_test, class_col = "class", type = c("rf", "svm")) {
  type <- match.arg(type)
  feats <- setdiff(names(df_test), class_col)
  lv <- levels(df_test[[class_col]])
  nd <- as.data.frame(df_test[, feats, drop = FALSE])
  ok <- stats::complete.cases(nd)
  if (!any(ok)) {
    stop("No test rows with complete predictor values.", call. = FALSE)
  }
  if (any(!ok)) {
    message("Dropped ", sum(!ok), " test row(s) with NA predictors before confusion matrix.")
  }
  dsub <- df_test[ok, , drop = FALSE]
  nd <- nd[ok, , drop = FALSE]
  ref <- factor(as.character(dsub[[class_col]]), levels = lv)
  if (type == "rf") {
    p <- stats::predict(model, newdata = nd)
  } else {
    p <- stats::predict(model, as.matrix(nd), type = "class")
  }
  if (is.matrix(p)) {
    p <- apply(p, 1L, function(z) names(which.max(z)))
  }
  if (is.factor(p)) {
    pc <- as.character(p)
  } else if (is.numeric(p)) {
    pc <- lv[as.integer(p)]
  } else {
    pc <- as.character(p)
  }
  if (length(pc) != length(ref)) {
    stop("Internal: pred length ", length(pc), " != ref length ", length(ref), call. = FALSE)
  }
  list(pred = factor(pc, levels = lv), ref = ref)
}

# --------------------------------------------------------------------------- #
# Accuracy: OA, Kappa, PA, UA (caret confusionMatrix)
# --------------------------------------------------------------------------- #

confusion_metrics_df <- function(pred, ref, model_name, sensor_name) {
  if (length(pred) != length(ref)) {
    stop(
      "confusion_metrics_df: pred (n=", length(pred), ") and ref (n=", length(ref), ") must match.",
      call. = FALSE
    )
  }
  cm <- caret::confusionMatrix(pred, ref)
  oa <- as.numeric(cm$overall["Accuracy"])
  kap <- as.numeric(cm$overall["Kappa"])
  # caret:: $byClass shape varies (matrix vs vector) across versions — derive PA/UA from the table.
  # Rows = Predicted, Columns = Reference (caret docs).
  tbl <- as.matrix(cm$table)
  cn_ref <- colnames(tbl)
  cn_pred <- rownames(tbl)
  n <- length(cn_ref)
  pa_vec <- rep(NA_real_, n)
  ua_vec <- rep(NA_real_, n)
  names(pa_vec) <- cn_ref
  names(ua_vec) <- cn_ref
  for (k in seq_len(n)) {
    ref_lab <- cn_ref[k]
    ic <- match(ref_lab, cn_ref)[1]
    ir <- match(ref_lab, cn_pred)[1]
    if (is.na(ic) || is.na(ir)) next
    ctot <- sum(tbl[, ic])
    rtot <- sum(tbl[ir, ])
    if (ctot > 0) pa_vec[k] <- tbl[ir, ic] / ctot
    if (rtot > 0) ua_vec[k] <- tbl[ir, ic] / rtot
  }
  cls_names <- cn_ref
  out <- data.frame(
    model = model_name,
    sensor = sensor_name,
    metric = c("OA", "Kappa", paste0("PA_", cls_names), paste0("UA_", cls_names)),
    value = c(oa, kap, unname(pa_vec), unname(ua_vec)),
    stringsAsFactors = FALSE
  )
  list(table = as.data.frame(cm$table), long = out, cm = cm)
}

# --------------------------------------------------------------------------- #
# Variable importance (RF only)
# --------------------------------------------------------------------------- #

rf_importance_df <- function(rf_model, sensor_name) {
  imp <- randomForest::importance(rf_model)
  imp <- as.data.frame(imp)
  imp$feature <- rownames(imp)
  imp$sensor <- sensor_name
  imp[order(-imp$MeanDecreaseGini), ]
}

# --------------------------------------------------------------------------- #
# Raster prediction & post-classification
# --------------------------------------------------------------------------- #

align_stack_for_model <- function(r, feature_names) {
  miss <- setdiff(feature_names, names(r))
  if (length(miss)) stop("Raster missing layers: ", paste(miss, collapse = ", "), call. = FALSE)
  r[[feature_names]]
}

get_phase3_prediction_cfg <- function(cfg) {
  cl <- cfg$classification %??% list()
  pr <- cl$prediction %??% list()
  list(
    memory_safe = if (is.null(pr$memory_safe)) TRUE else isTRUE(pr$memory_safe),
    tiles_x = as.integer(pr$tiles_x %??% 2L),
    tiles_y = as.integer(pr$tiles_y %??% 2L),
    skip_focal = isTRUE(pr$skip_focal),
    memfrac = {
      v <- pr$memfrac
      if (is.null(v) || !is.numeric(v) || length(v) != 1L || !is.finite(v)) NA_real_ else as.numeric(v)
    },
    cleanup_temp = if (is.null(pr$cleanup_temp)) TRUE else isTRUE(pr$cleanup_temp)
  )
}

grid_tile_extents <- function(r, nx, ny) {
  e <- as.vector(terra::ext(r))
  xr <- seq(e[1], e[2], length.out = as.integer(nx) + 1L)
  yr <- seq(e[3], e[4], length.out = as.integer(ny) + 1L)
  lst <- vector("list", as.integer(nx) * as.integer(ny))
  k <- 0L
  for (i in seq_len(as.integer(nx))) {
    for (j in seq_len(as.integer(ny))) {
      k <- k + 1L
      lst[[k]] <- terra::ext(xr[i], xr[i + 1L], yr[j], yr[j + 1L])
    }
  }
  lst
}

gdal_wopt_compress <- function() {
  c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
}

predict_rf_raster_to_file <- function(r, rf_model, feature_names, filename, overwrite = TRUE) {
  rs <- align_stack_for_model(r, feature_names)
  terra::predict(
    rs, rf_model, type = "response", na.rm = TRUE,
    filename = filename, overwrite = overwrite,
    wopt = list(datatype = "FLT4S", gdal = gdal_wopt_compress())
  )
  invisible(filename)
}

predict_svm_raster_tiled_files <- function(r, svm_model, feature_names, tile_dir, tiles_x = 2L, tiles_y = 2L) {
  rs <- align_stack_for_model(r, feature_names)
  exl <- grid_tile_extents(rs, tiles_x, tiles_y)
  paths <- character(0)
  for (i in seq_along(exl)) {
    message("  SVM tile ", i, "/", length(exl))
    rc <- suppressWarnings(tryCatch(
      terra::crop(rs, exl[[i]], snap = "near"),
      error = function(e) NULL
    ))
    if (is.null(rc) || terra::ncell(rc) < 1L) next
    pr <- predict_svm_raster(rc, svm_model, feature_names)
    tp <- file.path(tile_dir, sprintf("svm_tile_%02d.tif", length(paths) + 1L))
    terra::writeRaster(
      pr, tp, overwrite = TRUE, datatype = "FLT4S",
      gdal = gdal_wopt_compress()
    )
    paths <- c(paths, tp)
  }
  if (!length(paths)) {
    stop("SVM tiling produced no predictions (empty AOI or crop failure).", call. = FALSE)
  }
  if (length(paths) == 1L) return(paths[1])
  vrtf <- file.path(tile_dir, "svm_merged.vrt")
  terra::vrt(paths, filename = vrtf, overwrite = TRUE)
  vrtf
}

predict_rf_raster <- function(r, rf_model, feature_names) {
  rs <- align_stack_for_model(r, feature_names)
  terra::predict(rs, rf_model, type = "response", na.rm = TRUE)
}

predict_svm_raster <- function(r, svm_model, feature_names) {
  rs <- align_stack_for_model(r, feature_names)
  v <- terra::values(rs, mat = TRUE, na.rm = FALSE)
  out <- rep(NA_real_, nrow(v))
  ok <- stats::complete.cases(v)
  if (any(ok)) {
    pr <- stats::predict(svm_model, v[ok, , drop = FALSE], type = "class")
    out[ok] <- as.numeric(pr)
  }
  rr <- rs[[1]] * NA
  terra::values(rr) <- out
  names(rr) <- "class"
  rr
}

focal_majority_3x3 <- function(r_cls) {
  # Integer classes; promote to float (use S3 as.numeric — not terra::as.numeric)
  r_cls <- as.numeric(r_cls)
  maj <- function(x, ...) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_real_)
    tb <- sort(table(x), decreasing = TRUE)
    v <- names(tb)[1L]
    if (!nzchar(v)) return(NA_real_)
    as.double(v)
  }
  # terra >= 1.9 passes extra args (e.g. na.rm) into fun — use ... ; w = 3 is odd 3×3 window
  z <- terra::focal(r_cls, w = 3, fun = maj, pad = TRUE)
  terra::round(z, digits = 0)
}

focal_majority_3x3_to_file <- function(r_cls, filename, overwrite = TRUE) {
  r_cls <- as.numeric(r_cls)
  maj <- function(x, ...) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_real_)
    tb <- sort(table(x), decreasing = TRUE)
    v <- names(tb)[1L]
    if (!nzchar(v)) return(NA_real_)
    as.double(v)
  }
  terra::focal(
    r_cls, w = 3, fun = maj, pad = TRUE,
    filename = filename, overwrite = overwrite,
    wopt = list(datatype = "FLT4S", gdal = gdal_wopt_compress())
  )
  terra::rast(filename)
}

load_norm_params <- function(cfg) {
  p <- file.path(cfg$paths$processed_dir, "normalization_params.csv")
  if (!file.exists(p)) {
    warning("No normalization_params.csv — skipping spectral rules.", call. = FALSE)
    return(NULL)
  }
  utils::read.csv(p, stringsAsFactors = FALSE)
}

denormalize_layer <- function(r_norm, year, layer_name, params) {
  if (is.null(params)) return(NULL)
  row <- params[params$year == year & params$layer == layer_name, , drop = FALSE]
  if (nrow(row) != 1L) return(NULL)
  mn <- row$min[1]
  mx <- row$max[1]
  r_norm * (mx - mn) + mn
}

apply_spectral_rules <- function(cls, feats_norm, year, cfg,
                                 mndwi_water = 0.3, ndvi_veg = 0.6) {
  params <- load_norm_params(cfg)
  ndvi <- denormalize_layer(feats_norm[["NDVI"]], year, "NDVI", params)
  mndwi <- denormalize_layer(feats_norm[["MNDWI"]], year, "MNDWI", params)
  if (is.null(ndvi) || is.null(mndwi)) return(cls)
  out <- cls
  out <- terra::ifel(mndwi > mndwi_water, 4L, out)
  out <- terra::ifel(ndvi > ndvi_veg, 2L, out)
  out
}

post_classify_raster <- function(cls, feats_norm, year, cfg, focal_outfile = NULL) {
  z <- if (!is.null(focal_outfile)) {
    focal_majority_3x3_to_file(cls, focal_outfile, overwrite = TRUE)
  } else {
    focal_majority_3x3(cls)
  }
  z <- terra::round(z)
  z <- terra::clamp(z, lower = 1, upper = 4)
  apply_spectral_rules(z, feats_norm, year, cfg)
}

# --------------------------------------------------------------------------- #
# Urban / class area statistics (km²)
# --------------------------------------------------------------------------- #

class_area_km2 <- function(r_lulc) {
  px_m2 <- prod(terra::res(r_lulc))
  f <- terra::freq(r_lulc, digits = 6)
  if (is.null(f) || nrow(f) == 0) {
    return(data.frame(class = integer(), km2 = numeric()))
  }
  if ("layer" %in% names(f)) f <- f[f$layer == 1L, , drop = FALSE]
  km2 <- f$count * px_m2 / 1e6
  data.frame(class = as.integer(f$value), km2 = km2)
}

urban_area_stats_row <- function(year, r_lulc) {
  areas <- class_area_km2(r_lulc)
  tot <- sum(areas$km2, na.rm = TRUE)
  u <- areas$km2[match(1L, areas$class)]
  v <- areas$km2[match(2L, areas$class)]
  b <- areas$km2[match(3L, areas$class)]
  w <- areas$km2[match(4L, areas$class)]
  data.frame(
    Year = year,
    Urban_km2 = u %??% NA_real_,
    Vegetation_km2 = v %??% NA_real_,
    BareSoil_km2 = b %??% NA_real_,
    Water_km2 = w %??% NA_real_,
    Total_km2 = tot,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------- #
# Train both sensors + save artifacts
# --------------------------------------------------------------------------- #

#' @param training_path Path to `sf` file (gpkg/shp) with `class` column 1–4.
#' @param tune_svm If TRUE, runs e1071::tune (slow).
#' @param svm_cost,svm_gamma Used when `tune_svm` is FALSE.
train_phase3_models <- function(
    cfg = load_config(),
    training_path = NULL,
    class_col = "class",
    p_train = 0.7,
    seed = 42L,
    tune_svm = FALSE,
    svm_cost = 1,
    svm_gamma = NULL,
    min_per_class = 100L) {
  proc <- cfg$paths$processed_dir
  dir.create(proc, showWarnings = FALSE, recursive = TRUE)
  if (is.null(training_path)) training_path <- path_training_default(cfg)
  if (!file.exists(training_path)) {
    at <- get_auto_training_cfg(cfg)
    if (isTRUE(at$enable)) {
      message("No training file; generating samples via classification.auto_training …")
      auto_build_training_samples(cfg, out_path = training_path)
    } else {
      stop(
        "No training file: ", training_path,
        "\nSet classification.auto_training.enable: true (default) for WorldCover/GLOBE samples, ",
        "or run collect_training_mapedit() / add your own sf with `class` 1–4.",
        call. = FALSE
      )
    }
  }

  samples <- sf::st_read(training_path, quiet = TRUE)
  ys <- sentinel_training_year(cfg)
  yl <- landsat_training_year(cfg)
  fs <- path_features(ys, cfg)
  fl <- path_features(yl, cfg)
  if (!file.exists(fs)) stop("Missing Sentinel features: ", fs, call. = FALSE)
  if (!file.exists(fl)) stop("Missing Landsat features: ", fl, call. = FALSE)

  rs <- terra::rast(fs)
  rl <- terra::rast(fl)
  dfs <- extract_training_table(samples, rs, class_col)
  dfl <- extract_training_table(samples, rl, class_col)
  check_min_samples(dfs, "class", min_per_class)
  check_min_samples(dfl, "class", min_per_class)

  sp <- stratified_train_test_paired(dfs, dfl, "class", p_train, seed)

  # --- Sentinel RF + SVM ---
  rf_s <- train_rf_lulc(sp$train_a, "class")
  svm_s <- train_svm_lulc(sp$train_a, "class", cost = svm_cost, gamma = svm_gamma, tune = tune_svm)

  ev_rf_s <- predict_test_aligned(rf_s, sp$test_a, "class", "rf")
  ev_sv_s <- predict_test_aligned(svm_s, sp$test_a, "class", "svm")
  m_rf_s <- confusion_metrics_df(ev_rf_s$pred, ev_rf_s$ref, "RandomForest", "Sentinel-2")
  m_sv_s <- confusion_metrics_df(ev_sv_s$pred, ev_sv_s$ref, "SVM", "Sentinel-2")

  # --- Landsat RF + SVM ---
  rf_l <- train_rf_lulc(sp$train_b, "class")
  svm_l <- train_svm_lulc(sp$train_b, "class", cost = svm_cost, gamma = svm_gamma, tune = tune_svm)

  ev_rf_l <- predict_test_aligned(rf_l, sp$test_b, "class", "rf")
  ev_sv_l <- predict_test_aligned(svm_l, sp$test_b, "class", "svm")
  m_rf_l <- confusion_metrics_df(ev_rf_l$pred, ev_rf_l$ref, "RandomForest", "Landsat")
  m_sv_l <- confusion_metrics_df(ev_sv_l$pred, ev_sv_l$ref, "SVM", "Landsat")

  acc <- rbind(m_rf_s$long, m_sv_s$long, m_rf_l$long, m_sv_l$long)
  utils::write.csv(acc, file.path(proc, "accuracy_report.csv"), row.names = FALSE)

  imp <- rbind(
    rf_importance_df(rf_s, "Sentinel-2"),
    rf_importance_df(rf_l, "Landsat")
  )
  utils::write.csv(imp, file.path(proc, "variable_importance.csv"), row.names = FALSE)

  saveRDS(
    list(model = rf_s, feature_names = names(rs), sensor = "sentinel", train_year = ys),
    file.path(proc, "rf_model_sentinel.rds")
  )
  saveRDS(
    list(model = svm_s, feature_names = names(rs), sensor = "sentinel", train_year = ys),
    file.path(proc, "svm_model_sentinel.rds")
  )
  saveRDS(
    list(model = rf_l, feature_names = names(rl), sensor = "landsat", train_year = yl),
    file.path(proc, "rf_model_landsat.rds")
  )
  saveRDS(
    list(model = svm_l, feature_names = names(rl), sensor = "landsat", train_year = yl),
    file.path(proc, "svm_model_landsat.rds")
  )

  message("Wrote models + accuracy_report.csv + variable_importance.csv to ", proc)
  invisible(list(
    rf_sentinel = rf_s, svm_sentinel = svm_s,
    rf_landsat = rf_l, svm_landsat = svm_l,
    accuracy = acc
  ))
}

# --------------------------------------------------------------------------- #
# Predict all years (choose best model by Kappa on test — stored in accuracy table)
# --------------------------------------------------------------------------- #

pick_best_model_name <- function(acc_csv_path, sensor) {
  acc <- utils::read.csv(acc_csv_path, stringsAsFactors = FALSE)
  sub <- acc[acc$sensor == sensor & acc$metric == "Kappa", , drop = FALSE]
  if (!nrow(sub)) return("RandomForest")
  wm <- which.max(sub$value)
  sub$model[wm]
}

predict_phase3_all_years <- function(
    cfg = load_config(),
    model_choice = c("auto", "rf", "svm"),
    post_process = TRUE,
    memory_safe = NULL,
    years = NULL,
    temp_root = NULL,
    cleanup_temp = NULL) {
  proc <- cfg$paths$processed_dir
  acc_path <- file.path(proc, "accuracy_report.csv")
  if (!file.exists(acc_path)) stop("Run train_phase3_models() first.", call. = FALSE)

  pcfg <- get_phase3_prediction_cfg(cfg)
  if (is.null(memory_safe)) memory_safe <- pcfg$memory_safe
  if (is.null(cleanup_temp)) cleanup_temp <- pcfg$cleanup_temp

  mf_saved <- NULL
  if (isTRUE(memory_safe) && is.finite(pcfg$memfrac)) {
    mf_saved <- tryCatch(terra::terraOptions("memfrac"), error = function(e) NULL)
    terra::terraOptions(memfrac = pcfg$memfrac)
    on.exit(
      {
        if (!is.null(mf_saved)) {
          try(terra::terraOptions(memfrac = mf_saved), silent = TRUE)
        }
      },
      add = FALSE
    )
  }

  model_choice <- match.arg(model_choice)
  years_all <- sort(unique(phase3_years(cfg)))
  if (is.null(years)) {
    years <- years_all
  } else {
    years <- sort(unique(as.integer(years)))
    years <- intersect(years, years_all)
    if (!length(years)) {
      stop("No valid years to predict (check years vs phase3_years(cfg)).", call. = FALSE)
    }
  }

  tmp_top <- temp_root %??% tempdir()
  merge_urban_stats <- !is.null(years)
  stats_rows <- list()

  for (yr in years) {
    fp <- path_features(yr, cfg)
    if (!file.exists(fp)) {
      warning("Skip year ", yr, " — missing ", fp)
      next
    }
    message("Phase 3 predict ", yr, if (isTRUE(memory_safe)) " (memory-safe path)" else "")
    r <- terra::rast(fp)
    sens <- sensor_from_stack(r)

    mc <- model_choice
    if (mc == "auto") {
      mc <- pick_best_model_name(acc_path, if (sens == "sentinel") "Sentinel-2" else "Landsat")
      mc <- if (grepl("Random", mc, ignore.case = TRUE)) "rf" else "svm"
    }

    rds <- file.path(
      proc,
      if (sens == "sentinel") {
        if (mc == "rf") "rf_model_sentinel.rds" else "svm_model_sentinel.rds"
      } else {
        if (mc == "rf") "rf_model_landsat.rds" else "svm_model_landsat.rds"
      }
    )
    pack <- readRDS(rds)
    fn <- pack$feature_names

    td <- NULL
    td2 <- NULL
    if (isTRUE(memory_safe)) {
      td <- tempfile(pattern = paste0("p3_", yr, "_"), tmpdir = tmp_top)
      dir.create(td, showWarnings = FALSE, recursive = TRUE)
    }

    if (isTRUE(memory_safe)) {
      if (mc == "rf") {
        pred_tif <- file.path(td, "pred_rf.tif")
        predict_rf_raster_to_file(r, pack$model, fn, pred_tif)
        cls <- terra::rast(pred_tif)
      } else {
        pred_src <- predict_svm_raster_tiled_files(
          r, pack$model, fn, td, pcfg$tiles_x, pcfg$tiles_y
        )
        cls <- terra::rast(pred_src)
      }
    } else {
      cls <- if (mc == "rf") {
        predict_rf_raster(r, pack$model, fn)
      } else {
        predict_svm_raster(r, pack$model, fn)
      }
    }

    names(cls) <- "class"
    cls <- terra::round(cls)
    cls <- terra::clamp(cls, lower = 1, upper = 4)

    if (post_process) {
      if (isTRUE(memory_safe) && isTRUE(pcfg$skip_focal)) {
        cls <- apply_spectral_rules(cls, r, yr, cfg)
      } else if (isTRUE(memory_safe) && !isTRUE(pcfg$skip_focal)) {
        td2 <- tempfile(pattern = paste0("p3f_", yr, "_"), tmpdir = tmp_top)
        dir.create(td2, showWarnings = FALSE, recursive = TRUE)
        foc_path <- file.path(td2, "focal.tif")
        cls <- post_classify_raster(cls, r, yr, cfg, focal_outfile = foc_path)
      } else {
        cls <- post_classify_raster(cls, r, yr, cfg)
      }
      cls <- terra::round(cls)
      cls <- terra::clamp(cls, lower = 1, upper = 4)
    }

    out <- path_lulc(yr, cfg)
    terra::writeRaster(
      cls,
      out,
      overwrite = TRUE,
      datatype = "INT1U",
      gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
    )
    message("Wrote ", out, " (", sens, ", ", mc, ")")

    if (isTRUE(cleanup_temp)) {
      if (!is.null(td) && dir.exists(td)) unlink(td, recursive = TRUE)
      if (!is.null(td2) && dir.exists(td2)) unlink(td2, recursive = TRUE)
    }

    stats_rows[[length(stats_rows) + 1L]] <- urban_area_stats_row(yr, terra::rast(out))
    if (isTRUE(memory_safe)) gc(verbose = FALSE)
  }

  if (!length(stats_rows)) {
    warning("Phase 3 prediction: no rasters written (missing feature GeoTIFFs or no valid years?).")
    return(invisible(NULL))
  }

  tab <- do.call(rbind, stats_rows)
  out_stats <- file.path(proc, "urban_area_stats.csv")
  if (isTRUE(merge_urban_stats) && file.exists(out_stats)) {
    old <- utils::read.csv(out_stats, stringsAsFactors = FALSE)
    old <- old[!old$Year %in% years, , drop = FALSE]
    tab <- rbind(old, tab)
    tab <- tab[order(tab$Year), ]
  }
  utils::write.csv(tab, out_stats, row.names = FALSE)
  message("Wrote ", out_stats)
  invisible(tab)
}

#' Predict a single year (same options as `predict_phase3_all_years`, `years` fixed).
predict_phase3_year <- function(
    year,
    cfg = load_config(),
    model_choice = c("auto", "rf", "svm"),
    post_process = TRUE,
    memory_safe = NULL,
    temp_root = NULL,
    cleanup_temp = NULL) {
  predict_phase3_all_years(
    cfg = cfg,
    model_choice = model_choice,
    post_process = post_process,
    memory_safe = memory_safe,
    years = as.integer(year),
    temp_root = temp_root,
    cleanup_temp = cleanup_temp
  )
}

#' Full Phase 3: train (if models missing) + optional predict
#'
#' @param predict If FALSE, only training (or skip if models exist); use `predict_phase3_all_years()` later.
#' @param years Integer vector of years to predict; NULL = all phase-3 years.
#' @param memory_safe NULL = use `classification.prediction` in config (default TRUE). Uses disk-backed RF
#'   predict, tiled SVM + VRT, and focal-to-file to reduce RAM spikes.
run_phase3 <- function(
    cfg = load_config(),
    training_path = NULL,
    tune_svm = FALSE,
    model_choice = c("auto", "rf", "svm"),
    skip_train_if_models_exist = TRUE,
    predict = TRUE,
    years = NULL,
    memory_safe = NULL,
    temp_root = NULL,
    cleanup_temp = NULL) {
  with_pipeline_phase("Phase 3 — LULC training and/or prediction (can be slow)", {
  proc <- cfg$paths$processed_dir
  has_rf_s <- file.exists(file.path(proc, "rf_model_sentinel.rds"))
  if (!has_rf_s || !skip_train_if_models_exist) {
    pipeline_log("Phase 3: training models …")
    train_phase3_models(cfg, training_path = training_path, tune_svm = tune_svm)
  } else {
    pipeline_log("Phase 3: using existing models in ", proc, " (skip_train_if_models_exist=FALSE to retrain).")
  }
  if (isTRUE(predict)) {
    pipeline_log("Phase 3: predicting LULC for configured years …")
    predict_phase3_all_years(
      cfg,
      model_choice = model_choice,
      memory_safe = memory_safe,
      years = years,
      temp_root = temp_root,
      cleanup_temp = cleanup_temp
    )
  } else {
    pipeline_log("Phase 3: skipping prediction (predict=FALSE).")
    message("Skipping prediction (predict=FALSE). Run predict_phase3_all_years() or predict_phase3_year() when ready.")
    invisible(NULL)
  }
  })
}

if (interactive()) {
  message(
    "Phase 3: run_phase3() — training + prediction. ",
    "Use predict=FALSE to train only; predict_phase3_year(2020) for one year; ",
    "classification.prediction in study_area.yml for memory_safe / tiles / skip_focal."
  )
}
