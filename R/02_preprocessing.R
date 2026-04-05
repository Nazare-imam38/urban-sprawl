#' Phase 2 — Cloud masking, reflectance scaling, spectral indices, min–max normalization
#'
#' Inputs: Phase 1 GeoTIFFs (`data/raw/sentinel2/S2_composite_*.tif`, `data/raw/landsat/LS_composite_*.tif`).
#' Outputs: `data/processed/features_{year}.tif` (UTM), `data/processed/normalization_params.csv`.
#'
#' Sentinel-2: SCL mask → scale /1e4 → indices → 11 layers (6 SR + 5 indices).
#' Landsat: QA_PIXEL bits 1,3,4,5 → scale ×0.0000275−0.2 → indices → 11 layers (6 SR + 5 indices).
#'
#' Working directory = project root.

source("R/00_setup.R")

EPS <- 1e-6

# SCL: keep 4–7; mask 0–3, 8–11 (per spec)
SCL_KEEP <- c(4L, 5L, 6L, 7L)

# QA_PIXEL: mask if bits 1 (dilated cloud), 3 (cloud), 4 (cloud shadow), 5 (snow) set
QA_MASK_BITS <- bitwOr(bitwOr(2L, 8L), bitwOr(16L, 32L)) # 58

.safe_div <- function(num, den) num / (den + EPS)

# --------------------------------------------------------------------------- #
# Load Phase 1 rasters
# --------------------------------------------------------------------------- #

path_s2_composite <- function(year, cfg) {
  file.path(cfg$paths$raw_dir, "sentinel2", sprintf("S2_composite_%d.tif", year))
}

path_ls_composite <- function(year, cfg) {
  file.path(cfg$paths$raw_dir, "landsat", sprintf("LS_composite_%d.tif", year))
}

load_s2_phase1 <- function(year, cfg) {
  p <- path_s2_composite(year, cfg)
  if (!file.exists(p)) stop("Missing Phase 1 file: ", p, call. = FALSE)
  r <- terra::rast(p)
  exp <- c("B02", "B03", "B04", "B08", "B11", "B12", "SCL")
  if (!identical(names(r), exp)) {
    if (terra::nlyr(r) != length(exp)) {
      stop("S2 stack for ", year, " must have 7 layers; got ", terra::nlyr(r), call. = FALSE)
    }
    names(r) <- exp
    warning("S2 layer names were reassigned to: ", paste(exp, collapse = ", "))
  }
  r
}

load_ls_phase1 <- function(year, cfg) {
  p <- path_ls_composite(year, cfg)
  if (!file.exists(p)) stop("Missing Phase 1 file: ", p, call. = FALSE)
  r <- terra::rast(p)
  exp <- c("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL")
  if (!identical(names(r), exp)) {
    if (terra::nlyr(r) != length(exp)) {
      stop("Landsat stack for ", year, " must have 7 layers; got ", terra::nlyr(r), call. = FALSE)
    }
    names(r) <- exp
    warning("Landsat layer names were reassigned to: ", paste(exp, collapse = ", "))
  }
  r
}

# --------------------------------------------------------------------------- #
# Cloud / quality masking
# --------------------------------------------------------------------------- #

mask_from_scl <- function(scl) {
  terra::ifel(scl %in% SCL_KEEP, 1L, NA_integer_)
}

landsat_qa_mask <- function(qa) {
  terra::app(qa, fun = function(x, ...) {
    qi <- suppressWarnings(as.integer(round(x)))
    bad <- is.na(qi) | (bitwAnd(qi, QA_MASK_BITS) > 0L)
    ifelse(bad, NA_real_, 1)
  })
}

mask_sentinel2 <- function(r) {
  scl <- r[["SCL"]]
  m <- mask_from_scl(scl)
  refl <- r[[c("B02", "B03", "B04", "B08", "B11", "B12")]]
  terra::mask(refl, m)
}

mask_landsat <- function(r) {
  qa <- r[["QA_PIXEL"]]
  m <- landsat_qa_mask(qa)
  sr <- r[[c("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7")]]
  terra::mask(sr, m)
}

# --------------------------------------------------------------------------- #
# Reflectance scaling
# --------------------------------------------------------------------------- #

scale_s2_reflectance <- function(r_refl) {
  x <- r_refl / 10000
  terra::clamp(x, lower = 0, upper = 1, values = TRUE)
}

scale_landsat_reflectance <- function(r_sr) {
  x <- r_sr * 0.0000275 - 0.2
  terra::clamp(x, lower = 0, upper = 1, values = TRUE)
}

# --------------------------------------------------------------------------- #
# Spectral indices (scaled reflectance)
# --------------------------------------------------------------------------- #

indices_from_s2 <- function(r) {
  B02 <- r[["B02"]]
  B03 <- r[["B03"]]
  B04 <- r[["B04"]]
  B08 <- r[["B08"]]
  B11 <- r[["B11"]]
  NDVI <- .safe_div(B08 - B04, B08 + B04)
  NDBI <- .safe_div(B11 - B08, B11 + B08)
  MNDWI <- .safe_div(B03 - B11, B03 + B11)
  NBI <- .safe_div(B04 * B11, B08)
  den_bsi <- (B11 + B04) + (B08 + B02)
  BSI <- .safe_div((B11 + B04) - (B08 + B02), den_bsi)
  list(NDVI = NDVI, NDBI = NDBI, MNDWI = MNDWI, NBI = NBI, BSI = BSI)
}

indices_from_landsat <- function(r) {
  B2 <- r[["SR_B2"]]
  B3 <- r[["SR_B3"]]
  B4 <- r[["SR_B4"]]
  B5 <- r[["SR_B5"]]
  B6 <- r[["SR_B6"]]
  NDVI <- .safe_div(B5 - B4, B5 + B4)
  NDBI <- .safe_div(B6 - B5, B6 + B5)
  MNDWI <- .safe_div(B3 - B6, B3 + B6)
  NBI <- .safe_div(B4 * B6, B5)
  den_bsi <- (B6 + B4) - (B5 + B2)
  den_bsi2 <- (B6 + B4) + (B5 + B2)
  BSI <- .safe_div(den_bsi, den_bsi2)
  list(NDVI = NDVI, NDBI = NDBI, MNDWI = MNDWI, NBI = NBI, BSI = BSI)
}

# --------------------------------------------------------------------------- #
# Feature stacks (reflectance + indices)
# --------------------------------------------------------------------------- #

build_feature_stack_s2 <- function(r_scaled) {
  idx <- indices_from_s2(r_scaled)
  nm_refl <- c("B02", "B03", "B04", "B08", "B11", "B12")
  refl <- r_scaled[[nm_refl]]
  terra::rast(c(
    list(refl),
    list(idx$NDVI, idx$NDBI, idx$MNDWI, idx$NBI, idx$BSI)
  ))
}

build_feature_stack_landsat <- function(r_scaled) {
  idx <- indices_from_landsat(r_scaled)
  nm_sr <- c("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7")
  sr <- r_scaled[[nm_sr]]
  terra::rast(c(
    list(sr),
    list(idx$NDVI, idx$NDBI, idx$MNDWI, idx$NBI, idx$BSI)
  ))
}

feature_names_s2 <- function() {
  c("B02", "B03", "B04", "B08", "B11", "B12", "NDVI", "NDBI", "MNDWI", "NBI", "BSI")
}

feature_names_landsat <- function() {
  c("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "NDVI", "NDBI", "MNDWI", "NBI", "BSI")
}

# --------------------------------------------------------------------------- #
# Reproject to target CRS (UTM)
# --------------------------------------------------------------------------- #

project_to_study_crs <- function(r, crs_target) {
  ct <- terra::crs(r)
  if (is.na(ct) || ct == "") stop("Raster has no CRS", call. = FALSE)
  if (terra::same.crs(r, crs_target)) return(r)
  terra::project(r, crs_target, method = "bilinear", gdal = TRUE)
}

# --------------------------------------------------------------------------- #
# Min–max normalization (per layer); record min/max before normalize
# --------------------------------------------------------------------------- #

minmax_normalize_stack <- function(r, year, sensor, params_accum) {
  nms <- names(r)
  if (is.null(nms) || any(nzchar(nms) == FALSE)) {
    nms <- paste0("L", seq_len(terra::nlyr(r)))
    names(r) <- nms
  }
  mins <- terra::global(r, fun = "min", na.rm = TRUE)[, 1]
  maxs <- terra::global(r, fun = "max", na.rm = TRUE)[, 1]
  out_layers <- vector("list", terra::nlyr(r))
  for (i in seq_len(terra::nlyr(r))) {
    mn <- mins[i]
    mx <- maxs[i]
    ri <- r[[i]]
    if (!is.finite(mn) || !is.finite(mx) || abs(mx - mn) < EPS) {
      out_layers[[i]] <- ri * 0
      warning("Layer ", nms[i], " (", year, "): degenerate min/max; filled with 0.")
    } else {
      out_layers[[i]] <- (ri - mn) / (mx - mn)
    }
  }
  rout <- terra::rast(out_layers)
  names(rout) <- nms
  rows <- data.frame(
    year = year,
    sensor = sensor,
    layer = nms,
    min = mins,
    max = maxs,
    stringsAsFactors = FALSE
  )
  params_accum <- rbind(params_accum, rows)
  list(r = rout, params = params_accum)
}

# --------------------------------------------------------------------------- #
# Quality checks (§2.3)
# --------------------------------------------------------------------------- #

phase2_quality_checks <- function(r_features, r_ndvi, r_ndbi, label = "") {
  chk <- list(label = label)
  # terra >= 1.9 forwards na.rm into fun(); accept ... so custom fun is valid
  chk$na_frac_layers <- terra::global(r_features, fun = function(x, ...) mean(is.na(x)), na.rm = FALSE)
  mm_ndvi <- terra::minmax(r_ndvi)
  mm_ndbi <- terra::minmax(r_ndbi)
  chk$ndvi_minmax <- mm_ndvi
  chk$ndbi_minmax <- mm_ndbi
  refl_min <- terra::global(r_features[[1:min(6L, terra::nlyr(r_features))]], fun = "min", na.rm = TRUE)
  chk$refl_mins <- refl_min
  geom_ok <- TRUE
  if (terra::nlyr(r_features) > 1) {
    for (j in 2:terra::nlyr(r_features)) {
      if (!terra::compareGeom(r_features[[1]], r_features[[j]], stopOnError = FALSE)) {
        geom_ok <- FALSE
        break
      }
    }
  }
  chk$compareGeom_all <- geom_ok
  chk
}

print_phase2_checks <- function(chk) {
  message("--- Quality checks: ", chk$label, " ---")
  message("NDVI min/max: ", paste(round(chk$ndvi_minmax, 4), collapse = ", "))
  message("NDBI min/max: ", paste(round(chk$ndbi_minmax, 4), collapse = ", "))
  message("compareGeom all layers: ", chk$compareGeom_all)
  invisible(chk)
}

# --------------------------------------------------------------------------- #
# End-to-end per year
# --------------------------------------------------------------------------- #

process_sentinel2_year_phase2 <- function(year, cfg, crs_target, params_accum, quality_check = FALSE) {
  r <- load_s2_phase1(year, cfg)
  r_m <- mask_sentinel2(r)
  r_s <- scale_s2_reflectance(r_m)
  feat <- build_feature_stack_s2(r_s)
  names(feat) <- feature_names_s2()
  feat <- project_to_study_crs(feat, crs_target)
  if (quality_check) {
    i_ndvi <- which(names(feat) == "NDVI")
    i_ndbi <- which(names(feat) == "NDBI")
    chk <- phase2_quality_checks(feat, feat[[i_ndvi]], feat[[i_ndbi]], paste0("S2 ", year))
    print_phase2_checks(chk)
  }
  minmax_normalize_stack(feat, year, "Sentinel-2", params_accum)
}

process_landsat_year_phase2 <- function(year, cfg, crs_target, params_accum, quality_check = FALSE) {
  r <- load_ls_phase1(year, cfg)
  r_m <- mask_landsat(r)
  r_s <- scale_landsat_reflectance(r_m)
  feat <- build_feature_stack_landsat(r_s)
  names(feat) <- feature_names_landsat()
  feat <- project_to_study_crs(feat, crs_target)
  if (quality_check) {
    i_ndvi <- which(names(feat) == "NDVI")
    i_ndbi <- which(names(feat) == "NDBI")
    chk <- phase2_quality_checks(feat, feat[[i_ndvi]], feat[[i_ndbi]], paste0("Landsat ", year))
    print_phase2_checks(chk)
  }
  minmax_normalize_stack(feat, year, "Landsat", params_accum)
}

# --------------------------------------------------------------------------- #
# Run Phase 2 for all configured years
# --------------------------------------------------------------------------- #

#' @param quality_check If TRUE, print NDVI/NDBI minmax and geometry checks
#' @param overwrite If FALSE, skip years whose `features_{year}.tif` already exists
run_phase2 <- function(cfg = load_config(), quality_check = FALSE, overwrite = FALSE) {
  with_pipeline_phase("Phase 2 — preprocessing (features rasters per year)", {
  proc <- cfg$paths$processed_dir
  dir.create(proc, recursive = TRUE, showWarnings = FALSE)
  crs_t <- cfg$study_area$crs_projected
  csv_path <- file.path(proc, "normalization_params.csv")

  pc <- cfg$planetary_computer
  if (is.null(pc)) stop("config: planetary_computer block missing (years for Phase 2).", call. = FALSE)

  params_accum <- if (file.exists(csv_path)) {
    utils::read.csv(csv_path, stringsAsFactors = FALSE)
  } else {
    data.frame(
      year = integer(),
      sensor = character(),
      layer = character(),
      min = numeric(),
      max = numeric(),
      stringsAsFactors = FALSE
    )
  }

  for (yr in unlist(pc$sentinel2$years, use.names = FALSE)) {
    yr <- as.integer(yr)
    outf <- file.path(proc, sprintf("features_%d.tif", yr))
    if (file.exists(outf) && !overwrite) {
      pipeline_log("Phase 2: skip existing ", outf)
      next
    }
    params_accum <- params_accum[params_accum$year != yr, , drop = FALSE]
    if (!file.exists(path_s2_composite(yr, cfg))) {
      warning("No Phase 1 Sentinel-2 for ", yr, "; skipping.")
      next
    }
    pipeline_log("Phase 2: Sentinel-2 ", yr, " (compute + write) …")
    res <- process_sentinel2_year_phase2(yr, cfg, crs_t, params_accum, quality_check = quality_check)
    params_accum <- res$params
    terra::writeRaster(
      res$r,
      outf,
      overwrite = TRUE,
      gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
    )
    message("Wrote ", outf)
  }

  for (row in seq_along(pc$landsat$acquisitions)) {
    acq <- pc$landsat$acquisitions[[row]]
    yr <- as.integer(acq$year)
    outf <- file.path(proc, sprintf("features_%d.tif", yr))
    if (file.exists(outf) && !overwrite) {
      pipeline_log("Phase 2: skip existing ", outf)
      next
    }
    params_accum <- params_accum[params_accum$year != yr, , drop = FALSE]
    if (!file.exists(path_ls_composite(yr, cfg))) {
      warning("No Phase 1 Landsat for ", yr, "; skipping.")
      next
    }
    pipeline_log("Phase 2: Landsat ", yr, " (compute + write) …")
    res <- process_landsat_year_phase2(yr, cfg, crs_t, params_accum, quality_check = quality_check)
    params_accum <- res$params
    terra::writeRaster(
      res$r,
      outf,
      overwrite = TRUE,
      gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
    )
    message("Wrote ", outf)
  }

  utils::write.csv(params_accum, csv_path, row.names = FALSE)
  pipeline_log("Phase 2: wrote ", csv_path, " (", nrow(params_accum), " rows)")
  invisible(params_accum)
  })
}

#' Quick maps for §2.4 — NDVI, NDBI, MNDWI (normalized features raster)
plot_phase2_sanity_maps <- function(year, cfg = load_config()) {
  f <- file.path(cfg$paths$processed_dir, sprintf("features_%d.tif", year))
  if (!file.exists(f)) stop("Run Phase 2 first: ", f, call. = FALSE)
  r <- terra::rast(f)
  nm <- names(r)
  par(mfrow = c(1, 3))
  for (n in c("NDVI", "NDBI", "MNDWI")) {
    i <- match(n, nm)
    if (!is.na(i)) plot(r[[i]], main = paste(year, n))
  }
  invisible(r)
}

if (interactive()) {
  message("Phase 2 loaded. Run: run_phase2(quality_check = TRUE)")
}
