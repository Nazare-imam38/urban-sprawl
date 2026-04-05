#' Phase 1 — Data acquisition (Microsoft Planetary Computer)
#'
#' STAC search → sign URLs → stream COGs with `/vsicurl/` → crop/mask AOI → stack → GeoTIFF.
#' Sentinel-2: one best scene per year (lowest cloud), Jan–Mar; 20m bands resampled to 10m.
#' Landsat C2 L2: one scene per year at native 30m (PC assets: blue, green, red, nir08, …).
#' GlobeLand30: merge user-downloaded tiles, crop, project, save (manual download).
#'
#' Working directory = project root. Requires: rstac, terra, sf, httr, jsonlite, yaml.

source("R/00_setup.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

suppressPackageStartupMessages({
  req <- c("rstac", "httr", "jsonlite")
  miss <- req[!vapply(req, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(miss)) {
    stop("Install Phase 1 packages: ", paste(miss, collapse = ", "),
         "\ninstall.packages(c(\"rstac\",\"httr\",\"jsonlite\"))", call. = FALSE)
  }
  library(rstac)
  library(httr)
  library(jsonlite)
})

# --------------------------------------------------------------------------- #
# Config helpers
# --------------------------------------------------------------------------- #

get_pc_cfg <- function(cfg) {
  pc <- cfg$planetary_computer
  if (is.null(pc)) stop("config: add planetary_computer section to study_area.yml", call. = FALSE)
  pc
}

aoi_bbox_wgs84 <- function(cfg) {
  sa <- cfg$study_area
  c(sa$xmin, sa$ymin, sa$xmax, sa$ymax)
}

aoi_sf_polygon <- function(cfg) {
  sa <- cfg$study_area
  bb <- sf::st_bbox(
    c(xmin = sa$xmin, ymin = sa$ymin, xmax = sa$xmax, ymax = sa$ymax),
    crs = sa$crs_wgs84
  )
  sf::st_as_sf(sf::st_as_sfc(bb))
}

season_datetime_range <- function(year, month_start, month_end) {
  start <- sprintf("%d-%02d-01", year, month_start)
  mdays <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L)
  ld <- mdays[month_end]
  if (month_end == 2L) {
    leap <- (year %% 4L == 0L && year %% 100L != 0L) || (year %% 400L == 0L)
    if (leap) ld <- 29L
  }
  end <- sprintf("%d-%02d-%02d", year, month_end, ld)
  paste(start, end, sep = "/")
}

# --------------------------------------------------------------------------- #
# Planetary Computer: sign href (required for download)
# --------------------------------------------------------------------------- #

sign_planetary_href <- function(href, sign_endpoint) {
  if (!nzchar(href)) stop("Empty asset href", call. = FALSE)
  # R >= 4.5: URLencode() no longer has useUTF8= (RFC 3986–style encoding by default)
  q <- paste0(sign_endpoint, "?href=", utils::URLencode(href, reserved = TRUE))
  r <- httr::GET(q)
  if (httr::status_code(r) >= 400) {
    stop("Sign failed (", httr::status_code(r), "): ", substr(httr::content(r, "text", encoding = "UTF-8"), 1, 200),
         call. = FALSE)
  }
  txt <- httr::content(r, "text", encoding = "UTF-8")
  j <- jsonlite::fromJSON(txt)
  out <- j$href
  if (is.null(out) || !nzchar(out)) stop("Sign response missing href", call. = FALSE)
  out
}

rast_from_signed_href <- function(href, sign_endpoint) {
  signed <- sign_planetary_href(href, sign_endpoint)
  u <- paste0("/vsicurl/", signed)
  terra::rast(u)
}

item_cloud_cover <- function(feature) {
  p <- feature$properties
  if (is.null(p)) return(NA_real_)
  v <- p[["eo:cloud_cover"]]
  if (is.null(v)) v <- p$`eo:cloud_cover`
  land <- p[["landsat:cloud_cover_land"]]
  if (is.null(land)) land <- p$`landsat:cloud_cover_land`
  nv <- if (is.null(v) || !length(v)) NA_real_ else suppressWarnings(as.numeric(v)[1])
  nl <- if (is.null(land) || !length(land)) NA_real_ else suppressWarnings(as.numeric(land)[1])
  both <- c(nv, nl)
  both <- both[is.finite(both)]
  if (!length(both)) NA_real_ else min(both)
}

item_platform <- function(feature) {
  p <- feature$properties
  if (is.null(p) || is.null(p$platform)) return(NA_character_)
  as.character(p$platform)[1]
}

normalize_platform <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("_", "-", x, fixed = TRUE)
}

pick_best_feature <- function(features, platform_want = NULL, cloud_max = Inf) {
  if (!length(features)) return(NULL)
  rows <- lapply(features, function(f) {
    data.frame(
      cloud = item_cloud_cover(f),
      platform = item_platform(f),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  want <- if (is.null(platform_want) || is.na(platform_want)) NA_character_ else normalize_platform(platform_want)
  ok <- is.na(tab$cloud) | tab$cloud <= cloud_max
  if (!is.na(want)) {
    plat <- vapply(tab$platform, function(z) {
      if (is.na(z)) return(NA_character_)
      normalize_platform(z)
    }, character(1))
    ok <- ok & !is.na(plat) & plat == want
  }
  if (!any(ok)) return(NULL)
  cl <- tab$cloud
  cl[is.na(cl)] <- Inf
  ii <- which(ok)
  idx <- ii[order(cl[ii])[1]]
  features[[idx]]
}

stac_search_features <- function(stac_url, collection, bbox, datetime, limit = 500) {
  s <- rstac::stac(stac_url)
  req <- rstac::stac_search(
    s,
    collections = collection,
    bbox = bbox,
    datetime = datetime,
    limit = as.integer(limit)
  )
  doc <- rstac::get_request(req)
  feats <- doc$features
  if (is.null(feats)) {
    warning("STAC response has no features list; check bbox/datetime.")
    return(list())
  }
  feats
}

asset_href <- function(feature, key) {
  a <- feature$assets[[key]]
  if (is.null(a)) return(NA_character_)
  h <- a$href
  if (is.null(h)) NA_character_ else as.character(h)[1]
}

crop_mask_to_aoi <- function(r, aoi_sf_wgs84) {
  crs_r <- terra::crs(r)
  if (is.na(crs_r) || crs_r == "") stop("Raster has no CRS", call. = FALSE)
  aoi_p <- sf::st_transform(aoi_sf_wgs84, crs_r)
  v <- terra::vect(aoi_p)
  x <- terra::crop(r, v)
  terra::mask(x, v)
}

# --------------------------------------------------------------------------- #
# Sentinel-2: one scene per year → 7 bands @ 10m, GeoTIFF
# --------------------------------------------------------------------------- #

download_sentinel2_year <- function(
    year,
    cfg,
    overwrite = FALSE) {
  raw <- file.path(cfg$paths$raw_dir, "sentinel2")
  dir.create(raw, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(raw, sprintf("S2_composite_%d.tif", year))
  if (file.exists(out_path) && !overwrite) {
    message("Exists, loading: ", out_path)
    return(terra::rast(out_path))
  }

  pc <- get_pc_cfg(cfg)
  s2 <- pc$sentinel2
  bbox <- aoi_bbox_wgs84(cfg)
  aoi <- aoi_sf_polygon(cfg)
  sign_ep <- pc$sign_endpoint

  ms <- s2$season_month_start
  me <- s2$season_month_end
  dt <- season_datetime_range(year, ms, me)

  feats <- stac_search_features(pc$stac_url, s2$collection, bbox, dt, limit = s2$stac_limit)
  cloud_try <- as.numeric(c(s2$cloud_cover_max_percent, s2$cloud_cover_relaxed_percent))
  best <- NULL
  for (cm in cloud_try) {
    best <- pick_best_feature(feats, platform_want = NULL, cloud_max = cm)
    if (!is.null(best)) break
    message("No Sentinel-2 scenes with cloud <= ", cm, "% — trying next threshold or expand season.")
  }
  if (is.null(best)) {
    stop("No Sentinel-2 items for ", year, " in Jan–Mar. Expand dates in code or relax cloud.", call. = FALSE)
  }

  bands_10 <- s2$bands_10m
  bands_20 <- s2$bands_20m
  all_b <- c(bands_10, bands_20)

  layers <- vector("list", length(all_b))
  names(layers) <- all_b
  for (b in all_b) {
    href <- asset_href(best, b)
    if (is.na(href)) stop("Missing asset ", b, " on item ", best$id %||% "<no_id>", call. = FALSE)
    message("Streaming ", b, " …")
    layers[[b]] <- crop_mask_to_aoi(rast_from_signed_href(href, sign_ep), aoi)
  }

  ref <- layers[["B02"]]
  for (b in bands_20) {
    if (!terra::compareGeom(layers[[b]], ref, stopOnError = FALSE)) {
      layers[[b]] <- terra::resample(layers[[b]], ref, method = "bilinear")
    }
  }

  stk <- terra::rast(layers[all_b])
  names(stk) <- all_b
  terra::writeRaster(stk, out_path, overwrite = TRUE, gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES"))
  message("Wrote ", out_path)
  terra::rast(out_path)
}

# --------------------------------------------------------------------------- #
# Landsat C2 L2: 30m, PC band keys → SR_B* / QA_PIXEL names
# --------------------------------------------------------------------------- #

download_landsat_year <- function(
    year,
    platform,
    cfg,
    overwrite = FALSE) {
  raw <- file.path(cfg$paths$raw_dir, "landsat")
  dir.create(raw, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(raw, sprintf("LS_composite_%d.tif", year))
  if (file.exists(out_path) && !overwrite) {
    message("Exists, loading: ", out_path)
    return(terra::rast(out_path))
  }

  pc <- get_pc_cfg(cfg)
  ls <- pc$landsat
  bbox <- aoi_bbox_wgs84(cfg)
  aoi <- aoi_sf_polygon(cfg)
  sign_ep <- pc$sign_endpoint

  ms <- ls$season_month_start
  me <- ls$season_month_end
  dt <- season_datetime_range(year, ms, me)
  feats <- stac_search_features(pc$stac_url, ls$collection, bbox, dt, limit = ls$stac_limit)

  if (!length(feats)) {
    message(
      "No Landsat STAC items for ", year, " in months ", ms, "–", me,
      "; retrying full calendar year …"
    )
    dt <- season_datetime_range(year, 1L, 12L)
    feats <- stac_search_features(pc$stac_url, ls$collection, bbox, dt, limit = ls$stac_limit)
  }

  cloud_try <- unique(c(
    as.numeric(c(ls$cloud_cover_max_percent, ls$cloud_cover_relaxed_percent)),
    40, 60, 80, 100
  ))
  cloud_try <- cloud_try[is.finite(cloud_try)]
  cloud_try <- c(cloud_try, Inf)

  best <- NULL
  for (cm in cloud_try) {
    best <- pick_best_feature(feats, platform_want = platform, cloud_max = cm)
    if (!is.null(best)) break
  }

  if (is.null(best) && (as.integer(ms) != 1L || as.integer(me) != 12L)) {
    message("Landsat ", year, ": no scene in months ", ms, "–", me, "; searching full year …")
    dt_full <- season_datetime_range(year, 1L, 12L)
    feats_full <- stac_search_features(pc$stac_url, ls$collection, bbox, dt_full, limit = ls$stac_limit)
    if (length(feats_full)) {
      feats <- feats_full
      for (cm in cloud_try) {
        best <- pick_best_feature(feats, platform_want = platform, cloud_max = cm)
        if (!is.null(best)) break
      }
    }
  }

  if (is.null(best) && length(feats)) {
    plat_obs <- vapply(feats, function(f) normalize_platform(item_platform(f)), character(1))
    want_n <- normalize_platform(platform)
    any_plat <- any(!is.na(plat_obs) & plat_obs == want_n, na.rm = TRUE)
    if (!any_plat) {
      message(
        "No STAC item lists platform ", platform, " (", want_n, "); using lowest-cloud scene regardless of platform …"
      )
      for (cm in cloud_try) {
        best <- pick_best_feature(feats, platform_want = NULL, cloud_max = cm)
        if (!is.null(best)) break
      }
    }
  }
  if (is.null(best)) {
    stop(
      "No Landsat item for ", year, " (want platform ", platform, "). ",
      "Returned ", length(feats), " STAC item(s). ",
      "Try widening landsat.season_month_* in config/study_area.yml or increase stac_limit.",
      call. = FALSE
    )
  }

  bm <- ls$band_map
  keys <- names(bm)
  layers <- list()

  for (k in keys) {
    href <- asset_href(best, k)
    if (is.na(href)) stop("Missing Landsat asset key: ", k, call. = FALSE)
    out_nm <- as.character(bm[[k]])
    message("Streaming ", k, " → ", out_nm, " …")
    layers[[out_nm]] <- crop_mask_to_aoi(rast_from_signed_href(href, sign_ep), aoi)
  }

  nm_order <- vapply(keys, function(k) as.character(bm[[k]]), character(1))
  stk <- terra::rast(layers[nm_order])
  terra::writeRaster(stk, out_path, overwrite = TRUE, gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES"))
  message("Wrote ", out_path)
  terra::rast(out_path)
}

# --------------------------------------------------------------------------- #
# GlobeLand30 — merge manual tiles, crop, project, save
# --------------------------------------------------------------------------- #

prepare_globeland30 <- function(
    year,
    tile_paths,
    cfg,
    overwrite = FALSE) {
  if (!length(tile_paths)) stop("No GlobeLand30 tile paths provided.", call. = FALSE)
  raw <- file.path(cfg$paths$raw_dir, "globeland30")
  dir.create(raw, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(raw, sprintf("GL30_%d.tif", year))
  if (file.exists(out_path) && !overwrite) {
    message("Exists, loading: ", out_path)
    return(terra::rast(out_path))
  }

  rlist <- lapply(tile_paths, terra::rast)
  m <- if (length(rlist) == 1) rlist[[1]] else do.call(terra::merge, rlist)
  aoi <- aoi_sf_polygon(cfg)
  crs_t <- cfg$study_area$crs_projected
  aoi_t <- sf::st_transform(aoi, crs_t)
  v <- terra::vect(aoi_t)
  m2 <- terra::project(m, crs_t)
  x <- terra::crop(m2, v)
  x <- terra::mask(x, v)
  terra::writeRaster(x, out_path, overwrite = TRUE, gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES"))
  message("Wrote ", out_path)
  terra::rast(out_path)
}

# --------------------------------------------------------------------------- #
# Run full Phase 1
# --------------------------------------------------------------------------- #

run_phase1 <- function(cfg = load_config(), overwrite = FALSE) {
  with_pipeline_phase("Phase 1 — AOI datasets (Sentinel-2 + Landsat downloads)", {
    pc <- get_pc_cfg(cfg)

    for (yr in unlist(pc$sentinel2$years, use.names = FALSE)) {
      pipeline_log("Phase 1: Sentinel-2 year ", yr, " …")
      download_sentinel2_year(as.integer(yr), cfg, overwrite = overwrite)
    }

    acq <- pc$landsat$acquisitions
    for (i in seq_along(acq)) {
      row <- acq[[i]]
      pipeline_log("Phase 1: Landsat year ", row$year, " (", row$platform, ") …")
      download_landsat_year(row$year, row$platform, cfg, overwrite = overwrite)
    }

    message("GlobeLand30: place downloaded tiles on disk, then call prepare_globeland30(year, tile_paths, cfg).")
    invisible(TRUE)
  })
}

# Default execution when sourced interactively (comment out to import functions only)
if (identical(Sys.getenv("URBAN_SPRAWL_SKIP_PHASE1_RUN"), "")) {
  if (interactive()) {
    message("Phase 1 functions loaded. Run: run_phase1() or download_sentinel2_year(2018, load_config()).")
  }
}
