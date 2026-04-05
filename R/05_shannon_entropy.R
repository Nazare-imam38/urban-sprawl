#' Phase 5 — Shannon entropy, UEII, landscape metrics, hotspots, directional growth
#'
#' Inputs: Phase 3 `lulc_{year}.tif` in `data/processed/`. Centre + rings from `config$sprawl`.
#' Outputs: zone shapefiles, CSVs, `sprawl_entropy_map.tif`, `ueii_map.tif`, growth chart PNG.
#'
#' Working directory = project root.

source("R/00_setup.R")

# --------------------------------------------------------------------------- #
# Config & paths
# --------------------------------------------------------------------------- #

get_sprawl_cfg <- function(cfg) {
  s <- cfg$sprawl
  if (is.null(s)) stop("config: add `sprawl` block (centre_lon, centre_lat, rings).", call. = FALSE)
  s
}

path_lulc <- function(year, cfg) {
  # Relative to project root; callers can wrap with file.path(project_root(), ...)
  file.path(cfg$paths$processed_dir, sprintf("lulc_%d.tif", year))
}

phase_years_sorted <- function(cfg) {
  pc <- cfg$planetary_computer
  sort(unique(c(
    vapply(pc$landsat$acquisitions, function(z) as.integer(z$year), integer(1)),
    unlist(pc$sentinel2$years, use.names = FALSE)
  )))
}

study_aoi_polygon <- function(cfg) {
  sa <- cfg$study_area
  bb <- sf::st_bbox(
    c(xmin = sa$xmin, ymin = sa$ymin, xmax = sa$xmax, ymax = sa$ymax),
    crs = sa$crs_wgs84
  )
  sf::st_as_sf(sf::st_as_sfc(bb))
}

dissolve_by_id <- function(x, id_col = "zone_id") {
  ids <- unique(x[[id_col]])
  parts <- lapply(ids, function(id) {
    sub <- x[x[[id_col]] == id, ]
    if (nrow(sub) == 1L) return(sub)
    gcol <- attr(sub, "sf_column")
    one <- sub[1, ]
    one[[gcol]] <- sf::st_union(sf::st_geometry(sub))
    one
  })
  do.call(rbind, parts)
}

# --------------------------------------------------------------------------- #
# Concentric rings (0–5, 5–10, … km) in projected CRS
# --------------------------------------------------------------------------- #

city_centre_sf <- function(cfg) {
  s <- get_sprawl_cfg(cfg)
  p <- sf::st_point(c(s$centre_lon, s$centre_lat))
  sf::st_sf(geometry = sf::st_sfc(p, crs = cfg$study_area$crs_wgs84))
}

create_buffer_rings <- function(cfg) {
  s <- get_sprawl_cfg(cfg)
  crs_m <- cfg$study_area$crs_projected
  step_km <- s$ring_step_km %||% 5
  max_km <- s$ring_max_km %||% 30
  breaks_m <- seq(step_km, max_km, by = step_km) * 1000

  ctr <- sf::st_transform(city_centre_sf(cfg), crs_m)
  bufs <- lapply(breaks_m, function(d) sf::st_buffer(ctr, dist = d))
  bufs <- lapply(bufs, sf::st_make_valid)

  lst <- vector("list", length(breaks_m))
  g1 <- sf::st_geometry(bufs[[1]])
  lst[[1]] <- sf::st_sf(
    zone_id = 1L,
    inner_km = 0,
    outer_km = breaks_m[1] / 1000,
    geometry = g1,
    crs = sf::st_crs(crs_m)
  )
  for (i in seq_len(length(breaks_m))[-1]) {
    g <- sf::st_difference(sf::st_geometry(bufs[[i]]), sf::st_geometry(bufs[[i - 1]]))
    g <- sf::st_make_valid(g)
    lst[[i]] <- sf::st_sf(
      zone_id = i,
      inner_km = breaks_m[i - 1] / 1000,
      outer_km = breaks_m[i] / 1000,
      geometry = g,
      crs = sf::st_crs(crs_m)
    )
  }
  rings <- do.call(rbind, lst)

  aoi <- sf::st_transform(study_aoi_polygon(cfg), crs_m)
  rings <- sf::st_intersection(sf::st_make_valid(rings), sf::st_make_valid(aoi))
  rings <- sf::st_make_valid(rings)
  if (nrow(rings) == 0L) stop("No ring geometry after AOI clip — check centre and bbox.", call. = FALSE)
  dissolve_by_id(rings, "zone_id")
}

# --------------------------------------------------------------------------- #
# Urban counts per zone × year
# --------------------------------------------------------------------------- #

count_urban_per_zone <- function(lulc_rast, zones_sf) {
  results <- data.frame(
    zone_id = zones_sf$zone_id,
    urban_px = NA_integer_,
    total_px = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(zones_sf))) {
    zv <- terra::vect(zones_sf[i, ])
    # Ensure CRS matches raster before crop/mask (prevents "[crop] extents do not overlap" from CRS mismatch)
    if (!is.na(terra::crs(zv)) && !is.na(terra::crs(lulc_rast)) && terra::crs(zv) != terra::crs(lulc_rast)) {
      zv <- terra::project(zv, lulc_rast)
    }
    zc <- tryCatch(terra::crop(lulc_rast, zv), error = function(e) NULL)
    if (is.null(zc)) {
      results$urban_px[i] <- 0L
      results$total_px[i] <- 0L
      next
    }
    zm <- tryCatch(terra::mask(zc, zv), error = function(e) NULL)
    if (is.null(zm)) {
      results$urban_px[i] <- 0L
      results$total_px[i] <- 0L
      next
    }
    v <- terra::values(zm[[1]], mat = FALSE)
    results$urban_px[i] <- sum(v == 1L, na.rm = TRUE)
    results$total_px[i] <- sum(!is.na(v))
  }
  results
}

urban_area_km2_per_zone <- function(urban_px, lulc_rast) {
  px_km2 <- prod(terra::res(lulc_rast)) / 1e6
  urban_px * px_km2
}

zone_area_km2 <- function(zones_sf) {
  as.numeric(sf::st_area(zones_sf)) / 1e6
}

# --------------------------------------------------------------------------- #
# Shannon entropy (zonal urban proportions)
# --------------------------------------------------------------------------- #

shannon_entropy_from_counts <- function(urban_counts, n_zones) {
  total <- sum(urban_counts, na.rm = TRUE)
  if (total <= 0) {
    return(list(H = 0, H_rel = 0, Pi = rep(0, length(urban_counts))))
  }
  Pi <- urban_counts / total
  Pi_pos <- Pi[Pi > 0]
  H <- -sum(Pi_pos * log(Pi_pos))
  H_rel <- if (n_zones > 1L) H / log(n_zones) else 0
  list(H = H, H_rel = H_rel, Pi = Pi)
}

# --------------------------------------------------------------------------- #
# UEII: (UA_t2 - UA_t1) / (TLA * delta_years) * 100
# --------------------------------------------------------------------------- #

compute_ueii <- function(ua_t1, ua_t2, zone_area_km2, years_between) {
  if (years_between <= 0) return(rep(NA_real_, length(ua_t1)))
  (ua_t2 - ua_t1) / (zone_area_km2 * years_between) * 100
}

# --------------------------------------------------------------------------- #
# landscapemetrics (class 1 = Urban)
# --------------------------------------------------------------------------- #

lsm_urban_value <- function(lulc_rast, metric_fn) {
  x <- tryCatch(metric_fn(lulc_rast), error = function(e) NULL)
  if (is.null(x)) return(NA_real_)
  d <- as.data.frame(x)
  if (!nrow(d) || !"class" %in% names(d)) return(NA_real_)
  vv <- d$value[d$class == 1L]
  if (!length(vv)) NA_real_ else as.numeric(vv[[1]])
}

compute_landscape_metrics_year <- function(lulc_rast) {
  lulc_rast <- terra::as.int(terra::round(lulc_rast))
  c(
    PD = lsm_urban_value(lulc_rast, landscapemetrics::lsm_c_pd),
    LPI = lsm_urban_value(lulc_rast, landscapemetrics::lsm_c_lpi),
    ED = lsm_urban_value(lulc_rast, landscapemetrics::lsm_c_ed),
    AI = lsm_urban_value(lulc_rast, landscapemetrics::lsm_c_ai),
    MPA = lsm_urban_value(lulc_rast, landscapemetrics::lsm_c_area_mn)
  )
}

# --------------------------------------------------------------------------- #
# Eight directional sectors (bearings from North, clockwise)
# --------------------------------------------------------------------------- #

create_directional_sectors <- function(cfg) {
  s <- get_sprawl_cfg(cfg)
  crs_m <- cfg$study_area$crs_projected
  R <- s$sector_radius_m %||% 80000
  ctr <- sf::st_transform(city_centre_sf(cfg), crs_m)
  xy <- sf::st_coordinates(ctr)[1, , drop = FALSE]
  cx <- xy[1, 1]
  cy <- xy[1, 2]

  bearings <- seq(0, 315, by = 45)
  nm <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  polys <- vector("list", 8L)
  for (i in seq_len(8L)) {
    th1 <- bearings[i] * pi / 180
    th2 <- if (i < 8L) bearings[i + 1] * pi / 180 else 2 * pi
    arc <- seq(th1, th2, length.out = 40)
    ax <- cx + R * sin(arc)
    ay <- cy + R * cos(arc)
    coords <- rbind(c(cx, cy), cbind(ax, ay), c(cx, cy))
    polys[[i]] <- sf::st_polygon(list(coords))
  }
  sfc <- sf::st_sfc(polys, crs = crs_m)
  sec <- sf::st_sf(
    sector_id = 1L:8L,
    sector_name = nm,
    geometry = sfc,
    crs = sf::st_crs(crs_m)
  )
  aoi <- sf::st_transform(study_aoi_polygon(cfg), crs_m)
  sec2 <- sf::st_intersection(sf::st_make_valid(sec), sf::st_make_valid(aoi))
  sec2 <- sf::st_make_valid(sec2)
  dissolve_by_id(sec2, "sector_id")
}

count_urban_per_sector <- function(lulc_rast, sectors_sf) {
  out <- data.frame(
    sector_id = sectors_sf$sector_id,
    sector_name = sectors_sf$sector_name,
    urban_px = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(sectors_sf))) {
    zv <- terra::vect(sectors_sf[i, ])
    if (!is.na(terra::crs(zv)) && !is.na(terra::crs(lulc_rast)) && terra::crs(zv) != terra::crs(lulc_rast)) {
      zv <- terra::project(zv, lulc_rast)
    }
    zc <- tryCatch(terra::crop(lulc_rast, zv), error = function(e) NULL)
    if (is.null(zc)) {
      out$urban_px[i] <- 0L
      next
    }
    zm <- tryCatch(terra::mask(zc, zv), error = function(e) NULL)
    if (is.null(zm)) {
      out$urban_px[i] <- 0L
      next
    }
    v <- terra::values(zm[[1]], mat = FALSE)
    out$urban_px[i] <- sum(v == 1L, na.rm = TRUE)
  }
  px_km2 <- prod(terra::res(lulc_rast)) / 1e6
  out$urban_km2 <- out$urban_px * px_km2
  out
}

# --------------------------------------------------------------------------- #
# Rasterize zonal attribute (Pi, UEII) onto LULC grid
# --------------------------------------------------------------------------- #

rasterize_zone_field <- function(zones_sf, template_rast, field_name) {
  v <- terra::vect(zones_sf)
  if (!is.na(terra::crs(v)) && !is.na(terra::crs(template_rast)) && terra::crs(v) != terra::crs(template_rast)) {
    v <- terra::project(v, template_rast)
  }
  tryCatch(
    terra::rasterize(v, template_rast, field = field_name),
    error = function(e) template_rast
  )
}

# --------------------------------------------------------------------------- #
# Hotspots: UEII high + fringe distance
# --------------------------------------------------------------------------- #

identify_sprawl_hotspots <- function(zones_sf, ueii_df, fringe_km) {
  # ueii_df: zone_id, ueii (one period or max)
  u <- merge(
    sf::st_drop_geometry(zones_sf),
    ueii_df,
    by = "zone_id",
    all.x = TRUE
  )
  m <- mean(u$ueii, na.rm = TRUE)
  s <- stats::sd(u$ueii, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- 0
  hot <- !is.na(u$ueii) & u$ueii > (m + s) & u$inner_km > fringe_km
  idx <- which(hot)
  if (!length(idx)) {
    return(zones_sf[FALSE, ])
  }
  zones_sf[idx, ]
}

# --------------------------------------------------------------------------- #
# Main pipeline
# --------------------------------------------------------------------------- #

`%||%` <- function(x, y) if (is.null(x)) y else x

run_phase5 <- function(cfg = load_config()) {
  with_pipeline_phase("Phase 5 — Shannon entropy, UEII, sprawl metrics", {
  proc <- file.path(project_root(), cfg$paths$processed_dir)
  outd <- file.path(project_root(), cfg$paths$outputs_dir)
  dir.create(proc, showWarnings = FALSE, recursive = TRUE)
  dir.create(outd, showWarnings = FALSE, recursive = TRUE)

  s <- get_sprawl_cfg(cfg)
  fringe_km <- s$hotspot_fringe_inner_km %||% 10

  years <- phase_years_sorted(cfg)
  for (y in years) {
    lp <- file.path(project_root(), path_lulc(y, cfg))
    if (!file.exists(lp)) {
      stop("Missing ", lp, " — run Phase 3 first.", call. = FALSE)
    }
  }

  rings <- create_buffer_rings(cfg)
  sf::st_write(
    rings,
    file.path(proc, "buffer_rings.shp"),
    delete_dsn = TRUE,
    quiet = TRUE
  )

  sectors <- create_directional_sectors(cfg)
  sf::st_write(
    sectors,
    file.path(proc, "directional_sectors.shp"),
    delete_dsn = TRUE,
    quiet = TRUE
  )

  # --- Urban counts wide table ---
  wide <- sf::st_drop_geometry(rings)[, c("zone_id", "inner_km", "outer_km"), drop = FALSE]
  tmpl <- terra::rast(file.path(project_root(), path_lulc(years[length(years)], cfg)))
  # Ensure zones are in same CRS as rasters before terra crop/mask
  crs_r <- sf::st_crs(terra::crs(tmpl))
  rings <- sf::st_transform(rings, crs_r)
  sectors <- sf::st_transform(sectors, crs_r)
  wide$zone_area_km2 <- zone_area_km2(rings)
  n_zones <- nrow(rings)

  for (y in years) {
    lu <- terra::rast(file.path(project_root(), path_lulc(y, cfg)))
    lu <- terra::round(lu)
    cnt <- count_urban_per_zone(lu, rings)
    nm_px <- paste0("urban_px_", y)
    nm_km <- paste0("urban_km2_", y)
    wide[[nm_px]] <- cnt$urban_px[match(wide$zone_id, cnt$zone_id)]
    wide[[nm_km]] <- urban_area_km2_per_zone(wide[[nm_px]], lu)
  }
  utils::write.csv(wide, file.path(proc, "urban_counts_per_zone.csv"), row.names = FALSE)

  # --- Shannon entropy per year ---
  ent_rows <- list()
  pi_layers <- list()
  for (y in years) {
    col_px <- paste0("urban_px_", y)
    uc <- wide[[col_px]]
    se <- shannon_entropy_from_counts(uc, n_zones = n_zones)
    ent_rows[[length(ent_rows) + 1L]] <- data.frame(
      year = y,
      H = se$H,
      H_rel = se$H_rel,
      n_zones = n_zones,
      stringsAsFactors = FALSE
    )
    rings_pi <- rings
    rings_pi$Pi <- se$Pi[match(rings_pi$zone_id, wide$zone_id)]
    pi_layers[[as.character(y)]] <- rasterize_zone_field(rings_pi, tmpl, "Pi")
  }
  ent_df <- do.call(rbind, ent_rows)
  utils::write.csv(ent_df, file.path(proc, "shannon_entropy_results.csv"), row.names = FALSE)

  stk_pi <- terra::rast(pi_layers)
  names(stk_pi) <- paste0("Pi_", years)
  terra::writeRaster(
    stk_pi,
    file.path(outd, "sprawl_entropy_map.tif"),
    overwrite = TRUE,
    gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
  )

  # --- UEII per zone per transition ---
  trans <- cfg$cnn$transitions
  ueii_list <- list()
  ueii_maps <- list()
  for (t in trans) {
    y1 <- as.integer(t$pre)
    y2 <- as.integer(t$post)
    c1 <- paste0("urban_km2_", y1)
    c2 <- paste0("urban_km2_", y2)
    if (!all(c(c1, c2) %in% names(wide))) next
    dy <- y2 - y1
    ua1 <- wide[[c1]]
    ua2 <- wide[[c2]]
    ue <- compute_ueii(ua1, ua2, wide$zone_area_km2, dy)
    ueii_list[[length(ueii_list) + 1L]] <- data.frame(
      zone_id = wide$zone_id,
      pre_year = y1,
      post_year = y2,
      delta_years = dy,
      ueii = ue,
      stringsAsFactors = FALSE
    )
    rz <- rings
    rz$ueii <- ue[match(rz$zone_id, wide$zone_id)]
    nm <- paste0("UEII_", y1, "_", y2)
    ueii_maps[[nm]] <- rasterize_zone_field(rz, tmpl, "ueii")
  }
  if (length(ueii_list)) {
    ueii_df <- do.call(rbind, ueii_list)
    utils::write.csv(ueii_df, file.path(proc, "ueii_results.csv"), row.names = FALSE)
    stk_ue <- terra::rast(ueii_maps)
    terra::writeRaster(
      stk_ue,
      file.path(outd, "ueii_map.tif"),
      overwrite = TRUE,
      gdal = c("BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "TILED=YES")
    )
  }

  # Hotspots: use latest transition in table
  if (length(ueii_list)) {
    last_ue <- ueii_list[[length(ueii_list)]]
    hot <- identify_sprawl_hotspots(rings, last_ue[, c("zone_id", "ueii")], fringe_km = fringe_km)
    if (nrow(hot)) {
      shp_path <- file.path(proc, "sprawl_hotspots.shp")
      ok <- tryCatch(
        {
          sf::st_write(hot, shp_path, delete_dsn = TRUE, quiet = TRUE)
          TRUE
        },
        error = function(e) FALSE
      )
      if (!ok) {
        gpkg_path <- file.path(proc, "sprawl_hotspots.gpkg")
        sf::st_write(hot, gpkg_path, layer = "sprawl_hotspots", delete_layer = TRUE, quiet = TRUE)
        message("Wrote hotspots to ", gpkg_path, " (shapefile write failed).")
      }
    } else {
      message("No sprawl hotspot zones met criteria (UEII / fringe).")
    }
  }

  # --- Landscape metrics ---
  lm_rows <- list()
  for (y in years) {
    lu <- terra::rast(path_lulc(y, cfg))
    lu <- terra::round(lu)
    m <- compute_landscape_metrics_year(lu)
    lm_rows[[length(lm_rows) + 1L]] <- data.frame(
      year = y,
      PD = m[["PD"]],
      LPI = m[["LPI"]],
      ED = m[["ED"]],
      AI = m[["AI"]],
      MPA = m[["MPA"]],
      stringsAsFactors = FALSE
    )
  }
  utils::write.csv(do.call(rbind, lm_rows), file.path(proc, "landscape_metrics.csv"), row.names = FALSE)

  # --- Directional growth chart ---
  y0 <- min(years)
  y1 <- max(years)
  lu0 <- terra::rast(path_lulc(y0, cfg)) |> terra::round()
  lu1 <- terra::rast(path_lulc(y1, cfg)) |> terra::round()
  s0 <- count_urban_per_sector(lu0, sectors)
  s1 <- count_urban_per_sector(lu1, sectors)
  dg <- data.frame(
    sector = s0$sector_name,
    urban_km2_t0 = s0$urban_km2,
    urban_km2_t1 = s1$urban_km2,
    growth_pct = ifelse(
      s0$urban_km2 > 0,
      (s1$urban_km2 - s0$urban_km2) / s0$urban_km2 * 100,
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
  dg$growth_pct[is.infinite(dg$growth_pct)] <- NA_real_

  sector_levels <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  p <- ggplot2::ggplot(dg, ggplot2::aes(x = factor(sector, levels = sector_levels), y = growth_pct)) +
    ggplot2::geom_col(fill = "steelblue4") +
    ggplot2::labs(
      title = sprintf("Urban area growth by direction (%d → %d)", y0, y1),
      x = "Sector (from centre)",
      y = "Growth (%)",
      caption = "Based on LULC urban class (1) counts per wedge"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  ggplot2::ggsave(
    file.path(outd, "directional_growth_chart.png"),
    p,
    width = 8,
    height = 5,
    dpi = 150
  )

  pipeline_log("Phase 5 complete. See ", proc, " and ", outd)
  invisible(list(
    urban_wide = wide,
    entropy = ent_df,
    ueii = if (length(ueii_list)) do.call(rbind, ueii_list) else NULL,
    directional = dg
  ))
  })
}

if (interactive()) {
  message("Phase 5: run_phase5() after Phase 3 LULC rasters exist.")
}
