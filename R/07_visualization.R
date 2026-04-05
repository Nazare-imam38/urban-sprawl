#' Phase 7 — Publication-style maps, charts, WebGIS, optional R Markdown report
#'
#' Reads outputs from Phases 3–6. Writes under `outputs/maps/`, `outputs/charts/`,
#' `outputs/urban_sprawl_dashboard.html`. Optional: **ggcorrplot**, **leaflet.extras**,
#' **htmlwidgets**, **rmarkdown** (install if missing — see `run_phase7()` message).
#'
#' Working directory = project root.

source("R/00_setup.R")

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------------------------------------------------------- #
# Style
# --------------------------------------------------------------------------- #

# Codes 0–4 (numeric). Leaflet uses numeric domain; tmap v4 needs palette names = factor class labels.
lulc_class_labels <- c("Other/NoData", "Urban", "Vegetation", "Bare Soil", "Water")

lulc_hex <- c(
  "#9E9E9E", # 0 = Other / NoData
  "#E8392A", # 1 = Urban
  "#2ECC71", # 2 = Vegetation
  "#F5CBA7", # 3 = Bare soil
  "#3498DB"  # 4 = Water
)
names(lulc_hex) <- as.character(0:4)

# tmap categorical rasters use class names from levels(); palette names must match those labels.
lulc_hex_tmap <- stats::setNames(as.character(lulc_hex), lulc_class_labels)

# Leaflet: fixed order for domain 0:4 (do not rely on as.vector() order of a named vector).
lulc_palette_leaflet <- as.character(lulc_hex[as.character(0:4)])

# tmap v4: legend palette must match **active** levels only (downsampling can drop classes).
lulc_palette_active <- function(r) {
  lv <- terra::levels(r)[[1L]]
  if (is.null(lv) || !nrow(lv) || !("class" %in% names(lv))) {
    return(lulc_hex_tmap)
  }
  cls <- as.character(lv$class)
  pal <- lulc_hex_tmap[cls]
  pal <- pal[!is.na(pal)]
  if (!length(pal)) lulc_hex_tmap else pal
}

prep_lulc_for_tmap <- function(lu) {
  # Some rasters can contain 0/255/other codes (nodata or artefacts). Force to 0–4 and set factor levels
  # so tmap doesn't mis-handle labels/palette lengths and crash.
  r <- terra::round(lu)
  r <- terra::ifel(is.na(r), NA_integer_, r)
  r <- terra::ifel(r %in% 0:4, r, 0L)
  r <- terra::as.factor(r)
  lev <- data.frame(
    ID = 0:4,
    class = lulc_class_labels,
    stringsAsFactors = FALSE
  )
  levels(r) <- lev
  tryCatch(terra::droplevels(r), error = function(e) r)
}

change_hex_tmap <- c("New urban" = "#E74C3C")

prep_new_urban_overlay <- function(ch) {
  r <- terra::ifel(ch == 1L, 1L, NA_integer_)
  vals <- terra::values(r, mat = FALSE)
  if (!any(vals == 1L, na.rm = TRUE)) {
    return(NULL)
  }
  r <- terra::as.factor(r)
  levels(r) <- data.frame(ID = 1L, class = "New urban", stringsAsFactors = FALSE)
  names(r) <- "new_urban"
  r
}

theme_sprawl <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray40", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = 8, color = "gray60")
    )
}

# --------------------------------------------------------------------------- #
# Config helpers
# --------------------------------------------------------------------------- #

get_phase7_cfg <- function(cfg) {
  p <- cfg$phase7 %||% list()
  p$map_dpi <- as.numeric(p$map_dpi %||% 300)
  p$chart_dpi <- as.numeric(p$chart_dpi %||% 300)
  p$leaflet_aggregate_fact <- as.integer(p$leaflet_aggregate_fact %||% 4L)
  p
}

# All Phase 7 I/O must be anchored to project_root(), not relative getwd() (RStudio / Knit often change wd).
r_proj_root <- function() {
  normalizePath(project_root(), winslash = "/", mustWork = FALSE)
}

resolved_out_dir <- function(cfg) {
  od <- cfg$paths$outputs_dir %||% "outputs"
  if (nzchar(od) && (startsWith(od, "/") || grepl("^[A-Za-z]:[\\\\/]", od))) {
    normalizePath(od, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(r_proj_root(), od), winslash = "/", mustWork = FALSE)
  }
}

resolved_proc_dir <- function(cfg) {
  pd <- cfg$paths$processed_dir %||% "data/processed"
  if (nzchar(pd) && (startsWith(pd, "/") || grepl("^[A-Za-z]:[\\\\/]", pd))) {
    normalizePath(pd, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(r_proj_root(), pd), winslash = "/", mustWork = FALSE)
  }
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

path_lulc <- function(year, cfg) {
  file.path(resolved_proc_dir(cfg), sprintf("lulc_%d.tif", year))
}

path_change_bin <- function(y1, y2, cfg) {
  file.path(resolved_out_dir(cfg), sprintf("change_map_%d_%d.tif", y1, y2))
}

ensure_dirs <- function(cfg) {
  out <- resolved_out_dir(cfg)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out, "maps"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out, "charts"), showWarnings = FALSE, recursive = TRUE)
  invisible(NULL)
}

phase7_try_ggsave <- function(filename, plot, width, height, dpi) {
  tryCatch(
    {
      ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi)
      if (!file.exists(filename)) {
        warning("ggsave finished but file is missing: ", filename, call. = FALSE)
      } else {
        message("Wrote ", filename)
      }
      invisible(filename)
    },
    error = function(e) {
      warning(
        "Chart save failed (", basename(filename), "): ", conditionMessage(e),
        call. = FALSE
      )
      invisible(NULL)
    }
  )
}

# --------------------------------------------------------------------------- #
# tmap — LULC panel
# --------------------------------------------------------------------------- #

save_lulc_panel_map <- function(cfg, dpi = 300) {
  years <- phase_years_sorted(cfg)
  miss <- years[!vapply(years, function(y) file.exists(path_lulc(y, cfg)), logical(1))]
  if (length(miss)) {
    warning("LULC panel: missing rasters for years: ", paste(miss, collapse = ", "), call. = FALSE)
    years <- setdiff(years, miss)
  }
  if (length(years) < 2L) {
    warning("LULC panel: not enough years; skip.", call. = FALSE)
    return(invisible(NULL))
  }

  tmap::tmap_mode("plot")
  maps <- vector("list", length(years))
  for (i in seq_along(years)) {
    r <- prep_lulc_for_tmap(terra::rast(path_lulc(years[i], cfg)))
    maps[[i]] <- tmap::tm_shape(r) +
      tmap::tm_raster(
        style = "cat",
        palette = lulc_palette_active(r),
        title = ""
      ) +
      tmap::tm_layout(
        title = as.character(years[i]),
        title.position = c("center", "top"),
        title.size = 1.2,
        legend.show = (i == length(years)),
        frame = FALSE
      ) +
      tmap::tm_scale_bar(position = c("left", "bottom"), text.size = 0.6) +
      tmap::tm_compass(position = c("right", "top"), size = 1.2)
  }
  ncol <- min(3L, length(years))
  nrow <- ceiling(length(years) / ncol)
  panel <- do.call(tmap::tmap_arrange, c(maps, list(ncol = ncol, nrow = nrow)))
  out <- file.path(resolved_out_dir(cfg), "maps", "lulc_panel_map.png")
  tmap::tmap_save(panel, filename = out, width = 14, height = max(6, 3.5 * nrow), dpi = dpi)
  message("Wrote ", out)
  invisible(out)
}

# --------------------------------------------------------------------------- #
# tmap — change maps (pre LULC + new urban overlay)
# --------------------------------------------------------------------------- #

save_change_maps <- function(cfg, dpi = 300) {
  trans <- cfg$cnn$transitions
  if (is.null(trans) || !length(trans)) {
    warning("No cnn.transitions in config; skip change maps.", call. = FALSE)
    return(invisible(NULL))
  }
  tmap::tmap_mode("plot")
  for (t in trans) {
    y0 <- as.integer(t$pre)
    y1 <- as.integer(t$post)
    fp <- path_change_bin(y0, y1, cfg)
    lp <- path_lulc(y0, cfg)
    if (!file.exists(fp) || !file.exists(lp)) {
      warning("Skip change map ", y0, "→", y1, " — missing file.", call. = FALSE)
      next
    }
    lulc_pre <- prep_lulc_for_tmap(terra::rast(lp))
    ch <- terra::rast(fp)
    ch <- terra::resample(ch, lulc_pre, method = "near")
    overlay <- prep_new_urban_overlay(ch)

    m <- tmap::tm_shape(lulc_pre) +
      tmap::tm_raster(
        style = "cat",
        palette = lulc_palette_active(lulc_pre),
        title = "Land cover (pre)"
      )
    if (!is.null(overlay)) {
      m <- m +
        tmap::tm_shape(overlay) +
        tmap::tm_raster(
          style = "cat",
          palette = change_hex_tmap,
          title = "Change"
        )
    }
    m <- m +
      tmap::tm_layout(
        title = paste0("Urban expansion: ", y0, " → ", y1),
        title.position = c("center", "top"),
        legend.outside = TRUE,
        frame = FALSE
      ) +
      tmap::tm_scale_bar(position = c("left", "bottom"), text.size = 0.6) +
      tmap::tm_compass(position = c("right", "top"), size = 1.2)

    out <- file.path(
      resolved_out_dir(cfg), "maps",
      sprintf("change_map_%d_%d.png", y0, y1)
    )
    tmap::tmap_save(m, filename = out, width = 9, height = 8, dpi = dpi)
    message("Wrote ", out)
  }
  invisible(NULL)
}

# --------------------------------------------------------------------------- #
# Hotspot context map
# --------------------------------------------------------------------------- #

save_hotspot_map <- function(cfg, dpi = 300) {
  y_last <- max(phase_years_sorted(cfg))
  lu_path <- path_lulc(y_last, cfg)
  hot_path <- file.path(resolved_proc_dir(cfg), "sprawl_hotspots.shp")
  if (!file.exists(lu_path)) {
    warning("Hotspot map: missing ", lu_path, call. = FALSE)
    return(invisible(NULL))
  }
  lu <- prep_lulc_for_tmap(terra::rast(lu_path))
  tmap::tmap_mode("plot")
  m <- tmap::tm_shape(lu) +
    tmap::tm_raster(
      style = "cat",
      palette = lulc_palette_active(lu),
      title = ""
    ) +
    tmap::tm_layout(
      title = paste0("LULC ", y_last, " & sprawl hotspots"),
      title.position = c("center", "top"),
      frame = FALSE,
      legend.outside = TRUE
    ) +
    tmap::tm_scale_bar(position = c("left", "bottom"), text.size = 0.6) +
    tmap::tm_compass(position = c("right", "top"), size = 1.2)

  if (file.exists(hot_path)) {
    hot <- sf::st_read(hot_path, quiet = TRUE)
    hot <- sf::st_transform(hot, terra::crs(lu))
    m <- m +
      tmap::tm_shape(hot) +
      tmap::tm_fill(col = "#E74C3C", alpha = 0.28) +
      tmap::tm_borders(col = "#C0392B", lwd = 2)
  }

  out <- file.path(resolved_out_dir(cfg), "maps", "sprawl_hotspot_map.png")
  tmap::tmap_save(m, filename = out, width = 10, height = 9, dpi = dpi)
  message("Wrote ", out)
  invisible(out)
}

# --------------------------------------------------------------------------- #
# ggplot2 charts
# --------------------------------------------------------------------------- #

save_urban_growth_chart <- function(cfg, dpi = 300) {
  path <- file.path(resolved_proc_dir(cfg), "urban_area_stats.csv")
  if (!file.exists(path)) {
    warning("Missing urban_area_stats.csv — run Phase 3.", call. = FALSE)
    return(invisible(NULL))
  }
  u <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("Year", "Urban_km2") %in% names(u))) {
    warning("urban_area_stats.csv missing Year / Urban_km2.", call. = FALSE)
    return(invisible(NULL))
  }
  u <- u[order(u$Year), , drop = FALSE]

  p <- ggplot2::ggplot(u, ggplot2::aes(x = Year, y = Urban_km2)) +
    ggplot2::geom_area(fill = "#E74C3C", alpha = 0.15) +
    ggplot2::geom_line(color = "#E74C3C", linewidth = 1.2) +
    ggplot2::geom_point(color = "#C0392B", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(round(Urban_km2, 1), " km²")),
      vjust = -1.1,
      size = 3.2,
      color = "#C0392B"
    ) +
    ggplot2::labs(
      title = "Urban area expansion",
      subtitle = paste("Study area:", cfg$study_area$name),
      x = "Year",
      y = "Urban area (km²)",
      caption = "Source: classified LULC (Landsat & Sentinel-2)"
    ) +
    ggplot2::scale_x_continuous(breaks = unique(u$Year)) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    theme_sprawl()

  out <- file.path(resolved_out_dir(cfg), "charts", "urban_growth_trajectory.png")
  phase7_try_ggsave(out, p, width = 10, height = 6, dpi = dpi)
  invisible(out)
}

save_entropy_chart <- function(cfg, dpi = 300) {
  path <- file.path(resolved_proc_dir(cfg), "shannon_entropy_results.csv")
  if (!file.exists(path)) {
    warning("Missing shannon_entropy_results.csv — run Phase 5.", call. = FALSE)
    return(invisible(NULL))
  }
  e <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"H_rel" %in% names(e)) {
    warning("shannon_entropy_results.csv missing H_rel.", call. = FALSE)
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(e, ggplot2::aes(x = year, y = H_rel)) +
    ggplot2::geom_line(color = "#8E44AD", linewidth = 1.2) +
    ggplot2::geom_point(color = "#6C3483", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(H_rel, 3)),
      vjust = -1.1,
      size = 3.2
    ) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed", color = "#E74C3C", linewidth = 0.6) +
    ggplot2::annotate(
      "text",
      x = min(e$year, na.rm = TRUE),
      y = 0.52,
      label = "Moderate sprawl (H_rel = 0.5)",
      size = 3,
      color = "gray40",
      hjust = 0
    ) +
    ggplot2::annotate(
      "text",
      x = min(e$year, na.rm = TRUE),
      y = 0.72,
      label = "High sprawl (H_rel = 0.7)",
      size = 3,
      color = "#E74C3C",
      hjust = 0
    ) +
    ggplot2::labs(
      title = "Urban sprawl intensity (Shannon entropy)",
      subtitle = "Higher H_rel = more dispersed urban across buffer zones",
      x = "Year",
      y = "Relative Shannon entropy (H / ln(n))",
      caption = "n = number of concentric zones (Phase 5)"
    ) +
    ggplot2::scale_x_continuous(breaks = unique(e$year)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    theme_sprawl()

  out <- file.path(resolved_out_dir(cfg), "charts", "shannon_entropy_chart.png")
  phase7_try_ggsave(out, p, width = 10, height = 6, dpi = dpi)
  invisible(out)
}

save_driver_chart <- function(cfg, dpi = 300) {
  path <- file.path(resolved_out_dir(cfg), "logistic_regression_results.csv")
  if (!file.exists(path)) {
    warning("Missing logistic_regression_results.csv — run Phase 6.", call. = FALSE)
    return(invisible(NULL))
  }
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("term", "estimate", "std_error") %in% names(d))) {
    warning("logistic_regression_results.csv missing term/estimate/std_error.", call. = FALSE)
    return(invisible(NULL))
  }
  d <- d[d$term != "(Intercept)", , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))

  d$coefficient <- d$estimate
  d$ci_lower <- d$estimate - 1.96 * d$std_error
  d$ci_upper <- d$estimate + 1.96 * d$std_error
  d$direction <- ifelse(d$coefficient > 0, "Increases P(growth)", "Decreases P(growth)")
  d$predictor <- d$term

  d$predictor <- stats::reorder(d$predictor, abs(d$coefficient))
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = predictor, y = coefficient, fill = direction)
  ) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2,
      color = "gray30"
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.6, color = "gray25") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "Increases P(growth)" = "#E74C3C",
        "Decreases P(growth)" = "#3498DB"
      )
    ) +
    ggplot2::labs(
      title = "Spatial drivers of urban growth (logistic regression)",
      subtitle = "Coefficients with approximate 95% Wald intervals",
      x = NULL,
      y = "Coefficient (log-odds scale)",
      fill = NULL,
      caption = "Phase 6 sample-based model; positive = higher growth odds"
    ) +
    theme_sprawl() +
    ggplot2::theme(legend.position = "bottom")

  out <- file.path(resolved_out_dir(cfg), "charts", "driver_importance_chart.png")
  phase7_try_ggsave(out, p, width = 10, height = 7, dpi = dpi)
  invisible(out)
}

save_landscape_metrics_chart <- function(cfg, dpi = 300) {
  path <- file.path(resolved_proc_dir(cfg), "landscape_metrics.csv")
  if (!file.exists(path)) {
    warning("Missing landscape_metrics.csv — run Phase 5.", call. = FALSE)
    return(invisible(NULL))
  }
  lm <- utils::read.csv(path, stringsAsFactors = FALSE)
  long <- rbind(
    data.frame(year = lm$year, metric = "PD", value = lm$PD),
    data.frame(year = lm$year, metric = "LPI", value = lm$LPI),
    data.frame(year = lm$year, metric = "ED", value = lm$ED)
  )
  long <- long[stats::complete.cases(long), , drop = FALSE]

  p <- ggplot2::ggplot(long, ggplot2::aes(x = year, y = value, color = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = "Landscape structure metrics (urban class)",
      subtitle = "Patch density (PD), largest patch index (LPI), edge density (ED)",
      x = "Year",
      y = "Value",
      caption = "landscapemetrics on LULC class 1 (Urban)"
    ) +
    theme_sprawl() +
    ggplot2::theme(legend.position = "none")

  out <- file.path(resolved_out_dir(cfg), "charts", "landscape_metrics_chart.png")
  phase7_try_ggsave(out, p, width = 11, height = 4.5, dpi = dpi)
  invisible(out)
}

save_correlation_heatmap <- function(cfg, dpi = 300) {
  path <- file.path(resolved_out_dir(cfg), "correlation_matrix.csv")
  if (!file.exists(path)) {
    warning("Missing correlation_matrix.csv — run Phase 6.", call. = FALSE)
    return(invisible(NULL))
  }
  if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
    message("Install ggcorrplot for correlation heatmap: install.packages(\"ggcorrplot\")")
    return(invisible(NULL))
  }

  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"variable" %in% names(raw)) {
    warning("correlation_matrix.csv missing `variable` column.", call. = FALSE)
    return(invisible(NULL))
  }
  rn <- raw$variable
  mat <- as.matrix(raw[, setdiff(names(raw), "variable"), drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- rn

  p <- ggcorrplot::ggcorrplot(
    mat,
    method = "circle",
    type = "lower",
    lab = TRUE,
    lab_size = 3,
    colors = c("#3498DB", "white", "#E74C3C"),
    title = "Correlation matrix — sprawl drivers",
    ggtheme = theme_sprawl()
  )

  out <- file.path(resolved_out_dir(cfg), "charts", "correlation_heatmap.png")
  phase7_try_ggsave(out, p, width = 9, height = 8, dpi = dpi)
  invisible(out)
}

# --------------------------------------------------------------------------- #
# OSM roads (leaflet)
# --------------------------------------------------------------------------- #

fetch_osm_roads_wgs84 <- function(cfg) {
  aoi <- study_aoi_polygon(cfg)
  bb <- sf::st_bbox(aoi)
  q <- osmdata::opq(bbox = bb) |>
    osmdata::add_osm_feature(
      key = "highway",
      value = c(
        "motorway", "trunk", "primary", "secondary",
        "tertiary", "residential"
      )
    )
  od <- osmdata::osmdata_sf(q)
  ln <- od$osm_lines
  if (is.null(ln) || nrow(ln) == 0) return(NULL)
  sf::st_transform(ln, 4326)
}

# --------------------------------------------------------------------------- #
# Leaflet dashboard
# --------------------------------------------------------------------------- #

aggregate_categorical_raster <- function(r, fact) {
  if (fact < 2L) return(r)
  tryCatch(
    terra::aggregate(r, fact = fact, fun = "modal", na.rm = TRUE),
    error = function(e) terra::aggregate(r, fact = fact, fun = stats::median, na.rm = TRUE)
  )
}

save_leaflet_dashboard <- function(cfg, fact = 4L) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    message("Install htmlwidgets for dashboard HTML: install.packages(\"htmlwidgets\")")
    return(invisible(NULL))
  }

  y_last <- max(phase_years_sorted(cfg))
  lu_path <- path_lulc(y_last, cfg)
  y0 <- min(phase_years_sorted(cfg))
  ch_path <- path_change_bin(y0, y_last, cfg)

  if (!file.exists(lu_path)) {
    warning("Leaflet: missing ", lu_path, call. = FALSE)
    return(invisible(NULL))
  }

  lu <- prep_lulc_for_tmap(terra::rast(lu_path))
  lu_w <- lu |>
    aggregate_categorical_raster(fact) |>
    terra::project("EPSG:4326", method = "near", gdal = TRUE)

  ch_w <- NULL
  if (file.exists(ch_path)) {
    ch <- terra::rast(ch_path)
    ch <- terra::resample(ch, lu, method = "near")
    ch_w <- ch |>
      aggregate_categorical_raster(fact) |>
      terra::project("EPSG:4326", method = "near", gdal = TRUE)
  }

  rings_path <- file.path(resolved_proc_dir(cfg), "buffer_rings.shp")
  rings_w <- if (file.exists(rings_path)) {
    sf::st_read(rings_path, quiet = TRUE) |> sf::st_transform(4326)
  } else {
    NULL
  }

  hot_path <- file.path(resolved_proc_dir(cfg), "sprawl_hotspots.shp")
  hot_w <- if (file.exists(hot_path)) {
    sf::st_read(hot_path, quiet = TRUE) |> sf::st_transform(4326)
  } else {
    NULL
  }

  ueii_path <- file.path(resolved_proc_dir(cfg), "ueii_results.csv")
  if (!is.null(hot_w) && file.exists(ueii_path)) {
    ue <- utils::read.csv(ueii_path, stringsAsFactors = FALSE)
    if (nrow(ue)) {
      ue <- ue[order(ue$post_year, ue$pre_year, decreasing = TRUE), , drop = FALSE]
      u1 <- ue[1, , drop = FALSE]
      ue_sub <- ue[
        ue$pre_year == u1$pre_year & ue$post_year == u1$post_year,
        c("zone_id", "ueii"),
        drop = FALSE
      ]
      hot_w <- merge(hot_w, ue_sub, by = "zone_id", all.x = TRUE)
    }
  }

  roads_w <- tryCatch(fetch_osm_roads_wgs84(cfg), error = function(e) NULL)

  pal_lulc <- leaflet::colorFactor(
    palette = lulc_palette_leaflet,
    domain = 0:4,
    na.color = "transparent"
  )

  m <- leaflet::leaflet(
    width = "100%",
    height = "calc(100vh - 58px)",
    options = leaflet::leafletOptions(preferCanvas = TRUE)
  ) |>
    leaflet::addProviderTiles("CartoDB.Positron", group = "Light basemap") |>
    leaflet::addProviderTiles("Esri.WorldImagery", group = "Satellite") |>
    leaflet::addRasterImage(
      lu_w,
      colors = pal_lulc,
      opacity = 0.72,
      group = paste0("LULC ", y_last),
      project = FALSE
    )

  if (!is.null(ch_w)) {
    pal_ch <- leaflet::colorFactor(
      c("transparent", "#E74C3C"),
      domain = 0:1,
      na.color = "transparent"
    )
    m <- m |>
      leaflet::addRasterImage(
        ch_w,
        colors = pal_ch,
        opacity = 0.78,
        group = paste0("New urban ", y0, "–", y_last),
        project = FALSE
      )
  }

  if (!is.null(rings_w)) {
    lbl <- if ("inner_km" %in% names(rings_w) && "outer_km" %in% names(rings_w)) {
      sprintf(
        "Zone %s<br>Ring: %s–%s km",
        rings_w$zone_id,
        rings_w$inner_km,
        rings_w$outer_km
      )
    } else {
      sprintf("Zone %s", rings_w$zone_id)
    }
    m <- m |>
      leaflet::addPolygons(
        data = rings_w,
        fillColor = "transparent",
        color = "#8E44AD",
        weight = 2,
        dashArray = "4,4",
        opacity = 0.85,
        fillOpacity = 0,
        popup = lbl,
        group = "Buffer rings"
      )
  }

  if (!is.null(hot_w)) {
    pop <- if ("ueii" %in% names(hot_w)) {
      sprintf(
        "Sprawl hotspot (zone %s)<br>UEII: %s %% / yr",
        hot_w$zone_id,
        ifelse(is.na(hot_w$ueii), "—", format(round(hot_w$ueii, 2), nsmall = 2))
      )
    } else {
      sprintf("Sprawl hotspot (zone %s)", hot_w$zone_id)
    }
    m <- m |>
      leaflet::addPolygons(
        data = hot_w,
        fillColor = "#E74C3C",
        fillOpacity = 0.28,
        color = "#C0392B",
        weight = 2,
        popup = pop,
        group = "Sprawl hotspots"
      )
  }

  if (!is.null(roads_w) && nrow(roads_w)) {
    m <- m |>
      leaflet::addPolylines(
        data = roads_w,
        color = "#2C3E50",
        weight = 1,
        opacity = 0.55,
        group = "Road network"
      )
  }

  overlay_groups <- c(
    paste0("LULC ", y_last),
    if (!is.null(ch_w)) paste0("New urban ", y0, "–", y_last),
    if (!is.null(rings_w)) "Buffer rings",
    if (!is.null(hot_w)) "Sprawl hotspots",
    if (!is.null(roads_w) && nrow(roads_w)) "Road network"
  )

  m <- m |>
    leaflet::addLegend(
      position = "bottomright",
      colors = lulc_palette_leaflet,
      labels = lulc_class_labels,
      title = paste0("Land cover (", y_last, ")"),
      opacity = 0.9
    ) |>
    leaflet::addLayersControl(
      baseGroups = c("Light basemap", "Satellite"),
      overlayGroups = overlay_groups,
      options = leaflet::layersControlOptions(collapsed = FALSE)
    ) |>
    leaflet::addScaleBar(position = "bottomleft")

  if (requireNamespace("leaflet.extras", quietly = TRUE)) {
    m <- leaflet.extras::addFullscreenControl(m)
  }

  e <- terra::ext(lu_w)
  ex <- c(terra::xmin(e), terra::xmax(e), terra::ymin(e), terra::ymax(e))
  m <- m |>
    leaflet::fitBounds(ex[1], ex[3], ex[2], ex[4])

  out_html <- file.path(resolved_out_dir(cfg), "urban_sprawl_dashboard.html")
  dash_title <- paste0(cfg$study_area$name %||% "Study area", " — Urban sprawl explorer")

  if (requireNamespace("htmltools", quietly = TRUE)) {
    study_nm <- htmltools::htmlEscape(cfg$study_area$name %||% "Study area")
    dash_css <- paste0(
      "html,body{margin:0;padding:0;height:100%;font-family:system-ui,-apple-system,'Segoe UI',Roboto,sans-serif;}",
      "body{padding-top:58px;box-sizing:border-box;background:#eceff1;}",
      "#urban-dash-topbar{position:fixed;top:0;left:0;right:0;height:58px;padding:6px 1.1rem 8px;",
      "display:flex;flex-direction:column;justify-content:center;",
      "background:linear-gradient(135deg,#1a237e 0%,#3949ab 100%);color:#fff;",
      "z-index:20000;box-shadow:0 2px 14px rgba(0,0,0,.22);}",
      "#urban-dash-topbar .t1{margin:0;font-size:1.05rem;font-weight:600;line-height:1.2;}",
      "#urban-dash-topbar .t2{margin:3px 0 0;font-size:.78rem;opacity:.93;line-height:1.35;max-width:56rem;}",
      ".leaflet-control-layers{border-radius:12px!important;box-shadow:0 4px 22px rgba(0,0,0,.14)!important;",
      "border:1px solid rgba(0,0,0,.06)!important;}",
      ".leaflet-control-layers-expanded{padding:12px 14px!important;}",
      ".leaflet-control-layers-base label,.leaflet-control-layers-overlays label{font-size:13px;margin-bottom:5px;}",
      ".info.legend{border-radius:10px!important;box-shadow:0 2px 12px rgba(0,0,0,.1)!important;padding:8px 10px!important;}"
    )
    topbar <- htmltools::tags$div(
      id = "urban-dash-topbar",
      htmltools::tags$p(class = "t1", htmltools::HTML(paste0(study_nm, " — Urban sprawl explorer"))),
      htmltools::tags$p(
        class = "t2",
        "Use the layer control (top-left) for LULC, new urban (CNN), buffer rings, hotspots, and roads. Choose Light or Satellite basemap below."
      )
    )
    m <- htmlwidgets::prependContent(
      m,
      htmltools::tagList(
        htmltools::tags$style(htmltools::HTML(dash_css)),
        topbar
      )
    )
  }

  htmlwidgets::saveWidget(
    m,
    file = out_html,
    selfcontained = TRUE,
    title = dash_title
  )
  message("Wrote ", out_html)
  invisible(out_html)
}

# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

#' @param render_report If TRUE, renders `docs/urban_sprawl_report.Rmd` to **PDF** (primary) and **HTML**.
#'   PDF needs a LaTeX engine (e.g. `install.packages("tinytex"); tinytex::install_tinytex()`). If PDF fails, HTML is still attempted.
run_phase7 <- function(cfg = load_config(), render_report = TRUE) {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Install scales: install.packages(\"scales\")", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2 for Phase 7 charts: install.packages(\"ggplot2\")", call. = FALSE)
  }
  pr <- normalizePath(project_root(), winslash = "/", mustWork = FALSE)
  cfg_chk <- file.path(pr, "config", "study_area.yml")
  if (!file.exists(cfg_chk)) {
    stop(
      "Phase 7 needs the project root as working directory (missing ",
      cfg_chk,
      "). Use setwd() to your Urban Sprawl folder or open UrbanSprawl.Rproj, then run run_phase7().",
      call. = FALSE
    )
  }

  with_pipeline_phase("Phase 7 — maps, charts, dashboard (and optional report)", {
  p7 <- get_phase7_cfg(cfg)
  ensure_dirs(cfg)

  pipeline_log("Phase 7: LULC panel map …")
  save_lulc_panel_map(cfg, dpi = p7$map_dpi)
  pipeline_log("Phase 7: change maps …")
  save_change_maps(cfg, dpi = p7$map_dpi)
  pipeline_log("Phase 7: hotspot map …")
  save_hotspot_map(cfg, dpi = p7$map_dpi)

  pipeline_log("Phase 7: charts …")
  save_urban_growth_chart(cfg, dpi = p7$chart_dpi)
  save_entropy_chart(cfg, dpi = p7$chart_dpi)
  save_landscape_metrics_chart(cfg, dpi = p7$chart_dpi)
  save_driver_chart(cfg, dpi = p7$chart_dpi)
  save_correlation_heatmap(cfg, dpi = p7$chart_dpi)

  n_chart_png <- length(Sys.glob(file.path(resolved_out_dir(cfg), "charts", "*.png")))
  n_map_png <- length(Sys.glob(file.path(resolved_out_dir(cfg), "maps", "*.png")))
  pipeline_log("Phase 7: wrote ", n_chart_png, " chart PNG(s) and ", n_map_png, " map PNG(s) under ", resolved_out_dir(cfg))

  pipeline_log("Phase 7: Leaflet dashboard …")
  save_leaflet_dashboard(cfg, fact = p7$leaflet_aggregate_fact)

  if (isTRUE(render_report)) {
    rmd <- file.path(project_root(), "docs", "urban_sprawl_report.Rmd")
    if (!file.exists(rmd)) {
      message("Report template missing: ", rmd)
    } else if (!requireNamespace("rmarkdown", quietly = TRUE)) {
      message("Install rmarkdown to render the report: install.packages(\"rmarkdown\")")
    } else {
      root_n <- normalizePath(project_root(), winslash = "/", mustWork = FALSE)
      out_sub <- cfg$paths$outputs_dir %||% "outputs"
      out_is_abs <- nzchar(out_sub) && (startsWith(out_sub, "/") || grepl("^[A-Za-z]:[\\\\/]", out_sub))
      outputs_root <- if (out_is_abs) {
        normalizePath(out_sub, winslash = "/", mustWork = FALSE)
      } else {
        normalizePath(file.path(root_n, out_sub), winslash = "/", mustWork = FALSE)
      }
      rep_params <- list(
        project_root = root_n,
        study_name = cfg$study_area$name %||% "Study area",
        outputs_dir = out_sub,
        outputs_root = outputs_root
      )
      out_dir <- resolved_out_dir(cfg)

      pipeline_log("Phase 7: research report → PDF (requires LaTeX / TinyTeX) …")
      knit_root <- normalizePath(rep_params$project_root, winslash = "/", mustWork = FALSE)
      out_pdf <- tryCatch(
        rmarkdown::render(
          rmd,
          output_file = "urban_sprawl_report.pdf",
          output_dir = out_dir,
          params = rep_params,
          output_format = "pdf_document",
          knit_root_dir = knit_root,
          quiet = TRUE
        ),
        error = function(e) {
          warning(
            "PDF report failed: ", conditionMessage(e),
            "\nInstall TinyTeX: install.packages(\"tinytex\"); tinytex::install_tinytex()",
            call. = FALSE
          )
          NULL
        }
      )
      if (!is.null(out_pdf)) pipeline_log("Wrote PDF report: ", normalizePath(out_pdf))

      pipeline_log("Phase 7: research report → HTML …")
      out_html <- tryCatch(
        rmarkdown::render(
          rmd,
          output_file = "urban_sprawl_report.html",
          output_dir = out_dir,
          params = rep_params,
          output_format = "html_document",
          knit_root_dir = knit_root,
          quiet = TRUE
        ),
        error = function(e) {
          warning("HTML report failed: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (!is.null(out_html)) pipeline_log("Wrote HTML report: ", normalizePath(out_html))
    }
  }

  pipeline_log("Phase 7 complete. Maps: outputs/maps/ ; charts: outputs/charts/")
  invisible(TRUE)
  })
}

if (interactive()) {
  message("Phase 7: run_phase7() after Phases 3–6 (partial inputs OK; missing files are skipped).")
}
