# ggplot fallbacks for docs/urban_sprawl_report.Rmd when outputs/charts/*.png are missing.
# Sourced into knitr::knit_global() from the report setup chunk.

theme_sprawl_knit <- function() {
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

knit_urban_growth_plot <- function(root, study_name) {
  path <- normalizePath(file.path(root, "data", "processed", "urban_area_stats.csv"),
    winslash = "/", mustWork = FALSE
  )
  if (!file.exists(path)) return(invisible(NULL))
  u <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("Year", "Urban_km2") %in% names(u))) return(invisible(NULL))
  u <- u[order(u$Year), , drop = FALSE]
  nm <- if (is.null(study_name) || !nzchar(as.character(study_name)[1L])) {
    "Study area"
  } else {
    as.character(study_name)[1L]
  }
  ggplot2::ggplot(u, ggplot2::aes(x = Year, y = Urban_km2)) +
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
      subtitle = paste("Study area:", nm),
      x = "Year",
      y = "Urban area (km²)",
      caption = "Rendered at knit time (no Phase 7 PNG); source: urban_area_stats.csv"
    ) +
    ggplot2::scale_x_continuous(breaks = unique(u$Year)) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    theme_sprawl_knit()
}

knit_entropy_plot <- function(root) {
  path <- normalizePath(file.path(root, "data", "processed", "shannon_entropy_results.csv"),
    winslash = "/", mustWork = FALSE
  )
  if (!file.exists(path)) return(invisible(NULL))
  e <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"H_rel" %in% names(e)) return(invisible(NULL))
  ggplot2::ggplot(e, ggplot2::aes(x = year, y = H_rel)) +
    ggplot2::geom_line(color = "#8E44AD", linewidth = 1.2) +
    ggplot2::geom_point(color = "#6C3483", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(H_rel, 3)), vjust = -1.1, size = 3.2) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed", color = "#E74C3C", linewidth = 0.6) +
    ggplot2::labs(
      title = "Urban sprawl intensity (Shannon entropy)",
      subtitle = "Higher H_rel = more dispersed urban across buffer zones",
      x = "Year",
      y = "Relative Shannon entropy (H / ln(n))",
      caption = "Rendered at knit time; source: shannon_entropy_results.csv"
    ) +
    ggplot2::scale_x_continuous(breaks = unique(e$year)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    theme_sprawl_knit()
}

knit_landscape_plot <- function(root) {
  path <- normalizePath(file.path(root, "data", "processed", "landscape_metrics.csv"),
    winslash = "/", mustWork = FALSE
  )
  if (!file.exists(path)) return(invisible(NULL))
  lm <- utils::read.csv(path, stringsAsFactors = FALSE)
  long <- rbind(
    data.frame(year = lm$year, metric = "PD", value = lm$PD),
    data.frame(year = lm$year, metric = "LPI", value = lm$LPI),
    data.frame(year = lm$year, metric = "ED", value = lm$ED)
  )
  long <- long[stats::complete.cases(long), , drop = FALSE]
  if (!nrow(long)) return(invisible(NULL))
  ggplot2::ggplot(long, ggplot2::aes(x = year, y = value, color = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = "Landscape structure metrics (urban class)",
      subtitle = "Patch density (PD), largest patch index (LPI), edge density (ED)",
      x = "Year",
      y = "Value",
      caption = "Rendered at knit time; source: landscape_metrics.csv"
    ) +
    theme_sprawl_knit() +
    ggplot2::theme(legend.position = "none")
}

knit_drivers_plot <- function(oroot) {
  path <- normalizePath(file.path(oroot, "logistic_regression_results.csv"),
    winslash = "/", mustWork = FALSE
  )
  if (!file.exists(path)) return(invisible(NULL))
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("term", "estimate", "std_error") %in% names(d))) return(invisible(NULL))
  d <- d[d$term != "(Intercept)", , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  d$coefficient <- d$estimate
  d$ci_lower <- d$estimate - 1.96 * d$std_error
  d$ci_upper <- d$estimate + 1.96 * d$std_error
  d$direction <- ifelse(d$coefficient > 0, "Increases P(growth)", "Decreases P(growth)")
  d$predictor <- stats::reorder(d$term, abs(d$coefficient))
  ggplot2::ggplot(d, ggplot2::aes(x = predictor, y = coefficient, fill = direction)) +
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
      caption = "Rendered at knit time; source: logistic_regression_results.csv"
    ) +
    theme_sprawl_knit() +
    ggplot2::theme(legend.position = "bottom")
}

knit_corr_plot <- function(oroot) {
  if (!requireNamespace("ggcorrplot", quietly = TRUE)) return(invisible(NULL))
  path <- normalizePath(file.path(oroot, "correlation_matrix.csv"),
    winslash = "/", mustWork = FALSE
  )
  if (!file.exists(path)) return(invisible(NULL))
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"variable" %in% names(raw)) return(invisible(NULL))
  rn <- raw$variable
  mat <- as.matrix(raw[, setdiff(names(raw), "variable"), drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- rn
  ggcorrplot::ggcorrplot(
    mat,
    method = "circle",
    type = "lower",
    lab = TRUE,
    lab_size = 3,
    colors = c("#3498DB", "white", "#E74C3C"),
    title = "Correlation matrix — sprawl drivers",
    ggtheme = theme_sprawl_knit()
  )
}
