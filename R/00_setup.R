#' Phase 0 — Environment and paths
#'
#' Phase 1 uses Microsoft Planetary Computer (rstac + httr + jsonlite). Optional: rgee.

suppressPackageStartupMessages({
  pkgs <- c(
    "terra", "sf", "rstac", "httr", "jsonlite", "yaml",
    "randomForest", "e1071", "caret", "torch", "landscapemetrics",
    "mapedit", "osmdata", "tmap", "ggplot2", "leaflet", "scales"
  )
  to_install <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
})

project_root <- function() {
  getwd()
}

load_config <- function(path = NULL) {
  if (is.null(path)) {
    path <- file.path(project_root(), "config", "study_area.yml")
  }
  yaml::read_yaml(path)
}

# --- Terminal progress (stderr + flush: visible in RStudio & Rscript; not buffered like some message() sinks)

pipeline_log <- function(...) {
  line <- paste0(
    format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "),
    paste0(..., collapse = ""),
    "\n"
  )
  cat(line, file = stderr())
  flush(stderr())
}

pipeline_begin <- function(title) {
  pipeline_log("")
  pipeline_log("============================================================")
  pipeline_log("  >>> START: ", title)
  pipeline_log("  >>> (Long pauses are normal while rasters download or compute.)")
  pipeline_log("============================================================")
}

pipeline_ok <- function(title) {
  pipeline_log("============================================================")
  pipeline_log("  <<< FINISHED OK: ", title)
  pipeline_log("============================================================")
  pipeline_log("")
}

pipeline_stop_banner <- function(title) {
  pipeline_log("============================================================")
  pipeline_log("  !!! STOPPED: ", title)
  pipeline_log("  !!! See error message just above (or run traceback() in R).")
  pipeline_log("============================================================")
}

#' Run a phase body with visible start/end banners and immediate warnings.
with_pipeline_phase <- function(title, expr) {
  pipeline_begin(title)
  succeeded <- FALSE
  owarn <- getOption("warn")
  on.exit(
    {
      options(warn = owarn)
      if (!succeeded) pipeline_stop_banner(title)
    },
    add = TRUE
  )
  options(warn = 1)
  out <- eval(substitute(expr), parent.frame())
  succeeded <- TRUE
  pipeline_ok(title)
  out
}

init_gee <- function() {
  if (!requireNamespace("rgee", quietly = TRUE)) stop("Install rgee for optional GEE workflows")
  rgee::ee_Initialize()
}
