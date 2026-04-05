# One-time: install CRAN dependencies for the urban sprawl pipeline.
#
# Windows tips:
# - Close all other R / RStudio sessions if you see "Permission denied" or 00LOCK on a package.
# - Rtools is only needed to *compile* packages from source; CRAN binaries usually work without it.
# - Large downloads: this script raises the timeout. If a package still fails, re-run or install it alone:
#     install.packages("terra", repos = "https://cloud.r-project.org")

options(timeout = max(900, getOption("timeout")))

repos <- "https://cloud.r-project.org"

# Large / common timeout failures — install first, one-by-one if needed
heavy <- c("terra", "sf", "torch")
for (p in heavy) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " ...")
    tryCatch(
      install.packages(p, repos = repos, type = "binary"),
      error = function(e) message("Failed: ", p, " — try again: install.packages('", p, "')")
    )
  }
}

# Everything else (leaflet.extras omitted: often not yet built for bleeding-edge R — Phase 7 works without it)
pkgs <- c(
  "rstac", "httr", "jsonlite", "yaml",
  "randomForest", "e1071", "caret",
  "landscapemetrics", "mapedit", "osmdata", "tmap", "ggplot2",
  "leaflet", "stars", "scales",
  "ggcorrplot", "htmlwidgets", "rmarkdown", "base64enc"
)
miss <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(miss)) {
  install.packages(miss, repos = repos, type = "binary")
}

# Optional: fullscreen control on the Leaflet dashboard (skip if not on CRAN for your R version)
if (!requireNamespace("leaflet.extras", quietly = TRUE)) {
  message(
    "Optional package leaflet.extras not installed (CRAN may not have a build for this R version). ",
    "Phase 7 dashboard still works; only the fullscreen button is skipped.\n",
    "Try later: install.packages('leaflet.extras', repos = 'https://cloud.r-project.org')"
  )
}

message(
  "Done.\n",
  "Phase 4 CNN: library(torch); if needed: torch::install_torch()\n",
  "Phase 1: source(\"R/01_aoi_datasets.R\"); run_phase1()\n",
  "PDF research report (Phase 7): install.packages(\"tinytex\"); tinytex::install_tinytex()\n",
  "Optional GEE: install.packages(\"rgee\")"
)
