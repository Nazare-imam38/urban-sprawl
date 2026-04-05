#' Phase 4 — Siamese CNN urban change detection (torch)
#'
#' Change reference: (lulc_pre != 1) & (lulc_post == 1). Patches 64×64, stride 32,
#' label = majority of centre 8×8 on change_ref. Spatial 4×4 block split (11 train / 3 val / 2 test).
#' Imbalance: subsample unchanged + BCEWithLogitsLoss(pos_weight). Full-raster inference: batched torch
#' forward + cached feature matrices; default infer stride 32 (set 16 in config for denser maps).
#' Post: MMU + mask water (LULC post == 4).
#'
#' Requires **torch**. Working directory = project root.

source("R/00_setup.R")

if (!requireNamespace("torch", quietly = TRUE)) {
  stop("Install torch: install.packages(\"torch\"); then torch::install_torch() if needed.", call. = FALSE)
}

# --------------------------------------------------------------------------- #
# Paths & config
# --------------------------------------------------------------------------- #

path_features <- function(year, cfg) {
  file.path(cfg$paths$processed_dir, sprintf("features_%d.tif", year))
}

path_lulc <- function(year, cfg) {
  file.path(cfg$paths$processed_dir, sprintf("lulc_%d.tif", year))
}

path_change_ref <- function(y1, y2, cfg) {
  file.path(cfg$paths$processed_dir, sprintf("change_ref_%d_%d.tif", y1, y2))
}

path_change_prob <- function(y1, y2, cfg) {
  file.path(cfg$paths$outputs_dir, sprintf("change_prob_%d_%d.tif", y1, y2))
}

path_change_bin <- function(y1, y2, cfg) {
  file.path(cfg$paths$outputs_dir, sprintf("change_map_%d_%d.tif", y1, y2))
}

path_cnn_checkpoint <- function(cfg) {
  file.path(cfg$paths$outputs_dir, "cnn_model_best.pt")
}

path_cnn_report <- function(cfg) {
  file.path(cfg$paths$outputs_dir, "cnn_accuracy_report.csv")
}

path_patches_dir <- function(cfg) {
  file.path(cfg$paths$processed_dir, "patches")
}

get_cnn_cfg <- function(cfg) {
  cnn <- cfg$cnn
  if (is.null(cnn)) stop("config: add `cnn` block to study_area.yml", call. = FALSE)
  cnn
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------------------------------------------------------- #
# Canonical 11-channel order (Sentinel names) for both sensors
# --------------------------------------------------------------------------- #

CANON_FEATURE_NAMES <- c(
  "B02", "B03", "B04", "B08", "B11", "B12",
  "NDVI", "NDBI", "MNDWI", "NBI", "BSI"
)

LANDSAT_TO_CANON <- c(
  SR_B2 = "B02", SR_B3 = "B03", SR_B4 = "B04", SR_B5 = "B08",
  SR_B6 = "B11", SR_B7 = "B12",
  NDVI = "NDVI", NDBI = "NDBI", MNDWI = "MNDWI", NBI = "NBI", BSI = "BSI"
)

canonicalize_features <- function(r) {
  nm <- names(r)
  if ("SR_B2" %in% nm) {
    lst <- vector("list", length(CANON_FEATURE_NAMES))
    names(lst) <- CANON_FEATURE_NAMES
    for (i in seq_along(LANDSAT_TO_CANON)) {
      src <- names(LANDSAT_TO_CANON)[i]
      dst <- LANDSAT_TO_CANON[i]
      lst[[dst]] <- r[[src]]
    }
    return(terra::rast(lst))
  }
  miss <- setdiff(CANON_FEATURE_NAMES, nm)
  if (length(miss)) stop("Feature stack missing layers: ", paste(miss, collapse = ", "), call. = FALSE)
  r[[CANON_FEATURE_NAMES]]
}

align_to_template <- function(r, tmpl, method = "bilinear") {
  if (terra::compareGeom(r, tmpl, stopOnError = FALSE, rowcol = TRUE, ext = TRUE, crs = TRUE)) {
    return(r)
  }
  terra::resample(r, tmpl, method = method)
}

align_lulc <- function(lul, tmpl) {
  align_to_template(lul, tmpl, method = "near")
}

#' @param round_lulc If TRUE, round to integers 1–4 before change logic
build_change_ref <- function(lulc_pre, lulc_post, round_lulc = TRUE) {
  a <- lulc_pre
  b <- lulc_post
  if (round_lulc) {
    a <- terra::round(a)
    b <- terra::round(b)
  }
  ref <- (a != 1L) & (b == 1L)
  ref <- terra::ifel(ref, 1L, 0L)
  names(ref) <- "change"
  ref
}

# --------------------------------------------------------------------------- #
# Spatial blocks (4×4) → train / val / test (11 / 3 / 2)
# --------------------------------------------------------------------------- #

block_id_from_rc <- function(row_i, col_j, nr, nc) {
  br <- pmin(4L, pmax(1L, ceiling(4 * row_i / nr)))
  bc <- pmin(4L, pmax(1L, ceiling(4 * col_j / nc)))
  as.integer((br - 1L) * 4L + bc)
}

assign_blocks_to_splits <- function(seed = 42L) {
  set.seed(seed)
  lab <- sample(c(rep("train", 11L), rep("val", 3L), rep("test", 2L)))
  # Named list: reliable [[ "1" ]] … [[ "16" ]] (named atomic vectors can be picky in edge cases)
  stats::setNames(as.list(lab), as.character(1L:16L))
}

#' Seeds for 4×4 spatial train/val/test assignment (retry if change pixels sit only in val/test blocks)
phase4_spatial_split_seeds <- function(cfg) {
  cnn <- get_cnn_cfg(cfg)
  base <- as.integer(cnn$spatial_split_seed %||% 42L)
  fb <- cnn$spatial_split_seeds_fallback
  if (length(fb)) {
    return(unique(as.integer(c(base, unlist(fb, use.names = FALSE)))))
  }
  unique(c(
    base,
    99L, 123L, 7L, 2024L, 11L, 777L, 31415L, 2018L, 2010L, 2000L,
    555L, 13L, 21L, 500L, 1024L, 3L, 17L, 19L, 23L, 29L
  ))
}

#' All blocks labeled train — only for patch extraction; real splits applied via [assign_patches_to_spatial_splits]
phase4_dummy_train_splits <- function() {
  stats::setNames(as.list(rep("train", 16L)), as.character(1L:16L))
}

patch_center_rc <- function(r0, c0, ps) {
  c(r0 + (ps %/% 2L), c0 + (ps %/% 2L))
}

# --------------------------------------------------------------------------- #
# Patch extraction (arrays: channel × H × W)
# --------------------------------------------------------------------------- #

raster_to_matrix <- function(r) {
  stopifnot(inherits(r, "SpatRaster"), terra::nlyr(r) == 1L)
  # terra::values(..., mat=TRUE) is ncell × nlyr (columns = layers), NOT nrow × ncol.
  # Reshaping that with matrix(..., byrow=FALSE) mis-orders row-major cell layout.
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  m <- tryCatch(
    terra::as.matrix(r, wide = TRUE),
    error = function(e) NULL
  )
  if (!is.matrix(m) || nrow(m) != nr || ncol(m) != nc) {
    v <- as.vector(terra::values(r, mat = FALSE))
    m <- matrix(v, nrow = nr, ncol = nc, byrow = TRUE)
  }
  m
}

#' One full read per layer (call once per stack, then use [extract_patch_chw_from_mats])
feature_stack_as_matrices <- function(r_stack) {
  lapply(seq_len(terra::nlyr(r_stack)), function(L) raster_to_matrix(r_stack[[L]]))
}

extract_patch_chw_from_mats <- function(mats, r0, c0, ps) {
  m1 <- mats[[1L]]
  nr <- nrow(m1)
  nc <- ncol(m1)
  if (r0 < 1L || c0 < 1L || r0 + ps - 1L > nr || c0 + ps - 1L > nc) {
    return(NULL)
  }
  nl <- length(mats)
  a <- array(NA_real_, dim = c(nl, ps, ps))
  for (L in seq_len(nl)) {
    m <- mats[[L]]
    a[L, , ] <- m[r0:(r0 + ps - 1L), c0:(c0 + ps - 1L)]
  }
  a
}

extract_patch_chw <- function(r_stack, r0, c0, ps) {
  extract_patch_chw_from_mats(feature_stack_as_matrices(r_stack), r0, c0, ps)
}

centre_label_from_ref_mat <- function(ref_mat, r0, c0, ps, centre = 8L, thr = 0.5) {
  o <- (ps - centre) %/% 2L
  rr <- (r0 + o):(r0 + o + centre - 1L)
  cc <- (c0 + o):(c0 + o + centre - 1L)
  v <- ref_mat[rr, cc]
  mean(v, na.rm = TRUE) >= thr
}

centre_label_from_ref <- function(ref, r0, c0, ps, centre = 8L, thr = 0.5) {
  centre_label_from_ref_mat(raster_to_matrix(ref), r0, c0, ps, centre = centre, thr = thr)
}

extract_patches_transition <- function(
    feat_pre, feat_post, ref,
    ps, stride, block_splits, nr, nc, centre = 8L, centre_thr = 0.5) {
  message("Phase 4: loading feature rasters into memory for patch extraction …")
  mats_pre <- feature_stack_as_matrices(feat_pre)
  mats_post <- feature_stack_as_matrices(feat_post)
  ref_mat <- raster_to_matrix(ref)

  nrp <- nrow(mats_pre[[1L]])
  ncp <- ncol(mats_pre[[1L]])
  nrs <- nrow(mats_post[[1L]])
  ncs <- ncol(mats_post[[1L]])
  nrr <- nrow(ref_mat)
  ncr <- ncol(ref_mat)
  if (nrp != nr || ncp != nc || nrs != nr || ncs != nc || nrr != nr || ncr != nc) {
    stop(
      "Feature/ref matrix dims do not match template (", nr, "×", nc, "). ",
      "pre ", nrp, "×", ncp, ", post ", nrs, "×", ncs, ", ref ", nrr, "×", ncr,
      ". Re-run Phase 2–3 alignment or check rasters.",
      call. = FALSE
    )
  }

  pre_l <- list()
  post_l <- list()
  y_l <- list()
  sp_l <- list()
  r0_l <- list()
  c0_l <- list()
  n_iter <- 0L
  n_skip_spl <- 0L
  n_skip_bounds <- 0L
  n_skip_na <- 0L

  for (r0 in seq(1L, nr - ps + 1L, by = stride)) {
    for (c0 in seq(1L, nc - ps + 1L, by = stride)) {
      n_iter <- n_iter + 1L
      cc <- patch_center_rc(r0, c0, ps)
      bid <- as.character(block_id_from_rc(cc[1], cc[2], nr, nc))
      spl <- block_splits[[bid]]
      if (is.null(spl)) {
        n_skip_spl <- n_skip_spl + 1L
        next
      }

      pa <- extract_patch_chw_from_mats(mats_pre, r0, c0, ps)
      pb <- extract_patch_chw_from_mats(mats_post, r0, c0, ps)
      if (is.null(pa) || is.null(pb)) {
        n_skip_bounds <- n_skip_bounds + 1L
        next
      }
      if (any(is.na(pa)) || any(is.na(pb))) {
        n_skip_na <- n_skip_na + 1L
        next
      }

      lab <- centre_label_from_ref_mat(ref_mat, r0, c0, ps, centre = centre, thr = centre_thr)
      pre_l[[length(pre_l) + 1L]] <- pa
      post_l[[length(post_l) + 1L]] <- pb
      y_l[[length(y_l) + 1L]] <- lab
      sp_l[[length(sp_l) + 1L]] <- spl
      r0_l[[length(r0_l) + 1L]] <- r0
      c0_l[[length(c0_l) + 1L]] <- c0
    }
  }

  if (!length(pre_l)) {
    stop(
      "No patches extracted. Template ", nr, "×", nc, ", patch ", ps, ", stride ", stride,
      ". Grid cells visited: ", n_iter,
      ". Skipped — no split name: ", n_skip_spl, ", out of bounds: ", n_skip_bounds,
      ", NA in features: ", n_skip_na,
      ". If `no split name` ≈ all cells, block_splits lookup failed.",
      call. = FALSE
    )
  }

  n <- length(pre_l)
  ch <- dim(pre_l[[1]])[1]
  pre_arr <- array(NA_real_, dim = c(n, ch, ps, ps))
  post_arr <- array(NA_real_, dim = c(n, ch, ps, ps))
  yv <- logical(n)
  for (i in seq_len(n)) {
    pre_arr[i, , , ] <- pre_l[[i]]
    post_arr[i, , , ] <- post_l[[i]]
    yv[i] <- y_l[[i]]
  }
  list(
    pre = pre_arr,
    post = post_arr,
    y = yv,
    split = unlist(sp_l),
    r0 = unlist(r0_l, use.names = FALSE),
    c0 = unlist(c0_l, use.names = FALSE)
  )
}

assign_patches_to_spatial_splits <- function(r0v, c0v, nr, nc, ps, block_splits) {
  vapply(seq_along(r0v), function(i) {
    cc <- patch_center_rc(r0v[i], c0v[i], ps)
    bid <- as.character(block_id_from_rc(cc[1], cc[2], nr, nc))
    spl <- block_splits[[bid]]
    if (is.null(spl)) {
      return(NA_character_)
    }
    as.character(spl[[1L]])
  }, character(1L))
}

balance_indices <- function(y, ratio_unc_to_chg = 3L, seed = 42L) {
  set.seed(seed)
  idx_chg <- which(y)
  idx_unc <- which(!y)
  if (!length(idx_chg)) {
    stop("No positive (change) patches — widen AOI, date pair, or stride.", call. = FALSE)
  }
  n_need_unc <- min(length(idx_unc), ratio_unc_to_chg * length(idx_chg))
  idx_unc_s <- sample(idx_unc, n_need_unc)
  sample(c(idx_chg, idx_unc_s))
}

split_by_spatial_label <- function(obj, split_name) {
  w <- obj$split == split_name
  list(
    pre = obj$pre[w, , , , drop = FALSE],
    post = obj$post[w, , , , drop = FALSE],
    y = obj$y[w]
  )
}

# --------------------------------------------------------------------------- #
# torch model — Siamese encoder + MLP head (logits)
# --------------------------------------------------------------------------- #

encoder_module <- function(in_ch = 11L) {
  torch::nn_module(
    classname = "change_encoder",
    initialize = function(in_channels) {
      self$c1 <- torch::nn_conv2d(in_channels, 32L, 3L, padding = 1L)
      self$bn1 <- torch::nn_batch_norm2d(32L)
      self$c2 <- torch::nn_conv2d(32L, 64L, 3L, padding = 1L)
      self$bn2 <- torch::nn_batch_norm2d(64L)
      self$c3 <- torch::nn_conv2d(64L, 128L, 3L, padding = 1L)
      self$bn3 <- torch::nn_batch_norm2d(128L)
    },
    forward = function(x) {
      x <- torch::nnf_relu(self$bn1(self$c1(x)))
      x <- torch::nnf_max_pool2d(x, 2L)
      x <- torch::nnf_relu(self$bn2(self$c2(x)))
      x <- torch::nnf_max_pool2d(x, 2L)
      x <- torch::nnf_relu(self$bn3(self$c3(x)))
      x <- torch::nnf_max_pool2d(x, 2L)
      torch::torch_flatten(x, start_dim = 2L)
    }
  )(in_ch)
}

siamese_change_cnn <- function(in_ch = 11L) {
  torch::nn_module(
    classname = "siamese_change_cnn",
    initialize = function(in_channels) {
      self$encoder <- encoder_module(in_channels)
      self$fc1 <- torch::nn_linear(16384L, 256L)
      self$do1 <- torch::nn_dropout(0.5)
      self$fc2 <- torch::nn_linear(256L, 64L)
      self$do2 <- torch::nn_dropout(0.3)
      self$fc3 <- torch::nn_linear(64L, 1L)
    },
    forward = function(x_pre, x_post) {
      z1 <- self$encoder(x_pre)
      z2 <- self$encoder(x_post)
      z <- torch::torch_cat(list(z1, z2), dim = 2L)
      z <- torch::nnf_relu(self$fc1(z))
      z <- self$do1(z)
      z <- torch::nnf_relu(self$fc2(z))
      z <- self$do2(z)
      self$fc3(z)
    }
  )(in_ch)
}

array_to_torch_nchw <- function(arr) {
  # arr: N x C x H x W
  torch::torch_tensor(arr, dtype = torch::torch_float())
}

# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #

binary_metrics <- function(pred_prob, truth_bool, thr = 0.5) {
  pred <- pred_prob >= thr
  truth <- truth_bool
  tp <- sum(pred & truth)
  fp <- sum(pred & !truth)
  fn <- sum(!pred & truth)
  tn <- sum(!pred & !truth)
  prec <- tp / (tp + fp + 1e-9)
  rec <- tp / (tp + fn + 1e-9)
  f1 <- 2 * prec * rec / (prec + rec + 1e-9)
  iou <- tp / (tp + fp + fn + 1e-9)
  c(F1 = f1, Precision = prec, Recall = rec, IoU = iou, TP = tp, FP = fp, FN = fn, TN = tn)
}

torch_probs <- function(model, pre_t, post_t) {
  model$eval()
  torch::with_no_grad({
    logits <- model(pre_t, post_t)
    torch::torch_sigmoid(logits)$squeeze(2)
  })
}

# --------------------------------------------------------------------------- #
# Training
# --------------------------------------------------------------------------- #

train_siamese_one_epoch <- function(model, optim, idx, pre_arr, post_arr, y_log, pw, batch_size, grad_clip, device) {
  model$train()
  n <- length(idx)
  perm <- sample.int(n)
  tot_loss <- 0
  nb <- 0L
  for (start in seq(1L, n, by = batch_size)) {
    end <- min(n, start + batch_size - 1L)
    b <- perm[start:end]
    pre_b <- array_to_torch_nchw(pre_arr[b, , , , drop = FALSE])$to(device = device)
    post_b <- array_to_torch_nchw(post_arr[b, , , , drop = FALSE])$to(device = device)
    y_b <- torch::torch_tensor(as.numeric(y_log[b]), dtype = torch::torch_float())$to(device = device)
    optim$zero_grad()
    logits <- model(pre_b, post_b)$squeeze(2)
    loss <- torch::nnf_binary_cross_entropy_with_logits(logits, y_b, pos_weight = pw)
    loss$backward()
    torch::nn_utils_clip_grad_norm_(model$parameters, grad_clip)
    optim$step()
    tot_loss <- tot_loss + as.numeric(loss$item())
    nb <- nb + 1L
  }
  tot_loss / max(1L, nb)
}

eval_siamese_loss_f1 <- function(model, idx, pre_arr, post_arr, y_log, pw, batch_size, device) {
  model$eval()
  n <- length(idx)
  if (!n) return(list(loss = NA_real_, probs = numeric()))
  tot_loss <- 0
  probs <- numeric(n)
  torch::with_no_grad({
    o <- 1L
    for (start in seq(1L, n, by = batch_size)) {
      end <- min(n, start + batch_size - 1L)
      b <- start:end
      pre_b <- array_to_torch_nchw(pre_arr[b, , , , drop = FALSE])$to(device = device)
      post_b <- array_to_torch_nchw(post_arr[b, , , , drop = FALSE])$to(device = device)
      y_b <- torch::torch_tensor(as.numeric(y_log[b]), dtype = torch::torch_float())$to(device = device)
      logits <- model(pre_b, post_b)$squeeze(2)
      loss <- torch::nnf_binary_cross_entropy_with_logits(logits, y_b, pos_weight = pw)
      tot_loss <- tot_loss + as.numeric(loss$item())
      probs[o:(o + length(b) - 1L)] <- as.numeric(torch::torch_sigmoid(logits)$cpu())
      o <- o + length(b)
    }
  })
  list(loss = tot_loss / ceiling(n / batch_size), probs = probs)
}

train_siamese_model <- function(
    train_obj, val_obj, cnn_cfg, checkpoint_path, seed = 42L, device = "cpu") {
  set.seed(seed)
  torch::torch_manual_seed(seed)

  ch <- dim(train_obj$pre)[2]
  ps <- dim(train_obj$pre)[3]
  stopifnot(ps == cnn_cfg$patch_size %||% 64L)

  if (!dim(train_obj$pre)[1]) stop("No training patches in train split.", call. = FALSE)
  if (!dim(val_obj$pre)[1]) stop("No validation patches — check spatial blocks / AOI.", call. = FALSE)

  idx_tr <- balance_indices(train_obj$y, cnn_cfg$balance_unchanged_to_changed %||% 3L, seed)
  pre_tr <- train_obj$pre[idx_tr, , , , drop = FALSE]
  post_tr <- train_obj$post[idx_tr, , , , drop = FALSE]
  y_tr <- train_obj$y[idx_tr]

  pre_va <- val_obj$pre
  post_va <- val_obj$post
  y_va <- val_obj$y

  n_pos <- sum(y_tr)
  n_neg <- sum(!y_tr)
  pw <- torch::torch_tensor(n_neg / max(1L, n_pos), dtype = torch::torch_float())$to(device = device)

  model <- siamese_change_cnn(in_ch = ch)
  model$to(device = device)
  optim <- torch::optim_adam(model$parameters, lr = cnn_cfg$lr %||% 0.001)
  bs <- as.integer(cnn_cfg$batch_size %||% 32)
  max_ep <- as.integer(cnn_cfg$epochs_max %||% 50)
  patience <- as.integer(cnn_cfg$early_stop_patience %||% 10)
  clip <- cnn_cfg$grad_clip_norm %||% 1.0

  best_f1 <- -1
  bad <- 0L
  lr <- cnn_cfg$lr %||% 0.001

  n_tr <- dim(pre_tr)[1]
  n_va <- dim(pre_va)[1]

  for (ep in seq_len(max_ep)) {
    loss_tr <- train_siamese_one_epoch(
      model, optim, seq_len(n_tr), pre_tr, post_tr, y_tr, pw, bs, clip, device
    )
    ev <- eval_siamese_loss_f1(
      model, seq_len(n_va), pre_va, post_va, y_va, pw, bs, device
    )
    m <- binary_metrics(ev$probs, y_va)
    message(sprintf(
      "Epoch %d | train loss %.4f | val loss %.4f | val F1 %.4f | IoU %.4f",
      ep, loss_tr, ev$loss, m[["F1"]], m[["IoU"]]
    ))

    if (m[["F1"]] > best_f1) {
      best_f1 <- m[["F1"]]
      bad <- 0L
      torch::torch_save(list(state_dict = model$state_dict(), in_channels = ch), checkpoint_path)
    } else {
      bad <- bad + 1L
    }

    if (bad >= patience) {
      message("Early stopping.")
      break
    }

    # manual ReduceLROnPlateau-style: halve if val loss stagnates 3 epochs
    if (bad > 0L && bad %% 3L == 0L) {
      lr <- lr * 0.5
      optim <- torch::optim_adam(model$parameters, lr = lr)
      message("Reduced LR to ", lr)
    }
  }

  invisible(model)
}


# --------------------------------------------------------------------------- #
# Load weights & full-raster sliding window (stride infer)
# --------------------------------------------------------------------------- #

load_siamese_from_checkpoint <- function(checkpoint_path, device = "cpu") {
  ck <- torch::torch_load(checkpoint_path)
  ch <- ck$in_channels
  model <- siamese_change_cnn(in_ch = ch)
  model$load_state_dict(ck$state_dict)
  model$to(device = device)
  model
}

predict_patch_prob <- function(model, pre_a, post_a, device = "cpu") {
  # pre_a, post_a: C x H x W (single patch)
  pt <- torch::torch_tensor(
    array(pre_a, dim = c(1L, dim(pre_a)[1], dim(pre_a)[2], dim(pre_a)[3])),
    dtype = torch::torch_float()
  )$to(device = device)
  po <- torch::torch_tensor(
    array(post_a, dim = c(1L, dim(post_a)[1], dim(post_a)[2], dim(post_a)[3])),
    dtype = torch::torch_float()
  )$to(device = device)
  model$eval()
  torch::with_no_grad({
    lg <- model(pt, po)
    as.numeric(torch::torch_sigmoid(lg)$cpu())
  })
}

#' Extent of raster cells (r1,c1)–(r2,c2) inclusive
ext_from_rowcol_block <- function(r, r1, c1, r2, c2) {
  cells <- c(
    terra::cellFromRowCol(r, r1, c1),
    terra::cellFromRowCol(r, r1, c2),
    terra::cellFromRowCol(r, r2, c1),
    terra::cellFromRowCol(r, r2, c2)
  )
  xy <- terra::xyFromCell(r, cells)
  terra::ext(min(xy[, 1]), max(xy[, 1]), min(xy[, 2]), max(xy[, 2]))
}

normalize_readvalues_matrix <- function(m, k, nc) {
  if (!is.matrix(m)) {
    return(NULL)
  }
  if (nrow(m) == k * nc && ncol(m) == 1L) {
    return(matrix(as.vector(m), nrow = k, ncol = nc, byrow = TRUE))
  }
  if (nrow(m) != k || ncol(m) != nc) {
    v <- as.vector(m)
    return(matrix(v, nrow = k, ncol = nc, byrow = TRUE))
  }
  m
}

read_layer_strip_matrix <- function(lay, r1, k, nc) {
  mv <- terra::readValues(lay, row = r1, nrows = k, col = 1L, ncols = nc, mat = TRUE)
  if (!is.matrix(mv)) {
    v <- as.vector(terra::readValues(lay, row = r1, nrows = k, col = 1L, ncols = nc, mat = FALSE))
    mv <- matrix(v, nrow = k, ncol = nc, byrow = TRUE)
  }
  normalize_readvalues_matrix(mv, k, nc)
}

#' Read rows r1:r2 into per-layer matrices (GeoTIFFs need readStart before readValues)
strip_to_mats_readvalues <- function(r_stack, r1, r2) {
  k <- r2 - r1 + 1L
  nc <- terra::ncol(r_stack)
  nl <- terra::nlyr(r_stack)
  lapply(seq_len(nl), function(L) {
    lay <- r_stack[[L]]
    strip_via_crop <- function() {
      ex <- ext_from_rowcol_block(r_stack, r1, 1L, r2, nc)
      cr <- terra::crop(lay, ex, snap = "near")
      refg <- terra::rast(ex, nrows = k, ncols = nc, crs = terra::crs(r_stack))
      if (terra::nrow(cr) != k || terra::ncol(cr) != nc) {
        cr <- terra::resample(cr, refg, method = "bilinear")
      }
      raster_to_matrix(cr)
    }
    opened <- tryCatch(
      {
        terra::readStart(lay)
        TRUE
      },
      error = function(e) FALSE
    )
    if (opened) {
      on.exit(try(terra::readStop(lay), silent = TRUE), add = TRUE)
    }
    tryCatch(
      read_layer_strip_matrix(lay, r1, k, nc),
      error = function(e) strip_via_crop()
    )
  })
}

sliding_window_change_prob <- function(
    feat_pre, feat_post, model, ps, stride, device = "cpu",
    infer_batch_size = 16L, progress_every = 2500L,
    row_band_steps = 8L) {
  nr <- terra::nrow(feat_post)
  nc <- terra::ncol(feat_post)
  nl <- terra::nlyr(feat_post)
  stopifnot(terra::nlyr(feat_pre) == nl)

  maps <- new.env(parent = emptyenv())
  maps$sum_m <- matrix(0, nrow = nr, ncol = nc)
  maps$cnt_m <- matrix(0L, nrow = nr, ncol = nc)

  bs <- as.integer(infer_batch_size)
  if (bs < 1L) bs <- 1L
  n_band <- max(1L, as.integer(row_band_steps))
  buf_pre <- array(NA_real_, dim = c(bs, nl, ps, ps))
  buf_post <- array(NA_real_, dim = c(bs, nl, ps, ps))
  buf_slots <- vector("list", bs)
  bi <- 0L
  n_done <- 0L

  flush_batch <- function() {
    if (bi < 1L) return()
    pre_t <- array_to_torch_nchw(buf_pre[1:bi, , , , drop = FALSE])$to(device = device)
    post_t <- array_to_torch_nchw(buf_post[1:bi, , , , drop = FALSE])$to(device = device)
    model$eval()
    probs <- torch::with_no_grad({
      lg <- model(pre_t, post_t)$squeeze(2)
      as.numeric(torch::torch_sigmoid(lg)$cpu())
    })
    for (k in seq_len(bi)) {
      p <- probs[k]
      sl <- buf_slots[[k]]
      maps$sum_m[sl$rr, sl$cc] <- maps$sum_m[sl$rr, sl$cc] + p
      maps$cnt_m[sl$rr, sl$cc] <- maps$cnt_m[sl$rr, sl$cc] + 1L
    }
    bi <<- 0L
  }

  rows <- seq(1L, nr - ps + 1L, by = stride)
  cols <- seq(1L, nc - ps + 1L, by = stride)
  n_grid <- length(rows) * length(cols)

  message(
    "Phase 4: full-raster inference (row bands ≈ ", n_band, " r0 steps, batched ", bs, ") …"
  )
  message(
    "Phase 4: inference grid ", length(rows), " × ", length(cols), " ≈ ",
    format(n_grid, big.mark = ","), " patches …"
  )

  for (ib in seq(1L, length(rows), by = n_band)) {
    i2 <- min(length(rows), ib + n_band - 1L)
    r0_band <- rows[ib:i2]
    R1 <- min(r0_band)
    R2 <- max(r0_band) + ps - 1L

    need_rows <- R2 - R1 + 1L
    message(
      "Phase 4: inference strip rows ", R1, ":", R2, " (band ", (ib - 1L) %/% n_band + 1L,
      " of ", ceiling(length(rows) / n_band), ") …"
    )
    mats_pre <- strip_to_mats_readvalues(feat_pre, R1, R2)
    mats_post <- strip_to_mats_readvalues(feat_post, R1, R2)
    nloc <- nrow(mats_pre[[1L]])
    if (nloc != need_rows) {
      stop(
        "Strip read nrow ", nloc, " != expected ", need_rows, " for rows ", R1, ":", R2,
        call. = FALSE
      )
    }

    for (r0 in r0_band) {
      r_loc <- r0 - R1 + 1L
      for (c0 in cols) {
        pa <- extract_patch_chw_from_mats(mats_pre, r_loc, c0, ps)
        pb <- extract_patch_chw_from_mats(mats_post, r_loc, c0, ps)
        if (is.null(pa) || is.null(pb)) next
        if (any(is.na(pa)) || any(is.na(pb))) next

        bi <- bi + 1L
        buf_pre[bi, , , ] <- pa
        buf_post[bi, , , ] <- pb
        buf_slots[[bi]] <- list(rr = r0:(r0 + ps - 1L), cc = c0:(c0 + ps - 1L))

        n_done <- n_done + 1L
        if (progress_every > 0L && n_done %% progress_every == 0L) {
          message(
            "Phase 4: inference progress ", format(n_done, big.mark = ","), " / ~",
            format(n_grid, big.mark = ","), " patches"
          )
        }

        if (bi >= bs) {
          flush_batch()
        }
      }
    }
    flush_batch()
    rm(mats_pre, mats_post)
  }

  sum_m <- maps$sum_m
  cnt_m <- maps$cnt_m
  prob_m <- sum_m / cnt_m
  prob_m[cnt_m == 0L] <- NA_real_
  tpl <- terra::rast(feat_post[[1]])
  terra::values(tpl) <- c(t(prob_m))
  names(tpl) <- "change_prob"
  tpl
}

postprocess_change_map <- function(prob_r, lulc_post, min_px = 9L, thr = 0.5) {
  bin <- terra::ifel(prob_r > thr, 1L, 0L)
  lp <- terra::round(lulc_post)
  bin <- terra::ifel(lp == 4L, 0L, bin)
  one <- terra::ifel(bin == 1L, 1L, NA_integer_)
  pat <- terra::patches(one)
  f <- terra::freq(pat)
  if (is.null(f) || !nrow(f)) return(bin)
  if ("layer" %in% names(f)) f <- f[f$layer == 1L, , drop = FALSE]
  small_ids <- f$value[f$count < min_px & !is.na(f$value)]
  if (length(small_ids)) {
    for (sid in small_ids) {
      bin <- terra::ifel(pat == sid, 0L, bin)
    }
  }
  bin
}

# --------------------------------------------------------------------------- #
# One transition end-to-end
# --------------------------------------------------------------------------- #

run_phase4_transition <- function(
    y_pre, y_post, cfg, device = "cpu",
    skip_train = FALSE,
    checkpoint_path = NULL) {
  cnn <- get_cnn_cfg(cfg)
  ps <- as.integer(cnn$patch_size %||% 64L)
  stride_tr <- as.integer(cnn$patch_stride_train %||% 32L)
  stride_inf <- as.integer(cnn$patch_stride_infer %||% 32L)
  infer_bs <- as.integer(cnn$infer_batch_size %||% 16L)
  infer_prog <- as.integer(cnn$infer_progress_every %||% 2500L)
  infer_row_band <- as.integer(cnn$infer_row_band_steps %||% 8L)
  centre <- as.integer(cnn$center_label_px %||% 8L)
  centre_thr <- as.numeric(cnn$center_label_threshold %||% 0.5)

  dir.create(cfg$paths$outputs_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(path_patches_dir(cfg), showWarnings = FALSE, recursive = TRUE)
  if (is.null(checkpoint_path)) {
    checkpoint_path <- path_cnn_checkpoint(cfg)
  }

  fp_pre <- path_features(y_pre, cfg)
  fp_post <- path_features(y_post, cfg)
  lp_pre <- path_lulc(y_pre, cfg)
  lp_post <- path_lulc(y_post, cfg)
  for (p in c(fp_pre, fp_post, lp_pre, lp_post)) {
    if (!file.exists(p)) stop("Missing input: ", p, call. = FALSE)
  }

  rfp <- canonicalize_features(terra::rast(fp_post))
  rpp <- rfp
  r_pre <- align_to_template(canonicalize_features(terra::rast(fp_pre)), rpp)
  l_pre <- align_lulc(terra::rast(lp_pre), rpp)
  l_post <- align_lulc(terra::rast(lp_post), rpp)

  ref <- build_change_ref(l_pre, l_post)
  terra::writeRaster(ref, path_change_ref(y_pre, y_post, cfg), overwrite = TRUE)

  nr <- terra::nrow(rpp)
  nc <- terra::ncol(rpp)

  core <- extract_patches_transition(
    r_pre, rpp, ref, ps, stride_tr, phase4_dummy_train_splits(), nr, nc,
    centre = centre, centre_thr = centre_thr
  )
  n_pos_all <- sum(core$y)
  if (n_pos_all < 1L) {
    stop(
      "No positive (change) patches in the patch grid. Your change pixels may be too sparse ",
      "for the current patch labeling rule. Try lowering `cnn.center_label_threshold` (e.g. 0.1) ",
      "and/or `cnn.center_label_px` (e.g. 4) and/or use smaller `cnn.patch_stride_train` (e.g. 16). ",
      "Also confirm change_ref has positives after alignment.",
      call. = FALSE
    )
  }

  seeds <- phase4_spatial_split_seeds(cfg)
  full <- NULL
  tr <- va <- te <- NULL
  split_seed_used <- NA_integer_
  for (s in seeds) {
    message("Phase 4: spatial train/val/test blocks, seed = ", s, " …")
    block_splits <- assign_blocks_to_splits(s)
    spl <- assign_patches_to_spatial_splits(core$r0, core$c0, nr, nc, ps, block_splits)
    if (anyNA(spl)) {
      stop("Internal error: patch centre outside block map (NA split).", call. = FALSE)
    }
    full_try <- list(pre = core$pre, post = core$post, y = core$y, split = spl)
    tr_try <- split_by_spatial_label(full_try, "train")
    va_try <- split_by_spatial_label(full_try, "val")
    te_try <- split_by_spatial_label(full_try, "test")

    if (sum(tr_try$y) < 1L) {
      message(
        "Phase 4: positive change patches exist (", n_pos_all,
        ") but none fall in train blocks for this seed; retrying …"
      )
      next
    }
    if (!dim(va_try$pre)[1]) {
      message("Phase 4: no validation patches for this seed; retrying …")
      next
    }

    full <- full_try
    tr <- tr_try
    va <- va_try
    te <- te_try
    split_seed_used <- as.integer(s)
    break
  }
  if (is.null(full)) {
    stop(
      "Could not assign spatial blocks so train has change patches and val is non-empty. ",
      "Set cnn$spatial_split_seed (or cnn$spatial_split_seeds_fallback) in study_area.yml, ",
      "or reduce patch_stride_train.",
      call. = FALSE
    )
  }
  if (split_seed_used != seeds[[1L]]) {
    message("Phase 4: using spatial split seed ", split_seed_used, " (auto-selected).")
  }

  tag <- sprintf("%d_%d", y_pre, y_post)
  pdir <- path_patches_dir(cfg)
  saveRDS(tr, file.path(pdir, paste0("train_patches_", tag, ".rds")))
  saveRDS(va, file.path(pdir, paste0("val_patches_", tag, ".rds")))
  saveRDS(te, file.path(pdir, paste0("test_patches_", tag, ".rds")))

  if (!skip_train) {
    train_siamese_model(tr, va, cnn, checkpoint_path, device = device)
  }

  model <- load_siamese_from_checkpoint(checkpoint_path, device = device)

  # Test metrics (spatial split may assign zero patches to "test" on small AOIs)
  n_te <- as.integer(dim(te$pre)[1])
  probs_te <- numeric(n_te)
  bs <- as.integer(cnn$batch_size %||% 32)
  if (n_te > 0L) {
    torch::with_no_grad({
      o <- 1L
      for (start in seq(1L, n_te, by = bs)) {
        end <- min(n_te, start + bs - 1L)
        b <- start:end
        pre_b <- array_to_torch_nchw(te$pre[b, , , , drop = FALSE])$to(device = device)
        post_b <- array_to_torch_nchw(te$post[b, , , , drop = FALSE])$to(device = device)
        lg <- model(pre_b, post_b)$squeeze(2)
        probs_te[o:(o + length(b) - 1L)] <- as.numeric(torch::torch_sigmoid(lg)$cpu())
        o <- o + length(b)
      }
    })
  } else {
    message("Phase 4: no patches in spatial \"test\" split (n_test=0); test F1/IoU are NA — full-raster output still written.")
  }
  m_test <- if (n_te > 0L) {
    binary_metrics(probs_te, te$y)
  } else {
    c(F1 = NA_real_, Precision = NA_real_, Recall = NA_real_, IoU = NA_real_, TP = 0, FP = 0, FN = 0, TN = 0)
  }

  prob_r <- tryCatch(
    sliding_window_change_prob(
      r_pre, rpp, model, ps, stride_inf, device,
      infer_batch_size = infer_bs, progress_every = infer_prog,
      row_band_steps = infer_row_band
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("bad_alloc|std::bad_alloc|cannot allocate", msg, ignore.case = TRUE)) {
        message(
          "Phase 4: inference OOM (", msg, "). Retrying batch=4, row_band_steps=4 …"
        )
        sliding_window_change_prob(
          r_pre, rpp, model, ps, stride_inf, device,
          infer_batch_size = 4L, progress_every = infer_prog, row_band_steps = 4L
        )
      } else {
        stop(e)
      }
    }
  )
  terra::writeRaster(prob_r, path_change_prob(y_pre, y_post, cfg), overwrite = TRUE)

  bin_r <- postprocess_change_map(
    prob_r, l_post,
    min_px = as.integer(cnn$min_connected_pixels %||% 9L)
  )
  names(bin_r) <- "change"
  terra::writeRaster(bin_r, path_change_bin(y_pre, y_post, cfg), overwrite = TRUE)

  data.frame(
    transition = tag,
    pre_year = y_pre,
    post_year = y_post,
    F1 = m_test[["F1"]],
    Precision = m_test[["Precision"]],
    Recall = m_test[["Recall"]],
    IoU = m_test[["IoU"]],
    n_train_patches = dim(tr$pre)[1],
    n_val_patches = dim(va$pre)[1],
    n_test_patches = n_te,
    stringsAsFactors = FALSE
  )
}

#' List Phase 4 transitions from config (index matches `transition_indices` in [run_phase4])
phase4_list_transitions <- function(cfg = load_config()) {
  cnn <- get_cnn_cfg(cfg)
  trans <- cnn$transitions
  if (!length(trans)) {
    return(data.frame(idx = integer(), pre = integer(), post = integer()))
  }
  data.frame(
    idx = seq_along(trans),
    pre = vapply(trans, function(t) as.integer(t$pre), integer(1)),
    post = vapply(trans, function(t) as.integer(t$post), integer(1)),
    stringsAsFactors = FALSE
  )
}

#' Run all transitions in config (separate checkpoint per pair: `cnn_model_best_PRE_POST.pt`)
#'
#' @param one_transition If non-NULL, integer vector `c(pre, post)` to run a single pair
#'   (ignores `transition_indices`).
#' @param transition_indices Integer vector: which entries of `config$cnn$transitions` to run
#'   (1-based). Use [phase4_list_transitions] to see indices. Example: `transition_indices = 1L`
#'   runs only the first year pair so each run stays shorter on large AOIs.
run_phase4 <- function(
    cfg = load_config(), device = "cpu", skip_train = FALSE,
    one_transition = NULL, transition_indices = NULL) {
  with_pipeline_phase("Phase 4 — CNN change detection (CPU/GPU; long runs possible)", {
  cnn <- get_cnn_cfg(cfg)
  trans <- cnn$transitions
  if (!length(trans)) stop("config$cnn$transitions is empty.", call. = FALSE)

  rows <- list()
  if (!is.null(one_transition)) {
    trans <- list(list(pre = one_transition[1], post = one_transition[2]))
  } else if (!is.null(transition_indices)) {
    ti <- as.integer(transition_indices)
    if (any(ti < 1L | ti > length(trans))) {
      stop(
        "transition_indices must be in 1:", length(trans), " (see phase4_list_transitions())",
        call. = FALSE
      )
    }
    trans <- trans[ti]
  }

  n_run <- length(trans)
  pipeline_log(
    "Phase 4: running ", n_run, " transition(s)",
    if (n_run < length(cnn$transitions)) " (subset)" else "",
    " …"
  )

  for (t in trans) {
    y1 <- as.integer(t$pre)
    y2 <- as.integer(t$post)
    pipeline_log("Phase 4: transition ", y1, " → ", y2, " (train/infer) …")
    ck <- file.path(cfg$paths$outputs_dir, sprintf("cnn_model_best_%d_%d.pt", y1, y2))
    row <- run_phase4_transition(
      y1, y2, cfg,
      device = device, skip_train = skip_train,
      checkpoint_path = ck
    )
    rows[[length(rows) + 1L]] <- row
  }

  tab_new <- do.call(rbind, rows)
  rep_path <- path_cnn_report(cfg)
  partial <- !is.null(one_transition) || !is.null(transition_indices)
  if (partial && file.exists(rep_path)) {
    old <- utils::read.csv(rep_path, stringsAsFactors = FALSE)
    if ("transition" %in% names(old) && nrow(old)) {
      keep <- !old$transition %in% tab_new$transition
      tab <- rbind(old[keep, , drop = FALSE], tab_new)
      tab <- tab[order(tab$pre_year, tab$post_year), , drop = FALSE]
    } else {
      tab <- tab_new
    }
  } else {
    tab <- tab_new
  }
  utils::write.csv(tab, rep_path, row.names = FALSE)
  pipeline_log("Phase 4: wrote ", rep_path)
  invisible(tab)
  })
}

if (interactive()) {
  message(
    "Phase 4: run_phase4(device=\"cpu\") after Phases 2–3; ",
    "for one year-pair use transition_indices=1L or one_transition=c(pre,post); ",
    "phase4_list_transitions() shows indices. torch + torch::install_torch() if needed."
  )
}

