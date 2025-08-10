################################################################################
# source("ANN_SM_CV_lag_comparison_DOWNLOAD_DEM_ELEVATION_SM_sca_FIX.R")
#
# Purpose:
#   Cross-validation comparing:
#     Configurations
#       0) SM_sca     : single-feature baseline (uses SMAP_st only)
#       1) NoLag      : base + elevation (downloaded from online DEM)
#       2) 365Lag     : base + elevation + all lag features
#       3) Composite  : base + elevation + single weighted lag
#     CV Types
#       - LOSO   : Leave-One-Station-Out
#       - Monthly: Leave-One-Month-Out (temporal CV)
#       - Random : Random K-Fold across observations
#   Metrics: R2, r, RMSE, bias, ubRMSE, KGE
#
# Notes / Fixes:
#   - Elevation is created on-the-fly by downloading from online DEM tiles
#     using the 'elevatr' package (AWS Terrain Tiles) and saved to CSV.
#   - Robust handling of NA/constant columns for the 365-lag configuration:
#       * We DO NOT drop rows requiring complete.cases over all 365 lags.
#       * We clean lags per fold: drop all-NA/zero-variance columns based
#         on TRAINING data only, and impute remaining NA with training medians.
#       * Config name bug fixed: ensure the 365Lag predictions are stored under
#         the key "365Lag" (not "Lag365"), so results aggregate correctly.
################################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest", repos = "https://cloud.r-project.org")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
  if (!requireNamespace("sp", quietly = TRUE)) install.packages("sp", repos = "https://cloud.r-project.org")
  if (!requireNamespace("elevatr", quietly = TRUE)) install.packages("elevatr", repos = "https://cloud.r-project.org")
  library(randomForest)
  library(dplyr)
  library(sp)
  library(elevatr)
})

options(timeout = max(600, getOption("timeout")))

# ========= METRICS (robust to zero variance) =========
calculate_comprehensive_metrics <- function(predictions, observations) {
  valid_idx <- which(is.finite(predictions) & is.finite(observations))
  if (length(valid_idx) < 3) {
    return(list(
      R2 = NA, RMSE = NA, bias = NA, MAE = NA, ubRMSE = NA,
      KGE = NA, r = NA, beta = NA, alpha = NA, n = length(valid_idx)
    ))
  }
  pred <- predictions[valid_idx]
  obs  <- observations[valid_idx]

  rmse <- sqrt(mean((pred - obs)^2))
  bias <- mean(pred - obs)
  mae  <- mean(abs(pred - obs))
  ubrmse <- sqrt(mean((pred - obs - bias)^2))

  # R2 robust to zero variance in obs
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- if (ss_tot == 0) NA_real_ else (1 - ss_res / ss_tot)

  # r, alpha, beta & KGE robust handling
  sdp <- sd(pred)
  sdo <- sd(obs)
  r    <- if (sdp == 0 || sdo == 0) NA_real_ else suppressWarnings(cor(pred, obs))
  beta <- if (mean(obs) == 0) NA_real_ else mean(pred) / mean(obs)
  alpha<- if (sdo == 0) NA_real_ else sdp / sdo
  kge  <- if (is.na(r) || is.na(beta) || is.na(alpha)) NA_real_
          else 1 - sqrt((r - 1)^2 + (beta - 1)^2 + (alpha - 1)^2)

  list(R2 = r2, RMSE = rmse, bias = bias, MAE = mae, ubRMSE = ubrmse,
       KGE = kge, r = r, beta = beta, alpha = alpha, n = length(valid_idx))
}

determine_winner_multi <- function(values_named, metric_name) {
  ok <- is.finite(values_named)
  if (!any(ok)) return(NA_character_)
  v <- values_named
  if (metric_name %in% c("R2","r","KGE")) {
    names(v)[which.max(v)]
  } else if (metric_name %in% c("RMSE","ubRMSE","MAE")) {
    names(v)[which.min(v)]
  } else if (metric_name == "bias") {
    names(v)[which.min(abs(v))]
  } else NA_character_
}

trySuppressWarnings <- function(expr) {
  suppressWarnings(try(expr, silent = TRUE))
}

# ========= DATE HELPERS (for Monthly CV) =========
.try_parse_dates <- function(x) {
  fmts <- c(
    "%Y-%m-%d", "%Y/%m/%d", "%Y.%m.%d",
    "%d-%m-%Y", "%d/%m/%Y", "%d.%m.%Y",
    "%m-%d-%Y", "%m/%d/%Y", "%m.%d.%Y",
    "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S",
    "%Y-%m-%dT%H:%M:%OS", "%Y-%m-%dT%H:%M:%OSZ"
  )
  out <- suppressWarnings(as.POSIXct(x, tz = "UTC", tryFormats = fmts))
  out
}

extract_month_vector <- function(df, year_fallback = NULL) {
  mn_names <- c("month","Month","MON","mon","month_num","Month_num","monthid","mth","MTH")
  for (nm in mn_names) {
    if (nm %in% names(df)) {
      m <- suppressWarnings(as.integer(df[[nm]]))
      if (any(is.finite(m))) {
        m[m < 1 | m > 12] <- NA_integer_
        return(m)
      }
    }
  }

  date_names <- c("date","Date","dates","datetime","Datetime","timestamp","time","Time","Date_UTC","date_utc")
  for (nm in date_names) {
    if (nm %in% names(df)) {
      dt <- df[[nm]]
      if (inherits(dt, "Date")) {
        return(as.integer(format(dt, "%m")))
      } else if (inherits(dt, "POSIXct") || inherits(dt, "POSIXt")) {
        return(as.integer(format(dt, "%m")))
      } else {
        dt2 <- .try_parse_dates(dt)
        if (any(!is.na(dt2))) {
          return(as.integer(format(dt2, "%m")))
        }
      }
    }
  }

  doy_names <- c("doy","DOY","doy_utc","day_of_year","dayofyear","dayno","DayNo","DayOfYear","daynum","jday","JDAY","julian","JULIAN")
  for (nm in doy_names) {
    if (nm %in% names(df)) {
      doy <- suppressWarnings(as.integer(df[[nm]]))
      if (any(is.finite(doy)) && !is.null(year_fallback)) {
        origin <- as.Date(paste0(year_fallback, "-01-01"))
        doy[doy < 1] <- 1L
        doy[doy > 366] <- 366L
        dt <- origin + (doy - 1L)
        return(as.integer(format(dt, "%m")))
      }
    }
  }

  rep(NA_integer_, nrow(df))
}

# ========= DEM ELEVATION (download + save) =========
create_elevation_from_online_dem <- function(data, z = 10, smooth_factor = 0.3, save_csv_path = NULL) {
  if (!all(c("site_name","lat","lon") %in% names(data))) {
    stop("create_elevation_from_online_dem: data must contain 'site_name', 'lat', 'lon'.")
  }

  cat("Downloading elevation from online DEM (elevatr/AWS)...\n")
  pts_df <- data.frame(x = data$lon, y = data$lat, row_id = seq_len(nrow(data)))
  sp::coordinates(pts_df) <- ~ x + y
  sp::proj4string(pts_df) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

  elev_sp <- trySuppressWarnings(elevatr::get_elev_point(pts_df,
                                                         prj = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"),
                                                         src = "aws", z = z))
  elev_first <- if (inherits(elev_sp, "try-error")) rep(NA_real_, nrow(data)) else elev_sp@data$elevation
  elev <- elev_first

  if (mean(!is.finite(elev)) > 0.25) {
    cat("High NA rate at z =", z, "- retrying with z = 8 ...\n")
    elev_sp2 <- trySuppressWarnings(elevatr::get_elev_point(pts_df,
                                                            prj = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"),
                                                            src = "aws", z = 8))
    elev2 <- if (inherits(elev_sp2, "try-error")) rep(NA_real_, nrow(data)) else elev_sp2@data$elevation
    elev[!is.finite(elev)] <- elev2[!is.finite(elev)]
  }

  elev_filled <- elev
  if (anyNA(elev_filled)) {
    na_idx <- which(!is.finite(elev_filled))
    if (length(na_idx)) {
      for (i in na_idx) {
        st <- data$site_name[i]
        st_vals <- elev_filled[data$site_name == st & is.finite(elev_filled)]
        if (length(st_vals) >= 1) elev_filled[i] <- mean(st_vals)
      }
    }
    if (any(!is.finite(elev_filled))) {
      elev_filled[!is.finite(elev_filled)] <- median(elev_filled[is.finite(elev_filled)], na.rm = TRUE)
    }
  }

  elev_smoothed <- numeric(length(elev_filled))
  for (i in seq_along(elev_filled)) {
    st <- data$site_name[i]
    st_vals <- elev_filled[data$site_name == st]
    elev_smoothed[i] <- smooth_factor * mean(st_vals, na.rm = TRUE) + (1 - smooth_factor) * elev_filled[i]
  }

  elev_norm <- elev_smoothed
  rng <- range(elev_norm, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    elev_norm[] <- 0.5
  } else {
    elev_norm <- (elev_norm - rng[1]) / (rng[2] - rng[1])
  }

  cat("Elevation created. Range (normalized):", round(min(elev_norm),3), "-", round(max(elev_norm),3), "\n")

  if (!is.null(save_csv_path)) {
    out_df <- data.frame(
      site_name = data$site_name,
      lat = data$lat,
      lon = data$lon,
      elev_first_try_m = elev_first,
      elev_filled_m = elev_filled,
      elev_smoothed_m = elev_smoothed,
      elev_normalized_0_1 = elev_norm,
      stringsAsFactors = FALSE
    )
    suppressWarnings(dir.create(dirname(save_csv_path), showWarnings = FALSE, recursive = TRUE))
    write.csv(out_df, save_csv_path, row.names = FALSE)
    cat("DEM elevation CSV saved to:", save_csv_path, "\n")
  }

  elev_norm
}

# ========= FEATURE CLEANING HELPERS (to prevent NA predictions) =========
sanitize_predictors <- function(train_df, test_df) {
  train <- as.data.frame(train_df, check.names = FALSE)
  test  <- as.data.frame(test_df,  check.names = FALSE)

  # Drop columns that are completely NA in TRAIN
  keep <- colSums(is.finite(as.matrix(train))) > 0
  if (any(!keep)) {
    train <- train[, keep, drop = FALSE]
    test  <- test[,  keep, drop = FALSE]
  }

  # Impute TRAIN medians (then use the same medians for TEST)
  med <- vapply(train, function(x) median(x[is.finite(x)], na.rm = TRUE), numeric(1))
  med[!is.finite(med)] <- 0

  for (j in seq_along(train)) {
    idx <- !is.finite(train[[j]])
    if (any(idx)) train[[j]][idx] <- med[j]
  }
  for (j in seq_along(test)) {
    idx <- !is.finite(test[[j]])
    if (any(idx)) test[[j]][idx] <- med[j]
  }

  # Drop zero-variance columns in TRAIN
  sds <- vapply(train, sd, numeric(1))
  keep2 <- is.finite(sds) & (sds > 0)
  if (any(!keep2)) {
    train <- train[, keep2, drop = FALSE]
    test  <- test[,  keep2, drop = FALSE]
  }

  # Ensure identical column order
  test <- test[, colnames(train), drop = FALSE]

  list(train = train, test = test, medians = med[keep][keep2[keep]])
}

clean_lags <- function(st_train, st_test) {
  # Ensure lag names
  if (is.null(colnames(st_train)) || any(colnames(st_train) == "")) {
    colnames(st_train) <- paste0("lag_", seq_len(ncol(st_train)))
  }
  if (is.null(colnames(st_test)) || any(colnames(st_test) == "")) {
    colnames(st_test) <- colnames(st_train)
  }

  # Align columns
  common <- intersect(colnames(st_train), colnames(st_test))
  st_train <- as.data.frame(st_train[, common, drop = FALSE], check.names = FALSE)
  st_test  <- as.data.frame(st_test[,  common, drop = FALSE], check.names = FALSE)

  # Clean using training medians and drop zero-variance columns
  sanitize_predictors(st_train, st_test)
}

# ========= MODELING HELPERS =========
rf_predict_configs <- function(train_x_base, test_x_base, y_train, st_delt_train, st_delt_test,
                               rf_params) {

  # ---- Clean lag matrix per fold to avoid NA predictions in 365Lag ----
  lag_clean <- clean_lags(st_delt_train, st_delt_test)
  st_delt_train_clean <- lag_clean$train
  st_delt_test_clean  <- lag_clean$test

  # 0) SM_sca: single feature model using only SMAP_st
  sm_col <- if ("SMAP_st" %in% colnames(train_x_base)) "SMAP_st" else colnames(train_x_base)[1]
  train_sm <- data.frame(SM_sca = train_x_base[, sm_col, drop = TRUE])
  test_sm  <- data.frame(SM_sca = test_x_base[,  sm_col, drop = TRUE])
  sm_clean <- sanitize_predictors(train_sm, test_sm)
  set.seed(123)
  model_sm <- randomForest(
    x = sm_clean$train, y = y_train,
    ntree = rf_params$ntree, nodesize = rf_params$nodesize, maxnodes = rf_params$maxnodes
  )
  pred_sm <- as.numeric(predict(model_sm, sm_clean$test))

  # 1) NoLag: base + elevation (base already finite; still sanitize for safety)
  nolag_clean <- sanitize_predictors(train_x_base, test_x_base)
  set.seed(123)
  model_nolag <- randomForest(
    x = nolag_clean$train, y = y_train,
    ntree = rf_params$ntree, nodesize = rf_params$nodesize, maxnodes = rf_params$maxnodes
  )
  pred_nolag <- as.numeric(predict(model_nolag, nolag_clean$test))

  # 2) 365Lag: base + elevation + ALL CLEANED lags
  train_365 <- cbind(nolag_clean$train, st_delt_train_clean)
  test_365  <- cbind(nolag_clean$test,  st_delt_test_clean)
  # Final sanitize in case concatenation created zero-variance cols
  step365 <- sanitize_predictors(train_365, test_365)
  set.seed(123)
  model_365 <- randomForest(
    x = step365$train, y = y_train,
    ntree = rf_params$ntree, nodesize = rf_params$nodesize, maxnodes = rf_params$maxnodes
  )
  pred_365 <- as.numeric(predict(model_365, step365$test))

  # 3) Composite: base + elevation + single weighted lag (weights from RF importance)
  # Use cleaned lags to compute importance weights
  train_for_imp <- cbind(nolag_clean$train, st_delt_train_clean)
  set.seed(123)
  rf_imp <- randomForest(
    x = train_for_imp, y = y_train,
    ntree = 30, importance = TRUE
  )
  imp <- randomForest::importance(rf_imp, type = 1)
  if (is.matrix(imp) || is.data.frame(imp)) {
    imp_vec <- imp[,1]; names(imp_vec) <- rownames(imp)
  } else {
    imp_vec <- imp
  }
  lag_names <- colnames(st_delt_train_clean)
  lag_imp <- imp_vec[lag_names]
  if (any(is.na(lag_imp)) || sum(lag_imp, na.rm = TRUE) == 0) {
    lag_imp[is.na(lag_imp)] <- 0
    if (sum(lag_imp) == 0) lag_imp <- rep(1, length(lag_imp))
  }
  weights <- as.numeric(lag_imp) / sum(lag_imp)

  comp_train <- as.matrix(st_delt_train_clean) %*% weights
  comp_test  <- as.matrix(st_delt_test_clean)  %*% weights

  train_comp <- cbind(nolag_clean$train, composite_lag = as.vector(comp_train))
  test_comp  <- cbind(nolag_clean$test,  composite_lag = as.vector(comp_test))
  comp_clean <- sanitize_predictors(train_comp, test_comp)

  set.seed(123)
  model_comp <- randomForest(
    x = comp_clean$train, y = y_train,
    ntree = rf_params$ntree, nodesize = rf_params$nodesize, maxnodes = rf_params$maxnodes
  )
  pred_comp <- as.numeric(predict(model_comp, comp_clean$test))

  # IMPORTANT: return the 365-lag predictions using the EXACT key "365Lag"
  out <- list(
    SM_sca = pred_sm,
    NoLag = pred_nolag,
    `365Lag` = pred_365,
    Composite = pred_comp
  )
  out
}

init_results_container <- function() {
  list(
    SM_sca   = list(predictions = numeric(), observations = numeric(), aux = character()),
    NoLag    = list(predictions = numeric(), observations = numeric(), aux = character()),
    `365Lag` = list(predictions = numeric(), observations = numeric(), aux = character()),
    Composite= list(predictions = numeric(), observations = numeric(), aux = character())
  )
}

accumulate_results <- function(results_list, preds_list, obs_vec, aux_vec) {
  # Guard against any misnamed keys (e.g., "Lag365")
  if (!is.null(preds_list[["Lag365"]]) && is.null(preds_list[["365Lag"]])) {
    preds_list[["365Lag"]] <- preds_list[["Lag365"]]
  }
  configs <- c("SM_sca","NoLag","365Lag","Composite")
  for (nm in configs) {
    p <- preds_list[[nm]]
    if (is.null(p)) {
      # If still missing, append NA of correct length to keep arrays aligned
      p <- rep(NA_real_, length(obs_vec))
    }
    results_list[[nm]]$predictions  <- c(results_list[[nm]]$predictions,  p)
    results_list[[nm]]$observations <- c(results_list[[nm]]$observations, obs_vec)
    results_list[[nm]]$aux          <- c(results_list[[nm]]$aux,          aux_vec)
  }
  results_list
}

finalize_and_save_results <- function(results_list, cv_label, output_dir, timestamp) {
  configurations <- c("SM_sca","NoLag","365Lag","Composite")

  overall_metrics <- data.frame()
  for (cfg in configurations) {
    m <- calculate_comprehensive_metrics(results_list[[cfg]]$predictions,
                                         results_list[[cfg]]$observations)
    overall_metrics <- rbind(overall_metrics, data.frame(
      Configuration = cfg,
      n_predictions = length(results_list[[cfg]]$predictions),
      R2 = m$R2, r = m$r, RMSE = m$RMSE, bias = m$bias, ubRMSE = m$ubRMSE, KGE = m$KGE
    ))
  }

  metrics <- c("R2","r","RMSE","bias","ubRMSE","KGE")
  winner_table <- data.frame()
  for (met in metrics) {
    vals <- setNames(
      sapply(configurations, function(cf) overall_metrics[overall_metrics$Configuration==cf, met]),
      configurations
    )
    winner_table <- rbind(
      winner_table,
      data.frame(Metric = met, t(vals), Winner = determine_winner_multi(vals, met), check.names = FALSE)
    )
  }

  suppressWarnings(dir.create(output_dir, showWarnings = FALSE, recursive = TRUE))
  write.csv(winner_table,    file.path(output_dir, paste0(cv_label, "_lag_comparison_DEM_DL_", timestamp, ".csv")), row.names = FALSE)
  write.csv(overall_metrics, file.path(output_dir, paste0(cv_label, "_overall_metrics_DEM_DL_", timestamp, ".csv")), row.names = FALSE)

  list(comparison_table = winner_table, overall_metrics = overall_metrics, pooled = results_list)
}

# ========= CV DRIVERS =========
perform_loso_cv <- function(data_all, st_delt, base_cols, rf_params,
                            max_stations = 15) {

  unique_stations <- unique(data_all$site_name)
  set.seed(456)
  if (length(unique_stations) > max_stations) {
    unique_stations <- sample(unique_stations, max_stations)
    cat("Sampling", max_stations, "stations for LOSO efficiency.\n")
  }
  n_stations <- length(unique_stations)
  results_list <- init_results_container()

  cat("Starting LOSO cross-validation...\nProgress: ")

  for (i in seq_len(n_stations)) {
    test_station <- unique_stations[i]
    test_idx  <- which(data_all$site_name == test_station)
    train_idx <- which(data_all$site_name != test_station)

    if (length(test_idx) < 3 || length(train_idx) < 50) { cat("x"); next }

    train_x_base <- data_all[train_idx, base_cols, drop = FALSE]
    test_x_base  <- data_all[test_idx,  base_cols, drop = FALSE]

    preds <- rf_predict_configs(
      train_x_base, test_x_base,
      y_train = data_all$W_obs[train_idx],
      st_delt_train = st_delt[train_idx, , drop = FALSE],
      st_delt_test  = st_delt[test_idx,  , drop = FALSE],
      rf_params = rf_params
    )

    obs <- data_all$W_obs[test_idx]
    aux <- rep(test_station, length(test_idx))
    results_list <- accumulate_results(results_list, preds, obs, aux)

    if (i %% 5 == 0) cat(i) else cat(".")
  }
  cat(" Complete!\n\n")
  results_list
}

perform_monthly_cv <- function(data_all, st_delt, base_cols, rf_params, month_vec) {
  if (length(month_vec) != nrow(data_all)) stop("Month vector length does not match data_all.")
  months_present <- sort(unique(month_vec[is.finite(month_vec)]))
  if (length(months_present) < 2) {
    cat("Monthly CV skipped: insufficient or missing month information.\n")
    return(init_results_container())
  }

  results_list <- init_results_container()
  cat("Starting Monthly cross-validation (Leave-One-Month-Out)...\nProgress: ")

  for (m in months_present) {
    test_idx  <- which(month_vec == m)
    train_idx <- which(is.finite(month_vec) & month_vec != m)

    if (length(test_idx) < 3 || length(train_idx) < 50) { cat("x"); next }

    train_x_base <- data_all[train_idx, base_cols, drop = FALSE]
    test_x_base  <- data_all[test_idx,  base_cols, drop = FALSE]

    preds <- rf_predict_configs(
      train_x_base, test_x_base,
      y_train = data_all$W_obs[train_idx],
      st_delt_train = st_delt[train_idx, , drop = FALSE],
      st_delt_test  = st_delt[test_idx,  , drop = FALSE],
      rf_params = rf_params
    )

    obs <- data_all$W_obs[test_idx]
    aux <- rep(paste0("Month_", m), length(test_idx))
    results_list <- accumulate_results(results_list, preds, obs, aux)

    cat(".")
  }
  cat(" Complete!\n\n")
  results_list
}

perform_random_cv <- function(data_all, st_delt, base_cols, rf_params, k_folds = 5, seed = 1234) {
  n <- nrow(data_all)
  if (n < k_folds * 10) {
    k_folds <- max(2, min(5, floor(n / 10)))
  }
  set.seed(seed)
  fold_id <- sample(rep(1:k_folds, length.out = n))
  results_list <- init_results_container()

  cat("Starting Random K-Fold cross-validation (k =", k_folds, ")...\nProgress: ")

  for (k in 1:k_folds) {
    test_idx  <- which(fold_id == k)
    train_idx <- which(fold_id != k)

    if (length(test_idx) < 3 || length(train_idx) < 50) { cat("x"); next }

    train_x_base <- data_all[train_idx, base_cols, drop = FALSE]
    test_x_base  <- data_all[test_idx,  base_cols, drop = FALSE]

    preds <- rf_predict_configs(
      train_x_base, test_x_base,
      y_train = data_all$W_obs[train_idx],
      st_delt_train = st_delt[train_idx, , drop = FALSE],
      st_delt_test  = st_delt[test_idx,  , drop = FALSE],
      rf_params = rf_params
    )

    obs <- data_all$W_obs[test_idx]
    aux <- rep(paste0("Fold_", k), length(test_idx))
    results_list <- accumulate_results(results_list, preds, obs, aux)

    cat(".")
  }
  cat(" Complete!\n\n")
  results_list
}

# ========= MAIN =========
lag_comparison_all_cv <- function(year = "2017",
                                  output_dir = ".",
                                  random_k_folds = 5,
                                  loso_max_stations = 15) {

  cat("\n===========================================\n")
  cat("LAG COMPARISON WITH MULTIPLE CVs (Elevation downloaded from DEM)\n")
  cat("===========================================\n")
  cat("Year:", year, "\n\n")

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  suppressWarnings(dir.create(output_dir, showWarnings = FALSE, recursive = TRUE))

  region <- "AL"
  alpha  <- 5
  data_path <- "/Users/ecomm77/Dropbox/Data/"

  # Load main data
  dff <- read.csv(paste0(data_path, "Results/", year,
                         "_OMmap_var_BM_del_0_PLR_Y_1_et2_absorption_Y_kNDVI__kMPDI__damp_N__Alaska_ORI__Alaskaall.csv"))
  dff$OM_new[] <- 0
  dff$OP_new[] <- 0
  df_original <- na.omit(dff)

  # Load lag features
  cat("Loading lag features...\n")
  load(paste0(data_path, "ANN_data/Alaska__OM_delt_", year, "_", alpha, ".Rdata"))
  st_delt <- cc
  if (nrow(st_delt) != nrow(dff)) {
    st_delt <- st_delt[1:nrow(dff), , drop = FALSE]
  }

  # Filter to valid VOD
  idx <- which(df_original$VOD < 2.3)

  # Base vars
  w_obs    <- df_original$W_obs[idx]
  et_inv   <- df_original$etr_inv[idx]
  m_clay   <- df_original$clay_inv[idx]
  OM       <- df_original$OM_[idx]
  SMAP_st  <- df_original$ST_obs[idx]
  site_nm  <- df_original$site_name[idx]
  lat_     <- df_original$lat[idx]
  lon_     <- df_original$lon[idx]
  st_delt  <- st_delt[idx, , drop = FALSE]

  # Month vector
  month_vec_full <- extract_month_vector(df_original, year_fallback = year)
  month_vec_full <- month_vec_full[idx]

  # Remove rows with missing essential fields BEFORE elevation download
  base_ok <- is.finite(w_obs) & is.finite(et_inv) & is.finite(m_clay) &
             is.finite(OM) & is.finite(SMAP_st) &
             !is.na(site_nm) & is.finite(lat_) & is.finite(lon_)
  w_obs   <- w_obs[base_ok]
  et_inv  <- et_inv[base_ok]
  m_clay  <- m_clay[base_ok]
  OM      <- OM[base_ok]
  SMAP_st <- SMAP_st[base_ok]
  site_nm <- site_nm[base_ok]
  lat_    <- lat_[base_ok]
  lon_    <- lon_[base_ok]
  st_delt <- st_delt[base_ok, , drop = FALSE]
  month_vec_full <- month_vec_full[base_ok]

  # Download DEM elevation and smooth + SAVE CSV
  dem_input_df <- data.frame(site_name = site_nm, lat = lat_, lon = lon_, stringsAsFactors = FALSE)
  dem_csv_path <- file.path(output_dir, paste0("DEM_elevation_", timestamp, ".csv"))
  elevation <- create_elevation_from_online_dem(dem_input_df, z = 10, smooth_factor = 0.3, save_csv_path = dem_csv_path)

  # Keep all rows; just require finite elevation (lags cleaned per fold)
  valid_rows <- is.finite(elevation)
  w_obs    <- w_obs[valid_rows]
  et_inv   <- et_inv[valid_rows]
  m_clay   <- m_clay[valid_rows]
  OM       <- OM[valid_rows]
  SMAP_st  <- SMAP_st[valid_rows]
  site_nm  <- site_nm[valid_rows]
  elevation <- elevation[valid_rows]
  st_delt  <- st_delt[valid_rows, , drop = FALSE]
  month_vec_full <- month_vec_full[valid_rows]

  # Base + elevation dataframe
  data_all <- data.frame(
    et_inv = et_inv,
    m_clay = m_clay,
    OM = OM,
    SMAP_st = SMAP_st,
    elevation = elevation,
    site_name = site_nm,
    W_obs = w_obs,
    check.names = FALSE
  )

  # Ensure lag column names exist
  if (is.null(colnames(st_delt)) || any(colnames(st_delt) == "")) {
    colnames(st_delt) <- paste0("lag_", seq_len(ncol(st_delt)))
  }

  cat("\nDataset summary:\n")
  cat("Total observations:", nrow(data_all), "\n")
  cat("Unique stations:", length(unique(data_all$site_name)), "\n")
  cat("Lag feature columns:", ncol(st_delt), "\n\n")

  # Features
  base_cols <- c("et_inv","m_clay","OM","SMAP_st","elevation")

  # RF params
  rf_params <- list(ntree = 40, nodesize = 15, maxnodes = 40)
  cat("Using RF params: ntree=40, nodesize=15, maxnodes=40\n\n")

  # ===== LOSO =====
  loso_results <- perform_loso_cv(
    data_all = data_all,
    st_delt = st_delt,
    base_cols = base_cols,
    rf_params = rf_params,
    max_stations = loso_max_stations
  )
  loso_out <- finalize_and_save_results(loso_results, cv_label = "LOSO", output_dir = output_dir, timestamp = timestamp)

  # ===== Monthly (Leave-One-Month-Out) =====
  monthly_results <- perform_monthly_cv(
    data_all = data_all,
    st_delt = st_delt,
    base_cols = base_cols,
    rf_params = rf_params,
    month_vec = month_vec_full
  )
  monthly_out <- finalize_and_save_results(monthly_results, cv_label = "MonthlyCV", output_dir = output_dir, timestamp = timestamp)

  # ===== Random K-Fold =====
  random_results <- perform_random_cv(
    data_all = data_all,
    st_delt = st_delt,
    base_cols = base_cols,
    rf_params = rf_params,
    k_folds = random_k_folds,
    seed = 1234
  )
  random_out <- finalize_and_save_results(random_results, cv_label = "RandomCV", output_dir = output_dir, timestamp = timestamp)

  invisible(list(
    LOSO = loso_out,
    Monthly = monthly_out,
    Random = random_out
  ))
}

# ========= RUN =========
cat("Starting Lag Comparison Analysis with LOSO, Monthly CV, and Random CV (DEM elevation downloaded; includes SM_sca baseline; NA-safe 365Lag)...\n")
all_cv_results <- lag_comparison_all_cv(
  year = "2017",
  output_dir = ".",       # change if needed
  random_k_folds = 5,     # change k for Random CV
  loso_max_stations = 15  # change LOSO station cap
)
cat("\nAll CV runs complete!\n")
