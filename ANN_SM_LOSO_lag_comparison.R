################################################################################
# source("ANN_SM_LOSO_lag_comparison.R")
#
# Purpose:
#   Leave-One-Station-Out (LOSO) cross-validation comparing:
#   1. No-lag features
#   2. All 365 lag features
#   3. Single composite lag feature
#   Comprehensive metrics: r, RMSE, bias, ubRMSE, KGE, R2
#
# Training configurations:
#   AL_only: Alaska-only training (25 sites spanning diverse soil conditions)
#   AL_US:   Alaska+US training (enhanced spatial transferability)
#
################################################################################

library(randomForest)
library(raster)
library(sp)
library(dplyr)
library(ncdf4)


# Function to calculate comprehensive performance metrics including ubRMSE
calculate_comprehensive_metrics <- function(predictions, observations) {
  valid_idx <- which(!is.na(predictions) & !is.na(observations))
  if(length(valid_idx) < 3) {
    return(list(
      R2 = NA, RMSE = NA, bias = NA, MAE = NA, ubRMSE = NA,
      KGE = NA, r = NA, beta = NA, alpha = NA,
      n = length(valid_idx)
    ))
  }
  
  pred <- predictions[valid_idx]
  obs <- observations[valid_idx]
  
  # Basic metrics
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- 1 - (ss_res / ss_tot)
  rmse <- sqrt(mean((pred - obs)^2))
  bias <- mean(pred - obs)
  mae <- mean(abs(pred - obs))
  
  # Unbiased RMSE (ubRMSE) - RMSE after removing bias
  ubrmse <- sqrt(mean((pred - obs - bias)^2))
  
  # Pearson correlation
  r <- cor(pred, obs)
  
  # KGE components
  # Beta - bias ratio
  beta <- mean(pred) / mean(obs)
  
  # Alpha - variability ratio
  alpha <- sd(pred) / sd(obs)
  
  # KGE calculation
  kge <- 1 - sqrt((r - 1)^2 + (beta - 1)^2 + (alpha - 1)^2)
  
  return(list(
    R2 = r2, 
    RMSE = rmse, 
    bias = bias, 
    MAE = mae,
    ubRMSE = ubrmse,
    KGE = kge, 
    r = r, 
    beta = beta, 
    alpha = alpha,
    n = length(valid_idx)
  ))
}

# Function to load real DEM elevation data
load_real_dem <- function(data, csv) {
  if (!file.exists(csv) || !("site_name"%in%names(data))) return(rep(0.5,nrow(data)))
  dem <- tryCatch(read.csv(csv), error=function(e) NULL)
  if (is.null(dem) || !all(c("site_name","elevation")%in%names(dem))) return(rep(0.5,nrow(data)))
  e <- numeric(nrow(data))
  for (i in seq_len(nrow(data))){
    k <- which(dem$site_name==data$site_name[i])[1]
    e[i] <- ifelse(is.finite(k), dem$elevation[k], 500)
  }
  mn <- min(e); mx <- max(e); if (mx>mn) (e-mn)/(mx-mn) else rep(0.5,length(e))
}

# Function to load bulk density at 0cm from NetCDF
load_bd_0cm <- function(data, nc_path){
  if (!("lat"%in%names(data) && "lon"%in%names(data))) return(rep(1.3,nrow(data)))
  if (!file.exists(nc_path)) return(rep(1.3,nrow(data)))
  nc <- nc_open(nc_path); on.exit(try(nc_close(nc),silent=TRUE))
  lat <- ncvar_get(nc,"latitude"); lon <- ncvar_get(nc,"longitude")
  vname <- intersect(c("BD_0cm","bd_0cm","bulk_density_0cm","BD"), names(nc$var))
  if (!length(vname)) return(rep(1.3,nrow(data)))
  bd <- ncvar_get(nc, vname[1]); out <- numeric(nrow(data))
  for (i in seq_len(nrow(data))){
    ii <- which.min(abs(lat - data$lat[i])); jj <- which.min(abs(lon - data$lon[i]))
    val <- suppressWarnings(bd[jj,ii]); out[i] <- ifelse(is.finite(val), val, 1.3)
  }
  out
}

# Function to add real elevation and bulk density features
add_real_features <- function(data) {
  # Load real DEM elevation data
  csv_path <- "/Users/ecomm77/Dropbox/Data/Alaska_real_DEM_features_20250810_112821.csv"
  data$elevation <- load_real_dem(data, csv_path)
  
  # Load bulk density at 0cm
  nc_path <- "/Users/ecomm77/Dropbox/Data/OM_world/Chang/Hydraul_Param_SoilGrids_Schaap_sl1_Py.nc"
  data$bd_0cm <- load_bd_0cm(data, nc_path)
  
  return(data)
}

# Function to determine winner for each metric
determine_winner <- function(nolag, lag365, composite, metric_name) {
  values <- c(nolag, lag365, composite)
  names(values) <- c("NoLag", "365Lag", "Composite")
  
  # For these metrics, higher is better
  if(metric_name %in% c("R2", "r", "KGE")) {
    winner <- names(values)[which.max(values)]
  } 
  # For these metrics, lower is better (closer to 0 for bias)
  else if(metric_name %in% c("RMSE", "ubRMSE", "MAE")) {
    winner <- names(values)[which.min(values)]
  }
  # For bias, closer to 0 is better
  else if(metric_name == "bias") {
    winner <- names(values)[which.min(abs(values))]
  }
  else {
    winner <- "Unknown"
  }
  
  return(winner)
}

# Function to normalize dataframe columns for multi-region support
normalize_df_columns <- function(df, region_tag) {
  map1 <- function(df, target, cands, default = NA) {
    if (!(target %in% names(df))) {
      hit <- cands[cands %in% names(df)]
      if (length(hit) > 0) df[[target]] <- df[[hit[1]]]
    }
    if (!(target %in% names(df))) df[[target]] <- default
    df
  }
  num <- function(x) suppressWarnings(as.numeric(x))
  
  df <- map1(df, "W_obs", c("W_obs", "W_insitu", "W", "SM_insitu"))
  df <- map1(df, "OM_", c("OM_", "OM", "OM_pct"))
  df <- map1(df, "clay_inv", c("clay_inv", "clay", "clay_frac", "CLAY", "Clay"))
  df <- map1(df, "etr_inv", c("etr_inv", "etr", "ETR"))
  df <- map1(df, "ST_obs", c("ST_obs", "ST", "temp", "Ts"))
  df <- map1(df, "VOD", c("VOD", "tau", "VOD_L2"), 0.5)
  df <- map1(df, "W_sca", c("W_sca", "W_smap", "SMAP_sca"), 0.25)
  df <- map1(df, "site_name", c("site_name", "site", "station", "station_name"))
  df <- map1(df, "lat", c("lat", "latitude", "Latitude"), ifelse(region_tag == "AL", 64, 39))
  df <- map1(df, "lon", c("lon", "longitude", "Longitude"), ifelse(region_tag == "AL", -149, -98))
  
  if (all(is.na(df$site_name)) && "id_site" %in% names(df)) {
    df$site_name <- paste0(region_tag, "_Site_", df$id_site)
  }
  miss <- which(is.na(df$site_name) | df$site_name == "")
  if (length(miss) > 0) {
    df$site_name[miss] <- paste0(region_tag, "_", round(num(df$lat[miss]), 3), "_", round(num(df$lon[miss]), 3))
  }
  
  for (nm in c("W_obs", "OM_", "clay_inv", "etr_inv", "ST_obs", "VOD", "W_sca", "lat", "lon")) {
    if (nm %in% names(df)) df[[nm]] <- num(df[[nm]])
  }
  df$region <- region_tag
  df
}

# Main LOSO lag comparison function
loso_lag_comparison <- function(year = "2017", training_mode = "AL_only") {
  
  cat("\n===========================================\n")
  cat("LOSO LAG COMPARISON ANALYSIS\n")
  cat("===========================================\n")
  cat("Comparing: No-Lag vs 365-Lag vs Composite-Lag\n")
  cat("Year:", year, "\n")
  cat("Training mode:", training_mode, "\n\n")
  
  if (training_mode == "AL_only") {
    cat("=== ALASKA-ONLY TRAINING CONFIGURATION ===\n")
    cat("Training on Alaska sites spanning diverse soil conditions\n")
    cat("25 sites: lowland tundra to upland mineral soils\n\n")
    regions <- c("AL")
  } else if (training_mode == "AL_US") {
    cat("=== ALASKA+US TRAINING CONFIGURATION ===\n") 
    cat("Training on Alaska + USCRN/SCAN for enhanced spatial transferability\n")
    cat("Combining permafrost-affected and temperate mineral soils\n\n")
    regions <- c("AL", "US")
  } else {
    stop("Invalid training_mode. Use 'AL_only' or 'AL_US'")
  }
  
  # Setup
  alpha <- 5
  data_path <- "/Users/ecomm77/Dropbox/Data/"
  
  # Load data for each region
  cat("=== LOADING DATA ===\n")
  all_data <- list()
  all_lags <- list()
  
  for (reg in regions) {
    cat(sprintf("Loading %s region...\n", reg))
    if (reg == "US") {
      df_raw <- read.csv(paste0(data_path, "Results/", year, "_OMmap_var_BM_del_100_PLR_Y_1_et_absorption_Y___USCRN_ORI_USall.csv"))
      load(paste0(data_path, "ANN_data/US__OM_delt_", year, "_", alpha, ".Rdata"))
      cat("  Loaded USCRN/SCAN stations\n")
    } else {
      df_raw <- read.csv(paste0(data_path, "Results/", year, "_OMmap_var_BM_del_0_PLR_Y_1_et2_absorption_Y_kNDVI__kMPDI__damp_N__Alaska_ORI__Alaskaall.csv"))
      load(paste0(data_path, "ANN_data/Alaska__OM_delt_", year, "_", alpha, ".Rdata"))
      cat("  Loaded Alaska sites\n")
    }
    
    df_std <- normalize_df_columns(df_raw, reg)
    keep <- c("W_obs", "OM_", "clay_inv", "etr_inv", "ST_obs", "VOD", "W_sca", "site_name", "region", "lat", "lon")
    for (c in keep) if (!(c %in% names(df_std))) df_std[[c]] <- NA
    df_keep <- df_std[, keep]
    
    n <- min(nrow(df_keep), nrow(cc))
    all_data[[reg]] <- df_keep[seq_len(n), ]
    all_lags[[reg]] <- cc[seq_len(n), ]
    
    cat(sprintf("  %s: %d observations, %d unique sites\n", 
                reg, nrow(all_data[[reg]]), length(unique(all_data[[reg]]$site_name))))
  }
  
  # Combine data
  dfc <- do.call(rbind, all_data)
  st_delt <- do.call(rbind, all_lags)
  dfc$site_uid <- paste(dfc$region, dfc$site_name, sep = "_")
  
  cat(sprintf("\nCombined data: %d observations from %d regions\n", 
              nrow(dfc), length(regions)))
  print(table(dfc$region))
  
  # Apply VOD filter and process data
  idx <- which(is.finite(dfc$VOD) & dfc$VOD < 2.3)
  dfc <- dfc[idx, , drop = FALSE]
  st_delt <- st_delt[idx, , drop = FALSE]
  
  # Remove NA values
  valid_rows <- complete.cases(st_delt) & 
                with(dfc, is.finite(W_obs) & is.finite(etr_inv) & is.finite(clay_inv) & 
                     is.finite(OM_) & is.finite(ST_obs) & is.finite(lat) & is.finite(lon))
  
  dfc <- dfc[valid_rows, , drop = FALSE]
  st_delt <- st_delt[valid_rows, , drop = FALSE]
  
  # Extract variables for compatibility with existing code
  w_obs <- dfc$W_obs
  et_inv <- dfc$etr_inv
  m_clay <- dfc$clay_inv
  OM <- dfc$OM_
  SMAP_st <- dfc$ST_obs
  lat_ <- dfc$lat
  lon_ <- dfc$lon
  site_name_ <- dfc$site_uid  # Use unique site identifier
  
  # Create data frame with coordinates
  data_with_coords <- data.frame(
    et_inv = et_inv,
    m_clay = m_clay,
    OM = OM,
    SMAP_st = SMAP_st,
    lat = lat_,
    lon = lon_,
    site_name = site_name_,
    W_obs = w_obs
  )
  
  # Add real elevation and bulk density features
  cat("Adding real elevation and bulk density features...\n")
  enhanced_data <- add_real_features(data_with_coords)
  
  # Use ALL 365 lag features (no sampling)
  st_delt_sampled <- st_delt  # Use all lag features, not sampled
  
  # Get unique stations
  unique_stations <- unique(enhanced_data$site_name)
  
  # Use ALL stations for proper LOSO spatial validation
  n_stations <- length(unique_stations)
  
  cat("\nDataset summary:\n")
  cat("Total observations:", nrow(enhanced_data), "\n")
  cat("Unique stations:", n_stations, "\n")
  cat("Training regions:", paste(regions, collapse = ", "), "\n")
  if (training_mode == "AL_US") {
    region_counts <- table(dfc$region)
    for (reg in names(region_counts)) {
      sites_in_reg <- length(unique(dfc$site_uid[dfc$region == reg]))
      cat(sprintf("  %s: %d observations, %d sites\n", reg, region_counts[reg], sites_in_reg))
    }
  }
  cat("Full lag features:", ncol(st_delt), "\n")
  cat("Lag features (for 365-lag):", ncol(st_delt_sampled), "\n")
  cat("Using ALL", n_stations, "stations for proper LOSO spatial validation\n\n")
  
  # Define features
  base_features_cols <- c("et_inv", "m_clay", "OM", "SMAP_st")
#  additional_features_cols <- c("elevation", "bd_0cm")  # Real elevation and bulk density features
  additional_features_cols <- c("elevation")  # Real elevation and bulk density features  
  enhanced_features_cols <- c(base_features_cols, additional_features_cols)
  
  # Set parameters optimized for spatial generalization
  rf_params <- list(ntree = 100, nodesize = 5, maxnodes = NULL)
  cat("Using parameters optimized for spatial generalization: ntree=100, nodesize=5, maxnodes=NULL\n\n")
  
  # Initialize results storage for three configurations
  configurations <- c("NoLag", "365Lag", "Composite")
  results_list <- list()
  
  for(config in configurations) {
    results_list[[config]] <- list(
      predictions = numeric(),
      observations = numeric(),
      stations = character()
    )
  }
  
  # Progress indicator
  cat("Starting LOSO cross-validation for all three configurations...\n")
  cat("Progress: ")
  
  # Perform LOSO cross-validation
  for(i in 1:n_stations) {
    test_station <- unique_stations[i]
    
    # Create train/test splits
    test_idx <- which(enhanced_data$site_name == test_station)
    train_idx <- which(enhanced_data$site_name != test_station)
    
    if(length(test_idx) < 3 || length(train_idx) < 50) {
      cat("x")
      next
    }
    
    ## Configuration 1: No-Lag (Enhanced features only)
    train_features_nolag <- enhanced_data[train_idx, enhanced_features_cols]
    test_features_nolag <- enhanced_data[test_idx, enhanced_features_cols]
    
    set.seed(123)
    if(is.null(rf_params$maxnodes)) {
      model_nolag <- randomForest(
        x = train_features_nolag, 
        y = enhanced_data$W_obs[train_idx],
        ntree = rf_params$ntree,
        nodesize = rf_params$nodesize
      )
    } else {
      model_nolag <- randomForest(
        x = train_features_nolag, 
        y = enhanced_data$W_obs[train_idx],
        ntree = rf_params$ntree,
        nodesize = rf_params$nodesize,
        maxnodes = rf_params$maxnodes
      )
    }
    
    pred_nolag <- predict(model_nolag, test_features_nolag)
    
    ## Configuration 2: 365-Lag (Enhanced + all sampled lag features)
    train_features_365 <- cbind(enhanced_data[train_idx, enhanced_features_cols],
                               st_delt_sampled[train_idx, ])
    test_features_365 <- cbind(enhanced_data[test_idx, enhanced_features_cols],
                              st_delt_sampled[test_idx, ])
    
    set.seed(123)
    if(is.null(rf_params$maxnodes)) {
      model_365 <- randomForest(
        x = train_features_365, 
        y = enhanced_data$W_obs[train_idx],
        ntree = rf_params$ntree,
        nodesize = rf_params$nodesize
      )
    } else {
      model_365 <- randomForest(
        x = train_features_365, 
        y = enhanced_data$W_obs[train_idx],
        ntree = rf_params$ntree,
        nodesize = rf_params$nodesize,
        maxnodes = rf_params$maxnodes
      )
    }
    
    pred_365 <- predict(model_365, test_features_365)
    
    ## Configuration 3: Composite Lag (Enhanced + single composite lag)
    # Calculate composite lag using training data
    train_features_for_importance <- cbind(
      enhanced_data[train_idx, base_features_cols], 
      st_delt_sampled[train_idx, ]
    )
    
    set.seed(123)
    rf_importance <- randomForest(
      x = train_features_for_importance, 
      y = enhanced_data$W_obs[train_idx],
      ntree = 30, 
      importance = TRUE
    )
    
    importance_values <- importance(rf_importance, type = 1)
    lag_start_col <- length(base_features_cols) + 1
    lag_importance <- importance_values[lag_start_col:nrow(importance_values), 1]
    weights <- lag_importance / sum(lag_importance)
    
    # Create composite lag
    composite_lag_train <- as.matrix(st_delt_sampled[train_idx, ]) %*% weights
    composite_lag_test <- as.matrix(st_delt_sampled[test_idx, ]) %*% weights
    
    train_features_composite <- cbind(
      enhanced_data[train_idx, enhanced_features_cols],
      composite_lag = as.vector(composite_lag_train)
    )
    
    test_features_composite <- cbind(
      enhanced_data[test_idx, enhanced_features_cols],
      composite_lag = as.vector(composite_lag_test)
    )
    
    set.seed(123)
    if(is.null(rf_params$maxnodes)) {
      model_composite <- randomForest(
        x = train_features_composite, 
        y = enhanced_data$W_obs[train_idx],
        ntree = rf_params$ntree,
        nodesize = rf_params$nodesize
      )
    } else {
      model_composite <- randomForest(
        x = train_features_composite, 
        y = enhanced_data$W_obs[train_idx],
        ntree = rf_params$ntree,
        nodesize = rf_params$nodesize,
        maxnodes = rf_params$maxnodes
      )
    }
    
    pred_composite <- predict(model_composite, test_features_composite)
    
    # Store results
    obs <- enhanced_data$W_obs[test_idx]
    station_names <- rep(test_station, length(test_idx))
    
    results_list[["NoLag"]]$predictions <- c(results_list[["NoLag"]]$predictions, pred_nolag)
    results_list[["NoLag"]]$observations <- c(results_list[["NoLag"]]$observations, obs)
    results_list[["NoLag"]]$stations <- c(results_list[["NoLag"]]$stations, station_names)
    
    results_list[["365Lag"]]$predictions <- c(results_list[["365Lag"]]$predictions, pred_365)
    results_list[["365Lag"]]$observations <- c(results_list[["365Lag"]]$observations, obs)
    results_list[["365Lag"]]$stations <- c(results_list[["365Lag"]]$stations, station_names)
    
    results_list[["Composite"]]$predictions <- c(results_list[["Composite"]]$predictions, pred_composite)
    results_list[["Composite"]]$observations <- c(results_list[["Composite"]]$observations, obs)
    results_list[["Composite"]]$stations <- c(results_list[["Composite"]]$stations, station_names)
    
    # Progress indicator
    if(i %% 5 == 0) cat(i) else cat(".")
  }
  
  cat(" Complete!\n\n")
  
  # Calculate overall metrics for each configuration
  overall_metrics <- data.frame()
  
  for(config in configurations) {
    metrics <- calculate_comprehensive_metrics(
      results_list[[config]]$predictions,
      results_list[[config]]$observations
    )
    
    overall_row <- data.frame(
      Configuration = config,
      n_predictions = length(results_list[[config]]$predictions),
      R2 = metrics$R2,
      r = metrics$r,
      RMSE = metrics$RMSE,
      bias = metrics$bias,
      ubRMSE = metrics$ubRMSE,
      KGE = metrics$KGE
    )
    
    overall_metrics <- rbind(overall_metrics, overall_row)
  }
  
  # Create winner table
  winner_table <- data.frame(
    Metric = c("R2", "r", "RMSE", "bias", "ubRMSE", "KGE"),
    NoLag = c(overall_metrics$R2[1], overall_metrics$r[1], overall_metrics$RMSE[1], 
              overall_metrics$bias[1], overall_metrics$ubRMSE[1], overall_metrics$KGE[1]),
    Lag365 = c(overall_metrics$R2[2], overall_metrics$r[2], overall_metrics$RMSE[2], 
               overall_metrics$bias[2], overall_metrics$ubRMSE[2], overall_metrics$KGE[2]),
    Composite = c(overall_metrics$R2[3], overall_metrics$r[3], overall_metrics$RMSE[3], 
                  overall_metrics$bias[3], overall_metrics$ubRMSE[3], overall_metrics$KGE[3]),
    Winner = character(6)
  )
  
  # Determine winners
  for(i in 1:nrow(winner_table)) {
    metric_name <- winner_table$Metric[i]
    winner_table$Winner[i] <- determine_winner(
      winner_table$NoLag[i], 
      winner_table$Lag365[i], 
      winner_table$Composite[i], 
      metric_name
    )
  }
  
  # Display results
  cat("===========================================\n")
  cat("LOSO LAG COMPARISON RESULTS\n")
  cat("===========================================\n\n")
  
  cat("Overall Performance Comparison:\n")
  cat("-------------------------------\n")
  print(overall_metrics, row.names = FALSE, digits = 4)
  
  cat("\n\nDetailed Comparison Table:\n")
  cat("==========================\n")
  print(winner_table, row.names = FALSE, digits = 4)
  
  # Count wins
  win_counts <- table(winner_table$Winner)
  cat("\n\nWin Summary:\n")
  cat("============\n")
  for(config in names(win_counts)) {
    cat(sprintf("%s: %d wins\n", config, win_counts[config]))
  }
  
  # Determine overall best
  best_config <- names(win_counts)[which.max(win_counts)]
  cat(sprintf("\nOverall Best Configuration: %s\n", best_config))
  
  # Key insights
  cat("\n===========================================\n")
  cat("KEY INSIGHTS\n")
  cat("===========================================\n")
  
  # R2 comparison
  r2_values <- winner_table$NoLag[1:3]  # R2, r are first two
  names(r2_values) <- c("NoLag", "365Lag", "Composite")
  best_r2 <- names(r2_values)[which.max(r2_values)]
  
  cat(sprintf("Best R² performance: %s (R² = %.3f)\n", best_r2, max(r2_values)))
  
  # RMSE comparison
  rmse_values <- c(winner_table$NoLag[3], winner_table$Lag365[3], winner_table$Composite[3])
  names(rmse_values) <- c("NoLag", "365Lag", "Composite")
  best_rmse <- names(rmse_values)[which.min(rmse_values)]
  
  cat(sprintf("Best RMSE performance: %s (RMSE = %.3f)\n", best_rmse, min(rmse_values)))
  
  # KGE comparison
  kge_values <- c(winner_table$NoLag[6], winner_table$Lag365[6], winner_table$Composite[6])
  names(kge_values) <- c("NoLag", "365Lag", "Composite")
  best_kge <- names(kge_values)[which.max(kge_values)]
  
  cat(sprintf("Best KGE performance: %s (KGE = %.3f)\n", best_kge, max(kge_values)))
  
  if(best_r2 == best_rmse && best_rmse == best_kge) {
    cat(sprintf("\n✓ Clear winner across all major metrics: %s\n", best_r2))
  } else {
    cat("\n± Mixed results - different configurations excel in different metrics\n")
  }
  
  # Lag feature effectiveness
  nolag_r2 <- winner_table$NoLag[1]
  composite_r2 <- winner_table$Composite[1]
  lag365_r2 <- winner_table$Lag365[1]
  
  cat(sprintf("\nLag Feature Impact on R²:\n"))
  cat(sprintf("  No-Lag → Composite: %+.3f improvement\n", composite_r2 - nolag_r2))
  cat(sprintf("  No-Lag → 365-Lag: %+.3f improvement\n", lag365_r2 - nolag_r2))
  cat(sprintf("  Composite → 365-Lag: %+.3f improvement\n", lag365_r2 - composite_r2))
  
  # Save results
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  mode_suffix <- if(training_mode == "AL_only") "ALonly" else "AL_US"
  write.csv(winner_table, paste0("LOSO_lag_comparison_", mode_suffix, "_", timestamp, ".csv"), row.names = FALSE)
  write.csv(overall_metrics, paste0("LOSO_overall_metrics_", mode_suffix, "_", timestamp, ".csv"), row.names = FALSE)
  
  cat(sprintf("\nResults saved to:\n"))
  cat(sprintf("  - LOSO_lag_comparison_%s_%s.csv\n", mode_suffix, timestamp))
  cat(sprintf("  - LOSO_overall_metrics_%s_%s.csv\n", mode_suffix, timestamp))
  
  return(list(
    comparison_table = winner_table,
    overall_metrics = overall_metrics,
    detailed_results = results_list
  ))
}

# ===========================================
# RUN LOSO LAG COMPARISON
# ===========================================

# Configuration: Choose training mode
# "AL_only" = Alaska-only training (25 sites spanning diverse soil conditions)
# "AL_US"   = Alaska+US training (enhanced spatial transferability)
training_mode <- "AL_US"  # Change this to "AL_only" for Alaska-only training

cat("Starting LOSO Lag Comparison Analysis...\n")
comparison_results <- loso_lag_comparison(year = "2017", training_mode = training_mode)
cat("\nLOSO Lag Comparison complete!\n")