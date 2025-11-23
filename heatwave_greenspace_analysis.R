################################################################################
# Greenspace Landscape Patterns and Heatwave-related Mortality Analysis
# 
# This script implements the two-stage analytical framework described in:
# "The Overlooked Role of Green Space Landscape Patterns in Regulating 
# Heatwave-related Mortality Risk"
#
# Stage 1: City-specific associations using distributed lag non-linear models
# Stage 2: Meta-analysis to pool estimates and assess effect modification
#
# Author: Cheng He
# Date: 11/23/2025
################################################################################

# Load required packages -------------------------------------------------------
library(mgcv)          # Generalized Additive Models
library(dlnm)          # Distributed Lag Non-linear Models
library(mvmeta)        # Multivariate Meta-analysis
library(splines)       # Spline functions
library(dplyr)         # Data manipulation

################################################################################
# STAGE 1: CITY-SPECIFIC ANALYSIS FUNCTION
################################################################################

analyze_heatwave_mortality <- function(data_input, 
                                       lag_days = 10,
                                       df_time = 4,
                                       df_humidity = 4) {
  #' City-specific heatwave-mortality association using DLNM
  #' 
  #' @param data_input Data frame with columns: code_city, death, heatwave, 
  #'                   dow (day of week), rh (relative humidity), time
  #' @param lag_days Maximum lag days (default: 10)
  #' @param df_time Degrees of freedom per year for time trend (default: 4)
  #' @param df_humidity Degrees of freedom for humidity spline (default: 4)
  #' @return Pooled relative risk with 95% CI for heatwave exposure
  
  # Get list of cities
  regions <- as.character(unique(data_input$code_city))
  
  # Split data by city
  data_list <- lapply(regions, function(x) data_input[data_input$code_city == x, ])
  names(data_list) <- regions
  
  # Calculate temperature range bounds for crossbasis
  ranges <- t(sapply(data_list, function(x) range(x$heatwave, na.rm = TRUE)))
  bound <- colMeans(ranges)
  
  # Initialize storage for coefficients and covariance matrices
  n_cities <- length(data_list)
  coef_all <- matrix(NA, n_cities, 4, 
                     dimnames = list(regions, paste0("b", seq(4))))
  vcov_all <- vector("list", n_cities)
  names(vcov_all) <- regions
  
  # Define crossbasis parameters
  # Variable dimension: B-spline with degree 2 and 4 df
  varknots <- equalknots(bound, fun = "bs", degree = 2, df = 4)
  
  # Lag dimension: Natural spline with 4 df (3 internal knots in log scale)
  lagknots <- logknots(lag_days, df = 3, int = TRUE)
  
  argvar <- list(fun = "bs", degree = 2, knots = varknots, bound = bound)
  arglag <- list(fun = "ns", knots = lagknots)
  
  # Suppress warnings during model fitting
  options(warn = -1)
  
  # Loop through each city
  for (i in seq_along(data_list)) {
    
    city_data <- data_list[[i]]
    
    # Create crossbasis for heatwave exposure
    cb_heatwave <- crossbasis(city_data$heatwave, 
                              lag = lag_days,
                              argvar = argvar, 
                              arglag = arglag)
    
    # Fit GAM with quasi-Poisson distribution
    # Model specification:
    # - Heatwave exposure with distributed lag
    # - Day of week (factor)
    # - Relative humidity (natural spline with df_humidity df)
    # - Long-term trend and seasonality (natural spline with df_time df per year)
    model_gam <- gam(death ~ cb_heatwave + 
                       as.factor(dow) + 
                       ns(rh, df = df_humidity) + 
                       ns(time, df = df_time * 3),
                     family = quasipoisson(), 
                     data = city_data)
    
    # Extract cross-reduced estimates (cumulative effect over lag)
    cr_all <- crossreduce(cb_heatwave, model_gam)
    
    # Store coefficients and variance-covariance matrix
    coef_all[i, ] <- coef(cr_all)
    vcov_all[[i]] <- vcov(cr_all)
  }
  
  # Restore warning options
  options(warn = 0)
  
  ################################################################################
  # STAGE 2: META-ANALYSIS
  ################################################################################
  
  # Multivariate random-effects meta-analysis using REML
  meta_model <- mvmeta(coef_all ~ 1, vcov_all, method = "reml")
  
  # Predict overall cumulative relative risk
  xvar <- seq(bound[1], bound[2], by = 1)
  bvar <- do.call("onebasis", c(list(x = xvar), attr(cb_heatwave, "argvar")))
  
  pred_overall <- crosspred(bvar, 
                            coef = coef(meta_model), 
                            vcov = vcov(meta_model),
                            model.link = "log", 
                            by = 0.1,
                            from = bound[1], 
                            to = bound[2],
                            cen = 0,      # Reference value: no heatwave
                            cumul = TRUE) # Cumulative effect
  
  # Extract RR and 95% CI for heatwave vs. no heatwave (at x = 1)
  results <- with(pred_overall, 
                  c(RR = allRRfit["1"], 
                    Lower = allRRlow["1"], 
                    Upper = allRRhigh["1"]))
  
  return(round(results, 3))
}

################################################################################
# FUNCTION TO CALCULATE P-VALUE FOR DIFFERENCE BETWEEN TWO RRs
################################################################################

calculate_pvalue_diff <- function(RR1, RR1_lower, RR1_upper,
                                   RR2, RR2_lower, RR2_upper) {
  #' Calculate p-value for difference between two relative risks
  #' 
  #' @param RR1 Relative risk for group 1
  #' @param RR1_lower Lower 95% CI for group 1
  #' @param RR1_upper Upper 95% CI for group 1
  #' @param RR2 Relative risk for group 2
  #' @param RR2_lower Lower 95% CI for group 2
  #' @param RR2_upper Upper 95% CI for group 2
  #' @return P-value from two-sided Z-test
  
  # Calculate log(RR) and standard errors
  logRR1 <- log(RR1)
  SE1 <- (log(RR1_upper) - log(RR1_lower)) / (2 * 1.96)
  
  logRR2 <- log(RR2)
  SE2 <- (log(RR2_upper) - log(RR2_lower)) / (2 * 1.96)
  
  # Calculate difference and its standard error
  diff_est <- logRR1 - logRR2
  diff_se <- sqrt(SE1^2 + SE2^2)
  
  # Calculate Z-score and two-sided p-value
  z_score <- diff_est / diff_se
  pvalue <- 2 * pnorm(-abs(z_score))
  
  return(pvalue)
}

################################################################################
# MAIN ANALYSIS WORKFLOW
################################################################################

# Load data --------------------------------------------------------------------
# Mortality and meteorological data
data_main <- read.csv("ccches_heatwave.csv", header = TRUE)

# Select warm season (May-October)
data_main <- subset(data_main, season == "warm")

# Rename columns to match expected format
colnames(data_main)[colnames(data_main) == "total"] <- "death"
colnames(data_main)[colnames(data_main) == "temp"] <- "tmean"

# Load city-level metadata with landscape metrics
meta_data <- read.csv("city_metadata_greenspace.csv", header = TRUE)

# Define landscape metrics to analyze
landscape_metrics <- c("PLAND",   # Percentage of Landscape
                       "PD",      # Patch Density
                       "LPI",     # Largest Patch Index
                       "LSI",     # Landscape Shape Index
                       "SPLIT")   # Splitting Index

# Define heatwave definitions --------------------------------------------------
# Four definitions based on different intensity thresholds and durations
heatwave_definitions <- list(
  "92.5%_2day" = "hw_92.5_2d",
  "92.5%_3day" = "hw_92.5_3d",
  "95%_2day" = "hw_95_2d",
  "95%_3day" = "hw_95_3d"
)

# Initialize results storage
results_all <- list()

# Loop through each heatwave definition ----------------------------------------
for (hw_name in names(heatwave_definitions)) {
  
  cat("\n========================================\n")
  cat("Analyzing:", hw_name, "\n")
  cat("========================================\n")
  
  hw_column <- heatwave_definitions[[hw_name]]
  
  # Prepare data for this heatwave definition
  data_hw <- data_main
  colnames(data_hw)[colnames(data_hw) == hw_column] <- "heatwave"
  
  # Filter cities with at least 3 heatwave days
  data_hw <- data_hw %>%
    group_by(code_city) %>%
    filter(sum(heatwave, na.rm = TRUE) > 3) %>%
    ungroup()
  
  # Initialize results data frame for this heatwave definition
  results_df <- data.frame(
    Heatwave_Definition = character(),
    Landscape_Metric = character(),
    Category = character(),
    RR = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    P_value = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each landscape metric -----------------------------------------
  for (metric in landscape_metrics) {
    
    cat("\n  Analyzing metric:", metric, "\n")
    
    # Loop through three categories (low, medium, high)
    for (category in 1:3) {
      
      category_label <- c("Low", "Medium", "High")[category]
      cat("    Category:", category_label, "\n")
      
      # Select cities in this category
      cities_subset <- subset(meta_data, get(metric) == category)$code_city
      
      # Filter data for these cities
      data_subset <- subset(data_hw, code_city %in% cities_subset)
      
      # Check if subset has sufficient data
      if (nrow(data_subset) == 0) {
        cat("      Warning: No data available for this category\n")
        results_temp <- data.frame(
          Heatwave_Definition = hw_name,
          Landscape_Metric = metric,
          Category = category_label,
          RR = NA,
          Lower_CI = NA,
          Upper_CI = NA,
          P_value = NA
        )
      } else {
        # Perform analysis
        analysis_results <- analyze_heatwave_mortality(data_subset)
        
        results_temp <- data.frame(
          Heatwave_Definition = hw_name,
          Landscape_Metric = metric,
          Category = category_label,
          RR = analysis_results["RR"],
          Lower_CI = analysis_results["Lower"],
          Upper_CI = analysis_results["Upper"],
          P_value = NA  # Will be calculated later
        )
      }
      
      # Append to results
      results_df <- rbind(results_df, results_temp)
    }
    
    # Calculate p-values comparing Low vs High categories ---------------------
    metric_data <- subset(results_df, Landscape_Metric == metric)
    
    low_data <- subset(metric_data, Category == "Low")
    high_data <- subset(metric_data, Category == "High")
    
    # Calculate p-value for Low vs High comparison
    if (nrow(low_data) == 1 && nrow(high_data) == 1 && 
        !is.na(low_data$RR) && !is.na(high_data$RR)) {
      
      p_low_high <- calculate_pvalue_diff(
        low_data$RR, low_data$Lower_CI, low_data$Upper_CI,
        high_data$RR, high_data$Lower_CI, high_data$Upper_CI
      )
      
      # Update p-values in results
      results_df$P_value[results_df$Landscape_Metric == metric & 
                          results_df$Category == "Low"] <- "Reference"
      results_df$P_value[results_df$Landscape_Metric == metric & 
                          results_df$Category == "Medium"] <- NA
      results_df$P_value[results_df$Landscape_Metric == metric & 
                          results_df$Category == "High"] <- round(p_low_high, 4)
    }
  }
  
  # Store results for this heatwave definition
  results_all[[hw_name]] <- results_df
}

# Combine all results ----------------------------------------------------------
final_results <- do.call(rbind, results_all)
rownames(final_results) <- NULL

# Save results -----------------------------------------------------------------
write.csv(final_results, 
          "/mnt/user-data/outputs/heatwave_greenspace_results_all.csv",
          row.names = FALSE)

cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("Results saved to: heatwave_greenspace_results_all.csv\n")
cat("========================================\n")

################################################################################
# SENSITIVITY ANALYSES (Optional)
################################################################################

# 1. Using daily maximum temperature instead of mean temperature
# 2. Using population-weighted temperature data
# 3. Adjusting for PM2.5 concentrations
# 4. Alternative degrees of freedom for time and humidity
# 5. Using alternative greenspace dataset (Global Forest Change)

# Example: Sensitivity analysis with different df for time variable
sensitivity_df_time <- analyze_heatwave_mortality(data_hw, df_time = 6)

# Example: Sensitivity analysis with different df for humidity
sensitivity_df_humidity <- analyze_heatwave_mortality(data_hw, df_humidity = 6)

################################################################################
# AGE-STRATIFIED ANALYSIS (Optional)
################################################################################

# Analyze separately for age groups: 5-64 years and ≥65 years
# Replace 'death' column with age-specific mortality counts

# data_hw_young <- data_hw
# colnames(data_hw_young)[colnames(data_hw_young) == "death_5_64"] <- "death"
# results_young <- analyze_heatwave_mortality(data_hw_young)

# data_hw_old <- data_hw
# colnames(data_hw_old)[colnames(data_hw_old) == "death_65plus"] <- "death"
# results_old <- analyze_heatwave_mortality(data_hw_old)

################################################################################
# CLIMATE ZONE-STRATIFIED ANALYSIS (Optional)
################################################################################

# Analyze separately for different climate zones
# Filter cities by Köppen climate classification

# climate_zones <- c("Temperate", "Continental", "Arid")
# for (zone in climate_zones) {
#   cities_zone <- subset(meta_data, climate_type == zone)$code_city
#   data_zone <- subset(data_hw, code_city %in% cities_zone)
#   results_zone <- analyze_heatwave_mortality(data_zone)
# }

################################################################################
# END OF SCRIPT
################################################################################
