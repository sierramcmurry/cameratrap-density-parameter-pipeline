# =============================================================================
# 03_ACTIVITY.R
# Chapter 1: Camera Trap Density Estimation
# =============================================================================
#
# PURPOSE: Estimate activity levels (proportion of day active) for each species
#          using circular kernel density estimation (fitact from activity package).
#
# METHOD:  Uses first detection time of each sequence to estimate the 
#          probability density of activity across the 24-hour cycle.
#
# INPUT:   01_data_focal_species.csv from Step 1
# OUTPUT:  Activity estimates by species
#
# =============================================================================

# Load configuration
source("00_config.R")

library(tidyverse)
library(lubridate)
library(activity)

cat("\n")
cat("=============================================================\n")
cat("STEP 03: ACTIVITY ESTIMATION\n")
cat("Using circular kernel density (fitact)\n")
cat("=============================================================\n\n")

# =============================================================================
# 1. LOAD CLEAN DATA
# =============================================================================

cat("Loading data from Step 1...\n")

input_file <- paste0(OUTPUT_DIRS$processed, "01_data_focal_species.csv")
data <- read.csv(input_file, stringsAsFactors = FALSE)
data$timestamp_clean <- as.POSIXct(data$timestamp_clean, tz = "UTC")

cat("  - Loaded", nrow(data), "records\n")
cat("  - Sequences:", n_distinct(data$sequence_id_use), "\n")
cat("  - Species:", n_distinct(data$common_name_clean), "\n\n")

# =============================================================================
# 2. EXTRACT FIRST DETECTION TIME PER SEQUENCE
# =============================================================================

cat("Extracting first detection time per sequence...\n")

activity_data <- data %>%
  group_by(sequence_id_use, deployment_id_clean, common_name_clean) %>%
  summarise(
    first_detection = min(timestamp_clean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hour         = hour(first_detection),
    minute       = minute(first_detection),
    second       = second(first_detection),
    time_decimal = hour + minute/60 + second/3600,
    time_radians = (time_decimal / 24) * 2 * pi
  )

cat("  - Sequences for activity analysis:", nrow(activity_data), "\n")
cat("  - Time range:",
    sprintf("%02d:%02d", min(activity_data$hour), 0), "to",
    sprintf("%02d:%02d", max(activity_data$hour), 59), "\n\n")

# =============================================================================
# 2b. RADIAN VALIDATION
#
#     fitact expects values strictly in [0, 2π]. Values outside this range
#     cause a warning and unreliable activity estimates. Root causes:
#       - Timestamp parsed incorrectly (e.g. hour > 23)
#       - NA or negative time components slipping through
#       - DST or timezone edge cases producing out-of-range decimals
#
#     Fix: flag and remove any sequences with time_radians outside [0, 2π].
#     Report which deployments are affected so they can be investigated.
# =============================================================================

cat("Validating radian values...\n")

valid_range   <- activity_data$time_radians >= 0 &
  activity_data$time_radians <= 2 * pi &
  !is.na(activity_data$time_radians) &
  is.finite(activity_data$time_radians)

n_invalid <- sum(!valid_range)

if(n_invalid > 0) {
  cat("  WARNING:", n_invalid, "sequences have invalid radian values -- removed\n\n")
  
  # Report which deployments and species are affected
  invalid_records <- activity_data[!valid_range, ]
  cat("  Affected sequences:\n")
  invalid_records %>%
    select(sequence_id_use, deployment_id_clean, common_name_clean,
           first_detection, time_decimal, time_radians) %>%
    print(n = 20)
  
  # Flag by deployment for summary
  cat("\n  Affected deployments:\n")
  invalid_records %>%
    group_by(deployment_id_clean, common_name_clean) %>%
    summarise(n_invalid = n(), .groups = "drop") %>%
    print(n = 20)
  
  cat("\n  Possible causes:\n")
  cat("    - hour > 23 or hour < 0: timestamp parsed in wrong timezone\n")
  cat("    - time_decimal > 24 or < 0: DST or UTC offset issue\n")
  cat("    - Check raw timestamps for affected sequences above\n\n")
  
  # Remove invalid records
  activity_data <- activity_data[valid_range, ]
  cat("  Retained:", nrow(activity_data), "valid sequences\n\n")
  
} else {
  cat("  All", nrow(activity_data), "sequences have valid radian values [0, 2π]\n\n")
}

# =============================================================================
# 3. CALCULATE ACTIVITY BY SPECIES
# =============================================================================

cat("Calculating activity levels by species...\n")
cat("  - Minimum sample size:", THRESHOLDS$min_activity_n, "sequences\n\n")

calc_activity <- function(time_radians, min_n = 15, reps = 100) {
  if(length(time_radians) < min_n) {
    return(data.frame(n = length(time_radians),
                      activity = NA, se = NA, lcl = NA, ucl = NA))
  }
  
  fit <- tryCatch({
    fitact(time_radians, reps = reps, sample = "data", show = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if(is.null(fit)) {
    return(data.frame(n = length(time_radians),
                      activity = NA, se = NA, lcl = NA, ucl = NA))
  }
  
  data.frame(
    n        = length(time_radians),
    activity = as.numeric(fit@act[1]),
    se       = as.numeric(fit@act[2]),
    lcl      = as.numeric(fit@act[3]),
    ucl      = as.numeric(fit@act[4])
  )
}

activity_by_species <- activity_data %>%
  group_by(common_name_clean) %>%
  summarise(
    result = list(calc_activity(time_radians,
                                min_n = THRESHOLDS$min_activity_n,
                                reps = 100)),
    .groups = "drop"
  ) %>%
  unnest(result) %>%
  mutate(
    activity_pct = activity * 100,
    active_hours = activity * 24
  ) %>%
  arrange(desc(n))

cat("=== ACTIVITY RESULTS ===\n")
print(activity_by_species %>%
        select(common_name_clean, n, activity, se, lcl, ucl,
               activity_pct, active_hours),
      n = 20)

# =============================================================================
# 4. ACTIVITY BY DEPLOYMENT
# =============================================================================

cat("\nCalculating activity by deployment (for species with enough data)...\n")

abundant_species <- activity_by_species %>%
  filter(n >= 50) %>%
  pull(common_name_clean)

cat("  - Species with n ≥ 50:", paste(abundant_species, collapse = ", "), "\n")

activity_by_deployment <- activity_data %>%
  filter(common_name_clean %in% abundant_species) %>%
  group_by(deployment_id_clean, common_name_clean) %>%
  summarise(
    n_seqs = n(),
    result = list(calc_activity(time_radians, min_n = 10, reps = 50)),
    .groups = "drop"
  ) %>%
  unnest(result) %>%
  rename(n_detections = n_seqs)

cat("  - Deployment-level estimates:", nrow(activity_by_deployment), "\n\n")

# =============================================================================
# 5. SUMMARY STATISTICS
# =============================================================================

cat("=============================================================\n")
cat("ACTIVITY SUMMARY\n")
cat("=============================================================\n")

valid_species <- activity_by_species %>% filter(!is.na(activity))

cat("Species with valid estimates:", nrow(valid_species), "of",
    nrow(activity_by_species), "\n\n")

cat("Activity ranges:\n")
cat("  - Most active:",
    valid_species$common_name_clean[which.max(valid_species$activity)],
    sprintf("(%.1f%%)\n", max(valid_species$activity_pct, na.rm = TRUE)))
cat("  - Least active:",
    valid_species$common_name_clean[which.min(valid_species$activity)],
    sprintf("(%.1f%%)\n", min(valid_species$activity_pct, na.rm = TRUE)))
cat("  - Mean activity:",
    sprintf("%.1f%%\n", mean(valid_species$activity_pct, na.rm = TRUE)))

cat("=============================================================\n\n")

# =============================================================================
# 6. SAVE RESULTS
# =============================================================================

cat("Saving results...\n")

write.csv(activity_by_species,
          paste0(OUTPUT_DIRS$processed, "03_activity_by_species.csv"),
          row.names = FALSE)
cat("  - Species-level:", paste0(OUTPUT_DIRS$processed, "03_activity_by_species.csv"), "\n")

write.csv(activity_by_deployment,
          paste0(OUTPUT_DIRS$processed, "03_activity_by_deployment.csv"),
          row.names = FALSE)
cat("  - Deployment-level:", paste0(OUTPUT_DIRS$processed, "03_activity_by_deployment.csv"), "\n")

write.csv(activity_data,
          paste0(OUTPUT_DIRS$processed, "03_activity_detection_times.csv"),
          row.names = FALSE)
cat("  - Detection times:", paste0(OUTPUT_DIRS$processed, "03_activity_detection_times.csv"), "\n")

saveRDS(list(
  activity_by_species    = activity_by_species,
  activity_by_deployment = activity_by_deployment,
  activity_data          = activity_data
), paste0(OUTPUT_DIRS$processed, "03_activity_results.rds"))
cat("  - All results (RDS):", paste0(OUTPUT_DIRS$processed, "03_activity_results.rds"), "\n")

# =============================================================================
# 7. CALCULATE DAILY MOVEMENT RATE
# =============================================================================

cat("\nCalculating daily movement rate...\n")

sbd_results <- readRDS(paste0(OUTPUT_DIRS$processed, "02_sbd_results.rds"))
sbd_speeds  <- sbd_results$sbd_results

daily_movement <- sbd_speeds %>%
  select(common_name_clean, sbd_mean, sbd_se, sbd_lcl, sbd_ucl, n_speed = n) %>%
  left_join(
    activity_by_species %>%
      select(common_name_clean, activity, activity_se = se,
             activity_lcl = lcl, activity_ucl = ucl, n_activity = n),
    by = "common_name_clean"
  ) %>%
  mutate(
    dmr_km_day = sbd_mean * activity * 86400 / 1000,
    dmr_se     = dmr_km_day * sqrt((sbd_se / sbd_mean)^2 +
                                     (activity_se / activity)^2),
    dmr_lcl    = dmr_km_day - 1.96 * dmr_se,
    dmr_ucl    = dmr_km_day + 1.96 * dmr_se,
    dmr_m_day  = sbd_mean * activity * 86400
  ) %>%
  filter(!is.na(dmr_km_day)) %>%
  arrange(desc(n_speed))

cat("\n=== DAILY MOVEMENT RATE RESULTS ===\n")
print(daily_movement %>%
        select(common_name_clean, sbd_mean, activity,
               dmr_km_day, dmr_se, dmr_lcl, dmr_ucl),
      n = 15)

cat("\n=== DAILY MOVEMENT SUMMARY ===\n")
cat("  Highest:", daily_movement$common_name_clean[which.max(daily_movement$dmr_km_day)],
    sprintf("(%.1f km/day)\n", max(daily_movement$dmr_km_day, na.rm = TRUE)))
cat("  Lowest:", daily_movement$common_name_clean[which.min(daily_movement$dmr_km_day)],
    sprintf("(%.1f km/day)\n", min(daily_movement$dmr_km_day, na.rm = TRUE)))

# =============================================================================
# 8. SAVE DAILY MOVEMENT RESULTS
# =============================================================================

write.csv(daily_movement,
          paste0(OUTPUT_DIRS$processed, "03_daily_movement_rate.csv"),
          row.names = FALSE)
cat("\n  - Daily movement rate:", paste0(OUTPUT_DIRS$processed, "03_daily_movement_rate.csv"), "\n")

saveRDS(list(
  activity_by_species    = activity_by_species,
  activity_by_deployment = activity_by_deployment,
  activity_data          = activity_data,
  daily_movement         = daily_movement
), paste0(OUTPUT_DIRS$processed, "03_activity_results.rds"))
cat("  - Updated RDS with daily movement\n")

cat("\n")
cat("=============================================================\n")
cat("STEP 03 COMPLETE\n")
cat("=============================================================\n")
cat("\nNext: Run 04_staying_time.R\n\n")

# =============================================================================
# OBJECTS AVAILABLE FOR NEXT SCRIPT
# =============================================================================
# activity_by_species    - Activity estimates per species
# activity_by_deployment - Activity by deployment (abundant species only)
# activity_data          - Raw detection times for figures (validated)
# daily_movement         - DMR combining speed and activity