# =============================================================================
# 04_STAYING_TIME.R - CORRECTED
# Chapter 1: Camera Trap Density Estimation
# =============================================================================
#
# PURPOSE: Calculate staying time (time in view) for each species using
#          bounding box interpolation method.
#
# METHOD:  Based on original Western parameters script:
#          1. Parse bounding boxes and compute edge flags
#          2. Create 1-second timeline for each sequence
#          3. Interpolate gaps using edge detection logic
#          4. Count PRESENCE-SECONDS (time any animal in view)
#
# CORRECTIONS:
#   - Uses PRESENCE-SECONDS (sum of seconds with any detection)
#   - NOT animal-seconds (which would multiply by group size)
#   - Collapses burst photos: takes max bboxes per image, then max per timestamp
#   - 3 deer for 10 seconds = 10 presence-seconds (the group passage time)
#
# INPUT:   01_data_focal_species.csv from Step 1
# OUTPUT:  Staying time estimates by species and deployment
#
# =============================================================================

# Load configuration
source("00_config.R")

library(tidyverse)
library(lubridate)

cat("\n")
cat("=============================================================\n")
cat("STEP 04: STAYING TIME ESTIMATION (CORRECTED)\n")
cat("Using bounding box interpolation method - ANIMAL-SECONDS\n")
cat("=============================================================\n\n")

# =============================================================================
# 1. LOAD CLEAN DATA
# =============================================================================

cat("Loading data from Step 1...\n")

input_file <- paste0(OUTPUT_DIRS$processed, "01_data_for_density.csv")
data <- read.csv(input_file, stringsAsFactors = FALSE)
data$timestamp_clean <- as.POSIXct(data$timestamp_clean, tz = "UTC")

cat("  - Loaded", nrow(data), "records\n")
cat("  - Sequences:", n_distinct(data$sequence_id_use), "\n")
cat("  - Species:", n_distinct(data$common_name_clean), "\n\n")

# =============================================================================
# 2. HELPER FUNCTIONS (from original Western parameters script)
# =============================================================================

# Compute edge flag from parsed bbox coordinates
# UPDATED: uses xmin/xmax column names
compute_edge_flag <- function(xmin, xmax, edge_fraction = 1/6) {
  # Check if animal is on left or right edge of frame
  left <- xmin < edge_fraction
  right <- xmax > (1 - edge_fraction)
  
  case_when(
    left & !right ~ "left",
    right & !left ~ "right",
    TRUE ~ NA_character_
  )
}

# Create 1-second timeline for a sequence
create_1s_timeline <- function(df, time_col = "timestamp_clean") {
  df <- df[order(df[[time_col]]), ]
  df <- df[!duplicated(df[[time_col]]), ]
  
  min_time <- min(df[[time_col]], na.rm = TRUE)
  max_time <- max(df[[time_col]], na.rm = TRUE)
  
  # Create 1-second intervals
  times <- seq(min_time, max_time, by = "1 sec")
  
  # Merge with original data
  timeline <- data.frame(timestamp_clean = times)
  merged <- merge(timeline, df, by = "timestamp_clean", all.x = TRUE)
  
  return(merged)
}

# Interpolate gaps using edge detection logic
interpolate_gaps_with_edge <- function(df, rep_group_size,
                                       box_col = "num_bboxes",
                                       edge_col = "edge_flag") {
  df <- df[order(df$timestamp_clean), ]
  box_vals <- df[[box_col]]
  edge_vals <- df[[edge_col]]
  
  idx <- which(!is.na(box_vals))
  
  if (length(idx) < 2) {
    df[[box_col]][is.na(box_vals)] <- 0
    return(df)
  }
  
  for (i in seq_len(length(idx) - 1)) {
    s_i <- idx[i]
    e_i <- idx[i + 1]
    gap <- e_i - s_i - 1
    
    if (gap <= 0) next
    
    b_s <- box_vals[s_i]
    b_e <- box_vals[e_i]
    e_s <- edge_vals[s_i]
    e_e <- edge_vals[e_i]
    
    # Short gaps with detections on both ends: assume animal stayed
    if (gap <= 10 && b_s > 0 && b_e > 0) {
      box_vals[(s_i + 1):(e_i - 1)] <- rep_group_size
    }
    # Same edge entry/exit: assume animal left
    else if (!is.na(e_s) && !is.na(e_e) && e_s == e_e) {
      box_vals[(s_i + 1):(e_i - 1)] <- 0
    }
    # Otherwise: interpolate
    else {
      iv <- round(seq(b_s, b_e, length.out = gap + 2)[-c(1, gap + 2)])
      iv[iv < 1] <- 1
      box_vals[(s_i + 1):(e_i - 1)] <- iv
    }
  }
  
  df[[box_col]] <- ifelse(is.na(box_vals), 0, box_vals)
  return(df)
}

# =============================================================================
# 3. PREPARE DATA WITH EDGE FLAGS AND BBOX COUNTS
# =============================================================================

cat("Preparing data with edge flags and bbox counts...\n")

# CORRECTED: Handle burst photos by collapsing to max bboxes per timestamp
# Step 1: Count bboxes per image, compute edge flag
images_df <- data %>%
  mutate(
    edge_flag = compute_edge_flag(xmin, xmax, edge_fraction = 1/6)
  ) %>%
  # Group by image to count bboxes per image
  group_by(sequence_id_use, deployment_id_clean, common_name_clean, timestamp_clean, image_id) %>%
  summarise(
    bboxes_in_image = n(),
    edge_flag = first(na.omit(edge_flag)),
    .groups = "drop"
  ) %>%
  # Step 2: For each timestamp, take max bboxes (collapses burst photos)
  group_by(sequence_id_use, deployment_id_clean, common_name_clean, timestamp_clean) %>%
  summarise(
    num_bboxes = max(bboxes_in_image),  # Max animals seen in any single image at this timestamp
    edge_flag = first(na.omit(edge_flag)),
    .groups = "drop"
  )

cat("  - Prepared", nrow(images_df), "timestamp records\n")
cat("  - Unique sequences:", n_distinct(images_df$sequence_id_use), "\n\n")

# =============================================================================
# 4. EXPAND TIMELINE AND INTERPOLATE GAPS
# =============================================================================

cat("Creating 1-second timelines and interpolating gaps...\n")
cat("  (This may take a few minutes...)\n")

# Process each sequence
processed <- images_df %>%
  group_by(deployment_id_clean, common_name_clean, sequence_id_use) %>%
  group_modify(~ {
    # Create 1-second timeline
    tl <- create_1s_timeline(.x, time_col = "timestamp_clean")
    
    # Calculate representative group size for interpolation
    rep_size <- if (sum(.x$num_bboxes > 0, na.rm = TRUE) > 0) {
      round(mean(.x$num_bboxes[.x$num_bboxes > 0], na.rm = TRUE))
    } else {
      0
    }
    
    # Interpolate gaps
    interpolate_gaps_with_edge(tl, rep_group_size = rep_size)
  }) %>%
  ungroup()

cat("  - Processed timeline records:", nrow(processed), "\n\n")

# =============================================================================
# 5. CALCULATE STAYING TIME PER SEQUENCE
# =============================================================================

cat("Calculating staying time per sequence...\n")

# =============================================================================
# CORRECTED: Staying time = PRESENCE-SECONDS (time any animal present)
# NOT: sum(num_bboxes) which would be animal-seconds
#
# Example: 3 deer for 10 seconds = 10 presence-seconds (not 30)
# This matches the REST formula where s = time per event (group passage)
# =============================================================================
stay_seq <- processed %>%
  group_by(deployment_id_clean, common_name_clean, sequence_id_use) %>%
  summarise(
    staying_time = sum(num_bboxes > 0, na.rm = TRUE),  # PRESENCE-SECONDS
    animal_seconds = sum(num_bboxes, na.rm = TRUE),     # Keep for reference
    .groups = "drop"
  )

cat("  - Calculated staying time for", nrow(stay_seq), "sequences\n\n")

# =============================================================================
# 6. SUMMARIZE BY SPECIES
# =============================================================================

cat("Summarizing staying time by species...\n")

staying_by_species <- stay_seq %>%
  group_by(common_name_clean) %>%
  summarise(
    sample_size = n(),
    mean_staying_time = mean(staying_time, na.rm = TRUE),
    median_staying_time = median(staying_time, na.rm = TRUE),
    sd_staying_time = sd(staying_time, na.rm = TRUE),
    q025 = quantile(staying_time, 0.025, na.rm = TRUE),
    q975 = quantile(staying_time, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sample_size))

cat("\n=== STAYING TIME RESULTS (presence-seconds) ===\n")
print(staying_by_species, n = 20)

# =============================================================================
# 7. STAYING TIME BY DEPLOYMENT
# =============================================================================

cat("\nCalculating staying time by deployment...\n")

# Only for species with good sample sizes
abundant_species <- staying_by_species %>%
  filter(sample_size >= 20) %>%
  pull(common_name_clean)

cat("  - Species with n ≥ 20:", length(abundant_species), "\n")

staying_by_deployment <- stay_seq %>%
  filter(common_name_clean %in% abundant_species) %>%
  group_by(deployment_id_clean, common_name_clean) %>%
  summarise(
    sample_size = n(),
    mean_staying_time = mean(staying_time, na.rm = TRUE),
    q025 = quantile(staying_time, 0.025, na.rm = TRUE),
    q975 = quantile(staying_time, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(sample_size >= 3)

cat("  - Deployment-level estimates:", nrow(staying_by_deployment), "\n\n")

# =============================================================================
# 8. SUMMARY STATISTICS
# =============================================================================

cat("=============================================================\n")
cat("STAYING TIME SUMMARY (PRESENCE-SECONDS)\n")
cat("=============================================================\n")

valid_species <- staying_by_species %>% filter(sample_size >= 10)

cat("Species with ≥10 sequences:", nrow(valid_species), "\n\n")

cat("Staying time ranges (animal-seconds):\n")
cat("  - Longest:", 
    valid_species$common_name_clean[which.max(valid_species$mean_staying_time)],
    sprintf("(%.1f)\n", max(valid_species$mean_staying_time, na.rm = TRUE)))
cat("  - Shortest:", 
    valid_species$common_name_clean[which.min(valid_species$mean_staying_time)],
    sprintf("(%.1f)\n", min(valid_species$mean_staying_time, na.rm = TRUE)))
cat("  - Overall mean:", sprintf("%.1f\n", mean(valid_species$mean_staying_time, na.rm = TRUE)))

cat("=============================================================\n\n")

# =============================================================================
# 9. SAVE RESULTS
# =============================================================================

cat("Saving results...\n")

# Staying time by species (main results)
output_species <- paste0(OUTPUT_DIRS$processed, "04_staying_time_by_species.csv")
write.csv(staying_by_species, output_species, row.names = FALSE)
cat("  - Species-level:", output_species, "\n")

# Staying time by deployment
output_deployment <- paste0(OUTPUT_DIRS$processed, "04_staying_time_by_deployment.csv")
write.csv(staying_by_deployment, output_deployment, row.names = FALSE)
cat("  - Deployment-level:", output_deployment, "\n")

# Sequence-level data (for figures/diagnostics)
output_seq <- paste0(OUTPUT_DIRS$processed, "04_staying_time_sequences.csv")
write.csv(stay_seq, output_seq, row.names = FALSE)
cat("  - Sequence-level:", output_seq, "\n")

# Simplified version for density model
staying_for_density <- staying_by_species %>%
  select(species = common_name_clean, 
         staying_time = mean_staying_time,
         n = sample_size)

output_density <- paste0(OUTPUT_DIRS$processed, "04_staying_times_for_density.csv")
write.csv(staying_for_density, output_density, row.names = FALSE)
cat("  - For density model:", output_density, "\n")

# Save as RDS for next scripts
output_rds <- paste0(OUTPUT_DIRS$processed, "04_staying_time_results.rds")
saveRDS(list(
  staying_by_species = staying_by_species,
  staying_by_deployment = staying_by_deployment,
  stay_seq = stay_seq
), output_rds)
cat("  - All results (RDS):", output_rds, "\n")

cat("\n")
cat("=============================================================\n")
cat("STEP 04 COMPLETE\n")
cat("=============================================================\n")
cat("\nNext: Run 05_edd.R\n\n")

# =============================================================================
# OBJECTS AVAILABLE FOR NEXT SCRIPT
# =============================================================================
# staying_by_species     - Staying time estimates per species
# staying_by_deployment  - Staying time by deployment
# stay_seq               - Sequence-level staying time data