# =============================================================================
# 01_LOAD_AND_JOIN.R
# Chapter 1: Camera Trap Density Estimation
# 
# PURPOSE: Load and combine AI coordinate data, apply basic cleaning,
#          flag edge cases at frame and sequence level
# INPUT: Raw CSV files from Mohammad's AI pipeline
# OUTPUT: Combined, cleaned dataset ready for analysis
# =============================================================================

# Load configuration
source("00_config.R")

cat("\n")
cat("=============================================================\n")
cat("STEP 01: LOAD, COMBINE, AND FILTER DATA\n")
cat("=============================================================\n\n")

# =============================================================================
# 1. LOAD RAW DATA
# =============================================================================

cat("Loading AI coordinate data...\n")

browning <- read.csv(COORD_FILES$browning_elite, stringsAsFactors = FALSE)
cat("  Browning Recon Force Elite:", nrow(browning), "rows\n")

reconyx <- read.csv(COORD_FILES$reconyx, stringsAsFactors = FALSE)
cat("  Reconyx Hyperfire:", nrow(reconyx), "rows\n")

# Combine both camera types
data_raw <- bind_rows(browning, reconyx)
cat("  Combined total:", nrow(data_raw), "rows\n\n")

# =============================================================================
# 2. BASIC CLEANING
# =============================================================================

cat("Applying basic cleaning...\n")

data_clean <- data_raw %>%
  # Standardize column names for consistency with old code
  rename(
    common_name_clean = common_name,
    deployment_id_clean = deployment_id,
    sequence_id_use = sequence_id
  ) %>%
  # Parse timestamp
  mutate(
    timestamp_clean = case_when(
      # Format 1: "2024-09-23 19:53:25"
      !is.na(as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) ~
        as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      # Format 2: "2024-08-26T05:17:00Z"
      !is.na(as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")) ~
        as.POSIXct(timestamp, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      TRUE ~ as.POSIXct(NA)
    ),
    # Extract hour for activity analysis
    hour = hour(timestamp_clean) + minute(timestamp_clean)/60
  ) %>%
  # Remove rows with missing critical data
  filter(
    !is.na(world_z),
    !is.na(world_x),
    !is.na(timestamp_clean)
  )

cat("  After removing NAs:", nrow(data_clean), "rows\n")

# =============================================================================
# 3. USE AI BOUNDING BOXES (pre-normalized by Mohammad's pipeline)
# =============================================================================
#
# The x1_n, y1_n, x2_n, y2_n columns contain MegaDetector bounding boxes
# normalized to 0-1 by Mohammad's pipeline (pixel coords / image dimensions).
# These replace the Wildlife Insights bounding_boxes JSON column.
# =============================================================================

cat("Loading AI bounding boxes (pre-normalized)...\n")

data_clean <- data_clean %>%
  rename(xmin = x1_n, ymin = y1_n, xmax = x2_n, ymax = y2_n) %>%
  mutate(
    # Clamp to [0, 1] — MegaDetector boxes can extend past image edges
    xmin = pmax(0, pmin(1, xmin)),
    ymin = pmax(0, pmin(1, ymin)),
    xmax = pmax(0, pmin(1, xmax)),
    ymax = pmax(0, pmin(1, ymax)),
    # Derived bbox metrics
    bbox_width = xmax - xmin,
    bbox_height = ymax - ymin,
    bbox_area_pct = bbox_width * bbox_height * 100
  )

cat("  Bbox area range:", round(range(data_clean$bbox_area_pct, na.rm = TRUE), 2), "%\n")
cat("  Bbox value ranges - xmin:", round(range(data_clean$xmin, na.rm = TRUE), 3),
    "xmax:", round(range(data_clean$xmax, na.rm = TRUE), 3), "\n\n")
# =============================================================================
# 4. FLAG EDGE CASES (FRAME LEVEL)
# =============================================================================

cat("Flagging frame-level edge cases...\n")

data_clean <- data_clean %>%
  mutate(
    # Flag edge cases
    flag_bbox_large = bbox_area_pct > THRESHOLDS$max_bbox_pct,
    flag_too_close = world_z < THRESHOLDS$min_distance,
    flag_too_far = world_z > THRESHOLDS$max_distance
  )

# Report edge case counts
cat("\n=== FRAME-LEVEL EDGE CASES ===\n")
cat("  Bbox >", THRESHOLDS$max_bbox_pct, "%:", 
    sum(data_clean$flag_bbox_large, na.rm = TRUE), 
    "(", round(100 * mean(data_clean$flag_bbox_large, na.rm = TRUE), 2), "%)\n")
cat("  Too close (<", THRESHOLDS$min_distance, "m):", 
    sum(data_clean$flag_too_close, na.rm = TRUE),
    "(", round(100 * mean(data_clean$flag_too_close, na.rm = TRUE), 2), "%)\n")
cat("  Too far (>", THRESHOLDS$max_distance, "m):", 
    sum(data_clean$flag_too_far, na.rm = TRUE),
    "(", round(100 * mean(data_clean$flag_too_far, na.rm = TRUE), 2), "%)\n")

# =============================================================================
# 5. CALCULATE FRAME-TO-FRAME SPEED & FLAG IMPOSSIBLE SPEEDS
# =============================================================================

cat("\nCalculating frame-to-frame movements...\n")

# Extract frame number from filename for correct within-sequence ordering.
# Timestamps alone are unreliable -- multiple frames share the same second
# (e.g. all frames in sequence 9486720 are "16:37"). Sorting by timestamp
# only gives arbitrary order within a second, which scrambles the path and
# inflates both path distance and straight-line displacement.
# Filenames encode true capture order (e.g. IMG_0330.JPG < IMG_0332.JPG).
# We extract the numeric portion and use it as a tiebreaker after timestamp.

data_clean <- data_clean %>%
  mutate(
    frame_number = as.numeric(str_extract(filename, "\\d+(?=\\.[A-Za-z]{3,4}$)"))
  )

n_missing_frame_num <- sum(is.na(data_clean$frame_number))
if(n_missing_frame_num > 0) {
  cat("  Warning:", n_missing_frame_num,
      "rows could not extract frame number from filename\n")
  cat("  These will sort by timestamp only -- check filename format\n")
}

data_clean <- data_clean %>%
  arrange(sequence_id_use, timestamp_clean, frame_number) %>%
  group_by(sequence_id_use) %>%
  mutate(
    # Previous frame coordinates
    prev_x = lag(world_x),
    prev_z = lag(world_z),
    prev_time = lag(timestamp_clean),
    
    # Frame-to-frame distance
    frame_distance = sqrt((world_x - prev_x)^2 + (world_z - prev_z)^2),
    
    # Frame-to-frame time
    frame_time = as.numeric(difftime(timestamp_clean, prev_time, units = "secs")),
    
    # Frame-to-frame speed
    frame_speed = ifelse(frame_time > 0, frame_distance / frame_time, NA),
    
    # Flag impossible speeds
    flag_impossible_speed = !is.na(frame_speed) & frame_speed > THRESHOLDS$max_frame_speed
  ) %>%
  ungroup()

n_impossible <- sum(data_clean$flag_impossible_speed, na.rm = TRUE)
cat("  Flagged", n_impossible, "frame transitions with speed >", 
    THRESHOLDS$max_frame_speed, "m/s\n")

if(n_impossible > 0) {
  # Show which species have impossible speeds
  cat("  By species:\n")
  impossible_by_species <- data_clean %>%
    filter(flag_impossible_speed) %>%
    count(common_name_clean, sort = TRUE)
  print(impossible_by_species, n = 10)
}

# =============================================================================
# 6. COMBINED FRAME-LEVEL FLAG
# =============================================================================

data_clean <- data_clean %>%
  mutate(
    edge_case_flag = flag_bbox_large | flag_too_close | flag_too_far | flag_impossible_speed
  )

cat("\n=== TOTAL FRAME-LEVEL FLAGS ===\n")
cat("  Total flagged:", 
    sum(data_clean$edge_case_flag, na.rm = TRUE),
    "(", round(100 * mean(data_clean$edge_case_flag, na.rm = TRUE), 2), "%)\n")

# =============================================================================
# 7. SEQUENCE-LEVEL FLAGS
# =============================================================================

cat("\nCreating sequence-level quality flags...\n")

sequence_flags <- data_clean %>%
  group_by(sequence_id_use) %>%
  summarise(
    n_frames = n(),
    has_bbox_issue = any(flag_bbox_large, na.rm = TRUE),
    has_close_detection = any(flag_too_close, na.rm = TRUE),
    has_far_detection = any(flag_too_far, na.rm = TRUE),
    has_impossible_speed = any(flag_impossible_speed, na.rm = TRUE),
    pct_frames_flagged = mean(edge_case_flag, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    # Overall sequence quality flag
    seq_any_issue = has_bbox_issue | has_close_detection | has_far_detection | has_impossible_speed
  )

cat("=== SEQUENCE-LEVEL FLAGS ===\n")
cat("  Total sequences:", nrow(sequence_flags), "\n")
cat("  Sequences with any issue:", sum(sequence_flags$seq_any_issue), 
    "(", round(100 * mean(sequence_flags$seq_any_issue), 1), "%)\n")
cat("    - Bbox issues:", sum(sequence_flags$has_bbox_issue), "\n")
cat("    - Close detections:", sum(sequence_flags$has_close_detection), "\n")
cat("    - Far detections:", sum(sequence_flags$has_far_detection), "\n")
cat("    - Impossible speeds:", sum(sequence_flags$has_impossible_speed), "\n")

# Join sequence flags back to main data
data_clean <- data_clean %>%
  left_join(
    sequence_flags %>% select(sequence_id_use, seq_any_issue, pct_frames_flagged),
    by = "sequence_id_use"
  )

# =============================================================================
# 8. FILTER TO FOCAL SPECIES
# =============================================================================

cat("\nFiltering to focal species...\n")

data_focal <- data_clean %>%
  filter(common_name_clean %in% FOCAL_SPECIES)

cat("  Focal species detections:", nrow(data_focal), "\n")

# Species summary
species_summary <- data_focal %>%
  group_by(common_name_clean) %>%
  summarise(
    n_detections = n(),
    n_sequences = n_distinct(sequence_id_use),
    n_deployments = n_distinct(deployment_id_clean),
    n_flagged = sum(edge_case_flag, na.rm = TRUE),
    pct_flagged = round(100 * mean(edge_case_flag, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_detections))

cat("\n=== FOCAL SPECIES SUMMARY ===\n")
print(species_summary, n = 25)

# =============================================================================
# 9. CREATE ANALYSIS-SPECIFIC DATASETS
# =============================================================================

cat("\nCreating analysis-specific datasets...\n")

# Report excluded deployments
cat("  Excluding", length(EXCLUDE_DEPLOYMENTS), "deployments with calibration issues\n")
cat("    (", paste(EXCLUDE_DEPLOYMENTS, collapse = ", "), ")\n")

n_excluded <- sum(data_focal$deployment_id_clean %in% EXCLUDE_DEPLOYMENTS)
cat("    Detections excluded:", n_excluded, "\n")

# -----------------------------------------------------------------------------
# For DENSITY/EDD analysis:
#   - Exclude calibration-problem deployments
#   - Filter to max distance (26m)
#   - Set close detections to 0.5m
#   - Remove impossible speeds
# -----------------------------------------------------------------------------
data_density <- data_focal %>%
  filter(
    !deployment_id_clean %in% EXCLUDE_DEPLOYMENTS_ALL_COORDS,
    world_z <= THRESHOLDS$max_distance,
    !flag_impossible_speed
  ) %>%
  mutate(
    world_z_adj = case_when(
      flag_bbox_large | flag_too_close ~ 0.5,
      TRUE ~ world_z
    )
  )

cat("  Density dataset:", nrow(data_density), "rows\n")
cat("    (calibration issues excluded, <=26m, close->0.5m, impossible speeds removed)\n")

# -----------------------------------------------------------------------------
# For SPEED/SBD analysis:
#   Sequence-level exclusions (entire sequence dropped if any violation):
#
#   FIX 1 - bbox: If ANY frame in a sequence has bbox_area_pct > 80%, drop
#   the whole sequence. The animal was too close for valid coordinates.
#
#   FIX 2 - group size: Drop any sequence where group_size > 1 according to
#   the Wildlife Insights sequences metadata. Cannot track individual movement
#   when multiple animals are present -- the coordinates jump between
#   individuals across timestamps, inflating path distance.
# -----------------------------------------------------------------------------

# Identify sequences to exclude for speed analysis
seqs_with_large_bbox <- data_focal %>%
  filter(flag_bbox_large) %>%
  distinct(sequence_id_use) %>%
  pull(sequence_id_use)

# Load group_size from sequences metadata and exclude multi-animal sequences
sequences_meta <- read.csv(WI_FILES$sequences, stringsAsFactors = FALSE) %>%
  select(sequence_id, group_size) %>%
  mutate(sequence_id = as.character(sequence_id))

seqs_with_multi_animal <- sequences_meta %>%
  filter(group_size > 1) %>%
  pull(sequence_id)

cat("\n  Speed dataset sequence-level exclusions:\n")
cat("    Sequences with any bbox >", THRESHOLDS$max_bbox_pct, "%:",
    length(seqs_with_large_bbox), "\n")
cat("    Sequences with group_size > 1:",
    length(seqs_with_multi_animal), "\n")
cat("    (overlap between these two groups is allowed)\n")

data_speed <- data_focal %>%
  filter(
    !deployment_id_clean %in% EXCLUDE_DEPLOYMENTS_ALL_COORDS,
    !sequence_id_use %in% seqs_with_large_bbox,
    !sequence_id_use %in% seqs_with_multi_animal,
    !edge_case_flag | is.na(edge_case_flag)
  )

cat("  Speed dataset:", nrow(data_speed), "rows\n")
cat("    Sequences remaining:", n_distinct(data_speed$sequence_id_use), "\n")
cat("    (calibration deployments excluded, bbox sequences excluded,\n")
cat("     group size > 1 excluded, frame-level edge cases removed)\n")

# -----------------------------------------------------------------------------
# For ACTIVITY analysis:
#   - Keep ALL deployments (only uses timestamps, not coordinates)
# -----------------------------------------------------------------------------
data_activity <- data_focal

cat("  Activity dataset:", nrow(data_activity), "rows\n")
cat("    (all deployments included - activity only uses timestamps)\n")

# =============================================================================
# 10. SAVE PROCESSED DATA
# =============================================================================

cat("\nSaving processed data...\n")

# Save full data with flags (all detections)
write.csv(data_focal, 
          paste0(OUTPUT_DIRS$processed, "01_data_focal_species.csv"),
          row.names = FALSE)
cat("  Saved: 01_data_focal_species.csv\n")

# Save density-ready data
write.csv(data_density,
          paste0(OUTPUT_DIRS$processed, "01_data_for_density.csv"),
          row.names = FALSE)
cat("  Saved: 01_data_for_density.csv\n")

# Save speed-ready data
write.csv(data_speed,
          paste0(OUTPUT_DIRS$processed, "01_data_for_speed.csv"),
          row.names = FALSE)
cat("  Saved: 01_data_for_speed.csv\n")

# Save activity-ready data (same as focal but explicit)
write.csv(data_activity,
          paste0(OUTPUT_DIRS$processed, "01_data_for_activity.csv"),
          row.names = FALSE)
cat("  Saved: 01_data_for_activity.csv\n")

# Save sequence flags
write.csv(sequence_flags,
          paste0(OUTPUT_DIRS$processed, "01_sequence_flags.csv"),
          row.names = FALSE)
cat("  Saved: 01_sequence_flags.csv\n")

# Save species summary
write.csv(species_summary,
          paste0(OUTPUT_DIRS$tables, "01_species_summary.csv"),
          row.names = FALSE)
cat("  Saved: 01_species_summary.csv\n")

# =============================================================================
# 11. FILTERING SUMMARY TABLE
# =============================================================================

cat("\n")
cat("=============================================================\n")
cat("FILTERING SUMMARY\n")
cat("=============================================================\n")

summary_table <- data.frame(
  Metric = c(
    "Raw detections loaded",
    "After removing NAs",
    "Focal species detections",
    "Frame-level flags:",
    "  - Bbox >80%",
    paste0("  - Too close (<", THRESHOLDS$min_distance, "m)"),
    paste0("  - Too far (>", THRESHOLDS$max_distance, "m)"),
    "  - Impossible speed (>20 m/s)",
    "  - Total frames flagged",
    "Sequence-level:",
    "  - Total sequences",
    "  - Sequences with any issue",
    "Deployment exclusions:",
    "  - Calibration problems",
    "  - Detections excluded",
    "Speed dataset exclusions (sequence-level):",
    "  - Sequences with any bbox >80%",
    "  - Sequences with group_size > 1",
    "Output datasets:",
    "  - For density analysis",
    "  - For speed analysis",
    "  - For activity analysis"
  ),
  Count = c(
    nrow(data_raw),
    nrow(data_clean),
    nrow(data_focal),
    NA,
    sum(data_focal$flag_bbox_large, na.rm = TRUE),
    sum(data_focal$flag_too_close, na.rm = TRUE),
    sum(data_focal$flag_too_far, na.rm = TRUE),
    sum(data_focal$flag_impossible_speed, na.rm = TRUE),
    sum(data_focal$edge_case_flag, na.rm = TRUE),
    NA,
    n_distinct(data_focal$sequence_id_use),
    sum(sequence_flags$seq_any_issue[
      sequence_flags$sequence_id_use %in% data_focal$sequence_id_use]),
    NA,
    length(EXCLUDE_DEPLOYMENTS),
    n_excluded,
    NA,
    length(seqs_with_large_bbox),
    length(seqs_with_multi_animal),
    NA,
    nrow(data_density),
    nrow(data_speed),
    nrow(data_activity)
  ),
  Percent = c(
    "100%",
    paste0(round(100 * nrow(data_clean) / nrow(data_raw), 1), "%"),
    paste0(round(100 * nrow(data_focal) / nrow(data_clean), 1), "%"),
    "",
    paste0(round(100 * mean(data_focal$flag_bbox_large, na.rm = TRUE), 2), "%"),
    paste0(round(100 * mean(data_focal$flag_too_close, na.rm = TRUE), 2), "%"),
    paste0(round(100 * mean(data_focal$flag_too_far, na.rm = TRUE), 2), "%"),
    paste0(round(100 * mean(data_focal$flag_impossible_speed, na.rm = TRUE), 2), "%"),
    paste0(round(100 * mean(data_focal$edge_case_flag, na.rm = TRUE), 2), "%"),
    "",
    "100%",
    paste0(round(100 * mean(sequence_flags$seq_any_issue[
      sequence_flags$sequence_id_use %in% data_focal$sequence_id_use]), 1), "%"),
    "",
    "",
    paste0(round(100 * n_excluded / nrow(data_focal), 1), "%"),
    "",
    paste0(round(100 * length(seqs_with_large_bbox) /
                   n_distinct(data_focal$sequence_id_use), 1), "% of sequences"),
    paste0(round(100 * length(seqs_with_multi_animal) /
                   n_distinct(data_focal$sequence_id_use), 1), "% of sequences"),
    "",
    paste0(round(100 * nrow(data_density) / nrow(data_focal), 1), "% of focal"),
    paste0(round(100 * nrow(data_speed) / nrow(data_focal), 1), "% of focal"),
    "100% of focal"
  )
)

print(summary_table, row.names = FALSE)

cat("=============================================================\n\n")

# =============================================================================
# 12. CLEANUP & FINAL MESSAGE
# =============================================================================

# Clean up intermediate objects
rm(browning, reconyx, data_raw, bbox_parsed, bbox_df)

cat("STEP 01 COMPLETE\n")
cat("=============================================================\n")
cat("\nData objects available:\n")
cat("  data_clean      - all detections with flags\n")
cat("  data_focal      - focal species only (with flags)\n")
cat("  data_density    - for density/EDD (calibration issues excluded, <=26m)\n")
cat("  data_speed      - for SBD (calibration deployments excluded,\n")
cat("                    bbox sequences excluded, group size > 1 excluded)\n")
cat("  data_activity   - for activity (all deployments)\n")
cat("  sequence_flags  - sequence-level quality summary\n")
cat("  species_summary - species counts and flag rates\n")
cat("\nNext: Run 02_sbd_speed.R\n\n")