# =============================================================================
# CHAPTER 1 CONFIGURATION
# Camera Trap Density Estimation - Master Settings
# =============================================================================
# 
# PURPOSE: This file holds all settings, file paths, and thresholds in ONE place.
# Every other script starts with: source("00_config.R")
# Change something here → it updates everywhere automatically.
# =============================================================================

library(dplyr)
library(lubridate)
library(readr)
library(stringr)

# =============================================================================
# FILE PATHS
# =============================================================================

# Data directory (up 3 levels, then into dissertation_main)
DATA_DIR <- "../../../dissertation_main/01_data/"
# Results directory (up 2 levels, then into 03_results)
RESULTS_DIR <- "../../03_results/"
# Code directory (current location)
CODE_DIR <- "./"

# AI coordinate data files (from Mohammad's pipeline) - UPDATED FEB 2026
COORD_FILES <- list(
  browning_elite = paste0(DATA_DIR, "1_raw/AI/results-Browning Recon Force Elite-full-merged-normalized.csv"),
  reconyx = paste0(DATA_DIR, "1_raw/AI/results-Reconyx Hyperfire 1-full-merged-normalized.csv")
)

# Depth calibration data
DEPTH_FILE <- paste0(DATA_DIR, "1_raw/AI/depth_mmm.csv")

# Wildlife Insights files
WI_FILES <- list(
  images      = paste0(DATA_DIR, "1_raw/Camera/images_2006449.csv"),
  sequences   = paste0(DATA_DIR, "1_raw/Camera/sequences.csv"),
  deployments = paste0(DATA_DIR, "1_raw/Camera/deployments_allmysites.csv")
)

# Output directories (all within ch1/)
OUTPUT_DIRS <- list(
  processed    = paste0(RESULTS_DIR, "processed/"),
  model_output = paste0(RESULTS_DIR, "model_output/"),
  figures      = paste0(RESULTS_DIR, "final_figures/"),
  tables       = paste0(RESULTS_DIR, "tables/")
)

# =============================================================================
# FOCAL SPECIES
# =============================================================================

FOCAL_SPECIES <- c(
  # Eastern species
  "White-tailed Deer", "Wild Turkey", "Northern Raccoon", "Eastern Gray Squirrel",
  "Red Squirrel", "Coyote", "Virginia Opossum", "Eastern Fox Squirrel",
  "Eastern Cottontail", "Eastern Chipmunk", "Snowshoe Hare",
  "Bobcat", "Red Fox", "Grey Fox", "Fisher", "American Black Bear",
  # Western species  
  "Mule Deer", "Elk", "Puma", "Black-tailed Deer", "Pronghorn", "Moose"
)

# =============================================================================
# CAMERA SPECIFICATIONS
# =============================================================================

CAMERA_SPECS <- list(
  frame_rate            = 1.4,   # frames per second
  trigger_speed         = 0.5,   # seconds from detection to first image
  detection_width       = 8.1,   # meters - FOV width at typical distance
  typical_detection_dist = 10    # meters
)

# =============================================================================
# ANALYSIS THRESHOLDS
# =============================================================================

THRESHOLDS <- list(
  # Distance filters
  min_distance     = 1,    # meters - flag detections closer than this
  max_distance     = 26,   # meters - UPDATE TO 50 when new calibration data arrives
  
  # Bounding box filter
  max_bbox_pct     = 80,   # % of frame - flag if bbox larger than this
  
  # Coordinate anchoring (for SBD speed estimation)
  anchoring_threshold = 0.10,  # meters - movements smaller than this set to 0
  
  # Speed filters
  max_frame_speed    = 20,  # m/s - biologically impossible above this
  max_sequence_speed = 50,  # m/s - flag sequences with impossibly high avg speed
  
  # Minimum sample sizes
  min_activity_n  = 15,  # minimum detections for activity estimation
  min_sbd_n       = 10,  # minimum sequences for SBD estimation
  min_density_n   = 5    # minimum deployments for density estimation
)

# =============================================================================
# DEPLOYMENT EXCLUSIONS
# =============================================================================

# Deployments with calibration issues (max dist <=3m AND range <0.5m).
# AI depth coordinates are clamped/unreliable.
# Excluded from: speed (01_data_for_speed.csv) and EDD (05_EDD_estimation.R)
# Kept for: activity (timestamps only) and staying time
EXCLUDE_DEPLOYMENTS <- c("A24", "A26", "B13", "C16", "C3", "D14", "D2", "F1")

# Deployments with good imagery but bad AI coordinate outputs.
# Images are fine so timestamps are valid -- keep for activity and staying time.
# Excluded from: speed (01_data_for_speed.csv) and density/EDD (01_data_for_density.csv)
EXCLUDE_DEPLOYMENTS_COORDS <- c("B18")

# Combined exclusion for any analysis using AI coordinates
EXCLUDE_DEPLOYMENTS_ALL_COORDS <- c(EXCLUDE_DEPLOYMENTS, EXCLUDE_DEPLOYMENTS_COORDS)

# =============================================================================
# SBD CONFIGURATION (per Rowcliffe 2016 minimal filtering)
# =============================================================================

SBD_CONFIG <- list(
  bootstrap_reps = 1000,
  min_frames     = 2,
  min_sequences  = 10,
  max_tortuosity = 3
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Check if required files exist
check_data_files <- function() {
  cat("Checking data files...\n")
  for (name in names(COORD_FILES)) {
    if (file.exists(COORD_FILES[[name]])) {
      cat("  \u2713", name, "\n")
    } else {
      cat("  \u2717", name, "- FILE NOT FOUND:", COORD_FILES[[name]], "\n")
    }
  }
  for (name in names(WI_FILES)) {
    if (file.exists(WI_FILES[[name]])) {
      cat("  \u2713", name, "\n")
    } else {
      cat("  \u2717", name, "- FILE NOT FOUND:", WI_FILES[[name]], "\n")
    }
  }
}

# Create output directories if they don't exist
setup_output_dirs <- function() {
  for (dir in OUTPUT_DIRS) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created:", dir, "\n")
    }
  }
}

# =============================================================================
# STARTUP MESSAGE
# =============================================================================

cat("\n")
cat("=============================================================\n")
cat("CHAPTER 1 CONFIGURATION LOADED\n")
cat("=============================================================\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Focal species:", length(FOCAL_SPECIES), "\n")
cat("Max distance threshold:", THRESHOLDS$max_distance, "m\n")
cat("Anchoring threshold:", THRESHOLDS$anchoring_threshold, "m\n")
cat("Deployments excluded (calibration):",
    paste(EXCLUDE_DEPLOYMENTS, collapse = ", "), "\n")
cat("Deployments excluded (bad AI coords):",
    paste(EXCLUDE_DEPLOYMENTS_COORDS, collapse = ", "), "\n")
cat("=============================================================\n\n")

# Uncomment these to run checks on load:
check_data_files()
setup_output_dirs()