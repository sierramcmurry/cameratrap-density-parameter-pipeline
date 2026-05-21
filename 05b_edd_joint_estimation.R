# =============================================================================
# 05_EDD_ESTIMATION.R
# Chapter 1: Camera Trap Density Estimation
# =============================================================================
#
# PURPOSE: Estimate Effective Detection Distance (EDD) using NIMBLE model.
#          Uses SBD speeds, activity, and staying time from previous steps.
#
# KEY DECISIONS:
#
#   1. TRIGGER PICTURES FOR DETECTION HISTORY
#      Distance sampling requires one independent detection distance per
#      detection event. We use the trigger picture (first timestamp of each
#      sequence) — the moment the animal first entered the viewshed.
#      Using all burst frames would inflate detection counts and bias the
#      distance histogram toward longer distances.
#
#   2. DPT SCENE DEPTH AS PIXEL_DISTANCE COVARIATE
#      pixel_distance = scale(log(scene_depth)) from depth_mmm.csv.
#      This replaces mean detection distance, which was circular (using
#      detection distance to predict detection range r -> EDD). DPT scene
#      depth is a genuine environmental covariate reflecting habitat openness
#      that independently predicts camera detection range.
#
# INPUT:   Parameter estimates from Steps 2-4, clean data from Step 1
# OUTPUT:  05_edd_results.rds
#
# =============================================================================

source("00_config.R")

library(tidyverse)
library(lubridate)
library(nimble)
library(MCMCvis)
library(coda)

cat("\n=============================================================\n")
cat("STEP 05: EDD ESTIMATION\n")
cat("=============================================================\n\n")

# =============================================================================
# 1. EDD APPROXIMATION FUNCTION
# =============================================================================

EDD_approx_logmix <- nimbleFunction(run = function(r       = double(0),
                                                   shape_d = double(0),
                                                   shape_e = double(0),
                                                   increment = double(0),
                                                   B       = double(0)) {
  points        <- (1:(B/increment)) * increment
  area          <- (points + increment/2)^2 - (points - increment/2)^2
  relative_area <- area / sum(area)
  detprob       <- exp(-(points^2) / (2 * r^2)) *
    (1 / (1 + exp(shape_d * (shape_e - points))))
  E_captured    <- detprob * relative_area
  EDD           <- B * sqrt(sum(E_captured))
  return(EDD)
  returnType(double(0))
})

# =============================================================================
# 2. NIMBLE MODEL CODE
# =============================================================================

model_code <- nimbleCode({
  
  for (i in 1:n_deployments) {
    L[i, 1:n_bins] ~ dmulti(Pi_ij_c[i, 1:n_bins], N_i[i])
    N_i[i] ~ dbinom(sum_Pi[i], n[i])
    n[i]   ~ dpois(mu[i])
    
    log(mu[i]) <- B_0_mu + beta_mu * pixel_distance[i]
    log(r[i])  <- B_0_r  + beta_r  * pixel_distance[i]
    
    E[i]     <- EDD_approx_logmix(r[i], shape_d, shape_e, increment = 0.01, B)
    A_est[i] <- 0.5 * E[i]^2 * viewshed_angle[i]
    
    for (j in 1:n_bins) {
      P_ij[i, j]  <- exp(-(d_j[j]^2) / (2 * r[i]^2))
      Pi_ij[i, j] <- psi[j] * P_ij[i, j]
    }
    Pi_ij_c[i, 1:n_bins] <- Pi_ij[i, 1:n_bins] / sum(Pi_ij[i, 1:n_bins])
    sum_Pi[i]             <- sum(Pi_ij[i, 1:n_bins])
  }
  
  for (i in 1:n_deployments) {
    y_plot[i] ~ dpois(mu_plot[i])
    log(mu_plot[i]) <- log(D_plot[i]) + log(A_est[i]) +
      log(p_active * T[i]) - log(s_mean)
    
    y_move[i] ~ dpois(mu_move[i])
    log(mu_move[i]) <- log(D_move[i]) +
      log(E[i] * (2 + viewshed_angle[i]) * p_active * T[i]) +
      log(v) - log(3.14)
    
    D_plot[i] <- exp(B_0_D_p + beta_d_p * pixel_distance[i]) / 1e6
    D_move[i] <- exp(B_0_D_m + beta_d_m * pixel_distance[i]) / 1e6
  }
  
  D_plot_mean <- exp(B_0_D_p + beta_d_p * mean_pixel_dist)
  D_move_mean <- exp(B_0_D_m + beta_d_m * mean_pixel_dist)
  
  # Priors
  B_0_mu   ~ dnorm(0, 0.01)
  beta_mu  ~ dnorm(0, 0.01)
  B_0_r    ~ dnorm(0, 0.01)
  beta_r   ~ dnorm(0, 0.01)
  B_0_D_p  ~ dnorm(0, 0.01)
  B_0_D_m  ~ dnorm(0, 0.01)
  beta_d_p ~ dnorm(0, 0.01)
  beta_d_m ~ dnorm(0, 0.01)
  shape_d  ~ dunif(0, 10)
  shape_e  ~ dunif(0, 10)
})

# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

count_breaks_distance <- function(dist_breaks, dist) {
  result <- numeric(length(dist_breaks) - 1)
  for (i in seq_along(dist)) {
    w <- which(dist_breaks > dist[i])
    bin <- if (length(w) == 0) length(dist_breaks) - 1 else min(w) - 1
    result[bin] <- result[bin] + 1
  }
  return(result)
}

construct_dethist_array <- function(datlist, dist_breaks) {
  y_array <- matrix(0, nrow = nrow(datlist$deployments_df),
                    ncol = length(dist_breaks) - 1)
  for (i in 1:nrow(datlist$deployments_df)) {
    this_dets    <- datlist$individual_dets_df %>%
      filter(deployment_id == datlist$deployments_df$deployment_id[i])
    y_array[i, ] <- count_breaks_distance(dist_breaks, this_dets$distance)
  }
  return(y_array)
}

parse_angle <- function(angle_str) {
  if (is.na(angle_str) || angle_str == "" || angle_str == "not available") return(30)
  num <- as.numeric(str_extract(angle_str, "[0-9.]+"))
  if (is.na(num)) return(30)
  return(num)
}

# =============================================================================
# 4. LOAD DATA AND PARAMETERS
# =============================================================================

cat("Loading data and parameter estimates...\n")

data <- read.csv(paste0(OUTPUT_DIRS$processed, "01_data_for_density.csv"),
                 stringsAsFactors = FALSE)
data$timestamp_clean <- as.POSIXct(data$timestamp_clean, tz = "UTC")

sbd_results      <- readRDS(paste0(OUTPUT_DIRS$processed, "02_sbd_results.rds"))
activity_results <- readRDS(paste0(OUTPUT_DIRS$processed, "03_activity_results.rds"))
staying_results  <- readRDS(paste0(OUTPUT_DIRS$processed, "04_staying_time_results.rds"))

species_params <- sbd_results$sbd_results %>%
  select(species = common_name_clean, sbd_speed = sbd_mean, n_speed = n) %>%
  left_join(
    activity_results$activity_by_species %>%
      select(species = common_name_clean, activity_level = activity, n_activity = n),
    by = "species"
  ) %>%
  left_join(
    staying_results$staying_by_species %>%
      select(species = common_name_clean, staying_time = mean_staying_time,
             n_staying = sample_size),
    by = "species"
  ) %>%
  filter(!is.na(sbd_speed) & !is.na(activity_level) & !is.na(staying_time))

cat("  - Species with complete parameters:", nrow(species_params), "\n")
print(species_params)
cat("\n")

# =============================================================================
# 5. LOAD TECH SPECS AND DEPLOYMENTS
# =============================================================================

cat("Loading tech specs and deployment data...\n")

tech_specs_sub <- read.csv(paste0(DATA_DIR, "1_raw/AI/tech_specs.csv"),
                           stringsAsFactors = FALSE) %>%
  select(camera_name, angle) %>%
  mutate(
    angle_numeric     = sapply(angle, parse_angle),
    camera_name_clean = toupper(str_trim(camera_name))
  )

deployments <- read.csv(WI_FILES$deployments, stringsAsFactors = FALSE) %>%
  select(deployment_id, start_date, end_date, camera_name) %>%
  mutate(
    deployment_id = toupper(str_trim(deployment_id)),
    start_date    = as.Date(start_date),
    end_date      = as.Date(end_date)
  ) %>%
  distinct(deployment_id, .keep_all = TRUE)

cat("  - Tech specs:", nrow(tech_specs_sub), "camera types\n")
cat("  - Deployments:", nrow(deployments), "\n")

# =============================================================================
# 5b. DPT SCENE DEPTH COVARIATE
# =============================================================================

DEPTH_FILE <- paste0(DATA_DIR, "1_raw/AI/depth_mmm.csv")
if (file.exists(DEPTH_FILE)) {
  depth_cal <- read.csv(DEPTH_FILE, stringsAsFactors = FALSE) %>%
    select(deployment_id, scene_depth = mean)
  cat("  - DPT scene depth loaded:", nrow(depth_cal), "deployments\n\n")
} else {
  depth_cal <- NULL
  cat("  WARNING: depth_mmm.csv not found — falling back to mean detection\n")
  cat("  distance (circular — flag in methods if used)\n\n")
}

# =============================================================================
# 6. SURVEY EFFORT
# =============================================================================

survey_effort <- deployments %>%
  mutate(
    effort_secs = as.numeric(difftime(end_date, start_date, units = "secs")),
    effort_days = as.numeric(difftime(end_date, start_date, units = "days"))
  ) %>%
  filter(effort_secs > 0)

cat("Survey effort:", nrow(survey_effort), "deployments with valid effort\n\n")

# =============================================================================
# 7. EXTRACT TRIGGER PICTURES
# =============================================================================
# Use the first frame of each sequence (trigger picture) as the detection
# distance for the detection history. Avoids inflating counts from burst
# photography and matches the independent-detection assumption of distance
# sampling.

cat("Extracting trigger pictures...\n")

data <- data %>%
  mutate(
    deployment_id     = toupper(str_trim(deployment_id_clean)),
    camera_name_clean = toupper(str_trim(camera_name))
  ) %>%
  left_join(deployments %>% select(deployment_id, start_date), by = "deployment_id") %>%
  filter(as.Date(timestamp_clean) != start_date)   # remove calibration frames

cat("  - Records after removing start date:", nrow(data), "\n")

# Viewshed angles per deployment
viewshed_angles <- data %>%
  select(deployment_id, camera_name_clean) %>%
  distinct() %>%
  left_join(tech_specs_sub %>% select(camera_name_clean, angle_numeric),
            by = "camera_name_clean") %>%
  mutate(
    viewshed_angle_deg = ifelse(is.na(angle_numeric), 30, angle_numeric),
    viewshed_angle     = (viewshed_angle_deg / 360) * (2 * pi)
  )

# Sort by frame number — consistent with 02_sbd_speed.R
data <- data %>%
  mutate(frame_number = as.numeric(str_extract(filename, "\\d+"))) %>%
  arrange(sequence_id_use, frame_number)

# 1-point-per-second: first frame per unique timestamp
data_collapsed <- data %>%
  group_by(sequence_id_use, timestamp_clean) %>%
  slice(1) %>%
  ungroup()

# Trigger picture = first timestamp of each sequence
trigger_pictures <- data_collapsed %>%
  group_by(sequence_id_use) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(distance = as.numeric(world_z))

cat("  - Total trigger pictures:", nrow(trigger_pictures), "\n\n")
stopifnot(n_distinct(trigger_pictures$sequence_id_use) == nrow(trigger_pictures))

# =============================================================================
# 8. RUN EDD ANALYSIS FOR ALL SPECIES
# =============================================================================

all_species           <- species_params$species
all_edd_by_deployment <- list()

cat("Starting EDD analysis for", length(all_species), "species...\n\n")

for (SPECIES_TO_ANALYZE in all_species) {
  
  cat("\n##############################################\n")
  cat("ANALYZING:", SPECIES_TO_ANALYZE, "\n")
  cat("##############################################\n\n")
  
  params   <- species_params[species_params$species == SPECIES_TO_ANALYZE, ]
  s_mean   <- params$staying_time
  p_active <- params$activity_level
  v_sbd    <- params$sbd_speed
  
  cat("SBD Speed:    ", round(v_sbd, 3), "m/s\n")
  cat("Activity:     ", round(p_active, 3), "\n")
  cat("Staying Time: ", round(s_mean, 1), "s\n")
  
  data_species <- trigger_pictures %>%
    filter(common_name_clean == SPECIES_TO_ANALYZE)
  
  if (nrow(data_species) == 0) { cat("No data — SKIPPING\n"); next }
  cat("Trigger detections:", nrow(data_species), "\n")
  
  datlist <- list(
    deployments_df     = data_species %>% select(deployment_id) %>% distinct(),
    individual_dets_df = data_species %>% select(deployment_id, distance)
  )
  
  n_deployments <- nrow(datlist$deployments_df)
  cat("Deployments:", n_deployments, "\n")
  if (n_deployments < 3) { cat("Too few deployments — SKIPPING\n"); next }
  
  # Distance bins
  break_width_m  <- 3
  max_distance_m <- ceiling(max(data_species$distance, na.rm = TRUE) /
                              break_width_m) * break_width_m
  dist_breaks_m  <- seq(0, max_distance_m, by = break_width_m)
  dist_midpoints <- dist_breaks_m[-length(dist_breaks_m)] + break_width_m / 2
  
  total_area <- pi * max_distance_m^2
  psi <- sapply(1:(length(dist_breaks_m) - 1), function(i)
    (dist_breaks_m[i+1]^2 - dist_breaks_m[i]^2) / total_area)
  
  dethist_array <- construct_dethist_array(datlist, dist_breaks_m)
  L_matrix      <- matrix(0, nrow = nrow(dethist_array), ncol = length(dist_midpoints))
  for (i in seq_len(nrow(L_matrix))) L_matrix[i, ] <- as.vector(dethist_array[i, ])
  
  # Build detection_data with pixel_distance covariate
  # DPT scene depth preferred; falls back to mean_det_dist if unavailable
  detection_data <- datlist$deployments_df %>%
    left_join(survey_effort, by = "deployment_id") %>%
    left_join(
      data_species %>%
        group_by(deployment_id) %>%
        summarise(mean_det_dist = mean(distance, na.rm = TRUE), .groups = "drop"),
      by = "deployment_id"
    )
  
  if (!is.null(depth_cal)) {
    detection_data <- detection_data %>%
      left_join(depth_cal, by = "deployment_id") %>%
      mutate(
        scene_depth    = ifelse(is.na(scene_depth),
                                median(scene_depth, na.rm = TRUE), scene_depth),
        pixel_distance = as.numeric(scale(log(scene_depth))),
        pixel_distance = ifelse(is.na(pixel_distance), 0, pixel_distance)
      )
    cat("  pixel_distance: scale(log(scene_depth)) from DPT\n")
  } else {
    detection_data <- detection_data %>%
      mutate(
        pixel_distance = as.numeric(scale(log(mean_det_dist))),
        pixel_distance = ifelse(is.na(pixel_distance), 0, pixel_distance)
      )
    cat("  NOTE: pixel_distance using mean_det_dist fallback (circular)\n")
  }
  
  if (nrow(detection_data) != nrow(L_matrix)) { cat("Row mismatch — SKIPPING\n"); next }
  
  detection_data$N_i <- rowSums(L_matrix)
  
  density_data_aligned <- detection_data %>%
    select(deployment_id) %>%
    left_join(
      data_species %>%
        group_by(deployment_id) %>%
        summarise(y_plot = n_distinct(sequence_id_use),
                  y_move = n_distinct(sequence_id_use),
                  .groups = "drop"),
      by = "deployment_id"
    ) %>%
    mutate(y_plot = ifelse(is.na(y_plot), 0, y_plot),
           y_move = ifelse(is.na(y_move), 0, y_move))
  
  viewshed_aligned <- detection_data %>%
    left_join(viewshed_angles, by = "deployment_id") %>%
    mutate(viewshed_angle = ifelse(is.na(viewshed_angle),
                                   (30/360) * 2 * pi, viewshed_angle)) %>%
    pull(viewshed_angle)
  
  nimble_data <- list(
    L              = dethist_array,
    N_i            = detection_data$N_i,
    pixel_distance = detection_data$pixel_distance,
    d_j            = dist_midpoints,
    y_plot         = density_data_aligned$y_plot,
    y_move         = density_data_aligned$y_move,
    T              = detection_data$effort_secs,
    psi            = psi,
    viewshed_angle = viewshed_aligned
  )
  
  nimble_constants <- list(
    n_deployments   = nrow(detection_data),
    n_bins          = length(dist_midpoints),
    B               = max_distance_m,
    v               = v_sbd,
    mean_pixel_dist = mean(detection_data$pixel_distance),
    s_mean          = s_mean,
    p_active        = p_active
  )
  
  nimble_inits <- list(
    B_0_mu   = 5,
    B_0_r    = log(10),
    beta_mu  = 0.2,
    beta_r   = 0.5,
    B_0_D_p  = log(10),
    beta_d_p = 0,
    B_0_D_m  = log(10),
    beta_d_m = 0.2,
    shape_d  = 1,
    shape_e  = 1,
    n        = detection_data$N_i * 2
  )
  
  tryCatch({
    model         <- nimbleModel(code = model_code, data = nimble_data,
                                 constants = nimble_constants, inits = nimble_inits)
    compiled_mod  <- compileNimble(model)
    mcmc_conf     <- configureMCMC(model, monitors = c("E"))
    mcmc          <- buildMCMC(mcmc_conf)
    compiled_mcmc <- compileNimble(mcmc, project = model)
    
    cat("Running MCMC...\n")
    samples_raw <- runMCMC(compiled_mcmc,
                           niter   = 20000,
                           nburnin = 5000,
                           thin    = 10,
                           nchains = 2,
                           samplesAsCodaMCMC = TRUE)
    
    posterior_sum <- MCMCsummary(samples_raw)
    E_rows        <- grep("^E\\[", rownames(posterior_sum))
    
    if (length(E_rows) > 0) {
      E_values <- posterior_sum[E_rows, ]
      
      edd_by_deployment <- data.frame(
        deployment_id  = detection_data$deployment_id,
        Species        = SPECIES_TO_ANALYZE,
        SBD_Speed      = v_sbd,
        Activity       = p_active,
        Staying_Time   = s_mean,
        EDD_mean       = E_values[, "mean"],
        EDD_sd         = E_values[, "sd"],
        EDD_2.5        = E_values[, "2.5%"],
        EDD_97.5       = E_values[, "97.5%"],
        N_detections   = detection_data$N_i,
        pixel_distance = detection_data$pixel_distance,
        mean_det_dist  = detection_data$mean_det_dist,
        stringsAsFactors = FALSE
      )
      
      cat("\nEDD by deployment:", nrow(edd_by_deployment), "deployments\n")
      cat("EDD range:", round(min(edd_by_deployment$EDD_mean), 2), "-",
          round(max(edd_by_deployment$EDD_mean), 2), "m\n")
      cat("Mean EDD:", round(mean(E_values[, "mean"], na.rm = TRUE), 2), "m\n")
      
      all_edd_by_deployment[[SPECIES_TO_ANALYZE]] <- edd_by_deployment
      
    } else {
      cat("WARNING: No EDD values extracted — check model convergence\n")
    }
    
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

# =============================================================================
# 9. COMPILE AND SAVE RESULTS
# =============================================================================

cat("\n\n##############################################\n")
cat("COMPILING EDD RESULTS\n")
cat("##############################################\n\n")

if (length(all_edd_by_deployment) > 0) {
  
  edd_all_species <- do.call(rbind, all_edd_by_deployment)
  
  edd_summary <- edd_all_species %>%
    group_by(Species) %>%
    summarise(
      N_Deployments = n(),
      N_Detections  = sum(N_detections),
      SBD_Speed     = first(SBD_Speed),
      Activity      = first(Activity),
      Staying_Time  = first(Staying_Time),
      Mean_EDD      = round(mean(EDD_mean), 2),
      SD_EDD        = round(sd(EDD_mean),   2),
      Min_EDD       = round(min(EDD_mean),  2),
      Max_EDD       = round(max(EDD_mean),  2),
      .groups = "drop"
    )
  
  cat("=== EDD SUMMARY BY SPECIES ===\n")
  print(edd_summary)
  
  write.csv(edd_all_species,
            paste0(OUTPUT_DIRS$model_output, "05_EDD_by_deployment_all_species.csv"),
            row.names = FALSE)
  write.csv(edd_summary,
            paste0(OUTPUT_DIRS$model_output, "05_EDD_summary_by_species.csv"),
            row.names = FALSE)
  
  output_rds <- paste0(OUTPUT_DIRS$processed, "05_edd_results.rds")
  saveRDS(list(edd_all_species = edd_all_species,
               edd_summary     = edd_summary),
          output_rds)
  cat("\nResults saved to:", output_rds, "\n")
  
} else {
  cat("No EDD results — all species failed or were skipped\n")
}

cat("\n=============================================================\n")
cat("STEP 05 COMPLETE\n")
cat("=============================================================\n\n")