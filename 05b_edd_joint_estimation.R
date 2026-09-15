# =============================================================================
# 05b_EDD_JOINT_ESTIMATION.R
# Chapter 1: Camera Trap Density Estimation
# =============================================================================
#
# PURPOSE: Joint multi-species EDD estimation with partial pooling.
#
#   Decomposes the detection scale parameter:
#     log(σ_ik) = α_i + δ_k
#
#   where:
#     α_i = deployment-level random effect (site-level "visibility")
#     δ_k = species offset (how detectable is species k relative to average?)
#
#   Scene depth and camera model enter as covariates on the hyperprior:
#     α_i ~ Normal(γ_0 + γ_1 × pixel_distance_i + γ_cam[cam_model_i], τ_α)
#
#   Camera model is categorical (reference-level coding). The first model
#   (alphabetically) is the reference (γ_cam[1] = 0); additional models get
#   free offsets. This captures hardware differences in PIR sensitivity,
#   focal length, and FOV that affect detection range independently of habitat.
#
#   Identifiability constraint: sum-to-zero on δ_k
#     → α_i represents detection range for an "average" species at site i
#     → δ_k is each species' deviation from that average
#
#   The key benefit: data-rich species (deer, bear) at a camera inform the
#   site effect α_i, which then improves EDD estimates for rare species at
#   that same camera via partial pooling.
#
# DIAGNOSTIC: Compares prior vs posterior on γ_1 (scene depth effect) to
#             address whether the covariate is informative.
#
# INPUT:   Parameter estimates from Steps 2-4, clean data from Step 1
# OUTPUT:  05b_edd_joint_results.rds
#
# =============================================================================

source("00_config.R")

library(tidyverse)
library(lubridate)
library(nimble)
library(MCMCvis)
library(coda)

cat("\n=============================================================\n")
cat("STEP 05b: JOINT MULTI-SPECIES EDD ESTIMATION\n")
cat("=============================================================\n\n")

# =============================================================================
# 1. EDD APPROXIMATION FUNCTION (unchanged from 05)
# =============================================================================

EDD_approx_logmix <- nimbleFunction(run = function(r         = double(0),
                                                   shape_d   = double(0),
                                                   shape_e   = double(0),
                                                   increment = double(0),
                                                   B         = double(0)) {
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
# 2. JOINT NIMBLE MODEL
# =============================================================================
#
# Observation level: m = 1:M  (one entry per species × deployment with N>0)
# Deployment level:  i = 1:I  (all deployments)
# Species level:     k = 1:K  (all species)
#
# Each observation m maps to a deployment (deploy_idx[m]) and species
# (species_idx[m]).
#
# The multinomial models the SHAPE of the distance distribution for each
# species-deployment combo. The Poisson/Binomial layers model total encounter
# counts. Together they estimate the detection function.
#
# NOTE: We separate detection estimation from density estimation. This model
#       estimates deployment×species EDD only. Density is computed downstream
#       using EDD posteriors + activity/speed/staying time from Steps 2-4.
# =============================================================================

joint_model_code <- nimbleCode({
  
  # ---- Species offsets: sum-to-zero constraint ----
  # K-1 free parameters; Kth is deterministic = -sum(others)
  for (k in 1:(K_minus_1)) {
    delta_free[k] ~ dnorm(0, tau_delta)
    delta[k] <- delta_free[k]
  }
  delta[K] <- -sum(delta_free[1:K_minus_1])
  
  # ---- Camera model offsets: reference-level coding ----
  # cam_model_idx[i] = 1 is the reference (gamma_cam[1] = 0)
  # Additional camera models get free offsets
  gamma_cam[1] <- 0                          # reference camera model
  for (c in 2:n_cam_models) {
    gamma_cam[c] ~ dnorm(0, 0.01)            # offset vs reference model
  }
  
  # ---- Deployment-level random effects ----
  # Scene depth + camera model enter as covariates on the mean
  for (i in 1:n_deployments) {
    alpha[i] ~ dnorm(gamma_0 + gamma_1 * pixel_distance[i] + gamma_cam[cam_model_idx[i]],
                     tau_alpha)
  }
  
  # ---- Observation-level: detection model ----
  # Each observation m is one species at one deployment (with detections > 0)
  for (m in 1:M) {
    # Detection scale from deployment effect + species offset
    log(sigma[m]) <- alpha[deploy_idx[m]] + delta[species_idx[m]]
    
    # Detection probabilities across distance bins
    for (j in 1:n_bins) {
      P[m, j]  <- exp(-(d_j[j]^2) / (2 * sigma[m]^2))
      Pi[m, j] <- psi[j] * P[m, j]
    }
    Pi_c[m, 1:n_bins] <- Pi[m, 1:n_bins] / sum(Pi[m, 1:n_bins])
    sum_Pi[m]          <- sum(Pi[m, 1:n_bins])
    
    # Multinomial on binned distances
    L[m, 1:n_bins] ~ dmulti(Pi_c[m, 1:n_bins], N_obs[m])
    
    # Encounter count model
    N_obs[m] ~ dbinom(sum_Pi[m], n_avail[m])
    n_avail[m] ~ dpois(mu[m])
    log(mu[m]) <- mu_0 + mu_deploy[deploy_idx[m]] + mu_species[species_idx[m]]
    
    # EDD for this species × deployment
    E[m] <- EDD_approx_logmix(sigma[m], shape_d, shape_e, increment = 0.01, B)
  }
  
  # ---- Nuisance: encounter rate intercepts ----
  # These just need to be flexible enough to not constrain the detection model
  mu_0 ~ dnorm(0, 0.01)
  for (i in 1:n_deployments) {
    mu_deploy[i] ~ dnorm(0, tau_mu_deploy)
  }
  for (k in 1:K) {
    mu_species[k] ~ dnorm(0, tau_mu_species)
  }
  tau_mu_deploy  ~ dgamma(0.01, 0.01)
  tau_mu_species ~ dgamma(0.01, 0.01)
  
  # ---- Priors ----
  gamma_0   ~ dnorm(0, 0.01)    # intercept for deployment effect (precision=0.01 → SD=10)
  gamma_1   ~ dnorm(0, 0.01)    # scene depth slope — KEY diagnostic param
  # gamma_cam[c] priors defined above (reference-level coding)
  tau_alpha  ~ dgamma(0.01, 0.01)   # precision of deployment effects
  tau_delta  ~ dgamma(0.01, 0.01)   # precision of species offsets
  shape_d   ~ dunif(0, 10)      # logistic shoulder shape
  shape_e   ~ dunif(0, 10)      # logistic shoulder midpoint
  
  # ---- Derived quantities ----
  sd_alpha <- 1 / sqrt(tau_alpha)
  sd_delta <- 1 / sqrt(tau_delta)
})

# =============================================================================
# 3. HELPER FUNCTIONS (carried from 05)
# =============================================================================

count_breaks_distance <- function(dist_breaks, dist) {
  result <- numeric(length(dist_breaks) - 1)
  for (i in seq_along(dist)) {
    w   <- which(dist_breaks > dist[i])
    bin <- if (length(w) == 0) length(dist_breaks) - 1 else min(w) - 1
    result[bin] <- result[bin] + 1
  }
  return(result)
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

cat("  Species with complete parameters:", nrow(species_params), "\n")
print(species_params)

# =============================================================================
# 5. LOAD TECH SPECS, DEPLOYMENTS, AND SCENE DEPTH
# =============================================================================

cat("\nLoading deployment metadata...\n")

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

# Scene depth covariate
DEPTH_FILE <- paste0(DATA_DIR, "1_raw/AI/depth_mmm.csv")
if (file.exists(DEPTH_FILE)) {
  depth_cal <- read.csv(DEPTH_FILE, stringsAsFactors = FALSE) %>%
    select(deployment_id, scene_depth = mean)
  cat("  DPT scene depth loaded:", nrow(depth_cal), "deployments\n")
} else {
  stop("depth_mmm.csv not found — joint model requires DPT scene depth")
}

survey_effort <- deployments %>%
  mutate(
    effort_secs = as.numeric(difftime(end_date, start_date, units = "secs")),
    effort_days = as.numeric(difftime(end_date, start_date, units = "days"))
  ) %>%
  filter(effort_secs > 0)

# =============================================================================
# 6. BUILD JOINT TRIGGER PICTURE DATASET (ALL SPECIES)
# =============================================================================
# Key change from 05: instead of looping over species and subsetting, we
# build one combined dataset and index by species and deployment.

cat("\nBuilding joint trigger picture dataset...\n")

data <- data %>%
  mutate(
    deployment_id     = toupper(str_trim(deployment_id_clean)),
    camera_name_clean = toupper(str_trim(camera_name))
  ) %>%
  left_join(deployments %>% select(deployment_id, start_date), by = "deployment_id") %>%
  filter(as.Date(timestamp_clean) != start_date)   # remove calibration frames

# Sort by frame number, collapse to 1-point-per-second, take trigger picture
data <- data %>%
  mutate(frame_number = as.numeric(str_extract(filename, "\\d+"))) %>%
  arrange(sequence_id_use, frame_number)

data_collapsed <- data %>%
  group_by(sequence_id_use, timestamp_clean) %>%
  slice(1) %>%
  ungroup()

trigger_all <- data_collapsed %>%
  group_by(sequence_id_use) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(distance = as.numeric(world_z)) %>%
  filter(common_name_clean %in% species_params$species)

cat("  Total trigger pictures across all species:", nrow(trigger_all), "\n")
cat("  Species represented:", n_distinct(trigger_all$common_name_clean), "\n")
cat("  Deployments represented:", n_distinct(trigger_all$deployment_id), "\n")

# Species and deployment lookup tables
species_list    <- sort(unique(trigger_all$common_name_clean))
K               <- length(species_list)
species_lookup  <- setNames(1:K, species_list)

deploy_list     <- sort(unique(trigger_all$deployment_id))
n_deployments_total <- length(deploy_list)
deploy_lookup   <- setNames(1:n_deployments_total, deploy_list)

cat("\n  Species (K =", K, "):\n")
for (k in 1:K) cat("    ", k, ":", species_list[k], "\n")
cat("\n  Deployments (I =", n_deployments_total, ")\n")

# =============================================================================
# 7. COMMON DISTANCE BINS
# =============================================================================
# One bin structure shared across all species. B = global max detection distance.

break_width_m  <- 3
B_global       <- ceiling(max(trigger_all$distance, na.rm = TRUE) / break_width_m) * break_width_m
dist_breaks_m  <- seq(0, B_global, by = break_width_m)
dist_midpoints <- dist_breaks_m[-length(dist_breaks_m)] + break_width_m / 2
n_bins         <- length(dist_midpoints)

total_area <- pi * B_global^2
psi <- sapply(1:(length(dist_breaks_m) - 1), function(i)
  (dist_breaks_m[i+1]^2 - dist_breaks_m[i]^2) / total_area)

cat("  Max detection distance (B):", B_global, "m\n")
cat("  Number of distance bins:", n_bins, "\n\n")

# =============================================================================
# 8. BUILD OBSERVATION-LEVEL DATA (LONG FORMAT)
# =============================================================================
# Each "observation" m = one species × deployment combination with N_ik > 0.
# This replaces the per-species loop + separate L matrix from 05.

cat("Building observation-level detection histograms...\n")

# Identify all species × deployment combos with detections
obs_combos <- trigger_all %>%
  group_by(deployment_id, common_name_clean) %>%
  summarise(
    N_detections = n(),
    .groups      = "drop"
  ) %>%
  filter(N_detections > 0) %>%
  mutate(
    deploy_idx  = deploy_lookup[deployment_id],
    species_idx = species_lookup[common_name_clean]
  ) %>%
  arrange(deploy_idx, species_idx)

M <- nrow(obs_combos)
cat("  Observation-level entries (M):", M, "\n")
cat("  (= species × deployment combos with detections)\n\n")

# Build the detection histogram matrix L[m, j]
L_matrix <- matrix(0, nrow = M, ncol = n_bins)

for (m in 1:M) {
  this_dep <- obs_combos$deployment_id[m]
  this_sp  <- obs_combos$common_name_clean[m]
  
  dists <- trigger_all %>%
    filter(deployment_id == this_dep, common_name_clean == this_sp) %>%
    pull(distance)
  
  L_matrix[m, ] <- count_breaks_distance(dist_breaks_m, dists)
}

# Verify row sums match N_detections
stopifnot(all(rowSums(L_matrix) == obs_combos$N_detections))
cat("  Detection histogram verified: row sums match N_detections\n")

# =============================================================================
# 9. DEPLOYMENT-LEVEL COVARIATE
# =============================================================================

deploy_data <- tibble(deployment_id = deploy_list) %>%
  left_join(depth_cal, by = "deployment_id") %>%
  left_join(
    trigger_all %>%
      select(deployment_id, camera_name) %>%
      mutate(camera_name_clean = toupper(str_trim(camera_name))) %>%
      distinct(deployment_id, .keep_all = TRUE),
    by = "deployment_id"
  ) %>%
  mutate(
    scene_depth    = ifelse(is.na(scene_depth),
                            median(scene_depth, na.rm = TRUE), scene_depth),
    pixel_distance = as.numeric(scale(log(scene_depth))),
    pixel_distance = ifelse(is.na(pixel_distance), 0, pixel_distance)
  )

# Camera model categorical index
# Collapse to brand level (e.g., "BROWNING ELITE HP5" → "BROWNING")
# This captures the meaningful hardware differences (PIR sensor, focal length)
# without overfitting to rare model variants with few deployments.
deploy_data <- deploy_data %>%
  mutate(
    camera_brand = case_when(
      str_detect(camera_name_clean, "BROWNING") ~ "BROWNING",
      str_detect(camera_name_clean, "RECONYX")  ~ "RECONYX",
      TRUE ~ NA_character_  # flag unrecognized names for imputation below
    )
  )

# If camera_brand is missing (unrecognized name or NA), assign to most common brand
n_missing_brand <- sum(is.na(deploy_data$camera_brand))
if (n_missing_brand > 0) {
  most_common_brand <- names(sort(table(deploy_data$camera_brand), decreasing = TRUE))[1]
  bad_names <- unique(deploy_data$camera_name_clean[is.na(deploy_data$camera_brand)])
  cat("  WARNING:", n_missing_brand, "deployments with unrecognized camera name(s):",
      paste(bad_names, collapse = ", "), "\n")
  cat("           Assigning to most common brand:", most_common_brand, "\n")
  deploy_data$camera_brand[is.na(deploy_data$camera_brand)] <- most_common_brand
}

cam_model_list  <- sort(unique(deploy_data$camera_brand))
n_cam_models    <- length(cam_model_list)
cam_model_lookup <- setNames(1:n_cam_models, cam_model_list)

deploy_data <- deploy_data %>%
  mutate(cam_model_idx = cam_model_lookup[camera_brand])

cat("  pixel_distance: scale(log(scene_depth)) from DPT\n")
cat("  Range:", round(min(deploy_data$pixel_distance), 2), "to",
    round(max(deploy_data$pixel_distance), 2), "\n")
cat("  Camera brands (n =", n_cam_models, "):\n")
for (c in 1:n_cam_models) {
  n_deps <- sum(deploy_data$cam_model_idx == c)
  ref_tag <- ifelse(c == 1, " [REFERENCE]", "")
  cat("    ", c, ":", cam_model_list[c], "(", n_deps, "deployments)", ref_tag, "\n")
}
cat("\n")

# =============================================================================
# 10. DETECTION SUMMARY BY SPECIES × DEPLOYMENT
# =============================================================================

cat("Detection counts per species:\n")
trigger_all %>%
  group_by(common_name_clean) %>%
  summarise(
    n_triggers    = n(),
    n_deployments = n_distinct(deployment_id),
    .groups       = "drop"
  ) %>%
  arrange(desc(n_triggers)) %>%
  print(n = 20)
cat("\n")

# =============================================================================
# 11. ASSEMBLE NIMBLE DATA AND CONSTANTS
# =============================================================================

cat("Assembling NIMBLE inputs...\n\n")

nimble_data <- list(
  L              = L_matrix,
  N_obs          = obs_combos$N_detections,
  pixel_distance = deploy_data$pixel_distance,
  d_j            = dist_midpoints,
  psi            = psi
)

nimble_constants <- list(
  M              = M,
  n_deployments  = n_deployments_total,
  K              = K,
  K_minus_1      = K - 1,
  n_cam_models   = n_cam_models,
  cam_model_idx  = deploy_data$cam_model_idx,
  n_bins         = n_bins,
  B              = B_global,
  deploy_idx     = obs_combos$deploy_idx,
  species_idx    = obs_combos$species_idx
)

# Initial values
# NOTE: gamma_cam[1] is deterministic (=0, reference level), so don't init it.
#       Only init gamma_cam[2:n_cam_models] if n_cam_models > 1.
gamma_cam_init <- rep(0, n_cam_models)
nimble_inits <- list(
  gamma_0        = log(10),          # ~10m baseline detection range
  gamma_1        = 0.2,              # slight positive scene depth effect
  gamma_cam      = gamma_cam_init,   # camera model offsets (idx 1 ignored by model)
  tau_alpha      = 1,
  tau_delta      = 1,
  tau_mu_deploy  = 1,
  tau_mu_species = 1,
  mu_0           = 3,
  shape_d        = 1,
  shape_e        = 1,
  alpha          = rep(log(10), n_deployments_total),
  delta_free     = rep(0, K - 1),
  mu_deploy      = rep(0, n_deployments_total),
  mu_species     = rep(0, K),
  n_avail        = obs_combos$N_detections * 2
)

# =============================================================================
# 12. BUILD AND RUN MODEL
# =============================================================================

cat("Building NIMBLE model...\n")

# Parameters to monitor
monitors <- c(
  "E",                          # deployment×species EDD (primary output)
  "alpha",                      # deployment effects
  "delta",                      # species offsets (including derived Kth)
  "gamma_0", "gamma_1",         # hyperprior params (gamma_1 = scene depth diagnostic)
  "gamma_cam",                  # camera model offsets
  "tau_alpha", "tau_delta",     # precisions
  "sd_alpha", "sd_delta",       # standard deviations (derived)
  "sigma",                      # detection scale per observation
  "shape_d", "shape_e"          # logistic shoulder params
)

tryCatch({
  
  model         <- nimbleModel(code = joint_model_code, data = nimble_data,
                               constants = nimble_constants, inits = nimble_inits)
  compiled_mod  <- compileNimble(model)
  mcmc_conf     <- configureMCMC(model, monitors = monitors)
  mcmc          <- buildMCMC(mcmc_conf)
  compiled_mcmc <- compileNimble(mcmc, project = model)
  
  # ---- Quick diagnostic run first ----
  cat("\n--- DIAGNOSTIC RUN (short) ---\n")
  samples_diag <- runMCMC(compiled_mcmc,
                          niter   = 2000,
                          nburnin = 500,
                          thin    = 5,
                          nchains = 2,
                          samplesAsCodaMCMC = TRUE)
  
  diag_summary <- MCMCsummary(samples_diag, params = c("gamma_0", "gamma_1",
                                                        "gamma_cam",
                                                        "sd_alpha", "sd_delta"))
  cat("\nDiagnostic run — key parameters:\n")
  print(diag_summary)
  cat("\nCheck: are gamma_1 and sd_alpha reasonable? If so, proceed to full run.\n")
  cat("If gamma_1 is stuck at 0 or sd_alpha is huge, something is wrong.\n\n")
  
  # ---- Full production run ----
  cat("--- FULL PRODUCTION RUN ---\n")
  samples_raw <- runMCMC(compiled_mcmc,
                         niter   = 20000,
                         nburnin = 5000,
                         thin    = 10,
                         nchains = 2,
                         samplesAsCodaMCMC = TRUE)
  
  posterior_sum <- MCMCsummary(samples_raw)
  
  # =============================================================================
  # 13. SCENE DEPTH COVARIATE DIAGNOSTIC
  # =============================================================================
  # Roland's request: "report how useful the prior was"
  # Compare prior vs posterior for gamma_1
  
  cat("\n##############################################\n")
  cat("SCENE DEPTH COVARIATE DIAGNOSTIC\n")
  cat("##############################################\n\n")
  
  gamma1_summary <- MCMCsummary(samples_raw, params = "gamma_1")
  cat("gamma_1 (scene depth effect on deployment detection range):\n")
  print(gamma1_summary)
  
  # Extract posterior samples for gamma_1
  gamma1_samples <- MCMCchains(samples_raw, params = "gamma_1")
  
  # Prior: gamma_1 ~ Normal(0, precision=0.01) → SD = 10
  prior_sd   <- sqrt(1/0.01)
  prior_mean <- 0
  
  cat("\n  Prior:     Normal(", prior_mean, ",", prior_sd, ")\n")
  cat("  Posterior: mean =", round(gamma1_summary[, "mean"], 4),
      ", SD =", round(gamma1_summary[, "sd"], 4), "\n")
  cat("  95% CI:   [", round(gamma1_summary[, "2.5%"], 4), ",",
      round(gamma1_summary[, "97.5%"], 4), "]\n\n")
  
  # Is the posterior substantially narrower than the prior?
  precision_gain <- prior_sd / gamma1_summary[, "sd"]
  cat("  Precision gain (prior SD / posterior SD):", round(precision_gain, 1), "x\n")
  
  # Does 95% CI exclude zero?
  ci_excludes_zero <- (gamma1_summary[, "2.5%"] > 0) | (gamma1_summary[, "97.5%"] < 0)
  cat("  95% CI excludes zero:", ci_excludes_zero, "\n")
  
  if (ci_excludes_zero) {
    direction <- ifelse(gamma1_summary[, "mean"] > 0, "positive", "negative")
    cat("  → Scene depth has a", direction, "effect on detection range.\n")
    cat("    Cameras with greater scene depth (more open habitat) have",
        ifelse(direction == "positive", "LARGER", "SMALLER"), "detection range.\n")
  } else {
    cat("  → Scene depth effect is not clearly distinguishable from zero.\n")
    cat("    The covariate may have limited influence on detection range.\n")
  }
  
  # Also report deployment-level variance (how much variation does scene depth explain?)
  cat("\n  Deployment random effect SD (sd_alpha):\n")
  print(MCMCsummary(samples_raw, params = "sd_alpha"))
  cat("  Species offset SD (sd_delta):\n")
  print(MCMCsummary(samples_raw, params = "sd_delta"))
  
  # =============================================================================
  # 13b. CAMERA MODEL COVARIATE DIAGNOSTIC
  # =============================================================================
  
  cat("\n##############################################\n")
  cat("CAMERA MODEL COVARIATE DIAGNOSTIC\n")
  cat("##############################################\n\n")
  
  gamma_cam_summary <- MCMCsummary(samples_raw, params = "gamma_cam")
  cat("gamma_cam (camera model effects on deployment detection range):\n")
  cat("  Reference model (index 1):", cam_model_list[1], "= 0 (fixed)\n\n")
  
  for (c in 1:n_cam_models) {
    cat("  Camera model", c, "(", cam_model_list[c], "):")
    if (c == 1) {
      cat(" 0.000 [reference]\n")
    } else {
      cat(" mean =", round(gamma_cam_summary[c, "mean"], 4),
          ", 95% CI [", round(gamma_cam_summary[c, "2.5%"], 4), ",",
          round(gamma_cam_summary[c, "97.5%"], 4), "]\n")
      cam_ci_excludes_zero <- (gamma_cam_summary[c, "2.5%"] > 0) |
                              (gamma_cam_summary[c, "97.5%"] < 0)
      if (cam_ci_excludes_zero) {
        cam_dir <- ifelse(gamma_cam_summary[c, "mean"] > 0, "LARGER", "SMALLER")
        cat("    → Significantly different from reference;", cam_dir,
            "detection range\n")
      } else {
        cat("    → Not significantly different from reference\n")
      }
    }
  }
  cat("\n  On the exp scale, a gamma_cam offset of X means detection scale is\n")
  cat("  exp(X) times the reference model.\n")
  
  # =============================================================================
  # 14. EXTRACT EDD RESULTS
  # =============================================================================
  
  cat("\n##############################################\n")
  cat("EXTRACTING EDD RESULTS\n")
  cat("##############################################\n\n")
  
  E_rows <- grep("^E\\[", rownames(posterior_sum))
  
  if (length(E_rows) == M) {
    E_values <- posterior_sum[E_rows, ]
    
    edd_by_obs <- obs_combos %>%
      mutate(
        EDD_mean  = E_values[, "mean"],
        EDD_sd    = E_values[, "sd"],
        EDD_2.5   = E_values[, "2.5%"],
        EDD_97.5  = E_values[, "97.5%"]
      )
    
    cat("EDD estimates for", M, "species × deployment combinations\n\n")
    
    # Summary by species
    edd_species_summary <- edd_by_obs %>%
      group_by(common_name_clean) %>%
      summarise(
        N_Deployments = n(),
        N_Detections  = sum(N_detections),
        Mean_EDD      = round(mean(EDD_mean), 2),
        SD_EDD        = round(sd(EDD_mean), 2),
        Min_EDD       = round(min(EDD_mean), 2),
        Max_EDD       = round(max(EDD_mean), 2),
        .groups = "drop"
      ) %>%
      arrange(desc(N_Detections))
    
    cat("=== EDD SUMMARY BY SPECIES (joint model) ===\n")
    print(edd_species_summary, n = 20)
    
    # Species offsets (delta)
    cat("\n=== SPECIES OFFSETS (delta_k) ===\n")
    cat("Positive = detected at greater range than average; Negative = shorter range\n\n")
    delta_rows <- grep("^delta\\[", rownames(posterior_sum))
    delta_summary <- posterior_sum[delta_rows, ]
    delta_df <- data.frame(
      species = species_list,
      delta_mean = round(delta_summary[, "mean"], 3),
      delta_sd   = round(delta_summary[, "sd"], 3),
      delta_2.5  = round(delta_summary[, "2.5%"], 3),
      delta_97.5 = round(delta_summary[, "97.5%"], 3)
    ) %>% arrange(desc(delta_mean))
    print(delta_df)
    
  } else {
    cat("WARNING: Expected", M, "E values but found", length(E_rows), "\n")
    edd_by_obs         <- NULL
    edd_species_summary <- NULL
    delta_df           <- NULL
  }
  
  # =============================================================================
  # 15. SAVE RESULTS
  # =============================================================================
  
  cat("\n\nSaving results...\n")
  
  output_rds <- paste0(OUTPUT_DIRS$processed, "05b_edd_joint_results.rds")
  saveRDS(list(
    # EDD estimates
    edd_by_obs          = edd_by_obs,
    edd_species_summary = edd_species_summary,
    
    # Model structure info
    species_list   = species_list,
    species_lookup = species_lookup,
    deploy_list    = deploy_list,
    deploy_lookup  = deploy_lookup,
    obs_combos     = obs_combos,
    
    # Camera model info
    cam_model_list   = cam_model_list,
    cam_model_lookup = cam_model_lookup,
    deploy_cam_model = deploy_data %>% select(deployment_id, camera_brand, cam_model_idx),
    
    # Species offsets
    delta_summary = delta_df,
    
    # Covariate diagnostics
    gamma1_summary   = gamma1_summary,
    gamma_cam_summary = gamma_cam_summary,
    precision_gain   = precision_gain,
    ci_excludes_zero = ci_excludes_zero,
    
    # Full posterior summary
    posterior_summary = posterior_sum,
    
    # Raw MCMC samples (for downstream density estimation)
    mcmc_samples = samples_raw
  ), output_rds)
  
  cat("Results saved to:", output_rds, "\n")
  
  # Also save CSVs for quick inspection
  if (!is.null(edd_by_obs)) {
    write.csv(edd_by_obs,
              paste0(OUTPUT_DIRS$model_output, "05b_EDD_joint_by_obs.csv"),
              row.names = FALSE)
    write.csv(edd_species_summary,
              paste0(OUTPUT_DIRS$model_output, "05b_EDD_joint_species_summary.csv"),
              row.names = FALSE)
    write.csv(delta_df,
              paste0(OUTPUT_DIRS$model_output, "05b_species_offsets_delta.csv"),
              row.names = FALSE)
  }
  
}, error = function(e) {
  cat("\n!!! MODEL ERROR !!!\n")
  cat("Error:", conditionMessage(e), "\n\n")
  cat("If this failed, check:\n")
  cat("  1. Are there species-deployment combos with very few detections?\n")
  cat("     (Multinomial needs N >= 1 in each included combo)\n")
  cat("  2. Is B_global too large for some species? (sparse bins)\n")
  cat("  3. Are initial values reasonable?\n")
  cat("  4. Try running the diagnostic (niter=2000) first.\n")
})

cat("\n=============================================================\n")
cat("STEP 05b COMPLETE\n")
cat("=============================================================\n\n")
