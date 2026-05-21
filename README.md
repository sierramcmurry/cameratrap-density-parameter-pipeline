# Automated Parameter Estimation for Camera Trap Density Models

[![DOI](https://zenodo.org/badge/DOI/XXXX.svg)](https://doi.org/XXXX)

This repository contains the R code and example data for the automated parameter estimation pipeline described in:

> McMurry, S., Goldstein, B., Alyetam, M., & Kays, R. (in review). Automated Parameter Estimation for Camera Trap Density Models Using Computer Vision-Enhanced Distance Sampling. *Ecological Monographs*.

## Overview

Camera trap density models (REM, CTDS, REST, STE, TTE) require four parameters that are traditionally estimated through labor-intensive manual measurements:

1. **Effective detection distance (EDD)** — how far a camera reliably detects animals
2. **Movement speed** — how fast animals travel through detection zones
3. **Activity level** — what proportion of the day animals are active
4. **Staying time** — how long animals remain within the detection zone

This pipeline automates extraction of all four parameters directly from camera trap imagery using computer vision-derived coordinates, eliminating the need for manual annotation. The AI depth estimation pipeline (MegaDetector, SAM, DPT) that generates the input coordinates is available at [Alyetama/distance-estimation](https://github.com/Alyetama/distance-estimation). The key methodological contribution is a **joint multi-species hierarchical detection model** that estimates deployment-specific EDDs while borrowing strength across species through partial pooling.

## Pipeline Overview

The pipeline assumes you have already processed your camera trap images through the [AI depth estimation pipeline](https://github.com/Alyetama/distance-estimation) (MegaDetector + SAM + DPT) to generate per-frame world coordinates. Starting from those AI outputs and Wildlife Insights metadata, the R scripts estimate all density parameters:

```
AI coordinate outputs + WI metadata
        │
        ▼
  01_data_preparation.R     →  Clean, filter, join data
        │
        ▼
  02_sbd_speed.R            →  Movement speed (SBD framework)
        │
        ▼
  03_activity.R             →  Activity level (circular KDE)
        │
        ▼
  04_staying_time.R         →  Staying time (bounding box interpolation)
        │
        ▼
  05b_edd_joint_estimation.R →  Effective detection distance (joint hierarchical model)
```

## Repository Structure

```
├── README.md
├── LICENSE
├── R/
│   ├── 00_config.R                     # Configuration and file paths
│   ├── 01_data_preparation.R           # Data loading, joining, quality filtering
│   ├── 02_sbd_speed.R                  # SBD speed estimation with coordinate anchoring
│   ├── 03_activity.R                   # Activity level via circular kernel density
│   ├── 04_staying_time.R               # Staying time via bounding box interpolation
│   └── 05b_edd_joint_estimation.R      # Joint multi-species hierarchical EDD model (NIMBLE)
├── vignette/
│   ├── vignette_parameter_estimation.Rmd   # Worked example (start here)
│   └── example_vignette_data.csv           # Example dataset (28 deployments, 9 species)
└── data/
    └── README_data.md                  # Description of input data format
```

## Quick Start

### 1. Run the vignette

The fastest way to understand the pipeline is the worked example:

1. Clone this repository
2. Open `vignette/vignette_parameter_estimation.Rmd` in RStudio
3. Knit the document

The vignette walks through every step of the pipeline using an included example dataset of 28 camera deployments from montane forests in Montana, USA.

### 2. Run on your own data

To use the pipeline with your own camera trap data:

1. Process images through the [AI depth estimation pipeline](https://github.com/Alyetama/distance-estimation) (MegaDetector → SAM → DPT) to generate per-frame world coordinates
2. Export Wildlife Insights metadata (images, sequences, deployments)
3. Update paths in `R/00_config.R`
4. Run scripts 01 through 05 in sequence

## Input Data Requirements

The pipeline requires two main inputs:

### AI coordinate output (CSV)
One row per detected animal per frame, with columns:

| Column | Description |
|--------|------------|
| `deployment_id_clean` | Camera deployment identifier |
| `sequence_id_use` | Detection sequence identifier |
| `common_name_clean` | Species common name |
| `timestamp_clean` | Detection timestamp (POSIXct) |
| `world_x` | Lateral position relative to camera (meters) |
| `world_z` | Forward distance from camera (meters) |
| `x1, y1, x2, y2` | Bounding box pixel coordinates |
| `camera_name` | Camera model name |
| `detection_confidence` | AI detection confidence score |

### Deployment metadata (CSV)
One row per camera deployment:

| Column | Description |
|--------|------------|
| `deployment_id` | Camera deployment identifier |
| `start_date` | Deployment start (POSIXct) |
| `end_date` | Deployment end (POSIXct) |
| `camera_name` | Camera model |
| `latitude, longitude` | Deployment coordinates |

## Key Methods

### Size-Biased Distribution (SBD) Speed Estimation
Camera traps over-sample faster animals because they cross the detection zone more frequently. The SBD framework (Rowcliffe et al. 2016) corrects for this using the harmonic mean of observed speeds. We implement species-specific coordinate anchoring thresholds (determined via the Kneedle algorithm) to suppress AI coordinate jitter from posture changes.

### Joint Multi-Species Hierarchical EDD Model
The detection function scale parameter is decomposed as:

```
log(σ_ik) = α_i + δ_k
```

where `α_i` is a deployment-level random effect capturing site conditions and `δ_k` is a species-specific offset. A sum-to-zero constraint on species offsets ensures identifiability. This allows data-rich species to inform detection conditions at cameras where rare species have few observations.

### Staying Time via Bounding Box Interpolation
We construct a 1-second timeline for each detection sequence and interpolate gaps between trigger events using edge detection logic, distinguishing between animals remaining in the detection zone versus exiting and re-entering.

## R Package Dependencies

```r
install.packages(c(
  "tidyverse",
  "activity",
  "nimble",
  "MCMCvis",
  "coda",
  "patchwork",
  "scales"
))
```

## Citation

If you use this pipeline, please cite:

```
McMurry, S., Goldstein, B., Alyetam, M., & Kays, R. (in review). Automated Parameter 
Estimation for Camera Trap Density Models Using Computer Vision-Enhanced Distance 
Sampling. Ecological Monographs.
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

Sierra McMurry — Department of Forestry and Environmental Resources, North Carolina State University

For questions about the AI depth estimation pipeline, contact Mohammed Alyetam ([GitHub](https://github.com/Alyetama)).
