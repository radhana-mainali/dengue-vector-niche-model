library(terra)
library(sf)
library(dplyr)
library(ENMeval)

rm(list = ls())
setwd("/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling")

# 1. Load presence points (sf or data.frame) & environmental raster stack
pres <- st_read("Data_processed/Occurrences/Kraemer/both_clipped.shp")
env <- rast("Data_processed/Covariates/covariates_uncorrelated.tif")  

# 2. Define extent/region for background sampling (study area)
study_area <- st_read("Data_processed/Shapefiles/asia_kalaazar_study_area.shp")
study_mask   <- crop(mask(env[[1]], study_area), study_area)
study_raster <- study_mask  # single-layer for sampling

# 3. Create output directory for rasters
raster_dir <- "Analysis/MaxEnt_Output/prediction_rasters"
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

# 4. Parameters
n_bg_draws <- 5
n_splits <- 5
bg_size    <- 10000

# Determine number of cores to use
n_cores <- parallel::detectCores() - 1
cat(sprintf("ENMeval will use %d cores for parallel processing\n\n", n_cores))

all_runs <- list()

# 5. Loop through bg_draws and splits
cat("Starting model runs...\n")
start_time <- Sys.time()

for (i in seq_len(n_bg_draws)) {
  
  cat(sprintf("\n=== Background draw %d/%d ===\n", i, n_bg_draws))
  
  # Sample background points
  set.seed(i * 1000)  # Set seed for reproducible bg sampling
  bg_raw <- terra::spatSample(study_raster, size = bg_size, 
                              method = "random", na.rm = TRUE, xy = TRUE)
  bg_df <- as.data.frame(bg_raw)[,1:2]
  colnames(bg_df) <- c("longitude","latitude")
  
  # Presence coordinates
  pres_xy <- pres %>% st_coordinates() %>% as.data.frame() %>% 
    rename(longitude = X, latitude = Y)
  
  for (j in seq_len(n_splits)) {
    set.seed(i*1000 + j)
    
    cat(sprintf("  Running split %d/%d...", j, n_splits))
    
    # Run ENMevaluate with parallel = TRUE
    # ENMeval handles all parallel setup internally
    eval_res <- ENMevaluate(
      occs       = pres_xy,
      bg         = bg_df,
      envs       = env,
      algorithm  = "maxnet",
      partitions = "randomkfold",
      partition.settings = list(kfolds = 5),
      tune.args  = list(
        fc = c("L", "LQ"),
        rm = seq(0.5, 3, by = 0.5)
      ),
      raster.preds = TRUE,       # Generate prediction rasters
      parallel     = TRUE,       # ENMeval handles parallel backend internally
      numCores     = n_cores,    # Explicitly specify number of cores
      updateProgress = FALSE     # Disable progress bar for cleaner output
    )
    
    # Store the eval_res object
    run_name <- paste0("bg", i, "_split", j)
    all_runs[[run_name]] <- eval_res
    
    # Extract and save prediction raster for the best model
    best_idx <- which.min(eval_res@results$AICc)
    best_preds <- eval_res@predictions[[best_idx]]
    
    # Save the best model prediction raster
    writeRaster(best_preds, 
                filename = file.path(raster_dir, paste0(run_name, "_best_model.tif")),
                overwrite = TRUE)
    
    cat(sprintf(" Done (Best: fc=%s rm=%.1f, AICc=%.2f)\n", 
                eval_res@results$fc[best_idx], 
                eval_res@results$rm[best_idx],
                eval_res@results$AICc[best_idx]))
  }
}

end_time <- Sys.time()
cat(sprintf("\nAll runs completed in %.2f minutes\n", 
            as.numeric(difftime(end_time, start_time, units = "mins"))))

# Save all_runs object (contains models and metadata)
saveRDS(all_runs, 
        file = "Analysis/MaxEnt_Output/all_runs_maxnet_ENMeval.rds")

cat(sprintf("\nResults saved:\n"))
cat(sprintf("  Prediction rasters: %s\n", raster_dir))
cat(sprintf("  RDS object: Analysis/MaxEnt_Output/all_runs_maxnet_ENMeval.rds\n"))
