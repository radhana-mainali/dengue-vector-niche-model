library(terra)
library(dplyr)
library(ggplot2)

rm(list = ls())
setwd("/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling")

# Load the saved results
cat("Loading ENMeval results...\n")
all_runs <- readRDS("Analysis/MaxEnt_Output/all_runs_maxnet_ENMeval.rds")

cat(sprintf("Loaded %d model runs\n\n", length(all_runs)))

# =============================================================================
# 1. BASIC SUMMARY STATISTICS
# =============================================================================
cat("=== Summary Statistics ===\n")

# Best AICc for each run
best_aicc <- sapply(all_runs, function(x) min(x@results$AICc, na.rm = TRUE))

cat(sprintf("Mean best AICc: %.2f (SD=%.2f)\n", 
            mean(best_aicc, na.rm = TRUE), 
            sd(best_aicc, na.rm = TRUE)))
cat(sprintf("Range: %.2f - %.2f\n", 
            min(best_aicc, na.rm = TRUE), 
            max(best_aicc, na.rm = TRUE)))

# Best feature classes
best_fc <- sapply(all_runs, function(x) {
  best_idx <- which.min(x@results$AICc)
  as.character(x@results$fc[best_idx])
})

cat("\nBest model feature classes:\n")
print(table(best_fc))

# Best regularization multipliers
best_rm <- sapply(all_runs, function(x) {
  best_idx <- which.min(x@results$AICc)
  x@results$rm[best_idx]
})

cat("\nBest regularization multipliers:\n")
print(summary(best_rm))

# =============================================================================
# 2. DETAILED MODEL PERFORMANCE
# =============================================================================
cat("\n=== Model Performance Metrics ===\n")

# AUC values
best_auc <- sapply(all_runs, function(x) {
  best_idx <- which.min(x@results$AICc)
  x@results$auc.val.avg[best_idx]
})

cat(sprintf("Mean validation AUC: %.3f (SD=%.3f)\n", 
            mean(best_auc, na.rm = TRUE), 
            sd(best_auc, na.rm = TRUE)))
cat(sprintf("Range: %.3f - %.3f\n", 
            min(best_auc, na.rm = TRUE), 
            max(best_auc, na.rm = TRUE)))

# =============================================================================
# 3. CREATE SUMMARY DATA FRAME
# =============================================================================
cat("\n=== Creating Summary Table ===\n")

summary_df <- data.frame(
  run_name = names(all_runs),
  best_fc = best_fc,
  best_rm = best_rm,
  aicc = best_aicc,
  auc = best_auc,
  stringsAsFactors = FALSE
)

# Add bg_draw and split info
summary_df <- summary_df %>%
  mutate(
    bg_draw = as.integer(gsub("bg(\\d+)_split\\d+", "\\1", run_name)),
    split = as.integer(gsub("bg\\d+_split(\\d+)", "\\1", run_name))
  )

# Print first few rows
cat("\nFirst 10 rows of summary:\n")
print(head(summary_df, 10))

# Save summary table
write.csv(summary_df, 
          "Analysis/MaxEnt_Output/model_summary.csv", 
          row.names = FALSE)
cat("\nSummary table saved to: Analysis/MaxEnt_Output/model_summary.csv\n")

# =============================================================================
# 4. VISUALIZATIONS
# =============================================================================
cat("\n=== Creating Visualizations ===\n")

# Plot 1: AICc distribution
p1 <- ggplot(summary_df, aes(x = aicc)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7) +
  labs(title = "Distribution of Best Model AICc",
       x = "AICc", y = "Count") +
  theme_minimal()

ggsave("Analysis/MaxEnt_Output/aicc_distribution.png", p1, 
       width = 8, height = 6, dpi = 300)

# Plot 2: AUC distribution
p2 <- ggplot(summary_df, aes(x = auc)) +
  geom_histogram(bins = 20, color = "white", fill = "darkgreen", alpha = 0.7) +
  labs(title = "Distribution of Best Model AUC",
       x = "Validation AUC", y = "Count") +
  theme_minimal()

ggsave("Analysis/MaxEnt_Output/auc_distribution.png", p2, 
       width = 8, height = 6, dpi = 300)

# Plot 3: Feature class frequencies
fc_counts <- as.data.frame(table(best_fc))
colnames(fc_counts) <- c("Feature_Class", "Count")

p3 <- ggplot(fc_counts, aes(x = Feature_Class, y = Count)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  labs(title = "Best Model Feature Classes",
       x = "Feature Class", y = "Count") +
  theme_minimal()

ggsave("Analysis/MaxEnt_Output/feature_class_frequency.png", p3, 
       width = 8, height = 6, dpi = 300)

# Plot 4: Regularization multiplier distribution
p4 <- ggplot(summary_df, aes(x = factor(best_rm))) +
  geom_bar(fill = "purple", alpha = 0.7) +
  labs(title = "Best Model Regularization Multipliers",
       x = "Regularization Multiplier", y = "Count") +
  theme_minimal()

ggsave("Analysis/MaxEnt_Output/rm_frequency.png", p4, 
       width = 8, height = 6, dpi = 300)

cat("\nPlots saved to: Analysis/MaxEnt_Output/\n")

# =============================================================================
# 5. IDENTIFY BEST OVERALL MODEL
# =============================================================================
cat("\n=== Best Overall Model ===\n")

best_overall_idx <- which.min(summary_df$aicc)
best_model <- summary_df[best_overall_idx, ]

cat(sprintf("Best model: %s\n", best_model$run_name))
cat(sprintf("  Feature class: %s\n", best_model$best_fc))
cat(sprintf("  Regularization: %.1f\n", best_model$best_rm))
cat(sprintf("  AICc: %.2f\n", best_model$aicc))
cat(sprintf("  AUC: %.3f\n", best_model$auc))

# =============================================================================
# 6. LOAD AND EXAMINE BEST PREDICTION RASTER
# =============================================================================
cat("\n=== Best Prediction Raster ===\n")

best_raster_path <- file.path("Analysis/MaxEnt_Output/prediction_rasters",
                              paste0(best_model$run_name, "_best_model.tif"))

if (file.exists(best_raster_path)) {
  best_raster <- rast(best_raster_path)
  
  cat(sprintf("Loaded: %s\n", basename(best_raster_path)))
  cat(sprintf("Prediction range: %.4f - %.4f\n", 
              global(best_raster, "min", na.rm = TRUE)[1,1],
              global(best_raster, "max", na.rm = TRUE)[1,1]))
  
  # Save a simple plot
  png("Analysis/MaxEnt_Output/best_model_prediction_map.png", 
      width = 10, height = 8, units = "in", res = 600)
    plot(best_raster, col = rev(terrain.colors(100)), main = paste("Best Model Prediction"))
  dev.off()
  
  
  # save binary plot 
  m <- matrix(c(
    -Inf, 0.5, 0,
    0.5, Inf, 1
  ), ncol = 3, byrow = TRUE)
  bin_r <- classify(best_raster, m)

  bin_colors <- c("lightblue", "red")
  png("Analysis/MaxEnt_Output/best_model_binary_map.png",
      width = 10, height = 8, units = "in", res = 600)
    plot(bin_r, col = bin_colors, legend = FALSE, main = "Binary Prediction")
  dev.off()
  
  cat("\nPrediction map saved to: Analysis/MaxEnt_Output/best_model_prediction_map.png\n")
} 





# =============================================================================
# 7. FILTER HIGH-PERFORMING MODELS (AUC >= 0.85)
# =============================================================================
cat("\n=== High-Performing Models (AUC >= 0.85) ===\n")

# Filter models with AUC >= 0.85
high_auc_threshold <- 0.85
high_performers <- summary_df %>%
  filter(auc >= high_auc_threshold)

n_high <- nrow(high_performers)
cat(sprintf("Number of models with AUC >= %.2f: %d out of %d (%.1f%%)\n", 
            high_auc_threshold, 
            n_high, 
            nrow(summary_df),
            100 * n_high / nrow(summary_df)))



cat("\nSummary of high-performing models:\n")
cat(sprintf("  AUC range: %.3f - %.3f\n", 
            min(high_performers$auc), 
            max(high_performers$auc)))
cat(sprintf("  AICc range: %.2f - %.2f\n", 
            min(high_performers$aicc), 
            max(high_performers$aicc)))

# Save high-performing models list
write.csv(high_performers, 
          "Analysis/MaxEnt_Output/high_performing_models.csv", 
          row.names = FALSE)
cat("\nHigh-performing models saved to: Analysis/MaxEnt_Output/high_performing_models.csv\n")


# =============================================================================
# 8. COMPUTE ENSEMBLE MEAN AND CV FROM HIGH-PERFORMING MODELS
# =============================================================================
cat("\n=== Ensemble Predictions from High-Performing Models ===\n")

# Load all high-performing model rasters
raster_list <- list()

for (i in 1:nrow(high_performers)) {
  run_name <- high_performers$run_name[i]
  raster_path <- file.path("Analysis/MaxEnt_Output/prediction_rasters",
                           paste0(run_name, "_best_model.tif"))
  
  if (file.exists(raster_path)) {
    raster_list[[run_name]] <- rast(raster_path)
    cat(sprintf("  Loaded: %s (AUC=%.3f)\n", run_name, high_performers$auc[i]))
  } else {
    cat(sprintf("  Warning: Raster not found for %s\n", run_name))
  }
}

cat(sprintf("\nLoaded %d prediction rasters\n", length(raster_list)))

# Stack all rasters
pred_stack <- rast(raster_list)

# Compute mean prediction
cat("\nComputing ensemble mean...\n")
ensemble_mean <- app(pred_stack, fun = mean, na.rm = TRUE)
names(ensemble_mean) <- "ensemble_mean"

# Compute coefficient of variation (CV = SD / Mean)
cat("Computing coefficient of variation...\n")
ensemble_sd <- app(pred_stack, fun = sd, na.rm = TRUE)
ensemble_cv <- ensemble_sd / ensemble_mean
names(ensemble_cv) <- "ensemble_cv"

# Save ensemble rasters
writeRaster(ensemble_mean, 
            "Analysis/MaxEnt_Output/ensemble_mean_highAUC.tif", 
            overwrite = TRUE)
writeRaster(ensemble_cv, 
            "Analysis/MaxEnt_Output/ensemble_cv_highAUC.tif", 
            overwrite = TRUE)

cat("\nEnsemble rasters saved:\n")
cat("  Mean: Analysis/MaxEnt_Output/ensemble_mean_highAUC.tif\n")
cat("  CV: Analysis/MaxEnt_Output/ensemble_cv_highAUC.tif\n")

# Summary statistics
cat("\n--- Ensemble Mean Statistics ---\n")
cat(sprintf("  Range: %.4f - %.4f\n", 
            global(ensemble_mean, "min", na.rm = TRUE)[1,1],
            global(ensemble_mean, "max", na.rm = TRUE)[1,1]))
cat(sprintf("  Mean: %.4f\n", 
            global(ensemble_mean, "mean", na.rm = TRUE)[1,1]))

cat("\n--- Ensemble CV Statistics ---\n")
cat(sprintf("  Range: %.4f - %.4f\n", 
            global(ensemble_cv, "min", na.rm = TRUE)[1,1],
            global(ensemble_cv, "max", na.rm = TRUE)[1,1]))
cat(sprintf("  Mean CV: %.4f\n", 
            global(ensemble_cv, "mean", na.rm = TRUE)[1,1]))

# Create visualization plots
cat("\nCreating ensemble visualization plots...\n")

# Plot ensemble mean
png("Analysis/MaxEnt_Output/ensemble_mean_map.png", 
    width = 10, height = 8, units = "in", res = 600)
  plot(ensemble_mean, 
       main = sprintf("Ensemble Mean Prediction\n(n=%d models, AUC >= %.2f)", 
                      length(raster_list), high_auc_threshold),
       col = rev(terrain.colors(100)))
dev.off()


# save binary plot 
m
bin_ens <- classify(ensemble_mean, m)

png("Analysis/MaxEnt_Output/ensemble_binary_map.png",
    width = 10, height = 8, units = "in", res = 600)
  plot(bin_ens, col = bin_colors, legend = FALSE, main = "Ensemble Binary Prediction")
dev.off()    
    

# Plot ensemble CV
png("Analysis/MaxEnt_Output/ensemble_cv_map.png", 
    width = 10, height = 8, units = "in", res = 300)
  plot(ensemble_cv, 
       main = sprintf("Ensemble Coefficient of Variation\n(n=%d models, AUC >= %.2f)", 
                      length(raster_list), high_auc_threshold),
       col = heat.colors(100))
dev.off()

cat("\nEnsemble maps saved:\n")
cat("  Analysis/MaxEnt_Output/ensemble_mean_map.png\n")
cat("  Analysis/MaxEnt_Output/ensemble_cv_map.png\n")


# =============================================================================
# 9. UNCERTAINTY ANALYSIS
# =============================================================================
cat("\n=== Uncertainty Analysis ===\n")

# Classify CV into uncertainty categories
cv_vals <- values(ensemble_cv, na.rm = TRUE)
cv_quantiles <- quantile(cv_vals, probs = c(0.33, 0.67), na.rm = TRUE)

cat("CV Uncertainty Categories:\n")
cat(sprintf("  Low uncertainty (CV < %.3f): %.1f%% of area\n", 
            cv_quantiles[1],
            100 * sum(cv_vals < cv_quantiles[1], na.rm = TRUE) / length(cv_vals)))
cat(sprintf("  Medium uncertainty (%.3f <= CV < %.3f): %.1f%% of area\n", 
            cv_quantiles[1], cv_quantiles[2],
            100 * sum(cv_vals >= cv_quantiles[1] & cv_vals < cv_quantiles[2], na.rm = TRUE) / length(cv_vals)))
cat(sprintf("  High uncertainty (CV >= %.3f): %.1f%% of area\n", 
            cv_quantiles[2],
            100 * sum(cv_vals >= cv_quantiles[2], na.rm = TRUE) / length(cv_vals)))
    
  
    
    

    
# =============================================================================
# 10. VARIABLE IMPORTANCE FROM HIGH-PERFORMING MODELS
# =============================================================================
cat("\n=== Variable Importance Analysis ===\n")

# Extract variable importance from each high-performing model
var_importance_list <- list()

for (i in 1:nrow(high_performers)) {
  run_name <- high_performers$run_name[i]
  eval_obj <- all_runs[[run_name]]
  
  # Get the best model index for this run
  best_idx <- which.min(eval_obj@results$AICc)
  
  # Extract the model object
  best_model <- eval_obj@models[[best_idx]]
  
  # Get variable importance (permutation importance from maxnet)
  if (!is.null(best_model)) {
    # For maxnet models, we can extract beta coefficients as importance
    var_imp <- best_model$betas
    
    # Store in list
    var_importance_list[[run_name]] <- var_imp
    
    cat(sprintf("  Extracted variable importance for %s\n", run_name))
  }
}

    
# Combine all variable importance measures

# Get all unique variable names across models
all_vars <- unique(unlist(lapply(var_importance_list, names)))

# Create a matrix to store importance values
var_imp_matrix <- matrix(NA, nrow = length(all_vars), ncol = length(var_importance_list))
rownames(var_imp_matrix) <- all_vars
colnames(var_imp_matrix) <- names(var_importance_list)

# Fill the matrix
for (i in seq_along(var_importance_list)) {
  model_vars <- var_importance_list[[i]]
  for (var_name in names(model_vars)) {
    if (var_name %in% all_vars) {
      var_imp_matrix[var_name, i] <- abs(model_vars[var_name])
    }
  }
}

# Compute mean importance and standard deviation across models
mean_importance <- rowMeans(var_imp_matrix, na.rm = TRUE)
sd_importance <- apply(var_imp_matrix, 1, sd, na.rm = TRUE)

# Create summary data frame
var_imp_summary <- data.frame(
  variable = names(mean_importance),
  mean_importance = mean_importance,
  sd_importance = sd_importance,
  cv_importance = sd_importance / mean_importance
) %>%
  arrange(desc(mean_importance))

# Extract base variable names (remove feature class suffixes)
var_imp_summary$base_variable <- gsub("I\\(|\\)|\\^.*|\\:.*", "", var_imp_summary$variable)

# Aggregate by base variable
var_imp_aggregated <- var_imp_summary %>%
  group_by(base_variable) %>%
  summarise(
    total_importance = sum(mean_importance, na.rm = TRUE),
    mean_importance = mean(mean_importance, na.rm = TRUE),
    n_features = n()
  ) %>%
  arrange(desc(total_importance))


# Save variable importance
write.csv(var_imp_summary, 
          "Analysis/MaxEnt_Output/variable_importance_detailed.csv", 
          row.names = FALSE)
write.csv(var_imp_aggregated, 
          "Analysis/MaxEnt_Output/variable_importance_aggregated.csv", 
          row.names = FALSE)

cat("\nVariable importance tables saved:\n")
cat("  Analysis/MaxEnt_Output/variable_importance_detailed.csv\n")
cat("  Analysis/MaxEnt_Output/variable_importance_aggregated.csv\n")

# Create variable importance plot (aggregated)
top_n <- min(15, nrow(var_imp_aggregated))
top_vars <- var_imp_aggregated %>% head(top_n)

p_var <- ggplot(top_vars, aes(x = reorder(base_variable, total_importance), 
                              y = total_importance)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(title = sprintf("Variable Importance\n(All %d variables from %d high-performing models)", 
                       top_n, length(var_importance_list)),
       x = "Variable",
       y = "Total Importance (sum of absolute coefficients)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


ggsave("Analysis/MaxEnt_Output/variable_importance_plot.png", p_var, 
       width = 7, height = 8, dpi = 300)

cat("\nVariable importance plot saved:\n")
cat("  Analysis/MaxEnt_Output/variable_importance_plot.png\n")

# Create detailed plot showing variability across models
p_var_detailed <- ggplot(var_imp_summary, 
                         aes(x = reorder(variable, mean_importance), 
                             y = mean_importance)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  geom_errorbar(aes(ymin = pmax(0, mean_importance - sd_importance), 
                    ymax = mean_importance + sd_importance), 
                width = 0.3, alpha = 0.5) +
  coord_flip() +
  labs(title = sprintf("Variable Importance with Variability\n(All 22 features from %d high-performing models)", 
                       length(var_importance_list)),
       x = "Variable Feature",
       y = "Mean Importance ± SD") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))

ggsave("Analysis/MaxEnt_Output/variable_importance_detailed_plot.png", 
       p_var_detailed, 
       width = 7, height = 10, dpi = 300)

cat("  Analysis/MaxEnt_Output/variable_importance_detailed_plot.png\n")
      
    

# ==========================================================================
# 11. RESPONSE CURVES FOR TOP 6 VARIABLES
# ==========================================================================
cat("\n=== Response Curves for Top Variables ===\n")

# Get top 6 variables
top_6_vars <- var_imp_aggregated$base_variable[1:min(6, nrow(var_imp_aggregated))]
cat(sprintf("Creating response curves for top %d variables:\n", length(top_6_vars)))
cat(paste("  -", top_6_vars, collapse = "\n"))
cat("\n")

# Load environmental data
env <- rast("Data_processed/Covariates/covariates_uncorrelated.tif")

# Create response curves for each top variable
response_data_list <- list()

for (var_name in top_6_vars) {
  cat(sprintf("Processing %s...\n", var_name))
  
  # Check if variable exists in raster
  if (var_name %in% names(env)) {
    
    # Get the range of the variable
    var_raster <- env[[var_name]]
    var_values <- values(var_raster, na.rm = TRUE)
    var_range <- range(var_values, na.rm = TRUE)
    
    # Create sequence of values across the range
    var_seq <- seq(var_range[1], var_range[2], length.out = 100)
    
    # Store predictions from each high-performing model
    pred_matrix <- matrix(NA, nrow = length(var_seq), ncol = length(high_performers$run_name))
    
    for (i in 1:nrow(high_performers)) {
      run_name <- high_performers$run_name[i]
      eval_obj <- all_runs[[run_name]]
      best_idx <- which.min(eval_obj@results$AICc)
      best_model <- eval_obj@models[[best_idx]]
      
      if (!is.null(best_model)) {
        # Create prediction data frame with all variables at their mean
        # except the variable of interest
        pred_data <- as.data.frame(matrix(nrow = length(var_seq), 
                                          ncol = nlyr(env)))
        colnames(pred_data) <- names(env)
        
        # Set all variables to their mean
        for (col in names(env)) {
          pred_data[[col]] <- global(env[[col]], "mean", na.rm = TRUE)[1,1]
        }
        
        # Vary the target variable
        pred_data[[var_name]] <- var_seq
        
        # Predict using the model
        preds <- predict(best_model, pred_data, type = "cloglog")
        pred_matrix[, i] <- preds
      }
    }
    
    # Compute mean and SD across models
    mean_pred <- rowMeans(pred_matrix, na.rm = TRUE)
    sd_pred <- apply(pred_matrix, 1, sd, na.rm = TRUE)
    
    # Store in list
    response_data_list[[var_name]] <- data.frame(
      variable = var_name,
      value = var_seq,
      mean_response = mean_pred,
      sd_response = sd_pred,
      lower_ci = mean_pred - 1.96 * sd_pred / sqrt(ncol(pred_matrix)),
      upper_ci = mean_pred + 1.96 * sd_pred / sqrt(ncol(pred_matrix))
    )
    
  } else {
    cat(sprintf("  Warning: Variable %s not found in environmental data\n", var_name))
  }
}

# Combine all response data
response_df <- bind_rows(response_data_list)

# Save response curve data
write.csv(response_df, 
          "Analysis/MaxEnt_Output/response_curves_data.csv", 
          row.names = FALSE)
cat("\nResponse curve data saved to: Analysis/MaxEnt_Output/response_curves_data.csv\n")

# Create individual plots for each variable
cat("\nCreating response curve plots...\n")

for (var_name in top_6_vars) {
  if (var_name %in% response_df$variable) {
    var_data <- response_df %>% filter(variable == var_name)

    p_resp <- ggplot(var_data, aes(x = value, y = mean_response)) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
                  fill = "steelblue", alpha = 0.3) +
      geom_line(color = "darkblue", size = 1.2) +
      labs(title = paste("Response Curve:", var_name),
           subtitle = sprintf("Mean ± 95%% CI from %d models",
                              nrow(high_performers)),
           x = var_name,
           y = "Predicted Suitability") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 14))

    ggsave(paste0("Analysis/MaxEnt_Output/response_curve_", var_name, ".png"),
           p_resp, width = 8, height = 6, dpi = 300)

    cat(sprintf("  Saved: response_curve_%s.png\n", var_name))
  }
}

# Create combined plot with all response curves
p_combined <- ggplot(response_df, aes(x = value, y = mean_response)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = variable), 
              alpha = 0.2) +
  geom_line(aes(color = variable), size = 1) +
  facet_wrap(~variable, scales = "free_x", ncol = 3) +
  labs(title = sprintf("Response Curves for Top %d Variables", length(top_6_vars)),
       subtitle = sprintf("Mean ± 95%% CI from %d high-performing models", 
                          nrow(high_performers)),
       x = "Variable Value",
       y = "Predicted Suitability") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10))

ggsave("Analysis/MaxEnt_Output/response_curves_combined.png", 
       p_combined, width = 8, height = 8, dpi = 300)

cat("  Saved: response_curves_combined.png\n")


cat("\n=== Analysis Complete ===\n")    
