library(terra)
library(ggplot2)
library(reshape2)
library(terra)
library(sf)
library(caret)  # for findCorrelation
library(car)    # for VIF


rm(list = ls())

# ------------------------------------------------------
# Load WorldClim & ENVIREM rasters and combine
# ------------------------------------------------------

covdir <- "/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling/Data_processed/Covariates"

# Load both multiband rasters
wc <- rast(file.path(covdir, "worldclim30s_asia.tif"))
env <- rast(file.path(covdir, "envirem30s_asia.tif"))
pop <- rast(file.path(covdir, "humanpop30s_asia.tif"))
lc <- rast(file.path(covdir, "landcover30s_asia.tif"))



# custon function to check extent, resolution, origin, grid alignment, crs
check_match <- function(...){
  rasters <- list(...)
  if(length(rasters) < 2) stop("Need at least two rasters")
  ref <- rasters[[1]]
  # geometry attributes for reference
  ref_ext  <- ext(ref)
  ref_res  <- res(ref)
  ref_dim  <- dim(ref)[1:2]
  ref_orig <- origin(ref)
  ref_crs  <- crs(ref)
  
  for (i in 2:length(rasters)) {
    r <- rasters[[i]]
    cat("Comparing raster", i, "to reference\n")
    ok_ext  <- isTRUE(all.equal(ext(r),  ref_ext))
    ok_res  <- isTRUE(all.equal(res(r),  ref_res))
    ok_dim  <- all(dim(r)[1:2] == ref_dim)
    ok_orig <- isTRUE(all.equal(origin(r), ref_orig))
    ok_crs  <- identical(crs(r), ref_crs)
    
    cat("Summary of equality checks:\n")
    cat("  Extent identical?     ", ok_ext,  "\n")
    cat("  Resolution identical? ", ok_res,  "\n")
    cat("  Dimensions identical? ", ok_dim,  "\n")
    cat("  Origin identical?     ", ok_orig, "\n")
    cat("  CRS identical?        ", ok_crs,  "\n\n")
  }
}

check_match(wc[[1]], env[[1]])
check_match(wc[[1]], pop[[1]])
check_match(wc[[1]], lc[[1]])



# Combine into one stack
covars <- c(wc, env, pop, lc)

all.equal(ext(wc), ext(env), ext(pop), ext(lc))
# all.equal(res(wc), res(env), res(pop), res(lc))

# Resolution in arc-seconds
res_arcsec <- res(wc) * 3600
res_arcsec

# Resolution in meters (approximate at equator)
res_meters <- res(wc)[1] * 111320   # 1 degree ≈ 111.32 km
res_meters

covars
names(covars)
plot(covars)



# ------------------------------------------------------------------------
# Load study area polygon and presence points, and extract covariates
# ------------------------------------------------------------------------

# SpatVector of study/background polygon
study_area_sf <- st_read("/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling/Data_processed/Shapefiles/asia_kalaazar_study_area.shp")

# to use it with terra rasters, convert to SpatVector
study_v <- vect(study_area_sf)

plot(covars[[1]])
plot(study_v, add=TRUE)


# sf of presence points (lon/lat WGS84)
pres_pts1 <- st_read("/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling/Data_processed/Occurrences/Kraemer/aegypti_clipped.shp")
pres_pts2 <- st_read("/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling/Data_processed/Occurrences/Kraemer/albopictus_clipped.shp")
points(pres_pts1, col="red", pch=20, cex=0.5)
points(pres_pts2, col="red", pch=20, cex=0.5)

pres_pts1 <- pres_pts1 %>% mutate(STATUS = as.character(STATUS))
pres_pts2 <- pres_pts2 %>% mutate(STATUS = as.character(STATUS))
pres_pts <- bind_rows(pres_pts1, pres_pts2)
head(pres_pts)
dim(pres_pts)
# st_write(pres_pts, "/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling/Data_processed/Occurrences/Kraemer/both_clipped.shp", delete_dsn = TRUE)


# Extract covariate values at those pseudo-presence points
pres_df <- as.data.frame(extract(covars, pres_pts))
head(pres_df)
dim(pres_df)


# ------------------------------------------------------
# Sample background points and extract covariates
# ------------------------------------------------------

bg_pts <- spatSample(covars, size = 10000, method = "random", na.rm = TRUE, as.points = TRUE)  # or "random"
points(bg_pts, col="gray", pch=1, cex=0.5)

## Extract values
bg_df   <- as.data.frame(extract(covars, bg_pts))
head(bg_df)
dim(bg_df)


# ------------------------------------------------------
# Prepare predictor matrix (drop NAs, near-zero variance)
# ------------------------------------------------------

# merge presences and background points for colliinearity analyssi
X <- rbind(pres_df, bg_df); head(X)
X <- X[, -1, drop = FALSE]     # drop ID column
head(X)

X <- X[, sapply(X, function(z) sum(!is.na(z)) > 0)]        # drops a columns if it is entirely NA

# find predictors that:
# have very few unique values relative to the number of samples (e.g., almost constant everywhere).
# Or are extremely unbalanced (e.g., 99% of the data are one value, 1% another).
# Such variables are essentially uninformative for modeling and can cause problems in correlation or VIF calculations.
nzv <- nearZeroVar(X) # Identify columns with near-zero variance
if (length(nzv)) X <- X[, -nzv] # If any were found, drop them
head(X); dim(X)
X <- X[complete.cases(X), ]
head(X); dim(X)


# ------------------------------------------------------
# Correlation filter (Pearson) and visualization
# ------------------------------------------------------

cm <- cor(X, use = "pairwise.complete.obs", method = "pearson")
cm

# ----------------------
threshold <- 0.75
# ----------------------

cm_long <- melt(cm)
ggplot(cm_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(
    aes(label = round(value, 2),
        fontface = ifelse(abs(value) > threshold, "bold", "plain")),
    size = 2.5
  ) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text = element_text(size = 9)
  ) +
  labs(x = "", y = "", fill = "r")

# find predictors to drop that have high correlations
to_drop <- findCorrelation(cm, cutoff = threshold, names = TRUE, exact = TRUE); to_drop
length(to_drop)

# reduced predictors set
dim(X)
X_red <- X[, setdiff(colnames(X), to_drop)]
head(X_red)
dim(X_red)

# check one more time if anythis is above the threshold
cm_red <- cor(X_red, use = "pairwise.complete.obs", method = "pearson")
cm_red
cm_red_long <- melt(cm_red)
ggplot(cm_red_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(
    aes(label = round(value, 2),
        fontface = ifelse(abs(value) > threshold, "bold", "plain")),
    size = 2.5
  ) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text = element_text(size = 9)
  ) +
  labs(x = "", y = "", fill = "r")


to_drop_red <- findCorrelation(cm_red, cutoff = threshold, names = TRUE, exact = TRUE); to_drop_red
length(to_drop_red)

# so all looks good now


# ------------------------------------------------------
# VIF screen (iterative) on correlation-reduced set
# ------------------------------------------------------
vif_drop <- function(df, thr = 5){
  repeat {
    # fit dummy linear model; VIF depends only on predictor matrix
    fit <- lm(as.matrix(df[,1]) ~ ., data = as.data.frame(df))
    v  <- vif(fit)
    vmax <- max(v, na.rm = TRUE)
    if (vmax < thr) break
    worst <- names(which.max(v))
    message("Dropping variable with highest VIF: ", worst, " (VIF=", round(vmax, 2), ")")
    df <- df[, setdiff(colnames(df), worst)]
  }
  colnames(df)
}

# run VIF screening on scaled X_red
keep_names <- vif_drop(scale(X_red), thr = 10)
keep_names
length(keep_names)


# ------------------------------------------------------
# Subset raster stack with selected variables
# ------------------------------------------------------

covars
nlyr(covars)
names(covars)

covars_sel <- covars[[keep_names]]
covars_sel
names(covars_sel)

# ------------------------------------------------------
# Save for MaxEnt
# ------------------------------------------------------
writeRaster(covars_sel, "/Users/kpmainali/Dropbox/Documents/RESEARCH/DengueModeling/Data_processed/Covariates/covariates_uncorrelated.tif", overwrite = TRUE)
