# ==============================================================================
# Supplement 2: Clinical Phenotype Calculator with Projection Mapping
# Developed for Red Blood Cell Fatty Acid Profiling
# 
# Instructions:
# 1. Ensure "Anonymous_Phenotype_Core.RData" is in the same working directory.
# 2. Input dataframe should contain 29 required fatty acids and an 'ID' column.
# ==============================================================================

library(ggplot2)

predict_phenotype <- function(new_raw_df, core_file = "Anonymous_Phenotype_Core.RData") {
  
  # 1. Environment and File Check
  if(!file.exists(core_file)) {
    stop("❌ Error: Core data file 'Anonymous_Phenotype_Core.RData' not found in the working directory.")
  }
  load(core_file)
  
  # 2. Variable Standardization
  raw_names <- colnames(new_raw_df)
  raw_names <- gsub(":", "_", raw_names)
  raw_names <- gsub("^C", "c", raw_names)
  colnames(new_raw_df) <- tolower(raw_names)
  
  missing_vars <- setdiff(base_vars, colnames(new_raw_df))
  if(length(missing_vars) > 0) {
    stop(paste("❌ Missing required fatty acid variables:", paste(missing_vars, collapse=", ")))
  }
  
  # 3. Numeric Conversion and Percentage Formatting
  new_mat <- matrix(0, nrow = nrow(new_raw_df), ncol = length(base_vars))
  colnames(new_mat) <- base_vars
  for(v in base_vars) {
    val_raw <- as.character(new_raw_df[[v]])
    has_percent <- grepl("%", val_raw)
    suppressWarnings({
      num_vec <- as.numeric(gsub("%", "", val_raw))
    })
    new_mat[, v] <- ifelse(has_percent, num_vec / 100, num_vec)
  }
  
  # 4. Data Integrity Interceptor (NA/Missing Value Check)
  if(any(is.na(new_mat))) {
    na_rows <- which(rowSums(is.na(new_mat)) > 0)
    bad_ids <- if("id" %in% colnames(new_raw_df)) new_raw_df$id[na_rows] else na_rows
    stop(paste0(
      "\n❌ Data Integrity Alert!\n",
      "Missing values (NA) or non-numeric characters (e.g., 'ND', '<0.01') detected for the following Patient IDs:\n",
      paste(bad_ids, collapse = ", "), "\n",
      "PCA projection requires complete feature matrices. Please correct the input records."
    ))
  }
  
  # 5. Dual-Engine Zero Imputation
  # (Dynamically adapts to local batch LOD; falls back to baseline LOD for N=1)
  clean_mat <- new_mat
  for(v in base_vars) {
    x <- clean_mat[, v]
    is_zero <- (!is.na(x) & x == 0)
    
    if(any(is_zero)) {
      if(any(x > 0, na.rm = TRUE)) {
        fill_val <- min(x[x > 0], na.rm = TRUE) / 2
      } else {
        fill_val <- base_min_halves[[v]] 
      }
      clean_mat[is_zero, v] <- fill_val
    }
  }
  
  # 6. Centered Log-Ratio (CLR) Transformation
  geom_mean <- function(x) { exp(mean(log(x[!is.na(x)]))) } 
  clr_mat <- t(apply(clean_mat, 1, function(row) {
    g <- geom_mean(row)
    log(row / g)
  }))
  
  # 7. Absolute PCA Projection and Phenotype Assignment (99% Threshold)
  new_pca_scores <- scale(clr_mat, center = base_pca_center, scale = FALSE) %*% base_pca_rotation
  
  assign_logic <- function(point) {
    dists <- apply(base_centroids, 1, function(c) sqrt(sum((point - c)^2)))
    best <- which.min(dists)
    if(dists[best] > base_thresholds[best]) return("Unclassified (Outlier)")
    return(paste("Phenotype", best))
  }
  results <- apply(new_pca_scores, 1, assign_logic)
  
  # 8. Result Assembly
  all_pheno_levels <- c("Phenotype 1", "Phenotype 2", "Phenotype 3", 
                        "Phenotype 4", "Phenotype 5", "Unclassified (Outlier)")
  
  final_df <- data.frame(
    ID = if("id" %in% colnames(new_raw_df)) new_raw_df$id else 1:nrow(new_raw_df),
    Phenotype = factor(results, levels = all_pheno_levels)
  )
  final_df <- cbind(final_df, new_pca_scores[, 1:2, drop = FALSE])
  colnames(final_df)[3:4] <- c("PC1", "PC2")
  
  # 9. Visualization Engine (With Ghost Data for Complete Legend)
  cluster_colors <- c("1"="#B7B2D0", "2"="#E68B81", "3"="#7DA6C6", "4"="#84C3B7", "5"="#EAAA60")
  pheno_colors <- c("Phenotype 1"="#B7B2D0", "Phenotype 2"="#E68B81", 
                    "Phenotype 3"="#7DA6C6", "Phenotype 4"="#84C3B7", 
                    "Phenotype 5"="#EAAA60", "Unclassified (Outlier)"="black")
  
  ghost_df <- data.frame(
    PC1 = Inf, PC2 = Inf, # Moved to Inf to absolutely prevent rendering overlap
    Phenotype = factor(all_pheno_levels, levels = all_pheno_levels)
  )
  
  p <- ggplot() +
    geom_polygon(data = base_ellipses, aes(x = PC1, y = PC2, group = Cluster, fill = factor(Cluster)), 
                 alpha = 0.15, color = "gray50", linetype = "dashed") +
    scale_fill_manual(values = cluster_colors, name = "Reference Territories", drop = FALSE) +
    
    geom_point(data = ghost_df, aes(x = PC1, y = PC2, color = Phenotype), alpha = 0, size = 0.1) +
    geom_point(data = final_df, aes(x = PC1, y = PC2, color = Phenotype), size = 3, alpha = 0.8) +
    
    scale_color_manual(values = pheno_colors, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3, shape = 16))) +
    
    labs(x = "PC1", y = "PC2", title = "Phenotype Projection Map",
         subtitle = "New cohort projection onto baseline metabolic territories") +
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),
      
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      
      axis.title = element_text(size = 12),  
      axis.text = element_text(size = 12),  
      
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"), 
      legend.text = element_text(size = 12),               

      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(list(Data = final_df, Plot = p))
}

# ==============================================================================
# Example Usage (Uncomment to test):
# library(readxl)
# df_external <- read_excel("example_validation_cohort.xlsx")
# result <- predict_phenotype(df_external)
# print(result$Plot)
# head(result$Data)
# ==============================================================================