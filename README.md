# Erythrocyte-Fatty-Acid-Phenotype-Calculator
Supplement 2 for manuscript: "Erythrocyte Fatty Acid Phenotypes and Quality of Life Trajectories in Cancer Survivors"

Overview:
This archive contains the necessary R scripts, core reference parameters, and a sample dataset to execute the Phenotype Calculator described in the manuscript. 
This tool allows researchers and clinicians to project new, external red blood cell (RBC) fatty acid profiles onto the established baseline metabolic territories.

⚠️ CRITICAL REQUIREMENT (RBC vs. Plasma):
1. RBC Profiles Only (No Plasma/Serum): > This calculator is strictly calibrated for Red Blood Cell (RBC) membrane fatty acid profiles, which reflect long-term (approx. 120 days) dietary intake and metabolic homeostasis. It MUST NOT be applied to plasma or serum fatty acid data (which reflect short-term/recent meals), as this will result in severe misclassification.
2. Compositional Data Format (Fractional Scale):
Input data MUST be provided as relative abundance (decimals), as the calculator is strictly calibrated for a fractional scale (0.0–1.0). While the standard input format uses decimals (e.g., 0.015 for 1.5%), the internal engine features a smart interceptor that automatically detects and rescales data containing percentage signs (e.g., "1.5%" is rescaled to 0.015). However, it is STRICTLY PROHIBITED to input whole-number percentages without a percentage sign (e.g., entering 1.5 to represent 1.5%); doing so will cause a 100-fold scale error.

Contents of this Archive:
1. Phenotype_Calculator.R
   The executable R script containing the `predict_phenotype()` function.
2. Anonymous_Phenotype_Core.RData
   The serialized core parameter file containing the baseline PCA rotation matrix, cluster centroids, 99% distance rejection thresholds, and fallback zero-imputation limits.
3. example_validation_cohort.xlsx
   A desensitized, representative sample dataset (N=12) demonstrating various 
   phenotype assignments and outlier detection.

System Requirements:
- R (version 4.3.3 or higher recommended)
- Required R packages: `ggplot2` (for visualization), `readxl` (for reading the example file)

Instructions for Use:
1. Download all files from this repository into a single local directory (which will serve as your Working Directory).
2. Open R or RStudio and set your Working Directory to that downloaded folder.
   (e.g., `setwd("path/to/extracted/folder")`)
3. Load the calculator function by running:
   `source("Phenotype_Calculator.R")`
4. Load the required packages and run the example data:
   `library(readxl)`
   `df_external <- read_excel("example_validation_cohort.xlsx")`
   `result <- predict_phenotype(df_external)`
5. View the results:
   `head(result$Data)`  # Displays the classification results and PC coordinates
   `print(result$Plot)` # Renders the 2D Phenotype Projection Map

Input Data Format Requirements:
- The input dataframe must contain an 'ID' column.
- Smart Scaling for "%": The calculator features a smart interceptor for percentage signs. If your clinical data contains the "%" symbol (e.g., "1.5%"), the tool will automatically strip the character and convert the value to the required fractional scale (0.015) for accurate PCA projection.
- Numeric Consistency: If no percentage sign is present, the values must already be in decimal format (summing to approximately 1.0). Pure numeric inputs greater than 1.0 (without a "%" sign) will be treated as raw values and lead to erroneous results.
- Zero values are automatically imputed using a batch-adaptive algorithm to account for multi-center limits of detection (LOD).
- REQUIRED VARIABLES: The input data must contain the following 29 specific fatty acid variables (column names are case-insensitive and can use either colons or underscores, e.g., C16:0 or c16_0):
  1. C14:0 (myristic acid)
  2. C15:0 (pentadecanoic acid)
  3. C15:1 (pentadecenoic acid)
  4. C16:0 (palmitic acid)
  5. C16:1 (palmitoleic acid)
  6. C17:0 (margaric acid)
  7. C17:1 (heptadecenoic acid)
  8. C18:0 (stearic acid)
  9. C18:1n9c (oleic acid)
  10. C18:1n9t (elaidic acid)
  11. C18:2n6c (linoleic acid)
  12. C18:2n6t (linolelaidic acid)
  13. C18:3n3 (α-linolenic acid)
  14. C18:3n6 (γ-linolenic acid)
  15. C20:0 (arachidic acid)
  16. C20:1n9 (gondoic acid)
  17. C20:2 (eicosadienoic acid)
  18. C20:3n6 (dihomo-γ-linolenic acid)
  19. C20:4n6 (arachidonic acid)
  20. C20:5n3 (eicosapentaenoic acid)
  21. C22:0 (behenic acid)
  22. C22:1 (erucic acid)
  23. C22:4n6 (adrenic acid)
  24. C22:5n3 (docosapentaenoic acid)
  25. C22:5n6 (osbond acid)
  26. C22:6n3 (docosahexaenoic acid)
  27. C23:0 (tricosanoic acid)
  28. C24:0 (lignoceric acid)
  29. C24:1n9 (nervonic acid)

Methodological Details:
For a comprehensive explanation of the mathematical framework, please refer to "eMethods 3" in the Supplementary Appendix associated with the manuscript.
================================================================================
