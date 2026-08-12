# Whole-Blood-Fatty-Acid-Phenotype-Calculator

Supplement 2 for the manuscript: "Whole-blood fatty acid phenotypes and subsequent quality-of-life changes among cancer survivors"

## Overview

This repository contains the R script, core reference parameters, and example dataset required to run the Whole-Blood Fatty Acid Phenotype Calculator described in the manuscript. The tool supports reproducible phenotype assignment in external research cohorts by projecting new whole-blood fatty acid profiles into the fixed baseline phenotype space.

## Important input requirements

### Whole-blood profiles only

This calculator was developed using whole-blood fatty acid profiles measured from dried blood spot samples. It is intended for application to whole-blood fatty acid data generated using a comparable analytical workflow. It should not be applied directly to plasma, serum, or other analytical matrices, because matrix-specific differences in fatty acid composition may result in inappropriate phenotype assignments.

### Compositional data format

Input data should be provided as relative abundances on a fractional scale, with fatty acid values summing approximately to 1.0. For example, 1.5% should be entered as 0.015.

The calculator can automatically detect and rescale values containing percentage signs, such as "1.5%". Whole-number percentages without a percentage sign, such as 1.5 to represent 1.5%, should not be used. If such values are detected, users should rescale the data before running the calculator.

## Contents of this repository:
1. Phenotype_Calculator.R
   The executable R script containing the `predict_phenotype()` function.
2. Anonymous_Phenotype_Core.RData
   The serialized core parameter file containing the baseline PCA rotation matrix, cluster centroids, 99% distance rejection thresholds, and fallback zero-imputation limits.
3. example_validation_cohort.xlsx
   A de-identified example dataset (N=12) demonstrating phenotype assignment and outlier detection.

### System Requirements:
- R (version 4.3.3 or higher recommended)
- Required R packages: `ggplot2` (for visualization), `readxl` (for reading the example file)

### Instructions for Use:
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

### Input Data Format Requirements:
- An ID column is recommended. If no ID column is provided, row numbers will be used as profile identifiers.
- Percentage signs: The calculator automatically detects values containing percentage signs. If your input data contains the "%" symbol (e.g., "1.5%"), the tool will automatically strip the character and convert the value to the required fractional scale (0.015) for accurate PCA projection.
- Numeric Consistency: Whole-number percentages without a percentage sign, such as 1.5 to represent 1.5%, should not be used because they are interpreted as fractional-scale inputs.
- Zero values are automatically imputed using half the minimum nonzero value observed within the external dataset, with baseline half-minimum values used as fallback when no positive value is available for a required fatty acid.
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

### Methodological Details:
For a comprehensive explanation of the mathematical framework, please refer to "eMethods 3" in the Supplementary Appendix associated with the manuscript.
