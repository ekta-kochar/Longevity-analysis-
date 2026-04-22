# README for Data and Code Associated with:
Title of Manuscript: Mitochondrial genotypes affect thermal plasticity of longevity in Drosophila melanogaster
Authors: Ekta Kochar, Venkatesh Nagarajan-Radha, Rebecca E. Koch and Damian K Dowling 
Journal: Journal of Evolutionary Biology
Corresponding Author: Ekta Kochar (ekta.kochar.1502@gmail.com)  
Date: 2025-08-25  

---

## 1. Study Summary
This repository contains the datasets and analysis scripts associated with the manuscript "Mitochondrial genotypes affect thermal plasticity of longevity in Drosophila melanogaster" submitted to the Journal of Evolutionary Biology. In this study, we test how mitochondrial variation influences the longevity of flies under divergent thermal conditions. We use laboratory strains carrying different mitochondrial haplotypes expressed in an isogenic nuclear background and measure lifespan at two temperatures (18 °C and 28 °C). Our results show that mitochondrial variation affects longevity, though not always in the predicted direction, and we explicitly test both the mitochondrial climatic adaptation hypothesis and the mother’s curse hypothesis.

Data were collected by Ekta Kochar and Venkatesh Nagarajan-Radha, and the analysis code was written by Ekta Kochar.  


---

## 2. Repository Contents
├── data/
│ ├── Longevity_Data.csv # Original experimental data
│ ├── Longevity_MetaData.csv # # Sample information (genotype, sex, temperature, etc.)
├── scripts/
│ ├── Longevity_Rscript.R # Complete workflow: packages, data formatting, mixed-effects models, plots
  └── README.md # This file


---

## 3. Data Description
- Longevity_Data.csv: Dataset containing all variables ane measured trait, with one row per observation.  
- Longevity_MetaData.csv: Describes experimental design variables such as genotype, sex, and temperature.  


## 4. Code Description
- **Longevity_Rscript.R**:  
  This script contains the complete workflow used in the manuscript. It is organized into the following steps:  

  1. **Load packages**  
     Loads all required R packages:  
     - `lme4` / `lmerTest` (linear mixed-effects models)  
     - `carData`, `car` (ANOVA with type II/III sums of squares)  
     - `emmeans`, `multcomp` (post-hoc tests and contrasts)  
     - `effects`, `Hmisc` (effect plots, descriptive statistics)  
     - `ggplot2`, `ggbeeswarm` (visualizations)  

  2. **Set contrasts**  
     Specifies sum-to-zero contrasts (`contr.sum`) to ensure Type III ANOVA results are interpretable.  

  3. **Import and format data**  
     - Reads in the dataset (`Longevity_Data.csv`).  
     - Converts categorical predictors (e.g., `hapgroup`, `hapdup`, `location`, `block`, `bio`, `vial`, `temp`, `sex`) into factors.  

  4. **Fit linear mixed-effects model**  
     - Uses `lme4::lmer()` to model lifespan (`age`) as a function of:  
       - Fixed effects: haplogroup/haplotype, sex, temperature, block, and their interactions.  
       - Nested structure: haplotypes nested within haplogroups.  
       - Random effects:  
         - `vial` (random intercept, cohort-level variation).  
         - `bio` (strain ID, random intercepts and random slopes for sex, temp, and their interaction).  
     - Model fitted with REML estimation.  

  5. **Model outputs**  
     - Summarizes model results with `summary()`.  
     - Performs Type III ANOVA with `car::Anova()`.  
  
  6. **Visualization (later in the script)**  
     - Creates publication-quality plots of model results using ggplot2 and ggbeeswarm.  
---

## 5. Software Requirements
- R version: 4.3.1 (2023-06-16), running on macOS (aarch64-apple-darwin20)
- Key R packages:  
  - lme4  
  - lmerTest  
  - carData, car  
  - emmeans, multcomp  
  - effects, Hmisc  
  - ggplot2, ggbeeswarm  

---

## 6. Instructions
1. Install R (≥ 4.3.1) and required packages.  
2. Place `Longevity_Data.csv` in the `data/` folder.  
3. Open `Longevity_Rscript.R` and run from top to bottom.  

---

## 7. License & Availability
Data and code are released under the **CC-BY 4.0 license** (unless otherwise required by the journal).  

---

## 8. Contact
For questions about the data or analysis, please contact:  
Ekta Kochar (ekta.kochar.1502@gmail.com) 
