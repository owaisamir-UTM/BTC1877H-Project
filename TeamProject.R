# BTC1877H - Team Project
# Team 1
# Agata Wolochacz, Avery Fitzpatrick, Aya Finkelstein, Nina He, Owais Amir

################################################################################
# Package installation
################################################################################

library(dplyr)
library(janitor)
library(lubridate)
library(ggplot2)
library(survival)
library(glmnet)
library(table1)
library(tidyr)
library(broom)
library(pROC)
library(tree)
library(gt)

################################################################################
# Reading in data and copying
################################################################################

# Reading data in
data <- readxl::read_xlsx("transfusion data.xlsx") |> clean_names()

# Defining which subset of columns to keep
cols_to_keep <- c(
  "study_id_number", "or_date", "death_date", # Identifiers
  "age", "height", "weight", "gender_male",   # Demographics
  "type", "redo_lung_transplant", "intraoperative_ecls", "ecls_cpb", # Risk Factors
  "cystic_fibrosis", "idiopathic_pulmonary_hypertension",  # Diagnosis
  "las_score", "pre_hb", "pre_platelets", "pre_inr",       # Labs
  "total_24hr_rbc", "duration_of_icu_stay_days", "alive_12mths_yn", "massive_transfusion" # Outcomes
)

# Copying data to a new clean dataframe
clean_df <- data[ ,cols_to_keep]

# Cleaning up environment
rm(cols_to_keep)

################################################################################
# Cleaning up data types and setting factors
################################################################################

# Looking at the data structures
str(clean_df)

# Summary statistics of the data frame
summary(clean_df)

# Looking for any columns that have NA values in them
na_check <- apply(clean_df, 2, function(x) any(is.na(x)))

# Cleaning up data types
clean_df$death_date <- as.POSIXct(
  strptime(clean_df$death_date, format = "%d-%b-%Y"),
  tz = "UTC"
)

cols_to_clean <- c("gender_male", "redo_lung_transplant", "ecls_cpb", "cystic_fibrosis", "idiopathic_pulmonary_hypertension")
clean_df <- clean_df |> 
  mutate(across(all_of(cols_to_clean), ~ if_else(. == "FALSE", FALSE, TRUE)))

clean_df <- clean_df |> rename("gender" = "gender_male")
clean_df$gender <- factor(clean_df$gender,
                          levels = c(FALSE, TRUE),
                          labels = c("Female", "Male"))

clean_df$massive_transfusion <- if_else(clean_df$massive_transfusion == 1, TRUE, FALSE)

# Collapsing "Single Left Lung" and "Single Right Lung" as just "Single"
clean_df$type <- if_else(clean_df$type == "Bilateral", "Bilateral", "Single")
clean_df$type <- factor(clean_df$type, levels = c("Single", "Bilateral"))

clean_df$alive_12mths_yn <- if_else(clean_df$alive_12mths_yn == "Y", TRUE, FALSE)
clean_df <- clean_df |> rename("alive_12mths" = "alive_12mths_yn")

summary(clean_df)

# Cleaning up environment
rm(cols_to_clean, na_check)

################################################################################
# Creating composite variables
################################################################################

# Event variable where event (1= died, 0=alive/censored)
clean_df["event"] <- if_else(clean_df$alive_12mths == "TRUE", 0, 1)

# Calculating days alive post OR date up until censoring
clean_df["survival_time_raw"] <- as.numeric(clean_df$death_date - clean_df$or_date)
clean_df["survival_time"] <- if_else(clean_df$event == 1 & clean_df$survival_time_raw <= 365,
                                     clean_df$survival_time_raw, 365)

# Calculating BMI from weight and height
clean_df["bmi"] <- clean_df$weight/((clean_df$height/100)**2)

# Classifying a high risk patient 
clean_df["high_risk_patient"] <- if_else(clean_df$redo_lung_transplant == TRUE |
                                           clean_df$intraoperative_ecls == TRUE |
                                           clean_df$ecls_cpb == TRUE |
                                           clean_df$cystic_fibrosis == TRUE |
                                           clean_df$idiopathic_pulmonary_hypertension == TRUE, TRUE, FALSE)

# Classifying patients more likely to bleed due to coagulopathy
clean_df["baseline_coagulopathy"] <- if_else(clean_df$pre_inr > 1.5 | clean_df$pre_platelets < 150, TRUE, FALSE)

################################################################################
# Table 1 summary statistics
################################################################################

# Variable labels
label(clean_df$age) <- "Age (years)"
label(clean_df$height) <- "Height (cm)"
label(clean_df$weight) <- "Weight (kg)"
label(clean_df$bmi) <- "BMI (kg/m²)"
label(clean_df$gender) <- "Gender"
label(clean_df$type) <- "Transplant type"
label(clean_df$redo_lung_transplant) <- "Redo lung transplant"
label(clean_df$intraoperative_ecls) <- "Intraoperative ECLS"
label(clean_df$ecls_cpb) <- "CPB/ECLS use"
label(clean_df$cystic_fibrosis) <- "Cystic fibrosis"
label(clean_df$idiopathic_pulmonary_hypertension) <- "Idiopathic pulmonary HTN"
label(clean_df$las_score) <- "LAS score"
label(clean_df$pre_hb) <- "Pre-op Hb (g/L)"
label(clean_df$pre_platelets) <- "Pre-op platelets (×10⁹/L)"
label(clean_df$pre_inr) <- "Pre-op INR"
label(clean_df$total_24hr_rbc) <- "RBC transfusion 0–24h (units)"
label(clean_df$duration_of_icu_stay_days) <- "ICU LOS (days)"
label(clean_df$massive_transfusion) <- "Massive transfusion"
label(clean_df$high_risk_patient) <- "High-risk patient"
label(clean_df$baseline_coagulopathy) <- "Baseline coagulopathy"

# Generating Table1 Summary
table1(
  ~ age + height + weight + bmi +
    type +
    redo_lung_transplant + intraoperative_ecls + ecls_cpb +
    cystic_fibrosis + idiopathic_pulmonary_hypertension +
    las_score +
    pre_hb + pre_platelets + pre_inr +
    total_24hr_rbc +
    duration_of_icu_stay_days +
    massive_transfusion +
    high_risk_patient +
    baseline_coagulopathy
  | gender,
  data = clean_df
)

################################################################################
# Sub-Analysis 1.1: The "Screening" Model (Who needs ANY blood?)
################################################################################

# Create the outcome variable for this analysis: Any Transfusion
# 1 = patient received at least 1 unit of RBCs in first 24 hours
# 0 = patient received 0 units of RBCs
clean_df$any_transfusion <- as.numeric(clean_df$total_24hr_rbc > 0)

# Defining the predictor variables for the model
# These were selected based on clinical knowledge from literature review and available data 
predictor_vars <- c("age",                    # Patient age in years
                    "bmi",                    # Body Mass Index
                    "high_risk_patient",      # Composite high-risk indicator
                    "pre_hb",                 # Preoperative hemoglobin (g/L)
                    "baseline_coagulopathy",  # Coagulopathy indicator (INR>1.5 or platelets<150)
                    "type",                   # Transplant type (Single vs Bilateral)
                    "las_score")              # Lung Allocation Score (disease severity)

# Create a model 1 dataset with only predictor variables and outcome
model1_data <- clean_df[, c(predictor_vars, "any_transfusion")]

# Using complete case analysis, so removing rows with missing data
# LASSO regression requires complete data for all variables
model1_data <- model1_data[complete.cases(model1_data), ]

# Prepare predictor matrix (X) and outcome vector (Y) for glmnet
# model.matrix() converts factors to dummy variables and creates design matrix
# [,-1] removes the intercept column (glmnet adds its own intercept)
X <- model.matrix(any_transfusion ~ age + bmi + high_risk_patient + 
                    pre_hb + baseline_coagulopathy + type + las_score, 
                  data = model1_data)[, -1]

# Extract the binary outcome variable
Y <- model1_data$any_transfusion

# Set random seed for reproducibility of cross-validation results
set.seed(2025)

# Perform 10-fold cross-validated LASSO regression
# cv.glmnet() tests different penalty values (lambda) to find optimal model
# optimizing area under the curve (auc)
cv_lasso <- cv.glmnet(X, Y, 
                      family = "binomial", 
                      alpha = 1,
                      nfolds = 10,
                      type.measure = "auc")

# Plot cross-validation results showing model performance vs. lambda
# Save as PDF 
pdf("model1_lasso_cv_plot.pdf", width = 8, height = 6)
plot(cv_lasso)
title("LASSO Cross-Validation for Any Transfusion Model", line = 2.5)
dev.off() # to close the plot file 

# Fit final LASSO model using lambda.1se
# lambda.1se gives simpler model with fewer predictors, while still staying
# within the one standard error of the min cross-validation error 
# Chose simpler model to reduce risk of overfitting, but will confirm with auc
lasso_model <- glmnet(X, Y, 
                      family = "binomial", 
                      alpha = 1,
                      lambda = cv_lasso$lambda.1se)

# Extract coefficients from LASSO model
lasso_coefs <- coef(lasso_model)

# Identify which variables were selected (the non-zero coefficients)
selected_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
selected_coefs <- as.vector(lasso_coefs)[which(lasso_coefs != 0)]

# Print selected features and their LASSO coefficients
# Note: These are penalized coefficients, not the final odds ratios
for(i in 1:length(selected_vars)) {
  print(sprintf("%s: %.4f", selected_vars[i], selected_coefs[i]))
}

# Remove intercept from selected variables list
non_intercept_vars <- selected_vars[selected_vars != "(Intercept)"]

# Check if LASSO selected at least one predictor variable
if(length(non_intercept_vars) > 0) {
  
  # Build formula for standard logistic regression using LASSO-selected variables
  # To get Odds Ratios with CI
  formula_str <- "any_transfusion ~ high_risk_patient + pre_hb"
  
  # Fit standard (non-penalized) logistic regression with selected predictors
  # This gives us maximum likelihood estimates for odds ratios
  final_model <- glm(as.formula(formula_str), 
                     data = model1_data, 
                     family = binomial(link = "logit"))
  
  # Display full model summary including coefficients, SEs, z-values, and p-values
  print(summary(final_model))
  # exp() to converts log-odds to odds ratios, confint() for CIs
  or_table <- exp(cbind(OR = coef(final_model), confint(final_model)))
  print(round(or_table, 3))
  
  # Generate predicted probabilities for each patient
  model1_data$predicted_prob <- predict(final_model, type = "response")
  
  # Calculate Area Under the ROC Curve (AUC)
  # AUC measures the model's ability to discriminate between outcomes
  roc_obj <- roc(model1_data$any_transfusion, model1_data$predicted_prob, quiet = TRUE)
  auc_value <- auc(roc_obj)
  
  # Calculate 95% confidence interval for AUC
  auc_ci <- ci.auc(roc_obj)
  
  # Create and save ROC curve plot
  # This figure shows sensitivity vs. (1-specificity) across all thresholds
  pdf("model1_roc_curve_screening.pdf", width = 8, height = 8)
  plot(roc_obj, 
       main = "ROC Curve: Screening Model for Any Transfusion",
       col = "orchid", 
       lwd = 2,
       print.auc = TRUE,          # Print AUC value on plot
       print.auc.y = 0.4,         # Position of AUC text
       legacy.axes = TRUE)        
  abline(a = 0, b = 1, lty = 2, col = "gray")  # Add diagonal reference line (chance)
  dev.off()
  
  # Create results summary table for report
  results_summary <- data.frame(
    Variable = names(coef(final_model))[-1],          
    Coefficient = round(coef(final_model)[-1], 3),    
    OR = round(or_table[-1, 1], 3),                   
    CI_Lower = round(or_table[-1, 2], 3),             
    CI_Upper = round(or_table[-1, 3], 3),             
    P_Value = round(summary(final_model)$coefficients[-1, 4], 4)  
  )
  
  # Display results summary table
  print(results_summary)
  
  # Save results table to CSV file 
  write.csv(results_summary, "screening_model1_results.csv", row.names = FALSE)
  
} else {
  # Handle case where LASSO selected no variables (all coefficients = 0)
  # This can happen if lambda is too high or predictors have weak signals
  # Would try lambda.min for less aggressive shrinkage 
  print("LASSO selected no variables (all coefficients shrunk to zero).")
}

################################################################################
# Sub-Analysis 1.2: The "Crisis" Model (Who needs MASSIVE blood?)
################################################################################


################################################################################
# Sub-Analysis 1.3: The "Volume" Model (What drives amount?)
################################################################################


################################################################################
# Sub-Analysis 2.1: Survival Analysis (Mortality)
################################################################################


################################################################################
# Sub-Analysis 2.2: Resource Utilization (ICU LOS) 
################################################################################

