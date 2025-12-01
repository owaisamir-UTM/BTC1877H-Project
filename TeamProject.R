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
library(car)

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

# Flagging whether someone received a transfusion
clean_df["any_transfusion"] <- if_else(clean_df["total_24hr_rbc"] > 0, 1, 0)

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
label(clean_df$bmi) <- "BMI (kg/m²)"
label(clean_df$gender) <- "Gender"
label(clean_df$type) <- "Transplant type"
label(clean_df$redo_lung_transplant) <- "Redo lung transplant"
label(clean_df$intraoperative_ecls) <- "Intraoperative ECLS"
label(clean_df$ecls_cpb) <- "CPB/ECLS use"
label(clean_df$cystic_fibrosis) <- "Cystic fibrosis"
label(clean_df$idiopathic_pulmonary_hypertension) <- "Idiopathic pulmonary HTN"
label(clean_df$las_score) <- "LAS score"
label(clean_df$total_24hr_rbc) <- "RBC transfusion 0–24h (units)"
label(clean_df$duration_of_icu_stay_days) <- "ICU LOS (days)"
label(clean_df$massive_transfusion) <- "Massive transfusion"
label(clean_df$high_risk_patient) <- "High-risk patient"
label(clean_df$baseline_coagulopathy) <- "Baseline coagulopathy"

# Generating Table1 Summary
table1(
  ~ age + bmi +
    type +
    redo_lung_transplant + intraoperative_ecls + ecls_cpb +
    cystic_fibrosis + idiopathic_pulmonary_hypertension +
    las_score +
    total_24hr_rbc +
    duration_of_icu_stay_days +
    massive_transfusion +
    high_risk_patient +
    baseline_coagulopathy
  | factor(any_transfusion, levels = c(0,1), labels = c("No Transfusion", "Recieved Transfusion")),
  data = clean_df,
  caption = "<b>Table 1. Baseline Characteristics of the Study Cohort</b>",
)

################################################################################
# Sub-Analysis 1.1: The "Screening" Model (Who needs ANY blood?)
################################################################################

# Defining the predictor variables for the model
# These were selected based on clinical knowledge from literature review and available data 
predictor_vars <- c("age",                    # Patient age in years
                    "gender",                 # Biological sex
                    "bmi",                    # Body Mass Index
                    "high_risk_patient",      # Composite high-risk indicator
                    "pre_hb",                 # Preoperative hemoglobin (g/L)
                    "baseline_coagulopathy",  # Coagulopathy indicator (INR>1.5 or platelets<150)
                    "type",                   # Transplant type (Single vs Bilateral)
                    "las_score")              # Lung Allocation Score (disease severity)

# Create a model 1 dataset with only predictor variables and outcome
model1_data <- clean_df[, c(predictor_vars, "any_transfusion")]

# Using complete case analysis, so removing rows with missing data
model1_data <- model1_data[complete.cases(model1_data), ]

# Check sample size adequacy (EPV should be ≥10)
# Since we have 8 predictors need at least 80 (8x10) transfusion events total 
n_events <- sum(model1_data$any_transfusion)
n_predictors <- length(predictor_vars)
epv <- n_events / n_predictors
print(epv)

# Fit multivariable logistic regression model
screening_model <- glm(any_transfusion ~ age + gender + bmi + high_risk_patient + 
                         pre_hb + baseline_coagulopathy + type + las_score,
                       data = model1_data,
                       family = binomial(link = "logit"))

# Display full model summary including coefficients, SEs, z-values, and p-values
print(summary(screening_model))

# Check for multicollinearity using VIF
# Should all be below 5 
vif_values <- vif(screening_model)
print(vif_values)

# Calculate odds ratios with 95% confidence intervals
# exp() converts log-odds to odds ratios, confint() for CIs
or_table <- exp(cbind(OR = coef(screening_model), confint(screening_model)))
print(round(or_table, 3))

# Generate predicted probabilities for each patient
# This is used to calculate AUC, to measure overall model discrimination
model1_data$predicted_prob <- predict(screening_model, type = "response")

# Calculate Area Under the ROC Curve (AUC)
# AUC measures the model's ability to discriminate between transfused 
# vs non-transfused patients
roc_obj <- roc(model1_data$any_transfusion, model1_data$predicted_prob, quiet = TRUE)
auc_value <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj) # Calculate 95% confidence interval for AUC

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
  Variable = names(coef(screening_model))[-1],          
  Coefficient = round(coef(screening_model)[-1], 3),    
  OR = round(or_table[-1, 1], 3),                   
  CI_Lower = round(or_table[-1, 2], 3),             
  CI_Upper = round(or_table[-1, 3], 3),             
  P_Value = ifelse(summary(screening_model)$coefficients[-1, 4] < 0.001, 
                   "<0.001", 
                   round(summary(screening_model)$coefficients[-1, 4], 3))
)

# Display results summary table
print(results_summary)

# Save results table to CSV file 
write.csv(results_summary, "screening_model1_results.csv", row.names = FALSE)

################################################################################
# Sub-Analysis 1.2: The "Crisis" Model (Who needs MASSIVE blood?)
################################################################################

# Defining the logistic regression model for the Crisis Model (Sub-Analysis 1.2)
logreg_full_model <- glm(massive_transfusion ~ age + gender + type + bmi + high_risk_patient + 
                           pre_hb + baseline_coagulopathy + las_score, data = clean_df, family = binomial)

# Stepwise selection for the Crisis Model
logreg_stepwise <- stepAIC(logreg_full_model, direction = "both")

# Summary of the stepwise-selected Logistic Regression Model for the Crisis Model
summary(logreg_stepwise)


################################################################################
# Sub-Analysis 1.3: The "Volume" Model (What drives amount?)
################################################################################

# Linear regression for the Volume Model with all predictors
linreg_full_model <- lm(total_24hr_rbc ~ age + gender + type + bmi + high_risk_patient + 
                          pre_hb + baseline_coagulopathy + las_score, data = clean_df)

# Stepwise selection for the Volume Model
linreg_stepwise <- stepAIC(linreg_full_model, direction = "both")

# Summary of the stepwise-selected Linear Regression Model for the Volume Model
summary(linreg_stepwise)

################################################################################
# Sub-Analysis 2.1: Survival Analysis (Mortality)
################################################################################

                  # Kaplan-Meier estimator for the whole cohort
km_overall <- survfit(Surv(survival_time, event == 1) ~ 1, data = clean_df)

print(km_overall)
plot(km_overall, xlab = "Days from Transplant", ylab = "Survival", ylim = c(0.75, 1.0))

# 1-year survival probability for the whole cohort
summary(km_overall, time = 365)

# Stratified Kaplan-Meier curve by transfusion level

# Creating transfusion groups for visualization
clean_df$rbc_group <- cut(
  clean_df$total_24hr_rbc,
  breaks = c(-Inf, 0, 5, 10, Inf),
  labels = c("None (0)", "Low (1-5)", "Moderate (6-10)", "High (>10)")
)

# KM curve stratified by RBC transfusion group
km_by_rbc <- survfit(Surv(survival_time, event == 1) ~ rbc_group, data = clean_df)

print(km_by_rbc)

# Stratified KM plot
ggsurvplot(
  km_by_rbc,
  data = clean_df,
  pval = TRUE,
  xlab = "Days After Transplant",
  ylab = "Survival Probability",
  legend.title = "RBC Transfusion",
  legend.labs = c("None (0)", "Low (1-5)", "Moderate (6-10)", "High (>10)"),
  ylim = c(0.75, 1.0),
  palette = c("cyan", "darkorange", "green", "purple"),
  ggtheme = theme_bw()
)

# 1-year survival probability per transfusion group
summary(km_by_rbc, time = 365)

# Checking PH assumption visually before log-rank test

# Complementary log-log plot
plot(survfit(Surv(survival_time, event == 1) ~ rbc_group, data = clean_df), 
     fun = "cloglog",
     xlab = "log(Time in Days)",
     ylab = "log(-log(Survival))")
# Only 3 groups displayed because 0 deaths in High(>10) group

table(clean_df$rbc_group)
table(clean_df$rbc_group, clean_df$event)

# Log-rank test comparing survival between transfusion groups

survdiff(Surv(survival_time, event == 1) ~ rbc_group, data = clean_df)

# Cox Proportional Hazards model

# Filtering to complete cases (excluding missing LAS scores)
cox_df <- clean_df |> filter(!is.na(las_score))

# Fitting Cox model with RBC transfusion as primary predictor
cox_model <- coxph(
  Surv(survival_time, event == 1) ~ total_24hr_rbc + age + las_score + high_risk_patient + gender,
  data = cox_df
)

summary(cox_model)

# Testing proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)

# Schoenfeld residual plots for each variable
par(mfrow = c(2, 3))

plot(ph_test, var = 1, main = "RBC Transfusion")
plot(ph_test, var = 2, main = "Age")
plot(ph_test, var = 3, main = "LAS Score")
plot(ph_test, var = 4, main = "High Risk Patient")
plot(ph_test, var = 5, main = "Gender")

# Reset to single plot
par(mfrow = c(1, 1))

# Extracting results for reporting

# Hazard ratios with 95% confidence intervals
cox_results <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(
    HR = round(estimate, 2),
    CI_lower = round(conf.low, 2),
    CI_upper = round(conf.high, 2),
    p_value = round(p.value, 3)
  ) |>
  select(term, HR, CI_lower, CI_upper, p_value)

print(cox_results)

# C-index for model discrimination
c_index <- summary(cox_model)$concordance[1]
print(c_index)

# Visualizations for reporting

# Forest plot hazard data
forest_data <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(term = c("RBC Transfusion (per unit)", "Age (per year)", "LAS Score (per point)", "High Risk Patient", "Gender (Male)"))

# Plot
ggplot(forest_data, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
  ) +
  theme_bw()


# Descriptive statistics by survival status
clean_df |>
  filter(!is.na(las_score)) |>
  group_by(event) |>
  summarise(
    n = n(),
    mean_rbc = round(mean(total_24hr_rbc, na.rm = TRUE), 2),
    median_rbc = median(total_24hr_rbc, na.rm = TRUE),
    mean_age = round(mean(age, na.rm = TRUE), 2),
    mean_las = round(mean(las_score, na.rm = TRUE), 2),
    pct_high_risk = round(mean(high_risk_patient, na.rm = TRUE) * 100, 1),
    pct_male = round(mean(gender == "Male", na.rm = TRUE) * 100, 1)
  )

################################################################################
# Sub-Analysis 2.2: Resource Utilization (ICU LOS) 
################################################################################
# Multiple Linear Regression with Days in ICU as outcome, using all predictors

# Defining predictors
predictors <- c("total_24hr_rbc",
                "age",                   
                "gender",            
                "bmi",               
                "high_risk_patient", 
                "pre_hb",                
                "baseline_coagulopathy", 
                "type",                   
                "las_score")

# Subsetting data
ICU_df <- clean_df[, c(predictors, "duration_of_icu_stay_days")]

# Checking for missing values
colSums(is.na(ICU_df))

# Some missing values found, complete case analysis done 
ICU_df_comp <- na.omit(ICU_df)

# Creating formula
ICU_formula <- as.formula(
  paste("duration_of_icu_stay_days ~", paste(predictors, collapse = " + ")))

# Creating ICU outcome model
ICU_model <- lm(ICU_formula, data = ICU_df_comp)

# Investigating results
summary(ICU_model)

# Multicolinearity Check
vif(ICU_model)                  # all VIFs low
                                
# Obtaining 95% CIs
confint(ICU_model, level = 0.95)

# Making some plots 
plot(predict(ICU_model), ICU_df_comp$duration_of_icu_stay_days,
     xlab = "Predicted ICU LOS",
     ylab = "Observed ICU LOS",
     main = "Observed vs. Predicted ICU LOS")

abline(0, 1, col = "red", lwd = 2)

legend("right",
       legend = c("Perfect Fit Line"),
       col = "red",
       lty = 1,
       lwd = 2,
       bty = "n")

par(mfrow = c(1, 2))
hist(ICU_df_comp$duration_of_icu_stay_days, main = "Distribution of Length of of Stay in ICU",
     xlab = "Length of ICU stay (in days)", ylab = "Number of Patients")

hist(ICU_df_comp$total_24hr_rbc, main = "Distribution of RBC Units Transfused in 24hrs",
     xlab = "Units of RBC Transfused", ylab = "Number of Patients")








