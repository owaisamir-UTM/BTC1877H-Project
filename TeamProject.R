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

# Flagging whether someone received a transfusion
clean_df["recieved_transfusion"] <- if_else(clean_df["total_24hr_rbc"] > 0, TRUE, FALSE)

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

