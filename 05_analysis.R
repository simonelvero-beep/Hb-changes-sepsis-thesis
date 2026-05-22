# ================= 05_analysis.R ===================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")

# Packages ----------------------------------------------------------------
packages <- c("dplyr", "stringr", "lubridate", "lme4",
              "lmerTest", "broom.mixed", "ggplot2", "ggpubr",
              "patchwork", "moments")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Load data ---------------------------------------------------------------
sepsis_hb_48h     <- readRDS("data_clean/sepsis_hb_48h.rds")      # one row per episode
sepsis_hb_long    <- readRDS("data_clean/sepsis_hb_long.rds")     # repeated measures per episode
comparison_hb_48h <- readRDS("data_clean/comparison_hb_48h.rds")  # one row per episode

# Helper functions --------------------------------------------------------
format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

safe_spearman <- function(data, x, y = "hb_change") {
  df <- data |>
    dplyr::select(all_of(c(y, x))) |>
    filter(!is.na(.data[[y]]), !is.na(.data[[x]]))
  
  if (nrow(df) < 3) {
    return(data.frame(
      variable = x,
      rho      = NA_real_,
      p_value  = NA_real_,
      sig      = NA_character_
    ))
  }
  
  test <- suppressWarnings(cor.test(df[[y]], df[[x]], method = "spearman"))
  
  data.frame(
    variable = x,
    rho      = unname(test$estimate),
    p_value  = test$p.value,
    sig      = case_when(
      is.na(test$p.value)  ~ NA_character_,
      test$p.value < 0.001 ~ "***",
      test$p.value < 0.01  ~ "**",
      test$p.value < 0.05  ~ "*",
      TRUE                 ~ "ns"
    )
  )
}

get_cor_label <- function(cor_table, varname) {
  row <- cor_table |> filter(variable == varname)
  if (nrow(row) == 0 || is.na(row$rho[1])) return("rho = NA")
  paste0("rho = ", round(row$rho[1], 3), ", p ", format_p(row$p_value[1]))
}

# =========================================================================
# RESEARCH QUESTION 1
# Describe Hb change in first 48h in sepsis episodes
# =========================================================================
cat("=== RQ1: Hb change in sepsis episodes ===\n\n")
cat("Episodes:", nrow(sepsis_hb_48h), "\n")
cat("Unique patients:", n_distinct(sepsis_hb_48h$personal_id), "\n\n")

# Normality check ---------------------------------------------------------
cat("Normality check:\n")
print(shapiro.test(na.omit(sepsis_hb_48h$hb_change)))
cat("Skewness:", round(skewness(sepsis_hb_48h$hb_change, na.rm = TRUE), 3), "\n")
cat("Kurtosis:", round(kurtosis(sepsis_hb_48h$hb_change, na.rm = TRUE), 3), "\n")

# Wilcoxon signed rank - baseline vs last ---------------------------------
cat("\nWilcoxon signed rank - baseline vs last Hb:\n")
paired_hb <- sepsis_hb_48h |>
  filter(!is.na(hb_baseline), !is.na(hb_last))

print(wilcox.test(
  paired_hb$hb_baseline,
  paired_hb$hb_last,
  paired = TRUE
))

# Hb change by patient outcome --------------------------------------------
cat("\nHb change by patient outcome:\n")
outcome_data <- sepsis_hb_48h |>
  filter(!is.na(hb_change), !is.na(patient_outcome))

print(wilcox.test(hb_change ~ patient_outcome, data = outcome_data))

sepsis_hb_48h |>
  group_by(patient_outcome) |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# Hb change by measurement type -------------------------------------------
cat("\nHb change by measurement type:\n")
measurement_data <- sepsis_hb_48h |>
  filter(!is.na(hb_change), !is.na(measurement_type))

print(kruskal.test(hb_change ~ measurement_type, data = measurement_data))

cat("\nPost-hoc pairwise comparisons for measurement type (Bonferroni):\n")
print(pairwise.wilcox.test(
  measurement_data$hb_change,
  measurement_data$measurement_type,
  p.adjust.method = "bonferroni"
))

measurement_data |>
  group_by(measurement_type) |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# Correlation - Hb change vs ICU stay -------------------------------------
cat("\nCorrelation - Hb change vs ICU stay:\n")
icu_cor_data <- sepsis_hb_48h |>
  filter(!is.na(hb_change), !is.na(icu_hours))

cor_hb_icu <- cor.test(
  icu_cor_data$hb_change,
  icu_cor_data$icu_hours,
  method = "spearman"
)
print(cor_hb_icu)

# Sensitivity analysis - primary ------------------------------------------
cat("\nSensitivity analysis - primary episode-level dataset:\n")
sepsis_hb_48h |>
  summarise(
    n             = n(),
    n_valid       = sum(!is.na(hb_change)),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1)
  ) |>
  print()

# Sensitivity - exclude transfused episodes -------------------------------
cat("\nSensitivity - excluding transfused episodes (rbc_48h = 0):\n")
sepsis_hb_48h |>
  filter(rbc_48h == 0) |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1)
  ) |>
  print()

# Sensitivity - ICU stay >= 48h -------------------------------------------
cat("\nSensitivity - episodes with ICU stay >= 48h:\n")
sepsis_hb_48h |>
  filter(icu_hours >= 48) |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1)
  ) |>
  print()

# =========================================================================
# RESEARCH QUESTION 2
# Compare Hb change between sepsis vs neurotrauma/SAH
# =========================================================================
cat("\n=== RQ2: Comparison across diagnosis groups ===\n\n")

cat("Descriptive stats by group:\n")
comparison_hb_48h |>
  group_by(diagnosis_group) |>
  summarise(
    n               = n(),
    n_valid         = sum(!is.na(hb_change)),
    median_baseline = round(median(hb_baseline, na.rm = TRUE), 1),
    median_change   = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change      = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased   = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1),
    pct_dead        = round(mean(patient_outcome == "dead", na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) |>
  print(width = Inf)

# Kruskal-Wallis - overall difference between groups ----------------------
cat("\nKruskal-Wallis test:\n")
rq2_data <- comparison_hb_48h |>
  filter(!is.na(hb_change), !is.na(diagnosis_group))
rq2_kw <- kruskal.test(hb_change ~ diagnosis_group, data = rq2_data)
print(rq2_kw)

# Post-hoc pairwise Mann-Whitney ------------------------------------------
cat("\nPost-hoc pairwise comparisons (Bonferroni):\n")
rq2_pairwise <- pairwise.wilcox.test(
  rq2_data$hb_change,
  rq2_data$diagnosis_group,
  p.adjust.method = "bonferroni"
)
print(rq2_pairwise)

# Baseline Hb comparison --------------------------------------------------
cat("\nBaseline Hb comparison across groups:\n")
baseline_group_data <- comparison_hb_48h |>
  filter(!is.na(hb_baseline), !is.na(diagnosis_group))
print(kruskal.test(hb_baseline ~ diagnosis_group, data = baseline_group_data))

# Adjusted group comparison for hb_last -----------------------------------
cat("\nAdjusted group comparison for hb_last:\n")

ancova_data <- comparison_hb_48h |>
  mutate(
    diagnosis_group = factor(
      diagnosis_group,
      levels = c("Neurotrauma", "SAH", "Sepsis")
    )
  ) |>
  filter(
    !is.na(hb_last),
    !is.na(hb_baseline),
    !is.na(age),
    !is.na(sex),
    !is.na(hours_at_last),
    !is.na(diagnosis_group)
  )

cat("Adjusted model dataset size:\n")
print(ancova_data |> count(diagnosis_group))

ancova_model <- lm(
  hb_last ~ diagnosis_group + hb_baseline + age + sex + hours_at_last,
  data = ancova_data
)
print(summary(ancova_model))

# Explicit pairwise contrasts ---------------------------------------------
cat("\nExplicit pairwise contrasts:\n")

ancova_data_sah_ref <- ancova_data |>
  mutate(diagnosis_group = relevel(diagnosis_group, ref = "SAH"))

ancova_model_sah_ref <- lm(
  hb_last ~ diagnosis_group + hb_baseline + age + sex + hours_at_last,
  data = ancova_data_sah_ref
)

coef_main <- summary(ancova_model)$coefficients
coef_sah  <- summary(ancova_model_sah_ref)$coefficients

cat("Sepsis vs SAH (adjusted):\n")
print(round(coef_sah["diagnosis_groupSepsis", ], 4))

cat("\nSummary of adjusted pairwise comparisons:\n")
data.frame(
  Comparison = c("Sepsis vs Neurotrauma",
                 "SAH vs Neurotrauma",
                 "Sepsis vs SAH"),
  Beta = c(
    round(coef_main["diagnosis_groupSepsis", "Estimate"], 2),
    round(coef_main["diagnosis_groupSAH", "Estimate"], 2),
    round(coef_sah["diagnosis_groupSepsis", "Estimate"], 2)
  ),
  p_value = c(
    round(coef_main["diagnosis_groupSepsis", "Pr(>|t|)"], 4),
    round(coef_main["diagnosis_groupSAH", "Pr(>|t|)"], 4),
    round(coef_sah["diagnosis_groupSepsis", "Pr(>|t|)"], 4)
  )
) |>
  print()

# Crude vs adjusted by group ----------------------------------------------
cat("\nCrude vs adjusted Hb change by diagnosis group:\n")
ancova_data |>
  group_by(diagnosis_group) |>
  summarise(
    n                  = n(),
    mean_baseline      = round(mean(hb_baseline, na.rm = TRUE), 1),
    crude_change       = round(mean(hb_change, na.rm = TRUE), 1),
    mean_last          = round(mean(hb_last, na.rm = TRUE), 1),
    mean_hours_at_last = round(mean(hours_at_last, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# Sampling volume comparison ----------------------------------------------
cat("\nSampling volume comparison across groups:\n")
sampling_group_data <- comparison_hb_48h |>
  filter(!is.na(total_volume_ml), !is.na(diagnosis_group))

sampling_group_data |>
  group_by(diagnosis_group) |>
  summarise(
    median_volume = round(median(total_volume_ml, na.rm = TRUE), 1),
    iqr_volume    = round(IQR(total_volume_ml, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

rq2_sampling_kw <- kruskal.test(total_volume_ml ~ diagnosis_group, data = sampling_group_data)
print(rq2_sampling_kw)

# Visualization RQ2 -------------------------------------------------------
comparisons_rq2 <- list(
  c("Sepsis", "Neurotrauma"),
  c("Sepsis", "SAH"),
  c("Neurotrauma", "SAH")
)

p_rq2 <- comparison_hb_48h |>
  filter(!is.na(hb_change), !is.na(diagnosis_group)) |>
  ggplot(aes(x = diagnosis_group, y = hb_change, fill = diagnosis_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2) +
  stat_compare_means(
    comparisons     = comparisons_rq2,
    method          = "wilcox.test",
    p.adjust.method = "bonferroni",
    label           = "p.signif"
  ) +
  scale_fill_manual(values = c("Sepsis"      = "steelblue",
                               "Neurotrauma" = "darkorange",
                               "SAH"         = "forestgreen")) +
  labs(title = "Hb change by diagnosis group",
       x     = "Diagnosis group",
       y     = "Hb change (g/L)") +
  theme_minimal() +
  theme(legend.position = "none")

cat("\nPlot: Hb change by diagnosis group\n")
print(p_rq2)

# =========================================================================
# RESEARCH QUESTION 3
# Association between Hb change and inflammation/other variables
# =========================================================================
cat("\n=== RQ3: Association with inflammation and other variables ===\n\n")

# Spearman correlations ---------------------------------------------------
cat("Spearman correlations with Hb change:\n")
vars <- c("crp", "crp_24h", "lpk", "tpk", "mch", "mcv",
          "kreatinin", "bilirubin", "age", "hb_baseline",
          "fluid_balance_48h", "rbc_48h", "plasma_48h",
          "total_volume_ml", "n_bloodgas_draws", "n_labkemi_draws")

cor_results <- lapply(vars, function(v) safe_spearman(sepsis_hb_48h, v)) |>
  bind_rows() |>
  mutate(
    rho     = round(rho, 3),
    p_value = round(p_value, 4)
  )

print(cor_results)

# Crude univariable approximations ----------------------------------------
cat("\nCrude univariable effect-size approximations (NOT independent effects):\n")
lm_inflammation <- lm(hb_change ~ crp_24h,
                      data = sepsis_hb_48h, na.action = na.omit)
lm_dilution     <- lm(hb_change ~ fluid_balance_48h,
                      data = sepsis_hb_48h, na.action = na.omit)
lm_sampling     <- lm(hb_change ~ total_volume_ml,
                      data = sepsis_hb_48h, na.action = na.omit)

mean_crp_24h       <- mean(sepsis_hb_48h$crp_24h, na.rm = TRUE)
mean_fluid_balance <- mean(sepsis_hb_48h$fluid_balance_48h, na.rm = TRUE)
mean_sampling      <- mean(sepsis_hb_48h$total_volume_ml, na.rm = TRUE)

cat("Mean CRP 24h:", round(mean_crp_24h, 1), "mg/L\n")
cat("Crude univariable Hb association with inflammation:",
    round(coef(lm_inflammation)["crp_24h"] * mean_crp_24h, 1), "g/L\n\n")
cat("Mean fluid balance:", round(mean_fluid_balance, 1), "mL\n")
cat("Crude univariable Hb association with hemodilution:",
    round(coef(lm_dilution)["fluid_balance_48h"] * mean_fluid_balance, 1), "g/L\n\n")
cat("Mean sampling volume:", round(mean_sampling, 1), "mL\n")
cat("Crude univariable Hb association with blood sampling:",
    round(coef(lm_sampling)["total_volume_ml"] * mean_sampling, 1), "g/L\n\n")
cat("Observed median Hb change:",
    round(median(sepsis_hb_48h$hb_change, na.rm = TRUE), 1), "g/L\n")

# Dose-response - blood sampling tertiles ---------------------------------
cat("\nDose-response - sampling tertiles:\n")
sepsis_hb_48h <- sepsis_hb_48h |>
  mutate(
    sampling_tertile = ntile(total_volume_ml, 3),
    sampling_group   = factor(
      sampling_tertile,
      levels = c(1, 2, 3),
      labels = c("Low", "Medium", "High")
    )
  )

sampling_data <- sepsis_hb_48h |>
  filter(!is.na(sampling_group), !is.na(hb_change))

sampling_data |>
  group_by(sampling_group) |>
  summarise(
    n             = n(),
    mean_volume   = round(mean(total_volume_ml, na.rm = TRUE), 1),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) |>
  print()

sampling_kw <- kruskal.test(hb_change ~ sampling_group, data = sampling_data)
print(sampling_kw)

sampling_pairwise <- pairwise.wilcox.test(
  sampling_data$hb_change,
  sampling_data$sampling_group,
  p.adjust.method = "bonferroni"
)
print(sampling_pairwise)

# Confounding check -------------------------------------------------------
cat("\nConfounding check - inflammation and ICU stay across sampling groups:\n")
sampling_conf_data <- sepsis_hb_48h |>
  filter(!is.na(sampling_group))

sampling_conf_data |>
  group_by(sampling_group) |>
  summarise(
    median_crp     = round(median(crp,       na.rm = TRUE), 1),
    median_crp_24h = round(median(crp_24h,   na.rm = TRUE), 1),
    median_icu_hrs = round(median(icu_hours, na.rm = TRUE), 1),
    pct_dead       = round(mean(patient_outcome == "dead", na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) |>
  print()

cat("\nKruskal-Wallis - CRP across sampling groups:\n")
crp_kw <- kruskal.test(crp ~ sampling_group, data = sampling_conf_data)
print(crp_kw)

cat("\nKruskal-Wallis - ICU hours across sampling groups:\n")
icu_sampling_kw <- kruskal.test(icu_hours ~ sampling_group, data = sampling_conf_data)
print(icu_sampling_kw)

cat("\nCorrelation - sampling volume vs ICU hours:\n")
sampling_icu_cor_data <- sepsis_hb_48h |>
  filter(!is.na(total_volume_ml), !is.na(icu_hours))

sampling_icu_cor <- cor.test(
  sampling_icu_cor_data$total_volume_ml,
  sampling_icu_cor_data$icu_hours,
  method = "spearman"
)
print(sampling_icu_cor)

# Adjusted model ----------------------------------------------------------
cat("\nAdjusted model - sampling independent of ICU stay:\n")
model_sampling_adjusted <- lm(
  hb_change ~ total_volume_ml + icu_hours + crp_24h +
    fluid_balance_48h + rbc_48h + hb_baseline +
    age + sex,
  data      = sepsis_hb_48h,
  na.action = na.omit
)
print(summary(model_sampling_adjusted))

cat("\nBeta per mL sampled (adjusted):",
    round(coef(model_sampling_adjusted)["total_volume_ml"], 4), "g/L\n")

# =========================================================================
# PRIMARY LONGITUDINAL MODEL - lagged time-updated exposures
# =========================================================================
cat("\n--- Primary Longitudinal Model (lagged time-updated exposures) ---\n\n")

sepsis_hb_long_complete <- sepsis_hb_long |>
  filter(
    !is.na(crp_24h)                 &
      !is.na(lpk)                     &
      !is.na(tpk)                     &
      !is.na(kreatinin)               &
      !is.na(bilirubin)               &
      !is.na(hb_baseline)             &
      !is.na(cum_total_vol_ml_lag)    &
      !is.na(cum_fluid_balance_lag)   &
      !is.na(cum_rbc_lag)             &
      !is.na(cum_plasma_lag)          &
      !is.na(age)                     &
      !is.na(sex)                     &
      !is.na(hours_c)
  )

cat("Complete case dataset (lagged, centered time):\n")
cat("  Observations:", nrow(sepsis_hb_long_complete), "\n")
cat("  Episodes:    ", n_distinct(sepsis_hb_long_complete$episode_id), "\n")
cat("  Patients:    ", n_distinct(sepsis_hb_long_complete$personal_id), "\n\n")

# Base model --------------------------------------------------------------
model_base <- lmer(
  lab_result ~ hours_c +
    (1 + hours_c | episode_id),
  data    = sepsis_hb_long_complete,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

# Model without sampling --------------------------------------------------
model_no_sampling <- lmer(
  lab_result ~ hours_c +
    crp_24h * hours_c +
    lpk +
    cum_fluid_balance_lag +
    cum_rbc_lag +
    cum_plasma_lag +
    tpk + kreatinin +
    bilirubin + age + sex + hb_baseline +
    (1 + hours_c | episode_id),
  data    = sepsis_hb_long_complete,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

# Primary model -----------------------------------------------------------
model_primary <- lmer(
  lab_result ~ hours_c +
    crp_24h * hours_c +
    lpk +
    cum_fluid_balance_lag +
    cum_rbc_lag +
    cum_plasma_lag +
    cum_total_vol_ml_lag * hours_c +
    tpk + kreatinin +
    bilirubin + age + sex + hb_baseline +
    (1 + hours_c | episode_id),
  data    = sepsis_hb_long_complete,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

# Model comparisons -------------------------------------------------------
cat("Model comparison - base vs no_sampling vs primary:\n")
model_comp <- anova(model_base, model_no_sampling, model_primary)
print(model_comp)

# Fixed effects -----------------------------------------------------------
cat("\nFixed effects - primary model (lagged, centered time):\n")
lmm_results <- tidy(model_primary, effects = "fixed", conf.int = TRUE) |>
  mutate(across(where(is.numeric), ~ round(., 4)))
print(lmm_results, width = Inf)

cat("\nKey terms (main effects interpreted at hour 24):\n")
lmm_results |>
  filter(str_detect(term, "hours_c|cum_total|crp_24h")) |>
  print(width = Inf)

# Marginal effect of sampling at key time points --------------------------
cat("\nMarginal effect of 1 mL additional prior sampling at key time points:\n")
b_vol <- coef(summary(model_primary))["cum_total_vol_ml_lag", "Estimate"]
b_int <- tryCatch(
  coef(summary(model_primary))["hours_c:cum_total_vol_ml_lag", "Estimate"],
  error = function(e)
    coef(summary(model_primary))["cum_total_vol_ml_lag:hours_c", "Estimate"]
)
cat("  At hour  0:", round(b_vol + b_int * (0  - 24), 4), "g/L per mL\n")
cat("  At hour 24:", round(b_vol + b_int * (24 - 24), 4), "g/L per mL\n")
cat("  At hour 48:", round(b_vol + b_int * (48 - 24), 4), "g/L per mL\n")

# Per-draw effect on Hb ---------------------------------------------------
cat("\nPer-draw effect on Hb:\n")
lm_per_draw <- lm(
  hb_change ~ total_draws + icu_hours + crp_24h + hb_baseline + age + sex,
  data      = sepsis_hb_48h,
  na.action = na.omit
)
print(summary(lm_per_draw))

cat("\nEstimated Hb loss per blood draw:",
    round(coef(lm_per_draw)["total_draws"], 3), "g/L\n")
cat("At mean", round(mean(sepsis_hb_48h$total_draws, na.rm = TRUE), 1),
    "draws -> estimated total loss:",
    round(coef(lm_per_draw)["total_draws"] *
            mean(sepsis_hb_48h$total_draws, na.rm = TRUE), 1), "g/L\n")

# Forest plot - primary model fixed effects -------------------------------

# Edit these labels to whatever you want displayed on the y-axis
term_labels <- c(
  "hours_c"                       = "Time (hours, centered at 24h)",
  "crp_24h"                       = "CRP 24-48h (mg/L)",
  "lpk"                           = "Leukocytes (10⁹/L)",
  "cum_fluid_balance_lag"         = "Cumulative fluid balance (mL)",
  "cum_rbc_lag"                   = "Cumulative RBC transfusion (mL)",
  "cum_plasma_lag"                = "Cumulative plasma transfusion (mL)",
  "cum_total_vol_ml_lag"          = "Cumulative sampling volume (mL)",
  "tpk"                           = "Platelets (10⁹/L)",
  "kreatinin"                     = "Creatinine (µmol/L)",
  "bilirubin"                     = "Bilirubin (µmol/L)",
  "age"                           = "Age (years)",
  "sexmale"                       = "Sex (male)",
  "hb_baseline"                   = "Baseline Hb (g/L)",
  "hours_c:crp_24h"               = "Time × CRP interaction",
  "hours_c:cum_total_vol_ml_lag"  = "Time × Sampling volume interaction"
)

p_forest <- lmm_results |>
  filter(term != "(Intercept)") |>
  mutate(
    significant  = ifelse(p.value < 0.05, "Significant", "Not significant"),
    term         = factor(term, levels = rev(unique(term)))
  ) |>
  ggplot(aes(x = estimate, y = term, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Significant"     = "steelblue",
                                "Not significant" = "grey60")) +
  scale_y_discrete(labels = term_labels) +   # <-- applies custom labels
  labs(
    title    = "LMM Fixed Effects - Hb trajectory over 48h",
    subtitle = "Primary model: lagged exposures, time centered at 24h",
    x        = "Estimate (95% CI)",
    y        = "",
    color    = ""
  ) +
  theme_minimal()

cat("\nPlot: Forest plot for primary model\n")
print(p_forest)

# Observed vs predicted trajectory ----------------------------------------
sepsis_hb_long_complete$predicted <- predict(model_primary)

p_obs_pred <- ggplot(sepsis_hb_long_complete, aes(x = hours_since_admission)) +
  geom_line(aes(y = lab_result, group = episode_id),
            alpha = 0.1, color = "steelblue") +
  geom_line(aes(y = predicted, group = episode_id),
            alpha = 0.1, color = "red") +
  geom_smooth(aes(y = lab_result), color = "blue", se = TRUE) +
  geom_smooth(aes(y = predicted), color = "darkred", se = TRUE,
              linetype = "dashed") +
  labs(title    = "Hb trajectory over 48h - observed vs predicted",
       subtitle = "Blue = observed, red dashed = primary model predicted",
       x        = "Hours since ICU admission",
       y        = "Hemoglobin (g/L)") +
  theme_minimal()

cat("\nPlot: Observed vs predicted trajectory\n")
print(p_obs_pred)

# Mechanism plots ---------------------------------------------------------
p_inflam <- sepsis_hb_48h |>
  filter(!is.na(crp_24h), !is.na(hb_change)) |>
  ggplot(aes(x = crp_24h, y = hb_change)) +
  geom_point(alpha = 0.3, color = "tomato") +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  labs(title    = "Inflammation (crude association)",
       subtitle = get_cor_label(cor_results, "crp_24h"),
       x = "CRP at 24-48h (mg/L)", y = "Hb change (g/L)") +
  theme_minimal()

p_dilution <- sepsis_hb_48h |>
  filter(!is.na(fluid_balance_48h), !is.na(hb_change)) |>
  ggplot(aes(x = fluid_balance_48h, y = hb_change)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkblue", se = TRUE) +
  labs(title    = "Hemodilution (crude association)",
       subtitle = get_cor_label(cor_results, "fluid_balance_48h"),
       x = "Fluid balance 48h (mL)", y = "Hb change (g/L)") +
  theme_minimal()

p_sampling <- sepsis_hb_48h |>
  filter(!is.na(total_volume_ml), !is.na(hb_change)) |>
  ggplot(aes(x = total_volume_ml, y = hb_change)) +
  geom_point(alpha = 0.3, color = "forestgreen") +
  geom_smooth(method = "lm", color = "darkgreen", se = TRUE) +
  labs(title    = "Blood sampling (crude association)",
       subtitle = get_cor_label(cor_results, "total_volume_ml"),
       x = "Total sampling volume (mL)", y = "Hb change (g/L)") +
  theme_minimal()

p_dose <- ggplot(
  sepsis_hb_48h |> filter(!is.na(sampling_group), !is.na(hb_change)),
  aes(x = sampling_group, y = hb_change, fill = sampling_group)
) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2) +
  stat_compare_means(
    comparisons     = list(c("Low", "Medium"),
                           c("Low", "High"),
                           c("Medium", "High")),
    method          = "wilcox.test",
    p.adjust.method = "bonferroni",
    label           = "p.signif"
  ) +
  scale_fill_manual(values = c("Low"    = "steelblue",
                               "Medium" = "orange",
                               "High"   = "tomato")) +
  labs(title    = "Dose-response: sampling vs Hb decline",
       subtitle = paste0("Sampling tertiles, Kruskal-Wallis p ",
                         format_p(sampling_kw$p.value)),
       x        = "Sampling group",
       y        = "Hb change (g/L)") +
  theme_minimal() +
  theme(legend.position = "none")

p_traj <- sepsis_hb_long |>
  left_join(
    sepsis_hb_48h |> dplyr::select(episode_id, sampling_group),
    by = "episode_id"
  ) |>
  filter(!is.na(sampling_group), !is.na(lab_result),
         !is.na(hours_since_admission)) |>
  ggplot(aes(x = hours_since_admission, y = lab_result,
             color = sampling_group)) +
  geom_smooth(se = TRUE) +
  scale_color_manual(values = c("Low"    = "steelblue",
                                "Medium" = "orange",
                                "High"   = "tomato")) +
  labs(title = "Hb trajectory by sampling group",
       x     = "Hours since admission",
       y     = "Hemoglobin (g/L)",
       color = "Sampling group") +
  theme_minimal()

p_mechanisms <- (p_inflam + p_dilution + p_sampling) / (p_dose + p_traj)

cat("\nPlot: Mechanism plots\n")
print(p_mechanisms)

# =========================================================================
# RQ3 CONTRAST ANALYSIS
# Does blood sampling associate with Hb decline in neuro/SAH?
# =========================================================================
cat("\n=== RQ3 Contrast: Sampling vs Hb decline in neuro/SAH ===\n\n")

neuro_sah_rq3 <- comparison_hb_48h |>
  filter(
    diagnosis_group != "Sepsis",
    !is.na(hb_change),
    !is.na(total_volume_ml)
  )

cat("Neuro/SAH contrast dataset:\n")
cat("  Episodes:", nrow(neuro_sah_rq3), "\n")
print(neuro_sah_rq3 |> count(diagnosis_group))

cat("\nDescriptive summary - neuro/SAH:\n")
neuro_sah_rq3 |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1),
    median_volume = round(median(total_volume_ml, na.rm = TRUE), 1),
    iqr_volume    = round(IQR(total_volume_ml, na.rm = TRUE), 1)
  ) |>
  print()

cat("\nSpearman correlation - sampling vs Hb change (neuro/SAH):\n")
cor_neuro_sampling <- cor.test(
  neuro_sah_rq3$total_volume_ml,
  neuro_sah_rq3$hb_change,
  method = "spearman"
)
print(cor_neuro_sampling)

cat("\nSpearman correlations with Hb change - neuro/SAH group:\n")
vars_neuro <- c("crp", "crp_24h", "lpk", "tpk", "fluid_balance_48h",
                "rbc_48h", "total_volume_ml", "hb_baseline",
                "age", "icu_hours")

cor_results_neuro <- lapply(vars_neuro, function(v) {
  safe_spearman(neuro_sah_rq3, v)
}) |>
  bind_rows() |>
  mutate(rho = round(rho, 3), p_value = round(p_value, 4))

print(cor_results_neuro)

cat("\nAdjusted model - sampling in neuro/SAH:\n")
lm_neuro_sampling <- lm(
  hb_change ~ total_volume_ml + icu_hours + hb_baseline +
    crp_24h + fluid_balance_48h + rbc_48h + age + sex,
  data      = neuro_sah_rq3,
  na.action = na.omit
)
print(summary(lm_neuro_sampling))

lm_sepsis_comparable <- lm(
  hb_change ~ total_volume_ml + icu_hours + hb_baseline +
    crp_24h + fluid_balance_48h + rbc_48h + age + sex,
  data      = sepsis_hb_48h,
  na.action = na.omit
)

cat("\nSide-by-side: adjusted sampling effect in sepsis vs neuro/SAH:\n")
data.frame(
  group          = c("Sepsis", "Neurotrauma/SAH"),
  n              = c(nrow(model.frame(lm_sepsis_comparable)),
                     nrow(model.frame(lm_neuro_sampling))),
  beta_per_ml    = c(round(coef(lm_sepsis_comparable)["total_volume_ml"], 4),
                     round(coef(lm_neuro_sampling)["total_volume_ml"], 4)),
  p_value        = c(
    format_p(summary(lm_sepsis_comparable)$coefficients["total_volume_ml", "Pr(>|t|)"]),
    format_p(summary(lm_neuro_sampling)$coefficients["total_volume_ml", "Pr(>|t|)"])),
  mean_volume_ml = c(round(mean(sepsis_hb_48h$total_volume_ml, na.rm = TRUE), 1),
                     round(mean(neuro_sah_rq3$total_volume_ml, na.rm = TRUE), 1))
) |>
  print()

p_sampling_neuro <- neuro_sah_rq3 |>
  filter(!is.na(total_volume_ml), !is.na(hb_change)) |>
  ggplot(aes(x = total_volume_ml, y = hb_change, color = diagnosis_group)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", se = TRUE) +
  scale_color_manual(values = c("Neurotrauma" = "darkorange",
                                "SAH"         = "forestgreen")) +
  labs(title    = "Blood sampling vs Hb change - neuro/SAH group",
       subtitle = paste0("rho = ",
                         round(unname(cor_neuro_sampling$estimate), 3),
                         ", p ", format_p(cor_neuro_sampling$p.value)),
       x        = "Total sampling volume (mL)",
       y        = "Hb change (g/L)",
       color    = "Diagnosis group") +
  theme_minimal()

p_sampling_sepsis <- sepsis_hb_48h |>
  filter(!is.na(total_volume_ml), !is.na(hb_change)) |>
  ggplot(aes(x = total_volume_ml, y = hb_change)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkblue", se = TRUE) +
  labs(title    = "Blood sampling vs Hb change - sepsis",
       subtitle = get_cor_label(cor_results, "total_volume_ml"),
       x        = "Total sampling volume (mL)",
       y        = "Hb change (g/L)") +
  theme_minimal()

p_contrast <- p_sampling_sepsis + p_sampling_neuro +
  plot_annotation(title = "Sampling vs Hb decline: sepsis vs neuro/SAH")

cat("\nPlot: Contrast - sampling vs Hb change by group\n")
print(p_contrast)

# =========================================================================
# RQ3 EXTENSION - Does diagnosis group modify the sampling-Hb association?
# =========================================================================
cat("\n=== RQ3 Extension: Diagnosis x Sampling interaction ===\n\n")

interaction_data <- comparison_hb_48h |>
  mutate(
    diagnosis_group = factor(diagnosis_group,
                             levels = c("Sepsis", "Neurotrauma", "SAH"))
  ) |>
  filter(
    !is.na(hb_change),     !is.na(total_volume_ml),
    !is.na(hb_baseline),   !is.na(icu_hours),
    !is.na(age),           !is.na(sex),
    !is.na(fluid_balance_48h), !is.na(rbc_48h),
    !is.na(crp_24h)
  )

cat("Dataset size:\n")
print(interaction_data |> count(diagnosis_group))
cat("\n")

model_diag_main <- lm(
  hb_change ~ total_volume_ml + diagnosis_group +
    icu_hours + hb_baseline + age + sex +
    fluid_balance_48h + rbc_48h + crp_24h,
  data = interaction_data
)

model_diag_interaction <- lm(
  hb_change ~ total_volume_ml * diagnosis_group +
    icu_hours + hb_baseline + age + sex +
    fluid_balance_48h + rbc_48h + crp_24h,
  data = interaction_data
)

cat("LRT - does diagnosis modify the sampling effect?\n")
print(anova(model_diag_main, model_diag_interaction))

cat("\nInteraction terms:\n")
summary(model_diag_interaction)$coefficients |>
  as.data.frame() |>
  tibble::rownames_to_column("term") |>
  filter(stringr::str_detect(term, "volume|diagnosis")) |>
  mutate(across(where(is.numeric), ~ round(., 4))) |>
  print()

cat("\nSampling effect per diagnosis group:\n")
coefs <- coef(model_diag_interaction)
data.frame(
  group          = c("Sepsis", "Neurotrauma", "SAH"),
  beta_per_ml    = round(c(
    coefs["total_volume_ml"],
    coefs["total_volume_ml"] + coefs["total_volume_ml:diagnosis_groupNeurotrauma"],
    coefs["total_volume_ml"] + coefs["total_volume_ml:diagnosis_groupSAH"]
  ), 4),
  mean_volume_ml = round(c(
    mean(interaction_data$total_volume_ml[interaction_data$diagnosis_group == "Sepsis"],      na.rm = TRUE),
    mean(interaction_data$total_volume_ml[interaction_data$diagnosis_group == "Neurotrauma"], na.rm = TRUE),
    mean(interaction_data$total_volume_ml[interaction_data$diagnosis_group == "SAH"],         na.rm = TRUE)
  ), 1)
) |>
  mutate(estimated_loss = round(beta_per_ml * mean_volume_ml, 1)) |>
  print()

p_interaction <- interaction_data |>
  ggplot(aes(x = total_volume_ml, y = hb_change, color = diagnosis_group)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("Sepsis"      = "steelblue",
                                "Neurotrauma" = "darkorange",
                                "SAH"         = "forestgreen")) +
  labs(
    title    = "Sampling volume vs Hb change by diagnosis group",
    subtitle = paste0("Interaction p = ",
                      format_p(anova(model_diag_main,
                                     model_diag_interaction)$`Pr(>F)`[2])),
    x        = "Total sampling volume (mL)",
    y        = "Hb change (g/L)",
    color    = "Diagnosis group"
  ) +
  theme_minimal()

cat("\nPlot: Sampling x diagnosis interaction\n")
print(p_interaction)

# =========================================================================
# RESULTS SUMMARY — plain numbers by RQ
# =========================================================================
rho_sampling   <- cor_results |> filter(variable == "total_volume_ml")  |> pull(rho)
p_sampling_cor <- cor_results |> filter(variable == "total_volume_ml")  |> pull(p_value)
rho_fluid      <- cor_results |> filter(variable == "fluid_balance_48h")|> pull(rho)
p_fluid_cor    <- cor_results |> filter(variable == "fluid_balance_48h")|> pull(p_value)
rho_crp24      <- cor_results |> filter(variable == "crp_24h")          |> pull(rho)
p_crp24_cor    <- cor_results |> filter(variable == "crp_24h")          |> pull(p_value)
rho_baseline   <- cor_results |> filter(variable == "hb_baseline")      |> pull(rho)
p_baseline_cor <- cor_results |> filter(variable == "hb_baseline")      |> pull(p_value)

p_sampling_adj <- summary(model_sampling_adjusted)$coefficients["total_volume_ml", "Pr(>|t|)"]
p_add_sampling <- tail(na.omit(model_comp$`Pr(>Chisq)`), 1)
p_interaction  <- anova(model_diag_main, model_diag_interaction)$`Pr(>F)`[2]

cat("\n=== RESULTS SUMMARY ===\n\n")

# RQ1 ---------------------------------------------------------------------
cat("--- RQ1: Hb change in sepsis (n =", nrow(sepsis_hb_48h), "episodes) ---\n")
cat("Median Hb change (g/L):     ", round(median(sepsis_hb_48h$hb_change, na.rm = TRUE), 1), "\n")
cat("IQR (g/L):                  ", round(IQR(sepsis_hb_48h$hb_change, na.rm = TRUE), 1), "\n")
cat("% with Hb decrease:         ", round(mean(sepsis_hb_48h$hb_change < 0, na.rm = TRUE) * 100, 1), "%\n")
cat("Wilcoxon signed rank p:     <0.001\n")
cat("Sensitivity - no transfusion:\n")
sepsis_hb_48h |>
  filter(rbc_48h == 0) |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1)
  ) |> print()
cat("Sensitivity - ICU >= 48h:\n")
sepsis_hb_48h |>
  filter(icu_hours >= 48) |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1)
  ) |> print()

# RQ2 ---------------------------------------------------------------------
cat("\n--- RQ2: Group comparison ---\n")
comparison_hb_48h |>
  group_by(diagnosis_group) |>
  summarise(
    n              = n(),
    median_baseline = round(median(hb_baseline, na.rm = TRUE), 1),
    median_change  = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change     = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased  = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1),
    median_volume  = round(median(total_volume_ml, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print(width = Inf)

cat("Kruskal-Wallis Hb change p:            ", format_p(rq2_kw$p.value), "\n")
cat("ANCOVA Sepsis vs Neurotrauma: beta =",
    round(coef_main["diagnosis_groupSepsis", "Estimate"], 2),
    ", p", format_p(coef_main["diagnosis_groupSepsis", "Pr(>|t|)"]), "\n")
cat("ANCOVA Sepsis vs SAH:         beta =",
    round(coef_sah["diagnosis_groupSepsis", "Estimate"], 2),
    ", p", format_p(coef_sah["diagnosis_groupSepsis", "Pr(>|t|)"]), "\n")
cat("ANCOVA SAH vs Neurotrauma:    beta =",
    round(coef_main["diagnosis_groupSAH", "Estimate"], 2),
    ", p", format_p(coef_main["diagnosis_groupSAH", "Pr(>|t|)"]), "\n")
cat("Sampling volume KW p:                  ", format_p(rq2_sampling_kw$p.value), "\n")

# RQ3 ---------------------------------------------------------------------
cat("\n--- RQ3: Mechanisms ---\n")
cat("Spearman correlations with Hb change:\n")
cat("  hb_baseline:       rho =", round(rho_baseline, 3),  ", p", format_p(p_baseline_cor), "\n")
cat("  fluid_balance_48h: rho =", round(rho_fluid, 3),     ", p", format_p(p_fluid_cor), "\n")
cat("  crp_24h:           rho =", round(rho_crp24, 3),     ", p", format_p(p_crp24_cor), "\n")
cat("  total_volume_ml:   rho =", round(rho_sampling, 3),  ", p", format_p(p_sampling_cor), "\n\n")

cat("Dose-response (sampling tertiles):\n")
sampling_data |>
  group_by(sampling_group) |>
  summarise(
    n             = n(),
    mean_volume   = round(mean(total_volume_ml, na.rm = TRUE), 1),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    .groups = "drop"
  ) |> print()
cat("Kruskal-Wallis p:", format_p(sampling_kw$p.value), "\n\n")

cat("Multivariable model (adjusted):\n")
cat("  Sampling beta per mL:", round(coef(model_sampling_adjusted)["total_volume_ml"], 4),
    "g/L, p", format_p(p_sampling_adj), "\n\n")

cat("LMM primary model:\n")
cat("  cum_total_vol_ml_lag at hour 24:", round(b_vol, 4), "g/L per mL,",
    "p", format_p(lmm_results |> filter(term == "cum_total_vol_ml_lag") |> pull(p.value)), "\n")
cat("  LR test adding sampling:        p", format_p(p_add_sampling), "\n")
cat("  CRP x time interaction:         p",
    format_p(lmm_results |> filter(term == "hours_c:crp_24h") |> pull(p.value)), "\n\n")

cat("Per-draw effect:\n")
cat("  Beta per draw:", round(coef(lm_per_draw)["total_draws"], 3),
    "g/L, p", format_p(summary(lm_per_draw)$coefficients["total_draws", "Pr(>|t|)"]), "\n\n")

cat("Contrast analysis (neuro/SAH):\n")
cat("  Spearman rho:", round(unname(cor_neuro_sampling$estimate), 3),
    ", p", format_p(cor_neuro_sampling$p.value), "\n")
cat("  Adjusted beta per mL:",
    round(coef(lm_neuro_sampling)["total_volume_ml"], 4),
    "g/L, p",
    format_p(summary(lm_neuro_sampling)$coefficients["total_volume_ml", "Pr(>|t|)"]), "\n\n")

cat("Diagnosis x sampling interaction:\n")
cat("  F-test p:", format_p(p_interaction), "\n")
cat("  Sepsis beta:      ", round(coefs["total_volume_ml"], 4), "g/L per mL\n")
cat("  Neurotrauma beta: ", round(coefs["total_volume_ml"] +
                                    coefs["total_volume_ml:diagnosis_groupNeurotrauma"], 4), "g/L per mL\n")
cat("  SAH beta:         ", round(coefs["total_volume_ml"] +
                                    coefs["total_volume_ml:diagnosis_groupSAH"], 4), "g/L per mL\n")

# =========================================================================
# EXPORT ALL FIGURES AS PNG
# =========================================================================

# ── RQ1 figures ───────────────────────────────────────────────────────────

# Combined RQ1 overview
ggsave("figures/fig_rq1_combined.png",
       plot   = (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8),
       width  = 14, height = 20, dpi = 300, bg = "white")
cat("Saved: figures/fig_rq1_combined.png\n")

# ── RQ2 figures ───────────────────────────────────────────────────────────

ggsave("figures/fig_rq2_hb_change.png",
       plot   = p_rq2,
       width  = 8, height = 6, dpi = 300, bg = "white")
cat("Saved: figures/fig_rq2_hb_change.png\n")

# ── RQ3 figures ───────────────────────────────────────────────────────────

ggsave("figures/fig_rq3_mechanisms.png",
       plot   = p_mechanisms,
       width  = 14, height = 10, dpi = 300, bg = "white")
cat("Saved: figures/fig_rq3_mechanisms.png\n")

ggsave("figures/fig_rq3_forest.png",
       plot   = p_forest,
       width  = 10, height = 7, dpi = 300, bg = "white")
cat("Saved: figures/fig_rq3_forest.png\n")

ggsave("figures/fig_rq3_obs_pred.png",
       plot   = p_obs_pred,
       width  = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: figures/fig_rq3_obs_pred.png\n")

# ── Contrast figures ──────────────────────────────────────────────────────

ggsave("figures/fig_rq3_contrast.png",
       plot   = p_contrast,
       width  = 12, height = 6, dpi = 300, bg = "white")
cat("Saved: figures/fig_rq3_contrast.png\n")

ggsave("figures/fig_rq3_interaction.png",
       plot   = p_interaction,
       width  = 8, height = 6, dpi = 300, bg = "white")
cat("Saved: figures/fig_rq3_interaction.png\n")

cat("\nAll figures exported to figures/\n")

cat("\n✓ Analysis complete\n")
cat("Pipeline complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")