# ================== 04_descriptive_stat.R =========================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")

# Packages ----------------------------------------------------------------
packages <- c("dplyr", "ggplot2", "patchwork", "moments",
              "tidyr", "gtsummary", "stringr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Load merged data --------------------------------------------------------
sepsis_patients   <- readRDS("data_clean/sepsis_patients.rds")      # episode-level
sepsis_hb         <- readRDS("data_clean/sepsis_hb.rds")            # episode-level Hb on ICU
sepsis_hb_48h     <- readRDS("data_clean/sepsis_hb_48h.rds")        # one row per episode
sepsis_hb_long    <- readRDS("data_clean/sepsis_hb_long.rds")       # longitudinal, 0-48h
comparison_hb_48h <- readRDS("data_clean/comparison_hb_48h.rds")    # episode-level comparison



# =========================================================================
# RESEARCH QUESTION 1 - Hb change in sepsis episodes first 48h
# =========================================================================
cat("=== RQ1: Hb change in sepsis episodes ===\n\n")

cat("Sepsis dataset size:\n")
cat("  Episodes:        ", nrow(sepsis_hb_48h), "\n")
cat("  Unique patients: ", n_distinct(sepsis_hb_48h$personal_id), "\n\n")

# Summary table -----------------------------------------------------------
cat("\nSummary - Hb change:\n")
sepsis_hb_48h |>
  summarise(
    n_episodes        = n(),
    n_valid_48h       = sum(!is.na(hb_48h)),
    n_missing_48h     = sum(is.na(hb_48h)),
    mean_baseline     = round(mean(hb_baseline, na.rm = TRUE), 1),
    sd_baseline       = round(sd(hb_baseline, na.rm = TRUE), 1),
    median_baseline   = round(median(hb_baseline, na.rm = TRUE), 1),
    iqr_baseline      = round(IQR(hb_baseline, na.rm = TRUE), 1),
    mean_last         = round(mean(hb_last, na.rm = TRUE), 1),
    sd_last           = round(sd(hb_last, na.rm = TRUE), 1),
    median_last       = round(median(hb_last, na.rm = TRUE), 1),
    iqr_last          = round(IQR(hb_last, na.rm = TRUE), 1),
    median_change     = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change        = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased     = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1),
    mean_change_pct   = round(mean(hb_change_pct, na.rm = TRUE), 1),
    median_change_pct = round(median(hb_change_pct, na.rm = TRUE), 1),
    iqr_change_pct    = round(IQR(hb_change_pct, na.rm = TRUE), 1)
  ) |>
  print(width = Inf)

# Outcome distribution - analysis cohort (n = 859)
cat("\nEpisode outcome distribution (analysis cohort):\n")
sepsis_hb_48h |>
  count(patient_outcome) |>
  mutate(pct = round(n / sum(n) * 100, 1)) |>
  print()

# Subgroup by measurement type --------------------------------------------
cat("\nHb change by measurement type:\n")
sepsis_hb_48h |>
  group_by(measurement_type) |>
  summarise(
    n             = n(),
    median_change = round(median(hb_change, na.rm = TRUE), 1),
    iqr_change    = round(IQR(hb_change, na.rm = TRUE), 1),
    pct_decreased = round(mean(hb_change < 0, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) |>
  print()

# Missing 48h by ICU stay -------------------------------------------------
cat("\nMissing 48h measurement by ICU stay:\n")
sepsis_hb_48h |>
  mutate(missing_48h = is.na(hb_48h)) |>
  group_by(missing_48h) |>
  summarise(
    n                = n(),
    median_icu_hours = round(median(icu_hours, na.rm = TRUE), 1),
    min_icu_hours    = round(min(icu_hours, na.rm = TRUE), 1),
    max_icu_hours    = round(max(icu_hours, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# Missing data summary ----------------------------------------------------
cat("\nMissing data (%):\n")
sepsis_hb_48h |>
  summarise(across(everything(), ~ round(mean(is.na(.)) * 100, 1))) |>
  tidyr::pivot_longer(
    everything(),
    names_to  = "variable",
    values_to = "pct_missing"
  ) |>
  filter(pct_missing > 0) |>
  arrange(desc(pct_missing)) |>
  print()

# Time series summary -----------------------------------------------------
cat("\nHb by 6-hour intervals:\n")
sepsis_hb_long |>
  mutate(hour_bin = floor(hours_since_admission / 6) * 6) |>
  group_by(hour_bin) |>
  summarise(
    n         = n(),
    mean_hb   = round(mean(lab_result, na.rm = TRUE), 1),
    sd_hb     = round(sd(lab_result, na.rm = TRUE), 1),
    median_hb = round(median(lab_result, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# Normality check ---------------------------------------------------------
cat("\nNormality check:\n")
print(shapiro.test(na.omit(sepsis_hb_48h$hb_change)))
cat("Skewness:", round(skewness(sepsis_hb_48h$hb_change, na.rm = TRUE), 3), "\n")
cat("Kurtosis:", round(kurtosis(sepsis_hb_48h$hb_change, na.rm = TRUE), 3), "\n")

# =========================================================================
# Visualizations RQ1
# =========================================================================

# Plot 1 - Distribution of Hb change Sepsis cohort
p1 <- ggplot(sepsis_hb_48h, aes(x = hb_change)) +
  geom_histogram(binwidth = 2, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of Hb change - Sepsis cohort",
    x = "Hb change (g/L)",
    y = "Number of episodes"
  ) +
  theme_minimal()

# Plot 2 - QQ plot
p2 <- ggplot(sepsis_hb_48h, aes(sample = hb_change)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = "QQ-plot Hb change - Sepsis cohort",
    x = "Theoretical quantiles",
    y = "Observed quantiles"
  ) +
  theme_minimal()

# Plot 3 - Boxplot baseline vs last measurement
p3 <- sepsis_hb_48h |>
  tidyr::pivot_longer(
    cols      = c(hb_baseline, hb_last),
    names_to  = "timepoint",
    values_to = "hb"
  ) |>
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("hb_baseline", "hb_last"),
      labels = c("Baseline", "Last measurement")
    )
  ) |>
  ggplot(aes(x = timepoint, y = hb, fill = timepoint)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(
    title = "Hb at baseline vs last measurement - Sepsis cohort",
    x = "Timepoint",
    y = "Hemoglobin (g/L)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 4 - Hb trajectory over time (episode-level)
p4 <- sepsis_hb_long |>
  filter(lab_result >= 40, lab_result <= 200) |>  # remove extreme Hb values
  ggplot(aes(
    x = hours_since_admission,
    y = lab_result,
    group = episode_id
  )) +
  geom_line(alpha = 0.2, color = "steelblue") +
  geom_smooth(aes(group = 1), color = "red", se = TRUE) +
  labs(
    title = "Hb trajectory first 48h - Sepsis cohort",
    x = "Hours since ICU admission",
    y = "Hemoglobin (g/L)"
  ) +
  theme_minimal()

# Plot 5 - Hb change by patient outcome
p5 <- ggplot(
  sepsis_hb_48h,
  aes(x = patient_outcome, y = hb_change, fill = patient_outcome)
) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_fill_manual(values = c("alive" = "steelblue", "dead" = "tomato")) +
  labs(
    title = "Hb change by patient outcome - Sepsis cohort",
    x = "Outcome",
    y = "Hb change (g/L)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 6 - Hb change by measurement type
p6 <- sepsis_hb_48h |>
  filter(!is.na(hb_change)) |>
  ggplot(aes(x = measurement_type, y = hb_change, fill = measurement_type)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(
    title = "Hb change by measurement type",
    x = "",
    y = "Hb change (g/L)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

# Plot 7 - Hb trajectory by 6h bins
p7 <- sepsis_hb_long |>
  filter(hours_since_admission < 48) |>  # exclude the sparse hour 48 bin
  mutate(hour_bin = floor(hours_since_admission / 6) * 6) |>
  group_by(hour_bin) |>
  summarise(
    mean_hb = mean(lab_result, na.rm = TRUE),
    sd_hb   = sd(lab_result, na.rm = TRUE),
    n       = n(),
    .groups = "drop"
  ) |>
  mutate(
    se      = sd_hb / sqrt(n),
    ci_low  = mean_hb - 1.96 * se,
    ci_high = mean_hb + 1.96 * se
  ) |>
  ggplot(aes(x = hour_bin, y = mean_hb)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high),
              alpha = 0.2, fill = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  labs(
    title    = "Mean Hb trajectory by 6h intervals - Sepsis cohort",
    subtitle = "Mean ± 95% CI",
    x        = "Hours since ICU admission",
    y        = "Hemoglobin (g/L)"
  ) +
  theme_minimal()

# Plot 8 - Age distribution
p8 <- ggplot(sepsis_hb_48h, aes(x = age, fill = sex)) +
  geom_histogram(binwidth = 5, position = "dodge", color = "white") +
  scale_fill_manual(values = c("male" = "steelblue", "female" = "tomato")) +
  labs(
    title = "Age distribution by sex - Sepsis cohort",
    x = "Age (years)",
    y = "Count",
    fill = "Sex"
  ) +
  theme_minimal()

cat("\nDisplaying RQ1 plots...\n")
print((p1 + p2) / (p3 + p4) / (p5 + p6))
print(p7)
print(p8)
library(patchwork)

# Save all RQ1 figures — no titles or subtitles
ggsave("figures/fig1_distribution.png",
       plot = p1 + p2,
       width = 12, height = 5, dpi = 300, bg = "white")

ggsave("figures/fig2_baseline_trajectory.png",
       plot = p3 + p4,
       width = 12, height = 5, dpi = 300, bg = "white")

ggsave("figures/fig3_subgroups.png",
       plot = p5 + p6,
       width = 12, height = 5, dpi = 300, bg = "white")

ggsave("figures/fig4_mean_trajectory.png",
       plot = p7,
       width = 8, height = 5, dpi = 300, bg = "white")

ggsave("figures/fig5_age_sex.png",
       plot = p8,
       width = 8, height = 5, dpi = 300, bg = "white")

cat("All RQ1 figures saved to figures/\n")

# Combined RQ1 figure — all plots in one PNG
rq1_combined <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8)

ggsave("figures/fig_rq1_combined.png",
       plot   = rq1_combined,
       width  = 14,
       height = 20,
       dpi    = 300,
       bg     = "white")

cat("Saved: figures/fig_rq1_combined.png\n")

# =========================================================================
# RESEARCH QUESTION 2 - Comparison sepsis vs neurotrauma/SAH
# =========================================================================
cat("\n=== RQ2: Comparison across diagnosis groups ===\n\n")

cat("Comparison dataset size:\n")
comparison_hb_48h |>
  count(diagnosis_group) |>
  print()

# =========================================================================
# Visualizations RQ2
# =========================================================================

# Plot 9 - Hb change by diagnosis group
p9 <- comparison_hb_48h |>
  filter(!is.na(hb_change)) |>
  ggplot(
    aes(x = diagnosis_group, y = hb_change, fill = diagnosis_group)
  ) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2) +
  scale_fill_manual(values = c(
    "Sepsis"      = "steelblue",
    "Neurotrauma" = "darkorange",
    "SAH"         = "forestgreen"
  )) +
  labs(
    title = "Hb change by diagnosis group",
    x = "Diagnosis group",
    y = "Hb change (g/L)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 10 - Baseline Hb by group
p10 <- comparison_hb_48h |>
  filter(!is.na(hb_baseline)) |>
  ggplot(
    aes(x = diagnosis_group, y = hb_baseline, fill = diagnosis_group)
  ) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2) +
  scale_fill_manual(values = c(
    "Sepsis"      = "steelblue",
    "Neurotrauma" = "darkorange",
    "SAH"         = "forestgreen"
  )) +
  labs(
    title = "Baseline Hb by diagnosis group",
    x = "Diagnosis group",
    y = "Baseline Hb (g/L)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 11 - Sampling volume by group
p11 <- comparison_hb_48h |>
  filter(!is.na(total_volume_ml)) |>
  ggplot(aes(x = diagnosis_group, y = total_volume_ml, fill = diagnosis_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2) +
  scale_fill_manual(values = c(
    "Sepsis"      = "steelblue",
    "Neurotrauma" = "darkorange",
    "SAH"         = "forestgreen"
  )) +
  labs(
    title = "Blood sampling volume by diagnosis group",
    x = "Diagnosis group",
    y = "Total sampling volume (mL)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

cat("\nDisplaying RQ2 plots...\n")
print(p9 / (p10 + p11))

# ===Export RQ2 plots as PNG

# Combined RQ2 figure
ggsave("figures/fig_rq2_combined.png",
       plot   = p9 / (p10 + p11),
       width  = 10,
       height = 12,
       dpi    = 300,
       bg     = "white")
cat("Saved: figures/fig_rq2_combined.png\n")

# Individual plots if needed separately
ggsave("figures/fig_rq2_hb_change.png",
       plot   = p9,
       width  = 8,
       height = 6,
       dpi    = 300,
       bg     = "white")
cat("Saved: figures/fig_rq2_hb_change.png\n")

ggsave("figures/fig_rq2_baseline_hb.png",
       plot   = p10,
       width  = 8,
       height = 6,
       dpi    = 300,
       bg     = "white")
cat("Saved: figures/fig_rq2_baseline_hb.png\n")

ggsave("figures/fig_rq2_sampling_volume.png",
       plot   = p11,
       width  = 8,
       height = 6,
       dpi    = 300,
       bg     = "white")
cat("Saved: figures/fig_rq2_sampling_volume.png\n")

# =========================================================================
# COMBINED DESCRIPTIVE TABLE 1. Patient and episode characterisics by diagnosis group
# =========================================================================

# Prepare table dataset ---------------------------------------------------
table_data <- comparison_hb_48h |>
  mutate(
    diagnosis_group = factor(diagnosis_group,
                             levels = c("Sepsis", "Neurotrauma", "SAH")),
    patient_outcome = factor(patient_outcome,
                             levels = c("alive", "dead"),
                             labels = c("Alive", "Dead")),
    sex = factor(sex,
                 levels = c("male", "female"),
                 labels = c("Male", "Female"))
  )

# Vector-safe version of format_p for gtsummary
format_p_vec <- function(x) {
  sapply(x, function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("<0.001")
    sprintf("%.3f", p)
  })
}

table_combined <- table_data |>
  dplyr::select(
    diagnosis_group,
    age, sex,
    icu_hours, patient_outcome,
    hb_baseline, hb_change, hb_change_pct,
    crp, crp_24h, lpk, tpk, kreatinin, bilirubin,
    fluid_balance_48h, rbc_48h, plasma_48h, total_volume_ml, total_draws
  ) |>
  tbl_summary(
    by = diagnosis_group,
    label = list(
      age               ~ "Age (years)",
      sex               ~ "Sex",
      icu_hours         ~ "ICU stay (hours)",
      patient_outcome   ~ "Patient outcome",
      hb_baseline       ~ "Baseline Hb (g/L)",
      hb_change         ~ "Hb change (g/L)",
      hb_change_pct     ~ "Hb change (%)",
      crp               ~ "CRP admission (mg/L)",
      crp_24h           ~ "CRP 24-48h (mg/L)",
      lpk               ~ "Leukocytes (10⁹/L)",
      tpk               ~ "Platelets (10⁹/L)",
      kreatinin         ~ "Creatinine (µmol/L)",
      bilirubin         ~ "Bilirubin (µmol/L)",
      fluid_balance_48h ~ "Fluid balance 48h (mL)",
      rbc_48h           ~ "RBC transfusion (mL)",
      plasma_48h        ~ "Plasma transfusion (mL)",
      total_volume_ml   ~ "Blood sampling volume (mL)",
      total_draws       ~ "Total blood draws (n)"
    ),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing      = "ifany",
    missing_text = "n = Missing value"
  ) |>
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = format_p_vec
  ) |>
  add_overall(last = FALSE) |>
  bold_labels() |>
  modify_header(
    label  ~ "**Variable**",
    stat_0 ~ "**Overall**\nn = {N}",
    stat_1 ~ "**Sepsis**\nn = {n}",
    stat_2 ~ "**Neurotrauma**\nn = {n}",
    stat_3 ~ "**SAH**\nn = {n}"
  )

print(table_combined)

# Save with caption added via gt
table_combined |>
  as_gt() |>
  gt::tab_header(
    title   = "Table 1.Patient and episode characteristics by diagnosis group"
  ) |>
  gt::gtsave("figures/table1_combined.docx")

cat("Table saved to figures/table1_combined.docx\n")