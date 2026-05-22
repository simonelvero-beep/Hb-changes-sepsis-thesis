# =============== 03a_merge_sepsis_hb.R ====================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")

# Packages ----------------------------------------------------------------
packages <- c("dplyr", "stringr", "lubridate", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Load clean data ---------------------------------------------------------
civakohort_clean <- readRDS("data_clean/civakohort_clean.rds")
blodgashb_clean  <- readRDS("data_clean/blodgashb_clean.rds")

# Step 1 - Filter sepsis and septic shock episodes ------------------------
# Diagnos codes accoringd to ICD10-SE
sepsis_episodes <- civakohort_clean |>
  filter(
    if_any(
      c(ivadiagnos_1, ivadiagnos_2, ivadiagnos_3,
        diagnos_1, diagnos_2, diagnos_3),
      ~ str_detect(., "^A40|^A41|^R572|^R651")
    )
  )

# Remove readmissions within 30 days, but KEEP episode_id
sepsis_patients <- sepsis_episodes |>
  group_by(personal_id) |>
  arrange(icu_admission, .by_group = TRUE) |>
  mutate(
    days_since_last = as.numeric(difftime(
      icu_admission,
      lag(icu_admission),
      units = "days"
    ))
  ) |>
  filter(is.na(days_since_last) | days_since_last >= 30) |>
  ungroup() |>
  dplyr::select(
    episode_id, personal_id, icu_admission, icu_discharge, patient_outcome
  )

cat("Sepsis patients:\n")
cat("  Episodes kept:   ", nrow(sepsis_patients), "\n")
cat("  Unique patients: ", n_distinct(sepsis_patients$personal_id), "\n")

# Step 2 - Link Hb measurements to the correct ICU episode ----------------
# blodgashb_clean does not necessarily contain episode_id, so we assign it
# by matching on personal_id and requiring analysis_time to fall within
# the ICU admission-discharge interval.

sepsis_hb <- blodgashb_clean |>
  filter(parameter_name %in% c(
    "aB--Hemoglobin(PNA)",
    "MAN aB--Hemoglobin(PNA)"
  )) |>
  inner_join(
    sepsis_patients,
    by = "personal_id",
    relationship = "many-to-many"
  ) |>
  filter(
    !is.na(analysis_time),
    analysis_time >= icu_admission,
    analysis_time <= icu_discharge
  ) |>
  mutate(
    hours_since_admission = as.numeric(difftime(
      analysis_time, icu_admission, units = "hours"
    ))
  ) |>
  arrange(episode_id, analysis_time)

cat("  Episodes with Hb data: ", n_distinct(sepsis_hb$episode_id), "\n")
cat("  Hb measurements on ICU:", nrow(sepsis_hb), "\n")

# Optional safety check: one Hb row should belong to max one episode
dup_matches <- sepsis_hb |>
  count(personal_id, analysis_time, lab_result) |>
  filter(n > 1)

cat("  Potential duplicate episode matches:", nrow(dup_matches), "\n")

# Step 3 - Calculate baseline, last Hb <=60h, and Hb nearest 48h ----------
sepsis_hb_48h <- sepsis_hb |>
  group_by(episode_id, personal_id, icu_admission, icu_discharge, patient_outcome) |>
  arrange(analysis_time, .by_group = TRUE) |>
  summarise(
    hb_baseline = {
      idx_0_60 <- which(hours_since_admission >= 0 & hours_since_admission <= 60)
      if (length(idx_0_60) > 0) lab_result[idx_0_60[1]] else NA_real_
    },
    
    baseline_time = {
      idx_0_60 <- which(hours_since_admission >= 0 & hours_since_admission <= 60)
      if (length(idx_0_60) > 0) analysis_time[idx_0_60[1]] else as.POSIXct(NA)
    },
    
    hb_last = {
      idx_0_60 <- which(hours_since_admission >= 0 & hours_since_admission <= 60)
      if (length(idx_0_60) > 0) lab_result[idx_0_60[length(idx_0_60)]] else NA_real_
    },
    
    last_measurement_time = {
      idx_0_60 <- which(hours_since_admission >= 0 & hours_since_admission <= 60)
      if (length(idx_0_60) > 0) analysis_time[idx_0_60[length(idx_0_60)]] else as.POSIXct(NA)
    },
    
    hours_at_last = {
      idx_0_60 <- which(hours_since_admission >= 0 & hours_since_admission <= 60)
      if (length(idx_0_60) > 0) hours_since_admission[idx_0_60[length(idx_0_60)]] else NA_real_
    },
    
    hb_48h = {
      idx_48 <- which(hours_since_admission >= 36 & hours_since_admission <= 60)
      if (length(idx_48) > 0) {
        idx_best <- idx_48[which.min(abs(hours_since_admission[idx_48] - 48))]
        lab_result[idx_best]
      } else {
        NA_real_
      }
    },
    
    hb_48h_time = {
      idx_48 <- which(hours_since_admission >= 36 & hours_since_admission <= 60)
      if (length(idx_48) > 0) {
        idx_best <- idx_48[which.min(abs(hours_since_admission[idx_48] - 48))]
        analysis_time[idx_best]
      } else {
        as.POSIXct(NA)
      }
    },
    
    n_measurements_0_60h = sum(hours_since_admission >= 0 & hours_since_admission <= 60, na.rm = TRUE),
    n_measurements_total = n(),
    .groups = "drop"
  ) |>
  mutate(
    hb_change     = hb_last - hb_baseline,
    hb_change_pct = ((hb_last - hb_baseline) / hb_baseline) * 100,
    icu_hours     = as.numeric(difftime(icu_discharge, icu_admission, units = "hours")),
    measurement_type = case_when(
      !is.na(hb_48h) ~ "36-60h window",
      icu_hours < 36 & !is.na(hb_last) ~ "discharged before 36h",
      icu_hours >= 36 & !is.na(hb_last) ~ "no 36-60h Hb",
      TRUE ~ "no Hb within 60h"
    )
  )

# Step 4 - Add age and sex ------------------------------------------------
sepsis_hb_48h <- sepsis_hb_48h |>
  mutate(
    birthdate = as.Date(str_sub(personal_id, 1, 8), format = "%Y%m%d"),
    age       = round(as.numeric(difftime(icu_admission, birthdate, units = "days")) / 365.25, 1),
    sex_digit = as.numeric(str_sub(personal_id, 11, 11)),
    sex       = factor(ifelse(sex_digit %% 2 == 1, "male", "female"))
  ) |>
  dplyr::select(-birthdate, -sex_digit)

# Step 5 - Create longitudinal dataset at EPISODE level -------------------
sepsis_hb_long <- sepsis_hb |>
  filter(hours_since_admission >= 0, hours_since_admission <= 48) |>
  mutate(
    birthdate = as.Date(str_sub(personal_id, 1, 8), format = "%Y%m%d"),
    age       = round(as.numeric(difftime(icu_admission, birthdate, units = "days")) / 365.25, 1),
    sex_digit = as.numeric(str_sub(personal_id, 11, 11)),
    sex       = factor(ifelse(sex_digit %% 2 == 1, "male", "female"))
  ) |>
  dplyr::select(-birthdate, -sex_digit)

# Verify ------------------------------------------------------------------
cat("\nsepsis_hb_48h measurement types:\n")
print(sepsis_hb_48h |> count(measurement_type))

cat("\nAge and sex distribution:\n")
print(sepsis_hb_48h |> count(sex))
cat("Mean age:", round(mean(sepsis_hb_48h$age, na.rm = TRUE), 1), "\n")

cat("\nEpisodes in sepsis_hb_48h:", n_distinct(sepsis_hb_48h$episode_id), "\n")
cat("Episodes in sepsis_hb_long:", n_distinct(sepsis_hb_long$episode_id), "\n")

# Save --------------------------------------------------------------------
saveRDS(sepsis_patients, "data_clean/sepsis_patients.rds")
saveRDS(sepsis_hb,       "data_clean/sepsis_hb.rds")
saveRDS(sepsis_hb_48h,   "data_clean/sepsis_hb_48h.rds")
saveRDS(sepsis_hb_long,  "data_clean/sepsis_hb_long.rds")