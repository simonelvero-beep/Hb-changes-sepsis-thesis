# ============== 03b_merge_comparison.R =====================
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

# Definitions, needs to be updated in 03c aswell!!!
ml_per_draw_bloodgas <- 1.5
ml_per_draw_labkemi <- 10

diag_cols <- c(
  "ivadiagnos_1", "ivadiagnos_2", "ivadiagnos_3",
  "diagnos_1", "diagnos_2", "diagnos_3"
)

has_code <- function(x, pattern) {
  str_detect(
    coalesce(as.character(x), ""),
    regex(pattern, ignore_case = TRUE)
  )
}
# Load clean data ---------------------------------------------------------
civakohort_clean <- readRDS("data_clean/civakohort_clean.rds")
nivakohort_clean <- readRDS("data_clean/nivakohort_clean.rds")
blodgashb_clean  <- readRDS("data_clean/blodgashb_clean.rds")
labkemi_sampling <- readRDS("data_clean/labkemi_sampling.rds")
sepsis_hb_48h    <- readRDS("data_clean/sepsis_hb_48h.rds")

# Step 1 - Filter neurotrauma and SAH from civakohort ---------------------
neuro_sah_civa <- civakohort_clean |>
  filter(
    if_any(
      all_of(diag_cols),
      ~ has_code(.x, "^(S06|I60)")
    )
  ) |>
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
  mutate(cohort = "civa") |>
  mutate(
    diagnosis_group = case_when(
      if_any(all_of(diag_cols), ~ has_code(.x, "^I60")) ~ "SAH",
      if_any(all_of(diag_cols), ~ has_code(.x, "^S06")) ~ "Neurotrauma",
      TRUE ~ "Other"
    )
  ) |>
  dplyr::select(
    episode_id, cohort, personal_id, icu_admission, icu_discharge,
    patient_outcome, diagnosis_group,
    ivadiagnos_1, ivadiagnos_2, ivadiagnos_3,
    diagnos_1, diagnos_2, diagnos_3
  )

# Step 2 - Filter neurotrauma and SAH from nivakohort ---------------------
neuro_sah_niva <- nivakohort_clean |>
  filter(
    if_any(
      all_of(diag_cols),
      ~ has_code(.x, "^(S06|I60)")
    )
  ) |>
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
  mutate(cohort = "niva") |>
  mutate(
    diagnosis_group = case_when(
      if_any(all_of(diag_cols), ~ has_code(.x, "^I60")) ~ "SAH",
      if_any(all_of(diag_cols), ~ has_code(.x, "^S06")) ~ "Neurotrauma",
      TRUE ~ "Other"
    )
  ) |>
  dplyr::select(
    episode_id, cohort, personal_id, icu_admission, icu_discharge,
    patient_outcome, diagnosis_group,
    ivadiagnos_1, ivadiagnos_2, ivadiagnos_3,
    diagnos_1, diagnos_2, diagnos_3
  )

# Step 3 - Combine and keep only classified episodes ----------------------
neuro_sah_patients <- bind_rows(neuro_sah_civa, neuro_sah_niva) |>
  filter(diagnosis_group %in% c("SAH", "Neurotrauma"))|>

# exclude any episodes with sepsis codes
  filter(
    !if_any(
      all_of(diag_cols),
      ~ has_code(.x, "^(A40|A41|R572|R651)")
    )
  )

cat("Neuro/SAH patients by group and cohort:\n")
print(neuro_sah_patients |> count(diagnosis_group, cohort))

# Step 4 - Link Hb measurements to the correct ICU episode ----------------
# blodgashb_clean does not have reliable episode_id, so map by
# personal_id + time within ICU stay
neuro_sah_hb <- blodgashb_clean |>
  filter(parameter_name %in% c(
    "aB--Hemoglobin(PNA)",
    "MAN aB--Hemoglobin(PNA)"
  )) |>
  inner_join(
    neuro_sah_patients |>
      dplyr::select(
        episode_id, cohort, personal_id,
        icu_admission, icu_discharge,
        patient_outcome, diagnosis_group
      ),
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
  arrange(cohort, episode_id, analysis_time)

cat("\nEpisodes with Hb data:",
    n_distinct(paste(neuro_sah_hb$cohort, neuro_sah_hb$episode_id)), "\n")
cat("Hb measurements on ICU:", nrow(neuro_sah_hb), "\n")

dup_matches <- neuro_sah_hb |>
  count(personal_id, analysis_time, lab_result) |>
  filter(n > 1)

cat("Potential duplicate episode matches:", nrow(dup_matches), "\n")

# Step 5 - Calculate baseline and 48h Hb at episode level -----------------
neuro_sah_hb_48h <- neuro_sah_hb |>
  group_by(
    cohort, episode_id, personal_id,
    icu_admission, icu_discharge,
    patient_outcome, diagnosis_group
  ) |>
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
  ) |>
  mutate(
    birthdate = as.Date(str_sub(personal_id, 1, 8), format = "%Y%m%d"),
    age       = round(as.numeric(difftime(icu_admission, birthdate, units = "days")) / 365.25, 1),
    sex_digit = as.numeric(str_sub(personal_id, 11, 11)),
    sex       = factor(ifelse(sex_digit %% 2 == 1, "male", "female"))
  ) |>
  dplyr::select(-birthdate, -sex_digit)

# Step 6 - Create combined comparison dataset -----------------------------
comparison_hb_48h <- bind_rows(
  sepsis_hb_48h |>
    mutate(
      diagnosis_group = "Sepsis",
      cohort = "civa"
    ),
  neuro_sah_hb_48h
)

cat("\nComparison dataset by diagnosis group:\n")
print(comparison_hb_48h |> count(diagnosis_group))

# Step 7 - Calculate blood sampling for neuro/SAH at episode level --------

# Blood gas: still needs time-based mapping from personal_id to episode
neuro_sah_hb_sampling <- blodgashb_clean |>
  inner_join(
    neuro_sah_patients |>
      dplyr::select(episode_id, cohort, personal_id, icu_admission, icu_discharge),
    by = "personal_id",
    relationship = "many-to-many"
  ) |>
  filter(
    !is.na(analysis_time),
    analysis_time >= icu_admission,
    analysis_time <= icu_discharge
  ) |>
  mutate(
    hours_from_admission = as.numeric(difftime(
      analysis_time, icu_admission, units = "hours"
    ))
  ) |>
  filter(hours_from_admission >= 0, hours_from_admission <= 48) |>
  group_by(cohort, episode_id, analysis_time) |>
  slice(1) |>
  ungroup() |>
  group_by(cohort, episode_id) |>
  summarise(
    n_bloodgas_draws = n(),
    vol_bloodgas_ml  = n() * ml_per_draw_bloodgas,
    .groups = "drop"
  )

# Lab chemistry: labkemi_sampling already has correct episode_id
# Join by episode_id + personal_id only to add cohort/admission metadata
neuro_sah_labkemi_sampling <- labkemi_sampling |>
  inner_join(
    neuro_sah_patients |>
      dplyr::select(episode_id, personal_id, cohort, icu_admission, icu_discharge),
    by = c("episode_id", "personal_id")
  ) |>
  filter(
    !is.na(sample_date),
    sample_date >= icu_admission,
    sample_date <= icu_discharge
  ) |>
  mutate(
    hours_from_admission = as.numeric(difftime(
      sample_date, icu_admission, units = "hours"
    ))
  ) |>
  filter(hours_from_admission >= 0, hours_from_admission <= 48) |>
  arrange(cohort, episode_id, sample_date) |>
  group_by(cohort, episode_id) |>
  mutate(
    time_diff_min = as.numeric(difftime(
      sample_date, lag(sample_date), units = "mins"
    )),
    new_draw = is.na(time_diff_min) | time_diff_min > 0,
    draw_id  = cumsum(new_draw)
  ) |>
  summarise(
    n_labkemi_draws = n_distinct(draw_id),
    vol_labkemi_ml  = n_distinct(draw_id) * ml_per_draw_labkemi,
    .groups = "drop"
  )

# Combine blood gas and lab chemistry sampling for neuro/SAH
neuro_sah_sampling <- neuro_sah_hb_sampling |>
  full_join(neuro_sah_labkemi_sampling, by = c("cohort", "episode_id")) |>
  mutate(
    n_bloodgas_draws = replace_na(n_bloodgas_draws, 0),
    vol_bloodgas_ml  = replace_na(vol_bloodgas_ml, 0),
    n_labkemi_draws  = replace_na(n_labkemi_draws, 0),
    vol_labkemi_ml   = replace_na(vol_labkemi_ml, 0),
    total_draws      = n_bloodgas_draws + n_labkemi_draws,
    total_volume_ml  = vol_bloodgas_ml + vol_labkemi_ml
  ) |>
  dplyr::select(
    cohort, episode_id, total_volume_ml, total_draws,
    n_bloodgas_draws, n_labkemi_draws,
    vol_bloodgas_ml, vol_labkemi_ml
  )

# Get sepsis sampling - may be absent until 03c has run
sepsis_sampling <- sepsis_hb_48h |>
  mutate(cohort = "civa") |>
  dplyr::select(any_of(c(
    "cohort", "episode_id", "total_volume_ml", "total_draws",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  )))

# Combine sepsis and neuro/SAH sampling into one dataset
all_sampling <- bind_rows(sepsis_sampling, neuro_sah_sampling)

# Step 8 - Join sampling with comparison dataset --------------------------
comparison_hb_48h <- comparison_hb_48h |>
  dplyr::select(-any_of(c(
    "total_volume_ml", "total_draws",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  ))) |>
  left_join(all_sampling, by = c("cohort", "episode_id"))

# Verify sampling coverage by group
cat("\nSampling coverage by diagnosis group:\n")
comparison_hb_48h |>
  group_by(diagnosis_group) |>
  summarise(
    n               = n(),
    n_with_sampling = sum(!is.na(total_volume_ml)),
    median_volume   = round(median(total_volume_ml, na.rm = TRUE), 1),
    iqr_volume      = round(IQR(total_volume_ml, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# Save --------------------------------------------------------------------
saveRDS(neuro_sah_patients, "data_clean/neuro_sah_patients.rds")
saveRDS(neuro_sah_hb,       "data_clean/neuro_sah_hb.rds")
saveRDS(neuro_sah_hb_48h,   "data_clean/neuro_sah_hb_48h.rds")
saveRDS(comparison_hb_48h,  "data_clean/comparison_hb_48h.rds")