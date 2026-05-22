# ========== 03c_merge_inflammation.R =============================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")
# MUST run after 03a and 03b

# Packages ----------------------------------------------------------------
# NOTE: uses join_by() with inequality joins -> requires dplyr >= 1.1
packages <- c("dplyr", "stringr", "lubridate", "tidyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Definitions, needs to be updated in 03b aswell!!!
ml_per_draw_bloodgas <- 1.5
ml_per_draw_labkemi <- 10

# Load clean data ---------------------------------------------------------
sepsis_patients    <- readRDS("data_clean/sepsis_patients.rds")     # episode-level
sepsis_hb_48h      <- readRDS("data_clean/sepsis_hb_48h.rds")       # episode-level
sepsis_hb_long     <- readRDS("data_clean/sepsis_hb_long.rds")      # episode-level
labkemi_clean      <- readRDS("data_clean/labkemi_clean.rds")       # already has correct episode_id
labkemi_sampling   <- readRDS("data_clean/labkemi_sampling.rds")    # already has correct episode_id
vatskebalans_clean <- readRDS("data_clean/vatskebalans_clean.rds")
blodgashb_clean    <- readRDS("data_clean/blodgashb_clean.rds")

# Helper episode windows --------------------------------------------------
episodes_48h <- sepsis_patients |>
  mutate(
    followup_end = pmin(icu_discharge, icu_admission + hours(48))
  )

# =========================================================================
# STEP 1 - Link lab chemistry to sepsis episodes
# =========================================================================
# IMPORTANT:
# labkemi_clean already has correct episode_id -> use it directly
labkemi_sepsis <- labkemi_clean |>
  inner_join(
    sepsis_patients |>
      dplyr::select(episode_id, icu_admission, icu_discharge),
    by = "episode_id"
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
  )

# =========================================================================
# STEP 2 - Extract lab values at admission (0-24h) by EPISODE
# =========================================================================
lab_params <- c(
  "P-CRP", "B-Trombocyter", "B-Leukocyter",
  "P-Kreatinin", "P-Bilirubin", "(B)Erc-MCH", "(B)Erc-MCV"
)

sepsis_labs_admission <- labkemi_sepsis |>
  filter(
    parameter_name %in% lab_params,
    hours_from_admission >= 0,
    hours_from_admission <= 24
  ) |>
  group_by(episode_id, parameter_name) |>
  arrange(hours_from_admission, .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  dplyr::select(episode_id, parameter_name, lab_result) |>
  tidyr::pivot_wider(names_from = parameter_name, values_from = lab_result) |>
  rename(
    crp       = "P-CRP",
    tpk       = "B-Trombocyter",
    lpk       = "B-Leukocyter",
    kreatinin = "P-Kreatinin",
    bilirubin = "P-Bilirubin",
    mch       = "(B)Erc-MCH",
    mcv       = "(B)Erc-MCV"
  )

cat("Lab values at admission:\n")
cat("  Episodes with any lab:", nrow(sepsis_labs_admission), "\n")

# =========================================================================
# STEP 3 - Extract CRP at 24-48h by EPISODE
# =========================================================================
crp_24_48h <- labkemi_sepsis |>
  filter(
    parameter_name == "P-CRP",
    hours_from_admission >= 24,
    hours_from_admission <= 48
  ) |>
  group_by(episode_id) |>
  arrange(hours_from_admission, .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  dplyr::select(episode_id, crp_24h = lab_result)

cat("CRP at 24-48h:\n")
cat("  Episodes:", nrow(crp_24_48h), "\n")
cat("  Median:", round(median(crp_24_48h$crp_24h, na.rm = TRUE), 1), "\n")

# =========================================================================
# STEP 4 - Calculate fluid balance for first 48h by EPISODE
# =========================================================================
# IMPORTANT:
# Fluid balance does not have episode_id, so it still needs time-based
# mapping from personal_id to episode.
# Since balance is available as 24h periods, overlap with first 48h is an approximation.

sepsis_fluid_48h <- vatskebalans_clean |>
  inner_join(
    episodes_48h |>
      dplyr::select(episode_id, personal_id, icu_admission, icu_discharge, followup_end),
    by = "personal_id",
    relationship = "many-to-many"
  ) |>
  mutate(
    period_start = day_break_date,
    period_end   = day_break_date + hours(24)
  ) |>
  filter(
    !is.na(period_start),
    period_end > icu_admission,
    period_start < followup_end
  ) |>
  group_by(episode_id) |>
  summarise(
    fluid_balance_48h   = sum(fluid_balance_total,   na.rm = TRUE),
    colloid_balance_48h = sum(colloid_balance_total, na.rm = TRUE),
    total_balance_48h   = sum(total_balance,         na.rm = TRUE),
    rbc_48h             = sum(red_blood_cells,       na.rm = TRUE),
    colloid_48h         = sum(colloid,               na.rm = TRUE),
    platelets_48h       = sum(platelets,             na.rm = TRUE),
    plasma_48h          = sum(plasma,                na.rm = TRUE),
    n_days              = n(),
    .groups = "drop"
  )

cat("Fluid balance coverage:\n")
cat("  Episodes:", nrow(sepsis_fluid_48h), "\n")
cat(
  "  Median fluid balance:",
  round(median(sepsis_fluid_48h$fluid_balance_48h, na.rm = TRUE), 1),
  "mL\n"
)

# =========================================================================
# STEP 5 - Join static variables to sepsis_hb_48h and sepsis_hb_long
# =========================================================================
sepsis_hb_48h <- sepsis_hb_48h |>
  dplyr::select(-any_of(c(
    "crp", "lpk", "tpk", "kreatinin", "bilirubin",
    "mch", "mcv", "crp_24h",
    "fluid_balance_48h", "colloid_balance_48h", "total_balance_48h",
    "rbc_48h", "colloid_48h", "platelets_48h", "plasma_48h", "n_days"
  ))) |>
  left_join(sepsis_labs_admission, by = "episode_id") |>
  left_join(crp_24_48h,           by = "episode_id") |>
  left_join(sepsis_fluid_48h,     by = "episode_id")

sepsis_hb_long <- sepsis_hb_long |>
  dplyr::select(-any_of(c(
    "crp", "lpk", "tpk", "kreatinin", "bilirubin",
    "mch", "mcv", "crp_24h",
    "fluid_balance_48h", "colloid_balance_48h", "total_balance_48h",
    "rbc_48h", "colloid_48h", "platelets_48h", "plasma_48h", "n_days",
    "hb_baseline"
  ))) |>
  left_join(
    sepsis_hb_48h |>
      dplyr::select(
        episode_id, crp, lpk, tpk, mch, mcv,
        kreatinin, bilirubin, hb_baseline,
        crp_24h, fluid_balance_48h,
        colloid_balance_48h, total_balance_48h,
        rbc_48h, colloid_48h, platelets_48h, plasma_48h
      ),
    by = "episode_id"
  )

cat("\nLab coverage in sepsis_hb_48h:\n")
sepsis_hb_48h |>
  summarise(
    n_crp     = sum(!is.na(crp)),
    n_crp_24h = sum(!is.na(crp_24h)),
    n_lpk     = sum(!is.na(lpk)),
    n_fluid   = sum(!is.na(fluid_balance_48h)),
    n_rbc     = sum(!is.na(rbc_48h))
  ) |>
  print()

# =========================================================================
# STEP 6 - Calculate 48h blood sampling volume by EPISODE
# =========================================================================

# Blood gas draws: each unique timestamp = one draw
# blodgashb_clean does not have reliable episode_id -> still map by personal_id + time
hb_sampling <- blodgashb_clean |>
  inner_join(
    episodes_48h |>
      dplyr::select(episode_id, personal_id, icu_admission, followup_end),
    by = "personal_id",
    relationship = "many-to-many"
  ) |>
  filter(
    !is.na(analysis_time),
    analysis_time >= icu_admission,
    analysis_time <= followup_end
  ) |>
  group_by(episode_id, analysis_time) |>
  slice(1) |>
  ungroup() |>
  group_by(episode_id) |>
  summarise(
    n_bloodgas_draws = n(),
    vol_bloodgas_ml  = n() * ml_per_draw_bloodgas,
    .groups = "drop"
  )

# Lab chemistry draws: labkemi_sampling already has correct episode_id
labkemi_sampling_grouped <- labkemi_sampling |>
  inner_join(
    episodes_48h |>
      dplyr::select(episode_id, icu_admission, followup_end),
    by = "episode_id"
  ) |>
  filter(
    !is.na(sample_date),
    sample_date >= icu_admission,
    sample_date <= followup_end
  ) |>
  arrange(episode_id, sample_date) |>
  group_by(episode_id) |>
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

blood_sampling_48h <- hb_sampling |>
  full_join(labkemi_sampling_grouped, by = "episode_id") |>
  mutate(
    n_bloodgas_draws = replace_na(n_bloodgas_draws, 0),
    vol_bloodgas_ml  = replace_na(vol_bloodgas_ml, 0),
    n_labkemi_draws  = replace_na(n_labkemi_draws, 0),
    vol_labkemi_ml   = replace_na(vol_labkemi_ml, 0),
    total_draws      = n_bloodgas_draws + n_labkemi_draws,
    total_volume_ml  = vol_bloodgas_ml + vol_labkemi_ml
  )

cat("\nBlood sampling coverage:\n")
sepsis_hb_48h |>
  dplyr::select(-any_of(c(
    "total_draws", "total_volume_ml",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  ))) |>
  left_join(blood_sampling_48h, by = "episode_id") |>
  summarise(
    n_total_draws  = sum(!is.na(total_draws)),
    mean_draws     = round(mean(total_draws, na.rm = TRUE), 1),
    mean_volume_ml = round(mean(total_volume_ml, na.rm = TRUE), 1),
    max_volume_ml  = round(max(total_volume_ml, na.rm = TRUE), 1)
  ) |>
  print()

# Join blood sampling to sepsis datasets ----------------------------------
sepsis_hb_48h <- sepsis_hb_48h |>
  dplyr::select(-any_of(c(
    "total_draws", "total_volume_ml",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  ))) |>
  left_join(blood_sampling_48h, by = "episode_id")

sepsis_hb_long <- sepsis_hb_long |>
  dplyr::select(-any_of(c(
    "total_draws", "total_volume_ml",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  ))) |>
  left_join(blood_sampling_48h, by = "episode_id")

# =========================================================================
# STEP 7 - Update comparison dataset with complete sepsis sampling
# =========================================================================
comparison_hb_48h <- readRDS("data_clean/comparison_hb_48h.rds")

sepsis_sampling_updated <- sepsis_hb_48h |>
  mutate(cohort = "civa") |>
  dplyr::select(any_of(c(
    "cohort", "episode_id",
    "total_volume_ml", "total_draws",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  )))

comparison_non_sepsis_sampling <- comparison_hb_48h |>
  filter(diagnosis_group != "Sepsis") |>
  dplyr::select(any_of(c(
    "cohort", "episode_id",
    "total_volume_ml", "total_draws",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  )))

comparison_hb_48h <- comparison_hb_48h |>
  dplyr::select(-any_of(c(
    "total_volume_ml", "total_draws",
    "n_bloodgas_draws", "n_labkemi_draws",
    "vol_bloodgas_ml", "vol_labkemi_ml"
  ))) |>
  left_join(
    bind_rows(sepsis_sampling_updated, comparison_non_sepsis_sampling),
    by = c("cohort", "episode_id")
  )

cat("\nUpdated sampling coverage by diagnosis group:\n")
comparison_hb_48h |>
  group_by(diagnosis_group) |>
  summarise(
    n               = n(),
    n_with_sampling = sum(!is.na(total_volume_ml)),
    median_volume   = round(median(total_volume_ml, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

# =========================================================================
# STEP 8 - Create time-updated cumulative variables by EPISODE
# =========================================================================

# 8a. Cumulative blood gas sampling events --------------------------------
# blodgashb_clean still needs time-based mapping
blodgashb_timed <- blodgashb_clean |>
  inner_join(
    episodes_48h |>
      dplyr::select(episode_id, personal_id, icu_admission, followup_end),
    by = "personal_id",
    relationship = "many-to-many"
  ) |>
  filter(
    !is.na(analysis_time),
    analysis_time >= icu_admission,
    analysis_time <= followup_end
  ) |>
  group_by(episode_id, analysis_time) |>
  slice(1) |>
  ungroup() |>
  arrange(episode_id, analysis_time) |>
  group_by(episode_id) |>
  mutate(
    cum_bloodgas_draws = row_number(),
    cum_bloodgas_vol   = cum_bloodgas_draws * ml_per_draw_bloodgas
  ) |>
  ungroup() |>
  rename(event_time = analysis_time) |>
  dplyr::select(episode_id, event_time, cum_bloodgas_draws, cum_bloodgas_vol)

# 8b. Cumulative lab chemistry sampling events ----------------------------
# labkemi_sampling already has correct episode_id
labkemi_sampling_timed <- labkemi_sampling |>
  inner_join(
    episodes_48h |>
      dplyr::select(episode_id, icu_admission, followup_end),
    by = "episode_id"
  ) |>
  filter(
    !is.na(sample_date),
    sample_date >= icu_admission,
    sample_date <= followup_end
  ) |>
  arrange(episode_id, sample_date) |>
  group_by(episode_id) |>
  mutate(
    time_diff_min = as.numeric(difftime(
      sample_date, lag(sample_date), units = "mins"
    )),
    new_draw = is.na(time_diff_min) | time_diff_min > 5,
    draw_id  = cumsum(new_draw)
  ) |>
  group_by(episode_id, draw_id) |>
  slice(1) |>
  ungroup() |>
  arrange(episode_id, sample_date) |>
  group_by(episode_id) |>
  mutate(
    cum_labkemi_draws = row_number(),
    cum_labkemi_vol   = cum_labkemi_draws * ml_per_draw_labkemi
  ) |>
  ungroup() |>
  rename(event_time = sample_date) |>
  dplyr::select(episode_id, event_time, cum_labkemi_draws, cum_labkemi_vol)

# 8c. Cumulative fluid/transfusion events ---------------------------------
# Fluid balance still needs time-based mapping from personal_id to episode
fluid_timed <- vatskebalans_clean |>
  inner_join(
    episodes_48h |>
      dplyr::select(episode_id, personal_id, icu_admission, followup_end),
    by = "personal_id",
    relationship = "many-to-many"
  ) |>
  mutate(
    period_start = day_break_date,
    period_end   = day_break_date + hours(24)
  ) |>
  filter(
    !is.na(period_start),
    period_end > icu_admission,
    period_start < followup_end
  ) |>
  mutate(
    event_time = pmin(period_end, followup_end)
  ) |>
  arrange(episode_id, event_time) |>
  group_by(episode_id) |>
  mutate(
    cum_fluid_balance = cumsum(replace_na(fluid_balance_total, 0)),
    cum_rbc           = cumsum(replace_na(red_blood_cells, 0)),
    cum_plasma        = cumsum(replace_na(plasma, 0))
  ) |>
  ungroup() |>
  dplyr::select(episode_id, event_time, cum_fluid_balance, cum_rbc, cum_plasma)

# =========================================================================
# STEP 9 - Align cumulative exposures to each Hb measurement
# =========================================================================
# Instead of exact timestamp matching, we use the most recent cumulative
# value at or before each Hb measurement time within the same episode.

hb_base <- sepsis_hb_long |>
  mutate(row_id = row_number())

# Blood gas cumulative exposure before Hb time
bg_aligned <- hb_base |>
  dplyr::select(row_id, episode_id, analysis_time) |>
  left_join(
    blodgashb_timed,
    by = join_by(episode_id, analysis_time >= event_time),
    relationship = "many-to-many"
  ) |>
  arrange(row_id, event_time) |>
  group_by(row_id) |>
  slice_tail(n = 1) |>
  ungroup() |>
  dplyr::select(row_id, cum_bloodgas_draws, cum_bloodgas_vol)

# Lab chemistry cumulative exposure before Hb time
lab_aligned <- hb_base |>
  dplyr::select(row_id, episode_id, analysis_time) |>
  left_join(
    labkemi_sampling_timed,
    by = join_by(episode_id, analysis_time >= event_time),
    relationship = "many-to-many"
  ) |>
  arrange(row_id, event_time) |>
  group_by(row_id) |>
  slice_tail(n = 1) |>
  ungroup() |>
  dplyr::select(row_id, cum_labkemi_draws, cum_labkemi_vol)

# Fluid/transfusion cumulative exposure before Hb time
fluid_aligned <- hb_base |>
  dplyr::select(row_id, episode_id, analysis_time) |>
  left_join(
    fluid_timed,
    by = join_by(episode_id, analysis_time >= event_time),
    relationship = "many-to-many"
  ) |>
  arrange(row_id, event_time) |>
  group_by(row_id) |>
  slice_tail(n = 1) |>
  ungroup() |>
  dplyr::select(row_id, cum_fluid_balance, cum_rbc, cum_plasma)

# Join aligned cumulative variables back to Hb-long -----------------------
sepsis_hb_long <- hb_base |>
  dplyr::select(-any_of(c(
    "cum_bloodgas_draws", "cum_bloodgas_vol",
    "cum_labkemi_draws", "cum_labkemi_vol",
    "cum_fluid_balance", "cum_rbc", "cum_plasma",
    "cum_total_draws", "cum_total_vol_ml",
    "cum_total_vol_ml_lag", "cum_fluid_balance_lag",
    "cum_rbc_lag", "cum_plasma_lag",
    "hours_c", "row_id"
  ))) |>
  mutate(row_id = row_number()) |>
  left_join(bg_aligned,    by = "row_id") |>
  left_join(lab_aligned,   by = "row_id") |>
  left_join(fluid_aligned, by = "row_id") |>
  group_by(episode_id) |>
  arrange(analysis_time, .by_group = TRUE) |>
  mutate(
    cum_bloodgas_draws = replace_na(cum_bloodgas_draws, 0),
    cum_bloodgas_vol   = replace_na(cum_bloodgas_vol, 0),
    cum_labkemi_draws  = replace_na(cum_labkemi_draws, 0),
    cum_labkemi_vol    = replace_na(cum_labkemi_vol, 0),
    cum_fluid_balance  = replace_na(cum_fluid_balance, 0),
    cum_rbc            = replace_na(cum_rbc, 0),
    cum_plasma         = replace_na(cum_plasma, 0),
    
    cum_total_draws  = cum_bloodgas_draws + cum_labkemi_draws,
    cum_total_vol_ml = cum_bloodgas_vol   + cum_labkemi_vol,
    
    # Lagged cumulative values: exposure up to PREVIOUS Hb measurement
    cum_total_vol_ml_lag   = lag(cum_total_vol_ml, default = 0),
    cum_fluid_balance_lag  = lag(cum_fluid_balance, default = 0),
    cum_rbc_lag            = lag(cum_rbc, default = 0),
    cum_plasma_lag         = lag(cum_plasma, default = 0),
    
    # Center time at 24h
    hours_c = hours_since_admission - 24
  ) |>
  ungroup() |>
  dplyr::select(-row_id)

cat("\nTime-updated variables coverage in sepsis_hb_long:\n")
sepsis_hb_long |>
  summarise(
    n_obs                  = n(),
    n_with_cum_sample      = sum(!is.na(cum_total_vol_ml)),
    n_with_cum_sample_lag  = sum(!is.na(cum_total_vol_ml_lag)),
    n_with_cum_fluid       = sum(cum_fluid_balance > 0),
    mean_cum_vol           = round(mean(cum_total_vol_ml, na.rm = TRUE), 1),
    mean_cum_vol_lag       = round(mean(cum_total_vol_ml_lag, na.rm = TRUE), 1),
    max_cum_vol            = round(max(cum_total_vol_ml, na.rm = TRUE), 1)
  ) |>
  print()

# =========================================================================
# STEP 10 - Join lab and fluid variables to neuro/SAH group
# =========================================================================
# This enables contrast analyses in RQ3: does sampling associate with Hb
# decline in a group with less inflammation and hemodilution than sepsis?
#
# NOTE: labkemi_clean uses episode_id directly (same as for sepsis).
# vatskebalans_clean still requires time-based mapping from personal_id.

neuro_sah_patients  <- readRDS("data_clean/neuro_sah_patients.rds")
neuro_sah_hb_48h    <- readRDS("data_clean/neuro_sah_hb_48h.rds")

neuro_sah_episodes_48h <- neuro_sah_patients |>
  mutate(
    followup_end = pmin(icu_discharge, icu_admission + hours(48))
  )

# Step 10a - Link lab chemistry to neuro/SAH episodes ---------------------
labkemi_neuro <- labkemi_clean |>
  inner_join(
    neuro_sah_patients |>
      dplyr::select(episode_id, icu_admission, icu_discharge),
    by = "episode_id"
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
  )

cat("\nNeuro/SAH lab linkage:\n")
cat("  Lab rows linked:", nrow(labkemi_neuro), "\n")
cat("  Episodes covered:", n_distinct(labkemi_neuro$episode_id), "\n")

# Step 10b - Admission labs (0-24h) for neuro/SAH -------------------------
neuro_sah_labs_admission <- labkemi_neuro |>
  filter(
    parameter_name %in% lab_params,
    hours_from_admission >= 0,
    hours_from_admission <= 24
  ) |>
  group_by(episode_id, parameter_name) |>
  arrange(hours_from_admission, .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  dplyr::select(episode_id, parameter_name, lab_result) |>
  tidyr::pivot_wider(names_from = parameter_name, values_from = lab_result) |>
  rename(
    crp       = "P-CRP",
    tpk       = "B-Trombocyter",
    lpk       = "B-Leukocyter",
    kreatinin = "P-Kreatinin",
    bilirubin = "P-Bilirubin",
    mch       = "(B)Erc-MCH",
    mcv       = "(B)Erc-MCV"
  )

cat("  Neuro/SAH episodes with admission labs:", nrow(neuro_sah_labs_admission), "\n")

# Step 10c - CRP at 24-48h for neuro/SAH ---------------------------------
crp_24_48h_neuro <- labkemi_neuro |>
  filter(
    parameter_name == "P-CRP",
    hours_from_admission >= 24,
    hours_from_admission <= 48
  ) |>
  group_by(episode_id) |>
  arrange(hours_from_admission, .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  dplyr::select(episode_id, crp_24h = lab_result)

cat("  Neuro/SAH episodes with CRP 24-48h:", nrow(crp_24_48h_neuro), "\n")

# Step 10d - Fluid balance 48h for neuro/SAH ------------------------------
# Still needs time-based mapping because vatskebalans_clean has no episode_id
neuro_sah_fluid_48h <- vatskebalans_clean |>
  inner_join(
    neuro_sah_episodes_48h |>
      dplyr::select(episode_id, personal_id, icu_admission, icu_discharge, followup_end),
    by = "personal_id",
    relationship = "many-to-many"
  ) |>
  mutate(
    period_start = day_break_date,
    period_end   = day_break_date + hours(24)
  ) |>
  filter(
    !is.na(period_start),
    period_end > icu_admission,
    period_start < followup_end
  ) |>
  group_by(episode_id) |>
  summarise(
    fluid_balance_48h   = sum(fluid_balance_total,   na.rm = TRUE),
    colloid_balance_48h = sum(colloid_balance_total, na.rm = TRUE),
    total_balance_48h   = sum(total_balance,         na.rm = TRUE),
    rbc_48h             = sum(red_blood_cells,       na.rm = TRUE),
    colloid_48h         = sum(colloid,               na.rm = TRUE),
    platelets_48h       = sum(platelets,             na.rm = TRUE),
    plasma_48h          = sum(plasma,                na.rm = TRUE),
    n_days              = n(),
    .groups = "drop"
  )

cat("  Neuro/SAH episodes with fluid balance:", nrow(neuro_sah_fluid_48h), "\n")

# Step 10e - Join lab/fluid onto neuro_sah_hb_48h -------------------------
neuro_sah_hb_48h <- neuro_sah_hb_48h |>
  dplyr::select(-any_of(c(
    "crp", "lpk", "tpk", "kreatinin", "bilirubin",
    "mch", "mcv", "crp_24h",
    "fluid_balance_48h", "colloid_balance_48h", "total_balance_48h",
    "rbc_48h", "colloid_48h", "platelets_48h", "plasma_48h", "n_days"
  ))) |>
  left_join(neuro_sah_labs_admission, by = "episode_id") |>
  left_join(crp_24_48h_neuro,         by = "episode_id") |>
  left_join(neuro_sah_fluid_48h,      by = "episode_id")

cat("\nLab coverage in neuro_sah_hb_48h:\n")
neuro_sah_hb_48h |>
  summarise(
    n_crp     = sum(!is.na(crp)),
    n_crp_24h = sum(!is.na(crp_24h)),
    n_lpk     = sum(!is.na(lpk)),
    n_fluid   = sum(!is.na(fluid_balance_48h)),
    n_rbc     = sum(!is.na(rbc_48h))
  ) |>
  print()

# Step 10f - Update comparison_hb_48h with lab/fluid for all groups -------
# Re-build from scratch to ensure clean column state
comparison_hb_48h <- comparison_hb_48h |>
  dplyr::select(-any_of(c(
    "crp", "lpk", "tpk", "kreatinin", "bilirubin",
    "mch", "mcv", "crp_24h",
    "fluid_balance_48h", "colloid_balance_48h", "total_balance_48h",
    "rbc_48h", "colloid_48h", "platelets_48h", "plasma_48h", "n_days"
  ))) |>
  left_join(
    bind_rows(
      # Sepsis: already has all lab/fluid vars in sepsis_hb_48h
      sepsis_hb_48h |>
        dplyr::select(episode_id, crp, lpk, tpk, mch, mcv,
                      kreatinin, bilirubin, crp_24h,
                      fluid_balance_48h, colloid_balance_48h, total_balance_48h,
                      rbc_48h, colloid_48h, platelets_48h, plasma_48h),
      # Neuro/SAH: now has lab/fluid vars joined above
      neuro_sah_hb_48h |>
        dplyr::select(episode_id, crp, lpk, tpk, mch, mcv,
                      kreatinin, bilirubin, crp_24h,
                      fluid_balance_48h, colloid_balance_48h, total_balance_48h,
                      rbc_48h, colloid_48h, platelets_48h, plasma_48h)
    ),
    by = "episode_id"
  )

cat("\nLab/fluid coverage in comparison_hb_48h by diagnosis group:\n")
comparison_hb_48h |>
  group_by(diagnosis_group) |>
  summarise(
    n             = n(),
    has_crp       = sum(!is.na(crp)),
    has_crp_24h   = sum(!is.na(crp_24h)),
    has_fluid     = sum(!is.na(fluid_balance_48h)),
    has_sampling  = sum(!is.na(total_volume_ml)),
    .groups = "drop"
  ) |>
  print()

# Save --------------------------------------------------------------------
saveRDS(comparison_hb_48h,  "data_clean/comparison_hb_48h.rds")
saveRDS(neuro_sah_hb_48h,   "data_clean/neuro_sah_hb_48h.rds")   # now enriched with labs
saveRDS(sepsis_hb_48h,      "data_clean/sepsis_hb_48h.rds")
saveRDS(sepsis_hb_long,     "data_clean/sepsis_hb_long.rds")