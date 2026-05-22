# ============== 02_clean_data.R ===========================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")

# Packages ----------------------------------------------------------------
library(tidyverse)
library(janitor)
library(lubridate)    # for ymd_hm()

# Create output folder ----------------------------------------------------
dir.create("data_clean", showWarnings = FALSE)

# Help functions ----------------------------------------------------------
trim_all_character <- function(df) {
  df |>
    mutate(across(where(is.character), str_trim))
}

empty_to_na <- function(df) {
  df |>
    mutate(across(where(is.character), ~na_if(., "")))
}

standardize_pnr <- function(pnr) {
  pnr <- str_trim(pnr)
  case_when(
    str_detect(pnr, "^\\d{12}$")       ~ pnr,
    str_detect(pnr, "^\\d{8}-\\d{4}$") ~ str_remove(pnr, "-"),
    str_detect(pnr, "^\\d{6}-\\d{4}$") ~ paste0(ifelse(as.numeric(str_sub(pnr, 1, 2)) > 24, "19", "20"), str_remove(pnr, "-")),
    str_detect(pnr, "^\\d{10}$")       ~ paste0(ifelse(as.numeric(str_sub(pnr, 1, 2)) > 24, "19", "20"), pnr),
    TRUE ~ NA_character_
  )
}
#===================Data cleanup===================================
#--------------------- CIVA Kohort ---------------------
civakohort_clean <- civakohort_raw |>
  clean_names() |>
  trim_all_character() |>
  empty_to_na() |>
  rename(personal_id     = "personnummer",
         icu_admission   = "inskriven",
         icu_discharge   = "utskriven",
         episode_id      = "id",
         patient_outcome = "vardresultat") |>
  mutate(
    personal_id     = standardize_pnr(personal_id),
    episode_id      = as.character(episode_id),
    icu_admission   = as.POSIXct(icu_admission, tz = "Europe/Stockholm"),
    icu_discharge   = as.POSIXct(icu_discharge, tz = "Europe/Stockholm"),
    patient_outcome = factor(patient_outcome,
                             levels = c("Levande", "Avliden"),
                             labels = c("alive", "dead"))
  ) |>
  dplyr::select(episode_id, personal_id, icu_admission, icu_discharge,
                patient_outcome, ivadiagnos_1, ivadiagnos_2, ivadiagnos_3,
                diagnos_1, diagnos_2, diagnos_3) |>
  filter(!is.na(personal_id))

# NIVA Kohort -------------------------------------------------------------
nivakohort_clean <- nivakohort_raw |>
  clean_names() |>
  trim_all_character() |>
  empty_to_na() |>
  rename(personal_id     = "personnummer",
         icu_admission   = "inskriven",
         icu_discharge   = "utskriven",
         episode_id      = "id",
         patient_outcome = "vardresultat") |>
  mutate(
    personal_id     = standardize_pnr(personal_id),
    episode_id      = as.character(episode_id),
    icu_admission   = ymd_hm(icu_admission, tz = "Europe/Stockholm"),
    icu_discharge   = ymd_hm(icu_discharge, tz = "Europe/Stockholm"),
    patient_outcome = factor(patient_outcome,
                             levels = c("Levande", "Avliden"),
                             labels = c("alive", "dead"))
  ) |>
  dplyr::select(episode_id, personal_id, icu_admission, icu_discharge,
                patient_outcome, ivadiagnos_1, ivadiagnos_2, ivadiagnos_3,
                diagnos_1, diagnos_2, diagnos_3) |>
  filter(!is.na(personal_id))

# Labkemi - lab values ----------------------------------------------------
labkemi_clean <- labkemi_raw |>
  clean_names() |>
  trim_all_character() |>
  empty_to_na() |>
  rename(personal_id    = "personnummer",
         lab_result     = "resultat",
         parameter_name = "analys",
         unit           = "mattenhet",
         sample_date    = "provtagningsdatum",
         episode_id     = "id") |>
  mutate(
    personal_id    = standardize_pnr(personal_id),
    episode_id     = as.character(episode_id),
    sample_date    = as.POSIXct(sample_date,
                                format = "%Y-%m-%d %H:%M:%S",
                                tz     = "Europe/Stockholm"),
    lab_result     = as.numeric(str_replace(lab_result, ",", ".")),
    unit           = as.factor(unit),
    parameter_name = as.factor(parameter_name)
  ) |>
  dplyr::select(episode_id, personal_id, sample_date,
                parameter_name, lab_result, unit) |>
  filter(!is.na(personal_id)) |>
  filter(!is.na(lab_result))

# Labkemi - blood sampling counts -----------------------------------------
labkemi_sampling <- labkemi_raw |>
  clean_names() |>
  trim_all_character() |>
  empty_to_na() |>
  rename(personal_id    = "personnummer",
         lab_result     = "resultat",
         parameter_name = "analys",
         sample_date    = "provtagningsdatum",
         episode_id     = "id") |>
  mutate(
    personal_id = standardize_pnr(personal_id),
    episode_id  = as.character(episode_id),
    sample_date = as.POSIXct(sample_date,
                             format = "%Y-%m-%d %H:%M:%S",
                             tz     = "Europe/Stockholm")
  ) |>
  filter(!lab_result %in% c("Ej taget", "Provtagn fel",
                            "Aggregat", "Se kommentar", "Se kom,")) |>
  filter(!is.na(personal_id)) |>
  filter(!is.na(sample_date)) |>
  dplyr::select(episode_id, personal_id, sample_date, parameter_name)

# Vatskebalans ------------------------------------------------------------
vatskebalans_clean <- vatskebalans_raw |>
  clean_names() |>
  trim_all_character() |>
  empty_to_na() |>
  mutate(hospital_number = as.character(hospital_number)) |>
  rename(personal_id = "hospital_number") |>
  mutate(
    personal_id = standardize_pnr(personal_id),
    across(where(is.character) & !c(personal_id, day_break_date), as.numeric)
  ) |>
  dplyr::select(
    personal_id,
    day_break_date,
    fluid_balance_total  = vatskebalans_in_ut_totalt,
    colloid_balance_total = blod_kolloidbalans_in_ut_totalt,
    total_balance        = totalbalans,
    red_blood_cells      = erytrocyter,
    colloid              = kolloid,
    platelets            = trombocyter,
    plasma               = plasma       # ← removed trailing comma
  ) |>
  filter(!is.na(personal_id))

# Blood gas Hb (renamed from hbvatskebalans) ------------------------------
blodgashb_clean <- blodgashb_raw |>     # ← renamed from hbvatskebalans_raw
  clean_names() |>
  mutate(
    time = as.POSIXct(as.numeric(time) * 86400,
                      origin = "1899-12-30",
                      tz     = "UTC"),
    time = force_tz(time, tzone = "Europe/Stockholm")
  ) |>
  trim_all_character() |>
  empty_to_na() |>
  mutate(hospital_number = as.character(hospital_number)) |>
  rename(personal_id   = "hospital_number",
         analysis_time = "time",
         lab_result    = "value",
         unit          = "unit_name") |>
  mutate(
    personal_id = standardize_pnr(personal_id),
    lab_result  = as.numeric(lab_result)
  ) |>
  filter(!is.na(personal_id))

# Verify ------------------------------------------------------------------
cat("Clean data row counts:\n")
cat("  civakohort:   ", nrow(civakohort_clean),   "\n")
cat("  nivakohort:   ", nrow(nivakohort_clean),   "\n")
cat("  labkemi:      ", nrow(labkemi_clean),       "\n")
cat("  labkemi_samp: ", nrow(labkemi_sampling),    "\n")
cat("  vatskebalans: ", nrow(vatskebalans_clean),  "\n")
cat("  blodgashb:    ", nrow(blodgashb_clean),     "\n")

# Save --------------------------------------------------------------------
saveRDS(civakohort_clean,    "data_clean/civakohort_clean.rds")
saveRDS(nivakohort_clean,    "data_clean/nivakohort_clean.rds")
saveRDS(labkemi_clean,       "data_clean/labkemi_clean.rds")
saveRDS(labkemi_sampling,    "data_clean/labkemi_sampling.rds")   # ← added
saveRDS(vatskebalans_clean,  "data_clean/vatskebalans_clean.rds")
saveRDS(blodgashb_clean,     "data_clean/blodgashb_clean.rds")    # ← renamed