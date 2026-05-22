# ================== 02b_volume_assumption.R ============================
# Empirically calculates the labkemi volume per draw assumption
# based on tube types present at each sampling occasion.
#
# Run this BEFORE setting ml_per_draw_labkemi in 03b and 03c.
# Output: recommended ml_per_draw_labkemi value with justification.
# ========================================================================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")

library(dplyr)

# Load data ---------------------------------------------------------------
labkemi_clean   <- readRDS("data_clean/labkemi_clean.rds")
sepsis_patients <- readRDS("data_clean/sepsis_patients.rds")

# Tube volume assumptions (mL) --------------------------------------------
# retrieved from sampling instructions doc.plus region Uppsala
vol_edta   <- 4.0   # purple tube — CBC, differential
vol_sst    <- 4.5   # mintgreen tube — chemistry, CRP, bilirubin
vol_citrat <- 3.5   # blue tube - citrate — coagulation (PK, APTT)

# Tube classification -----------------------------------------------------
tube_classification <- labkemi_clean |>
  inner_join(
    sepsis_patients |> select(episode_id, icu_admission, icu_discharge),
    by = "episode_id"
  ) |>
  filter(sample_date >= icu_admission, sample_date <= icu_discharge) |>
  mutate(
    tube = case_when(
      parameter_name %in% c(
        "B-Trombocyter", "B-Leukocyter", "B-Hemoglobin[Hb]",
        "(B)Erc-MCV", "(B)Erc-MCH", "(B)Erc-MCHC", "B-Retikulocyter"
      ) ~ "EDTA",
      parameter_name %in% c(
        "P-Kreatinin", "P-CRP", "P-Bilirubin", "P-Bilirubin, konj."
      ) ~ "SST",
      parameter_name %in% c(
        "P-Protrombinkomplex[PK]", "P-APT tid", "B-Trombocyter Citrat"
      ) ~ "Citrat",
      TRUE ~ "Other"
    )
  )

# Calculate tube presence per sampling occasion ---------------------------
tube_per_occasion <- tube_classification |>
  group_by(episode_id, sample_date) |>
  summarise(
    has_edta   = any(tube == "EDTA"),
    has_sst    = any(tube == "SST"),
    has_citrat = any(tube == "Citrat"),
    has_other  = any(tube == "Other"),
    n_params   = n(),
    .groups = "drop"
  )

# Summary -----------------------------------------------------------------
volume_summary <- tube_per_occasion |>
  summarise(
    n_occasions        = n(),
    pct_edta           = round(mean(has_edta)   * 100, 1),
    pct_sst            = round(mean(has_sst)    * 100, 1),
    pct_citrat         = round(mean(has_citrat) * 100, 1),
    pct_other          = round(mean(has_other)  * 100, 1),
    median_params      = round(median(n_params), 1),
    mean_params        = round(mean(n_params), 1),
    mean_vol_per_draw  = round(mean(
      has_edta * vol_edta + has_sst * vol_sst + has_citrat * vol_citrat
    ), 1)
  )

cat("=== Empirical volume assumption ===\n\n")
cat("Tube volume assumptions used:\n")
cat("  EDTA (purple):  ", vol_edta,   "mL\n")
cat("  SST (yellow):   ", vol_sst,    "mL\n")
cat("  Citrat (blue):  ", vol_citrat, "mL\n\n")

cat("Tube presence per sampling occasion (sepsis episodes, 0-48h):\n")
cat("  Total occasions:      ", volume_summary$n_occasions, "\n")
cat("  With EDTA tube:       ", volume_summary$pct_edta,   "%\n")
cat("  With SST tube:        ", volume_summary$pct_sst,    "%\n")
cat("  With Citrat tube:     ", volume_summary$pct_citrat, "%\n")
cat("  Median params/draw:   ", volume_summary$median_params, "\n")
cat("  Mean params/draw:     ", volume_summary$mean_params, "\n\n")

cat(">>> Recommended ml_per_draw_labkemi:",
    volume_summary$mean_vol_per_draw, "mL <<<\n\n")

cat("Use this value in 03b_merge_comparison.R and 03c_merge_inflammation.R:\n")
cat("  ml_per_draw_labkemi <-", volume_summary$mean_vol_per_draw, "\n\n")