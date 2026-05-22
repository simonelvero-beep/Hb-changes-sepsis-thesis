# ============== run_all.R ================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")

cat("Starting analysis pipeline:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Step 1/8: Importing data...\n")
source("01_import_script.R")
cat("✓ Import complete\n\n")

cat("Step 2/8: Cleaning data...\n")
source("02_clean_data.R")
cat("✓ Clean complete\n\n")

cat("Step 3/8: Deriving labkemi volume assumption...\n")
source("02b_volume_assumption.R")
cat("✓ Volume assumption complete\n\n")

cat("Step 4/8: Merging sepsis Hb data...\n")
source("03a_merge_sepsis_hb.R")
cat("✓ Sepsis Hb merge complete\n\n")

cat("Step 5/8: Merging comparison data...\n")
source("03b_merge_comparison.R")
cat("✓ Comparison merge complete\n\n")

cat("Step 6/8: Merging inflammation data...\n")
source("03c_merge_inflammation.R")
cat("✓ Inflammation merge complete\n\n")

cat("Step 7/8: Descriptive statistics...\n")
source("04_descriptive_stat.R")
cat("✓ Descriptive stats complete\n\n")

cat("Step 8/8: Analysis...\n")
source("05_analysis.R")
cat("✓ Analysis complete\n\n")

cat("Pipeline complete:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
