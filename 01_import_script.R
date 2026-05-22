# =========== 01_import_script.R =================
# Working directory should be set to the project root.
# If using an RStudio Project (.Rproj), this is set automatically when you open it.
# Otherwise, uncomment and edit the line below:
# setwd("path/to/your/project")

# Packages ----------------------------------------------------------------
library(tidyverse)
library(readxl)

# Read raw data -----------------------------------------------------------

# CIVA Kohort
civakohort_raw <- read_tsv("data/CIVAkohort.csv",
                           locale         = locale(encoding = "Latin1"),
                           show_col_types = FALSE)

# NIVA Kohort
nivakohort_raw <- read.table("data/NIVAkohort.csv",
                             sep              = "\t",
                             header           = TRUE,
                             quote            = "\"",
                             stringsAsFactors = FALSE,
                             fileEncoding     = "latin1")

# Lab chemistry
labkemi_raw <- read.csv("data/labkemi.csv",
                        sep              = ";",
                        fileEncoding     = "UTF-8-BOM",
                        stringsAsFactors = FALSE)

# Fluid balance
vatskebalans_raw <- read_excel("data/vatskebalans.xlsx", sheet = 1)

# Blood gas Hb 
blodgashb_raw <- read_excel("data/vatskebalans.xlsx", sheet = 2)

# Validation checks -------------------------------------------------------
cat("Raw data loaded:\n")
cat("  civakohort:   ", nrow(civakohort_raw),  "rows\n")
cat("  nivakohort:   ", nrow(nivakohort_raw),  "rows\n")
cat("  labkemi:      ", nrow(labkemi_raw),      "rows\n")
cat("  vatskebalans: ", nrow(vatskebalans_raw), "rows\n")
cat("  blodgashb:    ", nrow(blodgashb_raw),    "rows\n")
