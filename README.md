# Hb-changes-sepsis-thesis

R code accompanying the medical thesis:

**Early Hemoglobin Decline in Sepsis: An Observational Cohort Study of Inflammation and Blood Sampling**

Simon Elverö, Uppsala University, *Läkarprogrammet, självständigt arbete (30 hp)*, 2026.
Supervisor: Miklós Lipcsey.

---

## About this repository

This repository contains the R code used for data processing, statistical analysis, and visualisation in the above thesis. The study is a retrospective observational cohort analysis of haemoglobin trajectories during the first 48–60 hours of ICU admission in patients with sepsis, with neurotrauma and subarachnoid haemorrhage patients as a comparison cohort.

**No patient data are included in this repository.** All scripts read from local data files that are not shared. The code is provided for transparency and reproducibility of the analytical pipeline.

The full console output of the pipeline is included as `ANALYSRESULTAT.txt`.

---

## Pipeline

Scripts run in numerical order:

| Script | Purpose |
|---|---|
| `01_import_script.R` | Imports the five raw data sources (CIVA, NIVA, Labkemi, Vatskebalans, BlodgasHb) and applies basic type coercion. |
| `02_clean_data.R` | Variable cleaning, harmonisation, timestamp standardisation to Stockholm time, derivation of age and sex from the Swedish personal identifier. |
| `02b_volume_assumption.R` | Derives the per-occasion clinical-chemistry sampling volume (≈10 mL) from the observed tube-type composition in the cleaned dataset. |
| `03a_merge_sepsis_hb.R` | Episode-level merging for the sepsis Hb-trajectory dataset (research question 1). |
| `03b_merge_comparison.R` | Episode-level merging for the three-group diagnosis comparison (research question 2). |
| `03c_merge_inflammation.R` | Episode-level merging for the inflammation, fluid balance, and sampling exposures, including lagged time-updated exposures (research question 3). |
| `04_descriptive_stat.R` | Table 1 and descriptive analyses. |
| `05_analysis.R` | Research-question analyses, multivariable regressions, and the primary linear mixed-effects model. |

`ANALYSRESULTAT.txt` is the captured console output from the full pipeline.

---

## Data sources

The five raw datasets are extracted from the routinely-collected ICU databases at Uppsala University Hospital:

- **CIVA** and **NIVA** — ICU registry data from the Central and Neuro intensive care units, providing episode identifiers, admission and discharge times, up to six ICD-coded diagnoses per episode, and patient outcome.
- **Labkemi** — clinical chemistry results with sampling time, sample type, and episode identifiers.
- **Vatskebalans** — daily fluid balance (07:00–07:00), including transfusion volumes.
- **BlodgasHb** — point-of-care blood-gas haemoglobin measurements; linked to ICU episodes using personal identifier and timestamp.

Data are not included in this repository. Researchers wishing to reproduce the analyses with new data would need equivalent extracts and would need to adapt the file paths in `01_import_script.R`.

---

## Cohort definitions

- **Sepsis cohort:** CIVA episodes with ICD-10-SE codes A40, A41, R57.2, or R65.1.
- **Comparison cohort:** CIVA and NIVA episodes with ICD-10-SE codes S06 (neurotrauma) or I60 (subarachnoid haemorrhage).
- Readmissions within 30 days of a previous ICU admission were excluded.

Final analytic cohorts: 859 sepsis episodes, 522 neurotrauma episodes, 336 SAH episodes.

---

## Software and packages

- R version 4.5.2
- Main packages: `dplyr`, `tidyr`, `gtsummary`, `rstatix`, `broom`, `ggplot2`, `ggpubr`, `patchwork`, `lme4`, `lmerTest`, `car`

---

## Ethics and data handling

This work was conducted as a quality improvement collaboration with the ICU department at Uppsala University Hospital. Under Swedish regulations, internal quality improvement work based on a healthcare provider's own data does not require approval from the Swedish Ethical Review Authority (Etikprövningsmyndigheten). All data were handled in accordance with the General Data Protection Regulation on secure servers on hospital-based computers. No patient-identifiable information appears in this repository.

---

## Citation

If you refer to this code, please cite the thesis:

> Elverö, S. (2026). *Early Hemoglobin Decline in Sepsis: An Observational Cohort Study of Inflammation and Blood Sampling.* Medical thesis (självständigt arbete 30 hp), Uppsala University.

---

## Contact

Simon Elverö — [your email, if you want to include one]
