# Family History of Cardiovascular Disease, Cardiovascular Health, and CVD Mortality

## Overview

This repository provides the R code used in the following study:

> **Family history of cardiovascular disease, cardiovascular health, and CVD mortality: a population-based cohort study of young and middle-aged Korean adults**  
> Sujeong Shin†\*, Sojin Kim†, Yoosoo Chang, Eunju Sung, Jae-Heon Kang\*  
> † These authors contributed equally to this work.  
> Corresponding authors: Sujeong Shin (sj1115.shin@samsung.com), Jae-Heon Kang (jaeheon.kang@samsung.com)

---

## Abstract

Family history of cardiovascular disease (CVD) and lifestyle behaviours represent complementary dimensions of cardiovascular risk, yet their joint contribution to CVD mortality remains limited. We aimed to investigate their joint association with CVD mortality in a large cohort of young and middle-aged Korean adults.

We analyzed **476,965 participants** in the Kangbuk Samsung Health Study (2002–2022). Lifestyle status was assessed using a six-component **modified Cardiovascular Health (mCVH) score** and categorized as high or low based on tertiles. Participants were classified into four groups by mCVH category and family history status. CVD mortality was ascertained via the national death registry. Hazard ratios (HRs) and 95% confidence intervals (CIs) were estimated using Cox proportional hazards models.

During a median follow-up of 5.72 years, 258 CVD deaths occurred. Both a low mCVH score and a positive family history were independently associated with elevated CVD mortality risk. Compared with those with a high mCVH score and no family history, participants with a low mCVH score and a positive family history had the highest risk (HR 2.54; 95% CI, 1.33–4.87). No significant multiplicative interaction was observed (p for interaction = 0.641).

**Keywords:** cardiovascular mortality; family history; healthy lifestyle; cardiovascular health score; cohort study; Korean adults

---

## Study Design

| Item | Details |
|---|---|
| Design | Population-based retrospective cohort study |
| Population | 476,965 Korean adults (Kangbuk Samsung Health Study, 2002–2022) |
| Exposure | Modified Cardiovascular Health (mCVH) score (6-component: blood pressure, blood glucose, blood lipids, BMI, smoking, physical activity) × family history of CVD |
| Outcome | CVD mortality (ICD-10 codes I00–I99) via national death registry |
| Statistical method | Cox proportional hazards models with interaction terms |
| Median follow-up | 5.72 years |

---

## Key Results

| Group | HR (95% CI) |
|---|---|
| High mCVH, No family history (reference) | 1.00 |
| High mCVH, Positive family history | — |
| Low mCVH, No family history | — |
| Low mCVH, Positive family history | 2.54 (1.33–4.87) |

Fully adjusted model accounting for age, sex, examination center, year, education, LDL cholesterol, history of hypertension, and history of dyslipidemia.  
p for interaction = 0.641 (no significant multiplicative interaction).

---

## Repository Structure

```
├── CVD_analysis.R      # Main analysis: preprocessing, mCVH score, Cox regression, adjusted survival curves
├── CVD_Baseline.R      # Baseline characteristics table
├── CVD_subgroup.R      # Subgroup analyses (smoking, BMI, BP, HbA1c, lipids)
├── KMcurve_fig2.R      # Figure 2: Adjusted survival curves (whole cohort)
├── KMcurve_fig3.R      # Figure 3: Adjusted survival curves (subgroup panels)
└── README.md
```

### Scripts

- **`CVD_analysis.R`**: Data preprocessing, exclusion criteria, mCVH score construction (LE6: blood pressure, HbA1c, blood lipids, nicotine, BMI, physical activity), Cox proportional hazards models with interaction terms (mCVH group × family history), and generation of adjusted survival curve data.
- **`CVD_Baseline.R`**: Descriptive statistics and baseline characteristics table by mCVH group.
- **`CVD_subgroup.R`**: Subgroup analyses stratified by smoking status, BMI, blood pressure, HbA1c level, and non-HDL cholesterol.
- **`KMcurve_fig2.R`**: Visualization of adjusted cumulative event curves for the four groups (whole cohort; Figure 2).
- **`KMcurve_fig3.R`**: Visualization of adjusted cumulative event curves across subgroup panels (Figure 3).

---

## Requirements

```r
install.packages(c("readr", "dplyr", "tidyr", "lubridate", "stringr",
                   "survival", "survminer", "rms",
                   "VIM", "mice",
                   "ggplot2", "patchwork", "readxl", "grid"))
```

| Package | Purpose |
|---|---|
| `dplyr`, `tidyr`, `readr`, `lubridate`, `stringr` | Data manipulation |
| `survival`, `survminer` | Cox proportional hazard models and survival curves |
| `rms` | Restricted cubic spline analyses |
| `VIM`, `mice` | Missing data visualization and multiple imputation |
| `ggplot2`, `patchwork`, `grid` | Visualization |
| `readxl` | Reading Excel-based results for figure generation |

---

## Data Availability

The dataset used in this study is from the **Kangbuk Samsung Health Study** and is not publicly available due to privacy and ethical restrictions. Researchers interested in data access should contact the corresponding authors:

**Sujeong Shin, MD, MPH**  
Department of Family Medicine, Kangbuk Samsung Hospital, Sungkyunkwan University School of Medicine  
29, Saemunan-ro, Jongno-gu, Seoul 03181, Republic of Korea  
Email: sj1115.shin@samsung.com

**Jae-Heon Kang, MD, PhD**  
Department of Family Medicine, Kangbuk Samsung Hospital, Sungkyunkwan University School of Medicine  
29, Saemunan-ro, Jongno-gu, Seoul 03181, Republic of Korea  
Email: jaeheon.kang@samsung.com

---

## Authors

| Name | ORCID | Affiliation |
|---|---|---|
| Sujeong Shin, MD, MPH\* | [0000-0002-7661-0935](https://orcid.org/0000-0002-7661-0935) | Department of Family Medicine, Kangbuk Samsung Hospital, Sungkyunkwan University School of Medicine |
| Sojin Kim, PhD | [0009-0006-5023-3189](https://orcid.org/0009-0006-5023-3189) | Human-Centered Artificial Intelligence Research Institute, Ewha Womans University |
| Yoosoo Chang, MD, PhD | — | Center for Cohort Studies & Department of Occupational and Environmental Medicine, Kangbuk Samsung Hospital |
| Eunju Sung, MD, PhD | — | Department of Family Medicine, Kangbuk Samsung Hospital |
| Jae-Heon Kang, MD, PhD\* | — | Department of Family Medicine, Kangbuk Samsung Hospital |

\* Corresponding authors  
† S.S. and S.K. contributed equally to this work.

---

## Author Contributions

S.S. and S.K. conceptualized the study and designed the methodology. S.K. performed formal statistical analysis and prepared the original draft. All authors reviewed and approved the final manuscript.

---

## Citation

> Shin S†, Kim S†, Chang Y, Sung E, Kang JH. Family history of cardiovascular disease, cardiovascular health, and CVD mortality: a population-based cohort study of young and middle-aged Korean adults. (under review).
