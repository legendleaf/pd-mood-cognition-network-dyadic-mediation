# Parkinson's Disease: Mood/Cognition, Quality of Life and Caregiver Burden

Analysis code for *Mood and Cognitive Symptoms Are Associated With Health-Related Quality of Life and Caregiver Burden in Parkinson's Disease: A Cross-Sectional Study*. This version includes the R1 analyses for 165 patients and 53 linked patient-caregiver pairs.

Participant data are not included. Access to de-identified study data requires approval from the corresponding author and remains subject to consent and privacy restrictions. The code can be shared separately for anonymous review.

## Run the analyses

Run from the repository root. Keep approved CSVs and workbooks in a separate private directory.

```r
Sys.setenv(
  PD_DATA_DIR = "/path/to/approved/private/csvs",
  PD_OUTPUT_DIR = "/path/to/local/results",
  PD_ITEM_WORKBOOK = "/path/to/finalized/patient_items.xlsx"
)
source("script/Step1_NetworkAnalysis.R")
source("script/Step2_Advanced_Analysis.R")
source("script/Step3_ledd_medication_analysis_serial_mediation.R")
source("script/Step4_PatientSensitivity.R")
source("script/Step5_TableChecks.R")
source("script/Step6_NetworkSensitivity.R")
```

Each step sources `Step0_DataPreparation.R`. It reads `analysis_patients.csv` and, when present, `analysis_pairs.csv`, checks unique IDs and patient-pair agreement, and computes derived scores in memory. Output must be outside the input directory. Files read by each step are checked for changes before the step finishes.

The finalized patient workbook is needed with the original CSVs. PDQ-39 items use codes 1-5; subtracting 1 gives scores 0-4. Their sum must match the workbook total. The loader corrects the documented constant social-domain offset in memory and computes SI as the unweighted mean of all eight domain scores. The raw item total and its 0-100 transform remain separate quantities. No CSV or workbook is overwritten.

The workbook may be omitted with approved R1 CSVs whose domain scores have already been checked. A legacy item-weighted SI is rejected without the workbook. The six-domain mean excludes emotional well-being and cognition; it is a sensitivity index, not a standard PDQ-39 SI.

Walking/standing and balance responses come from the workbook or an external `PD_FUNCTION_CSV` with `ID`, `independent_walk` and `balance_obstacle` (1=yes, 2=no). Without these inputs, Step 2 warns that the two functional-covariate models were omitted. These responses are not used to infer Hoehn and Yahr stages.

## Analyses and outputs

| Script | Analysis | Manuscript location |
| --- | --- | --- |
| Step 1 | 16-node network, strength, three-community bridge strength, case-dropping stability | Figure 1; Figure S1 |
| Step 2 | SES comparison; four patient-centered joint-outcome models; constructed differences | Figure 2; Table 3; Tables S8-S9 |
| Step 3 | Adjusted cross-sectional serial model; VIFs for all three equations | Figure 3; Figure S2; Table 4 |
| Step 4 | Full-cohort SI and six-domain regressions with HC3 uncertainty | Tables S11A-S11B |
| Step 5 | Distributions, functional descriptions, item completeness and WOQ count checks | Tables S2A-S2B, S4, S10 |
| Step 6 | Network intervals, mood/cognition contrasts and 14-node sensitivity | Tables S5-S7 |
| ItemReliability.py | Item-level internal consistency | Table S3 |

Outputs are analysis tables and figure components, not Word tables or final Adobe layouts. Original plot colors are retained. Figure S3's submission counts come from the documented recruitment summary and cannot be reconstructed from final CSVs alone; this code does not draw the flow diagram.

Network intervals use 1,000 participant resamples; case-dropping uses 500 subsets. Communities follow questionnaire source: seven PDQ-39, eight NMSQuest and one WOQ-9 node. Bridge strength is unnormalized cross-community connectivity. The reduced network removes PDQ emotional well-being and cognition. Centrality contrasts are exploratory and unadjusted. The stability plot's horizontal axis shows cases retained.

SES outputs distinguish the group strengths from the `cor_auto` networks and the strengths used internally by NCT. The NCT statistic and permutation p value belong to the latter; the two sets of strengths must not be substituted for each other. Neither comparison establishes equivalence from a nonsignificant result.

Dyadic intervals use 2,000 resamples; the serial model uses 5,000. Continuous variables are standardized, sex remains binary, and the dyadic models allow residual covariance between patient SI and caregiver ZBI. These are patient-centered joint-outcome models, not reciprocal APIMs. Serial products describe cross-sectional indirect associations. Proportions mediated are no longer reported.

Constructed differences compare caregiver DASS total with the patient PDQ communication/emotional-well-being mean or NMS mood/cognition count. Absolute-difference and absolute-residual terms enter after both linear components. Their coefficients, intervals and nested-model F tests are reported separately. Figure 2 lines show unadjusted trends; its displayed increments are component-adjusted tests. Bootstrap samples retain full-sample standardized scores and constructed indices, including the residual-generating regression: 2,000 draws for absolute differences and 1,000 for absolute residuals.

Patient regressions adjust for age, sex and disease duration without WOQ-9 or LEDD. HC3 uncertainty uses two-sided t intervals with residual degrees of freedom. Holm adjustment covers the two mood/cognition tests only.

## Item reliability

Set `PD_CAREGIVER_WORKBOOK` to the finalized caregiver workbook, then run in a shell with the same environment variables:

```text
python script/ItemReliability.py
```

This script also reads `analysis_caregivers.csv` and requires NumPy, pandas and openpyxl. Whole-participant percentile intervals use 2,000 draws. PDQ/NMS/ZBI/DASS retain the original NumPy seed and domain order (20260912); WOQ uses seed 20260930. WOQ code 1 means a symptom is present and improves after medication; all other recorded categories are nonpositive. Retained totals are checked against item counts. Binary alpha equals KR-20; internal consistency does not establish diagnostic validity or unidimensionality.

## Software and reproducibility

Required R packages: `dplyr`, `readr`, `ggplot2`, `qgraph`, `bootnet`, `networktools`, `NetworkComparisonTest`, `lavaan`, `boot`, `patchwork`, `sandwich` and `readxl` (for workbooks). Install them before running. Optional `PD_R_LIBRARY` supplies local library directories separated by the platform path separator.

Each step records its actual `sessionInfo()`. Initial analyses and revision checks used different R/package environments. A fresh bootstrap can differ slightly from an archived interval; point estimates and model specifications are checked separately.

Optional private archives can reproduce archived intervals:

- `PD_NETWORK_ARCHIVE`: original Step 1 RData. Node order, correlations, graph estimates and draw counts must match.
- `PD_MODEL_ARCHIVE_DIR`: revised model RDS files named `duration_only`, `age_sex_duration`, `add_balance_difficulty`, `add_inability_walk_stand` and `serial_corrected_si`. Standardized data and fresh coefficients must match before reuse.

These archives contain participant-level information and must remain private. Without them, models are refitted. Revision SEMs and the reduced-network bootstrap use seed 20260912; original network/SES and discrepancy reruns use seed 2024. Archived intervals preserve their original random streams.

Resample-count overrides `PD_NETWORK_BOOTSTRAP`, `PD_CASE_BOOTSTRAP`, `PD_NCT_PERMUTATIONS`, `PD_DYAD_BOOTSTRAP`, `PD_SERIAL_BOOTSTRAP`, `PD_DIFFERENCE_BOOTSTRAP`, `PD_RESIDUAL_BOOTSTRAP`, `PD_REDUCED_BOOTSTRAP` and `PD_RELIABILITY_BOOTSTRAP` are for testing. `PD_POINT_ONLY=1` disables SEM bootstrap uncertainty. Test outputs must not replace full-resample manuscript results.

## Sharing the code

Inputs, models, tables, figures and logs are ignored by Git. RDS/RData objects and scatter plots can contain individual observations. Do not upload generated outputs without a separate privacy review.

```text
python tools/release_code.py --output /path/to/anonymous-code.zip
```

The exporter uses an explicit file list, checks common identity/path/secret patterns and omits `.git`, data and outputs. ZIP metadata are normalized. This separates files from Git authorship and remote-account details; it does not anonymize the existing repository or its history. Automated checks cannot certify removal of every possible identifier.
