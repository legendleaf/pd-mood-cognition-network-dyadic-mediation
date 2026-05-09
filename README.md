# Mood/Cognition Symptoms in Parkinson's Disease

R code and de-identified data for the manuscript *Mood/Cognition Symptoms in Parkinson's Disease: Associations With Patient Quality of Life and Caregiver Burden*. See the manuscript for study details and references.

## Files

```
scripts/
  Step1_NetworkAnalysis.R                            # 16-node EBICglasso symptom network
  Step2_Advanced_Analysis.R                          # NCT, APIM, dyadic mismatch
  Step3_ledd_medication_analysis_serial_mediation.R  # serial mediation
data/
  analysis_patients.csv     # 165 patients x 51 columns
  analysis_caregivers.csv   # 53 caregivers x 6 columns
  analysis_pairs.csv        # 53 dyads, patient + caregiver columns joined
```

Each participant is identified by a sequential numeric `ID`. The patient file contains demographics, socioeconomic indicators, lifestyle variables, disease duration, LEDD, EQ-5D, the eight PDQ-39 domains plus the Summary Index, the eight NMSQuest domains plus total, the WOQ-9 total score, and the AD8. The caregiver file contains the ZBI total and the three DASS-21 subscales. The dyadic file is the patient and caregiver tables joined on `ID` for the 53 patients with a returning caregiver.

## Analysis steps

1. **Step 1 — Symptom network.** Estimates the 16-node EBICglasso network (PDQ-39, NMSQuest, WOQ-9), computes strength and bridge centrality, and assesses stability with non-parametric (1,000) and case-dropping (500) bootstraps. Generates the panels of Figure 1, Figure S1, and the centrality values reported in Table 2.

2. **Step 2 — Dyadic and secondary analyses.** Runs (i) the Network Comparison Test for SES moderation on a 10-node core network, (ii) a patient-centered actor–partner model linking patient mood/cognition symptoms to patient HRQoL and caregiver burden, and (iii) caregiver–patient discrepancy analyses (incremental ΔR², directional Spearman ρ, residualized sensitivity). Generates Figure 2 and Table 3.

3. **Step 3 — Serial mediation.** Fits the serial mediation LEDD → WOQ-9 → NMS mood/cognition → PDQ-39 Summary Index in `lavaan` with 5,000 percentile bootstrap resamples, controlling for age, sex, and disease duration. Generates Figure 3, Figure S2, and Table 4.

## Run order

Run from the repository root. Step 2 reads `step1_results.RData` written by Step 1.

```r
source("scripts/Step1_NetworkAnalysis.R")
source("scripts/Step2_Advanced_Analysis.R")
source("scripts/Step3_ledd_medication_analysis_serial_mediation.R")
```

Figures and tables are written under `results_network_analysis/`.

## Software

R 4.5.1 with `tidyverse`, `qgraph`, `bootnet`, `networktools`, `NetworkComparisonTest`, `lavaan`, `boot`, `car`, `showtext`, `patchwork`.

## Data access

De-identified analysis data are provided here. Raw participant-level questionnaire responses are available from the corresponding author on reasonable request.
