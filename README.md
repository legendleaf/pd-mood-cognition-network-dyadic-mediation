# Mood/Cognition Symptoms in Parkinson's Disease

This repository contains the R analysis code for the manuscript *Mood/Cognition Symptoms in Parkinson's Disease: Associations With Patient Quality of Life and Caregiver Burden*.

Participant-level datasets are not included in this repository because of participant privacy considerations. De-identified datasets generated and/or analysed during the current study are available from the corresponding author upon reasonable request.

## Files

```text
script/Step1_NetworkAnalysis.R
  Estimates the 16-node EBICglasso symptom network, computes strength and bridge centrality, and assesses network stability.

script/Step2_Advanced_Analysis.R
  Runs SES network comparison, patient-centered dyadic modelling, and caregiver-patient discrepancy analyses.

script/Step3_ledd_medication_analysis_serial_mediation.R
  Fits the serial mediation model linking LEDD, WOQ-9, NMS mood/cognition, and PDQ-39 Summary Index.

data/README.md
  Explains that participant-level datasets are not included in the public repository.
```

## Data access and expected local inputs

The analysis scripts expect the following CSV files in the repository root when run locally:

```text
analysis_patients.csv
analysis_pairs.csv
```

These files are not distributed with the public code repository. After data access is approved by the corresponding author, place the approved de-identified CSV files in the repository root using the filenames above. The separate `analysis_caregivers.csv` file is not required to run the current scripts.

## Analysis steps

Run the scripts from the repository root in this order:

```r
source("script/Step1_NetworkAnalysis.R")
source("script/Step2_Advanced_Analysis.R")
source("script/Step3_ledd_medication_analysis_serial_mediation.R")
```

Step 1 writes intermediate objects under `results_network_analysis/data_processed/`, which are used by Step 2. Figures and summary tables are written under `results_network_analysis/`.

## Software

Analyses were conducted in R 4.5.1. The scripts use:

```text
tidyverse
qgraph
bootnet
networktools
NetworkComparisonTest
lavaan
boot
car
ggplot2
showtext
sysfonts
patchwork
```

## Data availability statement

The analysis code is publicly available at GitHub (https://github.com/legendleaf/pd-mood-cognition-network-dyadic-mediation). The datasets generated and/or analysed during the current study are not publicly available due to participant privacy considerations, but are available from the corresponding author on reasonable request.
