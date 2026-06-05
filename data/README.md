# Data

Participant-level datasets are not included in this public repository because of participant privacy considerations.

After data access is approved by the corresponding author, place the approved de-identified CSV files in the repository root using the filenames expected by the analysis scripts:

```text
analysis_patients.csv     # 165 patients x 51 columns
analysis_caregivers.csv   # 53 caregivers x 6 columns
analysis_pairs.csv        # 53 dyads, patient + caregiver columns joined
```

These local data files are ignored by `.gitignore` and should not be committed.

Each participant is identified by a sequential numeric `ID`. The patient file contains demographics, socioeconomic indicators, lifestyle variables, disease duration, LEDD, EQ-5D, the eight PDQ-39 domains plus the Summary Index, the eight NMSQuest domains plus total, the WOQ-9 total score, and the AD8. The caregiver file contains the ZBI total and the three DASS-21 subscales. The dyadic file is the patient and caregiver tables joined on `ID` for the 53 patients with a returning caregiver.
