# Data

Participant-level datasets are not included in this public repository because of participant privacy considerations.

After access is approved, keep the de-identified CSV files outside the repository and set `PD_DATA_DIR` to that directory. Expected filenames are:

```text
analysis_patients.csv     # 165 patients; columns depend on the approved export
analysis_caregivers.csv   # 53 caregivers x 6 columns
analysis_pairs.csv        # 53 dyads, patient + caregiver columns joined
```

Do not put participant files in this directory. Column counts differ between the original and R1 exports because derived PDQ variables were added. Required fields and scoring checks are described in the main README and Step0_DataPreparation.R.

IDs are linkage keys; replacing names with IDs does not by itself make a clinical dataset anonymous.

Each participant is identified by a sequential numeric `ID`. The patient file contains demographics, socioeconomic indicators, lifestyle variables, disease duration, LEDD, EQ-5D, the eight PDQ-39 domains plus the Summary Index, the eight NMSQuest domains plus total, the WOQ-9 total score, and the AD8. The caregiver file contains the ZBI total and the three DASS-21 subscales. The dyadic file is the patient and caregiver tables joined on `ID` for the 53 patients with a returning caregiver.
