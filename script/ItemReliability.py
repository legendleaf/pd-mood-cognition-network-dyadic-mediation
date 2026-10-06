"""Item-level internal consistency for Supplementary Table S3."""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

PDQ = {
    "mobility": range(1, 11), "adl": range(11, 17), "emotional": range(17, 23),
    "stigma": range(23, 27), "social": range(27, 30), "cognition": range(30, 34),
    "communication": range(34, 37), "bodily": range(37, 40),
}
NMS = {
    "gi": range(1, 8), "urinary": [8, 9], "sexual": [18, 19],
    "cardiovascular": [20, 21], "sleep": range(22, 27), "mood": [12, 13, 15, 16, 17],
    "perception": [14, 30], "misc": [10, 11, 27, 28, 29],
}
DASS = {"depression": [3, 5, 10, 13, 16, 17, 21], "anxiety": [2, 4, 7, 9, 15, 19, 20],
        "stress": [1, 6, 8, 11, 12, 14, 18]}


def workbook(path, ids):
    book = load_workbook(path, read_only=True, data_only=True)
    try:
        rows = list(book.worksheets[0].values)
    finally:
        book.close()
    frame = pd.DataFrame(rows[1:], columns=rows[0]).set_index("ID")
    if not frame.index.is_unique or frame.index.hasnans:
        raise ValueError("Workbook IDs must be unique and nonmissing.")
    return frame.loc[ids]


def alpha(x):
    k = x.shape[1]
    variance = x.sum(axis=1).var(ddof=1)
    return k/(k-1)*(1 - x.var(axis=0, ddof=1).sum()/variance) if variance > 0 else np.nan


def reliability(items, rng, reps):
    x = items.to_numpy(dtype=float)
    if not np.isfinite(x).all():
        raise ValueError("Retained item matrix contains missing or nonfinite entries.")
    draws = np.asarray([alpha(x[rng.integers(0, len(x), len(x))]) for _ in range(reps)])
    if not np.isfinite(draws).all():
        raise ValueError("Nonfinite alpha bootstrap draws require review.")
    low, high = np.quantile(draws, [.025, .975], method="linear")
    return dict(n=len(x), items=x.shape[1], alpha=alpha(x), lower=low, upper=high,
                bootstrap_requested=reps, bootstrap_valid=len(draws))


def main():
    data_dir = Path(os.environ.get("PD_DATA_DIR", "private")).resolve(strict=True)
    output = Path(os.environ.get("PD_OUTPUT_DIR", "results_network_analysis")).resolve()
    patient_book = Path(os.environ["PD_ITEM_WORKBOOK"]).resolve(strict=True)
    caregiver_book = Path(os.environ["PD_CAREGIVER_WORKBOOK"]).resolve(strict=True)
    if output == data_dir or data_dir in output.parents:
        raise ValueError("Output directory must be outside the input directory.")
    for book in (patient_book, caregiver_book):
        if output == book.parent or book.parent in output.parents:
            raise ValueError("Workbooks and outputs must use separate directories.")
    patient_csv = data_dir / "analysis_patients.csv"
    caregiver_csv = data_dir / "analysis_caregivers.csv"
    paths = [patient_book, caregiver_book, patient_csv, caregiver_csv]
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    patients = pd.read_csv(patient_csv).set_index("ID")
    carers = pd.read_csv(caregiver_csv).set_index("ID")
    if not patients.index.is_unique or not carers.index.is_unique:
        raise ValueError("Analysis IDs must be unique.")
    rawp, rawc = workbook(patient_book, patients.index), workbook(caregiver_book, carers.index)
    reps = int(os.environ.get("PD_RELIABILITY_BOOTSTRAP", "2000"))
    if reps < 2:
        raise ValueError("At least two bootstrap draws are required.")
    rng = np.random.default_rng(20260912)
    results = []
    for name, indices in PDQ.items():
        items = rawp[[f"PDQ39_{i}" for i in indices]].astype(float)
        if not items.isin(range(1, 6)).all().all():
            raise ValueError("PDQ items must use codes 1-5.")
        results.append(dict(scale=f"PDQ39_{name}", cohort="patients", **reliability(items, rng, reps)))
    for name, indices in NMS.items():
        items = rawp[[f"\u975e\u8fd0\u52a8_{i}" for i in indices]]
        if not items.isin([1, 2]).all().all():
            raise ValueError("NMS items must use codes 1/2.")
        binary = (items == 1).astype(float)
        if not np.allclose(binary.sum(axis=1), patients[f"NMS_{name}"]):
            raise ValueError("NMS domain does not match retained counts.")
        results.append(dict(scale=f"NMS_{name}", cohort="patients", **reliability(binary, rng, reps)))
    items = rawc[[f"ZBI_{i}" for i in range(1, 23)]].astype(float) - 1
    if not items.isin(range(5)).all().all() or not np.allclose(items.sum(axis=1), carers.ZBI_score):
        raise ValueError("ZBI items do not match retained totals.")
    results.append(dict(scale="ZBI_score", cohort="caregivers", **reliability(items, rng, reps)))
    for name, indices in DASS.items():
        items = rawc[[f"DASS_{i}" for i in indices]].astype(float) - 1
        if not items.isin(range(4)).all().all() or not np.allclose(items.sum(axis=1)*2, carers[f"DASS_{name}"]):
            raise ValueError("DASS items do not match retained scores.")
        results.append(dict(scale=f"DASS_{name}", cohort="caregivers", **reliability(items, rng, reps)))
    woq = rawp[[f"WOQ9_{i}" for i in range(1, 10)]]
    if not woq.isin(range(4)).all().all():
        raise ValueError("WOQ items must use the documented 0/1/2/3 categories.")
    positive = (woq == 1).astype(float)
    if not np.array_equal(positive.sum(axis=1).to_numpy(), patients.WOQ9_total.to_numpy()):
        raise ValueError("WOQ positive-item count does not match retained totals.")
    results.append(dict(scale="WOQ9_total", cohort="patients",
                        **reliability(positive, np.random.default_rng(20260930), reps)))
    if any(hashlib.sha256(p.read_bytes()).hexdigest() != hashes[p] for p in paths):
        raise RuntimeError("An input file changed during execution.")
    (output / "tables").mkdir(parents=True, exist_ok=True)
    (output / "logs").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(output / "tables/item_reliability.csv", index=False)
    (output / "logs/reliability_environment.json").write_text(json.dumps(
        dict(numpy=np.__version__, pandas=pd.__version__, bootstrap=reps,
             domain_seed=20260912, woq_seed=20260930, intervals="percentile, linear quantiles",
             inputs_unchanged=True), indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
