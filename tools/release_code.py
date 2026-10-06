"""Export a code-only copy for anonymous review."""
import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = (
    ".gitignore", "README.md", "data/README.md",
    "script/Step0_DataPreparation.R", "script/Step1_NetworkAnalysis.R",
    "script/Step2_Advanced_Analysis.R", "script/Step3_ledd_medication_analysis_serial_mediation.R",
    "script/Step4_PatientSensitivity.R", "script/Step5_TableChecks.R",
    "script/Step6_NetworkSensitivity.R", "script/ItemReliability.py",
    "tools/release_code.py", "tests/check_pipeline.R", "tests/test_release.py",
)
PATTERNS = {
    "email address": r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}",
    "absolute Windows path": r"\b[A-Za-z]:[\\/](?:Users|Research|TempData|PythonScripts)[\\/]",
    "home directory": r"/(?:Users|home)/[A-Za-z0-9_.-]+/",
    "repository account URL": r"https?://(?:www\.)?github\.com/[^\s/]+/[^\s)]+",
    "access token": r"(?:ghp_|github_pat_)[A-Za-z0-9_]{20,}|\bsk-[A-Za-z0-9_-]{20,}",
}


def findings(text, denied=()):
    result = [label for label, pattern in PATTERNS.items() if re.search(pattern, text, re.I)]
    result.extend("author-specified private term" for term in denied if term and term.casefold() in text.casefold())
    return result


def payloads(root=ROOT, denied=()):
    entries = {}
    for relative in FILES:
        path = root / relative
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Missing or linked release file: {relative}")
        if root.resolve() not in path.resolve().parents:
            raise ValueError("Release file resolves outside the repository.")
        text = path.read_text(encoding="utf-8-sig")
        problems = findings(text, denied)
        if problems:
            raise ValueError(f"Review {relative}: {', '.join(problems)}")
        entries[relative] = text.replace("\r\n", "\n").encode("utf-8")
    manifest = {name: dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
                for name, data in entries.items()}
    entries["MANIFEST.json"] = json.dumps(manifest, indent=2, sort_keys=True).encode("utf-8")
    return entries


def export(output, denied=(), root=ROOT):
    entries = payloads(root, denied)
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError("Output already exists; choose a new filename.")
    if output.suffix.lower() != ".zip":
        raise ValueError("Output must be a ZIP file.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in entries.items():
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED)
    return len(entries)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--deny", action="append", default=[], help="Additional private name or identifier to exclude")
    args = parser.parse_args()
    count = export(args.output, args.deny)
    print(f"Exported {count} code/document entries. Git history, data and results were excluded.")
