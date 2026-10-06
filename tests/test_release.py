"""Release checks use synthetic strings and temporary output files."""
import importlib.util
import hashlib
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("release", ROOT / "tools/release_code.py")
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class ReleaseTests(unittest.TestCase):
    def test_common_identifiers(self):
        email = "researcher" + "@" + "example.org"
        private_path = "/".join(["C:", "Users", "Example", "private.csv"])
        self.assertIn("email address", release.findings(email))
        self.assertIn("absolute Windows path", release.findings(private_path))
        self.assertTrue(release.findings("private study label", ["study label"]))

    def test_code_only_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "review.zip"
            release.export(path)
            with zipfile.ZipFile(path) as archive:
                self.assertEqual(set(archive.namelist()), set(release.FILES) | {"MANIFEST.json"})
                self.assertFalse(archive.comment)
                manifest = json.loads(archive.read("MANIFEST.json"))
                self.assertEqual(set(manifest), set(release.FILES))
                for name, record in manifest.items():
                    data = archive.read(name)
                    self.assertEqual(len(data), record["bytes"])
                    self.assertEqual(hashlib.sha256(data).hexdigest(), record["sha256"])
                for info in archive.infolist():
                    self.assertEqual(info.date_time, (1980, 1, 1, 0, 0, 0))
                    self.assertFalse(info.comment)
            with self.assertRaises(FileExistsError):
                release.export(path)


if __name__ == "__main__":
    unittest.main()
