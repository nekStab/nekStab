import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

# validate.py lives at validation/validate.py; expose the repo root so the
# `validation` package (and src/ for the source-text checks below) resolve from
# any CWD, matching the convention in the sibling validation/tests modules.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from validation import validate  # noqa: E402


class ValidateScriptTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.example_root = self.root / "example"
        self.example_root.mkdir()
        self.cpu = validate.CPUInfo(cpu_list="0", count=1)

        patcher = mock.patch.object(validate, "NEKSTAB_ROOT", self.root)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.addCleanup(self._tmp.cleanup)

    def test_validate_case_fails_if_declared_copy_output_missing(self):
        case_dir = self.example_root / "flip_flop" / "baseflow"
        dest_dir = self.example_root / "flip_flop" / "stability" / "direct_Floquet"
        case_dir.mkdir(parents=True)
        dest_dir.mkdir(parents=True)

        case = {
            "description": "copy regression",
            "dir": "flip_flop/baseflow",
            "casename": "2cyl",
            "requires": [],
            "output": None,
            "copies_to": [
                {"file": "BF_2cyl0.f00001", "dest_dir": "flip_flop/stability/direct_Floquet"}
            ],
        }

        with mock.patch.object(validate, "_compile_case", return_value=True), \
             mock.patch.object(validate, "_run_case", return_value=True), \
             mock.patch.object(validate, "check_logfile_success", return_value=True), \
             mock.patch.object(validate, "get_optimal_nprocs", return_value=1):
            result = validate.validate_case(
                "flipflop_bf", case, check_only=False, nprocs=1, cpu=self.cpu, dry_run=False
            )

        # validate_case returns (status, elapsed_seconds); assert on the status.
        self.assertEqual(result[0], "fail")

    def test_discover_compilable_cases_ignores_hidden_placeholder_usr(self):
        case_dir = self.example_root / "cylinder" / "modal"
        case_dir.mkdir(parents=True)
        (case_dir / "SIZE").write_text("dummy\n")
        (case_dir / ".usr").write_text("")
        (case_dir / "1cyl.usr").write_text("")

        real_glob = Path.glob

        def fake_glob(path_obj, pattern, *args, **kwargs):
            # Forward *args/**kwargs so rglob's internal case_sensitive= call
            # (Python 3.13+) passes through to the real implementation.
            if path_obj == case_dir and pattern == "*.usr":
                return [case_dir / ".usr", case_dir / "1cyl.usr"]
            return real_glob(path_obj, pattern, *args, **kwargs)

        with mock.patch("pathlib.Path.glob", new=fake_glob):
            cases = validate.discover_compilable_cases(self.example_root)

        self.assertEqual(cases, [(case_dir, "1cyl")])

    def test_validate_case_materializes_requires_from_external_data_root(self):
        case_dir = self.example_root / "cylinder" / "dns"
        case_dir.mkdir(parents=True)

        data_root = self.root / "external_data"
        external_file = data_root / "example" / "cylinder" / "dns" / "rst_1cyl0.f00001"
        external_file.parent.mkdir(parents=True)
        external_file.write_text("restart\n")

        case = {
            "description": "external data prereq",
            "dir": "cylinder/dns",
            "casename": "1cyl",
            "requires": ["rst_1cyl0.f00001"],
            "output": None,
        }

        with mock.patch.dict(os.environ, {"NEKSTAB_DATA_ROOT": str(data_root)}), \
             mock.patch.object(validate, "check_logfile_success", return_value=True), \
             mock.patch.object(validate, "get_optimal_nprocs", return_value=1):
            result = validate.validate_case(
                "cyl_dns", case, check_only=True, nprocs=1, cpu=self.cpu, dry_run=False
            )

        # validate_case returns (status, elapsed_seconds); assert on the status.
        self.assertEqual(result[0], "pass")
        self.assertTrue((case_dir / "rst_1cyl0.f00001").exists())

    def test_outpost_ks_uses_per_scalar_extents_for_passive_scalars(self):
        source = (REPO_ROOT / "src" / "eigensolvers.f90").read_text()

        self.assertIn("n_scalar = nx1*ny1*nz1*nelfld(m + 1)", source)
        self.assertIn("oks_fp_ct_s(1:n_scalar, m) = matmul(", source)
        self.assertNotIn("call copy(oks_qt_s(1, m, i), Q(i)%t(1, m), n_temp)", source)
        self.assertNotIn("oks_fp_ct_s(1:n_temp, m) = matmul(", source)


if __name__ == "__main__":
    unittest.main()
