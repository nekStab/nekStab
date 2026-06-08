"""Catalog model regression tests.

Asserts the contract callers of `validation.catalog` rely on:

- `py_compile` succeeds;
- module imports without side effects;
- `FLOW_FAMILIES`, `CASES`, and `CAPABILITY_MATRIX` are non-empty;
- the cylinder demonstration ladder is fully represented;
- `current_path` and `proposed_folder` are stored separately for every
  case (no filesystem-rename required to use the catalog);
- post-mesh-unification statuses on cylinder UPO/Floquet/OTD cases reflect
  `NEEDS_DNS_SEED` (the IC files were archived 2026-05-21);
- the four public lookup functions resolve.

Runs as `python3 validation/tests/test_catalog.py` (no pytest dep) or
via pytest.
"""
from __future__ import annotations

import py_compile
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CATALOG_PY = REPO_ROOT / "validation" / "catalog.py"
CATALOG_CHECK_PY = REPO_ROOT / "validation" / "catalog_check.py"

# Ensure the validation/ dir is importable as a package-of-one.
sys.path.insert(0, str(REPO_ROOT))

from validation import catalog  # noqa: E402


REQUIRED_CYLINDER_CASES = (
    "cylinder/dns",
    "cylinder/ci_test",
    "cylinder/baseflow/sfd",
    "cylinder/baseflow/sfd_dyn",
    "cylinder/baseflow/sfd_dyn_oifs",
    "cylinder/baseflow/boostconv",
    "cylinder/baseflow/newton",
    "cylinder/baseflow/newton_dyn",
    "cylinder/baseflow/newton_dyn_temp",
    "cylinder/baseflow/newton_upo",
    "cylinder/stability/direct",
    "cylinder/stability/direct_Floquet",
    "cylinder/stability/adjoint",
    "cylinder/stability/adjoint_Floquet",
    "cylinder/stability/animate_modes",
    "cylinder/stability/animate_modes_with_UPO",
    "cylinder/postproc/sensitivity_budget_wavemaker",
    "cylinder/postproc/steady_force_sensitivity",
    "cylinder/otd",
    "cylinder/modal",
    "cylinder/RANS",
)

POST_MESH_UNIFICATION_NEEDS_SEED = (
    "cylinder/baseflow/newton_dyn_temp",
    "cylinder/baseflow/newton_upo",
    "cylinder/stability/direct_Floquet",
    "cylinder/stability/adjoint_Floquet",
    "cylinder/stability/animate_modes_with_UPO",
    "cylinder/otd",
)


def test_catalog_compiles() -> None:
    py_compile.compile(str(CATALOG_PY), doraise=True)


def test_nonempty() -> None:
    assert len(catalog.FLOW_FAMILIES) >= 10, (
        f"FLOW_FAMILIES has only {len(catalog.FLOW_FAMILIES)}"
    )
    assert len(catalog.CASES) >= 25, (
        f"CASES has only {len(catalog.CASES)}"
    )
    assert len(catalog.CAPABILITY_MATRIX) >= 50, (
        f"CAPABILITY_MATRIX has only {len(catalog.CAPABILITY_MATRIX)}"
    )


def test_cylinder_ladder_complete() -> None:
    missing = [c for c in REQUIRED_CYLINDER_CASES if c not in catalog.CASES]
    assert not missing, f"missing cylinder cases: {missing}"


def test_paths_separated() -> None:
    """current_path and proposed_folder are stored as separate fields.

    Before the 2026-05-21 cylinder migration, they had to differ — the
    catalog distinguished "today's path" from "future flat-layout path".
    After the migration, many cases have current_path == proposed_folder
    because today's path *is* the flat-layout target. That's the desired
    end-state. We only require both fields are non-empty (and
    current_path resolves, which catalog_check.py enforces separately).
    """
    for case_id, case in catalog.CASES.items():
        assert case.current_path is not None, f"{case_id}: current_path is None"
        assert case.proposed_folder, f"{case_id}: proposed_folder empty"


def test_post_unification_status() -> None:
    for case_id in POST_MESH_UNIFICATION_NEEDS_SEED:
        case = catalog.get_case(case_id)
        assert case.status == catalog.CaseStatus.NEEDS_DNS_SEED, (
            f"{case_id}: expected NEEDS_DNS_SEED, got {case.status.value}"
        )
        assert case.blockers, (
            f"{case_id}: needs a blocker pointing at the archived 1996-mesh seed"
        )


def test_thermal_target_re30() -> None:
    case = catalog.get_case("cylinder/baseflow/newton_dyn_temp")
    assert case.target_parameters.Re == 30.0, (
        f"thermal target Re is {case.target_parameters.Re}, expected 30"
    )
    assert case.seed_parameters.seed_path == "BFre25t_1cyl0.f00001"


def test_capability_lookups() -> None:
    assert (catalog.get_capability("cylinder_re100", "baseflow")
            == catalog.CapabilityState.IMPLEMENTED)
    assert (catalog.get_capability("slot_FST", "dns")
            == catalog.CapabilityState.IMPLEMENTED)
    assert (catalog.get_capability("nonexistent_family", "dns")
            == catalog.CapabilityState.NOT_APPLICABLE)


def test_public_api() -> None:
    families = catalog.get_flow_families()
    assert families and all(
        isinstance(f, catalog.FlowFamily) for f in families
    )
    all_cases = catalog.get_cases()
    assert all_cases and all(
        isinstance(c, catalog.MethodCase) for c in all_cases
    )
    cyl_cases = catalog.get_cases("cylinder_re100")
    assert len(cyl_cases) >= 18, (
        f"cylinder_re100 filter returned only {len(cyl_cases)}"
    )
    case = catalog.get_case("cylinder/dns")
    assert case.case_id == "cylinder/dns"


def test_catalog_check_lint_passes() -> None:
    """Run the standalone catalog lint as an end-to-end gate."""
    assert CATALOG_CHECK_PY.is_file(), f"missing: {CATALOG_CHECK_PY}"
    result = subprocess.run(
        [sys.executable, str(CATALOG_CHECK_PY)],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, (
        f"catalog_check exited {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert "catalog lint passed" in result.stdout


def test_settings_lint_catches_deliberate_breakage() -> None:
    """Prove the settings/capability checks fail loudly
    when given broken inputs (not via mutating the real catalog).
    """
    from validation import catalog_check
    families = list(catalog.FLOW_FAMILIES.values())

    # 1. static-cylinder non-RANS at Re=200 without a blocker -> must fail
    broken_re = catalog.MethodCase(
        case_id="cylinder/test_bad_re",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder/baseflow/newton",
        proposed_folder="999_test_bad_re", sort_prefix="999",
        label="test", mode_name="test", legacy_uparam01="2.0",
        expected_behavior="should be Re=100",
        status=catalog.CaseStatus.VALIDATED,
        compute_class=catalog.ComputeClass.BOUNDED_LOCAL_VALIDATION,
        target_parameters=catalog.TargetParameters(Re=200.0, model="DNS"),
    )
    errs = catalog_check.lint(families, [broken_re], {}, "1.0")
    assert any("deviates from canonical Re=100" in e for e in errs), errs

    # 2. RANS adjoint -> must fail
    broken_rans = catalog.MethodCase(
        case_id="cylinder/test_bad_rans_adj",
        flow_family="cylinder_re100", method_lane="adjoint",
        current_path="example/cylinder/baseflow/newton",
        proposed_folder="999_test_bad_rans_adj", sort_prefix="999",
        label="test", mode_name="test", legacy_uparam01="3.2",
        expected_behavior="adjoint on RANS — not allowed",
        status=catalog.CaseStatus.VALIDATED,
        compute_class=catalog.ComputeClass.BOUNDED_LOCAL_VALIDATION,
        target_parameters=catalog.TargetParameters(Re=100.0, model="RANS"),
    )
    errs = catalog_check.lint(families, [broken_rans], {}, "1.0")
    assert any("no adjoint formulation" in e for e in errs), errs

    # 3. seed_Re != target_Re without seed_status -> must fail
    broken_seed = catalog.MethodCase(
        case_id="cylinder/test_bad_seed",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder/baseflow/newton",
        proposed_folder="999_test_bad_seed", sort_prefix="999",
        label="test", mode_name="test", legacy_uparam01="2.0",
        expected_behavior="seed Re differs but no seed_status",
        status=catalog.CaseStatus.VALIDATED,
        compute_class=catalog.ComputeClass.BOUNDED_LOCAL_VALIDATION,
        target_parameters=catalog.TargetParameters(Re=100.0, model="DNS"),
        seed_parameters=catalog.SeedParameters(seed_Re=50.0),
    )
    errs = catalog_check.lint(families, [broken_seed], {}, "1.0")
    assert any("seed_parameters.seed_status is empty" in e for e in errs), errs


def test_artifact_paths_under_validation_figures() -> None:
    saw_validation_figure = False
    for case in catalog.CASES.values():
        for art in case.artifacts:
            if "validation/figures/" in art.path:
                saw_validation_figure = True
                break
    assert saw_validation_figure, (
        "expected at least one artifact pointing at validation/figures/"
    )


def main() -> int:
    failures: list[str] = []
    for fn in (
        test_catalog_compiles,
        test_nonempty,
        test_cylinder_ladder_complete,
        test_paths_separated,
        test_post_unification_status,
        test_thermal_target_re30,
        test_capability_lookups,
        test_public_api,
        test_catalog_check_lint_passes,
        test_settings_lint_catches_deliberate_breakage,
        test_artifact_paths_under_validation_figures,
    ):
        try:
            fn()
            print(f"  ok  {fn.__name__}")
        except AssertionError as exc:
            print(f"  FAIL {fn.__name__}: {exc}")
            failures.append(fn.__name__)
        except Exception as exc:
            print(f"  ERROR {fn.__name__}: {exc}")
            failures.append(fn.__name__)
    if failures:
        print(f"\n{len(failures)} test(s) failed: {failures}")
        return 1
    print("\nall catalog checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
