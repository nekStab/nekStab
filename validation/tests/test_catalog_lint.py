"""Tests for validation.catalog_check."""
from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation import catalog, catalog_check  # noqa: E402


def run_lint(cases):
    return catalog_check.lint(
        catalog.families(),
        cases,
        catalog.CAPABILITY_MATRIX,
        catalog.schema_version(),
    )


def test_catalog_passes_lint():
    assert run_lint(catalog.all_cases()) == []


def test_duplicate_prefix_is_reported():
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/dns")
    duplicate = replace(
        base,
        case_id="synthetic/duplicate_prefix",
        proposed_folder="example/cylinder_re100/000_dns",
    )
    errors = run_lint(cases + [duplicate])
    assert errors
    assert any("duplicate" in error.lower() for error in errors)


def test_missing_current_path_is_reported():
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/dns")
    missing = replace(
        base,
        case_id="synthetic/missing_current_path",
        sort_prefix="998",
        current_path="/nonexistent/path/that/does/not/exist",
        proposed_folder="example/cylinder/998_synthetic_missing_current_path",
    )
    errors = run_lint(cases + [missing])
    assert errors
    assert any("current_path" in error or "missing" in error for error in errors)


def test_malformed_proposed_folder_is_reported():
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/dns")
    malformed = replace(
        base,
        case_id="synthetic/malformed_proposed_folder",
        sort_prefix="999",
        proposed_folder="bad_no_prefix",
    )
    errors = run_lint(cases + [malformed])
    assert errors
    assert any("proposed_folder" in error for error in errors)



def test_dangling_prerequisite_is_reported():
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/dns")
    bad_prereq = replace(
        base,
        case_id="synthetic/dangling_prereq",
        sort_prefix="991",
        current_path="",
        proposed_folder="example/cylinder/991_synthetic_dangling_prereq",
        prerequisites=("nonexistent/case/that/does/not/exist",),
    )
    errors = run_lint(cases + [bad_prereq])
    assert errors, "expected at least one error for dangling prerequisite"
    assert any(
        ("prerequisite" in e.lower() or "unknown" in e.lower()) for e in errors
    ), f"expected prerequisite/unknown error, got: {errors}"


def test_no_dangling_prerequisites_in_real_catalog():
    errors = run_lint(catalog.all_cases())
    dangling = [
        e for e in errors
        if "prerequisite" in e.lower() and "unknown" in e.lower()
    ]
    assert not dangling, f"real catalog has dangling prerequisites: {dangling}"


def test_cycle_detection():
    """A synthetic 3-node cycle A->B->C->A must be reported."""
    base = catalog.get_case("cylinder/dns")
    case_a = replace(
        base,
        case_id="synthetic/cycle_a",
        sort_prefix="995",
        current_path="",
        proposed_folder="example/cylinder/995_synthetic_cycle_a",
        prerequisites=("synthetic/cycle_c",),
    )
    case_b = replace(
        base,
        case_id="synthetic/cycle_b",
        sort_prefix="996",
        current_path="",
        proposed_folder="example/cylinder/996_synthetic_cycle_b",
        prerequisites=("synthetic/cycle_a",),
    )
    case_c = replace(
        base,
        case_id="synthetic/cycle_c",
        sort_prefix="997",
        current_path="",
        proposed_folder="example/cylinder/997_synthetic_cycle_c",
        prerequisites=("synthetic/cycle_b",),
    )
    errors = run_lint([case_a, case_b, case_c])
    assert any("cycle" in e.lower() for e in errors), (
        f"expected cycle error, got: {errors}"
    )


# ---------------------------------------------------------------------------
# Numbering consistency and uparam preservation invariants
# ---------------------------------------------------------------------------

def test_sort_prefix_uparam_family_passes_on_real_catalog():
    """All real catalog cases pass the sort_prefix/uparam family rule."""
    errors = run_lint(catalog.all_cases())
    numbering_errors = [e for e in errors if "does not match legacy_uparam01" in e]
    assert not numbering_errors, f"real catalog has uparam family mismatches: {numbering_errors}"


def test_sort_prefix_uparam_family_violation():
    """sort_prefix in the 200-range with uparam01='1.1' (hundreds mismatch) must be reported."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/baseflow/sfd")
    bad = replace(
        base,
        case_id="synthetic/bad_uparam_family",
        # sort_prefix hundreds=2, but legacy_uparam01="1.1" → expected hundreds=1
        sort_prefix="210",
        legacy_uparam01="1.1",
        proposed_folder="example/cylinder/210_synthetic_bad_uparam_family",
    )
    errors = run_lint(cases + [bad])
    assert any("does not match legacy_uparam01" in e for e in errors), (
        f"expected uparam-family mismatch error, got: {errors}"
    )


def test_uparam_not_null_solver_lane_passes_on_real_catalog():
    """All real solver-lane cases have non-null legacy_uparam01."""
    errors = run_lint(catalog.all_cases())
    null_errors = [e for e in errors if "requires a legacy_uparam01" in e]
    assert not null_errors, f"real catalog has null uparam on solver lanes: {null_errors}"


def test_uparam_not_null_solver_lane_violation():
    """A baseflow-lane case with legacy_uparam01=None must be reported."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/baseflow/sfd")
    bad = replace(
        base,
        case_id="synthetic/null_uparam_baseflow",
        sort_prefix="115",
        legacy_uparam01=None,
        proposed_folder="example/cylinder/115_synthetic_null_uparam_baseflow",
    )
    errors = run_lint(cases + [bad])
    assert any("requires a legacy_uparam01" in e for e in errors), (
        f"expected uparam-null error, got: {errors}"
    )


def test_no_duplicate_prefix_variant_passes_on_real_catalog():
    """No real catalog cases share (flow_family, sort_prefix, variant_suffix)."""
    errors = run_lint(catalog.all_cases())
    dup_errors = [e for e in errors if "duplicates" in e and "sort_prefix" in e]
    assert not dup_errors, f"real catalog has duplicate prefix-variant pairs: {dup_errors}"


def test_no_duplicate_prefix_variant_violation():
    """Two cases with the same proposed_folder (same variant suffix) in the same family must be reported."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/baseflow/sfd")
    dup = replace(
        base,
        case_id="synthetic/dup_variant",
        sort_prefix="110",
        # Identical proposed_folder → same (flow_family, sort_prefix, variant_suffix)
        proposed_folder="example/cylinder_re100/110_baseflow_sfd/akervik",
    )
    errors = run_lint(cases + [dup])
    assert any("duplicates" in e and "sort_prefix" in e for e in errors), (
        f"expected duplicate prefix-variant error, got: {errors}"
    )


def test_variant_subdir_allowed_stage_passes_on_real_catalog():
    """No real catalog cases have variant subdirs in disallowed stage ranges."""
    errors = run_lint(catalog.all_cases())
    variant_errors = [e for e in errors if "variant subdir" in e and "outside the allowed spine" in e]
    assert not variant_errors, f"real catalog has variant-subdir violations: {variant_errors}"


def test_variant_subdir_allowed_stage_violation():
    """A proposed_folder with 4 segments under a 600-series prefix must be reported."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/modal")
    bad = replace(
        base,
        case_id="synthetic/bad_variant_subdir",
        sort_prefix="600",
        legacy_uparam01=None,
        # 4-segment path: example / cylinder / 600_cyl_modal / pod_variant
        proposed_folder="example/cylinder/600_cyl_modal/pod_variant",
    )
    errors = run_lint(cases + [bad])
    assert any("variant subdir" in e and "outside the allowed spine" in e for e in errors), (
        f"expected variant-subdir stage error, got: {errors}"
    )



# ---------------------------------------------------------------------------
# ArtifactRole / Artifact.role schema enforcement
# ---------------------------------------------------------------------------

def test_artifact_role_enum_complete():
    """EvidenceRole (== ArtifactRole) must contain the expected member names."""
    required_members = {
        "REFERENCE_IMAGE", "DNS_TIME_HISTORY", "RESIDUAL_HISTORY",
        "SPECTRUM", "EIGENVALUE_TABLE", "MODE_SHAPE",
        "ANIMATION", "CHECKPOINT_STATS", "PROVENANCE", "LOG", "BLOCKING_NOTE",
    }
    actual_members = {m.name for m in catalog.EvidenceRole}
    missing = required_members - actual_members
    assert not missing, f"EvidenceRole is missing members: {missing}"
    # ArtifactRole alias must be importable and identical
    assert catalog.ArtifactRole is catalog.EvidenceRole


def test_real_catalog_artifact_roles_pass_lint():
    """All artifacts in the real catalog carry valid EvidenceRole values."""
    errors = run_lint(catalog.all_cases())
    role_errors = [e for e in errors if "artifact role" in e and "not in EvidenceRole" in e]
    assert not role_errors, f"real catalog has invalid artifact roles: {role_errors}"


def test_artifact_role_violation_fails_lint():
    """An Artifact whose .role is a raw string (not enum) must be caught by lint."""
    import types
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/dns")

    # Build a fake artifact with a sentinel string role to bypass frozen-dataclass typing
    bad_artifact = catalog.Artifact.__new__(catalog.Artifact)
    object.__setattr__(bad_artifact, "role", "not_an_enum_value")
    object.__setattr__(bad_artifact, "path", "fake/path.txt")
    object.__setattr__(bad_artifact, "freshness", None)
    object.__setattr__(bad_artifact, "provenance", None)

    # Inject the bad artifact into a copy of an existing case
    from dataclasses import replace
    bad_case = replace(
        base,
        case_id="synthetic/bad_artifact_role",
        sort_prefix="993",
        current_path="",
        proposed_folder="example/cylinder/993_synthetic_bad_artifact_role",
        artifacts=(bad_artifact,),
    )
    errors = run_lint(cases + [bad_case])
    assert any("artifact role" in e and "not in EvidenceRole" in e for e in errors), (
        f"expected artifact-role enum error, got: {errors}"
    )


# ---------------------------------------------------------------------------
# Evidence diversity checks
# ---------------------------------------------------------------------------

def test_validated_cases_have_figure_and_history():
    """All upgraded VALIDATED cases in the real catalog must have diverse artifacts."""
    errors = run_lint(catalog.all_cases())
    diversity_errors = [e for e in errors if "VALIDATED case has no" in e]
    assert not diversity_errors, f"real catalog has diversity failures: {diversity_errors}"


def test_diversity_check_missing_history_fails():
    """A VALIDATED case with only a FIGURE artifact must fail the diversity check."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/dns")
    bad = replace(
        base,
        case_id="synthetic/no_history",
        sort_prefix="994",
        current_path="",
        proposed_folder="example/cylinder/994_synthetic_no_history",
        artifacts=(catalog.Artifact(role=catalog.EvidenceRole.REFERENCE_IMAGE,
                                     path="validation/figures/cylinder_dns_Re50.png",
                                     freshness="current"),),
        status=catalog.CaseStatus.VALIDATED,
    )
    errors = run_lint(cases + [bad])
    assert any("VALIDATED case has no HISTORY" in e for e in errors), (
        f"expected history-diversity error, got: {errors}"
    )


def test_diversity_check_missing_figure_fails():
    """A VALIDATED case with only a HISTORY artifact must fail the diversity check."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/dns")
    bad = replace(
        base,
        case_id="synthetic/no_figure",
        sort_prefix="992",
        current_path="",
        proposed_folder="example/cylinder/992_synthetic_no_figure",
        artifacts=(catalog.Artifact(role=catalog.EvidenceRole.DNS_TIME_HISTORY,
                                     path="example/cylinder/dns/1cyl.his",
                                     freshness="current"),),
        status=catalog.CaseStatus.VALIDATED,
    )
    errors = run_lint(cases + [bad])
    assert any("VALIDATED case has no FIGURE" in e for e in errors), (
        f"expected figure-diversity error, got: {errors}"
    )


# ---------------------------------------------------------------------------
# Evidence diversity checks
# ---------------------------------------------------------------------------

def test_validated_cases_have_figure_and_history():
    """All non-postproc VALIDATED cases must have FIGURE + HISTORY artifacts."""
    errors = run_lint(catalog.all_cases())
    diversity_errors = [e for e in errors if "VALIDATED case has no" in e]
    assert not diversity_errors, f"real catalog has diversity failures: {diversity_errors}"


def test_diversity_check_missing_history_fails():
    """A VALIDATED baseflow case with only a FIGURE artifact must fail the diversity check."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/baseflow/sfd")
    bad = replace(
        base,
        case_id="synthetic/no_history",
        sort_prefix="994",
        current_path="",
        proposed_folder="example/cylinder/994_synthetic_no_history",
        # Only a figure, no history
        artifacts=(catalog.Artifact(
            role=catalog.EvidenceRole.REFERENCE_IMAGE,
            path="validation/figures/cylinder_dns_Re50.png",
            freshness="current",
        ),),
        status=catalog.CaseStatus.VALIDATED,
    )
    errors = run_lint(cases + [bad])
    assert any("VALIDATED case has no HISTORY" in e for e in errors), (
        f"expected history-diversity error, got: {errors}"
    )


def test_diversity_check_missing_figure_fails():
    """A VALIDATED baseflow case with only a HISTORY artifact must fail the diversity check."""
    cases = list(catalog.all_cases())
    base = catalog.get_case("cylinder/baseflow/sfd")
    bad = replace(
        base,
        case_id="synthetic/no_figure",
        sort_prefix="992",
        current_path="",
        proposed_folder="example/cylinder/992_synthetic_no_figure",
        # Only history, no figure
        artifacts=(catalog.Artifact(
            role=catalog.EvidenceRole.RESIDUAL_HISTORY,
            path="example/cylinder/baseflow/sfd/residu.dat",
            freshness="current",
        ),),
        status=catalog.CaseStatus.VALIDATED,
    )
    errors = run_lint(cases + [bad])
    assert any("VALIDATED case has no FIGURE" in e for e in errors), (
        f"expected figure-diversity error, got: {errors}"
    )


def main() -> None:
    test_catalog_passes_lint()
    test_duplicate_prefix_is_reported()
    test_missing_current_path_is_reported()
    test_malformed_proposed_folder_is_reported()
    print("catalog lint tests passed")


if __name__ == "__main__":
    main()
