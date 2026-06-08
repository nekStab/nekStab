"""Smoke tests for validation/catalog.py.

Run as: python3 validation/tests/test_catalog_smoke.py
Or:     python3 -m pytest validation/tests/test_catalog_smoke.py -v
"""
import sys
import os

# Allow running as a plain script without installing the package
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from validation import catalog  # noqa: E402


def test_import():
    """catalog imports cleanly with no side effects."""
    assert catalog is not None


def test_families_nonempty():
    """families() returns at least 13 entries (spine contract: 13 geometries)."""
    fams = catalog.families()
    assert len(fams) >= 13, (
        f"Expected >= 13 families, got {len(fams)}: {[f.family_id for f in fams]}"
    )


def test_cylinder_re100_exists():
    """family('cylinder_re100') is registered in the catalog."""
    fam = catalog.family("cylinder_re100")
    assert fam is not None
    assert fam.family_id == "cylinder_re100"


def test_cylinder_re100_has_cases():
    """static_cylinder (the data store for cylinder_re100 cases) has >= 10 method cases."""
    # cylinder_re100 is a public-name alias; actual cases use flow_family='static_cylinder'
    cyl_cases = catalog.cases_for("cylinder_re100")
    assert len(cyl_cases) >= 10, (
        f"Expected >= 10 cylinder cases, got {len(cyl_cases)}: "
        f"{[c.case_id for c in cyl_cases]}"
    )


def test_method_cases_fields():
    """Every MethodCase has nonempty sort_prefix, nonempty label, and non-None compute_class."""
    for case in catalog.all_cases():
        assert case.sort_prefix, f"Empty sort_prefix on {case.case_id}"
        assert case.label, f"Empty label on {case.case_id}"
        assert case.compute_class is not None, f"None compute_class on {case.case_id}"


def test_by_legacy_uparam_sfd():
    """by_legacy_uparam(1.1) returns >= 4 SFD cases (3 existing dirs + 1 planned stub)."""
    sfd_cases = catalog.by_legacy_uparam(1.1)
    assert len(sfd_cases) >= 4, (
        f"Expected >= 4 cases for uparam=1.1, got {len(sfd_cases)}: "
        f"{[c.case_id for c in sfd_cases]}"
    )


def test_family_can_re1m_adjoint_false():
    """family_can('cylinder_re1m', 'adjoint') == False (RANS-only; no adjoint lane)."""
    result = catalog.family_can("cylinder_re1m", "adjoint")
    assert result is False, f"Expected False, got {result!r}"


def test_family_can_cylinder_re100_sfd_true():
    """family_can('cylinder_re100', 'sfd') == True (alias resolves to static_cylinder)."""
    result = catalog.family_can("cylinder_re100", "sfd")
    assert result is True, (
        f"Expected True for cylinder_re100/sfd, got {result!r}. "
        f"Check that CAPABILITY_MATRIX has ('static_cylinder', 'sfd') = implemented."
    )


def test_schema_version():
    """schema_version() returns '0.1.0'."""
    v = catalog.schema_version()
    assert v == "0.1.0", f"Got {v!r}"


if __name__ == "__main__":
    test_import()
    print("PASS: import")
    test_families_nonempty()
    print("PASS: families nonempty (>= 13)")
    test_cylinder_re100_exists()
    print("PASS: cylinder_re100 family exists")
    test_cylinder_re100_has_cases()
    print("PASS: cylinder_re100 cases >= 10")
    test_method_cases_fields()
    print("PASS: method case fields")
    test_by_legacy_uparam_sfd()
    print("PASS: by_legacy_uparam(1.1) >= 4")
    test_family_can_re1m_adjoint_false()
    print("PASS: family_can cylinder_re1m adjoint == False")
    test_family_can_cylinder_re100_sfd_true()
    print("PASS: family_can cylinder_re100 sfd == True")
    test_schema_version()
    print("PASS: schema_version == '0.1.0'")
    print("\nAll smoke tests passed.")
