"""Tests for scripts/evidence_report.py — collect_findings() function."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

# Make repo root importable
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from validation.catalog import (
    Artifact,
    CaseStatus,
    ComputeClass,
    EvidenceRequirement,
    EvidenceRole,
    MethodCase,
    SeedParameters,
    TargetParameters,
)

EVIDENCE_REPORT = REPO_ROOT / "scripts" / "evidence_report.py"


def _make_case(
    *,
    case_id: str = "test/case",
    flow_family: str = "cylinder_re100",
    status: CaseStatus = CaseStatus.LOCAL_SMOKE,
    evidence_requirements: tuple = (),
    artifacts: tuple = (),
    evidence_prereqs: tuple = (),
) -> MethodCase:
    """Build a minimal MethodCase for testing."""
    return MethodCase(
        case_id=case_id,
        flow_family=flow_family,
        method_lane="dns",
        current_path="example/cylinder/dns",
        proposed_folder="example/cylinder_re100/000_dns",
        sort_prefix="000",
        label="Test Case",
        mode_name=None,
        legacy_uparam01=None,
        expected_behavior="test",
        status=status,
        compute_class=ComputeClass.QUICK_LOCAL_SMOKE,
        target_parameters=TargetParameters(Re=100.0),
        seed_parameters=SeedParameters(),
        evidence_requirements=evidence_requirements,
        artifacts=artifacts,
        evidence_prereqs=evidence_prereqs,
    )


# ---------------------------------------------------------------------------
# Import collect_findings from the script
# ---------------------------------------------------------------------------
from evidence_report import collect_findings


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_missing_required_no_artifact(tmp_path: Path) -> None:
    """Required EvidenceRequirement with no matching Artifact -> MISSING."""
    req = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE, required=True)
    case = _make_case(evidence_requirements=(req,), artifacts=())
    findings = collect_findings([case], catalog_root=tmp_path)
    assert len(findings) == 1
    assert findings[0]["issue"] == "MISSING"
    assert findings[0]["role"] == EvidenceRole.REFERENCE_IMAGE.value
    assert findings[0]["artifact_path"] is None


def test_contradictory_validated_file_missing(tmp_path: Path) -> None:
    """VALIDATED case with required artifact declared but file not on disk -> CONTRADICTORY."""
    req = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE, required=True)
    art = Artifact(
        role=EvidenceRole.REFERENCE_IMAGE,
        path="validation/figures/no_such_file.png",
        freshness="current",
    )
    case = _make_case(
        status=CaseStatus.VALIDATED,
        evidence_requirements=(req,),
        artifacts=(art,),
    )
    findings = collect_findings([case], catalog_root=tmp_path)
    assert len(findings) == 1
    assert findings[0]["issue"] == "CONTRADICTORY"


def test_missing_non_validated_file_absent(tmp_path: Path) -> None:
    """Non-VALIDATED case with required artifact declared but file missing -> MISSING."""
    req = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE, required=True)
    art = Artifact(
        role=EvidenceRole.REFERENCE_IMAGE,
        path="validation/figures/no_such_file.png",
        freshness="current",
    )
    case = _make_case(
        status=CaseStatus.LOCAL_SMOKE,
        evidence_requirements=(req,),
        artifacts=(art,),
    )
    findings = collect_findings([case], catalog_root=tmp_path)
    assert len(findings) == 1
    assert findings[0]["issue"] == "MISSING"


def test_stale_by_freshness(tmp_path: Path) -> None:
    """Artifact with freshness='stale' and file on disk -> STALE."""
    # Create the artifact file
    rel_path = "validation/figures/stale_file.png"
    full_path = tmp_path / rel_path
    full_path.parent.mkdir(parents=True, exist_ok=True)
    full_path.write_bytes(b"fake png")

    req = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE, required=True)
    art = Artifact(role=EvidenceRole.REFERENCE_IMAGE, path=rel_path, freshness="stale")
    case = _make_case(
        status=CaseStatus.LOCAL_SMOKE,
        evidence_requirements=(req,),
        artifacts=(art,),
    )
    findings = collect_findings([case], catalog_root=tmp_path)
    assert len(findings) == 1
    assert findings[0]["issue"] == "STALE"


def test_stale_by_age(tmp_path: Path) -> None:
    """Artifact with freshness='current' but mtime 200 days ago -> STALE."""
    rel_path = "validation/figures/old_file.png"
    full_path = tmp_path / rel_path
    full_path.parent.mkdir(parents=True, exist_ok=True)
    full_path.write_bytes(b"fake png")

    # Set mtime to 200 days ago
    old_time = time.time() - 200 * 86400
    os.utime(full_path, (old_time, old_time))

    req = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE, required=True)
    art = Artifact(role=EvidenceRole.REFERENCE_IMAGE, path=rel_path, freshness="current")
    case = _make_case(
        status=CaseStatus.LOCAL_SMOKE,
        evidence_requirements=(req,),
        artifacts=(art,),
    )
    findings = collect_findings([case], catalog_root=tmp_path, stale_days=90)
    assert len(findings) == 1
    assert findings[0]["issue"] == "STALE"


def test_clean_no_findings(tmp_path: Path) -> None:
    """Recent file with freshness='current' -> no findings."""
    rel_path = "validation/figures/fresh_file.png"
    full_path = tmp_path / rel_path
    full_path.parent.mkdir(parents=True, exist_ok=True)
    full_path.write_bytes(b"fake png")
    # mtime is now (default), well within 90 days

    req = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE, required=True)
    art = Artifact(role=EvidenceRole.REFERENCE_IMAGE, path=rel_path, freshness="current")
    case = _make_case(
        status=CaseStatus.LOCAL_SMOKE,
        evidence_requirements=(req,),
        artifacts=(art,),
    )
    findings = collect_findings([case], catalog_root=tmp_path, stale_days=90)
    assert findings == []


def test_optional_requirement_missing_no_finding(tmp_path: Path) -> None:
    """Optional requirement (required=False) with no artifact -> no finding."""
    req = EvidenceRequirement(role=EvidenceRole.ANIMATION, required=False)
    case = _make_case(evidence_requirements=(req,), artifacts=())
    findings = collect_findings([case], catalog_root=tmp_path)
    assert findings == []


def test_json_output(tmp_path: Path) -> None:
    """--json flag produces valid JSON list."""
    result = subprocess.run(
        [sys.executable, str(EVIDENCE_REPORT), "--json", "--catalog-root", str(REPO_ROOT)],
        capture_output=True,
        text=True,
    )
    assert result.returncode in (0, 1), f"stderr={result.stderr!r}"
    data = json.loads(result.stdout)
    assert isinstance(data, list)
    # Each element must have the required keys
    for item in data:
        for key in ("case_id", "flow_family", "role", "issue", "artifact_path"):
            assert key in item, f"Missing key {key!r} in {item}"


def test_evidence_prereq_contradictory(tmp_path: Path) -> None:
    """VALIDATED case with evidence_prereq glob matching no files -> CONTRADICTORY."""
    case = _make_case(
        status=CaseStatus.VALIDATED,
        evidence_requirements=(),
        artifacts=(),
        evidence_prereqs=("nonexistent/glob/**/*.png",),
    )
    findings = collect_findings([case], catalog_root=tmp_path)
    assert any(f["issue"] == "CONTRADICTORY" for f in findings)
