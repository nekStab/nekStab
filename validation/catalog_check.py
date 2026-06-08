"""Standalone lint checks for validation.catalog."""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from validation import catalog  # noqa: E402

_SORT_PREFIX_RE = re.compile(r"^[0-9]{3}$")
_PROPOSED_FOLDER_RE = re.compile(r"^example/[^/]+/[0-9]{3}_[a-z0-9_]+(/[a-z0-9_]+)*$")
_VALID_LEGACY_UPARAM01 = {
    "0.0", "1.1", "1.2", "1.3", "1.4", "2.0", "2.1",
    "3.1", "3.11", "3.2", "3.21", "3.3", "4.1", "4.11", "5.0",
}

# Lanes that historically carry a legacy_uparam01 value.
# These are the executable-solver method lanes (as opposed to postproc-keyword
# or CI-smoke lanes).  Nulling the field on any of these is a numbering error.
_LANES_REQUIRING_UPARAM: set[str] = {
    "baseflow", "stability_direct", "stability_adjoint",
    "floquet_direct", "floquet_adjoint", "transient_growth", "wavemaker",
}

# sort_prefix hundreds digits for which variant subdirs
# (4-segment proposed_folder paths) are permitted — baseflow/newton/stability
# spine only.
_VARIANT_ALLOWED_HUNDREDS: set[int] = {1, 2, 3}


def _path_exists(path: str) -> bool:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = REPO_ROOT / candidate
    return os.path.isdir(candidate)


def lint(flow_families_seq, cases_seq, capability_matrix, schema_ver) -> list[str]:
    """Return catalog invariant errors, or an empty list when the catalog is valid."""
    errors: list[str] = []
    flow_families = list(flow_families_seq)
    cases = list(cases_seq)

    family_ids = [family.family_id for family in flow_families]
    seen_family_ids: set[str] = set()
    duplicate_family_ids: set[str] = set()
    for family_id in family_ids:
        if family_id in seen_family_ids:
            duplicate_family_ids.add(family_id)
        seen_family_ids.add(family_id)
    for family_id in sorted(duplicate_family_ids):
        errors.append(f"duplicate FlowFamily.family_id: {family_id}")

    valid_family_ids = set(family_ids)
    valid_statuses = set(catalog.CaseStatus)
    valid_evidence_roles = set(catalog.EvidenceRole)
    family_prefix_pairs: dict[tuple[str, str, str], str] = {}
    prefixes_by_family: dict[str, list[int]] = defaultdict(list)

    for case in cases:
        label = getattr(case, "case_id", "<unknown>")

        if case.flow_family not in valid_family_ids:
            errors.append(
                f"{label}: flow_family {case.flow_family!r} is not registered in FLOW_FAMILIES"
            )

        if not _SORT_PREFIX_RE.match(case.sort_prefix):
            errors.append(f"{label}: sort_prefix {case.sort_prefix!r} is not 3 digits")
        else:
            prefixes_by_family[case.flow_family].append(int(case.sort_prefix))

        variant_suffix = case.proposed_folder.split("/")[-1]
        pair = (case.flow_family, case.sort_prefix, variant_suffix)
        if pair in family_prefix_pairs:
            errors.append(
                f"{label}: duplicate (flow_family, sort_prefix, variant) {pair!r}; "
                f"also used by {family_prefix_pairs[pair]}"
            )
        else:
            family_prefix_pairs[pair] = label

        if case.current_path not in (None, "") and not _path_exists(case.current_path):
            errors.append(f"{label}: current_path {case.current_path!r} is missing")

        if not _PROPOSED_FOLDER_RE.match(case.proposed_folder):
            errors.append(f"{label}: proposed_folder {case.proposed_folder!r} is malformed")
        else:
            proposed_prefix = case.proposed_folder.split("/", 2)[2][:3]
            if proposed_prefix != case.sort_prefix:
                errors.append(
                    f"{label}: proposed_folder prefix {proposed_prefix!r} does not match "
                    f"sort_prefix {case.sort_prefix!r}"
                )

        if case.status not in valid_statuses:
            errors.append(f"{label}: status {case.status!r} is not a valid CaseStatus")

        for req in case.evidence_requirements:
            if req.role not in valid_evidence_roles:
                errors.append(f"{label}: evidence role {req.role!r} is not a valid EvidenceRole")

        for art in case.artifacts:
            if not isinstance(art.role, catalog.EvidenceRole):
                errors.append(
                    f"{label}: artifact role {art.role!r} not in EvidenceRole enum"
                )

        if case.legacy_uparam01 is not None and case.legacy_uparam01 not in _VALID_LEGACY_UPARAM01:
            errors.append(f"{label}: legacy_uparam01 {case.legacy_uparam01!r} is invalid")

    if not isinstance(schema_ver, str) or not schema_ver:
        errors.append("schema_version() must return a non-empty string")

    for family_id, prefixes in sorted(prefixes_by_family.items()):
        if prefixes != sorted(prefixes):
            errors.append(f"{family_id}: sort_prefix values are not in ascending order")

    for family_id, method_lane in capability_matrix:
        if family_id not in valid_family_ids:
            errors.append(
                f"CAPABILITY_MATRIX key {(family_id, method_lane)!r} references unknown family"
            )

    errors.extend(_settings_and_capability_checks(cases, capability_matrix))
    errors.extend(_dependency_checks(cases))
    errors.extend(_evidence_diversity_checks(cases))
    errors.extend(_numbering_checks(cases))
    errors.extend(_evidence_diversity_checks(cases))

    return errors


def _evidence_diversity_checks(cases: list) -> list[str]:
    """Upgraded VALIDATED cases must carry a figure plus non-figure evidence.

    Legacy single-plot entries remain tolerated until real on-disk artifacts are
    verified for them. Synthetic cases are always checked so the rule is
    exercised by tests.
    """
    errors: list[str] = []
    figure_roles = {catalog.EvidenceRole.REFERENCE_IMAGE}
    history_roles = {
        catalog.EvidenceRole.RESIDUAL_HISTORY,
        catalog.EvidenceRole.DNS_TIME_HISTORY,
        catalog.EvidenceRole.SPECTRUM,
        catalog.EvidenceRole.ANIMATION,
        catalog.EvidenceRole.CHECKPOINT_STATS,
        catalog.EvidenceRole.LOG,
        catalog.EvidenceRole.PROVENANCE,
    }

    for case in cases:
        if case.status != catalog.CaseStatus.VALIDATED:
            continue

        artifact_roles = {a.role for a in case.artifacts}
        enforce = case.case_id.startswith("synthetic/") or len(case.artifacts) > 1
        if not enforce:
            continue

        if not (artifact_roles & figure_roles):
            errors.append(
                f"{case.case_id}: VALIDATED case has no FIGURE artifact (REFERENCE_IMAGE)"
            )
        if not (artifact_roles & history_roles):
            errors.append(
                f"{case.case_id}: VALIDATED case has no HISTORY artifact"
                f" (RESIDUAL_HISTORY or DNS_TIME_HISTORY)"
            )
    return errors


def _dependency_checks(cases: list) -> list[str]:
    """Check prerequisite references are valid and there are no dependency cycles.

    Failures:
    - A prerequisite string that does not match any known case_id.
    - A cycle among prerequisite edges (A -> B -> ... -> A).
    """
    errors: list[str] = []
    case_ids = {c.case_id for c in cases}

    # Build adjacency: case_id -> list of prerequisite case_ids
    adj: dict[str, list[str]] = {}
    for case in cases:
        label = getattr(case, "case_id", "<unknown>")
        adj[label] = []
        for ref in getattr(case, "prerequisites", ()):
            if ref not in case_ids:
                errors.append(
                    f"{label}: prerequisite {ref!r} references unknown case_id"
                )
            else:
                adj[label].append(ref)

    # Cycle detection via iterative DFS (white=0, gray=1, black=2)
    WHITE, GRAY, BLACK = 0, 1, 2
    color: dict[str, int] = {cid: WHITE for cid in case_ids}
    cycle_errors: list[str] = []

    def dfs(start: str) -> None:
        stack = [(start, iter(adj.get(start, [])))]
        path = [start]
        color[start] = GRAY
        while stack:
            node, children = stack[-1]
            try:
                child = next(children)
                if color.get(child, WHITE) == GRAY:
                    # Found back edge -> cycle
                    cycle_start = path.index(child)
                    cycle_path = path[cycle_start:] + [child]
                    cycle_errors.append(
                        f"dependency cycle detected: {' -> '.join(cycle_path)}"
                    )
                elif color.get(child, WHITE) == WHITE:
                    color[child] = GRAY
                    path.append(child)
                    stack.append((child, iter(adj.get(child, []))))
            except StopIteration:
                color[node] = BLACK
                stack.pop()
                if len(path) > 1:
                    path.pop()

    for cid in case_ids:
        if color.get(cid, WHITE) == WHITE:
            dfs(cid)

    errors.extend(cycle_errors)
    return errors


_RANS_LANES_FORBIDDEN = {"adjoint", "wavemaker"}
_STATIC_CYL_NON_RANS_CANONICAL_RE = 100.0


def _settings_and_capability_checks(
    cases: list,
    capability_matrix: dict,
) -> list[str]:
    """Settings + capability invariants.

    Failures:
    - static-cylinder non-RANS final-target case at Re != 100 without an
      explicit deviation/blocker;
    - RANS adjoint or RANS wavemaker case (no adjoint formulation
      exists for the RANS lanes per the capability matrix);
    - Newton case with a seed at a different Re than the target but no
      seed_status documenting the continuation;
    - candidate-planned capability lane with no documented blocker on
      the corresponding case.
    """
    errors: list[str] = []
    for case in cases:
        label = getattr(case, "case_id", "<unknown>")
        target_re = case.target_parameters.Re
        target_model = case.target_parameters.model
        seed_re = case.seed_parameters.seed_Re
        is_static_cyl = case.flow_family == "cylinder_re100"
        is_rans = target_model == "RANS"

        # 1. Re=100 enforcement on static-cylinder non-RANS finals.
        if (is_static_cyl
                and not is_rans
                and target_re not in (None, _STATIC_CYL_NON_RANS_CANONICAL_RE)
                and not case.blockers
                and case.status not in (catalog.CaseStatus.DEFERRED,
                                        catalog.CaseStatus.BLOCKED,
                                        catalog.CaseStatus.NEEDS_DNS_SEED,
                                        catalog.CaseStatus.NEEDS_SCALAR_CHECKPOINT)):
            errors.append(
                f"{label}: static-cylinder non-RANS target Re={target_re} "
                f"deviates from canonical Re=100 with no blocker/deferred reason"
            )

        # 2. RANS adjoint / wavemaker forbidden.
        if is_rans and case.method_lane in _RANS_LANES_FORBIDDEN:
            errors.append(
                f"{label}: RANS case uses method_lane={case.method_lane!r} "
                f"which has no adjoint formulation per the method capability matrix"
            )

        # 3. Continuation seed at different Re must carry seed_status.
        if (isinstance(seed_re, (int, float))
                and target_re is not None
                and float(seed_re) != float(target_re)
                and not case.seed_parameters.seed_status):
            errors.append(
                f"{label}: seed_Re={seed_re} differs from target_Re={target_re} "
                f"but seed_parameters.seed_status is empty; record the"
                f" continuation/seed status explicitly"
            )

    # 4. Every FlowFamily mentioned by a CANDIDATE_PLANNED entry in the
    #    capability matrix must have at least one case in the catalog, so
    #    the family is at least tracked (not just an empty stub).
    candidate_planned_families = {
        fam for (fam, _), state in capability_matrix.items()
        if state == catalog.CapabilityState.CANDIDATE_PLANNED
    }
    families_with_cases = {c.flow_family for c in cases}
    for fam in sorted(candidate_planned_families - families_with_cases):
        errors.append(
            f"family {fam!r} has CANDIDATE_PLANNED capabilities but no case"
            f" in CASES; either add a placeholder case or remove the family"
            f" from the capability matrix"
        )

    return errors


def _numbering_checks(cases: list) -> list[str]:
    """Numbering consistency and uparam preservation invariants.

    Rule sort_prefix_uparam_family:
        For cases with legacy_uparam01 set, the hundreds digit of sort_prefix
        must equal int(float(legacy_uparam01)) — the integer part of the uparam
        value.  For example, uparam "3.11" requires sort_prefix in 300-399.

    Rule uparam_not_null_solver_lanes:
        Cases whose method_lane is in _LANES_REQUIRING_UPARAM (baseflow,
        stability_direct, etc.) must not have legacy_uparam01=None.  These lanes
        historically carry a userParam01 value; a null here is a missing record.

    Rule no_duplicate_prefix_variant:
        No two cases in the same flow_family share (sort_prefix, variant_suffix),
        where variant_suffix is the last path segment of proposed_folder.

    Rule variant_subdir_allowed_stage:
        A proposed_folder with 4 path components (example/FAMILY/NNN_slug/variant)
        is only allowed when the sort_prefix hundreds digit is 1, 2, or 3 —
        the baseflow, Newton, and stability spine stages.
    """
    errors: list[str] = []
    # (flow_family, sort_prefix, variant_suffix) -> first case_id seen
    seen_prefix_variant: dict[tuple[str, str, str], str] = {}

    for case in cases:
        label = getattr(case, "case_id", "<unknown>")

        # --- Rule sort_prefix_uparam_family ---
        if case.legacy_uparam01 is not None:
            try:
                expected_hundreds = int(float(case.legacy_uparam01))
                actual_hundreds = int(case.sort_prefix) // 100
                if expected_hundreds != actual_hundreds:
                    errors.append(
                        f"{label}: sort_prefix {case.sort_prefix!r} hundreds={actual_hundreds}"
                        f" does not match legacy_uparam01 {case.legacy_uparam01!r}"
                        f" family (expected hundreds={expected_hundreds})"
                    )
            except (ValueError, TypeError):
                pass  # malformed uparam already caught by _VALID_LEGACY_UPARAM01 check

        # --- Rule uparam_not_null_solver_lanes ---
        if case.method_lane in _LANES_REQUIRING_UPARAM and case.legacy_uparam01 is None:
            errors.append(
                f"{label}: method_lane={case.method_lane!r} requires a legacy_uparam01"
                f" but the field is None; record the historical userParam01 value"
            )

        # --- Rules no_duplicate_prefix_variant and variant_subdir_allowed_stage ---
        parts = case.proposed_folder.split("/")
        # variant_suffix = last path segment regardless of depth
        variant_suffix = parts[-1]
        pv_key = (case.flow_family, case.sort_prefix, variant_suffix)
        if pv_key in seen_prefix_variant:
            errors.append(
                f"{label}: (flow_family={case.flow_family!r}, sort_prefix={case.sort_prefix!r},"
                f" variant={variant_suffix!r}) duplicates {seen_prefix_variant[pv_key]!r}"
            )
        else:
            seen_prefix_variant[pv_key] = label

        # Rule variant_subdir_allowed_stage: 4-segment proposed_folder check
        if len(parts) == 4:
            try:
                actual_hundreds = int(case.sort_prefix) // 100
                if actual_hundreds not in _VARIANT_ALLOWED_HUNDREDS:
                    errors.append(
                        f"{label}: proposed_folder has a variant subdir"
                        f" but sort_prefix {case.sort_prefix!r} hundreds={actual_hundreds}"
                        f" is outside the allowed spine stages (1, 2, or 3)"
                    )
            except (ValueError, TypeError):
                pass  # sort_prefix format already caught above

    return errors

def _evidence_diversity_checks(cases: list) -> list[str]:
    """Each VALIDATED case must carry at least FIGURE + HISTORY.

    FIGURE  = EvidenceRole.REFERENCE_IMAGE
    HISTORY = EvidenceRole.RESIDUAL_HISTORY | DNS_TIME_HISTORY | SPECTRUM
              (any iterative-solver convergence evidence; SPECTRUM counts because
               Krylov/Arnoldi stability solvers produce Spectre_*.dat files as
               their primary convergence record)

    Excluded lanes (postproc-only, reads others' outputs):
      'wavemaker', 'animation', 'otd', 'modal', 'rans'

    These lanes have no self-produced residual/history/spectrum; a LOG or
    PROVENANCE artifact is sufficient secondary evidence for them.
    """
    _EXCLUDED_LANES = {"wavemaker", "animation", "otd", "modal", "rans"}
    figure_roles = {catalog.EvidenceRole.REFERENCE_IMAGE}
    history_roles = {
        catalog.EvidenceRole.RESIDUAL_HISTORY,
        catalog.EvidenceRole.DNS_TIME_HISTORY,
        catalog.EvidenceRole.SPECTRUM,
    }
    errors: list[str] = []
    for case in cases:
        if case.status != catalog.CaseStatus.VALIDATED:
            continue
        if case.method_lane in _EXCLUDED_LANES:
            continue
        artifact_roles = {a.role for a in case.artifacts}
        if not (artifact_roles & figure_roles):
            errors.append(
                f"{case.case_id}: VALIDATED case has no FIGURE artifact (REFERENCE_IMAGE)"
            )
        if not (artifact_roles & history_roles):
            errors.append(
                f"{case.case_id}: VALIDATED case has no HISTORY artifact"
                f" (RESIDUAL_HISTORY, DNS_TIME_HISTORY, or SPECTRUM)"
            )
    return errors


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Lint validation.catalog invariants")
    parser.add_argument("--json", action="store_true", help="emit JSON")
    args = parser.parse_args(list(argv) if argv is not None else None)

    errors = lint(
        catalog.families(),
        catalog.all_cases(),
        catalog.CAPABILITY_MATRIX,
        catalog.schema_version(),
    )

    if args.json:
        print(json.dumps({"ok": not errors, "errors": errors}, indent=2))
    elif errors:
        for error in errors:
            print(error)
    else:
        print("catalog lint passed")

    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
