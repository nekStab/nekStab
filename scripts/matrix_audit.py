#!/usr/bin/env python3
"""matrix_audit.py — audit nekStab example/ vs the locked geometry × NNN matrix.

Usage:
    scripts/matrix_audit.py             # markdown table to stdout
    scripts/matrix_audit.py --json      # machine-readable
    scripts/matrix_audit.py --missing   # only the missing cells

Mirrors plan.md §1.3. Update the GEOMETRIES and STAGES dicts when those change.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Locked matrix — mirrors plan.md §1.2
# ---------------------------------------------------------------------------
STAGES: dict[str, dict[str, str]] = {
    "000": {"label": "DNS",                        "uparam01": "0"},
    "110": {"label": "baseflow SFD",               "uparam01": "1.1"},
    "120": {"label": "baseflow BoostConv",         "uparam01": "1.2"},
    "130": {"label": "baseflow DMT",               "uparam01": "1.3"},
    "140": {"label": "baseflow TDF",               "uparam01": "1.4"},
    "210": {"label": "baseflow Newton",            "uparam01": "2/2.1/2.2"},
    "310": {"label": "stability direct",           "uparam01": "3.1"},
    "311": {"label": "stability direct Floquet",   "uparam01": "3.11"},
    "320": {"label": "stability adjoint",          "uparam01": "3.2"},
    "321": {"label": "stability adjoint Floquet",  "uparam01": "3.21"},
    "330": {"label": "transient growth",           "uparam01": "3.3"},
    "411": {"label": "wavemaker / energy budget",  "uparam01": "4.1/4.11"},
    "500": {"label": "OTD",                        "uparam01": "5"},
    "600": {"label": "modal POD",                  "uparam01": "postproc"},
    "610": {"label": "modal DMD",                  "uparam01": "postproc"},
    "620": {"label": "modal SPOD",                 "uparam01": "postproc"},
}

# SFD variant subdirs to check under any 110_* directory
SFD_VARIANTS = ["akervik", "casacuberta", "dyn", "dyn_oifs"]

# ---------------------------------------------------------------------------
# Locked matrix — mirrors plan.md §1.3
# ---------------------------------------------------------------------------
_FULL = ["000", "110", "120", "130", "140", "210", "310", "311",
         "320", "321", "330", "411", "500", "600", "610", "620"]

GEOMETRIES: dict[str, list[str]] = {
    "cylinder_re100":   _FULL,
    "cylinder_re1m":    ["000", "310"],
    "naca0012":         _FULL,
    "thermosyphon":     ["000", "210", "310", "320", "411"],
    "flip_flop":        ["000", "210", "311", "321", "411"],
    "tpjet":            ["000", "210", "311", "321", "411"],
    "lid_driven":       _FULL,
    "cubic_cavity":     _FULL,
    "cubic_cavity_upo": ["000", "311"],
    "blasius":          _FULL,
    "poiseuille":       _FULL,
    "bfs":              ["000", "330"],
    "slot_FST":         ["000"],
}

# Disk-name aliases: plan.md geometry name → actual example/ subdir name(s) to try.
# First match wins. If none match, the cell is TODO.
DISK_ALIASES: dict[str, list[str]] = {
    "cylinder_re100": ["cylinder_re100", "cylinder"],
    "cylinder_re1m":  ["cylinder_re1m", "cylinder/RANS"],
    "poiseuille":     ["poiseuille", "poiseuille_OTD", "poiseuille_RANS"],
    "bfs":            ["bfs", "back_fstep"],
}


def find_repo_root(start: Path) -> Path:
    """Walk up from *start* until a .git directory is found."""
    here = start.resolve()
    for parent in [here, *here.parents]:
        if (parent / ".git").exists():
            return parent
    # Fallback: parent of scripts/
    return here.parent


def resolve_geom_dir(root: Path, geom: str) -> Optional[Path]:
    """Return the first existing directory for *geom* under root/example/."""
    candidates = DISK_ALIASES.get(geom, [geom])
    for alias in candidates:
        p = root / "example" / alias
        if p.is_dir():
            return p
    return None


def find_nnn_dir(geom_dir: Path, nnn: str) -> Optional[Path]:
    """Return the first directory matching <nnn>_* under *geom_dir*."""
    matches = sorted(geom_dir.glob(f"{nnn}_*"))
    for m in matches:
        if m.is_dir():
            return m
    return None


def check_sfd_variants(nnn_dir: Path) -> list[str]:
    """Return list of SFD variant subdirs present under *nnn_dir*."""
    found = []
    for v in SFD_VARIANTS:
        if (nnn_dir / v).is_dir():
            found.append(v)
    return found


# ---------------------------------------------------------------------------
# Cell record
# ---------------------------------------------------------------------------

def audit(root: Path) -> list[dict]:
    records: list[dict] = []

    for geom, stages in GEOMETRIES.items():
        geom_dir = resolve_geom_dir(root, geom)

        for nnn in stages:
            stage_label = STAGES[nnn]["label"]

            if geom_dir is None:
                records.append({
                    "geom": geom, "nnn": nnn, "stage": stage_label,
                    "status": "TODO", "path": None, "variants": [],
                })
                continue

            nnn_dir = find_nnn_dir(geom_dir, nnn)
            if nnn_dir is None:
                records.append({
                    "geom": geom, "nnn": nnn, "stage": stage_label,
                    "status": "TODO", "path": None, "variants": [],
                })
                continue

            rel = str(nnn_dir.relative_to(root))
            variants = check_sfd_variants(nnn_dir) if nnn == "110" else []
            records.append({
                "geom": geom, "nnn": nnn, "stage": stage_label,
                "status": "present", "path": rel, "variants": variants,
            })

    return records


# ---------------------------------------------------------------------------
# Output formatters
# ---------------------------------------------------------------------------

def _status_cell(rec: dict) -> str:
    if rec["status"] == "TODO":
        return "TODO"
    path = rec["path"] or ""
    variants = rec["variants"]
    if variants:
        return f"✓ {path} [{', '.join(variants)}]"
    return f"✓ {path}"


def fmt_markdown(records: list[dict]) -> str:
    header = "| Geometry | NNN | Stage | Status |\n"
    header += "|----------|-----|-------|--------|\n"
    rows = []
    for rec in records:
        rows.append(
            f"| {rec['geom']} | {rec['nnn']} | {rec['stage']} | {_status_cell(rec)} |"
        )
    return header + "\n".join(rows)


def fmt_missing(records: list[dict]) -> str:
    header = "| Geometry | NNN | Stage |\n"
    header += "|----------|-----|-------|\n"
    rows = []
    for rec in records:
        if rec["status"] == "TODO":
            rows.append(f"| {rec['geom']} | {rec['nnn']} | {rec['stage']} |")
    if not rows:
        return header + "*(none — all cells present)*"
    return header + "\n".join(rows)


def fmt_json(records: list[dict]) -> str:
    out = []
    for rec in records:
        out.append({
            "geom":    rec["geom"],
            "nnn":     rec["nnn"],
            "stage":   rec["stage"],
            "status":  rec["status"],
            "path":    rec["path"],
            "variants": rec["variants"],
        })
    return json.dumps(out, indent=2)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit nekStab example/ vs the locked geometry × NNN matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--json",    action="store_true", help="emit machine-readable JSON")
    group.add_argument("--missing", action="store_true", help="emit only missing (TODO) cells")
    args = parser.parse_args()

    root = find_repo_root(Path(__file__).parent)
    records = audit(root)

    if args.json:
        print(fmt_json(records))
    elif args.missing:
        print(fmt_missing(records))
    else:
        print(fmt_markdown(records))


if __name__ == "__main__":
    main()
