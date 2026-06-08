#!/usr/bin/env python3
"""Gather validation plot artifacts into validation/figures/.

Default mode is catalog-driven. Legacy mode keeps the original walk-based
collection behavior for plot.py/plot.png example directories.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent
NEKSTAB_ROOT = EXAMPLE_DIR.parent
OUTPUT_DIR = NEKSTAB_ROOT / "validation" / "figures"

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    from validation.catalog import all_cases, EvidenceRole
except ImportError as exc:
    all_cases = None
    EvidenceRole = None
    CATALOG_IMPORT_ERROR = exc
else:
    CATALOG_IMPORT_ERROR = None


SKIP_DIRS = {".venv", "__pycache__", "modal"}


def find_python():
    """Find a Python interpreter that has pymech installed."""
    candidates = [sys.executable]
    for venv in [Path.home() / ".venv", NEKSTAB_ROOT / ".venv", EXAMPLE_DIR / ".venv"]:
        python = venv / "bin" / "python"
        if python.exists():
            candidates.insert(0, str(python))

    for python in candidates:
        try:
            result = subprocess.run(
                [python, "-c", "import pymech"],
                capture_output=True,
                timeout=5,
            )
        except Exception:
            continue
        if result.returncode == 0:
            return python

    return sys.executable


PYTHON = find_python()


def find_par_file(case_dir):
    """Find the .par file in case_dir or nearest parent under example/."""
    directory = case_dir
    while directory != EXAMPLE_DIR.parent:
        par_files = list(directory.glob("*.par"))
        if par_files:
            return par_files[0]
        directory = directory.parent
    return None


def extract_params(par_file):
    """Extract Re and/or Ra from a .par file."""
    if par_file is None:
        return {}

    text = par_file.read_text()
    params = {}

    match = re.search(r"(?i)viscosity\s*=\s*([-.\deE]+)", text)
    if match:
        viscosity = float(match.group(1))
        if viscosity < 0:
            params["Re"] = int(abs(viscosity))

    match = re.search(r"(?i)userParam06\s*=\s*([-.\deE]+)", text)
    if match:
        rayleigh = float(match.group(1))
        if rayleigh > 0:
            params["Ra"] = int(rayleigh)

    if params.get("Ra") == 0:
        del params["Ra"]

    return params


def make_name(case_dir, params):
    """Build descriptive filename from directory path and parameters."""
    relative_path = case_dir.relative_to(EXAMPLE_DIR)
    parts = [part for part in relative_path.parts if part not in SKIP_DIRS]
    base = "_".join(parts)

    for key in ("Re", "Ra"):
        if key in params:
            tag = key + str(params[key])
            if tag.lower() not in base.lower():
                base += "_" + tag

    return base + ".png"


def copy_or_report(src, dest, list_only, label):
    if list_only:
        print("  " + label.ljust(55) + " <- " + str(src.relative_to(NEKSTAB_ROOT)))
        return

    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() != dest.resolve():
        shutil.copy2(src, dest)


def catalog_roles():
    return {
        EvidenceRole.REFERENCE_IMAGE,
        EvidenceRole.SPECTRUM,
        EvidenceRole.ANIMATION,
        EvidenceRole.MODE_SHAPE,
        EvidenceRole.DNS_TIME_HISTORY,
    }


def catalog_summary_label(role_counts):
    return (
        "Catalog-driven: "
        + str(role_counts.get(EvidenceRole.REFERENCE_IMAGE, 0))
        + " figures + "
        + str(role_counts.get(EvidenceRole.SPECTRUM, 0))
        + " spectra + "
        + str(role_counts.get(EvidenceRole.ANIMATION, 0))
        + " animations + "
        + str(role_counts.get(EvidenceRole.MODE_SHAPE, 0))
        + " mode shapes + "
        + str(role_counts.get(EvidenceRole.DNS_TIME_HISTORY, 0))
        + " dns-history collected"
    )


def collect_catalog(list_only=False):
    """Collect registered catalog artifacts into validation/figures/."""
    if all_cases is None:
        print("WARNING: validation.catalog not importable; falling back to legacy mode.")
        print("  Import error: " + str(CATALOG_IMPORT_ERROR))
        collect_legacy(list_only=list_only)
        return

    if not list_only:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    wanted_roles = catalog_roles()
    role_counts = {
        EvidenceRole.REFERENCE_IMAGE: 0,
        EvidenceRole.SPECTRUM: 0,
        EvidenceRole.ANIMATION: 0,
        EvidenceRole.MODE_SHAPE: 0,
        EvidenceRole.DNS_TIME_HISTORY: 0,
    }
    missing = []

    for case in all_cases():
        for artifact in case.artifacts:
            if artifact.role not in wanted_roles:
                continue

            src = NEKSTAB_ROOT / artifact.path
            dest_name = Path(artifact.path).name
            dest = OUTPUT_DIR / dest_name

            if not src.exists():
                missing.append(artifact.path)
                continue

            copy_or_report(src, dest, list_only, dest_name)
            role_counts[artifact.role] += 1

    print(catalog_summary_label(role_counts))
    print("Missing (artifact registered but file absent): " + str(len(missing)))
    for path in missing:
        print("  " + path)


def collect_legacy(generate=False, list_only=False):
    """Original walk-based collection logic."""
    plot_scripts = sorted(EXAMPLE_DIR.rglob("plot.py"))
    plot_scripts = [
        path for path in plot_scripts
        if not any(skip in path.parts for skip in SKIP_DIRS)
    ]

    if not list_only:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    collected = []
    skipped = []

    for script in plot_scripts:
        case_dir = script.parent
        plot_png = case_dir / "plot.png"

        if generate:
            print(
                "  Running " + str(script.relative_to(NEKSTAB_ROOT)) + " ...",
                end=" ",
                flush=True,
            )
            try:
                result = subprocess.run(
                    [PYTHON, str(script)],
                    capture_output=True,
                    text=True,
                    timeout=120,
                    cwd=str(case_dir),
                    env={**os.environ, "MPLBACKEND": "Agg"},
                )
            except subprocess.TimeoutExpired:
                print("TIMEOUT")
            except Exception as exc:
                print("ERROR (" + str(exc) + ")")
            else:
                if result.returncode == 0:
                    print("OK")
                else:
                    if result.stderr:
                        error = result.stderr.strip().split("\n")[-1]
                    else:
                        error = "unknown"
                    print("FAIL (" + error + ")")

        if not plot_png.exists():
            skipped.append(case_dir.relative_to(EXAMPLE_DIR))
            continue

        par_file = find_par_file(case_dir)
        params = extract_params(par_file)
        name = make_name(case_dir, params)

        if list_only:
            print("  " + name.ljust(55) + " <- " + str(case_dir.relative_to(EXAMPLE_DIR)) + "/")
        else:
            shutil.copy2(plot_png, OUTPUT_DIR / name)
            collected.append(name)

    if list_only:
        if skipped:
            print("\nNo plot.png yet (" + str(len(skipped)) + "):")
            for path in skipped:
                print("  " + str(path) + "/")
        return

    print("\nCollected " + str(len(collected)) + " plots -> " + str(OUTPUT_DIR.relative_to(NEKSTAB_ROOT)) + "/")
    for name in collected:
        print("  " + name)
    if skipped:
        print("\nSkipped " + str(len(skipped)) + " (no plot.png):")
        for path in skipped:
            print("  " + str(path) + "/")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Gather plot files into validation/figures/.",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--catalog", action="store_true", help="Use catalog-driven collection")
    mode.add_argument("--legacy", action="store_true", help="Use legacy walk-based collection")
    parser.add_argument(
        "--generate",
        action="store_true",
        help="Legacy only: run plot.py scripts before collecting",
    )
    parser.add_argument(
        "--list",
        "-l",
        action="store_true",
        help="Show what would be collected without copying",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        dest="list",
        help="Alias for --list",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.legacy:
        collect_legacy(generate=args.generate, list_only=args.list)
        return

    if args.generate:
        print("WARNING: --generate is legacy-only and is ignored in catalog mode.")

    collect_catalog(list_only=args.list)


if __name__ == "__main__":
    main()
