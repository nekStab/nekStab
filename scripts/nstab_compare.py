#!/usr/bin/env python3
"""nstab_compare.py — overlay residual decay across variant subdirs of a stage.

Usage:
    scripts/nstab_compare.py <stage_dir> [--no-field-diff] [--out <png_dir>]

Discovers variant subdirs of <stage_dir> (skips _-prefixed and dotdirs). For
each variant, reads a residual log (search order: residual.log, residu_sfd.dat,
residu_dmt.dat, residu_tdf.dat, residu_boost.dat). Expects 3 cols
(time, residual, rate) or 2 cols (time, residual — rate computed via finite
diff). Emits <stage_dir>/residual_decay.png (or --out dir).

Optional .compare.yaml in <stage_dir>:
    title:      str
    x_label:    str
    y_label:    str
    log_y:      bool (default true)
    field_diff: bool (default false; .fld reader is TBD)

TODO: .fld field-diff (needs pymech or in-house reader).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- YAML loader — prefer PyYAML; fall back to a minimal flat key:value parser ---
try:
    import yaml as _yaml

    def _load_yaml(path: Path) -> dict:
        with path.open() as fh:
            return _yaml.safe_load(fh) or {}

except ImportError:  # pragma: no cover

    def _load_yaml(path: Path) -> dict:  # type: ignore[misc]
        """Fallback: parse flat 'key: value' YAML (no nesting, no lists)."""
        result: dict = {}
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                continue
            k, _, v = line.partition(":")
            v = v.strip()
            if v.lower() in ("true", "yes"):
                v = True  # type: ignore[assignment]
            elif v.lower() in ("false", "no"):
                v = False  # type: ignore[assignment]
            result[k.strip()] = v
        return result


# --- Residual log search order ---
LOG_NAMES = [
    "residual.log",
    "residu_sfd.dat",
    "residu_dmt.dat",
    "residu_tdf.dat",
    "residu_boost.dat",
]


def _find_log(variant_dir: Path) -> Optional[Path]:
    for name in LOG_NAMES:
        p = variant_dir / name
        if p.is_file():
            return p
    return None


def _load_log(log_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (time, residual) arrays from a 2- or 3-column log."""
    try:
        data = np.loadtxt(log_path, comments=("#", "!"))
    except Exception as exc:
        raise ValueError(f"Cannot parse {log_path}: {exc}") from exc
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(
            f"{log_path}: need ≥2 columns (time, residual); got {data.shape[1]}"
        )
    t = data[:, 0]
    r = np.abs(data[:, 1])
    return t, r


# --- Main ---
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Overlay residual decay across variant subdirs of a nekStab stage.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("stage_dir", type=Path, help="Stage directory containing variant subdirs.")
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        metavar="PNG_DIR",
        help="Output directory for residual_decay.png (default: stage_dir).",
    )
    p.add_argument(
        "--no-field-diff",
        action="store_true",
        help="Suppress field diff even if .compare.yaml requests it.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    stage_dir: Path = args.stage_dir.resolve()

    if not stage_dir.is_dir():
        sys.exit(f"ERROR: stage_dir does not exist or is not a directory: {stage_dir}")

    # --- load optional config ---
    cfg: dict = {}
    yaml_path = stage_dir / ".compare.yaml"
    if yaml_path.is_file():
        cfg = _load_yaml(yaml_path)

    title: str = cfg.get("title", f"{stage_dir.name} — variant comparison")
    x_label: str = cfg.get("x_label", "Time")
    y_label: str = cfg.get("y_label", "Residual")
    log_y: bool = bool(cfg.get("log_y", True))
    field_diff: bool = bool(cfg.get("field_diff", False))

    if args.no_field_diff:
        field_diff = False

    if field_diff:
        print("NOTE: field_diff requested but .fld reader is TBD — skipping.")

    # --- discover variant subdirs ---
    variants = sorted(
        d for d in stage_dir.iterdir()
        if d.is_dir() and not d.name.startswith(("_", "."))
    )
    if not variants:
        sys.exit(f"ERROR: No variant subdirs found in {stage_dir}")

    # --- collect data ---
    found: list[tuple[str, np.ndarray, np.ndarray]] = []
    skipped: list[str] = []
    for vdir in variants:
        log = _find_log(vdir)
        if log is None:
            skipped.append(vdir.name)
            continue
        try:
            t, r = _load_log(log)
        except ValueError as exc:
            print(f"WARNING: {exc} — skipping {vdir.name}")
            skipped.append(vdir.name)
            continue
        found.append((vdir.name, t, r))

    if skipped:
        print(f"Skipped variants (no log found): {', '.join(skipped)}")
    if not found:
        sys.exit("ERROR: No residual logs found in any variant subdir.")

    # --- plot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, t, r in found:
        ax.plot(t, r, label=name)

    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    fig.tight_layout()

    out_dir = args.out.resolve() if args.out else stage_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "residual_decay.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
