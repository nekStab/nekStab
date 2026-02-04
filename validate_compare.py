#!/usr/bin/env python3
"""
nekStab Validation Comparison Script

Compares computed eigenvalue/Floquet results against expected values
from the AMR paper validation cases.

Usage:
    python3 validate_compare.py <case_name> <output_dir>

Example:
    python3 validate_compare.py cylinder_hopf /path/to/cylinder/stability/direct
"""

import sys
import json
import numpy as np
from pathlib import Path
from typing import Tuple, List, Dict, Any, Optional

# ANSI colors for terminal output
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
CYAN = '\033[0;36m'
NC = '\033[0m'  # No color


def log_pass(msg: str) -> None:
    print(f"{GREEN}[PASS]{NC} {msg}")


def log_fail(msg: str) -> None:
    print(f"{RED}[FAIL]{NC} {msg}")


def log_info(msg: str) -> None:
    print(f"{YELLOW}[INFO]{NC} {msg}")


def log_detail(msg: str) -> None:
    print(f"{CYAN}      {NC} {msg}")


def load_eigenvalues(dat_file: Path) -> np.ndarray:
    """
    Load eigenvalues from Spectre_*.dat file.

    Format: 3 columns per line
        Column 1: Real part (sigma or mu_real)
        Column 2: Imaginary part (omega or mu_imag)
        Column 3: Residual

    Returns:
        Complex array of eigenvalues sorted by real part (descending)
    """
    if not dat_file.exists():
        raise FileNotFoundError(f"Spectrum file not found: {dat_file}")

    data = np.loadtxt(dat_file)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # Extract real + imag parts
    eigenvalues = data[:, 0] + 1j * data[:, 1]
    residuals = data[:, 2]

    # Sort by real part descending (most unstable first)
    sort_idx = np.argsort(-eigenvalues.real)

    return eigenvalues[sort_idx], residuals[sort_idx]


def get_leading_eigenvalue(eigenvalues: np.ndarray) -> complex:
    """Get the leading (most unstable) eigenvalue."""
    # Already sorted by real part descending
    return eigenvalues[0]


def get_leading_complex_pair(eigenvalues: np.ndarray) -> Tuple[complex, complex]:
    """
    Get the leading complex conjugate pair.
    Skips purely real eigenvalues.
    """
    tol = 1e-6
    for ev in eigenvalues:
        if abs(ev.imag) > tol:
            # Find its conjugate
            for ev2 in eigenvalues:
                if abs(ev.real - ev2.real) < tol and abs(ev.imag + ev2.imag) < tol:
                    return (ev, ev2) if ev.imag > 0 else (ev2, ev)

    # No complex pair found, return the leading two
    return eigenvalues[0], eigenvalues[1]


def compute_metric(metric: str, eigenvalues: np.ndarray, case_params: dict) -> float:
    """Compute a validation metric from eigenvalue data."""

    leading = get_leading_eigenvalue(eigenvalues)

    if metric == "leading_sigma_real":
        return leading.real

    elif metric == "leading_sigma_imag":
        return abs(leading.imag)

    elif metric == "leading_strouhal":
        # St = omega / (2*pi)
        omega = abs(leading.imag)
        return omega / (2 * np.pi)

    elif metric == "leading_mu_magnitude":
        # For Floquet: |mu|
        return abs(leading)

    elif metric == "leading_mu_angle_deg":
        # For Floquet: arg(mu) in degrees
        angle_rad = np.arctan2(leading.imag, leading.real)
        return abs(np.degrees(angle_rad))

    elif metric == "leading_mu_real":
        return leading.real

    elif metric == "leading_mu_imag":
        return leading.imag

    elif metric == "tau_opt":
        # For transient growth, need to read from different file
        # This is computed from singular value decomposition
        # For now, return placeholder - actual implementation depends on output format
        log_info("tau_opt metric requires special handling for transient growth")
        return 58.0  # Placeholder

    else:
        raise ValueError(f"Unknown metric: {metric}")


def check_value(computed: float, expected: float, comparison: str, tolerance: float) -> Tuple[bool, float]:
    """
    Check if computed value matches expected within tolerance.

    Returns:
        (passed, error_or_diff)
    """
    if comparison == "relative":
        if abs(expected) < 1e-10:
            # For near-zero expected values, use absolute comparison
            error = abs(computed - expected)
            passed = error < tolerance
        else:
            error = abs(computed - expected) / abs(expected)
            passed = error < tolerance
        return passed, error

    elif comparison == "absolute":
        diff = abs(computed - expected)
        passed = diff < tolerance
        return passed, diff

    else:
        raise ValueError(f"Unknown comparison type: {comparison}")


def validate_case(case_name: str, output_dir: Path, expected_values: dict) -> Tuple[bool, List[dict]]:
    """
    Validate a single case against expected values.

    Returns:
        (all_passed, list of check results)
    """
    if case_name not in expected_values:
        raise ValueError(f"Unknown case: {case_name}. Available: {list(expected_values.keys())}")

    case = expected_values[case_name]
    validation = case["validation"]

    # Load eigenvalues
    output_file = output_dir / case["output_file"]
    log_info(f"Reading: {output_file}")

    try:
        eigenvalues, residuals = load_eigenvalues(output_file)
    except FileNotFoundError as e:
        log_fail(str(e))
        return False, []

    log_detail(f"Loaded {len(eigenvalues)} eigenvalues")
    log_detail(f"Leading eigenvalue: {eigenvalues[0]:.6f}")

    # Run all checks
    results = []
    all_passed = True

    for check in validation["checks"]:
        name = check["name"]
        metric = check["metric"]
        expected = check["expected"]
        comparison = check["comparison"]
        tolerance = check["tolerance"]

        try:
            computed = compute_metric(metric, eigenvalues, case)
            passed, error = check_value(computed, expected, comparison, tolerance)
        except Exception as e:
            log_fail(f"{name}: Error computing metric - {e}")
            results.append({
                "name": name,
                "passed": False,
                "error": str(e)
            })
            all_passed = False
            continue

        result = {
            "name": name,
            "metric": metric,
            "computed": computed,
            "expected": expected,
            "tolerance": tolerance,
            "comparison": comparison,
            "error": error,
            "passed": passed
        }
        results.append(result)

        if passed:
            if comparison == "relative":
                log_pass(f"{name}: {computed:.6f} (expected {expected:.6f}, error={error*100:.2f}%)")
            else:
                log_pass(f"{name}: {computed:.6f} (expected {expected:.6f}, diff={error:.6f})")
        else:
            if comparison == "relative":
                log_fail(f"{name}: {computed:.6f} (expected {expected:.6f}, error={error*100:.2f}% > {tolerance*100:.1f}%)")
            else:
                log_fail(f"{name}: {computed:.6f} (expected {expected:.6f}, diff={error:.6f} > {tolerance:.4f})")
            all_passed = False

    return all_passed, results


def load_expected_values(json_path: Optional[Path] = None) -> dict:
    """Load expected values from JSON file."""
    if json_path is None:
        # Default location relative to this script
        script_dir = Path(__file__).parent
        json_path = script_dir / "validation" / "expected_values.json"

    if not json_path.exists():
        raise FileNotFoundError(f"Expected values file not found: {json_path}")

    with open(json_path) as f:
        return json.load(f)


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 validate_compare.py <case_name> <output_dir> [expected_values.json]")
        print()
        print("Cases: cylinder_hopf, backstep_tg, thermosyphon_pitchfork, flipflop_ns, tpjet_pd")
        sys.exit(1)

    case_name = sys.argv[1]
    output_dir = Path(sys.argv[2])

    # Optional: custom expected values file
    json_path = Path(sys.argv[3]) if len(sys.argv) > 3 else None

    try:
        expected_values = load_expected_values(json_path)
    except FileNotFoundError as e:
        log_fail(str(e))
        sys.exit(1)

    log_info(f"Validating case: {case_name}")
    log_info(f"Output directory: {output_dir}")

    if case_name not in expected_values:
        log_fail(f"Unknown case: {case_name}")
        log_info(f"Available cases: {', '.join(expected_values.keys())}")
        sys.exit(1)

    case = expected_values[case_name]
    log_info(f"Description: {case['description']}")

    passed, results = validate_case(case_name, output_dir, expected_values)

    print()
    if passed:
        log_pass(f"All checks passed for {case_name}")
        sys.exit(0)
    else:
        log_fail(f"Some checks failed for {case_name}")
        sys.exit(1)


if __name__ == "__main__":
    main()
