"""Stop rule of the Arnoldi pass on non-normal operators.

The rule is ks_leading_settled in src/ks_select.f90. NSTILL is set in
src/krylov_decomposition.f90.

A small residual bounds the error of a Ritz value only through the condition
number kappa of the eigenvalue. The error is up to kappa times the residual.
The forced jet of tpjet_re2000 has kappa = 4e5. Its adjoint run stopped at the
residual 1.1e-6 with a multiplier 5e-3 off. The run now also waits until the
leading values stand still for NSTILL steps.

The test builds random non-normal operators in Schur form. Each operator has
one leading multiplier of 1.16 with a large condition number, a cluster of
slow multipliers just inside the unit circle, and some coupling inside the
cluster. The matvec is exact or exact to 1e-9. Arnoldi runs in numpy. The
three leading Ritz values of every step go to the compiled ks_leading_settled.

For each operator the test compares the error of the leading value at the
step where a rule stops with the error after the whole basis. The whole basis
is what a run without early exit gives. A stop is worse when its error is
more than three times the error of the whole basis plus eigen_tol.

Run the file as a script to print the table for more operators:
    python3 validation/tests/test_ritz_stop_rule.py <trials> <seed> [<seed> ...]
The test needs gfortran and numpy. It skips when gfortran is missing.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"

TOL = 1.0e-6  # eigen_tol of the run
BASIS = 128  # k_dim of the run
MU = 1.16  # leading multiplier of every operator
WINDOWS = {"rest1": 1, "rest3": 3, "rest10": 10}  # steps of rest besides NSTILL
RULES = ("residual", "rest1", "rest3", "rest5", "rest10")

# A stop is worse than the whole basis when its error exceeds
# WORSE_FACTOR times the error of the whole basis plus WORSE_FLOOR.
WORSE_FACTOR = 3.0
WORSE_FLOOR = TOL

Trace = list[
    list[tuple[float, float, float]]
]  # per step: (re, im, residual) of 3 values


@dataclass(frozen=True)
class Stop:
    """The step where a rule stops and the error of the leading value there."""

    step: int
    error: float


@dataclass(frozen=True)
class Record:
    """One operator: its condition number and the stop of each rule."""

    kappa: float
    full: Stop
    stops: dict[str, Stop]

    def worse(self, rule: str) -> bool:
        """Return True when the rule stops with a larger error than the whole basis."""
        return self.stops[rule].error > WORSE_FACTOR * self.full.error + WORSE_FLOOR


def nstill_of_the_run() -> int:
    """Return NSTILL, the steps of rest that krylov_decomposition.f90 asks for."""
    text = (ROOT / "src" / "krylov_decomposition.f90").read_text()
    match = re.search(r"integer,\s*parameter\s*::\s*NSTILL\s*=\s*(\d+)", text)
    assert match is not None, "NSTILL not found in krylov_decomposition.f90"
    return int(match.group(1))


def build_driver(build: Path) -> Path:
    """Compile ks_select.f90 with the driver that prints the count of each step."""
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        raise FileNotFoundError("gfortran not found")
    flags = ["-fdefault-real-8", "-fdefault-double-8", "-ffree-line-length-none", "-O0"]
    sources = [
        ROOT / "src" / "argsort.f90",
        ROOT / "src" / "ks_select.f90",
        FORTRAN / "ks_settled_driver.f90",
    ]
    for src in sources:
        subprocess.run([gfortran, *flags, "-c", str(src)], cwd=build, check=True)
    exe = build / "driver"
    objects = [f"{src.stem}.o" for src in sources]
    subprocess.run([gfortran, *flags, "-o", str(exe), *objects], cwd=build, check=True)
    return exe


def make_operator(
    rng: np.random.Generator,
    n: int,
    kappa_target: float,
    coupling: float,
    width: float,
    centre: float,
) -> tuple[npt.NDArray[np.float64], float]:
    """Return a matrix of order n and the condition number of its leading multiplier.

    The Schur form has the multiplier MU at entry (0, 0). Row 0 couples it to
    all other entries. This spreads its left eigenvector and sets the
    condition number. A random orthogonal matrix hides the structure from
    Arnoldi.
    """
    n_cluster = n - 14
    sign = np.where(rng.random(n_cluster) < 0.5, 1.0, -1.0)
    cluster = (centre + width * (rng.random(n_cluster) - 0.5)) * sign
    rest = rng.uniform(0.2, 0.85, 13) * np.where(rng.random(13) < 0.5, 1.0, -1.0)
    schur = np.diag(np.concatenate(([MU], cluster, rest)))
    schur += np.triu(coupling * rng.standard_normal((n, n)) / np.sqrt(n), 1)
    schur[0, 1:] = kappa_target * 0.15 / np.sqrt(n) * rng.standard_normal(n - 1)
    schur[0, 0] = MU
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    # The right eigenvector of MU is e_0. The left one is the null vector of
    # (schur - MU I)^T with unit norm, so kappa = 1 / |left[0]|.
    left = np.linalg.svd((schur - MU * np.eye(n)).T)[2][-1]
    kappa = float(1.0 / abs(left[0]))
    return q @ schur @ q.T, kappa


def ritz_trace(
    matrix: npt.NDArray[np.float64], steps: int, noise: float, rng: np.random.Generator
) -> Trace:
    """Run Arnoldi and return, per step, the three leading Ritz entries.

    An entry is (real part, imaginary part, residual). The entries are sorted by
    decreasing modulus. The matvec adds a random error of relative size noise.
    """
    n = matrix.shape[0]
    basis = np.zeros((n, steps + 1))
    hess = np.zeros((steps + 1, steps))
    start = rng.standard_normal(n)
    basis[:, 0] = start / np.linalg.norm(start)
    trace: Trace = []
    for m in range(1, steps + 1):
        f = matrix @ basis[:, m - 1]
        if noise > 0:
            f = f + noise * np.linalg.norm(f) * rng.standard_normal(n) / np.sqrt(n)
        for _ in range(2):  # classical Gram-Schmidt, twice
            c = basis[:, :m].T @ f
            f = f - basis[:, :m] @ c
            hess[:m, m - 1] += c
        hess[m, m - 1] = np.linalg.norm(f)
        if hess[m, m - 1] < 1.0e-12:
            break
        basis[:, m] = f / hess[m, m - 1]
        w, y = np.linalg.eig(hess[:m, :m])
        res = np.abs(hess[m, m - 1] * y[m - 1, :]) / np.linalg.norm(y, axis=0)
        trace.append(
            [
                (float(w[i].real), float(w[i].imag), float(res[i]))
                for i in np.argsort(-np.abs(w))[:3]
            ]
        )
    return trace


def counts_of(driver: Path, trace: Trace, nstill: int) -> list[tuple[int, int]]:
    """Return (settled count, residual-only count) per step for one target."""
    lines = [f"1 {nstill} {TOL!r}"]
    for entries in trace:
        lines.append(str(len(entries)))
        lines += [f"{re!r} {im!r} {res!r}" for re, im, res in entries]
    run = subprocess.run(
        [str(driver)],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        check=True,
    )
    pairs = [tuple(int(x) for x in line.split()) for line in run.stdout.splitlines()]
    return [(p[0], p[1]) for p in pairs]


def stop_of(reached: list[bool], error: list[float]) -> Stop:
    """Return the first step with reached True, or the last step of the pass.

    The last step is where a pass ends when the rule never fires.
    """
    step = next((i + 1 for i, flag in enumerate(reached) if flag), len(reached))
    return Stop(step, error[step - 1])


def sweep(driver: Path, trials: int, seed: int) -> list[Record]:
    """Return one record per random operator."""
    rng = np.random.default_rng(seed)
    nstill = nstill_of_the_run()
    records = []
    for _ in range(trials):
        n = int(rng.choice([200, 300]))
        kappa_target = 10 ** rng.uniform(0, 6)
        coupling = float(rng.choice([0.0, 0.02, 0.1]))
        width = float(rng.choice([0.01, 0.05, 0.15]))
        centre = float(rng.choice([0.9, 0.95, 0.99]))
        noise = float(rng.choice([0.0, 1.0e-9]))
        matrix, kappa = make_operator(rng, n, kappa_target, coupling, width, centre)
        trace = ritz_trace(matrix, BASIS, noise, rng)
        error = [abs(complex(entries[0][0], entries[0][1]) - MU) for entries in trace]

        counts = counts_of(driver, trace, nstill)
        stops = {
            "residual": stop_of([c[1] >= 1 for c in counts], error),
            "rest5": stop_of([c[0] >= 1 for c in counts], error),
        }
        for rule, window in WINDOWS.items():
            by_window = counts_of(driver, trace, window)
            stops[rule] = stop_of([c[0] >= 1 for c in by_window], error)
        records.append(Record(kappa, Stop(len(trace), error[-1]), stops))
    return records


def table(records: list[Record]) -> str:
    """Return the summary of a sweep as text."""
    n = len(records)
    mean_full = np.mean([r.full.step for r in records])
    head = f"{n} operators, whole basis ends at mean step {mean_full:.1f}"
    rows = [head, "rule        worse than whole basis   error > 1e-4   mean stop step"]
    for rule in RULES:
        worse = sum(r.worse(rule) for r in records)
        large = sum(r.stops[rule].error > 1.0e-4 for r in records)
        mean_step = np.mean([r.stops[rule].step for r in records])
        share = 100.0 * worse / n
        rows.append(
            f"{rule:10s} {worse:6d} ({share:5.1f} %) {large:12d} {mean_step:16.1f}"
        )
    rows.append(
        f"whole basis: error > 1e-4 in {sum(r.full.error > 1.0e-4 for r in records)}"
    )
    return "\n".join(rows)


@pytest.fixture(scope="module")
def driver(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not found")
    return build_driver(tmp_path_factory.mktemp("ritz_stop_rule"))


@pytest.fixture(scope="module")
def records(driver: Path) -> list[Record]:
    return sweep(driver, trials=40, seed=2026)


def test_nstill_is_five() -> None:
    assert nstill_of_the_run() == 5


def test_the_sample_covers_strongly_non_normal_operators(records: list[Record]) -> None:
    assert sum(r.kappa > 1.0e4 for r in records) >= 8
    assert sum(r.kappa < 100.0 for r in records) >= 4


def test_the_residual_alone_stops_on_a_wrong_value_in_many_operators(
    records: list[Record],
) -> None:
    # Without this, the test would prove nothing: the residual would be enough.
    worse = sum(r.worse("residual") for r in records)
    assert worse >= 0.3 * len(records)


def test_the_stop_after_rest_is_as_accurate_as_the_whole_basis(
    records: list[Record],
) -> None:
    for r in records:
        assert r.stops["rest5"].error <= 10.0 * r.full.error + 1.0e-5, r


def test_a_longer_rest_never_stops_earlier(records: list[Record]) -> None:
    for r in records:
        steps = [r.stops[rule].step for rule in RULES]
        assert steps == sorted(steps), r


def test_waiting_costs_steps_but_not_the_whole_basis(records: list[Record]) -> None:
    mean_rest = np.mean([r.stops["rest5"].step for r in records])
    mean_residual = np.mean([r.stops["residual"].step for r in records])
    mean_full = np.mean([r.full.step for r in records])
    assert mean_rest - mean_residual <= 25.0
    assert mean_rest < 0.5 * mean_full


def main(argv: list[str]) -> None:
    trials = int(argv[1])
    seeds = [int(s) for s in argv[2:]] or [2026]
    with tempfile.TemporaryDirectory() as tmp:
        exe = build_driver(Path(tmp))
        records: list[Record] = []
        for seed in seeds:
            records += sweep(exe, trials, seed)
    print(table(records))


if __name__ == "__main__":
    main(sys.argv)
