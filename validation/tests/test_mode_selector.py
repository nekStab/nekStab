"""Mode selector compatibility: strings, flags and legacy userParam01 decimals.

src/mode_config.f90 is compiled with a stand-in for the Nek5000 bridge
(fortran/mode_selector_stub.f90) and run once per input by
fortran/mode_selector_driver.f90. The expected flags below are written by hand,
not read from the source, so a change to the decoder shows up as a failure.
Needs gfortran; the tests skip when it is missing.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FORTRAN = Path(__file__).resolve().parent / "fortran"

# userParam01 -> flags that must be true after resolution, and uparam(1) written back.
DECIMALS = {
    0: ({"ifDNS"}, 0.0),
    0.1: ({"ifLinDNS"}, 0.1),
    1.1: ({"ifSFD"}, 1.1),
    1.2: ({"ifBoostConv"}, 1.2),
    1.3: ({"ifDMT"}, 1.3),
    1.4: ({"ifTDF"}, 1.4),
    2.0: ({"isNewtonFP"}, 2.0),
    2.1: ({"isNewtonPO"}, 2.1),
    2.2: ({"isNewtonPO_T"}, 2.2),
    3.1: ({"isDirect"}, 3.1),
    3.11: ({"isFloquetDirect"}, 3.11),
    3.2: ({"isAdjoint"}, 3.2),
    3.21: ({"isFloquetAdjoint"}, 3.21),
    3.3: ({"isTransientGrowth"}, 3.3),
    3.31: ({"isFloquetTransientGrowth"}, 3.31),
    4: ({"ifEnergyBudget", "ifWavemaker", "ifBFSensitivity"}, 4.1),
    4.1: ({"ifEnergyBudget"}, 4.1),
    4.11: ({"ifEnergyBudget", "ifFloquet"}, 4.11),
    4.2: ({"ifWavemaker"}, 4.2),
    4.3: ({"ifBFSensitivity"}, 4.3),
    4.41: ({"ifForceSensReal"}, 4.41),
    4.42: ({"ifForceSensImag"}, 4.42),
    4.43: ({"ifDeltaForcing"}, 4.43),
    4.5: ({"ifAnimateMode"}, 4.5),
    4.51: ({"ifAnimateBFDeform"}, 4.51),
    4.52: ({"ifAnimateFloquet"}, 4.52),
    5: ({"ifotd"}, 5.0),
    6: ({"ifpod", "ifdmd", "ifspod"}, 6.1),
    6.1: ({"ifpod"}, 6.1),
    6.2: ({"ifdmd"}, 6.2),
    6.3: ({"ifspod"}, 6.3),
}

# nekstab_mode string -> flags that must be true (uparam(1) = 0 in the .par).
STRINGS = {
    "dns": {"ifDNS"},
    "DNS": {"ifDNS"},
    "linear_dns": {"ifLinDNS"},
    "sfd": {"ifSFD"},
    "boostconv": {"ifBoostConv"},
    "tdf": {"ifTDF"},
    "dmt": {"ifDMT"},
    "newton": {"isNewtonFP"},
    "newton_fp": {"isNewtonFP"},
    "upo": {"isNewtonPO"},
    "newton_po_t": {"isNewtonPO_T"},
    "direct": {"isDirect"},
    "floquet_direct": {"isFloquetDirect"},
    "adjoint": {"isAdjoint"},
    "floquet_adjoint": {"isFloquetAdjoint"},
    "transient_growth": {"isTransientGrowth"},
    "floquet_tg": {"isFloquetTransientGrowth"},
    "energy_budget": {"ifEnergyBudget"},
    "energy_budget_floquet": {"ifEnergyBudget", "ifFloquet"},
    "wavemaker": {"ifWavemaker"},
    "bf_sensitivity": {"ifBFSensitivity"},
    "force_sensitivity_real": {"ifForceSensReal"},
    "force_sensitivity_imag": {"ifForceSensImag"},
    "delta_forcing": {"ifDeltaForcing"},
    "animate_mode": {"ifAnimateMode"},
    "animate_bf_deform": {"ifAnimateBFDeform"},
    "animate_floquet": {"ifAnimateFloquet"},
    "otd": {"ifotd"},
    "pod": {"ifpod"},
    "dmd": {"ifdmd"},
    "spod": {"ifspod"},
}


@pytest.fixture(scope="module")
def driver(tmp_path_factory) -> Path:
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran not found")
    build = tmp_path_factory.mktemp("mode_selector")
    flags = ["-fdefault-real-8", "-fdefault-double-8", "-ffree-line-length-none", "-O0"]
    sources = [
        ROOT / "src" / "mode_codes.f90",
        FORTRAN / "mode_selector_stub.f90",
        ROOT / "src" / "mode_config.f90",
        FORTRAN / "mode_selector_driver.f90",
    ]
    for src in sources:
        subprocess.run([gfortran, *flags, "-c", str(src)], cwd=build, check=True)
    exe = build / "driver"
    subprocess.run(
        [gfortran, *flags, "-o", str(exe), *(f"{s.stem}.o" for s in sources)],
        cwd=build,
        check=True,
    )
    return exe


def resolve(exe: Path, up1: float, *args: str) -> tuple[int, set[str], float | None, str]:
    """Run the driver; return (exit status, true flags, uparam(1) after, output)."""
    run = subprocess.run([str(exe), repr(up1), *args], capture_output=True, text=True)
    out = run.stdout + run.stderr
    match = re.search(r"^(.*)\| uparam1=\s*(\S+)", out, re.M)
    if run.returncode != 0 or match is None:
        return run.returncode, set(), None, out
    return run.returncode, set(match.group(1).split()), float(match.group(2)), out


@pytest.mark.parametrize("up1", sorted(DECIMALS))
def test_legacy_decimal(driver, up1):
    flags, written_back = DECIMALS[up1]
    status, got, after, out = resolve(driver, up1)
    assert status == 0, out
    assert got == flags
    assert after == pytest.approx(written_back, abs=1e-6)


@pytest.mark.parametrize("name", sorted(STRINGS))
def test_string_selects_mode(driver, name):
    status, got, _, out = resolve(driver, 0, name)
    assert status == 0, out
    assert got == STRINGS[name]


def test_string_overrides_userparam(driver):
    status, got, after, out = resolve(driver, 3.1, "adjoint")
    assert status == 0, out
    assert got == {"isAdjoint"}
    assert after == pytest.approx(3.2, abs=1e-6)


def test_string_overrides_flags(driver):
    status, got, _, out = resolve(driver, 0, "sfd", "isDirect")
    assert status == 0, out
    assert got == {"ifSFD"}


def test_flags_select_mode(driver):
    status, got, after, out = resolve(driver, 0, "-", "isNewtonFP")
    assert status == 0, out
    assert got == {"isNewtonFP"}
    assert after == pytest.approx(2.0, abs=1e-6)


def test_floquet_flag_modifies_stability_flag(driver):
    status, got, after, out = resolve(driver, 0, "-", "isAdjoint", "ifFloquet")
    assert status == 0, out
    assert got == {"ifFloquet", "isFloquetAdjoint"}
    assert after == pytest.approx(3.21, abs=1e-6)


@pytest.mark.parametrize("up1", [3.15, 3.104, 7.3, 0.2, 1.0, 2.5, 4.4, 6.4, -1.0, 12.0])
def test_unknown_userparam_stops_the_run(driver, up1):
    # An unknown code must stop the run. It used to run a DNS with a warning.
    status, got, _, out = resolve(driver, up1)
    assert status == 3, out
    assert "ERROR: Unknown mode" in out
    assert not got


def test_unknown_string_stops_the_run(driver):
    status, _, _, out = resolve(driver, 0, "direct_typo")
    assert status == 3, out
    assert "Unknown nekstab_mode" in out


@pytest.mark.parametrize(
    "args",
    [("-", "isDirect", "isAdjoint"), ("-", "isDirect", "ifSFD"), ("-", "isDirect", "isTransientGrowth")],
)
def test_conflicting_flags_stop_the_run(driver, args):
    status, _, _, out = resolve(driver, 0, *args)
    assert status == 3, out
    assert "ERROR" in out


def test_orphan_floquet_flag_stops_the_run(driver):
    status, _, _, out = resolve(driver, 0, "-", "ifFloquet")
    assert status == 3, out
    assert "ifFloquet" in out


def test_every_example_userparam01_is_a_mode(driver):
    # Every numeric userParam01 in a shipped .par file must decode.
    seen = set()
    for par in sorted((ROOT / "example").rglob("*.par")):
        for line in par.read_text(errors="replace").splitlines():
            match = re.match(r"\s*userParam01\s*=\s*([-+0-9.eE]+)\s*(?:#.*)?$", line, re.I)
            if match:
                seen.add(float(match.group(1)))
    assert seen, "no userParam01 found in example/**/*.par"
    for up1 in sorted(seen):
        status, got, _, out = resolve(driver, up1)
        assert status == 0 and got, f"userParam01 = {up1}: {out}"
