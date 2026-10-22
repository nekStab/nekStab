"""High-level FST inflow driver: wavenumbers -> OSS modes -> one consolidated file.

In-memory replacement for the legacy two-stage, hundreds-of-tiny-files MATLAB
workflow.  Every mode is generated in memory (optionally in parallel) and written
to ONE consolidated file that the Fortran reads with a single open + one loop.

Two output formats (choose with `fmt`):

  fmt="dat"  -- formatted ASCII, human-checkable.  The shared y-grid is written
                ONCE (it is identical for every mode), then per mode:
                    # --- mode m ---
                    <omega> <gamma> <beta>
                    Ny rows of:  Re(u) Im(u)  Re(v) Im(v)  Re(w) Im(w)
                Header line:  <n_modes> <Ny> <Re>;  '#' lines are comments.

  fmt="bin"  -- little-endian stream binary (read in Fortran with access='stream'):
                    8 bytes magic 'NKSTFST2'
                    int32 numk, int32 nmodes, int32 Ny
                    float64 Re, okini, okfin, length, tu
                    float64 y[Ny]
                              (columns Re(u) Im(u) Re(v) Im(v) Re(w) Im(w))
                Much faster to parse than ASCII for large n_modes.

y rows are ascending (wall -> free stream); a reader may flip if it prefers the
opposite orientation.
"""

import struct
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .wavenumber import generate_wavenumbers
from .oss import solve_oss_mode

_MAGIC = b"NKSTFST2"  # 8-byte binary-format identifier


def _solve_one(task):
    """Module-level (picklable) worker for one OSS mode.

    task = (index, omega, gamma, beta, Re, Ny, Ly, seed_seq).  Each mode gets its
    own SeedSequence so results are reproducible AND independent of execution
    order -- safe under parallel evaluation.
    """
    index, omega, gamma, beta, Re, Ny, Ly, seed_seq = task
    rng = np.random.default_rng(seed_seq)
    y, U, V, W, diag = solve_oss_mode(omega, gamma, beta, Re, Ny, Ly, rng)
    return index, {"omega": float(omega), "gamma": float(gamma),
                   "beta": float(beta), "U": U, "V": V, "W": W, "y": y,
                   "newton_res": diag["newton_res"],
                   "spectral_div": diag["spectral_div"]}


def generate_fst_inflow(Re, Ny, Ly, numk, kmin, kmax, seed=None, jobs=1,
                        verbose=True, tu=0.010, length=None):
    """Generate all FST modes in memory (serially if jobs<=1, else in parallel).


    Returns a dict with keys: wavenumbers, y, modes, Re, Ny, Ly, numk, nmodes,
    okini, okfin, length, tu, max_newton_res, max_spectral_div.
    Each mode dict has omega, gamma, beta, U, V, W (complex (Ny,)), y, newton_res,
    spectral_div.
    """
    wavenumbers = generate_wavenumbers(numk, kmin, kmax, seed=seed)
    n_modes = wavenumbers.shape[0]

    # Independent, reproducible per-mode RNG streams.  spawn() guarantees the
    # streams are statistically independent and tied to `seed`, regardless of the
    # order in which the parallel workers finish.
    child_seeds = np.random.SeedSequence(seed).spawn(n_modes)
    tasks = [(i, float(w[0]), float(w[1]), float(w[2]), Re, Ny, Ly, child_seeds[i])
             for i, w in enumerate(wavenumbers)]

    results = [None] * n_modes
    if jobs and jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for done, (index, mode) in enumerate(ex.map(_solve_one, tasks), 1):
                results[index] = mode
                if verbose and (done % 10 == 0 or done == n_modes):
                    print(f"  {done}/{n_modes} modes done (parallel, jobs={jobs})")
    else:
        for task in tasks:
            index, mode = _solve_one(task)
            results[index] = mode
            if verbose and (index % 10 == 0 or index == n_modes - 1):
                print(f"  mode {index + 1}/{n_modes}: "
                      f"(w,g,b)=({mode['omega']:.3f},{mode['gamma']:.3f},"
                      f"{mode['beta']:.3f}) newton_res={mode['newton_res']:.1e}")

    y_shared = results[0]["y"]
    max_res = max(m["newton_res"] for m in results)
    max_div = max(m["spectral_div"] for m in results)
    
    # Handle length default
    if length is None:
        length = 1.80 / kmax
    
    if verbose:
        print(f"Generated {n_modes} modes; worst Newton res={max_res:.2e}, "
              f"worst spectral div={max_div:.2e}")

    return {"wavenumbers": wavenumbers, "y": y_shared, "modes": results,
            "Re": Re, "Ny": Ny, "Ly": Ly,
            "numk": int(numk), "nmodes": int(n_modes // numk),
            "okini": float(kmin), "okfin": float(kmax),
            "length": float(length), "tu": float(tu),
            "max_newton_res": max_res, "max_spectral_div": max_div}

def write_inflow_file(path, result, fmt="dat"):
    """Write the in-memory FST result to one consolidated file ('dat' or 'bin')."""
    # Contract: the header stores numk and nmodes, and the Fortran reader loops
    # over exactly numk*nmodes modes. Guard against any drift between that product
    # and the actual mode list, which would silently misalign the file stream.
    numk, nmodes, n_modes = result["numk"], result["nmodes"], len(result["modes"])
    if numk * nmodes != n_modes:
        raise ValueError(
            f"FST mode-count mismatch: numk*nmodes={numk}*{nmodes}={numk * nmodes} "
            f"!= {n_modes} generated modes; the consolidated file requires exactly "
            f"numk*nmodes modes.")
    if fmt == "bin":
        return _write_binary(path, result)
    if fmt == "dat":
        return _write_formatted(path, result)
    raise ValueError(f"unknown fmt {fmt!r} (use 'dat' or 'bin')")


def _write_formatted(path, result):
    modes = result["modes"]
    y = result["y"]
    n_modes = len(modes)
    Ny = len(y)
    Re = result["Re"]
    numk = result["numk"]
    nmodes = result["nmodes"]
    okini = result["okini"]
    okfin = result["okfin"]
    length = result["length"]
    tu = result["tu"]
    with open(path, "w") as f:
        f.write(f"# nekStab FST inflow -- {n_modes} modes, Ny={Ny}, Re={Re:g}. "
                f"Shared y-grid written ONCE below; per mode: 'omega gamma beta' "
                f"then Ny rows of [Re(u) Im(u) Re(v) Im(v) Re(w) Im(w)].\n")
        f.write(f"{numk:d} {nmodes:d} {Ny:d} {Re:.8e} {okini:.8e} {okfin:.8e} {length:.8e} {tu:.8e}\n")
        f.write("# y-grid (Ny values, ascending: wall -> free stream)\n")
        f.write(" ".join(f"{yy:.8e}" for yy in y) + "\n")
        for m, mode in enumerate(modes, start=1):
            f.write(f"# --- mode {m} ---\n")
            f.write(f"{mode['omega']:.8e} {mode['gamma']:.8e} {mode['beta']:.8e}\n")
            U, V, W = mode["U"], mode["V"], mode["W"]
            for j in range(Ny):
                f.write(f"{U[j].real:.6e} {U[j].imag:.6e} "
                        f"{V[j].real:.6e} {V[j].imag:.6e} "
                        f"{W[j].real:.6e} {W[j].imag:.6e}\n")
    return path


def _write_binary(path, result):
    modes = result["modes"]
    y = np.asarray(result["y"], dtype="<f8")
    n_modes = len(modes)
    Ny = int(len(y))
    Re = float(result["Re"])
    numk = int(result["numk"])
    nmodes = int(result["nmodes"])
    okini = float(result["okini"])
    okfin = float(result["okfin"])
    length = float(result["length"])
    tu = float(result["tu"])
    with open(path, "wb") as f:
        f.write(_MAGIC)
        f.write(struct.pack("<iii", numk, nmodes, Ny))
        f.write(struct.pack("<ddddd", Re, okini, okfin, length, tu))
        y.tofile(f)
        for mode in modes:
            f.write(struct.pack("<ddd", mode["omega"], mode["gamma"], mode["beta"]))
            U, V, W = mode["U"], mode["V"], mode["W"]
            block = np.empty((Ny, 6), dtype="<f8")
            block[:, 0] = U.real; block[:, 1] = U.imag
            block[:, 2] = V.real; block[:, 3] = V.imag
            block[:, 4] = W.real; block[:, 5] = W.imag
            block.tofile(f)
    return path


def read_inflow_binary(path):
    """Read back a binary inflow file (for round-trip verification / Python use).

    Returns (numk, nmodes, Ny, Re, okini, okfin, length, tu, y, modes) where modes is a list of
    (omega, gamma, beta, block[Ny,6]).
    """
    with open(path, "rb") as f:
        if f.read(8) != _MAGIC:
            raise ValueError("not a nekStab FST binary file (bad magic)")
        numk, nmodes, Ny = struct.unpack("<iii", f.read(12))
        Re, okini, okfin, length, tu = struct.unpack("<ddddd", f.read(40))
        y = np.fromfile(f, dtype="<f8", count=Ny)
        modes = []
        for _ in range(numk * nmodes):
            om, ga, be = struct.unpack("<ddd", f.read(24))
            block = np.fromfile(f, dtype="<f8", count=Ny * 6).reshape(Ny, 6)
            modes.append((om, ga, be, block))
    return numk, nmodes, Ny, Re, okini, okfin, length, tu, y, modes
