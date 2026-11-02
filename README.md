![nekStab banner](assets/nekstab-readme-banner-v4.jpg)

| OS | Arch | Compiler | MPI | Status |
|:---|:----:|:---------|:----|:------:|
| Ubuntu (latest runner) | x86_64 | gfortran 14 | OpenMPI / MPICH | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=main&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |
| Ubuntu (latest runner) | x86_64 | ifx 2025.2 | Intel MPI | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=main&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |
| macOS 26 | ARM64 | gfortran 14 | OpenMPI | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=main&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |

**nekStab** is a toolbox for global stability and bifurcation analysis using the spectral element solver [Nek5000](https://github.com/Nek5000/Nek5000). Released under BSD-3-Clause license.

## Features

### Steady-State & Periodic Orbit Computation
| Method | Description |
|--------|-------------|
| **SFD** | Selective Frequency Damping for unstable steady states |
| **BoostConv** | Residual acceleration for slow convergence |
| **TDF** | Time-Delayed Feedback for periodic orbits |
| **Newton-Krylov** | Quadratic convergence for fixed points, UPOs (natural frequency), and forced periodic orbits |

### Global Stability Analysis
| Analysis | Steady Flows | Time-Periodic (Floquet) |
|----------|:------------:|:-----------------------:|
| **Direct eigenmodes** | ✓ | ✓ |
| **Adjoint eigenmodes** | ✓ | ✓ |
| **Transient growth** | ✓ | ✓ |

### Sensitivity & Receptivity
- **Wavemaker** — structural sensitivity to feedback
- **Energy budget** — production, dissipation, transport decomposition (steady + Floquet)
- **Base flow sensitivity** — response to mean flow modifications
- **Steady force sensitivity** — optimal placement for passive control (imaginary part, 4.42, not validated in 2.0.0)
- **Mode animation** — reconstruct eigenmode dynamics over one period

### Modal Decomposition
| Method | Description |
|--------|-------------|
| **POD** | Proper Orthogonal Decomposition (energy-ranked modes) |
| **DMD** | Dynamic Mode Decomposition (frequency-ranked modes) |
| **SPOD** | Spectral POD (frequency-resolved coherent structures) |
| **DMT filter bank** | Online second-order filters on the DNS: the spatial structure of the flow at chosen frequencies |

### Advanced Capabilities
- **OTD modes** — real-time Lyapunov vectors for chaotic flows
- **Linearized DNS** — perturbation evolution around base flows
- **Scalar transport** — temperature and passive scalar stability
- **Sponge zones** — non-reflecting boundaries for open flows
- **Vortex identification** — λ₂, Q-criterion, Ω-criterion output
- **FFT module** — FFTW3/MKL integration for spectral analysis with Nek5000

### Why nekStab?
- **Matrix-free**: No Jacobian storage — scales to millions of DoFs
- **Krylov-based**: Time-stepper as linear operator — eigenvalues from snapshots
- **Spectral accuracy**: Leverages Nek5000's high-order elements
- **MPI parallel**: Efficient on laptops to supercomputers
- **Validated**: 80 example stages in 15 flow families, each with reference results (see [`EXAMPLES.md`](EXAMPLES.md)). CI builds one case, the cylinder DNS `example/cylinder_re100/000_dns/ci_test`, and runs it for 10 time steps as a smoke test; it does not run the stages.
- **FFT included**: Smart linking to FFTW3 (GCC) or MKL (Intel) — no manual setup

## Quick Reference

| Mode | String | Aliases | Description |
|:----:|--------|---------|-------------|
| 0 | `'dns'` | | Direct Numerical Simulation |
| 0.1 | `'linear_dns'` | `'lindns'` | Linearized DNS (perturbation) |
| 1.1 | `'sfd'` | | Selective Frequency Damping |
| 1.2 | `'boostconv'` | `'boost'` | BoostConv acceleration |
| 1.3 | `'dmt'` | | Dynamic Mode Tracking, base-flow stabilization (Queguineur et al. 2019) |
| 1.4 | `'tdf'` | | Time-Delayed Feedback |
| 2 | `'newton_fp'` | `'newton'` | Newton for fixed points |
| 2.1 | `'newton_po'` | `'upo'` | Newton for periodic orbits (unknown period) |
| 2.2 | `'newton_po_t'` | `'forced_upo'` | Newton for periodic orbits (forced period) |
| 3.1 | `'direct'` | | Direct stability eigenmodes |
| 3.11 | `'floquet_direct'` | | Floquet direct analysis |
| 3.2 | `'adjoint'` | | Adjoint stability eigenmodes |
| 3.21 | `'floquet_adjoint'` | | Floquet adjoint analysis |
| 3.3 | `'transient_growth'` | `'tg'` | Optimal perturbations |
| 3.31 | `'floquet_tg'` | | Floquet transient growth |
| 4.1 | `'energy_budget'` | | Kinetic energy budget |
| 4.11 | `'energy_budget_floquet'` | | Floquet energy budget |
| 4.2 | `'wavemaker'` | | Structural sensitivity |
| 4.3 | `'bf_sensitivity'` | | Base flow sensitivity |
| 4.41 | `'force_sensitivity_real'` | | Steady force sensitivity (real) |
| 4.42 | `'force_sensitivity_imag'` | | Steady force sensitivity (imag), not validated in 2.0.0 |
| 4.43 | `'delta_forcing'` | | Delta forcing response |
| 4.5 | `'animate_mode'` | `'animate'` | Mode animation |
| 4.51 | `'animate_bf_deform'` | | Base flow deformation animation |
| 4.52 | `'animate_floquet'` | | Floquet mode animation |
| 5 | `'otd'` | | Optimally Time-Dependent modes about a frozen base flow |
| 5.1 | `'otd_dns'` | | Optimally Time-Dependent modes about the base flow that the DNS advances |
| 6.1 | `'pod'` | | Proper Orthogonal Decomposition |
| 6.2 | `'dmd'` | | Dynamic Mode Decomposition |
| 6.3 | `'spod'` | | Spectral POD |
| 6.4 | `'dmt_bank'` | | DMT filter bank: filters on the running DNS at chosen frequencies |

See [DOC.md](DOC.md) for full parameter reference.

> **Large domains:** `residualTol` is not a fixed accuracy for stability runs. Nek5000 stops its solves on a residual averaged over the volume `V`, so on a unit-norm mode the error of a solver call grows with `sqrt(V)`: about 5e-7 at `residualTol = 1e-9` on a domain of 300 by 200. The eigensolver therefore tightens `residualTol` to `eigen_tol/(25*sqrt(V))` when that is smaller, and it prints the value at start-up. See [DOC.md](DOC.md), Troubleshooting.

## Examples

Ready-to-run cases live under `example/`, from canonical cylinder workflows to
thermal, airfoil, cavity, channel, and Floquet examples.
See [EXAMPLES.md](EXAMPLES.md) for the stage rationale, case catalog, and
public validation status.

## First Steps

**Prerequisites**

Linux (GCC)
```bash
sudo apt -y install build-essential gfortran libmpich-dev libopenblas-dev libfftw3-dev cmake m4 htop
```

macOS
```bash
brew install open-mpi gfortran fftw wget git cmake htop
```

**Installation**

```bash
git clone --depth=1 --branch v2.0.0 https://github.com/nekStab/nekStab.git
cd nekStab
./Nek5000setup.sh
```

Use `--depth=1`: a full clone also downloads the history of the repository, about 1 GB. To develop, fetch the history later with `git fetch --unshallow`.

The branch `main` is the default branch and holds the releases. Changes reach it by pull request.

The release was validated with Nek5000 commit `d083ef3`. `Nek5000setup.sh` clones the current Nek5000 `master`. The axisymmetric adjoint correction in `src/forcing.f90` is zero on a Nek5000 that already has the radius in its integration weights, so a newer `master` does not change the result of the adjoint runs.

To use the validated commit, run `git -C Nek5000 checkout d083ef3` in the `nekStab` folder (the setup script makes a full clone, so the commit is there), then rebuild the tools with `./maketools genmap genbox` in `Nek5000/tools` and the cases with `mks <case> --fresh`. To return to the default, `git -C Nek5000 checkout master`.

Use git for the install. A source archive from GitHub has no `.git` folder, and a build from it reports the version `unknown`.

| Command | What it does |
|---------|--------------|
| `git clone --depth=1 --branch v2.0.0 ...` | Downloads only the files of the release 2.0.0 (fastest, minimal download) |
| `cd nekStab` | Enter the nekStab directory |
| `./Nek5000setup.sh` | Downloads and configures Nek5000 inside `nekStab/Nek5000/` |

`./Nek5000setup.sh` asks four questions on a fresh clone (a fifth, about removing the folder, appears when `Nek5000/` already exists): install packages (default `n`), clone Nek5000 (default `y`), build `genmap` and `genbox` (default `y`), add the environment variables to your shell file (default `n`). The defaults work: press Enter four times, then load the environment with `source exports.sh` from the `nekStab` folder (`source ~/nekStab/exports.sh` if you cloned into your home directory; see below).

Its package list differs from the prerequisites above: it installs `libopenmpi-dev` (Linux) or `open-mpi` (macOS) and the X11 headers, and not MPICH, OpenBLAS, FFTW or `m4`. Use the prerequisites above (MPICH on Linux, OpenMPI on macOS) and keep the default `n` for the package question. Do not install both MPI libraries: on Linux the script's choice is OpenMPI and this README's is MPICH, and two installs can leave `mpiexec` and the compiler wrapper `mpif90` from different MPI libraries. On macOS install `open-mpi` and not `mpich`: the FFTW formula of Homebrew depends on `open-mpi` and conflicts with `mpich`.

> **Already have Nek5000?** Skip `./Nek5000setup.sh` and point `NEK_SOURCE_ROOT` to your existing installation (see below).

**Simplest way, no file to edit:** load the environment with one line.

```bash
source ~/nekStab/exports.sh
```

The path assumes that the clone is in your home directory; otherwise use the path of your clone.

It sets `NEKSTAB_SOURCE_ROOT` to the folder that holds `exports.sh`, keeps a
`NEK_SOURCE_ROOT` that you set before (an existing Nek5000), puts both `bin/` folders
on `PATH`, removes the stack limit and, in bash, turns on tab completion of `mks` (for zsh see
"Shell completion" below). Type it in
each new terminal, or add that same line once to the end of `~/.bashrc` (or
`~/.zshrc`) to have it always. Use `source`, not `bash exports.sh`: a script that runs
as a program cannot change the terminal that started it.

**Or set the variables yourself** in your shell config (`~/.bashrc` or `~/.zshrc`):
```bash
# nekStab path
export NEKSTAB_SOURCE_ROOT=$HOME/nekStab
export PATH=$NEKSTAB_SOURCE_ROOT/bin:$PATH

# Nek5000 path (adjust if using existing installation)
export NEK_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT/Nek5000  # ← or your existing path
export PATH=$NEK_SOURCE_ROOT/bin:$PATH

# Stack/core limits for large simulations
ulimit -s unlimited
ulimit -c unlimited
```

| Variable / Setting | Purpose |
|--------------------|---------|
| `NEKSTAB_SOURCE_ROOT` | Location of nekStab source (build scripts reference this) |
| `NEK_SOURCE_ROOT` | Location of Nek5000 (required by Nek5000 build system) |
| `PATH` additions | Makes `mks`, `nekbmpi`, `genmap`, etc. available anywhere |
| `ulimit -s unlimited` | Removes stack size limit — prevents crashes in large runs |
| `ulimit -c unlimited` | Enables core dumps for debugging crashes |

Then reload: `source ~/.bashrc` (or restart terminal).

### Shell completion

Enable tab completion for case names and common `mks` targets/options:

```bash
source "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash"
```

For zsh, load Bash completion first:

```bash
autoload -Uz bashcompinit && bashcompinit
source "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash"
```

## Compilation

Go to a given example folder and compile the code:
```bash
cd ~/nekStab/example/cylinder_re100/210_baseflow_newton/fp
mks 1cyl
```

**History points.** A case that calls `hpts` needs the history file `<case>.his`, which holds the probe list: the first line is the number of probes, then one line of coordinates (`x y z`) per probe. The file is gitignored because a run fills it with data. Each example stage that uses history points tracks the list as `<case>.probes` (`1cyl.probes` for the case `1cyl`; a few stages track a ready `<case>.his` instead), and `mks` copies a `.probes` file to `<case>.his` when that file does not exist. `mks` never overwrites an existing `.his`. To change the probes, edit `<case>.probes`, delete `<case>.his` and run `mks` again. If the `.usr` calls `hpts` and neither file exists, `mks` prints a warning that says what to write.

### Compiler Selection

nekStab supports multiple Fortran compilers. The build script auto-detects available compilers in this order: `ifx` → `ifort` → `gfortran`. The Intel compilers are on the path only after `source /opt/intel/oneapi/setvars.sh` (the default install folder of oneAPI). When they are on the path, `mks` picks `ifx` over `gfortran` unless `NEKSTAB_FC` is set. To force a specific compiler:

```bash
NEKSTAB_FC=ifx mks 1cyl    # Intel LLVM (recommended for Intel/AMD CPUs)
NEKSTAB_FC=ifort mks 1cyl  # Intel Classic (available on many HPC systems)
NEKSTAB_FC=gcc mks 1cyl    # GCC/gfortran
```

### Performance Considerations

| Compiler | Performance | Notes |
|----------|-------------|-------|
| **ifx** (Intel LLVM) | Fastest | Requires Intel oneAPI 2024+, uses MKL |
| **ifort** (Intel Classic) | Fast | Available on HPC systems, discontinued in oneAPI 2025 |
| **gfortran** (GCC) | Good | Universal, uses system BLAS/LAPACK |

gfortran 15.2 and ifx 2026.2 were tested and pass `cylinder_re100/310_stability_direct` (leading growth rate 0.1247875).

> **Note:** Intel discontinued `ifort` in oneAPI 2025, but it remains available on many HPC systems. New installations should use `ifx`.

The Intel compilers use MKL (Math Kernel Library) which provides highly optimized BLAS/LAPACK routines with runtime CPU dispatching.

### Debug Mode

For debugging crashes or numerical issues:
```bash
mks 1cyl --debug
```

| What it enables | Why it helps |
|-----------------|--------------|
| Debug symbols (`-g`) | Use `gdb` to inspect variables and step through code |
| Stack traces (`-fbacktrace`) | See exact file:line when a crash occurs |
| Compiler warnings (`-Wall`) | Catch uninitialized variables, type mismatches |

> **Tip:** When a crash happens in debug mode, the output shows which subroutine and line number failed — much easier than hunting through a release build.

## Running

After successful compilation:

```bash
nekbmpi 1cyl 8       # Run case "1cyl" on 8 MPI processes
tail -f logfile      # Monitor output in real-time
killall nek5000      # Stop the simulation
```

| Command | What it does |
|---------|--------------|
| `nekbmpi 1cyl 8` | Writes `SESSION.NAME`, then launches `nek5000` via `mpiexec` with 8 processes |
| `tail -f logfile` | Streams solver output — watch convergence, time steps, diagnostics |
| `killall nek5000` | Gracefully stops all running Nek5000 processes |

> **Tip:** Use `nekbmpi 1cyl 8 &` to run in background, then `tail -f logfile` in the same terminal.

**Rank count.** A run needs at least as many MPI ranks as `lpmin` in the `SIZE` file of the case. The array sizes of `SIZE` follow from `lelt = lelg/lpmin + 3`, so fewer ranks put more elements on one rank than `lelt` allows, and Nek5000 stops at start with a message like `requires lelt >= 532 for np = 4`. Stages differ: most cylinder stages have `lpmin = 8` and `cylinder_re100/412_postproc_steady_force_sensitivity` has `lpmin = 4`, so read `lpmin` in the `SIZE` file of the stage you run.

**Running without `nekbmpi`.** `nek5000` reads the case name and the case directory from the file `SESSION.NAME`. `nekbmpi` and `run.local.slurm` write it. A bare `mpiexec -np 8 ./nek5000` does not, and stops with `Error while reading SESSION.NAME`. Write it first, in the case folder:

```bash
printf "1cyl\n%s/\n" "$(pwd)" > SESSION.NAME    # line 1: case name, line 2: this folder
mpiexec -np 8 ./nek5000
```

A stage whose `.usr` calls `hpts` holds a probe list `<case>.probes`. Its history file `<case>.his` needs no extra step: `mks` writes it from that list when the stage has no `.his` file. A stage that ships a `.his` with data uses that file as the probe list, because `hpts` rewrites the file after the coordinates.

**Stages that need field files from an earlier stage.** Each stage ships its start fields and its main result: see [`example/REFERENCE_RESULTS.md`](example/REFERENCE_RESULTS.md). The other fields are written by a run. The cubic cavity stages ship one start field per case, in the stage that makes it. A stage whose `startFrom` or mode files are missing stops at start with `byte_read() :: fopen failure2!` and `Error reading restart header` (or `Could not open <file>`). The README of the stage lists the files under Prerequisites or Inputs and names the stage that produces them: run that stage first and copy the files. `run.local.slurm` of many stages does these copies, and the two-phase stages (DNS, then Newton or modal analysis) swap the `.par` file between the phases. A bare `mpiexec` does neither. Without Slurm, run the script as a shell script from the stage folder: `SLURM_NTASKS=8 bash run.local.slurm > logfile 2>&1`. `SLURM_NTASKS` is the rank count, at least the `lpmin` of the stage.

**A run overwrites the shipped results.** After a run, `git status` and `git diff --stat` list the files that changed. A field or a figure is binary, so git shows only that it changed. From the repository root, `python3 scripts/check_against_ref.py example/<stage>` compares the run with the quantities in `ref/reference.json`; it does not compare the other fields or the figures. `git restore -- <file>` brings a shipped file back. Name each changed field or figure: a restore of the whole stage folder also discards your own edits to `.par`, `.usr` and `plot.py`.

**Slurm memory.** `run.local.slurm` of the examples asks `--mem-per-cpu=2000` in 65 of 84 scripts, which is 16 GB for 8 tasks; the cubic cavity scripts ask 3000 or 4000. On a machine with less memory the job stays `PENDING` with reason `MaxMemPerLimit`, and `--mem` on the `sbatch` line is refused together with `--mem-per-cpu`. Lower the request on the `sbatch` line: `sbatch --mem-per-cpu=500 run.local.slurm` (500 worked for `cylinder_re100/310_stability_direct`). Check the queue with `squeue`.

A DNS that needs more time than one job may run is advanced to a target time through a chain of restarts: see [`auto_restart/README.md`](auto_restart/README.md).

For more information, see the [Documentation](DOC.md).

## Mode Selection

Three equivalent ways to select the operating mode. When more than one is set, the first row wins.

| Method | Where | Example |
|--------|-------|---------|
| String (preferred for new cases) | `.usr` | `nekstab_mode = 'floquet_adjoint'` |
| Flags | `.usr` | `isAdjoint = .true.` + `ifFloquet = .true.` |
| uparam (legacy decimal code, kept for existing cases) | `.par` | `userParam01 = 3.21` |

A `userParam01` that is not one of the codes in the table above, or an unknown
`nekstab_mode` string, stops the run with an error that names the value. Earlier
versions ran a DNS instead.

## Development

**nekStab** is maintained by [Ricardo Frantz](https://github.com/ricardofrantz).

See [FORTRAN_STYLE.md](FORTRAN_STYLE.md) for the current Fortran coding style conventions and guiding principles.

> **See also:** [LightKrylov](https://github.com/nekStab/LightKrylov) is a modern Fortran library providing abstract Krylov methods. [neklab](https://github.com/nekStab/neklab) is the next-generation bifurcation and stability analysis toolbox for Nek5000, built on LightKrylov. Both are actively developed by [Jean-Christophe Loiseau](https://github.com/loiseaujc) and [Simon Kern](https://github.com/Simkern).

Website: [https://nekstab.github.io/](https://nekstab.github.io/)

## Acknowledgments

Several capabilities in nekStab started from work shared by others:

- The free-stream turbulence inflow generation builds on original work by M. A. Bucci.
- The optimally time-dependent (OTD) modes capability grew out of groundwork contributed by [Simon Kern](https://github.com/Simkern) (KTH); the implementation follows Babaee & Sapsis (2016).
- The sponge/fringe damping was first prototyped against the [KTH Toolbox](https://github.com/KTH-Nek5000/KTH_Toolbox) and follows Nordström, Nordin & Henningson (1999) and Lundbladh et al. (1999).

The full list of authors and contributors is in [AUTHORS](AUTHORS).

## Citation

When using **nekStab**, please cite [Frantz et al. (2023)](https://doi.org/10.1115/1.4056808). See [CITATIONS.md](CITATIONS.md) for BibTeX entries and additional references.

## Disclaimer

This software is provided "as is", without warranty of any kind. To the extent
permitted by law, the authors and contributors are not liable for any damage, loss
or claim arising from its use or misuse. You are responsible for how you use it and
for following the laws and rules that apply to you. The full terms are in
[LICENSE](LICENSE).

This is research software. It has not been validated for engineering design,
certification or safety-critical use. Check its results independently before you
rely on them.
