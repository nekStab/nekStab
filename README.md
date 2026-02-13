```
                 __   _____  __          __
   ____   ___   / /__/ ___/ / /_ ____ _ / /_
  / __ \ / _ \ / //_/\__ \ / __// __ `// __ \
 / / / //  __// ,<  ___/ // /_ / /_/ // /_/ /
/_/ /_/ \___//_/|_|/____/ \__/ \__,_//_.___/
```

| OS | Arch | Compiler | MPI | Status |
|:---|:----:|:---------|:----|:------:|
| Ubuntu 24.04 | x86_64 | gfortran 14 | OpenMPI / MPICH | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=dev&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |
| Ubuntu 24.04 | x86_64 | ifort 2024 / ifx 2025 | Intel MPI | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=dev&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |
| macOS 26 | ARM64 | gfortran 14 | OpenMPI | [![CI](https://img.shields.io/github/actions/workflow/status/nekStab/nekStab/ci.yml?branch=dev&label=)](https://github.com/nekStab/nekStab/actions/workflows/ci.yml) |

**nekStab** is a toolbox for global stability and bifurcation analysis using the spectral element solver [Nek5000](https://github.com/Nek5000/Nek5000). Released under BSD-3-Clause license.

## Features

### Steady-State & Periodic Orbit Computation
| Method | Description |
|--------|-------------|
| **SFD** | Selective Frequency Damping for unstable steady states |
| **BoostConv** | Residual acceleration for slow convergence |
| **TDF** | Time-Delayed Feedback for periodic orbits |
| **Newton-Krylov** | Quadratic convergence for fixed points and UPOs |

### Global Stability Analysis
| Analysis | Steady Flows | Time-Periodic (Floquet) |
|----------|:------------:|:-----------------------:|
| **Direct eigenmodes** | ✓ | ✓ |
| **Adjoint eigenmodes** | ✓ | ✓ |
| **Transient growth** | ✓ | ✓ |

### Sensitivity & Receptivity
- **Wavemaker** — structural sensitivity to feedback
- **Base flow sensitivity** — response to mean flow modifications
- **Forcing response** — optimal and localized forcing analysis

### Modal Decomposition
| Method | Description |
|--------|-------------|
| **POD** | Proper Orthogonal Decomposition (energy-ranked modes) |
| **DMD** | Dynamic Mode Decomposition (frequency-ranked modes) |
| **SPOD** | Spectral POD (frequency-resolved coherent structures) |

### Advanced Capabilities
- **OTD modes** — real-time Lyapunov vectors for chaotic flows
- **Linearized DNS** — perturbation evolution around base flows
- **Scalar transport** — temperature and passive scalar stability
- **Sponge zones** — non-reflecting boundaries for open flows
- **Vortex identification** — λ₂, Q-criterion, Ω-criterion output
- **FFT module** — FFTW3/MKL integration for spectral analysis with Nek5000

### Why nekStab?
- **Matrix-free**: No Jacobian storage — scales to millions of DoFs
- **Krylov-based**: Time-stepper as linear operator — eigenvalues from snapshots, no math required
- **Spectral accuracy**: Leverages Nek5000's high-order elements
- **MPI parallel**: Efficient on laptops to supercomputers
- **Validated**: Benchmarked against canonical flows (cylinder, cavity, jets)
- **FFT included**: Smart linking to FFTW3 (GCC) or MKL (Intel) — no manual setup

## Quick Reference

| Mode | String | Description |
|:----:|--------|-------------|
| 0 | `'dns'` | Direct Numerical Simulation |
| 0.1 | `'linear_dns'` | Linearized DNS (perturbation) |
| 1.1 | `'sfd'` | Selective Frequency Damping |
| 1.2 | `'boostconv'` | BoostConv acceleration |
| 1.4 | `'tdf'` | Time-Delayed Feedback |
| 2.0 | `'newton_fp'` | Newton for fixed points |
| 2.1 | `'newton_po'` | Newton for periodic orbits |
| 3.1 | `'direct'` | Direct stability eigenmodes |
| 3.11 | `'floquet_direct'` | Floquet direct analysis |
| 3.2 | `'adjoint'` | Adjoint stability eigenmodes |
| 3.21 | `'floquet_adjoint'` | Floquet adjoint analysis |
| 3.3 | `'transient_growth'` | Optimal perturbations |
| 3.31 | `'floquet_tg'` | Floquet transient growth |
| 4.1 | `'energy_budget'` | Kinetic energy budget |
| 4.11 | `'energy_budget_floquet'` | Floquet energy budget |
| 4.2 | `'wavemaker'` | Structural sensitivity |
| 4.3 | `'bf_sensitivity'` | Base flow sensitivity |
| 5 | `'otd'` | Optimally Time-Dependent modes |
| 6.1 | `'pod'` | Proper Orthogonal Decomposition |
| 6.2 | `'dmd'` | Dynamic Mode Decomposition |
| 6.3 | `'spod'` | Spectral POD |

See [DOC.md](DOC.md) for full parameter reference (includes animation and forcing modes).

## Examples

Ready-to-run cases in `example/`:

| Case | Physics | Demonstrated Features |
|------|---------|----------------------|
| `cylinder/` | 2D/3D wake | DNS, Newton, stability, Floquet, wavemaker, OTD, POD/DMD/SPOD |
| `back_fstep/` | Separation | Convective instability, transient growth |
| `lid_driven/` | Confined | Steady bifurcations |
| `blasius/` | Boundary layer | Tollmien-Schlichting waves |
| `flip_flop/` | Side-by-side cylinders | Neimark-Sacker, Floquet |
| `thersyphon/` | Buoyancy-driven | Pitchfork + Hopf bifurcations |
| `torus/` | Curved pipe | Dean instability, 3D modes |
| `tpjet/` | Forced jet | Floquet period-doubling |
| `slot_FST/` | Free-stream turbulence | Synthetic inflow generation |

## First Steps

**Prerequisites**

Linux (GCC)
```bash
sudo apt -y install build-essential gfortran libmpich-dev libopenblas-dev libfftw3-dev cmake m4 htop
```

macOS
```bash
brew install mpich gfortran fftw wget git cmake htop
```

**Installation**

```bash
git clone --depth=1 https://github.com/nekStab/nekStab.git
cd nekStab
./Nek5000setup.sh
```

| Command | What it does |
|---------|--------------|
| `git clone --depth=1 ...` | Downloads only the latest version (fastest, minimal download) |
| `cd nekStab` | Enter the nekStab directory |
| `./Nek5000setup.sh` | Downloads and configures Nek5000 inside `nekStab/Nek5000/` |

> **Already have Nek5000?** Skip `./Nek5000setup.sh` and point `NEK_SOURCE_ROOT` to your existing installation (see below).

Add to your shell config (`~/.bashrc` or `~/.zshrc`):
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
cd ~/nekStab/example/cylinder/baseflow/newton
mks 1cyl
```

### Compiler Selection

nekStab supports multiple Fortran compilers. The build script auto-detects available compilers in this order: `ifx` → `ifort` → `gfortran`. To force a specific compiler:

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
nekbmpi 1cyl 4       # Run case "1cyl" on 4 MPI processes
tail -f logfile      # Monitor output in real-time
killall nek5000      # Stop the simulation
```

| Command | What it does |
|---------|--------------|
| `nekbmpi 1cyl 4` | Launches `nek5000` via `mpirun` with 4 processes (adjust to your CPU cores) |
| `tail -f logfile` | Streams solver output — watch convergence, time steps, diagnostics |
| `killall nek5000` | Gracefully stops all running Nek5000 processes |

> **Tip:** Use `nekbmpi 1cyl 4 &` to run in background, then `tail -f logfile` in the same terminal.

For more information, see the [Documentation](DOC.md).

## Mode Selection

Three equivalent ways to select operating mode:

| Method | Where | Example |
|--------|-------|---------|
| String | `.usr` | `nekstab_mode = 'floquet_adjoint'` |
| Flags | `.usr` | `isAdjoint = .true.` + `ifFloquet = .true.` |
| uparam | `.par` | `userParam01 = 3.21` |

## Development

**nekStab** is maintained by [Ricardo Frantz](https://github.com/ricardofrantz).

> **See also:** [LightKrylov](https://github.com/nekStab/LightKrylov) is a modern Fortran library providing abstract Krylov methods. [neklab](https://github.com/nekStab/neklab) is the next-generation bifurcation and stability analysis toolbox for Nek5000, built on LightKrylov. Both are actively developed by [Jean-Christophe Loiseau](https://github.com/loiseaujc) and [Simon Kern](https://github.com/Simkern).

Website: [https://nekstab.github.io/](https://nekstab.github.io/)

## Citation

When using **nekStab**, please cite [Frantz et al. (2023)](https://doi.org/10.1115/1.4056808). See [CITATIONS.md](CITATIONS.md) for BibTeX entries and additional references.
