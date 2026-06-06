"""nekStab Python tools.

A small package of pre/post-processing utilities for the nekStab spectral-element
stability solver.  Currently provides:

  nekstab.fst -- Free-Stream-Turbulence inflow generation (continuous-spectrum
                 Orr-Sommerfeld-Squire modes; Brandt et al. 2004).

Run the FST generator as a CLI:
    PYTHONPATH=python python -m nekstab.fst --help
or install editable and run from anywhere:
    uv pip install -e python/
    python -m nekstab.fst --help
"""

__version__ = "0.1.0"
