# `validation/figures/_stale/` — historical artifacts

Figures that were once referenced by the gallery but are no longer
catalog-referenced. They are preserved as historical evidence (so
nothing is lost) but are excluded from active gallery scanning.

## Why preserved (not deleted)

- The figure was the visual evidence for a prior operating point,
  parameter choice, or method that we've since revised. Deleting would
  lose the historical record of how the case used to look.
- The figure may be re-referenced if/when the corresponding case is
  re-introduced (e.g., a Re=2500 NACA 0012 run is regenerated).

## Convention

- `validation/build_gallery.py` does NOT glob into this directory; only
  `validation/figures/*.png` (top-level) is scanned.
- Catalog `Artifact` paths must not point at files under `_stale/`.
  If a stale file becomes relevant again, move it back to the top
  level and re-add the catalog reference.

## Current contents

- **`naca0012_Re2500.png`** — second operating point of the NACA 0012
  case. Catalog now targets Re=2000 only; no matching `.par` exists
  for Re=2500. Moved here 2026-05-21 when the family was renamed to
  `naca0012_re2000/`.
