# nekStab `example/_templates/` -- NNN stage templates

Drop-in starting points for new cases. Copy the relevant stage dir into a
geometry case, rename, edit the placeholders marked `# TEMPLATE:`, and set
`viscosity = -<Re>`.

```
cp -r example/_templates/310_stability_direct example/<geom>/310_stability_direct
cd example/<geom>/310_stability_direct
# rename case.par, fill SESSION.NAME, edit Re + sponge
```

---

## NNN stage table

| NNN | Stage | uparam01 | Subroutine | Notes |
|-----|-------|----------|------------|-------|
| 000 | DNS | 0 | (no-op) | cold start or restart |
| 110 | baseflow SFD | 1.1 | `SFD` (fixedp.f90) | variants: akervik, casacuberta, dyn, dyn_oifs |
| 120 | baseflow BoostConv | 1.2 | `BoostConv` (fixedp.f90) | Krylov acceleration |
| 130 | baseflow DMT | 1.3 | `dmt` (src/dmt.f90) | **TO PORT** -- Queguineur 2019 |
| 140 | baseflow TDF | 1.4 | `tdf` (fixedp.f90) | time-delayed feedback |
| 210 | baseflow Newton | 2 / 2.1 / 2.2 | `newton_krylov` | variants: fp, upo, upo_T |
| 310 | stability direct | 3.1 | Krylov-Schur | variants: schur, arnoldi, findiff |
| 311 | stability direct Floquet | 3.11 | (same family, periodic BF) | needs UPO from 210 |
| 320 | stability adjoint | 3.2 | Krylov-Schur adjoint | variants: schur, findiff |
| 321 | stability adjoint Floquet | 3.21 | (same family, periodic BF) | needs UPO from 210 |
| 330 | transient growth | 3.3 | Krylov-Schur (direct+adjoint) | |
| 411 | wavemaker / energy budget | 4.1 / 4.11 | `stability_energy_budget*` | variants: direct, floquet |
| 500 | OTD | 5 | OTD module | co-evolves with DNS |
| 600 | modal POD | postproc keyword | `modal_pod` | |
| 610 | modal DMD | postproc keyword | `modal_dmd` | |
| 620 | modal SPOD | postproc keyword | `modal_spod` | |

> Authoritative dispatch source: grep `uparam(1) ==` in `src/` and `src/mode_config.f90`.

---

## Variant subdirs -- canonical example: `110_baseflow_sfd/`

The `110_baseflow_sfd/` dir shows the variant subdir convention:

```
110_baseflow_sfd/
    case.par             <- base template (base -- variants override)
    akervik/case.par     <- classical SFD (Akervik 2006); Delta from leading eigenvalue
    casacuberta/case.par <- dynamic SFD (Casacuberta 2018); ifdyntol enabled
    dyn/case.par         <- dyn-tol SFD (userParam08 factor, userParam09 stride)
    dyn_oifs/case.par    <- dyn-tol + OIFS (targetCFL~5, extrapolation=OIFS)
```

Rule: physics (mesh, .usr) DRY at NNN level; strategy (parameters) explicit per subdir.

---

## Reserved userParam slots

| Slot | Convention |
|------|------------|
| 01 | mode dispatch |
| 02 | OTD mode count |
| 03 | schur_tgt (stability) |
| 04 | SFD Delta |
| 05 | SFD chi OR forcing St_D (TDF/forced cases) |
| 06 | Rayleigh / isothermal flag / jet amplitude (case-specific) |
| 07 | k_dim -- universal across all Krylov-using cases |
| 08 | sponge-left OR dyn-tol factor |
| 09 | sponge-right OR dyn-tol stride |
| 10 | sponge strength |
| 11 | DMT: omega_c (target frequency) |
| 12 | DMT: beta (filter width) |
| 13 | DMT: chi (feedback gain) |
| 14 | DMT: dmt_skp (steps before forcing) |
| 15 | DMT: dmt_tol (convergence tolerance) |
| 16-20 | reserved -- future stabilization methods |

---

## Template conventions

- `viscosity = -<Re>` encodes `1/Re` (Nek5000: negative => Reynolds number).
- `userParam08/09/10` hold sponge config for stability/postproc stages.
- `startFrom = BF_<name>0.f00001` for stability/postproc; `rst_<name>0.f00001` for DNS.
- `# TEMPLATE:` marks fields the user MUST edit before running.
- Variant subdirs share the parent mesh (.re2/.ma2/.usr/SIZE); `.par` toggles the variant.

---

## What is NOT in a template

- `<basename>.usr` -- boundary conditions, forcing, geometry (case-specific)
- `<basename>.box` / `.re2` / `.ma2` -- mesh files (geometry-specific)
- `SIZE` -- element counts, polynomial order (geometry-specific)
- `SESSION.NAME` -- generated per case; must match .par basename
