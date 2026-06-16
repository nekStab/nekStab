# Fortran Coding Style and Guiding Conventions for nekStab

This document is the reference for the coding style and conventions used in the nekStab Fortran sources. It aims for a codebase that stays consistent, readable, and maintainable.

## Source Format

- Modern free-form is the default for nekStab modules. Lines are at most **132 columns**.
- The Nek5000 interface files are **fixed-form** (code in columns 7+, continuation marker in column 6): `nekstab_nek_bridge.f90`, `main.f90`, `usr_wrappers.f90`, and the `NEKSTAB*` include. They compile in Nek5000's fixed-form `.usr` context, so their column layout and hand-wrapped continuations are part of the spec and must be preserved.
- No tabs. No trailing whitespace.
- Free-form modules indent **3 spaces per level**. Continuation uses `&`; aligned continuation in `public ::` lists and long argument lists is the house style (see `energy_budget.f90`).

## Lexical Conventions

- **Comments** use `!`. Legacy column-1 `c` / `C` / `*` comment markers are not used.
- **Keywords and intrinsics** are lowercase (`do`, `end if`, `sqrt`, `.true.`). User-defined names keep their own casing.
- **Relational operators** use the symbolic form: `>`, `>=`, `<`, `<=`, `==`, `/=` — not `.gt.`, `.ge.`, `.lt.`, `.le.`, `.eq.`, `.ne.`.
- One space around binary operators and after separators; the `::` in declarations is aligned within a declaration group.

## Module Structure

Free-form modules follow a fixed skeleton — the conventions below are already uniform across the sources, and keeping them explicit prevents the most common build breaks (a missing `public` line is a link failure, not a compile error).

- Every module declares `implicit none`.
- `private` by default; export the public API with an explicit `public ::` list (see the header template above, which mirrors that list in prose).
- All dummy arguments carry an `intent(in)`, `intent(out)`, or `intent(inout)`.
- Watch for name collisions with `include 'TOTAL'`: it pulls in many short names (e.g. `fw`), so a local variable of the same name causes cascading declaration errors. Use a distinct name (`fwrk`, not `fw`).

## Comment Style

### Module / File Level (Template)
Use a consistent header:

```fortran
!-----------------------------------------------------------------------
! filename.f90 — One-line title (keep it natural, not marketing)
!
! Purpose:
!   2–4 lines. What this actually does in the context of stability work.
!   Mention historical origin or tricky part if relevant.
!
! Public interface:
!   foo  — short, lowercase description
!   bar  — ...
!
! Dependencies:
!   list real ones
!
! See also:
!   otherfile.f90 — why it matters
!-----------------------------------------------------------------------
```

- Use an em dash `—` (not `--`) in the title line and in the top-level Public interface list.
- Structured blocks only where they genuinely help discoverability. Do not force identical sub-routine banners on every tiny helper.

### Inside Routines
Preferred — explanatory and useful:

```fortran
!  WHY: robust default; clean-converging cases accept alpha=1 on trial 0
!  and stay bit-identical to the old full-step behavior.

!  Note: if ifseed_* are all false, the 'useric' subroutine in the .usr
!  prescribes the initial seed. This is the historical path.

! Sponge parameters. This follows the Nordstrom/Nordin/Henningson
! smooth-mask formulation (1999) as used in the original fringe setup.
```

Short inline comments are acceptable and stay terse — do **not** expand them into full sentences:

```fortran
k_dim = 100          ! usual starting value; users override in .usr
ifnewton_backtrack = .true.   ! default on in v2.0; see backtracking section
```

### Section Dividers
Allow light variation — do not standardize to a single character:

```fortran
!  ── Log file management ──────────────────────────────────────────
!  User-settable mode flags (alternative to uparam encoding)
```

### Existing Comments
- Preserve explanatory blocks ("The old code did...", "WHY:", historical norms) verbatim — content and voice intact.
- Short inline comments with minor imperfections (brevity, terse phrasing) are intentional; do not mass-rewrite them.
- A homogenized header keeps at least one explanatory sentence or historical note.
- The comment marker is `!`; the comment text itself is never altered to fit the marker.

Reference files that already match this style: `energy_budget.f90`, `mode_config.f90`, `otd.f90`, `newton_krylov.f90` (log management block), `fixedp.f90` (terse close helpers with historical note).

## Close Lines (return and end subroutine / function / module)

### Core Rule
- Every implementation body closes with the **full name**:
  - `end subroutine <exact_name>`
  - `end function <exact_name>`
  - `end module <exact_name>`
- Never a bare `end subroutine` or lone `end` for actual procedures.

### Breathing Room
- Prefer exactly one blank line before `end subroutine` (etc.) after a terminating `end if` / `end do` (or for clean separation after the last statement).
- This is for readability, not mandatory if the last statement is already a clean close.

### Return Statements
- **Preserve all strategic and early `return` statements**, especially the classic Nek userchk dispatcher pattern.
- In the central dispatcher (`nekStab` in `main.f90`):
  - Keep explicit `return` after `nek_end` + `ifbfcv` for mode early exits.
  - Do **not** add a `return` immediately before the final `end subroutine nekStab` — let control fall off the if-ladder. This is the "skip the rest of the timestep" idiom.
- In normal subroutines: omit an explicit `return` right before the close unless it clarifies non-obvious control flow (e.g., after a set of early returns inside an if-ladder).
- Guard returns at the top of master-only / I/O subs:
  ```fortran
  if (nid /= 0) return
  ```
  Add a short comment only when the reason is not obvious from context:
  ```fortran
  if (nid /= 0) return   ! master rank only for file I/O and printing
  ```

### Visual Close
- The `!-------` comment that often follows the final `end subroutine` in `main.f90` is fine as a soft file-level separator.
- Inside modules, no need to decorate the `end subroutine` line itself.

## Numeric Literals (Double Precision)

**Rule**: In double-precision contexts (d0 suffix, `nekStab_dp`, `C_DOUBLE`, double arrays/params, real computations, time, eps, tols, eigenvalue math, etc.) use explicit full forms:

```fortran
B = 0.0d0
x1 = -b / 2.0d0 / a
if (d >= 0.0d0) then
...
ttime = 0.0d0
telapsed = ttime / 3600.0d0
eek = 0.50d0 / volvm1
PI = 3.141592653589793d0
TWOPI = 2.0d0 * PI
```

- Avoid ambiguous or bare forms in double code: `0.`, `0.d0`, `60.`, `2.`, `1.0/3.0`, `11./6.`, `>= 0.`, `1.-beta`, `1e-6`, `0.99e30` (etc.).
- Chained expressions and comparisons must be explicit.
- Keep genuine single-precision, uparam mode decimals (legacy single-float encoding for mode selection), format strings, and sentinels exactly as-is.
- Preserve any existing comments around literals (e.g., "historical norm").

This makes double intent explicit and prevents accidental single-precision promotion. The reason it matters here: Nek5000 compiles with `-r8`, so default `real` variables are already double, but a bare literal like `0.` is still a single-precision constant — it loses digits before being promoted. The explicit `d0` keeps the constant double from the start.