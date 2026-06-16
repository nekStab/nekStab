#!/usr/bin/env bash
# mkstage.sh — instantiate one numbered example stage from example/_templates/.
#
# Mechanical part of the "every case in place" campaign: copies the stage
# template .par, sets Re + startFrom + sponge, hardlinks the (read-only) mesh and
# copies the build infra from the geometry root, writes a run.local.slurm, and
# refreshes the README scaffold. Physics judgment (which stages apply, which
# startFrom) is the caller's; this script only does the file mechanics.
#
# Usage:
#   mkstage.sh <geom> <stage_rel> <case> <Re> <startFrom|cold> [ntasks] [spongeL,R,S] [time]
# Example:
#   mkstage.sh cubic_cavity_re1914 310_stability_direct/direct cav 1914 BF_cav0.f00001 16 0,0,0 12:00:00
set -euo pipefail

EX="$(cd "$(dirname "$0")/../example" && pwd)"
geom="$1"; stage="$2"; case="$3"; Re="$4"; startfrom="$5"
ntasks="${6:-8}"; sponge="${7:-0,0,0}"; walltime="${8:-12:00:00}"

root="$EX/$geom"
dst="$root/$stage"
[ -d "$root" ] || { echo "ERR: geom root missing: $root" >&2; exit 1; }

# template = first path component of stage (NNN_name); variant subpath kept if present
nnn="${stage%%/*}"; variant="${stage#"$nnn"}"; variant="${variant#/}"
tdir="$EX/_templates/$nnn${variant:+/$variant}"
[ -f "$tdir/case.par" ] || tdir="$EX/_templates/$nnn"      # fall back to base template
[ -f "$tdir/case.par" ] || { echo "ERR: no template for $stage ($tdir)" >&2; exit 1; }

mkdir -p "$dst"
cp "$tdir/case.par" "$dst/$case.par"

# --- edit .par: Re, startFrom, sponge ---
sed -i "s|^viscosity = .*|viscosity = -${Re}.0   # Re=${Re}|" "$dst/$case.par"
if [ "$startfrom" = "cold" ]; then
  sed -i "s|^startFrom = .*|# startFrom =            # cold start|" "$dst/$case.par"
else
  sed -i "s|^startFrom = .*|startFrom = ${startfrom}|" "$dst/$case.par"
fi
IFS=, read -r sL sR sS <<<"$sponge"
sed -i "s|^userParam08 = .*|userParam08 = ${sL}          # sponge left|"     "$dst/$case.par" 2>/dev/null || true
sed -i "s|^userParam09 = .*|userParam09 = ${sR}          # sponge right|"    "$dst/$case.par" 2>/dev/null || true
sed -i "s|^userParam10 = .*|userParam10 = ${sS}          # sponge strength|" "$dst/$case.par" 2>/dev/null || true

# --- mesh (hardlink, read-only) + infra (copy) ---
for m in "$case.re2" "$case.ma2" "$case.box"; do [ -f "$root/$m" ] && ln -f "$root/$m" "$dst/$m"; done
for f in SIZE "$case.usr" NEKSTAB.inc makefile makefile_usr.inc; do [ -f "$root/$f" ] && cp "$root/$f" "$dst/$f"; done
# seed field present in root and named as startFrom -> copy (so a run can read it; never hardlink, runs overwrite)
[ "$startfrom" != "cold" ] && [ -f "$root/$startfrom" ] && cp "$root/$startfrom" "$dst/$startfrom"

mem=2000; [ "$ntasks" -ge 16 ] && mem=3000
cat > "$dst/run.local.slurm" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${geom}-${nnn}
#SBATCH --partition=local
#SBATCH --nodes=1
#SBATCH --ntasks=${ntasks}
#SBATCH --mem-per-cpu=${mem}
#SBATCH --time=${walltime}
#SBATCH --output=logfile
#SBATCH --error=logfile

set -euo pipefail
printf "%s\n%s/\n" "${case}" "\$(pwd)" > SESSION.NAME
mpiexec -np "\$SLURM_NTASKS" ./nek5000
EOF

# --- README (concise, instantiated) ---
mode=$(grep -oE "userParam01 = [0-9.]+" "$dst/$case.par" | awk '{print $3}')
ic="$startfrom"; [ "$startfrom" = "cold" ] && ic="cold start (useric)"

# Guard / note for mode/uparam consistency (see src/mode_codes.f90 and mode_config.f90 for SSOT).
# Valid userParam01 values are defined centrally; templates and .par should match EXAMPLES.md catalog.
# If adding new mode, update mode_codes + mkstage comment + template README scaffold.

cat > "$dst/README.md" <<EOF
# ${geom} / ${stage}

**Stage**: ${nnn}${variant:+ (${variant})}
**uparam01**: ${mode}
**Re**: ${Re}
**Start from**: ${ic}
**Ranks**: ${ntasks}

Stage rationale is shared in the repo-root \`EXAMPLES.md\`. Run with
\`sbatch run.local.slurm\`. Sponge disabled (closed/internal geometry) where
applicable; tune endTime/k_dim for this operating point before the final run.
EOF

echo "OK: $geom/$stage"
grep -E "^startFrom|^userParam01|^viscosity|^userParam0[789]|^userParam10" "$dst/$case.par" | sed 's/^/   /'
if grep -q "TEMPLATE" "$dst/$case.par"; then echo "   WARN: residual TEMPLATE markers remain"; fi
