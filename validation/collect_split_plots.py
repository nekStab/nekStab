#!/usr/bin/env python3
"""collect_split_plots.py — one-shot collection of split per-panel PNGs/GIFs.

Copies each case's per-panel plot files to validation/figures/ with canonical names.
Skips files that don't exist (partial data is expected).
Also deletes legacy single-panel canonical PNGs that have been superseded.

Run from repo root:
    uv run python validation/collect_split_plots.py
"""
from __future__ import annotations
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FIGURES = REPO / "validation" / "figures"

# Map of (case_dir, [(local_filename, canonical_filename), ...])
MAP = [
    ("example/cylinder_re100/411_postproc_wavemaker", [
        ("plot_wavemaker.png", "cylinder_postproc_wavemaker_field_Re50.png"),
        ("plot_direct_mode.png", "cylinder_postproc_wavemaker_direct_mode_Re50.png"),
        ("plot_adjoint_mode.png", "cylinder_postproc_wavemaker_adjoint_mode_Re50.png"),
    ]),
    ("example/cylinder_re100/412_postproc_steady_force_sensitivity", [
        ("plot_real.png", "cylinder_postproc_steady_force_sensitivity_real_Re50.png"),
        ("plot_imag.png", "cylinder_postproc_steady_force_sensitivity_imag_Re50.png"),
    ]),
    ("example/back_fstep_re500/transient_growth", [
        ("plot_baseflow.png", "back_fstep_transient_growth_baseflow_Re500.png"),
        ("plot_optimal_perturbation.png", "back_fstep_transient_growth_optimal_perturbation_Re500.png"),
        ("plot_optimal_response.png", "back_fstep_transient_growth_optimal_response_Re500.png"),
        ("plot_envelope.png", "back_fstep_transient_growth_envelope_Re500.png"),
    ]),
    ("example/cylinder_re100/410_postproc_animate_modes", [
        ("plot.gif", "cylinder_stability_animate_modes_Re50.gif"),
    ]),
    ("example/thermosyphon_ra500/stability/direct", [
        ("plot_spectrum.png", "thermosyphon_stability_direct_spectrum_Ra500.png"),
        ("plot_baseflow.png", "thermosyphon_stability_direct_baseflow_Ra500.png"),
        ("plot_mode.png", "thermosyphon_stability_direct_mode_Ra500.png"),
    ]),
    ("example/flip_flop_re62/stability/direct_Floquet", [
        ("plot_spectrum.png", "flip_flop_stability_direct_Floquet_spectrum_Re62.png"),
        ("plot_upo_snapshot.png", "flip_flop_stability_direct_Floquet_upo_snapshot_Re62.png"),
        ("plot_mode.png", "flip_flop_stability_direct_Floquet_mode_Re62.png"),
    ]),
    ("example/tpjet_re2005/stability/direct_Floquet", [
        ("plot_spectrum.png", "tpjet_stability_direct_Floquet_spectrum_Re1900.png"),
        ("plot_baseflow.png", "tpjet_stability_direct_Floquet_baseflow_Re1900.png"),
        ("plot_mode.png", "tpjet_stability_direct_Floquet_mode_Re1900.png"),
    ]),
    ("example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn_oifs", [
        ("plot_ic.png",                  "cylinder_baseflow_sfd_casacuberta_dyn_oifs_ic_Re50.png"),
        ("plot_residual.png",            "cylinder_baseflow_sfd_casacuberta_dyn_oifs_residual_Re50.png"),
        ("plot_residual_vs_walltime.png","cylinder_baseflow_sfd_casacuberta_dyn_oifs_residual_vs_walltime_Re50.png"),
        ("plot_bf.png",                  "cylinder_baseflow_sfd_casacuberta_dyn_oifs_bf_Re50.png"),
    ]),
    ("example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs", [
        ("plot_ic.png",                  "cylinder_baseflow_sfd_akervik_dyn_oifs_ic_Re50.png"),
        ("plot_residual.png",            "cylinder_baseflow_sfd_akervik_dyn_oifs_residual_Re50.png"),
        ("plot_residual_vs_walltime.png","cylinder_baseflow_sfd_akervik_dyn_oifs_residual_vs_walltime_Re50.png"),
        ("plot_bf.png",                  "cylinder_baseflow_sfd_akervik_dyn_oifs_bf_Re50.png"),
    ]),
    ("example/cylinder_re100/600_modal_pod", [
        ("plot_snapshot.png",      "cylinder_modal_snapshot_Re100.png"),
        ("plot_pod_spectrum.png",  "cylinder_modal_pod_spectrum_Re100.png"),
        ("plot_pod_mode1.png",     "cylinder_modal_pod_mode1_Re100.png"),
        ("plot_pod_mode2.png",     "cylinder_modal_pod_mode2_Re100.png"),
        ("plot_dmd_spectrum.png",  "cylinder_modal_dmd_spectrum_Re100.png"),
        ("plot_dmd_mode1.png",     "cylinder_modal_dmd_mode1_Re100.png"),
        ("plot_spod_spectrum.png", "cylinder_modal_spod_spectrum_Re100.png"),
        ("plot_spod_mode1.png",    "cylinder_modal_spod_mode1_Re100.png"),
    ]),
    ("example/cylinder_re100/110_baseflow_sfd/akervik_dyn", [
        ("plot_convergence.png", "cylinder_baseflow_sfd_akervik_dyn_convergence_Re50.png"),
        ("plot_dynamic_tolerance.png", "cylinder_baseflow_sfd_akervik_dyn_dynamic_tolerance_Re50.png"),
        ("plot_baseflow.png", "cylinder_baseflow_sfd_akervik_dyn_baseflow_Re50.png"),
    ]),
    ("example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn", [
        ("plot_convergence.png", "cylinder_baseflow_sfd_casacuberta_dyn_convergence_Re50.png"),
        ("plot_dynamic_tolerance.png", "cylinder_baseflow_sfd_casacuberta_dyn_dynamic_tolerance_Re50.png"),
        ("plot_baseflow.png", "cylinder_baseflow_sfd_casacuberta_dyn_baseflow_Re50.png"),
    ]),
    ("example/poiseuille_re5k", [
        ("plot_lyapunov_exponents.png", "poiseuille_lyapunov_exponents_Re5000.png"),
        ("plot_otd_residuals.png", "poiseuille_otd_residuals_Re5000.png"),
        ("plot_otd_mode1.png", "poiseuille_otd_mode1_Re5000.png"),
        ("plot_otd_mode2.png", "poiseuille_otd_mode2_Re5000.png"),
    ]),
    ('example/cylinder_re100/000_dns', [
        ('plot_snapshot.png', 'cylinder_dns_snapshot_Re50.png'),
        ('plot_signal.png', 'cylinder_dns_signal_Re50.png'),
        ('plot_fft.png', 'cylinder_dns_fft_Re50.png'),
        ('plot_phase.png', 'cylinder_dns_phase_Re50.png'),
    ]),
    ('example/cylinder_re100/110_baseflow_sfd/akervik', [
        ('plot_ic.png', 'cylinder_baseflow_sfd_ic_Re50.png'),
        ('plot_residual.png', 'cylinder_baseflow_sfd_residual_Re50.png'),
        ('plot_bf.png', 'cylinder_baseflow_sfd_bf_Re50.png'),
    ]),
    ('example/cylinder_re100/110_baseflow_sfd/casacuberta', [
        ('plot_ic.png', 'cylinder_baseflow_sfd_casacuberta_ic_Re50.png'),
        ('plot_residual.png', 'cylinder_baseflow_sfd_casacuberta_residual_Re50.png'),
        ('plot_bf.png', 'cylinder_baseflow_sfd_casacuberta_bf_Re50.png'),
    ]),
    ('example/cylinder_re100/110_baseflow_sfd/akervik_dyn', [
        ('plot_ic.png', 'cylinder_baseflow_sfd_akervik_dyn_ic_Re50.png'),
        ('plot_residual.png', 'cylinder_baseflow_sfd_akervik_dyn_residual_Re50.png'),
        ('plot_dynamic_tolerance.png', 'cylinder_baseflow_sfd_akervik_dyn_dynamic_tolerance_Re50.png'),
        ('plot_bf.png', 'cylinder_baseflow_sfd_akervik_dyn_bf_Re50.png'),
    ]),
    ('example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn', [
        ('plot_ic.png', 'cylinder_baseflow_sfd_casacuberta_dyn_ic_Re50.png'),
        ('plot_residual.png', 'cylinder_baseflow_sfd_casacuberta_dyn_residual_Re50.png'),
        ('plot_dynamic_tolerance.png', 'cylinder_baseflow_sfd_casacuberta_dyn_dynamic_tolerance_Re50.png'),
        ('plot_bf.png', 'cylinder_baseflow_sfd_casacuberta_dyn_bf_Re50.png'),
    ]),
    ('example/cylinder_re100/120_baseflow_boostconv', [
        ('plot_ic.png', 'cylinder_baseflow_boostconv_ic_Re50.png'),
        ('plot_residual.png', 'cylinder_baseflow_boostconv_residual_Re50.png'),
        ('plot_bf.png', 'cylinder_baseflow_boostconv_bf_Re50.png'),
    ]),
    ('example/cylinder_re100/210_baseflow_newton/fp', [
        ('plot_ic.png', 'cylinder_baseflow_newton_ic_Re50.png'),
        ('plot_residual.png', 'cylinder_baseflow_newton_residual_Re50.png'),
        ('plot_bf.png', 'cylinder_baseflow_newton_bf_Re50.png'),
    ]),
    ('example/cylinder_re100/210_baseflow_newton/dyn', [
        ('plot_ic.png', 'cylinder_baseflow_newton_dyn_ic_Re50.png'),
        ('plot_residual.png', 'cylinder_baseflow_newton_dyn_residual_Re50.png'),
        ('plot_bf.png', 'cylinder_baseflow_newton_dyn_bf_Re50.png'),
    ]),
    ('example/cylinder_re180/210_baseflow_newton/upo', [
        ('plot_ic.png', 'cylinder_baseflow_newton_upo_ic_Re180.png'),
        ('plot_residual.png', 'cylinder_baseflow_newton_upo_residual_Re180.png'),
        ('plot_bf.png', 'cylinder_baseflow_newton_upo_bf_Re180.png'),
    ]),
    ('example/back_fstep_re500/baseflow', [
        ('plot_ic.png', 'back_fstep_baseflow_ic_Re500.png'),
        ('plot_residual.png', 'back_fstep_baseflow_residual_Re500.png'),
        ('plot_bf.png', 'back_fstep_baseflow_bf_Re500.png'),
    ]),
    ('example/flip_flop_re62/baseflow', [
        ('plot_ic.png', 'flip_flop_baseflow_ic_Re62.png'),
        ('plot_residual.png', 'flip_flop_baseflow_residual_Re62.png'),
        ('plot_bf.png', 'flip_flop_baseflow_bf_Re62.png'),
    ]),
    ('example/tpjet_re2005/baseflow/newton', [
        ('plot_ic.png', 'tpjet_baseflow_newton_ic_Re1900.png'),
        ('plot_residual.png', 'tpjet_baseflow_newton_residual_Re1900.png'),
        ('plot_bf.png', 'tpjet_baseflow_newton_bf_Re1900.png'),
    ]),
    ('example/tpjet_re2005/baseflow/tdf', [
        ('plot_ic.png', 'tpjet_baseflow_tdf_ic_Re2005.png'),
        ('plot_bf.png', 'tpjet_baseflow_tdf_bf_Re2005.png'),
    ]),
    ('example/naca0012_re2000', [
        ('plot_snapshot.png', 'naca0012_dns_snapshot_Re2000.png'),
        ('plot_signal.png', 'naca0012_dns_signal_Re2000.png'),
        ('plot_fft.png', 'naca0012_dns_fft_Re2000.png'),
        ('plot_phase.png', 'naca0012_dns_phase_Re2000.png'),
    ]),
    ('example/slot_fst_re495', [
        ('plot_signal.png', 'slot_FST_dns_signal_Re495.png'),
        ('plot_fft.png', 'slot_FST_dns_fft_Re495.png'),
        ('plot_phase.png', 'slot_FST_dns_phase_Re495.png'),
    ]),
    ('example/lid_driven_re3600', [
        ('plot_ic.png', 'lid_driven_baseflow_ic_Re3600.png'),
        ('plot_residual.png', 'lid_driven_baseflow_residual_Re3600.png'),
        ('plot_bf.png', 'lid_driven_baseflow_bf_Re3600.png'),
    ]),
    ('example/thermosyphon_ra500/baseflow', [
        ('plot_ic.png', 'thermosyphon_baseflow_ic_Ra500.png'),
        ('plot_residual.png', 'thermosyphon_baseflow_residual_Ra500.png'),
        ('plot_bf.png', 'thermosyphon_baseflow_bf_Ra500.png'),
    ]),
    ('example/cylinder_re30_thermal/210_baseflow_newton/dyn_temp', [
        ('plot_residual.png', 'cylinder_baseflow_newton_dyn_temp_residual_Re50.png'),
    ]),
    ('example/cylinder_re100/320_stability_adjoint', [
        ('plot_bf.png', 'cylinder_stability_adjoint_bf_Re50.png'),
        ('plot_spectrum.png', 'cylinder_stability_adjoint_spectrum_Re50.png'),
        ('plot_mode.png', 'cylinder_stability_adjoint_mode_Re50.png'),
    ]),
    ('example/cylinder_re180/311_stability_direct_floquet', [
        ('plot_spectrum.png', 'cylinder_stability_direct_Floquet_spectrum_Re180.png'),
        ('plot_mode.png', 'cylinder_stability_direct_Floquet_mode_Re180.png'),
        ('plot_bf.png', 'cylinder_stability_direct_Floquet_bf_Re180.png'),
    ]),
    ('example/cylinder_re180/321_stability_adjoint_floquet', [
        ('plot_spectrum.png', 'cylinder_stability_adjoint_Floquet_spectrum_Re180.png'),
        ('plot_mode.png', 'cylinder_stability_adjoint_Floquet_mode_Re180.png'),
        ('plot_bf.png', 'cylinder_stability_adjoint_Floquet_bf_Re180.png'),
    ]),
    ('example/cylinder_re180/410_postproc_animate_modes/upo', [
        ('plot.gif', 'cylinder_stability_animate_modes_with_UPO_Re180.gif'),
        ('plot.png', 'cylinder_stability_animate_modes_with_UPO_still_Re180.png'),
    ]),
    ('example/cylinder_re100/310_stability_direct/direct', [
        ('plot_bf.png',              'cylinder_stability_direct_bf_Re50.png'),
        ('plot_spectrum_circle.png', 'cylinder_stability_direct_spectrum_circle_Re50.png'),
        ('plot_spectrum_ns_log.png', 'cylinder_stability_direct_spectrum_ns_log_Re50.png'),
        ('plot_mode_vx.png',         'cylinder_stability_direct_mode_vx_Re50.png'),
        ('plot_mode_vy.png',         'cylinder_stability_direct_mode_vy_Re50.png'),
    ]),

]

# Legacy canonical PNGs that are superseded by the split per-panel files.
# These will be deleted if they exist.
LEGACY_PNGS = [
    "cylinder_postproc_sensitivity_budget_wavemaker_Re50.png",
    "cylinder_postproc_steady_force_sensitivity_Re50.png",
    "back_fstep_transient_growth_Re500.png",
    "cylinder_stability_animate_modes_Re50.png",  # replaced by .gif
    "cylinder_baseflow_newton_dyn_temp_Re50.png",
    "thermosyphon_baseflow_Ra500.png",
    "thermosyphon_stability_direct_Ra500.png",
    "flip_flop_stability_direct_Floquet_Re62.png",
    "tpjet_stability_direct_Floquet_Re1900.png",
    "cylinder_stability_direct_Floquet_Re180.png",   # superseded by _spectrum/_mode/_bf split
    "cylinder_stability_adjoint_Floquet_Re180.png",  # superseded by _spectrum/_mode/_bf split
    "cylinder_stability_animate_modes_with_UPO_Re100.png",  # stale Re100 name, now Re180 gif/still
    "lid_driven_Re3600.png",
    "cylinder_baseflow_sfd_akervik_dyn_Re50.png",
    "cylinder_baseflow_sfd_casacuberta_dyn_Re50.png",
    "cylinder_baseflow_sfd_akervik_dyn_oifs_Re50.png",
    "cylinder_baseflow_sfd_casacuberta_dyn_oifs_Re50.png",
    "poiseuille_Re5000.png",
]


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    n_copied = 0
    n_skipped = 0
    n_deleted = 0

    for case_dir_rel, file_map in MAP:
        case_dir = REPO / case_dir_rel
        for local_name, canonical_name in file_map:
            src = case_dir / local_name
            dst = FIGURES / canonical_name
            if src.exists():
                shutil.copy2(src, dst)
                print(f"COPY  {case_dir_rel}/{local_name} -> validation/figures/{canonical_name}")
                n_copied += 1
            else:
                print(f"SKIP  {case_dir_rel}/{local_name}  (not found)")
                n_skipped += 1

    for legacy_name in LEGACY_PNGS:
        legacy_path = FIGURES / legacy_name
        if legacy_path.exists():
            legacy_path.unlink()
            print(f"DEL   validation/figures/{legacy_name}")
            n_deleted += 1
        else:
            print(f"SKIP-DEL  validation/figures/{legacy_name}  (not present)")

    print()
    print(f"Summary: {n_copied} copied, {n_skipped} skipped, {n_deleted} deleted.")


if __name__ == "__main__":
    main()
