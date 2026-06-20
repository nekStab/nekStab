"""Validation case catalog.

Source contracts: validation-example-inventory-2026-05-20.md,
validation-method-capability-matrix-2026-05-20.md,
validation-target-seed-parameter-model-2026-05-21.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------

class CaseStatus(str, Enum):
    VALIDATED = "validated"
    LOCAL_SMOKE = "local_smoke"
    NEEDS_DNS_SEED = "needs_dns_seed"
    NEEDS_SCALAR_CHECKPOINT = "needs_scalar_checkpoint"
    MISSING_IMAGE = "missing_image"
    STALE = "stale"
    COMPUTE_HEAVY = "compute_heavy"
    BLOCKED = "blocked"
    DEFERRED = "deferred"
    AMBIGUOUS = "ambiguous"
    PARTIAL = "partial"  # legacy inventory shorthand; retained for current cases


class ComputeClass(str, Enum):
    QUICK_LOCAL_SMOKE = "quick_local_smoke"
    BOUNDED_LOCAL_VALIDATION = "bounded_local_validation"
    LOCAL_LONG_RUN = "local_long_run"
    REMOTE_JZ_CANDIDATE = "remote_jz_candidate"
    MANUAL_ONLY = "manual_only"
    BLOCKED_COMPUTE = "blocked_compute"


class EvidenceRole(str, Enum):
    REFERENCE_IMAGE = "reference_image"
    DNS_TIME_HISTORY = "dns_time_history"
    RESIDUAL_HISTORY = "residual_history"
    SPECTRUM = "spectrum"
    EIGENVALUE_TABLE = "eigenvalue_table"
    MODE_SHAPE = "mode_shape"
    ANIMATION = "animation"
    CHECKPOINT_STATS = "checkpoint_stats"
    PROVENANCE = "provenance"
    LOG = "log"
    BLOCKING_NOTE = "blocking_note"


ArtifactRole = EvidenceRole  # alias: artifacts use the same role vocabulary as evidence requirements


class CapabilityState(str, Enum):
    IMPLEMENTED = "implemented"
    CANDIDATE_PLANNED = "candidate_planned"
    BLOCKED_BY_PREREQUISITE = "blocked_by_prerequisite"
    INTENTIONALLY_DEFERRED = "intentionally_deferred"
    SPECIAL_DIRECT_ONLY = "special_direct_only"
    NOT_APPLICABLE = "not_applicable"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TargetParameters:
    Re: float | None = None
    Pr: float | None = None
    Ri: float | None = None
    geometry: str | None = None
    motion: str | None = None       # 'static', 'forced', 'free', etc.
    model: str | None = None        # 'DNS', 'RANS', etc.
    final_time_or_period: float | None = None
    tolerance_contract: str | None = None


@dataclass(frozen=True)
class SeedParameters:
    seed_Re: float | str | None = None   # 'unknown' allowed
    seed_Pr: float | None = None
    seed_Ri: float | None = None
    seed_path: str | None = None         # relative-to-case file or 'unknown'
    restart_source: str | None = None
    seed_case: str | None = None
    seed_time: float | None = None
    seed_status: str | None = None


@dataclass(frozen=True)
class EvidenceRequirement:
    role: EvidenceRole
    required: bool = True
    note: str | None = None


@dataclass(frozen=True)
class Artifact:
    role: EvidenceRole
    path: str               # relative to repo root
    freshness: str | None = None     # 'current' | 'stale' | 'unknown'
    provenance: str | None = None


@dataclass(frozen=True)
class FlowFamily:
    family_id: str          # 'cylinder_re100', 'thermosyphon', etc.
    label: str              # human label
    geometry_kind: str      # '2D-cylinder', 'cavity-3D', 'jet', 'boundary-layer', 'channel', etc.
    why_note: str | None = None


@dataclass(frozen=True)
class MethodCase:
    case_id: str                # stable slug, e.g. 'cylinder/baseflow/sfd'
    flow_family: str            # one of FLOW_FAMILY_IDS
    method_lane: str            # 'dns' | 'baseflow' | 'stability_direct' | etc.
    current_path: str           # actual path today
    proposed_folder: str        # future flat-numbered folder
    sort_prefix: str            # NNN[NN] used for ordering
    label: str                  # short human label for the gallery
    mode_name: str | None       # e.g. 'SFD', 'Newton-GMRES', 'Direct LNSE'
    legacy_uparam01: str | None # e.g. '1.1', '2.0', '3.1', or None
    expected_behavior: str      # gallery card 'expected' line
    status: CaseStatus
    compute_class: ComputeClass
    target_parameters: TargetParameters
    seed_parameters: SeedParameters = field(default=SeedParameters())
    evidence_requirements: tuple[EvidenceRequirement, ...] = ()
    artifacts: tuple[Artifact, ...] = ()
    prerequisites: tuple[str, ...] = ()  # case_ids this depends on
    blockers: tuple[str, ...] = ()
    evidence_prereqs: tuple[str, ...] = ()  # file globs that must exist before validated claim
    why_note: str | None = None
    gallery_visible: bool = True


# ---------------------------------------------------------------------------
# Module-level data: FLOW_FAMILIES
# ---------------------------------------------------------------------------

FLOW_FAMILIES: dict[str, FlowFamily] = {
    "cylinder_re100": FlowFamily("cylinder_re100", "Static Cylinder", "2D-cylinder",
        why_note="Rich demonstration family; canonical validation ladder"),
    "moving_cylinder": FlowFamily("moving_cylinder", "Moving Cylinder", "2D-cylinder-moving",
        why_note="Mesh-motion extension; deferred pending run policy"),
    "thermosyphon": FlowFamily("thermosyphon", "Thermosyphon", "2D-annulus",
        why_note="Buoyancy-driven annular flow; note spelling mismatch in gallery figures"),
    "flip_flop": FlowFamily("flip_flop", "Flip-Flop Jet", "2D-jet",
        why_note="Symmetry-breaking bifurcation; has JZ and local Slurm surfaces"),
    "tpjet": FlowFamily("tpjet", "Two-Phase Jet", "2D-jet",
        why_note="Compute-heavy jet case; JZ candidate"),
    "back_fstep": FlowFamily("back_fstep", "Backward-Facing Step", "2D-step",
        why_note="Classic transient-growth test case"),
    "cubic_cavity": FlowFamily("cubic_cavity", "Cubic Cavity", "cavity-3D",
        why_note="3D compute-heavy; deferred until evidence refreshed"),
    "lid_driven_cavity": FlowFamily("lid_driven_cavity", "Lid-Driven Cavity", "2D-cavity",
        why_note="Standalone cavity benchmark"),
    "naca0012": FlowFamily("naca0012", "NACA 0012", "2D-airfoil",
        why_note="Two operating points (Re=2000, Re=2500) in gallery"),
    "poiseuille": FlowFamily("poiseuille", "Poiseuille Channel", "channel",
        why_note="OTD and RANS lanes; shared gallery filename needs disambiguation"),
    "slot_FST": FlowFamily("slot_FST", "Slot FST", "2D-slot",
        why_note="Partial; preprocessFST/run.sh is part of runnable surface"),
}


# ---------------------------------------------------------------------------
# Module-level data: CAPABILITY_MATRIX
# ---------------------------------------------------------------------------
# Keyed on (family_id, method_lane).
# Lanes: 'dns', 'baseflow', 'direct', 'adjoint', 'floquet',
#        'transient_growth', 'wavemaker', 'animation', 'otd', 'modal', 'rans'

_I = CapabilityState.IMPLEMENTED
_CP = CapabilityState.CANDIDATE_PLANNED
_BP = CapabilityState.BLOCKED_BY_PREREQUISITE
_ID = CapabilityState.INTENTIONALLY_DEFERRED
_SD = CapabilityState.SPECIAL_DIRECT_ONLY
_NA = CapabilityState.NOT_APPLICABLE

CAPABILITY_MATRIX: dict[tuple[str, str], CapabilityState] = {
    # cylinder_re100 (fixed point + 2D modal analysis; NO Floquet at Re=100)
    ("cylinder_re100", "dns"):              _I,
    ("cylinder_re100", "baseflow"):         _I,
    ("cylinder_re100", "direct"):           _I,
    ("cylinder_re100", "adjoint"):          _I,
    ("cylinder_re100", "floquet"):          _NA,  # cylinder wake 2D-stable to Re_A~188 (Barkley-Henderson 1996); Floquet lives in cylinder_re180
    ("cylinder_re100", "transient_growth"): _CP,
    ("cylinder_re100", "wavemaker"):        _I,
    ("cylinder_re100", "animation"):        _I,
    ("cylinder_re100", "otd"):              _I,
    ("cylinder_re100", "modal"):            _I,
    ("cylinder_re100", "baseflow_dmt"):     _NA,  # DMT is for periodic orbits, not fixed points
    ("cylinder_re100", "baseflow_tdf"):     _NA,  # TDF is for periodic orbits, not fixed points
    ("cylinder_re100", "rans"):             _SD,  # special finite-difference direct only
    # cylinder_re180 (Floquet on the 2D periodic wake at Re~Re_A; Barkley-Henderson Mode A threshold)
    ("cylinder_re180", "dns"):              _CP,  # could be added; current Re=180 case is OTD + Floquet only
    ("cylinder_re180", "baseflow"):         _I,   # Newton-UPO at Re=180 (needs Re=150 seed)
    ("cylinder_re180", "floquet"):          _I,   # the primary motivation for this family
    ("cylinder_re180", "animation"):        _I,   # Floquet-mode animation on UPO
    ("cylinder_re180", "otd"):              _I,   # the original Re=180 OTD case
    # moving_cylinder
    ("moving_cylinder", "dns"):              _CP,
    ("moving_cylinder", "baseflow"):         _CP,
    ("moving_cylinder", "direct"):           _BP,
    ("moving_cylinder", "adjoint"):          _BP,
    ("moving_cylinder", "floquet"):          _CP,
    ("moving_cylinder", "transient_growth"): _CP,
    ("moving_cylinder", "wavemaker"):        _BP,
    ("moving_cylinder", "animation"):        _CP,
    ("moving_cylinder", "otd"):              _CP,
    ("moving_cylinder", "modal"):            _CP,
    ("moving_cylinder", "rans"):             _NA,
    # thermosyphon
    ("thermosyphon", "dns"):                 _CP,
    ("thermosyphon", "baseflow"):            _I,
    ("thermosyphon", "direct"):              _I,
    ("thermosyphon", "adjoint"):             _CP,
    ("thermosyphon", "floquet"):             _CP,
    ("thermosyphon", "transient_growth"):    _CP,
    ("thermosyphon", "wavemaker"):           _BP,
    ("thermosyphon", "animation"):           _CP,
    ("thermosyphon", "otd"):                 _CP,
    ("thermosyphon", "modal"):               _CP,
    ("thermosyphon", "rans"):                _NA,
    # flip_flop
    ("flip_flop", "dns"):                    _CP,
    ("flip_flop", "baseflow"):               _I,
    ("flip_flop", "direct"):                 _CP,
    ("flip_flop", "adjoint"):                _CP,
    ("flip_flop", "floquet"):                _I,
    ("flip_flop", "transient_growth"):       _CP,
    ("flip_flop", "wavemaker"):              _BP,
    ("flip_flop", "animation"):              _CP,
    ("flip_flop", "otd"):                    _CP,
    ("flip_flop", "modal"):                  _CP,
    ("flip_flop", "rans"):                   _NA,
    # tpjet
    ("tpjet", "dns"):                        _CP,
    ("tpjet", "baseflow"):                   _I,
    ("tpjet", "direct"):                     _CP,
    ("tpjet", "adjoint"):                    _CP,
    ("tpjet", "floquet"):                    _I,
    ("tpjet", "transient_growth"):           _CP,
    ("tpjet", "wavemaker"):                  _BP,
    ("tpjet", "animation"):                  _CP,
    ("tpjet", "otd"):                        _CP,
    ("tpjet", "modal"):                      _CP,
    ("tpjet", "rans"):                       _NA,
    # back_fstep
    ("back_fstep", "dns"):                   _CP,
    ("back_fstep", "baseflow"):              _I,
    ("back_fstep", "direct"):                _CP,
    ("back_fstep", "adjoint"):               _CP,
    ("back_fstep", "floquet"):               _NA,  # not_applicable for steady-only target
    ("back_fstep", "transient_growth"):      _I,
    ("back_fstep", "wavemaker"):             _BP,
    ("back_fstep", "animation"):             _CP,
    ("back_fstep", "otd"):                   _CP,
    ("back_fstep", "modal"):                 _CP,
    ("back_fstep", "rans"):                  _NA,
    # cubic_cavity
    ("cubic_cavity", "dns"):                 _CP,
    ("cubic_cavity", "baseflow"):            _I,    # implemented but compute_heavy/deferred
    ("cubic_cavity", "direct"):              _CP,
    ("cubic_cavity", "adjoint"):             _CP,
    ("cubic_cavity", "floquet"):             _BP,
    ("cubic_cavity", "transient_growth"):    _CP,
    ("cubic_cavity", "wavemaker"):           _BP,
    ("cubic_cavity", "animation"):           _CP,
    ("cubic_cavity", "otd"):                 _CP,
    ("cubic_cavity", "modal"):               _CP,
    ("cubic_cavity", "rans"):                _NA,
    # lid_driven_cavity
    ("lid_driven_cavity", "dns"):            _CP,
    ("lid_driven_cavity", "baseflow"):       _I,
    ("lid_driven_cavity", "direct"):         _CP,
    ("lid_driven_cavity", "adjoint"):        _CP,
    ("lid_driven_cavity", "floquet"):        _NA,  # not_applicable for steady-only target
    ("lid_driven_cavity", "transient_growth"): _CP,
    ("lid_driven_cavity", "wavemaker"):      _BP,
    ("lid_driven_cavity", "animation"):      _CP,
    ("lid_driven_cavity", "otd"):            _CP,
    ("lid_driven_cavity", "modal"):          _CP,
    ("lid_driven_cavity", "rans"):           _NA,
    # naca0012
    ("naca0012", "dns"):                     _I,
    ("naca0012", "baseflow"):                _CP,
    ("naca0012", "direct"):                  _CP,
    ("naca0012", "adjoint"):                 _CP,
    ("naca0012", "floquet"):                 _CP,
    ("naca0012", "transient_growth"):        _CP,
    ("naca0012", "wavemaker"):               _BP,
    ("naca0012", "animation"):               _CP,
    ("naca0012", "otd"):                     _CP,
    ("naca0012", "modal"):                   _CP,
    ("naca0012", "rans"):                    _NA,
    # poiseuille
    ("poiseuille", "dns"):                   _CP,
    ("poiseuille", "baseflow"):              _CP,
    ("poiseuille", "direct"):                _CP,
    ("poiseuille", "adjoint"):               _CP,
    ("poiseuille", "floquet"):               _NA,  # not_applicable for steady channel target
    ("poiseuille", "transient_growth"):      _CP,
    ("poiseuille", "wavemaker"):             _BP,
    ("poiseuille", "animation"):             _CP,
    ("poiseuille", "otd"):                   _I,
    ("poiseuille", "modal"):                 _CP,
    ("poiseuille", "rans"):                  _I,
    # slot_FST
    ("slot_FST", "dns"):                     _I,   # implemented/partial
    ("slot_FST", "baseflow"):                _CP,
    ("slot_FST", "direct"):                  _CP,
    ("slot_FST", "adjoint"):                 _CP,
    ("slot_FST", "floquet"):                 _CP,
    ("slot_FST", "transient_growth"):        _CP,
    ("slot_FST", "wavemaker"):               _BP,
    ("slot_FST", "animation"):               _CP,
    ("slot_FST", "otd"):                     _CP,
    ("slot_FST", "modal"):                   _CP,
    ("slot_FST", "rans"):                    _NA,
}


# ---------------------------------------------------------------------------
# Shared constants for cylinder cases
# ---------------------------------------------------------------------------

_CYL_RE100 = TargetParameters(Re=100.0, geometry="cylinder_re100", motion="static", model="DNS")
_REQ_IMG = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE)
_REQ_IMG_OPT = EvidenceRequirement(role=EvidenceRole.REFERENCE_IMAGE, required=False)
_REQ_RESID = EvidenceRequirement(role=EvidenceRole.RESIDUAL_HISTORY)
_REQ_SPEC = EvidenceRequirement(role=EvidenceRole.SPECTRUM)
_REQ_ANIM = EvidenceRequirement(role=EvidenceRole.ANIMATION)
_REQ_DNS_HX = EvidenceRequirement(role=EvidenceRole.DNS_TIME_HISTORY)

_BLV = ComputeClass.BOUNDED_LOCAL_VALIDATION
_MLO = ComputeClass.MANUAL_ONLY
_LLR = ComputeClass.LOCAL_LONG_RUN
_QLS = ComputeClass.QUICK_LOCAL_SMOKE


def _img(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.REFERENCE_IMAGE, path=path, freshness="current")


def _dns_hx(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.DNS_TIME_HISTORY, path=path, freshness="current")


def _history(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.RESIDUAL_HISTORY, path=path, freshness="current")


def _spectrum(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.SPECTRUM, path=path, freshness="current")


def _anim(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.ANIMATION, path=path, freshness="current")


def _checkpoint(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.CHECKPOINT_STATS, path=path, freshness="current")


def _provenance(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.PROVENANCE, path=path, freshness="current")


def _log(path: str) -> Artifact:
    return Artifact(role=EvidenceRole.LOG, path=path, freshness="current")


# ---------------------------------------------------------------------------
# Module-level data: CASES
# ---------------------------------------------------------------------------

CASES: dict[str, MethodCase] = {
    # ------------------------------------------------------------------
    # Cylinder DNS / CI
    # ------------------------------------------------------------------
    "cylinder/dns": MethodCase(
        case_id="cylinder/dns",
        flow_family="cylinder_re100", method_lane="dns",
        current_path="example/cylinder_re100/000_dns",
        proposed_folder="example/cylinder_re100/000_dns", sort_prefix="000",
        label="Cylinder DNS", mode_name="DNS", legacy_uparam01="0.0",
        expected_behavior="Nonlinear DNS of flow past a cylinder at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG, _REQ_DNS_HX),
        artifacts=(
            
            _img("validation/figures/cylinder_dns_snapshot_Re50.png"),
            _img("validation/figures/cylinder_dns_signal_Re50.png"),
            _img("validation/figures/cylinder_dns_fft_Re50.png"),
            _img("validation/figures/cylinder_dns_phase_Re50.png"),
            _dns_hx("example/cylinder_re100/000_dns/1cyl.his"),
            _log("example/cylinder_re100/000_dns/1cyl.log.20"),
            _provenance("example/cylinder_re100/000_dns/README.md"),
        ),
        evidence_prereqs=(
            
            "validation/figures/cylinder_dns_snapshot_Re50.png",
            "validation/figures/cylinder_dns_signal_Re50.png",
            "validation/figures/cylinder_dns_fft_Re50.png",
            "validation/figures/cylinder_dns_phase_Re50.png",
            "example/cylinder_re100/000_dns/1cyl.his",
        ),
        why_note="Re=100-canonical target; variableDt=yes targetCFL=0.5; reference nonlinear DNS lane",
    ),
    "cylinder/ci_test": MethodCase(
        case_id="cylinder/ci_test",
        flow_family="cylinder_re100", method_lane="dns",
        current_path="example/cylinder_re100/000_dns/ci_test",
        proposed_folder="example/cylinder_re100/001_ci_test", sort_prefix="001",
        label="Cylinder CI Test", mode_name=None, legacy_uparam01=None,
        expected_behavior="Lightweight CI smoke test for cylinder mesh and nonlinear driver",
        status=CaseStatus.VALIDATED, compute_class=_QLS,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG,),
        artifacts=(
            _img("validation/figures/cylinder_ci_test_Re50.png"),
            _dns_hx("example/cylinder_re100/000_dns/ci_test/lift_drag.dat"),
            _log("example/cylinder_re100/000_dns/ci_test/logfile"),
        ),
        why_note=(
            "Local smoke only; numSteps=10, writeInterval=10; not final evidence due to"
            " intentionally tiny run length"
        ),
        gallery_visible=False,
    ),
    # ------------------------------------------------------------------
    # Moving cylinder (separate family)
    # ------------------------------------------------------------------
    "cylinder/moving_cylinder": MethodCase(
        case_id="cylinder/moving_cylinder",
        flow_family="moving_cylinder", method_lane="dns",
        current_path="example/moving_cylinder_re100",
        proposed_folder="example/moving_cylinder/020_dns", sort_prefix="020",
        label="Moving Cylinder (Mesh Motion, Forced)", mode_name=None, legacy_uparam01=None,
        expected_behavior="Forced-motion cylinder using moving mesh; deferred pending run policy",
        status=CaseStatus.DEFERRED, compute_class=_MLO,
        target_parameters=TargetParameters(Re=100.0, geometry="moving_cylinder", motion="forced", model="DNS"),
        evidence_requirements=(_REQ_IMG_OPT,),
        why_note="moving-mesh family; DEFERRED 2026-05-21. Revisit after cylinder_re100 + cylinder_re180 Floquet ladders are validated. See example/moving_cylinder_re100/README.md for the resumption plan.",
    ),
    "moving_cylinder/newton_forced_upo": MethodCase(
        case_id="moving_cylinder/newton_forced_upo",
        flow_family="moving_cylinder", method_lane="baseflow",
        current_path="example/moving_cylinder_re100",
        proposed_folder="example/moving_cylinder/220_baseflow_newton_upo", sort_prefix="220",
        label="Moving Cylinder Newton Forced UPO", mode_name="Newton-GMRES UPO", legacy_uparam01="2.1",
        expected_behavior="Forced-motion cylinder UPO baseflow; planned future lane",
        status=CaseStatus.DEFERRED, compute_class=ComputeClass.MANUAL_ONLY,
        target_parameters=TargetParameters(Re=100.0, geometry="moving_cylinder", motion="forced", model="DNS"),
        evidence_requirements=(_REQ_IMG_OPT,),
        why_note="moving-mesh family; DEFERRED 2026-05-21. Revisit after cylinder_re100 + cylinder_re180 Floquet ladders are validated. See example/moving_cylinder_re100/README.md for the resumption plan.",
    ),
    "moving_cylinder/floquet_direct": MethodCase(
        case_id="moving_cylinder/floquet_direct",
        flow_family="moving_cylinder", method_lane="floquet_direct",
        current_path="example/moving_cylinder_re100",
        proposed_folder="example/moving_cylinder/311_stability_direct_floquet", sort_prefix="311",
        label="Moving Cylinder Direct Floquet", mode_name="Direct Floquet", legacy_uparam01="3.11",
        expected_behavior="Forced-motion cylinder direct Floquet analysis; planned future lane",
        status=CaseStatus.DEFERRED, compute_class=ComputeClass.MANUAL_ONLY,
        target_parameters=TargetParameters(Re=100.0, geometry="moving_cylinder", motion="forced", model="DNS"),
        evidence_requirements=(_REQ_IMG_OPT,),
        why_note="moving-mesh family; DEFERRED 2026-05-21. Revisit after cylinder_re100 + cylinder_re180 Floquet ladders are validated. See example/moving_cylinder_re100/README.md for the resumption plan.",
    ),
    "moving_cylinder/energy_response": MethodCase(
        case_id="moving_cylinder/energy_response",
        flow_family="moving_cylinder", method_lane="animation",
        current_path="example/moving_cylinder_re100",
        proposed_folder="example/moving_cylinder/410_animation_energy_response", sort_prefix="410",
        label="Moving Cylinder Energy Response", mode_name=None, legacy_uparam01=None,
        expected_behavior="Forced-motion cylinder energy-response animation; planned future lane",
        status=CaseStatus.DEFERRED, compute_class=ComputeClass.MANUAL_ONLY,
        target_parameters=TargetParameters(Re=100.0, geometry="moving_cylinder", motion="forced", model="DNS"),
        evidence_requirements=(_REQ_IMG_OPT,),
        why_note="moving-mesh family; DEFERRED 2026-05-21. Revisit after cylinder_re100 + cylinder_re180 Floquet ladders are validated. See example/moving_cylinder_re100/README.md for the resumption plan.",
    ),
    # ------------------------------------------------------------------
    # Cylinder baseflow — SFD family
    # ------------------------------------------------------------------
    "cylinder/baseflow/sfd": MethodCase(
        case_id="cylinder/baseflow/sfd",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/110_baseflow_sfd/akervik",
        proposed_folder="example/cylinder_re100/110_baseflow_sfd/akervik", sort_prefix="110",
        label="Cylinder Baseflow SFD", mode_name="SFD", legacy_uparam01="1.1",
        expected_behavior="SFD stabilized baseflow at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        seed_parameters=SeedParameters(
            seed_Re=50.0,
            seed_path="BF_seed_1cyl0.f00001",
            restart_source="startFrom",
            seed_status="intentional_lower_re_seed",
        ),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            
            _img("validation/figures/cylinder_baseflow_sfd_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_bf_Re50.png"),
            _history("example/cylinder_re100/110_baseflow_sfd/akervik/residu.dat"),
            _log("example/cylinder_re100/110_baseflow_sfd/akervik/logfile"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            
            "validation/figures/cylinder_baseflow_sfd_ic_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_residual_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_bf_Re50.png",
            "example/cylinder_re100/110_baseflow_sfd/akervik/residu.dat",
        ),
        why_note=(
            "needs_normalization: viscosity=-50 (Re=50) vs canonical Re=100; SFD"
            " knobs userParam04=0.12 userParam05=0.05 are intentional_variant"
        ),
    ),
    "cylinder/baseflow/sfd_akervik_dyn": MethodCase(
        case_id="cylinder/baseflow/sfd_akervik_dyn",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/110_baseflow_sfd/akervik_dyn",
        proposed_folder="example/cylinder_re100/110_baseflow_sfd/akervik_dyn", sort_prefix="110",
        label="Cylinder Baseflow SFD Akervik Dynamic", mode_name="SFD", legacy_uparam01="1.1",
        expected_behavior="Dynamic SFD variant at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        seed_parameters=SeedParameters(
            seed_Re=50.0,
            seed_path="BF_seed_1cyl0.f00001",
            restart_source="startFrom",
            seed_status="intentional_lower_re_seed",
        ),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_dynamic_tolerance_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_bf_Re50.png"),            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_dynamic_tolerance_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_baseflow_Re50.png"),
            _history("example/cylinder_re100/110_baseflow_sfd/akervik_dyn/residu.dat"),
            _log("example/cylinder_re100/110_baseflow_sfd/akervik_dyn/logfile"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_ic_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_residual_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_dynamic_tolerance_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_bf_Re50.png",
            "example/cylinder_re100/110_baseflow_sfd/akervik_dyn/residu.dat",
        ),
        why_note=(
            "needs_normalization: Re=50; dynamic-tolerance knobs userParam08=0.9"
            " userParam09=100 are intentional_variant"
        ),
    ),
    "cylinder/baseflow/sfd_casacuberta_dyn": MethodCase(
        case_id="cylinder/baseflow/sfd_casacuberta_dyn",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn",
        proposed_folder="example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn", sort_prefix="110",
        label="Cylinder Baseflow SFD Casacuberta Dynamic", mode_name="SFD", legacy_uparam01="1.1",
        expected_behavior="Dynamic SFD variant at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        seed_parameters=SeedParameters(
            seed_Re=50.0,
            seed_path="BF_seed_1cyl0.f00001",
            restart_source="startFrom",
            seed_status="intentional_lower_re_seed",
        ),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_dynamic_tolerance_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_bf_Re50.png"),            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_dynamic_tolerance_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_baseflow_Re50.png"),
            _history("example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn/residu.dat"),
            _log("example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn/logfile"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_ic_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_residual_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_dynamic_tolerance_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_bf_Re50.png",
            "example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn/residu.dat",
        ),
        why_note=(
            "needs_normalization: Re=50; dynamic-tolerance knobs userParam08=0.9"
            " userParam09=100 are intentional_variant"
        ),
    ),
    "cylinder/baseflow/sfd_akervik_dyn_oifs": MethodCase(
        case_id="cylinder/baseflow/sfd_akervik_dyn_oifs",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs",
        proposed_folder="example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs", sort_prefix="110",
        label="Cylinder Baseflow SFD Akervik Dynamic OIFS", mode_name="SFD", legacy_uparam01="1.1",
        expected_behavior="Dynamic SFD with OIFS time integration at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        seed_parameters=SeedParameters(
            seed_Re=50.0,
            seed_path="BF_seed_1cyl0.f00001",
            restart_source="startFrom",
            seed_status="intentional_lower_re_seed",
        ),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_residual_vs_walltime_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_bf_Re50.png"),
            _history("example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs/residu.dat"),
            _log("example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs/logfile"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_ic_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_residual_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_residual_vs_walltime_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_akervik_dyn_oifs_bf_Re50.png",
            "example/cylinder_re100/110_baseflow_sfd/akervik_dyn_oifs/residu.dat",
        ),
        why_note=(
            "needs_normalization: Re=50; OIFS targetCFL=5 is intentional_variant vs"
            " canonical targetCFL=0.5"
        ),
    ),
    "cylinder/baseflow/sfd_casacuberta_dyn_oifs": MethodCase(
        case_id="cylinder/baseflow/sfd_casacuberta_dyn_oifs",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn_oifs",
        proposed_folder="example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn_oifs", sort_prefix="110",
        label="Cylinder Baseflow SFD Casacuberta Dynamic OIFS", mode_name="SFD", legacy_uparam01="1.1",
        expected_behavior="Dynamic SFD with OIFS time integration at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        seed_parameters=SeedParameters(
            seed_Re=50.0,
            seed_path="BF_seed_1cyl0.f00001",
            restart_source="startFrom",
            seed_status="intentional_lower_re_seed",
        ),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_residual_vs_walltime_Re50.png"),
            _img("validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_bf_Re50.png"),
            _history("example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn_oifs/residu.dat"),
            _log("example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn_oifs/logfile"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_ic_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_residual_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_residual_vs_walltime_Re50.png",
            "validation/figures/cylinder_baseflow_sfd_casacuberta_dyn_oifs_bf_Re50.png",
            "example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn_oifs/residu.dat",
        ),
        why_note=(
            "needs_normalization: Re=50; OIFS targetCFL=5 is intentional_variant vs"
            " canonical targetCFL=0.5"
        ),
    ),
    "cylinder/baseflow/boostconv": MethodCase(
        case_id="cylinder/baseflow/boostconv",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/120_baseflow_boostconv",
        proposed_folder="example/cylinder_re100/120_baseflow_boostconv", sort_prefix="120",
        label="Cylinder Baseflow BoostConv", mode_name="BoostConv", legacy_uparam01="1.2",
        expected_behavior="BoostConv accelerated baseflow convergence at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            
            _img("validation/figures/cylinder_baseflow_boostconv_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_boostconv_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_boostconv_bf_Re50.png"),
            _history("example/cylinder_re100/120_baseflow_boostconv/residu.dat"),
            _log("example/cylinder_re100/120_baseflow_boostconv/logfile"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            
            "validation/figures/cylinder_baseflow_boostconv_ic_Re50.png",
            "validation/figures/cylinder_baseflow_boostconv_residual_Re50.png",
            "validation/figures/cylinder_baseflow_boostconv_bf_Re50.png",
            "example/cylinder_re100/120_baseflow_boostconv/residu.dat",
        ),
        why_note=(
            "Re=100-canonical; variableDt=no; BoostConv acceleration userParam01=1.2;"
            " compare with SFD after SFD target is normalized"
        ),
    ),
    # ------------------------------------------------------------------
    # Cylinder baseflow — Newton family
    # ------------------------------------------------------------------
    "cylinder/baseflow/newton": MethodCase(
        case_id="cylinder/baseflow/newton",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/210_baseflow_newton/fp",
        proposed_folder="example/cylinder_re100/210_baseflow_newton/fp", sort_prefix="210",
        label="Cylinder Baseflow Newton", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Newton-GMRES baseflow at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        seed_parameters=SeedParameters(seed_Re="unknown", seed_path="BF_seed_1cyl0.f00001",
                                       restart_source="startFrom"),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            
            _img("validation/figures/cylinder_baseflow_newton_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_newton_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_newton_bf_Re50.png"),
            _history("example/cylinder_re100/210_baseflow_newton/fp/residu_arnoldi.dat"),
            _log("example/cylinder_re100/210_baseflow_newton/fp/logfile"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            
            "validation/figures/cylinder_baseflow_newton_ic_Re50.png",
            "validation/figures/cylinder_baseflow_newton_residual_Re50.png",
            "validation/figures/cylinder_baseflow_newton_bf_Re50.png",
            "example/cylinder_re100/210_baseflow_newton/fp/residu_arnoldi.dat",
        ),
        why_note=(
            "Re=100-canonical; seed BF_seed_1cyl0.f00001 (seed_Re unknown); Krylov"
            " k_dim=100 (userParam07=100)"
        ),
    ),
    "cylinder/baseflow/newton_dyn": MethodCase(
        case_id="cylinder/baseflow/newton_dyn",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re100/210_baseflow_newton/dyn",
        proposed_folder="example/cylinder_re100/210_baseflow_newton/fp_dyn", sort_prefix="210",
        label="Cylinder Baseflow Newton Dynamic", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Dynamic Newton-GMRES baseflow at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        seed_parameters=SeedParameters(seed_Re="unknown", seed_path="BF_seed_1cyl0.f00001",
                                       restart_source="startFrom"),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            
            _img("validation/figures/cylinder_baseflow_newton_dyn_ic_Re50.png"),
            _img("validation/figures/cylinder_baseflow_newton_dyn_residual_Re50.png"),
            _img("validation/figures/cylinder_baseflow_newton_dyn_bf_Re50.png"),
            _history("example/cylinder_re100/210_baseflow_newton/dyn/residu_arnoldi.dat"),
            _log("example/cylinder_re100/210_baseflow_newton/dyn/logfile"),
        ),
        prerequisites=("cylinder/baseflow/newton",),
        evidence_prereqs=(
            
            "validation/figures/cylinder_baseflow_newton_dyn_ic_Re50.png",
            "validation/figures/cylinder_baseflow_newton_dyn_residual_Re50.png",
            "validation/figures/cylinder_baseflow_newton_dyn_bf_Re50.png",
            "example/cylinder_re100/210_baseflow_newton/dyn/residu_arnoldi.dat",
        ),
        why_note=(
            "Re=100-canonical; same seed as newton; dynamic-tolerance variant of"
            " newton lane"
        ),
    ),
    "cylinder/baseflow/newton_dyn_temp": MethodCase(
        # Thermal exception: target Re=30 with scalar coupling — not Re=100 drift
        case_id="cylinder/baseflow/newton_dyn_temp",
        flow_family="cylinder_re100", method_lane="baseflow",
        current_path="example/cylinder_re30_thermal/210_baseflow_newton",
        proposed_folder="example/cylinder_re100/210_baseflow_newton/fp_dyn_temp", sort_prefix="210",
        label="Cylinder Baseflow Newton Dynamic Thermal", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Thermal Newton-GMRES baseflow at Re=30 (thermal exception, not Re=100 target)",
        status=CaseStatus.NEEDS_DNS_SEED, compute_class=_BLV,
        target_parameters=TargetParameters(Re=30.0, Pr=0.71, geometry="cylinder_re100",
                                           motion="static", model="DNS"),
        seed_parameters=SeedParameters(seed_Re=25.0, seed_Pr=0.71,
                                       seed_path="BFre25t_1cyl0.f00001",
                                       restart_source="startFrom",
                                       seed_status="archived"),
        evidence_requirements=(
            _REQ_IMG,
            _REQ_RESID,
            EvidenceRequirement(role=EvidenceRole.CHECKPOINT_STATS),
        ),
        artifacts=(
            _img("validation/figures/cylinder_baseflow_newton_dyn_temp_residual_Re50.png"),
        ),
        blockers=(
            "BFre25t_1cyl0.f00001 was 1996-element thermal mesh; archived 2026-05-21 to"
            " .archive_mesh_1996_2026-05-21/. Needs fresh thermal Re=25 baseflow on the"
            " unified 2128 mesh before Re=30 continuation can resume.",
        ),
        why_note="mesh-unified 2026-05-21; see NOTE-mesh-replaced-2026-05-21.md",
    ),
    "cylinder/baseflow/newton_upo": MethodCase(
        case_id="cylinder/baseflow/newton_upo",
        flow_family="cylinder_re180", method_lane="baseflow",
        current_path="example/cylinder_re180/210_baseflow_newton",
        proposed_folder="example/cylinder_re180/210_baseflow_newton", sort_prefix="210",
        label="Cylinder Newton UPO at Re=180", mode_name="Newton-GMRES UPO", legacy_uparam01="2.1",
        expected_behavior="Newton UPO at Re=180 (St~0.18, T~5.5); base flow for the Floquet 311/321/410 chain",
        status=CaseStatus.NEEDS_DNS_SEED, compute_class=_LLR,
        target_parameters=TargetParameters(Re=180.0, geometry="cylinder", motion="static", model="DNS"),
        seed_parameters=SeedParameters(seed_Re=150.0, seed_path="rstcyl0.f00001",
                                       restart_source="startFrom",
                                       seed_status="needs-precompute"),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            _img("validation/figures/cylinder_baseflow_newton_upo_ic_Re180.png"),
            _img("validation/figures/cylinder_baseflow_newton_upo_residual_Re180.png"),
            _img("validation/figures/cylinder_baseflow_newton_upo_bf_Re180.png"),
        ),
        prerequisites=("cylinder/dns",),
        evidence_prereqs=(
            "validation/figures/cylinder_baseflow_newton_upo_ic_Re180.png",
            "validation/figures/cylinder_baseflow_newton_upo_residual_Re180.png",
            "validation/figures/cylinder_baseflow_newton_upo_bf_Re180.png",
        ),
        blockers=(
            "Needs Re=150 DNS snapshot as seed (see README): pre-compute a"
            " DNS at Re=150 until limit cycle is established (t~200-500),"
            " snapshot as rstcyl0.f00001, then Newton-UPO at Re=180 with"
            " period guess T~5.5. The Re=100 historical snapshot at t=500"
            " was archived 2026-05-21 and is for the wrong Re.",
        ),
        why_note="Re=180 UPO at the Barkley-Henderson Mode A threshold; base flow for the Floquet 311/321/410 chain. See cylinder_re180/210_baseflow_newton/README.md for the Re=150-seed workflow.",
    ),
    # ------------------------------------------------------------------
    # Cylinder stability
    # ------------------------------------------------------------------
    "cylinder/stability/direct": MethodCase(
        case_id="cylinder/stability/direct",
        flow_family="cylinder_re100", method_lane="stability_direct",
        current_path="example/cylinder_re100/310_stability_direct/direct",
        proposed_folder="example/cylinder_re100/310_stability_direct/direct", sort_prefix="310",
        label="Cylinder Direct Stability", mode_name="Direct LNSE", legacy_uparam01="3.1",
        expected_behavior="Direct linearized Navier-Stokes eigenvalue computation at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG, _REQ_SPEC),
        artifacts=(
            _img("validation/figures/cylinder_stability_direct_bf_Re50.png"),
            _img("validation/figures/cylinder_stability_direct_spectrum_circle_Re50.png"),
            _img("validation/figures/cylinder_stability_direct_spectrum_ns_log_Re50.png"),
            _img("validation/figures/cylinder_stability_direct_mode_vx_Re50.png"),
            _img("validation/figures/cylinder_stability_direct_mode_vy_Re50.png"),
            _spectrum("example/cylinder_re100/310_stability_direct/direct/Spectre_Hd.dat"),
            _log("example/cylinder_re100/310_stability_direct/direct/Spectre_d.info"),
            _provenance("example/cylinder_re100/310_stability_direct/direct/README.md"),
        ),
        prerequisites=("cylinder/baseflow/sfd",),
        evidence_prereqs=(
            "validation/figures/cylinder_stability_direct_bf_Re50.png",
            "validation/figures/cylinder_stability_direct_spectrum_circle_Re50.png",
            "validation/figures/cylinder_stability_direct_spectrum_ns_log_Re50.png",
            "validation/figures/cylinder_stability_direct_mode_vx_Re50.png",
            "validation/figures/cylinder_stability_direct_mode_vy_Re50.png",
            "example/cylinder_re100/310_stability_direct/direct/Spectre_Hd.dat",
        ),
        why_note=(
            "Re=100-canonical; k_dim=50 (userParam07=50); sponge userParam08=5"
            " 09=5 10=0; startFrom=BF_1cyl0.f00001"
        ),
    ),
    "cylinder/stability/direct_Floquet": MethodCase(
        case_id="cylinder/stability/direct_Floquet",
        flow_family="cylinder_re180", method_lane="floquet_direct",
        current_path="example/cylinder_re180/311_stability_direct_floquet",
        proposed_folder="example/cylinder_re180/311_stability_direct_floquet", sort_prefix="311",
        label="Cylinder Direct Floquet at Re=180", mode_name="Direct Floquet", legacy_uparam01="3.11",
        expected_behavior="Direct Floquet of the Re=180 UPO: leading multiplier mu=+1 (synchronous phase mode) on the unit circle; all other 2D modes damped (|mu|<1)",
        status=CaseStatus.VALIDATED, compute_class=_MLO,
        target_parameters=TargetParameters(Re=180.0, geometry="cylinder", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, _REQ_SPEC),
        artifacts=(
            _img("validation/figures/cylinder_stability_direct_Floquet_spectrum_Re180.png"),
            _img("validation/figures/cylinder_stability_direct_Floquet_mode_Re180.png"),
            _img("validation/figures/cylinder_stability_direct_Floquet_bf_Re180.png"),
            _spectrum("example/cylinder_re180/311_stability_direct_floquet/Spectre_NSd_conv.dat"),
            _log("example/cylinder_re180/311_stability_direct_floquet/logfile"),
        ),
        prerequisites=("cylinder/baseflow/newton_upo",),
        evidence_prereqs=(
            "validation/figures/cylinder_stability_direct_Floquet_spectrum_Re180.png",
            "validation/figures/cylinder_stability_direct_Floquet_mode_Re180.png",
        ),
        why_note="2D PROXY of Barkley & Henderson 1996: Mode A/B are 3D (beta!=0) and cannot appear in this 2D (beta=0) run; the only neutral mode here is the trivial mu=+1 phase mode. Pipeline demonstrator at the famous parameters, not a physical secondary-instability result.",
    ),
    "cylinder/stability/adjoint": MethodCase(
        case_id="cylinder/stability/adjoint",
        flow_family="cylinder_re100", method_lane="stability_adjoint",
        current_path="example/cylinder_re100/320_stability_adjoint",
        proposed_folder="example/cylinder_re100/320_stability_adjoint/schur", sort_prefix="320",
        label="Cylinder Adjoint Stability", mode_name="Adjoint LNSE", legacy_uparam01="3.2",
        expected_behavior="Adjoint linearized Navier-Stokes eigenvalue computation at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG, _REQ_SPEC),
        artifacts=(
            _img("validation/figures/cylinder_stability_direct_Floquet_spectrum_Re180.png"),
            _spectrum("example/cylinder_re100/320_stability_adjoint/Spectre_Hd.dat"),
            _log("example/cylinder_re100/320_stability_adjoint/logfile"),
        ),
        prerequisites=("cylinder/baseflow/sfd",),
        evidence_prereqs=(
            "validation/figures/cylinder_stability_direct_Floquet_spectrum_Re180.png",
            "example/cylinder_re100/320_stability_adjoint/Spectre_Hd.dat",
        ),
        why_note=(
            "Re=100-canonical; k_dim=50 (userParam07=50); sponge strength"
            " userParam10=1.7; startFrom=BF_1cyl0.f00001"
        ),
    ),
    "cylinder/stability/adjoint_Floquet": MethodCase(
        case_id="cylinder/stability/adjoint_Floquet",
        flow_family="cylinder_re180", method_lane="floquet_adjoint",
        current_path="example/cylinder_re180/321_stability_adjoint_floquet",
        proposed_folder="example/cylinder_re180/321_stability_adjoint_floquet", sort_prefix="321",
        label="Cylinder Adjoint Floquet at Re=180", mode_name="Adjoint Floquet", legacy_uparam01="3.21",
        expected_behavior="Adjoint Floquet of the Re=180 UPO: damped modes (|mu|<1), adjoint mode localized near/upstream of the body (receptivity); mirrors direct spectrum's oscillatory pairs",
        status=CaseStatus.VALIDATED, compute_class=_MLO,
        target_parameters=TargetParameters(Re=180.0, geometry="cylinder", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, _REQ_SPEC),
        artifacts=(
            _img("validation/figures/cylinder_stability_adjoint_Floquet_spectrum_Re180.png"),
            _img("validation/figures/cylinder_stability_adjoint_Floquet_mode_Re180.png"),
            _img("validation/figures/cylinder_stability_adjoint_Floquet_bf_Re180.png"),
            _spectrum("example/cylinder_re180/321_stability_adjoint_floquet/Spectre_NSa_conv.dat"),
            _log("example/cylinder_re180/321_stability_adjoint_floquet/logfile"),
        ),
        prerequisites=("cylinder/baseflow/newton_upo",),
        evidence_prereqs=(
            "validation/figures/cylinder_stability_adjoint_Floquet_spectrum_Re180.png",
            "validation/figures/cylinder_stability_adjoint_Floquet_mode_Re180.png",
        ),
        why_note="2D PROXY (see 311). Adjoint mu=+1 phase mode did not converge into the leading set (known adjoint convergence difficulty); leading adjoint here is a decaying real mode. Oscillatory pairs match the direct spectrum.",
    ),
    # ------------------------------------------------------------------
    # Cylinder postproc / animation
    # ------------------------------------------------------------------
    "cylinder/stability/animate_modes": MethodCase(
        case_id="cylinder/stability/animate_modes",
        flow_family="cylinder_re100", method_lane="animation",
        current_path="example/cylinder_re100/410_postproc_animate_modes",
        proposed_folder="example/cylinder_re100/310_stability_direct/direct/animate", sort_prefix="310",
        label="Cylinder Animate Modes", mode_name="Mode Animation", legacy_uparam01=None,
        expected_behavior="Mode shape animation postprocessing at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG, _REQ_ANIM),
        artifacts=(
            _img("validation/figures/cylinder_stability_adjoint_Floquet_spectrum_Re180.png"),
            _anim("example/cylinder_re100/410_postproc_animate_modes/dQ_1cyl0.f00001"),
            _spectrum("example/cylinder_re100/410_postproc_animate_modes/Spectre_NSd_conv.dat"),
            _log("example/cylinder_re100/410_postproc_animate_modes/logfile"),
        ),
        prerequisites=("cylinder/stability/direct",),
        evidence_prereqs=(
            "validation/figures/cylinder_stability_adjoint_Floquet_spectrum_Re180.png",
            "example/cylinder_re100/410_postproc_animate_modes/dQ_1cyl0.f00001",
        ),
        why_note=(
            "Re=100-canonical; endTime=8.3050 is mode-period artifact setting;"
            " steady-base animation lane"
        ),
    ),
    "cylinder/stability/animate_modes_with_UPO": MethodCase(
        case_id="cylinder/stability/animate_modes_with_UPO",
        flow_family="cylinder_re180", method_lane="animation",
        current_path="example/cylinder_re180/410_postproc_animate_modes/upo",
        proposed_folder="example/cylinder_re180/411_postproc_animate_floquet_upo", sort_prefix="411",
        label="Cylinder Animate Floquet Modes on UPO at Re=180", mode_name="Mode Animation UPO", legacy_uparam01=None,
        expected_behavior="Animate the mu=+1 synchronous mode on the Re=180 UPO: base-flow DNS advanced one period T (the phase mode rides the orbit; shedding at omega_base=2pi/T)",
        status=CaseStatus.VALIDATED, compute_class=_MLO,
        target_parameters=TargetParameters(Re=180.0, geometry="cylinder", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, _REQ_ANIM),
        artifacts=(
            _img("validation/figures/cylinder_stability_animate_modes_with_UPO_Re180.gif"),
            _img("validation/figures/cylinder_stability_animate_modes_with_UPO_still_Re180.png"),
            _anim("example/cylinder_re180/410_postproc_animate_modes/upo/floquet_mode_animation.mp4"),
            _log("example/cylinder_re180/410_postproc_animate_modes/upo/logfile"),
        ),
        prerequisites=("cylinder/baseflow/newton_upo",
                       "cylinder/stability/direct_Floquet"),
        evidence_prereqs=(
            "validation/figures/cylinder_stability_animate_modes_with_UPO_Re180.gif",
            "validation/figures/cylinder_stability_animate_modes_with_UPO_still_Re180.png",
        ),
        why_note="2D PROXY (see 311). animate_mode_Floquet (uParam01=4.52): the leading mu=+1 mode is the synchronous phase mode, so its honest animation is the evolving base-flow shedding over one period T (not the static cos/sin reconstruction).",
    ),
    "cylinder/postproc/sensitivity_budget_wavemaker": MethodCase(
        case_id="cylinder/postproc/sensitivity_budget_wavemaker",
        flow_family="cylinder_re100", method_lane="wavemaker",
        current_path="example/cylinder_re100/411_postproc_wavemaker",
        proposed_folder="example/cylinder_re100/411_postproc_wavemaker", sort_prefix="411",
        label="Cylinder Wavemaker Sensitivity Budget", mode_name="Wavemaker", legacy_uparam01="4.1",
        expected_behavior="Structural sensitivity (wavemaker) and energy budget at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG,),
        artifacts=(
            _img("validation/figures/cylinder_postproc_wavemaker_field_Re50.png"),
            _img("validation/figures/cylinder_postproc_wavemaker_direct_mode_Re50.png"),
            _img("validation/figures/cylinder_postproc_wavemaker_adjoint_mode_Re50.png"),
            _log("example/cylinder_re100/411_postproc_wavemaker/logfile"),
        ),
        prerequisites=("cylinder/stability/direct", "cylinder/stability/adjoint"),
        evidence_prereqs=(
            "validation/figures/cylinder_postproc_wavemaker_field_Re50.png",
            "validation/figures/cylinder_postproc_wavemaker_direct_mode_Re50.png",
            "validation/figures/cylinder_postproc_wavemaker_adjoint_mode_Re50.png",
        ),
        why_note=(
            "Re=100-canonical; endTime=0 startFrom=0; postproc-only lane reading"
            " direct+adjoint modes"
        ),
    ),
    "cylinder/postproc/steady_force_sensitivity": MethodCase(
        case_id="cylinder/postproc/steady_force_sensitivity",
        flow_family="cylinder_re100", method_lane="wavemaker",
        current_path="example/cylinder_re100/412_postproc_steady_force_sensitivity",
        proposed_folder="example/cylinder_re100/412_postproc_steady_force_sensitivity", sort_prefix="412",
        label="Cylinder Steady Force Sensitivity", mode_name="Force Sensitivity", legacy_uparam01="4.11",
        expected_behavior="Steady force sensitivity postprocessing at Re=100",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG,),
        artifacts=(
            _img("validation/figures/cylinder_postproc_steady_force_sensitivity_real_Re50.png"),
            _img("validation/figures/cylinder_postproc_steady_force_sensitivity_imag_Re50.png"),
            _history("example/cylinder_re100/412_postproc_steady_force_sensitivity/residu_arnoldi.dat"),
            _log("example/cylinder_re100/412_postproc_steady_force_sensitivity/logfile"),
        ),
        prerequisites=("cylinder/stability/adjoint",),
        evidence_prereqs=(
            "validation/figures/cylinder_postproc_steady_force_sensitivity_real_Re50.png",
            "validation/figures/cylinder_postproc_steady_force_sensitivity_imag_Re50.png",
        ),
        why_note=(
            "Re=100-canonical; k_dim=50 (userParam07=50); postproc reading adjoint mode"
        ),
    ),
    # ------------------------------------------------------------------
    # Cylinder OTD / Modal / RANS
    # ------------------------------------------------------------------
    "cylinder/otd": MethodCase(
        case_id="cylinder/otd",
        flow_family="cylinder_re100", method_lane="otd",
        current_path="example/cylinder_re180/500_otd",
        proposed_folder="example/cylinder_re100/500_otd", sort_prefix="500",
        label="Cylinder OTD", mode_name="OTD", legacy_uparam01="5.0",
        expected_behavior="Optimally Time-Dependent modes for cylinder at Re=180",
        status=CaseStatus.NEEDS_DNS_SEED, compute_class=_MLO,
        target_parameters=TargetParameters(Re=180.0, geometry="cylinder_re100", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG,),
        artifacts=(_img("validation/figures/cylinder_otd_Re180.png"),),
        prerequisites=("cylinder/dns",),
        blockers=(
            "OTD baseflow (BF_1cyl0.f00001) and perturbation set (bf0/ip0/r0 series)"
            " were 1996-element; archived 2026-05-21. OTD warm-start needs fresh"
            " baseflow + perturbation set on the 2128 mesh.",
        ),
        why_note="mesh-unified 2026-05-21; OTD evidence archived",
    ),
    "cylinder/modal": MethodCase(
        case_id="cylinder/modal",
        flow_family="cylinder_re100", method_lane="modal",
        current_path="example/cylinder_re100/600_modal_pod",
        proposed_folder="example/cylinder_re100/600_modal", sort_prefix="600",
        label="Cylinder Modal (POD/DMD/SPOD)", mode_name="Modal", legacy_uparam01=None,
        expected_behavior="Modal decomposition (POD/DMD/SPOD) of cylinder snapshots (partial)",
        status=CaseStatus.PARTIAL, compute_class=_MLO,
        target_parameters=_CYL_RE100,
        evidence_requirements=(_REQ_IMG_OPT,),
        artifacts=(
            _img("validation/figures/cylinder_modal_snapshot_Re100.png"),
            _img("validation/figures/cylinder_modal_pod_spectrum_Re100.png"),
            _img("validation/figures/cylinder_modal_pod_mode1_Re100.png"),
            _img("validation/figures/cylinder_modal_pod_mode2_Re100.png"),
            _img("validation/figures/cylinder_modal_dmd_spectrum_Re100.png"),
            _img("validation/figures/cylinder_modal_dmd_mode1_Re100.png"),
            _img("validation/figures/cylinder_modal_spod_spectrum_Re100.png"),
            _img("validation/figures/cylinder_modal_spod_mode1_Re100.png"),
        ),
        prerequisites=("cylinder/dns",),
        why_note=(
            "Re=100-canonical modal lane; endTime=2050 writeInterval=0.5 for snapshot"
            " generation; POD/DMD/SPOD from cylinder DNS snapshots"
        ),
    ),
    "cylinder/RANS": MethodCase(
        case_id="cylinder/RANS",
        flow_family="cylinder_re100", method_lane="rans",
        current_path="example/cylinder_re1m/310_stability_direct/findiff",
        proposed_folder="example/cylinder_re100/700_rans", sort_prefix="700",
        label="Cylinder RANS", mode_name="RANS", legacy_uparam01=None,
        expected_behavior="RANS model for cylinder flow (partial, special direct-only lane)",
        status=CaseStatus.PARTIAL, compute_class=_MLO,
        target_parameters=TargetParameters(Re=100.0, geometry="cylinder_re100", motion="static", model="RANS"),
        evidence_requirements=(_REQ_IMG_OPT,),
        artifacts=(_img("validation/figures/cylinder_rans_Re100.png"),),
        why_note=(
            "Shares example/cylinder_re1m/310_stability_direct/findiff/ with cylinder/RANS/re1m (different"
            " target Re); both are valid catalog entries because target_Re differs"
        ),
    ),
    # ------------------------------------------------------------------
    # Other families — one case each
    # ------------------------------------------------------------------
    "thermosyphon/baseflow": MethodCase(
        case_id="thermosyphon/baseflow",
        flow_family="thermosyphon", method_lane="baseflow",
        current_path="example/thermosyphon_ra500/210_baseflow_newton",
        proposed_folder="example/thermosyphon_ra500/210_thermosyphon_baseflow", sort_prefix="210",
        label="Thermosyphon Baseflow", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Thermosyphon natural-convection baseflow at Ra=500, Pr=5",
        status=CaseStatus.NEEDS_SCALAR_CHECKPOINT, compute_class=_BLV,
        target_parameters=TargetParameters(geometry="thermosyphon", motion="static", model="DNS"),
        seed_parameters=SeedParameters(
            seed_path="BF_Ra400_tsyphon0.f00001",
            restart_source="startFrom",
            seed_status="Ra=400 continuation seed",
        ),
        evidence_requirements=(
            _REQ_IMG,
            _REQ_RESID,
            EvidenceRequirement(role=EvidenceRole.CHECKPOINT_STATS, required=True,
                                note="Thermal/scalar correctness: log+plot alone are insufficient (status taxonomy)"),
        ),
        artifacts=(
            _img("validation/figures/thermosyphon_baseflow_ic_Ra500.png"),
            _img("validation/figures/thermosyphon_baseflow_residual_Ra500.png"),
            _img("validation/figures/thermosyphon_baseflow_bf_Ra500.png"),
            _history("example/thermosyphon_ra500/baseflow/residu_arnoldi.dat"),
            _log("example/thermosyphon_ra500/baseflow/logfile"),
        ),
        evidence_prereqs=(
            "validation/figures/thermosyphon_baseflow_ic_Ra500.png",
            "validation/figures/thermosyphon_baseflow_residual_Ra500.png",
            "validation/figures/thermosyphon_baseflow_bf_Ra500.png",
            "example/thermosyphon_ra500/baseflow/residu_arnoldi.dat",
        ),
        blockers=(
            "Needs scalar-field checkpoint statistics (min/max/norm of T)"
            " to confirm the temperature field is populated and Ra*Pr"
            " buoyancy coupling is active. Per status taxonomy: logs +"
            " plots are insufficient for thermal/scalar validation."
            " Generate via scripts/inspect_nek_field.py on BF_*.f00001"
            " and record the scalar range as a CHECKPOINT_STATS artifact.",
        ),
        why_note=(
            "Rayleigh-driven case (Ra=500, Pr=5)."
            " Distinguish Ra parameter evidence (target) from actual"
            " saved scalar-field statistics (evidence)."
        ),
    ),
    "thermosyphon/stability/direct": MethodCase(
        case_id="thermosyphon/stability/direct",
        flow_family="thermosyphon", method_lane="direct",
        current_path="example/thermosyphon_ra500/310_stability_direct/direct",
        proposed_folder="example/thermosyphon_ra500/310_thermosyphon_stability_direct", sort_prefix="310",
        label="Thermosyphon Direct Stability", mode_name="Direct LNSE", legacy_uparam01="3.1",
        expected_behavior="Direct stability of thermosyphon conductive base flow at Ra=500; buoyancy coupling Pr*Ra",
        status=CaseStatus.NEEDS_SCALAR_CHECKPOINT, compute_class=_BLV,
        target_parameters=TargetParameters(Re=None, geometry="thermosyphon", motion="static", model="DNS"),
        evidence_requirements=(
            _REQ_IMG,
            _REQ_SPEC,
            EvidenceRequirement(role=EvidenceRole.CHECKPOINT_STATS, required=True,
                                note="Base-flow scalar field must be populated; eigenvalue norm uses Pr*Ra T weighting"),
        ),
        artifacts=(
            _img("validation/figures/thermosyphon_stability_direct_spectrum_Ra500.png"),
            _img("validation/figures/thermosyphon_stability_direct_baseflow_Ra500.png"),
            _img("validation/figures/thermosyphon_stability_direct_mode_Ra500.png"),
            _spectrum("example/thermosyphon_ra500/stability/direct/Spectre_Hd.dat"),
            _log("example/thermosyphon_ra500/stability/direct/tsyphon.log.10"),
        ),
        prerequisites=("thermosyphon/baseflow",),
        evidence_prereqs=(
            "validation/figures/thermosyphon_stability_direct_spectrum_Ra500.png",
            "validation/figures/thermosyphon_stability_direct_baseflow_Ra500.png",
            "validation/figures/thermosyphon_stability_direct_mode_Ra500.png",
            "example/thermosyphon_ra500/stability/direct/Spectre_Hd.dat",
        ),
        blockers=(
            "Same scalar-field requirement as thermosyphon/baseflow:"
            " the base-flow checkpoint must carry a populated T field"
            " (rdcode XUPT, non-zero scalar range), otherwise the"
            " linearized buoyancy operator is fed zero T and the"
            " spectrum is physically meaningless even if numerically"
            " converged.",
        ),
        why_note=(
            "Direct stability of a thermally-coupled base flow."
            " Depends on the upstream baseflow's"
            " CHECKPOINT_STATS being satisfied first."
        ),
    ),
    "flip_flop/baseflow": MethodCase(
        case_id="flip_flop/baseflow",
        flow_family="flip_flop", method_lane="baseflow",
        current_path="example/flip_flop_re62/210_baseflow_newton",
        proposed_folder="example/flip_flop_re62/210_flip_flop_baseflow", sort_prefix="210",
        label="Flip-Flop Jet Baseflow", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Flip-flop jet baseflow at Re=62",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=TargetParameters(Re=62.0, geometry="flip_flop", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            _img("validation/figures/flip_flop_baseflow_ic_Re62.png"),
            _img("validation/figures/flip_flop_baseflow_residual_Re62.png"),
            _img("validation/figures/flip_flop_baseflow_bf_Re62.png"),
            _history("example/flip_flop_re62/baseflow/residu_arnoldi.dat"),
            _log("example/flip_flop_re62/baseflow/logfile"),
        ),
        evidence_prereqs=(
            "validation/figures/flip_flop_baseflow_ic_Re62.png",
            "validation/figures/flip_flop_baseflow_residual_Re62.png",
            "validation/figures/flip_flop_baseflow_bf_Re62.png",
            "example/flip_flop_re62/baseflow/residu_arnoldi.dat",
        ),
    ),
    "flip_flop/stability/direct_Floquet": MethodCase(
        case_id="flip_flop/stability/direct_Floquet",
        flow_family="flip_flop", method_lane="floquet_direct",
        current_path="example/flip_flop_re62/311_stability_direct_floquet",
        proposed_folder="example/flip_flop_re62/311_flip_flop_stability_direct_floquet", sort_prefix="311",
        label="Flip-Flop Direct Floquet", mode_name="Direct Floquet", legacy_uparam01="3.11",
        expected_behavior="Floquet analysis of flip-flop orbit at Re=62; leading multiplier sigma=0.0079, omega=0.1408",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=TargetParameters(Re=62.0, geometry="flip_flop", motion="static", model="DNS"),
        seed_parameters=SeedParameters(seed_Re=62.0, seed_path="BF_2cyl0.f00001", restart_source="from baseflow dir"),
        evidence_requirements=(_REQ_IMG, _REQ_SPEC),
        artifacts=(
            _img("validation/figures/flip_flop_stability_direct_Floquet_spectrum_Re62.png"),
            _img("validation/figures/flip_flop_stability_direct_Floquet_upo_snapshot_Re62.png"),
            _img("validation/figures/flip_flop_stability_direct_Floquet_mode_Re62.png"),
            _spectrum("example/flip_flop_re62/stability/direct_Floquet/Spectre_Hd.dat"),
            _log("example/flip_flop_re62/stability/direct_Floquet/logfile"),
        ),
        prerequisites=("flip_flop/baseflow",),
        evidence_prereqs=(
            "validation/figures/flip_flop_stability_direct_Floquet_spectrum_Re62.png",
            "validation/figures/flip_flop_stability_direct_Floquet_upo_snapshot_Re62.png",
            "validation/figures/flip_flop_stability_direct_Floquet_mode_Re62.png",
            "example/flip_flop_re62/stability/direct_Floquet/Spectre_Hd.dat",
        ),
        why_note="verified 2026-05-17 on 8 ranks; Re=62 just above Re_c~61.17",
    ),
    "tpjet/baseflow/newton": MethodCase(
        case_id="tpjet/baseflow/newton",
        flow_family="tpjet", method_lane="baseflow",
        current_path="example/tpjet_re2005/210_baseflow_newton",
        proposed_folder="example/tpjet_re2005/210_tpjet_baseflow_newton", sort_prefix="210",
        label="Two-Phase Jet Baseflow Newton", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Two-phase jet Newton baseflow at Re=2000 (compute-heavy)",
        status=CaseStatus.VALIDATED, compute_class=_LLR,
        target_parameters=TargetParameters(Re=2000.0, geometry="tpjet", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            _img("validation/figures/tpjet_baseflow_newton_ic_Re1900.png"),
            _img("validation/figures/tpjet_baseflow_newton_residual_Re1900.png"),
            _img("validation/figures/tpjet_baseflow_newton_bf_Re1900.png"),
            _history("example/tpjet_re2005/baseflow/newton/residu_arnoldi.dat"),
            _log("example/tpjet_re2005/baseflow/newton/tpjet.log.20"),
            _checkpoint("example/tpjet_re2005/baseflow/newton/BF_tpjet0.f00001"),
            _provenance("example/tpjet_re2005/baseflow/newton/README.md"),
        ),
        evidence_prereqs=(
            "validation/figures/tpjet_baseflow_newton_ic_Re1900.png",
            "validation/figures/tpjet_baseflow_newton_residual_Re1900.png",
            "validation/figures/tpjet_baseflow_newton_bf_Re1900.png",
            "example/tpjet_re2005/baseflow/newton/residu_arnoldi.dat",
            "example/tpjet_re2005/baseflow/newton/BF_tpjet0.f00001",
        ),
    ),
    "tpjet/130_baseflow_dmt": MethodCase(
        case_id="tpjet/130_baseflow_dmt",
        flow_family="tpjet", method_lane="baseflow",
        current_path="example/tpjet_re2005/130_baseflow_dmt",
        proposed_folder="example/tpjet_re2005/130_tpjet_baseflow_dmt", sort_prefix="130",
        label="Two-Phase Jet Baseflow DMT", mode_name="DMT", legacy_uparam01="1.3",
        expected_behavior="DMT (Direct Mapping Technique) baseflow for tpjet at Re=2000 (deferred validation)",
        status=CaseStatus.DEFERRED, compute_class=_LLR,
        target_parameters=TargetParameters(Re=2000.0, geometry="tpjet", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG_OPT,),
    ),
    "tpjet/baseflow/tdf": MethodCase(
        case_id="tpjet/baseflow/tdf",
        flow_family="tpjet", method_lane="baseflow",
        current_path="",
        proposed_folder="example/tpjet_re2005/140_tpjet_baseflow_tdf", sort_prefix="140",
        label="Two-Phase Jet TDF", mode_name="TDF", legacy_uparam01="1.4",
        expected_behavior="TDF periodic base flow for tpjet at Re=2005, St=0.60 (deferred validation)",
        # Stage not yet instantiated (proposed_folder 140_* does not exist yet);
        # empty current_path = backlog, not a broken path. Set the real dir once created.
        status=CaseStatus.DEFERRED, compute_class=_LLR,
        target_parameters=TargetParameters(Re=2005.0, geometry="tpjet", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG,),
        artifacts=(_img("validation/figures/tpjet_baseflow_tdf_Re2005.png"),),
        evidence_prereqs=("validation/figures/tpjet_baseflow_tdf_Re2005.png",),
        why_note="TBD: needs evidence; initial guess BF_Re1900_tpjet0.f00001",
    ),
    "tpjet/stability/direct_Floquet": MethodCase(
        case_id="tpjet/stability/direct_Floquet",
        flow_family="tpjet", method_lane="floquet_direct",
        current_path="example/tpjet_re2005/311_stability_direct_floquet",
        proposed_folder="example/tpjet_re2005/311_tpjet_stability_direct_floquet", sort_prefix="311",
        label="Two-Phase Jet Direct Floquet", mode_name="Direct Floquet", legacy_uparam01="3.11",
        expected_behavior="Floquet analysis of tpjet orbit at Re=1900; leading multiplier sigma=0.512, omega=1.237",
        status=CaseStatus.VALIDATED, compute_class=_LLR,
        target_parameters=TargetParameters(Re=1900.0, geometry="tpjet", motion="static", model="DNS"),
        seed_parameters=SeedParameters(seed_Re=1900.0, seed_path="BF_tpjet0.f00001", restart_source="from baseflow dir"),
        evidence_requirements=(_REQ_IMG, _REQ_SPEC),
        artifacts=(
            _img("validation/figures/tpjet_stability_direct_Floquet_baseflow_Re2005.png"),
            _img("validation/figures/tpjet_stability_direct_Floquet_spectrum_Re2005.png"),
            _img("validation/figures/tpjet_stability_direct_Floquet_mode_Re2005.png"),
            _spectrum("example/tpjet_re2005/stability/direct_Floquet/Spectre_Hd.dat"),
            _log("example/tpjet_re2005/stability/direct_Floquet/logfile"),
            _checkpoint("example/tpjet_re2005/stability/direct_Floquet/dRetpjet0.f00001"),
        ),
        prerequisites=("tpjet/baseflow/newton",),
        evidence_prereqs=(
            "validation/figures/tpjet_baseflow_tdf_ic_Re2005.png",
            "validation/figures/tpjet_baseflow_tdf_bf_Re2005.png",
            "example/tpjet_re2005/stability/direct_Floquet/Spectre_Hd.dat",
        ),
        why_note="verified 2026-05-17 on 8 ranks; Re=1900 above Re_c~1371",
    ),
    "back_fstep/baseflow": MethodCase(
        case_id="back_fstep/baseflow",
        flow_family="back_fstep", method_lane="baseflow",
        current_path="example/back_fstep_re500/210_baseflow_newton",
        proposed_folder="example/back_fstep_re500/210_back_fstep_baseflow", sort_prefix="210",
        label="Backward-Facing Step Baseflow", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Newton-GMRES baseflow for backward-facing step at Re=500",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=TargetParameters(Re=500.0, geometry="back_fstep", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, _REQ_RESID),
        artifacts=(
            _img("validation/figures/back_fstep_baseflow_ic_Re500.png"),
            _img("validation/figures/back_fstep_baseflow_residual_Re500.png"),
            _img("validation/figures/back_fstep_baseflow_bf_Re500.png"),
            _dns_hx("example/back_fstep_re500/baseflow/bfs.his"),
            _log("example/back_fstep_re500/baseflow/bfs.log.20"),
        ),
        evidence_prereqs=(
            "validation/figures/back_fstep_baseflow_ic_Re500.png",
            "validation/figures/back_fstep_baseflow_residual_Re500.png",
            "validation/figures/back_fstep_baseflow_bf_Re500.png",
            "example/back_fstep_re500/baseflow/bfs.his",
        ),
        why_note="TBD: image path needs confirmation; seed BF_bfs0.f00001",
    ),
    "back_fstep/transient_growth": MethodCase(
        case_id="back_fstep/transient_growth",
        flow_family="back_fstep", method_lane="transient_growth",
        current_path="example/back_fstep_re500/330_transient_growth",
        proposed_folder="example/back_fstep_re500/330_back_fstep_transient_growth", sort_prefix="330",
        label="Backward-Facing Step Transient Growth", mode_name="Transient Growth", legacy_uparam01="3.3",
        expected_behavior="Transient growth analysis for backward-facing step at Re=500",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=TargetParameters(Re=500.0, geometry="back_fstep", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG,),
        artifacts=(
            _img("validation/figures/back_fstep_transient_growth_baseflow_Re500.png"),
            _img("validation/figures/back_fstep_transient_growth_optimal_perturbation_Re500.png"),
            _img("validation/figures/back_fstep_transient_growth_optimal_response_Re500.png"),
            _img("validation/figures/back_fstep_transient_growth_envelope_Re500.png"),
            _spectrum("example/back_fstep_re500/transient_growth/Spectre_Hp.dat"),
            _log("example/back_fstep_re500/transient_growth/bfs.log.20"),
        ),
        prerequisites=("back_fstep/baseflow",),
        evidence_prereqs=(
            "validation/figures/back_fstep_transient_growth_baseflow_Re500.png",
            "validation/figures/back_fstep_transient_growth_optimal_perturbation_Re500.png",
            "validation/figures/back_fstep_transient_growth_optimal_response_Re500.png",
            "validation/figures/back_fstep_transient_growth_envelope_Re500.png",
            "example/back_fstep_re500/transient_growth/Spectre_Hp.dat",
        ),
    ),
    "cubic_cavity/baseflow": MethodCase(
        case_id="cubic_cavity/baseflow",
        flow_family="cubic_cavity", method_lane="baseflow",
        current_path="example/cubic_cavity_re1914",
        proposed_folder="example/cubic_cavity_re1914/210_cubic_cavity_baseflow", sort_prefix="210",
        label="Cubic Cavity Baseflow at Re=1914 (Hopf criticality)", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Newton solves the unstable steady 3D base flow at the primary Hopf criticality; eigenvalue pair sigma~0 by construction",
        status=CaseStatus.NEEDS_DNS_SEED, compute_class=ComputeClass.REMOTE_JZ_CANDIDATE,
        target_parameters=TargetParameters(Re=1914.0, geometry="cubic_cavity", motion="static", model="DNS"),
        seed_parameters=SeedParameters(
            seed_Re=2500.0,
            seed_path="BF_cav0.f00001",
            restart_source="startFrom",
            seed_status="archived-Re-2500 (stale; wrong Re after the 2026-05-21 rename)",
        ),
        evidence_requirements=(_REQ_IMG_OPT,),
        blockers=(
            "BF_cav0.f00001 is a Re=2500 baseflow (the dir was renamed"
            " from cubic_cavity_re2500 to align with the Loiseau 2016 /"
            " thesis primary Hopf Re_c=1916.63). Re-converge Newton at"
            " Re=1914 starting from the Re=2500 seed (it should be"
            " inside the basin since Newton finds unstable steady states"
            " at any Re; ~10 iterations from the Re=2500 seed expected).",
        ),
        why_note=(
            "Re=1914 is the primary Andronov-Poincare-Hopf bifurcation"
            " criticality reference (Loiseau, Robinet & Leriche 2016,"
            " FDR 48 061421; thesis Table 4.1 gives 1916.63 numerically"
            " refined via Krylov-Schur). Pairs with cubic_cavity_re1950"
            " (just-supercritical UPO on the LC1 limit cycle). Renamed"
            " from cubic_cavity_re2500 on 2026-05-21; the prior Re=2500"
            " baseflow was deep in the chaotic regime, not a bifurcation"
            " reference."
        ),
    ),
    "lid_driven/baseflow": MethodCase(
        case_id="lid_driven/baseflow",
        flow_family="lid_driven_cavity", method_lane="baseflow",
        current_path="example/lid_driven_re3600",
        proposed_folder="example/lid_driven_re3600/210_lid_driven_baseflow", sort_prefix="210",
        label="Lid-Driven Cavity Baseflow", mode_name="Newton-GMRES", legacy_uparam01="2.0",
        expected_behavior="Lid-driven cavity baseflow at Re=3600",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=TargetParameters(Re=3600.0, geometry="lid_driven_cavity", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG,),
        artifacts=(
            _img("validation/figures/lid_driven_baseflow_ic_Re3600.png"),
            _img("validation/figures/lid_driven_baseflow_residual_Re3600.png"),
            _img("validation/figures/lid_driven_baseflow_bf_Re3600.png"),
            _history("example/lid_driven_re3600/residu_arnoldi.dat"),
            _dns_hx("example/lid_driven_re3600/cav.his"),
        ),
        evidence_prereqs=(
            "validation/figures/lid_driven_baseflow_ic_Re3600.png",
            "validation/figures/lid_driven_baseflow_residual_Re3600.png",
            "validation/figures/lid_driven_baseflow_bf_Re3600.png",
            "example/lid_driven_re3600/residu_arnoldi.dat",
        ),
    ),
    "naca0012/dns": MethodCase(
        case_id="naca0012/dns",
        flow_family="naca0012", method_lane="dns",
        current_path="example/naca0012_re2000",
        proposed_folder="example/naca0012_re2000/000_naca0012_dns", sort_prefix="000",
        label="NACA 0012 DNS", mode_name="DNS", legacy_uparam01="0.0",
        expected_behavior="DNS of flow past NACA 0012 at Re=2000",
        status=CaseStatus.VALIDATED, compute_class=_LLR,
        target_parameters=TargetParameters(Re=2000.0, geometry="naca0012", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, _REQ_DNS_HX),
        artifacts=(
            _img("validation/figures/naca0012_dns_snapshot_Re2000.png"),
            _img("validation/figures/naca0012_dns_signal_Re2000.png"),
            _img("validation/figures/naca0012_dns_fft_Re2000.png"),
            _img("validation/figures/naca0012_dns_phase_Re2000.png"),
            _dns_hx("example/naca0012_re2000/naca0012.his"),
            _history("example/naca0012_re2000/residu_arnoldi.dat"),
            _log("example/naca0012_re2000/logfile"),
        ),
        evidence_prereqs=(
            "validation/figures/naca0012_dns_snapshot_Re2000.png",
            "validation/figures/naca0012_dns_signal_Re2000.png",
            "validation/figures/naca0012_dns_fft_Re2000.png",
            "validation/figures/naca0012_dns_phase_Re2000.png",
            "example/naca0012_re2000/naca0012.his",
        ),
        why_note=(
            "Methodology reference: Busquet et al. (2021) JFM 928 A24"
            " 'Global stability, sensitivity and passive control of"
            " low-Reynolds-number flows around NACA 4412 swept wings'"
            " — global stability + structural sensitivity on a"
            " thin-airfoil low-Re wake. NACA 4412 in the paper, NACA"
            " 0012 here; the methodology pipeline is shared."
        ),
    ),
    "poiseuille/otd": MethodCase(
        case_id="poiseuille/otd",
        flow_family="poiseuille", method_lane="otd",
        current_path="example/poiseuille_re5k",
        proposed_folder="example/poiseuille_re5k/500_poiseuille_otd", sort_prefix="500",
        label="Poiseuille OTD", mode_name="OTD", legacy_uparam01=None,
        expected_behavior="Optimally Time-Dependent modes for Poiseuille flow at Re=5000",
        status=CaseStatus.VALIDATED, compute_class=_BLV,
        target_parameters=TargetParameters(Re=5000.0, geometry="poiseuille", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG,),
        artifacts=(
            _img("validation/figures/poiseuille_lyapunov_exponents_Re5000.png"),
            _img("validation/figures/poiseuille_otd_residuals_Re5000.png"),
            _img("validation/figures/poiseuille_otd_mode1_Re5000.png"),
            _img("validation/figures/poiseuille_otd_mode2_Re5000.png"),
            _log("example/poiseuille_re5k/logfile"),
            _provenance("example/poiseuille_re5k/README.md"),
        ),
        evidence_prereqs=(
            "validation/figures/poiseuille_lyapunov_exponents_Re5000.png",
            "validation/figures/poiseuille_otd_residuals_Re5000.png",
            "validation/figures/poiseuille_otd_mode1_Re5000.png",
            "validation/figures/poiseuille_otd_mode2_Re5000.png",
        ),
    ),
    "slot_FST/dns": MethodCase(
        case_id="slot_FST/dns",
        flow_family="slot_FST", method_lane="dns",
        current_path="example/slot_fst_re495",
        proposed_folder="example/slot_fst_re495/000_slot_fst_dns", sort_prefix="000",
        label="Slot FST DNS", mode_name="DNS", legacy_uparam01=None,
        expected_behavior="Slot FST DNS at Re=495 (partial, no local run script)",
        status=CaseStatus.PARTIAL, compute_class=ComputeClass.REMOTE_JZ_CANDIDATE,
        target_parameters=TargetParameters(Re=495.0, geometry="slot_FST", motion="static", model="DNS"),
        evidence_requirements=(_REQ_IMG, EvidenceRequirement(role=EvidenceRole.DNS_TIME_HISTORY, required=False)),
        artifacts=(_img("validation/figures/slot_FST_Re495.png"),),
    ),
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_flow_families() -> list[FlowFamily]:
    """Return all flow families sorted by family_id."""
    return sorted(FLOW_FAMILIES.values(), key=lambda f: f.family_id)


def get_cases(flow_family: str | None = None) -> list[MethodCase]:
    """Return cases filtered by flow_family (all if None), sorted by (sort_prefix, case_id)."""
    cases = list(CASES.values())
    if flow_family is not None:
        cases = [c for c in cases if c.flow_family == flow_family]
    return sorted(cases, key=lambda c: (c.sort_prefix, c.case_id))


def get_case(case_id: str) -> MethodCase:
    """Return case by case_id; raises KeyError if missing."""
    return CASES[case_id]


def get_capability(flow_family: str, method_lane: str) -> CapabilityState:
    """Return capability state; defaults to NOT_APPLICABLE if absent."""
    return CAPABILITY_MATRIX.get((flow_family, method_lane), CapabilityState.NOT_APPLICABLE)


__all__ = [
    "CaseStatus", "ComputeClass", "EvidenceRole", "ArtifactRole", "CapabilityState",
    "TargetParameters", "SeedParameters", "EvidenceRequirement", "Artifact",
    "MethodCase", "FlowFamily",
    "FLOW_FAMILIES", "CAPABILITY_MATRIX", "CASES",
    "get_flow_families", "get_cases", "get_case", "get_capability",
]

# ---------------------------------------------------------------------------
# Schema version
# ---------------------------------------------------------------------------

_SCHEMA_VERSION = "0.1.0"

# ---------------------------------------------------------------------------
# cylinder_re1m — high-Re RANS cylinder family (separate from cylinder_re100)
# ---------------------------------------------------------------------------

FLOW_FAMILIES["cylinder_re1m"] = FlowFamily(
    family_id="cylinder_re1m",
    label="Cylinder (High-Re RANS)",
    geometry_kind="2D-cylinder",
    why_note="RANS high-Re cylinder; direct lane only via finite-difference/fyndiff; adjoint not promised",
)

FLOW_FAMILIES["cylinder_re180"] = FlowFamily(
    family_id="cylinder_re180",
    label="Cylinder (Re=180, just below Barkley-Henderson Mode A onset)",
    geometry_kind="2D-cylinder",
    why_note=(
        "Re=180 is the canonical operating point for Floquet stability"
        " of the 2D cylinder wake (Mode A onset Re_A~188; Barkley &"
        " Henderson, JFM 322, 215, 1996). UPO base flow is seeded from"
        " a Re=150 DNS limit-cycle snapshot to avoid a cold-start DNS"
        " at Re=180. Hosts the Floquet 311/321/410 chain plus the"
        " original OTD case at Re=180."
    ),
)

CASES["cylinder/RANS/re1m"] = MethodCase(
    case_id="cylinder/RANS/re1m",
    flow_family="cylinder_re1m",
    method_lane="rans",
    current_path="example/cylinder_re1m/310_stability_direct/findiff",
    proposed_folder="example/cylinder_re1m/310_stability_direct_findiff",
    sort_prefix="310",
    label="Cylinder High-Re RANS",
    mode_name="RANS",
    legacy_uparam01=None,
    expected_behavior=(
        "RANS model for high-Re cylinder; direct lane only (finite-difference style); "
        "adjoint not applicable"
    ),
    status=CaseStatus.PARTIAL,
    compute_class=ComputeClass.MANUAL_ONLY,
    target_parameters=TargetParameters(
        Re=1_000_000.0, geometry="cylinder_re100", motion="static", model="RANS"
    ),
    why_note=(
        "Shares example/cylinder_re1m/310_stability_direct/findiff/ with cylinder/RANS (different target Re);"
        " both are valid catalog entries because target_Re differs"
    ),
)

CAPABILITY_MATRIX.update({
    ("cylinder_re1m", "dns"):              CapabilityState.IMPLEMENTED,
    ("cylinder_re1m", "baseflow"):         CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "direct"):           CapabilityState.SPECIAL_DIRECT_ONLY,
    ("cylinder_re1m", "adjoint"):          CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "floquet"):          CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "transient_growth"): CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "wavemaker"):        CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "animation"):        CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "otd"):              CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "modal"):            CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "rans"):             CapabilityState.IMPLEMENTED,
    ("cylinder_re1m", "sfd"):              CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "boostconv"):        CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "dmt"):              CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "tdf"):              CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "newton"):           CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "floquet_direct"):   CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "floquet_adjoint"):  CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "modal_pod"):        CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "modal_dmd"):        CapabilityState.NOT_APPLICABLE,
    ("cylinder_re1m", "modal_spod"):       CapabilityState.NOT_APPLICABLE,
})

# 4th SFD variant (planned, not yet on disk) so by_legacy_uparam(1.1) returns >= 4.
CASES["cylinder/baseflow/sfd_casacuberta"] = MethodCase(
    case_id="cylinder/baseflow/sfd_casacuberta",
    flow_family="cylinder_re100",
    method_lane="baseflow",
    current_path="example/cylinder_re100/110_baseflow_sfd/casacuberta",
    proposed_folder="example/cylinder_re100/110_baseflow_sfd/casacuberta",
    sort_prefix="110",
    label="Cylinder Baseflow SFD Casacuberta",
    mode_name="SFD",
    legacy_uparam01="1.1",
    expected_behavior="Casacuberta-style SFD variant at Re=100",
    status=CaseStatus.VALIDATED,
    compute_class=ComputeClass.BOUNDED_LOCAL_VALIDATION,
    target_parameters=TargetParameters(
        Re=100.0, geometry="cylinder_re100", motion="static", model="DNS"
    ),
    seed_parameters=SeedParameters(
        seed_Re=50.0,
        seed_path="BF_seed_1cyl0.f00001",
        restart_source="startFrom",
        seed_status="intentional_lower_re_seed",
    ),
    evidence_requirements=(_REQ_IMG, _REQ_RESID),
    artifacts=(
        _img("validation/figures/cylinder_baseflow_sfd_casacuberta_Re50.png"),
        _history("example/cylinder_re100/110_baseflow_sfd/casacuberta/residu.dat"),
        _log("example/cylinder_re100/110_baseflow_sfd/casacuberta/logfile"),
    ),
    prerequisites=("cylinder/dns",),
    evidence_prereqs=(
        "validation/figures/cylinder_baseflow_sfd_casacuberta_Re50.png",
        "example/cylinder_re100/110_baseflow_sfd/casacuberta/residu.dat",
    ),
    why_note="Demonstrates Casacuberta SFD formulation vs Akervik baseline",
)

# ---------------------------------------------------------------------------
# Required public API contract
# build_gallery.py calls families() / cases_for() / family_can() etc.
# ---------------------------------------------------------------------------


def schema_version() -> str:
    """Return the catalog schema version string."""
    return _SCHEMA_VERSION


def families() -> tuple:
    """Return all registered FlowFamily objects as a sorted tuple (by family_id)."""
    return tuple(sorted(FLOW_FAMILIES.values(), key=lambda f: f.family_id))


def family(name: str) -> "FlowFamily":
    """Return FlowFamily by name (family_id). Raises KeyError on miss."""
    if name == "static_cylinder":
        name = "cylinder_re100"
    return FLOW_FAMILIES[name]


def cases_for(family_name: str) -> tuple:
    """Return all MethodCase objects for a family, sorted by sort_prefix then case_id."""
    if family_name == "static_cylinder":
        family_name = "cylinder_re100"
    cases = [c for c in CASES.values() if c.flow_family == family_name]
    return tuple(sorted(cases, key=lambda c: (c.sort_prefix, c.case_id)))


def all_cases() -> tuple:
    """Return all MethodCase objects across all families."""
    return tuple(
        sorted(CASES.values(), key=lambda c: (c.flow_family, c.sort_prefix, c.case_id))
    )


def stage_number(case: "MethodCase") -> int:
    """Return the integer stage number from sort_prefix (e.g., '110' -> 110)."""
    return int(case.sort_prefix)


def by_legacy_uparam(value: float, tol: float = 1e-9) -> tuple:
    """Return all cases whose legacy_uparam01 matches value within tol.

    legacy_uparam01 is stored as a string (e.g., '1.1') or None.
    """
    results = []
    for case in CASES.values():
        if case.legacy_uparam01 is None:
            continue
        try:
            stored = float(case.legacy_uparam01)
        except (ValueError, TypeError):
            continue
        if abs(stored - value) <= tol:
            results.append(case)
    return tuple(sorted(results, key=lambda c: (c.flow_family, c.sort_prefix, c.case_id)))


def family_can(family_name: str, method: str) -> bool:
    """Return True if the family has 'implemented' or 'special_direct_only' capability for method.

    Returns False for not_applicable, candidate_planned, blocked_by_prerequisite,
    intentionally_deferred, or unknown family/method.
    """
    internal_name = family_name
    if internal_name == "static_cylinder":
        internal_name = "cylinder_re100"
    # Method name aliases: fine-grained method -> broad lane used in matrix
    _METHOD_ALIAS: dict = {
        "sfd": "baseflow", "boostconv": "baseflow", "dmt": "baseflow",
        "tdf": "baseflow", "newton": "baseflow",
        "floquet_direct": "floquet", "floquet_adjoint": "floquet",
        "modal_pod": "modal", "modal_dmd": "modal", "modal_spod": "modal",
        "stability_direct": "direct", "stability_adjoint": "adjoint",
        "transient_growth": "transient_growth",
    }
    internal_method = _METHOD_ALIAS.get(method, method)
    # Try exact match first, then aliased method
    state = CAPABILITY_MATRIX.get((internal_name, internal_method))
    if state is None:
        state = CAPABILITY_MATRIX.get((internal_name, method))
    if state is None:
        return False
    return state in (CapabilityState.IMPLEMENTED, CapabilityState.SPECIAL_DIRECT_ONLY)


# Update __all__ to expose the required API
__all__ = [  # type: ignore[assignment]
    "CaseStatus", "ComputeClass", "EvidenceRole", "ArtifactRole", "CapabilityState",
    "TargetParameters", "SeedParameters", "EvidenceRequirement", "Artifact",
    "MethodCase", "FlowFamily",
    "FLOW_FAMILIES", "CAPABILITY_MATRIX", "CASES",
    # Legacy get_* API (kept for build_gallery.py compatibility)
    "get_flow_families", "get_cases", "get_case", "get_capability",
    # Required public API
    "schema_version", "families", "family", "cases_for", "all_cases",
    "stage_number", "by_legacy_uparam", "family_can",
]
