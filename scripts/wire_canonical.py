#!/usr/bin/env python3
"""Wire 9 done cases into catalog.py + collect_split_plots.py with canonical panel names.

Per-case maps:
  case_id -> [(local_panel_filename, canonical_root_with_re_tag), ...]

Updates:
  - catalog.py: replace each affected case's artifacts (single _img) with N _img() entries
                + updates evidence_prereqs to match
  - collect_split_plots.py: extends MAP with (local, canonical) tuples
"""
from pathlib import Path
import re

REPO = Path(__file__).resolve().parents[1]
CATALOG = REPO / 'validation' / 'catalog.py'
COLLECT = REPO / 'validation' / 'collect_split_plots.py'

# (case_id, case_dir relative, list of (local_filename, canonical_basename_with_Re_tag))
WIRE = [
    ('cylinder/dns', 'example/cylinder_re100/000_dns', [
        ('plot_snapshot.png', 'cylinder_dns_snapshot_Re50.png'),
        ('plot_signal.png',   'cylinder_dns_signal_Re50.png'),
        ('plot_fft.png',      'cylinder_dns_fft_Re50.png'),
        ('plot_phase.png',    'cylinder_dns_phase_Re50.png'),
    ]),
    ('cylinder/baseflow/sfd', 'example/cylinder_re100/110_baseflow_sfd/akervik', [
        ('plot_ic.png',        'cylinder_baseflow_sfd_ic_Re50.png'),
        ('plot_residual.png',  'cylinder_baseflow_sfd_residual_Re50.png'),
        ('plot_bf.png',        'cylinder_baseflow_sfd_bf_Re50.png'),
    ]),
    ('cylinder/baseflow/sfd_casacuberta', 'example/cylinder_re100/110_baseflow_sfd/casacuberta', [
        ('plot_ic.png',        'cylinder_baseflow_sfd_casacuberta_ic_Re50.png'),
        ('plot_residual.png',  'cylinder_baseflow_sfd_casacuberta_residual_Re50.png'),
        ('plot_bf.png',        'cylinder_baseflow_sfd_casacuberta_bf_Re50.png'),
    ]),
    ('cylinder/baseflow/sfd_akervik_dyn', 'example/cylinder_re100/110_baseflow_sfd/akervik_dyn', [
        ('plot_ic.png',                  'cylinder_baseflow_sfd_akervik_dyn_ic_Re50.png'),
        ('plot_residual.png',            'cylinder_baseflow_sfd_akervik_dyn_residual_Re50.png'),
        ('plot_dynamic_tolerance.png',   'cylinder_baseflow_sfd_akervik_dyn_dynamic_tolerance_Re50.png'),
        ('plot_bf.png',                  'cylinder_baseflow_sfd_akervik_dyn_bf_Re50.png'),
    ]),
    ('cylinder/baseflow/sfd_casacuberta_dyn', 'example/cylinder_re100/110_baseflow_sfd/casacuberta_dyn', [
        ('plot_ic.png',                  'cylinder_baseflow_sfd_casacuberta_dyn_ic_Re50.png'),
        ('plot_residual.png',            'cylinder_baseflow_sfd_casacuberta_dyn_residual_Re50.png'),
        ('plot_dynamic_tolerance.png',   'cylinder_baseflow_sfd_casacuberta_dyn_dynamic_tolerance_Re50.png'),
        ('plot_bf.png',                  'cylinder_baseflow_sfd_casacuberta_dyn_bf_Re50.png'),
    ]),
    ('cylinder/baseflow/boostconv', 'example/cylinder_re100/120_baseflow_boostconv', [
        ('plot_ic.png',        'cylinder_baseflow_boostconv_ic_Re50.png'),
        ('plot_residual.png',  'cylinder_baseflow_boostconv_residual_Re50.png'),
        ('plot_bf.png',        'cylinder_baseflow_boostconv_bf_Re50.png'),
    ]),
    ('cylinder/baseflow/newton', 'example/cylinder_re100/210_baseflow_newton/fp', [
        ('plot_ic.png',        'cylinder_baseflow_newton_ic_Re50.png'),
        ('plot_residual.png',  'cylinder_baseflow_newton_residual_Re50.png'),
        ('plot_bf.png',        'cylinder_baseflow_newton_bf_Re50.png'),
    ]),
    ('cylinder/baseflow/newton_dyn', 'example/cylinder_re100/210_baseflow_newton/dyn', [
        ('plot_ic.png',        'cylinder_baseflow_newton_dyn_ic_Re50.png'),
        ('plot_residual.png',  'cylinder_baseflow_newton_dyn_residual_Re50.png'),
        ('plot_bf.png',        'cylinder_baseflow_newton_dyn_bf_Re50.png'),
    ]),
    ('cylinder/baseflow/newton_upo', 'example/cylinder_re180/210_baseflow_newton/upo', [
        ('plot_ic.png',        'cylinder_baseflow_newton_upo_ic_Re180.png'),
        ('plot_residual.png',  'cylinder_baseflow_newton_upo_residual_Re180.png'),
        ('plot_bf.png',        'cylinder_baseflow_newton_upo_bf_Re180.png'),
    ]),
]

def update_catalog():
    """Replace the artifacts=(...) block for each affected case in catalog.py."""
    text = CATALOG.read_text()
    replacements = 0
    for case_id, _, panels in WIRE:
        # Build new _img list and evidence_prereqs list
        new_imgs = ',\n            '.join(
            f'_img("validation/figures/{can}")' for _, can in panels
        )
        new_prereq_paths = ',\n            '.join(
            f'"validation/figures/{can}"' for _, can in panels
        )
        # Find the MethodCase block for this case_id
        pat = rf'("{re.escape(case_id)}":\s*MethodCase\([^)]*?artifacts=\()(.*?)(\),\s*\n)'
        # use re.DOTALL to span lines; find the immediate _img(...) entries in artifacts=(
        m = re.search(
            rf'("{re.escape(case_id)}":\s*MethodCase\(.*?artifacts=\(\s*)(.*?)(\s*\),\s*\n)',
            text, re.DOTALL
        )
        if not m:
            print(f'  SKIP (case not found): {case_id}')
            continue
        head, body, tail = m.group(1), m.group(2), m.group(3)
        # Keep any non-_img() artifacts (e.g. _log, _history, _spectrum)
        keep_lines = []
        for line in body.split('\n'):
            stripped = line.strip().rstrip(',')
            if stripped and not stripped.startswith('_img('):
                keep_lines.append(line)
        # Build new body
        new_body_lines = [f'            _img("validation/figures/{can}"),' for _, can in panels]
        if keep_lines:
            new_body_lines.extend(line.rstrip().rstrip(',') + ',' for line in keep_lines)
        new_body = '\n'.join(new_body_lines)
        new_full = head + '\n' + new_body + tail
        text = text.replace(m.group(0), new_full, 1)
        replacements += 1
        print(f'  catalog updated: {case_id} -> {len(panels)} _img() entries')

        # Update evidence_prereqs for this case too
        pat_prereq = rf'("{re.escape(case_id)}":\s*MethodCase\(.*?evidence_prereqs=\(\s*)(.*?)(\s*\),\s*\n)'
        mp = re.search(pat_prereq, text, re.DOTALL)
        if mp:
            head_p, body_p, tail_p = mp.group(1), mp.group(2), mp.group(3)
            keep_p_lines = []
            for line in body_p.split('\n'):
                stripped = line.strip().rstrip(',')
                if stripped and 'validation/figures/' not in stripped:
                    keep_p_lines.append(line)
            new_p_lines = [f'            "validation/figures/{can}",' for _, can in panels]
            if keep_p_lines:
                new_p_lines.extend(line.rstrip().rstrip(',') + ',' for line in keep_p_lines)
            new_p_body = '\n'.join(new_p_lines)
            new_p_full = head_p + '\n' + new_p_body + tail_p
            text = text.replace(mp.group(0), new_p_full, 1)
            print(f'    evidence_prereqs updated too')

    CATALOG.write_text(text)
    print(f'\nTotal catalog edits: {replacements}')


def update_collect_map():
    """Insert per-case entries into collect_split_plots.py MAP."""
    text = COLLECT.read_text()
    # Find the closing `]` of the MAP list
    # Insert new entries right before `]`
    new_entries = []
    for _, case_dir, panels in WIRE:
        plist = ',\n        '.join(
            f'({local!r}, {can!r})' for local, can in panels
        )
        new_entries.append(f'    ({case_dir!r}, [\n        {plist},\n    ]),')
    new_block = '\n'.join(new_entries) + '\n'

    # Insert just before the closing `]\n` of MAP = [
    m = re.search(r'(MAP\s*=\s*\[\s*\n)(.*?)(\n\])', text, re.DOTALL)
    if not m:
        print('  ERROR: could not find MAP block in collect_split_plots.py')
        return
    head, body, tail = m.group(1), m.group(2), m.group(3)
    # Strip duplicate case_dir entries we already have
    existing_dirs = set()
    for case_dir, _, _ in [(c, None, None) for _, c, _ in WIRE]:
        existing_dirs.add(case_dir)
    body_lines = body.split('\n')
    kept_body_lines = []
    skip_until_close = False
    for line in body_lines:
        if any(f"({d!r}, [" in line for d in existing_dirs):
            skip_until_close = True
            continue
        if skip_until_close:
            if line.strip().endswith(']),'):
                skip_until_close = False
            continue
        kept_body_lines.append(line)
    new_body = '\n'.join(kept_body_lines).rstrip('\n') + '\n' + new_block.rstrip('\n')
    new_full = head + new_body + tail
    text = text.replace(m.group(0), new_full)
    COLLECT.write_text(text)
    print(f'collect MAP updated with {len(WIRE)} cases')


if __name__ == '__main__':
    update_catalog()
    update_collect_map()
