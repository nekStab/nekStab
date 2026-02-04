#!/usr/bin/env python3
"""
Compare benchmark results from GCC, ifort, and ifx compiler builds.
Plots velocity signals and computes error estimation.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

RESULTS_DIR = Path(__file__).parent / "benchmark_results"

def load_his_file(filepath):
    """Load Nek5000 .his file (history points)."""
    data = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
        # Skip header lines (first 2)
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    data.append([float(x) for x in parts[:4]])
                except ValueError:
                    continue
    return np.array(data) if data else None

def load_timing(filepath):
    """Load timing results from CSV (no header)."""
    timing = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 3:
                compiler = parts[0]
                nek_time = float(parts[1]) if parts[1] else 0
                wall_time = float(parts[2]) if parts[2] else 0
                timing[compiler] = {'nek_time': nek_time, 'wall_time': wall_time}
    return timing

def compute_errors(data1, data2):
    """Compute various error metrics between two datasets."""
    t1, t2 = data1[:, 0], data2[:, 0]
    t_min = max(t1.min(), t2.min())
    t_max = min(t1.max(), t2.max())
    dt = min(np.diff(t1).mean(), np.diff(t2).mean())
    t_common = np.arange(t_min, t_max, dt)

    errors = {}
    for col, name in enumerate(['vx', 'vy', 'vz'], start=1):
        v1 = np.interp(t_common, t1, data1[:, col])
        v2 = np.interp(t_common, t2, data2[:, col])
        diff = v1 - v2
        errors[name] = {
            'max_abs_error': np.max(np.abs(diff)),
            'mean_abs_error': np.mean(np.abs(diff)),
            'rms_error': np.sqrt(np.mean(diff**2)),
            'max_rel_error': np.max(np.abs(diff)) / (np.max(np.abs(v1)) + 1e-15),
        }
    return errors

def main():
    print("=" * 65)
    print("nekStab Benchmark: GCC vs ifort (Classic) vs ifx (LLVM)")
    print("=" * 65)

    timing_csv = RESULTS_DIR / "timing.csv"

    if not RESULTS_DIR.exists():
        print(f"Error: Results directory '{RESULTS_DIR}' not found.")
        print("Run benchmark.sh first to generate results.")
        sys.exit(1)

    # Load timing data
    timing = {}
    if timing_csv.exists():
        timing = load_timing(timing_csv)
        print("\n Performance Comparison:")
        print("-" * 50)
        print(f"  {'Compiler':<10} {'Nek5000 time':>14} {'Wall time':>12} {'Speedup':>10}")
        print("-" * 50)

        ref_time = timing.get('gcc', {}).get('nek_time', 1)
        for compiler in ['gcc', 'ifort', 'ifx']:
            if compiler in timing:
                t = timing[compiler]
                speedup = ref_time / t['nek_time'] if t['nek_time'] > 0 else 0
                print(f"  {compiler.upper():<10} {t['nek_time']:>12.3f}s {t['wall_time']:>10.3f}s {speedup:>9.2f}x")
        print("-" * 50)

    # Load velocity data
    datasets = {}
    compilers = ['gcc', 'ifort', 'ifx']

    for compiler in compilers:
        his_path = RESULTS_DIR / f"{compiler}.his"
        if his_path.exists():
            data = load_his_file(his_path)
            if data is not None:
                datasets[compiler] = data
                print(f"  Loaded {compiler}: {len(data)} timesteps")

    if len(datasets) < 2:
        print(f"\nWarning: Only found {len(datasets)} data files: {list(datasets.keys())}")

    # Compute pairwise errors
    pairs = [('gcc', 'ifort'), ('gcc', 'ifx'), ('ifort', 'ifx')]
    all_errors = {}

    print("\n Accuracy Comparison:")
    print("-" * 50)

    for c1, c2 in pairs:
        if c1 in datasets and c2 in datasets:
            errors = compute_errors(datasets[c1], datasets[c2])
            all_errors[(c1, c2)] = errors
            max_rel = max(e['max_rel_error'] for e in errors.values())
            print(f"  {c1.upper()} vs {c2.upper()}: max relative error = {max_rel:.2e}")

    # Create plots
    fig = plt.figure(figsize=(15, 10))

    colors = {'gcc': '#2ca02c', 'ifort': '#1f77b4', 'ifx': '#d62728'}
    labels = {'gcc': 'GCC 13.3', 'ifort': 'ifort 2021 (Classic)', 'ifx': 'ifx 2025 (LLVM)'}

    # Velocity plots (top row)
    for idx, (col, name) in enumerate([(1, 'vx'), (2, 'vy')], start=1):
        ax = fig.add_subplot(2, 3, idx)
        for compiler in ['gcc', 'ifort', 'ifx']:
            if compiler in datasets:
                data = datasets[compiler]
                ax.plot(data[:, 0], data[:, col],
                       color=colors[compiler],
                       label=labels[compiler],
                       alpha=0.8,
                       linewidth=1.5 if compiler == 'gcc' else 1.0)
        ax.set_xlabel('Time')
        ax.set_ylabel(name)
        ax.set_title(f'Velocity {name}')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)

    # Error plots (middle)
    ax = fig.add_subplot(2, 3, 3)
    if len(datasets) >= 2:
        ref_compiler = 'gcc' if 'gcc' in datasets else list(datasets.keys())[0]
        ref_data = datasets[ref_compiler]
        t_ref = ref_data[:, 0]

        for compiler in datasets:
            if compiler != ref_compiler:
                data = datasets[compiler]
                t_common = np.linspace(max(t_ref.min(), data[:, 0].min()),
                                      min(t_ref.max(), data[:, 0].max()), 500)
                v_ref = np.interp(t_common, t_ref, ref_data[:, 1])  # vx only
                v_comp = np.interp(t_common, data[:, 0], data[:, 1])
                diff = np.abs(v_ref - v_comp)
                ax.plot(t_common, diff, label=f'{labels[compiler]}',
                       color=colors[compiler], alpha=0.8)

        ax.set_xlabel('Time')
        ax.set_ylabel('|vx_gcc - vx_other|')
        ax.set_title('Velocity difference vs GCC')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')

    # Timing bar chart
    ax = fig.add_subplot(2, 3, 4)
    if timing:
        compilers_t = [c for c in ['gcc', 'ifort', 'ifx'] if c in timing]
        nek_times = [timing[c]['nek_time'] for c in compilers_t]
        wall_times = [timing[c]['wall_time'] for c in compilers_t]

        x = np.arange(len(compilers_t))
        width = 0.35

        bars1 = ax.bar(x - width/2, nek_times, width, label='Nek5000 time',
                      color=[colors[c] for c in compilers_t], alpha=0.7)
        bars2 = ax.bar(x + width/2, wall_times, width, label='Wall time',
                      color=[colors[c] for c in compilers_t], alpha=0.4)

        ax.set_ylabel('Time (s)')
        ax.set_title('Performance Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels([c.upper() for c in compilers_t])
        ax.legend()
        for bar, t in zip(bars1, nek_times):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                   f'{t:.3f}', ha='center', va='bottom', fontsize=9)

    # Speedup chart
    ax = fig.add_subplot(2, 3, 5)
    if timing and 'gcc' in timing:
        ref_time = timing['gcc']['nek_time']
        compilers_t = [c for c in ['gcc', 'ifort', 'ifx'] if c in timing]
        speedups = [ref_time / timing[c]['nek_time'] for c in compilers_t]

        bars = ax.bar(range(len(compilers_t)), speedups,
                     color=[colors[c] for c in compilers_t])
        ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
        ax.set_ylabel('Speedup vs GCC')
        ax.set_title('Relative Performance')
        ax.set_xticks(range(len(compilers_t)))
        ax.set_xticklabels([c.upper() for c in compilers_t])
        for bar, s in zip(bars, speedups):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                   f'{s:.2f}x', ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Summary text
    ax = fig.add_subplot(2, 3, 6)
    ax.axis('off')

    summary_lines = ["BENCHMARK SUMMARY", "-" * 30]

    if timing:
        fastest = min(timing.items(), key=lambda x: x[1]['nek_time'])
        summary_lines.append(f"Fastest: {fastest[0].upper()} ({fastest[1]['nek_time']:.3f}s)")

        if 'gcc' in timing:
            ref = timing['gcc']['nek_time']
            for c in ['ifort', 'ifx']:
                if c in timing:
                    speedup = ref / timing[c]['nek_time']
                    summary_lines.append(f"{c.upper()} vs GCC: {speedup:.2f}x")

    if all_errors:
        summary_lines.append("")
        summary_lines.append("Accuracy:")
        for (c1, c2), errs in all_errors.items():
            max_rel = max(e['max_rel_error'] for e in errs.values())
            summary_lines.append(f"  {c1}/{c2}: {max_rel:.2e}")

    ax.text(0.1, 0.9, '\n'.join(summary_lines), transform=ax.transAxes,
           fontsize=11, verticalalignment='top', fontfamily='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('nekStab Compiler Benchmark: GCC vs ifort vs ifx',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_file = RESULTS_DIR / "benchmark_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n Plot saved to: {output_file}")

    # Final summary
    print("\n" + "=" * 65)
    print("FINAL RESULTS")
    print("=" * 65)

    if timing:
        sorted_t = sorted(timing.items(), key=lambda x: x[1]['nek_time'])
        print("\nRanking (fastest to slowest):")
        for i, (compiler, t) in enumerate(sorted_t, 1):
            speedup = timing['gcc']['nek_time'] / t['nek_time'] if 'gcc' in timing else 1
            print(f"  {i}. {compiler.upper():<8} {t['nek_time']:.3f}s (Nek), {t['wall_time']:.3f}s (wall) - {speedup:.2f}x vs GCC")

    print("=" * 65)

if __name__ == "__main__":
    main()
