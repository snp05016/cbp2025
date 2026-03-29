import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ==========================================
# CONFIGURATION & CONSTANTS
# ==========================================

# Output Directory (can be overridden via CLI)
DEFAULT_OUTPUT_DIR = 'subsumption_analysis/figures_all_pairs'

# Plot Styling for High Readability
# Removed 'grid.axis' to fix KeyError
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 14,
    'axes.labelsize': 15,
    'axes.titlesize': 18,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'axes.grid': False,       # We will turn this on manually for Y-axis only
    'grid.alpha': 0.4,
    'grid.linestyle': '--'
})

# ==========================================
# DATA PROCESSING
# ==========================================

SKIP_CSV_NAMES = {
    'combined_predictor_results.csv',
    'pairwise_deltas.csv',
    'subsumption_summary.csv',
    'status.csv',
}


def _slugify(text: str) -> str:
    text = str(text).strip()
    text = re.sub(r'\s+', '_', text)
    text = re.sub(r'[^A-Za-z0-9._-]+', '_', text)
    text = re.sub(r'_+', '_', text)
    return text.strip('._-') or 'predictor'


def discover_predictor_csvs(results_dir: Path, recursive: bool) -> list[Path]:
    if recursive:
        candidates = sorted(p for p in results_dir.rglob('*.csv') if p.is_file())
    else:
        candidates = sorted(p for p in results_dir.glob('*.csv') if p.is_file())

    csvs: list[Path] = []
    for p in candidates:
        if p.name in SKIP_CSV_NAMES:
            continue
        csvs.append(p)
    return csvs


def load_predictor_table(csv_file: Path, predictor_name: str) -> pd.DataFrame | None:
    """Load a single predictor CSV and normalize to [trace_name, predictor_name, MPKI]."""
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"✗ Failed to read {csv_file}: {e}")
        return None

    if 'Run' not in df.columns or 'MPKI' not in df.columns:
        print(f"- Skipping {csv_file.name} (missing 'Run' and/or 'MPKI')")
        return None

    out = df[['Run', 'MPKI']].copy()
    out.rename(columns={'Run': 'trace_name'}, inplace=True)
    out['trace_name'] = out['trace_name'].astype(str).str.strip()
    out['MPKI'] = pd.to_numeric(out['MPKI'], errors='coerce')
    out['predictor_name'] = predictor_name
    out = out.dropna(subset=['trace_name', 'MPKI'])
    return out[['trace_name', 'predictor_name', 'MPKI']]


def compute_pairwise_deltas(long_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-trace MPKI deltas for every predictor pair.

    Delta definition: Δ(A→B) = MPKI_A - MPKI_B. Positive means B is better.

    Excludes traces where either MPKI is <= 0 (treated as invalid/degenerate).
    """
    predictors = sorted(long_df['predictor_name'].unique())
    deltas: list[pd.DataFrame] = []

    for i, pred_a in enumerate(predictors):
        data_a = long_df[long_df['predictor_name'] == pred_a][['trace_name', 'MPKI']].rename(
            columns={'MPKI': 'MPKI_A'}
        )
        for pred_b in predictors[i + 1:]:
            data_b = long_df[long_df['predictor_name'] == pred_b][['trace_name', 'MPKI']].rename(
                columns={'MPKI': 'MPKI_B'}
            )
            merged = pd.merge(data_a, data_b, on='trace_name', how='inner')
            valid = merged[(merged['MPKI_A'] > 0) & (merged['MPKI_B'] > 0)].copy()
            if valid.empty:
                continue
            valid['pred_a'] = pred_a
            valid['pred_b'] = pred_b
            valid['delta'] = valid['MPKI_A'] - valid['MPKI_B']
            deltas.append(valid[['pred_a', 'pred_b', 'trace_name', 'MPKI_A', 'MPKI_B', 'delta']])

    if not deltas:
        return pd.DataFrame(columns=['pred_a', 'pred_b', 'trace_name', 'MPKI_A', 'MPKI_B', 'delta'])

    return pd.concat(deltas, ignore_index=True)

# ==========================================
# PLOTTING FUNCTION (Violin Plot)
# ==========================================

def plot_pair_violin(deltas: np.ndarray, pred_a: str, pred_b: str, output_dir: Path) -> Path:
    """Generate a single violin plot for one predictor pair."""
    output_dir.mkdir(parents=True, exist_ok=True)

    deltas = deltas[np.isfinite(deltas)]
    median = float(np.median(deltas)) if len(deltas) else float('nan')

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.grid(True, axis='y', linestyle='--', alpha=0.4)

    parts = ax.violinplot([deltas], showmeans=False, showmedians=False, vert=True)
    for pc in parts['bodies']:
        pc.set_facecolor('#1f77b4')
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)

    for partname in ('cbars', 'cmins', 'cmaxes'):
        if partname in parts:
            v = parts[partname]
            v.set_edgecolor('black')
            v.set_linewidth(1)

    ax.scatter([1], [median], marker='o', color='white', s=60, zorder=3,
               edgecolors='black', linewidth=1.5, label='Median Δ')

    ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.8)

    ax.set_xticks([1])
    ax.set_xticklabels([f"{pred_a}\n↓\n{pred_b}"])
    ax.set_ylabel(r'$\Delta$ MPKI (Higher is better for 2nd predictor)')
    ax.set_title(f"Delta MPKI Distribution (Violin)\n{pred_a} vs {pred_b}")
    ax.legend(loc='upper right')
    plt.tight_layout()

    safe_a = _slugify(pred_a)
    safe_b = _slugify(pred_b)
    filepath = output_dir / f"delta_violin_{safe_a}_vs_{safe_b}.png"
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return filepath

# ==========================================
# MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate pairwise delta subsumption violin plots for all predictors'
    )
    parser.add_argument('--results-dir', type=str, default='results',
                        help='Directory containing predictor CSVs (default: results)')
    parser.add_argument('--recursive', action='store_true',
                        help='Recursively search for CSVs under --results-dir')
    parser.add_argument('--output-dir', type=str, default=DEFAULT_OUTPUT_DIR,
                        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Discovering CSVs under: {results_dir}")
    csv_files = discover_predictor_csvs(results_dir, recursive=args.recursive)
    if not csv_files:
        print(f"Error: no CSV files found under {results_dir}")
        return

    print(f"Found {len(csv_files)} CSV(s). Loading MPKI tables...")
    tables: list[pd.DataFrame] = []
    for csv_file in csv_files:
        predictor_name = csv_file.stem
        t = load_predictor_table(csv_file, predictor_name=predictor_name)
        if t is not None and not t.empty:
            tables.append(t)

    if not tables:
        print("Error: no usable predictor CSVs (need columns 'Run' and 'MPKI').")
        return

    long_df = pd.concat(tables, ignore_index=True)
    predictors = sorted(long_df['predictor_name'].unique())
    print(f"Loaded {len(long_df)} rows across {len(predictors)} predictor(s).")

    print("Computing pairwise deltas...")
    delta_df = compute_pairwise_deltas(long_df)
    if delta_df.empty:
        print("No valid pairwise deltas could be computed (check MPKI values / shared traces).")
        return

    delta_csv = output_dir / 'pairwise_deltas_all_predictors.csv'
    delta_df.to_csv(delta_csv, index=False)
    print(f"✓ Wrote delta table: {delta_csv}")

    print("Generating per-pair violin plots...")
    pairs_out_dir = output_dir / 'pairwise_violin_plots'
    pairs_out_dir.mkdir(parents=True, exist_ok=True)

    plot_count = 0
    for (pred_a, pred_b), grp in delta_df.groupby(['pred_a', 'pred_b'], sort=True):
        deltas = grp['delta'].to_numpy(dtype=float)
        if len(deltas) == 0:
            continue
        out = plot_pair_violin(deltas, pred_a=pred_a, pred_b=pred_b, output_dir=pairs_out_dir)
        plot_count += 1
        print(f"✓ {out}")

    print(f"Done. Generated {plot_count} plot(s) in: {pairs_out_dir}")

if __name__ == "__main__":
    main()