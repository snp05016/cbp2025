"""
Delta Curve Visualization for Branch Predictor Subsumption Analysis

Generates sorted delta curves (CDF-style line plots) showing MPKI improvements
across all traces for each predictor pair. Curves above zero indicate subsumption;
curves crossing zero indicate complementary behavior.

Usage:
    python delta_curve_visualization.py <combined_csv> [--output-dir <dir>]

Output:
    - PNG plots: one per predictor pair (sorted delta curves)
    - CSV files: delta values per pair
    - Summary statistics per curve
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path
from typing import Tuple, Dict


class DeltaCurveGenerator:
    """Generates sorted delta curves for predictor pairs."""
    
    def __init__(self, df: pd.DataFrame, output_dir: str = 'delta_curves'):
        """
        Initialize with predictor results.
        
        Args:
            df: DataFrame with columns [trace_name, predictor_name, MPKI]
            output_dir: Directory for output plots and CSVs
        """
        self.df = df.copy()
        self.predictors = sorted(df['predictor_name'].unique())
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _get_valid_deltas(self, pred_a: str, pred_b: str) -> Tuple[np.ndarray, int]:
        """
        Extract valid delta values for predictors A and B.
        
        Excludes traces where either predictor has zero MPKI.
        
        Args:
            pred_a: Name of predictor A
            pred_b: Name of predictor B
            
        Returns:
            Tuple of (sorted_deltas array, skipped_count)
        """
        # Get data for each predictor
        data_a = self.df[self.df['predictor_name'] == pred_a][['trace_name', 'MPKI']]
        data_b = self.df[self.df['predictor_name'] == pred_b][['trace_name', 'MPKI']]
        
        # Merge on trace name
        merged = pd.merge(data_a, data_b, on='trace_name', suffixes=('_A', '_B'))
        
        total_traces = len(merged)
        
        # Filter out traces where either predictor has zero MPKI
        valid = merged[(merged['MPKI_A'] > 0) & (merged['MPKI_B'] > 0)].copy()
        
        skipped = total_traces - len(valid)
        
        # Compute delta: Δ_A→B = MPKI_A - MPKI_B (positive = B better)
        deltas = (valid['MPKI_A'] - valid['MPKI_B']).values
        
        # Sort deltas for curve plotting
        sorted_deltas = np.sort(deltas)
        
        return sorted_deltas, skipped, valid[['trace_name', 'MPKI_A', 'MPKI_B']].assign(
            delta=deltas
        ).sort_values('delta')
    
    def _compute_statistics(self, deltas: np.ndarray) -> Dict:
        """
        Compute statistics for delta curve.
        
        Args:
            deltas: Array of delta values
            
        Returns:
            Dictionary with statistics
        """
        return {
            'mean': np.mean(deltas),
            'median': np.median(deltas),
            'std': np.std(deltas),
            'min': np.min(deltas),
            'max': np.max(deltas),
            'pct_positive': 100 * np.sum(deltas > 0) / len(deltas),
            'pct_negative': 100 * np.sum(deltas < 0) / len(deltas),
            'pct_zero': 100 * np.sum(deltas == 0) / len(deltas),
        }
    
    def _classify_curve(self, deltas: np.ndarray, stats: Dict) -> str:
        """
        Classify the subsumption relationship based on delta curve.
        
        Args:
            deltas: Array of delta values
            stats: Statistics dictionary
            
        Returns:
            Classification string
        """
        pct_pos = stats['pct_positive']
        mean_delta = stats['mean']
        min_delta = stats['min']
        
        if pct_pos >= 90 and min_delta > -0.5:
            return "Strong Subsumption"
        elif pct_pos >= 75 and min_delta > -1.0:
            return "Moderate Subsumption"
        elif 30 <= pct_pos <= 70:
            return "Complementary"
        elif pct_pos < 30 and mean_delta < 0:
            return "A Dominates B"
        else:
            return "Mixed/Unclear"
    
    def generate_curve(self, pred_a: str, pred_b: str) -> None:
        """
        Generate sorted delta curve plot for predictor pair.
        
        Args:
            pred_a: Name of predictor A (baseline)
            pred_b: Name of predictor B (comparison)
        """
        deltas, skipped, delta_df = self._get_valid_deltas(pred_a, pred_b)
        
        if len(deltas) == 0:
            print(f"⚠️  No valid traces for {pred_a} vs {pred_b}")
            return
        
        # Compute statistics
        stats = self._compute_statistics(deltas)
        classification = self._classify_curve(deltas, stats)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Plot sorted delta curve
        trace_indices = np.arange(len(deltas))
        ax.plot(trace_indices, deltas, linewidth=2.5, color='#2E86AB', label='Δ MPKI (A→B)')
        
        # Add zero line
        ax.axhline(y=0, color='red', linestyle='--', linewidth=2, label='Δ = 0 (No Change)', alpha=0.7)
        
        # Shade regions
        ax.fill_between(trace_indices, 0, deltas, where=(deltas > 0), 
                       alpha=0.2, color='green', label='B Better (Δ > 0)')
        ax.fill_between(trace_indices, 0, deltas, where=(deltas < 0), 
                       alpha=0.2, color='red', label='A Better (Δ < 0)')
        
        # Labels and title
        ax.set_xlabel('Traces Sorted by Δ MPKI', fontsize=12, fontweight='bold')
        ax.set_ylabel('Δ MPKI (MPKI_A - MPKI_B)', fontsize=12, fontweight='bold')
        ax.set_title(
            f'Delta Curve: {pred_a} vs {pred_b}\n'
            f'{classification} | Mean Δ = {stats["mean"]:.4f} | '
            f'{stats["pct_positive"]:.1f}% Better',
            fontsize=14, fontweight='bold'
        )
        
        # Grid and legend
        ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.8)
        ax.legend(loc='upper left', fontsize=11, framealpha=0.95)
        
        # Add statistics box
        stats_text = (
            f'Valid Traces: {len(deltas)} (Skipped: {skipped})\n'
            f'Mean: {stats["mean"]:+.4f} | Median: {stats["median"]:+.4f}\n'
            f'Std Dev: {stats["std"]:.4f}\n'
            f'Range: [{stats["min"]:+.4f}, {stats["max"]:+.4f}]\n'
            f'Positive: {stats["pct_positive"]:.1f}% | Negative: {stats["pct_negative"]:.1f}%'
        )
        ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
               family='monospace')
        
        plt.tight_layout()
        
        # Save plot
        plot_file = self.output_dir / f'delta_curve_{pred_a}_vs_{pred_b}.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {plot_file}")
        plt.close()
        
        # Save delta data
        csv_file = self.output_dir / f'deltas_{pred_a}_vs_{pred_b}.csv'
        delta_df.to_csv(csv_file, index=False)
        
        # Print summary
        print(f"\n{'='*80}")
        print(f"Delta Curve: {pred_a} vs {pred_b}")
        print(f"{'='*80}")
        print(f"Classification: {classification}")
        print(f"Valid Traces: {len(deltas)} (Skipped: {skipped} zero-MPKI)")
        print(f"\nStatistics:")
        print(f"  Mean Δ:     {stats['mean']:+.4f} MPKI")
        print(f"  Median Δ:   {stats['median']:+.4f} MPKI")
        print(f"  Std Dev:    {stats['std']:.4f} MPKI")
        print(f"  Range:      [{stats['min']:+.4f}, {stats['max']:+.4f}]")
        print(f"\nDistribution:")
        print(f"  {stats['pct_positive']:.1f}% traces: {pred_b} better (Δ > 0)")
        print(f"  {stats['pct_negative']:.1f}% traces: {pred_a} better (Δ < 0)")
        print(f"  {stats['pct_zero']:.1f}% traces: identical performance")
        print(f"\nInterpretation:")
        
        if classification == "Strong Subsumption":
            print(f"  {pred_b} strictly subsumes {pred_a}: improvements across nearly all traces")
            print(f"  with minimal regressions. {pred_b} is strictly superior.")
        elif classification == "Moderate Subsumption":
            print(f"  {pred_b} mostly subsumes {pred_a}: improvements on majority of traces")
            print(f"  with acceptable worst-case behavior. {pred_b} is generally superior.")
        elif classification == "Complementary":
            print(f"  {pred_a} and {pred_b} are complementary: each excels on different workloads.")
            print(f"  A hybrid or ensemble approach may be beneficial.")
        elif classification == "A Dominates B":
            print(f"  {pred_a} dominates {pred_b}: consistently better across all traces.")
        else:
            print(f"  Mixed behavior: no clear dominance pattern detected.")
    
    def generate_all_curves(self) -> None:
        """Generate delta curves for all predictor pairs."""
        pairs = []
        for i, pred_a in enumerate(self.predictors):
            for pred_b in self.predictors[i+1:]:
                pairs.append((pred_a, pred_b))
        
        print(f"\nGenerating {len(pairs)} delta curve(s)...\n")
        
        for pred_a, pred_b in pairs:
            self.generate_curve(pred_a, pred_b)
        
        print(f"\n{'='*80}")
        print(f"All delta curves saved to: {self.output_dir.absolute()}")
        print(f"{'='*80}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Generate sorted delta curves for predictor subsumption analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('input_csv', type=str,
                       help='Combined CSV file with all predictor results')
    parser.add_argument('--output-dir', type=str, default='delta_curves',
                       help='Output directory for plots and data (default: delta_curves)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input_csv).exists():
        print(f"Error: Input file '{args.input_csv}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Load data
    print(f"Loading combined predictor data from {args.input_csv}...")
    try:
        df = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded {len(df)} entries for {df['predictor_name'].nunique()} predictors")
    print(f"Predictors: {', '.join(sorted(df['predictor_name'].unique()))}\n")
    
    # Check required columns
    required_cols = ['trace_name', 'predictor_name', 'MPKI']
    missing = set(required_cols) - set(df.columns)
    if missing:
        print(f"Error: Missing required columns: {missing}", file=sys.stderr)
        sys.exit(1)
    
    # Generate curves
    generator = DeltaCurveGenerator(df, args.output_dir)
    generator.generate_all_curves()


if __name__ == '__main__':
    main()
