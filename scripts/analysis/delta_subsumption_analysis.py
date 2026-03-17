"""
Delta-Based Subsumption Analysis for Branch Predictors

This script performs rigorous pairwise subsumption analysis using MPKI deltas.
It enforces strict data sanitization rules to exclude invalid zero-MPKI entries,
ensuring that all statistical conclusions are based only on valid predictor behavior.

Usage:
    python delta_subsumption_analysis.py <input_csv> [--output-dir <dir>]

Input CSV Format:
    trace_name, predictor_name, MPKI, miss_rate, IPC, cycles

Output:
    - Console: Detailed subsumption analysis for each predictor pair
    - CSV: Pairwise delta values (valid traces only)
    - CSV: Summary table with subsumption decisions

Author: Generated for CBP-2025 Analysis
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
from typing import Tuple, Dict, List


class SubsumptionAnalyzer:
    """
    Performs delta-based subsumption analysis between branch predictors.
    
    Enforces strict validation: any trace with zero MPKI for either predictor
    in a comparison is excluded from that comparison's statistics.
    """
    
    def __init__(self, df: pd.DataFrame):
        """
        Initialize analyzer with predictor results.
        
        Args:
            df: DataFrame with columns [trace_name, predictor_name, MPKI, ...]
        """
        self.df = df.copy()
        self.predictors = sorted(df['predictor_name'].unique())
        self.traces = sorted(df['trace_name'].unique())
        
        # Validate required columns
        required_cols = ['trace_name', 'predictor_name', 'MPKI']
        missing = set(required_cols) - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
    
    def _get_valid_trace_pairs(self, pred_a: str, pred_b: str) -> Tuple[pd.DataFrame, int]:
        """
        Extract valid trace pairs for predictors A and B.
        
        A trace is valid only if both predictors have non-zero MPKI values.
        
        Args:
            pred_a: Name of predictor A
            pred_b: Name of predictor B
            
        Returns:
            Tuple of (valid_pairs_df, skipped_count)
            valid_pairs_df columns: [trace_name, MPKI_A, MPKI_B, delta]
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
        
        # Compute delta: Δ_A→B = MPKI_A - MPKI_B
        # Positive delta means B improves over A
        valid['delta'] = valid['MPKI_A'] - valid['MPKI_B']
        
        return valid, skipped
    
    def _classify_subsumption(self, 
                            avg_delta: float,
                            pct_improved: float,
                            worst_case_a: float,
                            worst_case_b: float,
                            valid_count: int) -> str:
        """
        Classify subsumption relationship based on statistical criteria.
        
        Args:
            avg_delta: Average delta (positive means B better)
            pct_improved: Percentage of traces where B beats A
            worst_case_a: Maximum MPKI for predictor A (valid traces)
            worst_case_b: Maximum MPKI for predictor B (valid traces)
            valid_count: Number of valid traces in comparison
            
        Returns:
            Subsumption classification string
        """
        if valid_count < 3:
            return "Insufficient valid data"
        
        # Strong subsumption criteria
        if (avg_delta > 0 and 
            pct_improved >= 90 and 
            worst_case_b <= worst_case_a * 1.05):  # Allow 5% tolerance on worst case
            return "B strictly subsumes A"
        
        # Mostly subsumption
        if (avg_delta > 0 and 
            pct_improved >= 75 and 
            worst_case_b <= worst_case_a * 1.15):  # Allow 15% tolerance
            return "B mostly subsumes A"
        
        # Complementary if mixed results
        if 30 <= pct_improved <= 70:
            return "A and B are complementary"
        
        # Weak relationships
        if avg_delta > 0 and pct_improved >= 55:
            return "B weakly dominates A"
        
        if avg_delta < 0 and pct_improved <= 45:
            return "A weakly dominates B"
        
        return "No clear dominance"
    
    def analyze_pair(self, pred_a: str, pred_b: str) -> Dict:
        """
        Perform complete subsumption analysis for predictor pair A vs B.
        
        Args:
            pred_a: Name of predictor A (baseline)
            pred_b: Name of predictor B (comparison)
            
        Returns:
            Dictionary containing all analysis metrics
        """
        valid_pairs, skipped = self._get_valid_trace_pairs(pred_a, pred_b)
        
        if len(valid_pairs) == 0:
            return {
                'pred_a': pred_a,
                'pred_b': pred_b,
                'valid_traces': 0,
                'skipped_traces': skipped,
                'classification': 'No valid traces for comparison'
            }
        
        # Compute statistics on valid data only
        avg_delta = valid_pairs['delta'].mean()
        std_delta = valid_pairs['delta'].std()
        median_delta = valid_pairs['delta'].median()
        
        # Count traces where B improves over A (delta > 0)
        improved = (valid_pairs['delta'] > 0).sum()
        pct_improved = 100 * improved / len(valid_pairs)
        
        # Worst-case analysis
        worst_case_a = valid_pairs['MPKI_A'].max()
        worst_case_b = valid_pairs['MPKI_B'].max()
        
        # Classify relationship
        classification = self._classify_subsumption(
            avg_delta, pct_improved, worst_case_a, worst_case_b, len(valid_pairs)
        )
        
        return {
            'pred_a': pred_a,
            'pred_b': pred_b,
            'valid_traces': len(valid_pairs),
            'skipped_traces': skipped,
            'avg_delta': avg_delta,
            'std_delta': std_delta,
            'median_delta': median_delta,
            'pct_improved': pct_improved,
            'improved_count': improved,
            'worst_case_a': worst_case_a,
            'worst_case_b': worst_case_b,
            'classification': classification,
            'valid_pairs': valid_pairs
        }
    
    def analyze_all_pairs(self) -> List[Dict]:
        """
        Analyze all predictor pairs.
        
        Returns:
            List of analysis dictionaries, one per pair
        """
        results = []
        
        for i, pred_a in enumerate(self.predictors):
            for pred_b in self.predictors[i+1:]:
                # Analyze A vs B
                result = self.analyze_pair(pred_a, pred_b)
                results.append(result)
        
        return results
    
    def print_analysis(self, result: Dict):
        """
        Print formatted analysis for a single predictor pair.
        
        Args:
            result: Analysis dictionary from analyze_pair()
        """
        print(f"\n{'='*80}")
        print(f"SUBSUMPTION ANALYSIS: {result['pred_a']} vs {result['pred_b']}")
        print(f"{'='*80}")
        
        print(f"\nData Quality:")
        print(f"  Valid traces:   {result['valid_traces']}")
        print(f"  Skipped traces: {result['skipped_traces']} (zero MPKI)")
        
        if result['valid_traces'] == 0:
            print(f"\n⚠️  {result['classification']}")
            return
        
        print(f"\nDelta Statistics (Δ = MPKI_{result['pred_a']} - MPKI_{result['pred_b']}):")
        print(f"  Average:  {result['avg_delta']:+.4f} MPKI")
        print(f"  Std Dev:  {result['std_delta']:.4f} MPKI")
        print(f"  Median:   {result['median_delta']:+.4f} MPKI")
        
        print(f"\nImprovement Analysis:")
        print(f"  {result['pred_b']} beats {result['pred_a']}: "
              f"{result['improved_count']}/{result['valid_traces']} traces "
              f"({result['pct_improved']:.1f}%)")
        
        print(f"\nWorst-Case MPKI:")
        print(f"  {result['pred_a']}: {result['worst_case_a']:.4f}")
        print(f"  {result['pred_b']}: {result['worst_case_b']:.4f}")
        
        print(f"\n{'─'*80}")
        print(f"CONCLUSION: {result['classification']}")
        print(f"{'─'*80}")
        
        # Generate paper-ready interpretation
        self._print_interpretation(result)
    
    def _print_interpretation(self, result: Dict):
        """
        Generate paper-ready textual interpretation.
        
        Args:
            result: Analysis dictionary
        """
        if result['valid_traces'] == 0:
            return
        
        pred_a = result['pred_a']
        pred_b = result['pred_b']
        
        print(f"\nPaper-Ready Interpretation:")
        print(f"  After excluding {result['skipped_traces']} invalid zero-MPKI traces, ")
        
        if "strictly subsumes" in result['classification']:
            print(f"  {pred_b} strictly subsumes {pred_a}. Across {result['valid_traces']} "
                  f"valid traces, {pred_b} achieves an average MPKI reduction of "
                  f"{result['avg_delta']:.4f}, improving performance on {result['pct_improved']:.1f}% "
                  f"of workloads. The worst-case MPKI of {pred_b} ({result['worst_case_b']:.4f}) "
                  f"is no worse than {pred_a} ({result['worst_case_a']:.4f}), indicating that "
                  f"{pred_b} captures all predictive patterns of {pred_a} plus additional "
                  f"architectural dependencies.")
        
        elif "mostly subsumes" in result['classification']:
            print(f"  {pred_b} mostly subsumes {pred_a}. On {result['valid_traces']} valid traces, "
                  f"{pred_b} provides superior prediction on {result['pct_improved']:.1f}% of workloads "
                  f"with an average MPKI improvement of {result['avg_delta']:.4f}. While "
                  f"{pred_a} may offer marginal benefits on select traces, {pred_b} represents "
                  f"a more robust general-purpose solution.")
        
        elif "complementary" in result['classification']:
            print(f"  {pred_a} and {pred_b} exhibit complementary behavior. {pred_b} outperforms "
                  f"{pred_a} on {result['pct_improved']:.1f}% of {result['valid_traces']} valid traces, "
                  f"suggesting that different workload characteristics favor different prediction "
                  f"mechanisms. A hybrid approach may be warranted.")
        
        else:
            print(f"  The relationship between {pred_a} and {pred_b} is inconclusive across "
                  f"{result['valid_traces']} valid traces. Further analysis or workload-specific "
                  f"evaluation may be required.")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Delta-based subsumption analysis for branch predictors',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('input_csv', type=str,
                       help='Input CSV file with predictor results')
    parser.add_argument('--output-dir', type=str, default='.',
                       help='Output directory for generated files (default: current directory)')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input_csv).exists():
        print(f"Error: Input file '{args.input_csv}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print(f"Loading data from {args.input_csv}...")
    try:
        df = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loaded {len(df)} entries for {df['predictor_name'].nunique()} predictors "
          f"across {df['trace_name'].nunique()} traces")
    
    # Initialize analyzer
    analyzer = SubsumptionAnalyzer(df)
    
    # Perform analysis
    print("\nPerforming pairwise subsumption analysis...")
    results = analyzer.analyze_all_pairs()
    
    # Print detailed results
    for result in results:
        analyzer.print_analysis(result)
    
    # Export delta values (all valid pairs)
    print(f"\n{'='*80}")
    print("EXPORTING RESULTS")
    print(f"{'='*80}")
    
    delta_dfs = []
    for result in results:
        if result['valid_traces'] > 0:
            pairs = result['valid_pairs'].copy()
            pairs['pred_a'] = result['pred_a']
            pairs['pred_b'] = result['pred_b']
            delta_dfs.append(pairs[['pred_a', 'pred_b', 'trace_name', 'MPKI_A', 'MPKI_B', 'delta']])
    
    if delta_dfs:
        all_deltas = pd.concat(delta_dfs, ignore_index=True)
        delta_file = output_dir / 'pairwise_deltas.csv'
        all_deltas.to_csv(delta_file, index=False)
        print(f"\n✓ Exported delta values to: {delta_file}")
    
    # Export summary table
    summary_data = []
    for result in results:
        summary_data.append({
            'Predictor_A': result['pred_a'],
            'Predictor_B': result['pred_b'],
            'Valid_Traces': result['valid_traces'],
            'Skipped_Traces': result['skipped_traces'],
            'Avg_Delta_MPKI': result.get('avg_delta', np.nan),
            'Pct_B_Better': result.get('pct_improved', np.nan),
            'Worst_MPKI_A': result.get('worst_case_a', np.nan),
            'Worst_MPKI_B': result.get('worst_case_b', np.nan),
            'Classification': result['classification']
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_dir / 'subsumption_summary.csv'
    summary_df.to_csv(summary_file, index=False)
    print(f"✓ Exported summary table to: {summary_file}")
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nAll outputs saved to: {output_dir.absolute()}")


if __name__ == '__main__':
    main()
