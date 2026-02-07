"""
Automated subsumption analysis runner for all predictor result files.

This script automatically discovers and processes all predictor CSV files
in the current directory, combines them, and runs delta-based subsumption analysis.

Usage:
    python run_subsumption_analysis.py [--output-dir <dir>]
"""

import pandas as pd
import argparse
import sys
from pathlib import Path

# Import the analyzer from the main script
from delta_subsumption_analysis import SubsumptionAnalyzer


def discover_csv_files() -> list:
    """
    Discover all predictor CSV files in the current directory.
    
    Returns:
        List of Path objects for relevant CSV files
    """
    current_dir = Path('.')
    csv_files = []
    
    # Files to process
    target_files = [
        'baseline.csv',
        'load-value-correlator-man-results.csv',
        'register-value-aware-toru-results.csv',
        'tage-sc-l-alberto-ros-results.csv',
        'tage-scl-andrez-seznec-results.csv'
    ]
    
    for filename in target_files:
        filepath = current_dir / filename
        if filepath.exists():
            csv_files.append(filepath)
            print(f"✓ Found: {filename}")
    
    return csv_files


def extract_predictor_name(filename: str) -> str:
    """
    Extract a clean predictor name from the CSV filename.
    
    Args:
        filename: CSV filename
        
    Returns:
        Clean predictor name
    """
    name_map = {
        'baseline.csv': 'Baseline',
        'load-value-correlator-man-results.csv': 'LVCP',
        'register-value-aware-toru-results.csv': 'RVA-Toru',
        'tage-sc-l-alberto-ros-results.csv': 'TAGE-SC-L',
        'tage-scl-andrez-seznec-results.csv': 'TAGE-SCL'
    }
    
    return name_map.get(filename, Path(filename).stem)


def load_and_combine_csvs(csv_files: list) -> pd.DataFrame:
    """
    Load all CSV files and combine them into a single DataFrame.
    
    Args:
        csv_files: List of Path objects for CSV files
        
    Returns:
        Combined DataFrame with predictor_name column added
    """
    dfs = []
    
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            
            # Strip whitespace from all string columns
            for col in df.columns:
                if df[col].dtype == object:
                    df[col] = df[col].astype(str).str.strip()
            
            # Strip column names
            df.columns = df.columns.str.strip()
            
            # Add predictor name if not present
            if 'predictor_name' not in df.columns:
                predictor_name = extract_predictor_name(csv_file.name)
                df['predictor_name'] = predictor_name
            
            dfs.append(df)
            print(f"  Loaded {len(df)} entries from {csv_file.name}")
            
        except Exception as e:
            print(f"  Error loading {csv_file.name}: {e}", file=sys.stderr)
            continue
    
    if not dfs:
        raise ValueError("No CSV files could be loaded successfully")
    
    combined = pd.concat(dfs, ignore_index=True)
    
    # Strip whitespace from combined column names
    combined.columns = combined.columns.str.strip()
    
    # Use Run column as the unique trace identifier
    if 'Run' in combined.columns:
        combined['trace_name'] = combined['Run'].astype(str).str.strip()
    else:
        raise ValueError("Required 'Run' column not found in CSV files")
    
    # Normalize remaining column names
    if 'MR' in combined.columns:
        combined.rename(columns={'MR': 'miss_rate'}, inplace=True)
    
    # Ensure MPKI is numeric
    combined['MPKI'] = pd.to_numeric(combined['MPKI'], errors='coerce')
    
    # Ensure required columns exist
    required_cols = ['trace_name', 'predictor_name', 'MPKI']
    
    missing = set(required_cols) - set(combined.columns)
    if missing:
        print(f"\nAvailable columns: {list(combined.columns)}")
        raise ValueError(f"Missing required columns after mapping: {missing}")
    
    # Select only the columns we need for analysis
    analysis_cols = ['trace_name', 'predictor_name', 'MPKI', 'miss_rate', 'IPC', 'Cycles']
    available_analysis_cols = [col for col in analysis_cols if col in combined.columns]
    
    result = combined[available_analysis_cols].copy()
    
    return result


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Automated subsumption analysis for all predictor results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--output-dir', type=str, default='subsumption_analysis',
                       help='Output directory for analysis results (default: subsumption_analysis)')
    
    args = parser.parse_args()
    
    print("="*80)
    print("AUTOMATED SUBSUMPTION ANALYSIS")
    print("="*80)
    
    # Discover CSV files
    print("\nDiscovering predictor CSV files...")
    csv_files = discover_csv_files()
    
    if not csv_files:
        print("\nError: No CSV files found in current directory", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nFound {len(csv_files)} CSV file(s)")
    
    # Load and combine data
    print("\nLoading and combining data...")
    try:
        combined_df = load_and_combine_csvs(csv_files)
    except Exception as e:
        print(f"\nError combining CSV files: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nCombined dataset:")
    print(f"  Total entries: {len(combined_df)}")
    print(f"  Predictors: {combined_df['predictor_name'].nunique()}")
    print(f"  Traces: {combined_df['trace_name'].nunique()}")
    print(f"  Predictor names: {', '.join(sorted(combined_df['predictor_name'].unique()))}")
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save combined dataset
    combined_file = output_dir / 'combined_predictor_results.csv'
    combined_df.to_csv(combined_file, index=False)
    print(f"\n✓ Saved combined dataset to: {combined_file}")
    
    # Initialize analyzer
    print("\n" + "="*80)
    print("PERFORMING SUBSUMPTION ANALYSIS")
    print("="*80)
    
    analyzer = SubsumptionAnalyzer(combined_df)
    
    # Perform analysis
    results = analyzer.analyze_all_pairs()
    
    # Print detailed results
    for result in results:
        analyzer.print_analysis(result)
    
    # Export delta values
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
            'Avg_Delta_MPKI': result.get('avg_delta', float('nan')),
            'Pct_B_Better': result.get('pct_improved', float('nan')),
            'Worst_MPKI_A': result.get('worst_case_a', float('nan')),
            'Worst_MPKI_B': result.get('worst_case_b', float('nan')),
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
    print(f"\nGenerated files:")
    print(f"  - combined_predictor_results.csv")
    print(f"  - pairwise_deltas.csv")
    print(f"  - subsumption_summary.csv")


if __name__ == '__main__':
    main()
