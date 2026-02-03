# Analyze Baseline Difficulty (Intrinsic Hardness) across benchmark categories
import pandas as pd
import numpy as np


def load_file(file_path):
    data = pd.read_csv(file_path)
    return data


def get_category(run_name):
    """Extract category from run name (e.g., 'int_11_trace' -> 'int')"""
    return run_name.split('_')[0]


def analyze_baseline_difficulty(data):
    """
    Analyze baseline difficulty metrics: MPKI, CycWPPKI, MR
    These metrics show intrinsic hardness of benchmarks.
    """
    # Add category column
    data['Category'] = data['Run'].apply(get_category)
    
    # Categories to analyze
    categories = ['int', 'fp', 'compress', 'infra', 'media', 'web']
    
    print("\n" + "="*130)
    print(f"{'BASELINE DIFFICULTY ANALYSIS (INTRINSIC HARDNESS)':^130}")
    print("="*130)
    
    print("\nMETRICS EXPLANATION:")
    print("-"*130)
    print("  • MPKI (Mispredictions Per Kilo Instructions) → Frequency of wrong predictions")
    print("  • CycWPPKI (Cycles Wrong-Path Per Kilo Instructions) → Cost of each mistake")
    print("  • MR (Misprediction Rate) → Confirms difficulty is not due to branch density")
    print("="*130)
    
    # Analyze each category
    for category in categories:
        cat_data = data[data['Category'] == category].copy()
        
        if len(cat_data) == 0:
            continue
        
        # Sort by MPKI (descending) - primary difficulty metric
        cat_data_sorted = cat_data.sort_values('MPKI', ascending=False)
        
        print(f"\n{'─'*130}")
        print(f"{category.upper()} BENCHMARKS (n={len(cat_data)}) - Ranked by MPKI (Hardest → Easiest)")
        print(f"{'─'*130}")
        print(f"{'Rank':<6} {'Benchmark':<30} {'MPKI':>12} {'CycWPPKI':>12} {'MR':>12} {'BrPerCyc':>12} "
              f"{'Cycles':>15} {'IPC':>8}")
        print("-"*130)
        
        for rank, (_, row) in enumerate(cat_data_sorted.iterrows(), 1):
            print(f"{rank:<6} {row['Run'].replace('_trace', ''):<30} "
                  f"{row['MPKI']:>12.4f} {row['CycWPPKI']:>12.4f} {row['MR']:>12} "
                  f"{row['BrPerCyc']:>12.4f} {row['Cycles']:>15,.0f} {row['IPC']:>8.4f}")
        
        # Calculate statistics for this category
        print("-"*130)
        print(f"{'CATEGORY STATS':<36} "
              f"{'Mean:':>6} {cat_data['MPKI'].mean():>5.4f} "
              f"{'Mean:':>6} {cat_data['CycWPPKI'].mean():>5.2f} "
              f"{'Mean:':>6} {cat_data['MR'].apply(lambda x: float(x.strip('%'))).mean():>5.2f}% "
              f"{'Mean:':>6} {cat_data['BrPerCyc'].mean():>5.4f}")
        
        def calculate_geomean(values):
            """Calculate geometric mean using logarithms, filtering out zeros and negatives"""
            filtered = values[values > 0]
            if len(filtered) == 0:
                return 0.0
            return np.exp(np.mean(np.log(filtered)))
        
        geomean_mpki = calculate_geomean(cat_data['MPKI'])
        geomean_cycwppki = calculate_geomean(cat_data['CycWPPKI'])
        mr_values = cat_data['MR'].apply(lambda x: float(x.strip('%')))
        geomean_mr = calculate_geomean(mr_values)
        geomean_brpercyc = calculate_geomean(cat_data['BrPerCyc'])
        
        print(f"{'':>36} "
              f"{'GMean:':>6} {geomean_mpki:>5.4f} "
              f"{'GMean:':>6} {geomean_cycwppki:>5.2f} "
              f"{'GMean:':>6} {geomean_mr:>5.2f}% "
              f"{'GMean:':>6} {geomean_brpercyc:>5.4f}")
        
        # Identify hardest benchmarks
        hardest = cat_data_sorted.iloc[0]
        easiest = cat_data_sorted.iloc[-1]
        print(f"\n  HARDEST: {hardest['Run'].replace('_trace', '')} (MPKI={hardest['MPKI']:.4f})")
        print(f"  EASIEST: {easiest['Run'].replace('_trace', '')} (MPKI={easiest['MPKI']:.4f})")
        print(f"  RANGE: {hardest['MPKI'] - easiest['MPKI']:.4f} MPKI difference")
    
    # Overall summary across all categories
    print(f"\n{'═'*130}")
    print(f"{'OVERALL SUMMARY (ALL CATEGORIES)':^130}")
    print(f"{'═'*130}")
    
    print(f"\n{'Category':<12} {'Count':>7} {'MPKI':>12} {'CycWPPKI':>12} {'MR (%)':>12} "
          f"{'BrPerCyc':>12} {'Avg Cycles':>15}")
    print("-"*130)
    
    summary_data = []
    for category in categories:
        cat_data = data[data['Category'] == category]
        
        if len(cat_data) == 0:
            continue
        
        def calculate_geomean(values):
            filtered = values[values > 0]
            if len(filtered) == 0:
                return 0.0
            return np.exp(np.mean(np.log(filtered)))
        
        geomean_mpki = calculate_geomean(cat_data['MPKI'])
        geomean_cycwppki = calculate_geomean(cat_data['CycWPPKI'])
        mr_values = cat_data['MR'].apply(lambda x: float(x.strip('%')))
        geomean_mr = calculate_geomean(mr_values)
        geomean_brpercyc = calculate_geomean(cat_data['BrPerCyc'])
        avg_cycles = cat_data['Cycles'].mean()
        
        summary_data.append({
            'Category': category.upper(),
            'Count': len(cat_data),
            'MPKI': geomean_mpki,
            'CycWPPKI': geomean_cycwppki,
            'MR': geomean_mr,
            'BrPerCyc': geomean_brpercyc,
            'Cycles': avg_cycles
        })
    
    # Sort by MPKI for overall comparison
    summary_df = pd.DataFrame(summary_data).sort_values('MPKI', ascending=False)
    
    for _, row in summary_df.iterrows():
        print(f"{row['Category']:<12} {row['Count']:>7.0f} {row['MPKI']:>12.4f} "
              f"{row['CycWPPKI']:>12.2f} {row['MR']:>12.2f} "
              f"{row['BrPerCyc']:>12.4f} {row['Cycles']:>15,.2f}")
    
    print("="*130)
    
    # Key insights
    print(f"\n{'KEY INSIGHTS':^130}")
    print("="*130)
    
    hardest_cat = summary_df.iloc[0]
    easiest_cat = summary_df.iloc[-1]
    
    print(f"\nHARDEST CATEGORY: {hardest_cat['Category']}")
    print(f"  → Geometric Mean MPKI: {hardest_cat['MPKI']:.4f}")
    print(f"  → Cost per mistake (CycWPPKI): {hardest_cat['CycWPPKI']:.2f}")
    print(f"  → Misprediction Rate: {hardest_cat['MR']:.2f}%")
    
    print(f"\nEASIEST CATEGORY: {easiest_cat['Category']}")
    print(f"  → Geometric Mean MPKI: {easiest_cat['MPKI']:.4f}")
    print(f"  → Cost per mistake (CycWPPKI): {easiest_cat['CycWPPKI']:.2f}")
    print(f"  → Misprediction Rate: {easiest_cat['MR']:.2f}%")
    
    print(f"\nDIFFICULTY RANGE:")
    print(f"  → MPKI varies from {easiest_cat['MPKI']:.4f} to {hardest_cat['MPKI']:.4f}")
    print(f"  → {hardest_cat['MPKI']/easiest_cat['MPKI']:.2f}x harder to predict")
    
    print("\n" + "="*130 + "\n")


def main():
    """Main function to analyze baseline difficulty."""
    
    # Load results
    file_path = 'results.csv'
    data = load_file(file_path)
    
    # Analyze
    analyze_baseline_difficulty(data)


if __name__ == "__main__":
    main()
