# Compare MPKI vs 50PercMPKI for benchmarks within same category
import pandas as pd
import matplotlib.pyplot as plt


def load_file(file_path):
    data = pd.read_csv(file_path)
    return data


def get_category(run_name):
    """Extract category from run name (e.g., 'int_11_trace' -> 'int')"""
    return run_name.split('_')[0]


def analyze_mpki_consistency(data, category, similarity_threshold=0.10, difference_threshold=0.20):
    """
    Analyze MPKI vs 50PercMPKI consistency within a benchmark category.
    
    Finds benchmarks where:
    1. MPKI ≈ 50PercMPKI (consistent behavior throughout execution)
    2. MPKI ≠ 50PercMPKI (phase behavior - different in first vs second half)
    
    Args:
        data: DataFrame with benchmark results
        category: Benchmark category ('int', 'fp', 'compress', 'infra', 'media', 'web')
        similarity_threshold: Max % difference to consider values similar (default: 10%)
        difference_threshold: Min % difference to consider values different (default: 20%)
    
    Returns:
        Dictionary with 'similar' and 'different' benchmark lists
    """
    # Filter data by category
    data['Category'] = data['Run'].apply(get_category)
    filtered_data = data[data['Category'] == category].copy()
    
    if len(filtered_data) == 0:
        print(f"Warning: No benchmarks found in category '{category}'")
        return {'similar': [], 'different': []}
    
    similar_benchmarks = []
    different_benchmarks = []
    
    for _, row in filtered_data.iterrows():
        run = row['Run']
        mpki = row['MPKI']
        mpki_50perc = row['50PercMPKI']
        
        # Calculate relative difference using max formula: |val₁ - val₂| / max(val₁, val₂)
        diff = abs(mpki - mpki_50perc)
        max_val = max(mpki, mpki_50perc)
        
        if max_val == 0:
            relative_diff = 0
        else:
            relative_diff = diff / max_val
        
        benchmark_info = {
            'run': run,
            'MPKI': mpki,
            '50PercMPKI': mpki_50perc,
            'abs_diff': diff,
            'relative_diff_percent': relative_diff * 100,
            'Cycles': row['Cycles'],
            '50PercCycles': row['50PercCycles'],
            'MR': row['MR'],
            '50PercMR': row['50PercMR']
        }
        
        # Categorize based on relative difference
        if relative_diff <= similarity_threshold:
            similar_benchmarks.append(benchmark_info)
        elif relative_diff >= difference_threshold:
            different_benchmarks.append(benchmark_info)
    
    # Sort by relative difference
    similar_benchmarks.sort(key=lambda x: x['relative_diff_percent'])
    different_benchmarks.sort(key=lambda x: x['relative_diff_percent'], reverse=True)
    
    return {
        'similar': similar_benchmarks,
        'different': different_benchmarks
    }


def print_analysis_results(results, category, similarity_threshold, difference_threshold):
    """Print detailed analysis results to CLI."""
    
    print("\n" + "="*100)
    print(f"MPKI CONSISTENCY ANALYSIS FOR {category.upper()} BENCHMARKS")
    print("="*100)
    print(f"\nCriteria:")
    print(f"  - SIMILAR (≈): Relative difference ≤ {similarity_threshold*100}% (consistent throughout execution)")
    print(f"  - DIFFERENT (≠): Relative difference ≥ {difference_threshold*100}% (phase behavior detected)")
    print(f"  - Formula: |MPKI - 50PercMPKI| / max(MPKI, 50PercMPKI)")
    print("-"*100)
    
    similar = results['similar']
    different = results['different']
    
    # Print SIMILAR benchmarks (MPKI ≈ 50PercMPKI)
    print(f"\n{'CONSISTENT BENCHMARKS (MPKI ≈ 50PercMPKI)':^100}")
    print("="*100)
    print(f"Found {len(similar)} benchmark(s) with consistent MPKI throughout execution:\n")
    
    if similar:
        print(f"{'Rank':<6} {'Benchmark':<25} {'MPKI':>10} {'50%MPKI':>10} {'Diff%':>10} {'MR':>12} {'50%MR':>12}")
        print("-"*100)
        for i, bench in enumerate(similar, 1):
            print(f"{i:<6} {bench['run'].replace('_trace', ''):<25} "
                  f"{bench['MPKI']:>10.4f} {bench['50PercMPKI']:>10.4f} "
                  f"{bench['relative_diff_percent']:>9.2f}% "
                  f"{bench['MR']:>12} {bench['50PercMR']:>12}")
        
        print(f"\n{'INSIGHT':>10}: These benchmarks show uniform misprediction behavior across the entire trace.")
        print(f"{'':>10}  No significant phase changes detected.")
    else:
        print(f"  No benchmarks found with MPKI difference ≤ {similarity_threshold*100}%")
    
    # Print DIFFERENT benchmarks (MPKI ≠ 50PercMPKI)
    print(f"\n\n{'PHASE-BEHAVIOR BENCHMARKS (MPKI ≠ 50PercMPKI)':^100}")
    print("="*100)
    print(f"Found {len(different)} benchmark(s) with significant phase behavior:\n")
    
    if different:
        print(f"{'Rank':<6} {'Benchmark':<25} {'MPKI':>10} {'50%MPKI':>10} {'Diff%':>10} {'MR':>12} {'50%MR':>12}")
        print("-"*100)
        for i, bench in enumerate(different, 1):
            trend = "↑ Harder" if bench['50PercMPKI'] > bench['MPKI'] else "↓ Easier"
            print(f"{i:<6} {bench['run'].replace('_trace', ''):<25} "
                  f"{bench['MPKI']:>10.4f} {bench['50PercMPKI']:>10.4f} "
                  f"{bench['relative_diff_percent']:>9.2f}% "
                  f"{bench['MR']:>12} {bench['50PercMR']:>12}  {trend}")
        
        print(f"\n{'INSIGHT':>10}: These benchmarks show distinct phase changes:")
        print(f"{'':>10}  - '↑ Harder' = Second half has MORE mispredictions (becomes harder to predict)")
        print(f"{'':>10}  - '↓ Easier' = Second half has FEWER mispredictions (becomes easier to predict)")
    else:
        print(f"  No benchmarks found with MPKI difference ≥ {difference_threshold*100}%")
    
    # Summary statistics
    print(f"\n\n{'SUMMARY STATISTICS':^100}")
    print("="*100)
    print(f"Total {category.upper()} benchmarks analyzed: {len(similar) + len(different) + (len(results.get('neither', [])) if 'neither' in results else 0)}")
    print(f"  - Consistent behavior (≈):     {len(similar):>3} benchmarks ({len(similar)/(len(similar)+len(different))*100 if (len(similar)+len(different))>0 else 0:.1f}%)")
    print(f"  - Phase behavior (≠):          {len(different):>3} benchmarks ({len(different)/(len(similar)+len(different))*100 if (len(similar)+len(different))>0 else 0:.1f}%)")
    print("="*100 + "\n")


def plot_comparison(results, category):
    """Create simple but effective visualization."""
    
    similar = results['similar']
    different = results['different']
    
    if not similar and not different:
        print(f"No data to plot for {category} category.")
        return
    
    # Combine for overall view
    all_benchmarks = similar + different
    
    if not all_benchmarks:
        return
    
    # Sort by absolute difference for better visualization
    all_benchmarks.sort(key=lambda x: x['relative_diff_percent'], reverse=True)
    
    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: MPKI vs 50PercMPKI comparison
    runs = [b['run'].replace('_trace', '') for b in all_benchmarks]
    mpki = [b['MPKI'] for b in all_benchmarks]
    mpki_50 = [b['50PercMPKI'] for b in all_benchmarks]
    
    x = range(len(runs))
    width = 0.35
    
    bars1 = ax1.bar([i - width/2 for i in x], mpki, width, label='MPKI (Full)', 
                    color='steelblue', alpha=0.8)
    bars2 = ax1.bar([i + width/2 for i in x], mpki_50, width, label='50PercMPKI (2nd Half)', 
                    color='coral', alpha=0.8)
    
    ax1.set_xlabel('Benchmark', fontsize=12, fontweight='bold')
    ax1.set_ylabel('MPKI', fontsize=12, fontweight='bold')
    ax1.set_title(f'{category.upper()} Benchmarks: MPKI vs 50PercMPKI Comparison', 
                  fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(runs, rotation=45, ha='right', fontsize=9)
    ax1.legend(fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Relative difference percentage
    diff_percents = [b['relative_diff_percent'] for b in all_benchmarks]
    colors = ['green' if d <= 10 else 'orange' if d <= 20 else 'red' for d in diff_percents]
    
    bars = ax2.bar(x, diff_percents, color=colors, alpha=0.7)
    ax2.axhline(y=10, color='green', linestyle='--', linewidth=1, label='10% (Similar threshold)')
    ax2.axhline(y=20, color='red', linestyle='--', linewidth=1, label='20% (Different threshold)')
    
    ax2.set_xlabel('Benchmark', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Relative Difference %', fontsize=12, fontweight='bold')
    ax2.set_title('MPKI Consistency: Lower = More Consistent Throughout Execution', 
                  fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(runs, rotation=45, ha='right', fontsize=9)
    ax2.legend(fontsize=9)
    ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    filename = f'mpki_consistency_{category}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Plot saved as: {filename}")


def main():
    """Main function to analyze MPKI consistency."""
    
    # Load results
    file_path = 'results.csv'
    data = load_file(file_path)
    
    # Configuration - CHANGE THESE VALUES AS NEEDED
    category = 'int'              # Category: 'int', 'fp', 'compress', 'infra', 'media', 'web'
    similarity_threshold = 0.10   # 10% - consider similar/consistent
    difference_threshold = 0.20   # 20% - consider different/phase behavior
    
    # Analyze
    results = analyze_mpki_consistency(data, category, similarity_threshold, difference_threshold)
    
    # Print detailed results to CLI
    print_analysis_results(results, category, similarity_threshold, difference_threshold)
    
    # Create visualization
    plot_comparison(results, category)
    
    return results


if __name__ == "__main__":
    results = main()
