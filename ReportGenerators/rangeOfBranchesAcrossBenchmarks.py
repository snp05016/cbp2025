# Compare Branch Density (BrPerCyc) across benchmark categories
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - no popups
import matplotlib.pyplot as plt
import numpy as np


def load_file(file_path):
    data = pd.read_csv(file_path)
    return data


def get_category(run_name):
    """Extract category from run name (e.g., 'int_11_trace' -> 'int')"""
    return run_name.split('_')[0]


def analyze_branch_density(data):
    """
    Analyze Branch Density (BrPerCyc) across different benchmark categories.
    
    Returns:
        DataFrame with statistics for each category
    """
    # Add category column
    data['Category'] = data['Run'].apply(get_category)
    
    # Group by category and calculate statistics
    categories = ['int', 'fp', 'compress', 'infra', 'media', 'web']
    results = []
    
    for category in categories:
        cat_data = data[data['Category'] == category]['BrPerCyc']
        
        if len(cat_data) == 0:
            continue
        
        stats = {
            'Category': category.upper(),
            'Count': len(cat_data),
            'Median': cat_data.median(),
            'Mean': cat_data.mean(),
            'Min': cat_data.min(),
            'Max': cat_data.max(),
            'Range': cat_data.max() - cat_data.min(),
            'Std': cat_data.std(),
            'Q1': cat_data.quantile(0.25),
            'Q3': cat_data.quantile(0.75),
            'IQR': cat_data.quantile(0.75) - cat_data.quantile(0.25)
        }
        results.append(stats)
    
    return pd.DataFrame(results), data


def print_analysis(stats_df, data):
    """Print detailed branch density analysis."""
    
    print("\n" + "="*120)
    print(f"{'BRANCH DENSITY ANALYSIS: BrPerCyc (Branches Per Cycle)':^120}")
    print("="*120)
    
    print("\n" + "WHAT THIS METRIC TELLS YOU:")
    print("-"*120)
    print("  BrPerCyc measures how control-flow heavy a workload is:")
    print("    • HIGH BrPerCyc   → Tight loops, simple branches, high branch density")
    print("    • LOW BrPerCyc    → Control-heavy workloads with complex branches, fewer branches executed per cycle")
    print("    • MODERATE values → Balanced control flow")
    print("\n" + "-"*120)
    
    # Sort by median for better readability
    stats_df = stats_df.sort_values('Median', ascending=False)
    
    print(f"\n{'STATISTICAL SUMMARY BY CATEGORY':^120}")
    print("="*120)
    print(f"{'Category':<12} {'Count':>7} {'Median':>10} {'Mean':>10} {'Min':>10} {'Max':>10} "
          f"{'Range':>10} {'Std':>10} {'Q1':>10} {'Q3':>10}")
    print("-"*120)
    
    for _, row in stats_df.iterrows():
        print(f"{row['Category']:<12} {row['Count']:>7.0f} {row['Median']:>10.4f} {row['Mean']:>10.4f} "
              f"{row['Min']:>10.4f} {row['Max']:>10.4f} {row['Range']:>10.4f} {row['Std']:>10.4f} "
              f"{row['Q1']:>10.4f} {row['Q3']:>10.4f}")
    
    print("-"*120)
    
    # Calculate and print geometric mean using formula: geomean = (x1 * x2 * ... * xn)^(1/n)
    # Or equivalently: geomean = exp(mean(log(x)))
    
    def calculate_geomean(values):
        """Calculate geometric mean using logarithms, filtering out zeros and negatives"""
        filtered = values[values > 0]
        if len(filtered) == 0:
            return 0.0
        return np.exp(np.mean(np.log(filtered)))
    
    geomean_median = calculate_geomean(stats_df['Median'])
    geomean_mean = calculate_geomean(stats_df['Mean'])
    geomean_min = calculate_geomean(stats_df['Min'])
    geomean_max = calculate_geomean(stats_df['Max'])
    geomean_range = calculate_geomean(stats_df['Range'])
    geomean_std = calculate_geomean(stats_df['Std'])
    geomean_q1 = calculate_geomean(stats_df['Q1'])
    geomean_q3 = calculate_geomean(stats_df['Q3'])
    total_count = stats_df['Count'].sum()
    
    print(f"{'GEOMEAN':<12} {total_count:>7.0f} {geomean_median:>10.4f} {geomean_mean:>10.4f} "
          f"{geomean_min:>10.4f} {geomean_max:>10.4f} {geomean_range:>10.4f} {geomean_std:>10.4f} "
          f"{geomean_q1:>10.4f} {geomean_q3:>10.4f}")
    
    print("="*120)
    
    # Interpretation
    print(f"\n{'INTERPRETATION':^120}")
    print("="*120)
    
    for _, row in stats_df.iterrows():
        cat = row['Category']
        median = row['Median']
        range_val = row['Range']
        
        # Classify based on median
        if median >= 0.6:
            density = "HIGH"
            interpretation = "Tight loops with high branch density - branches execute frequently"
        elif median >= 0.4:
            density = "MODERATE-HIGH"
            interpretation = "Balanced control flow with good branch execution rate"
        elif median >= 0.25:
            density = "MODERATE"
            interpretation = "Moderate control flow - mix of branches and other instructions"
        else:
            density = "LOW"
            interpretation = "Control-heavy with complex branches - fewer branches per cycle"
        
        # Assess consistency
        if range_val < 0.2:
            consistency = "CONSISTENT"
        elif range_val < 0.4:
            consistency = "MODERATE VARIATION"
        else:
            consistency = "HIGH VARIATION"
        
        print(f"\n{cat}:")
        print(f"  Density Level:  {density} (Median = {median:.4f})")
        print(f"  Consistency:    {consistency} (Range = {range_val:.4f})")
        print(f"  Interpretation: {interpretation}")
    
    print("\n" + "="*120)
    
    # Calculate overall geometric means across ALL benchmarks
    print(f"\n{'OVERALL GEOMETRIC MEANS (ALL CATEGORIES)':^120}")
    print("="*120)
    
    def calculate_geomean(values):
        """Calculate geometric mean using logarithms, filtering out zeros and negatives"""
        # Filter out zero and negative values
        filtered = values[values > 0]
        if len(filtered) == 0:
            return 0.0
        return np.exp(np.mean(np.log(filtered)))
    
    # Calculate geomeans for key metrics across all benchmarks
    all_brpercyc = data['BrPerCyc'].values
    all_numbr = data['NumBr'].values
    all_cycles = data['Cycles'].values
    all_mpki = data['MPKI'].values
    
    geomean_brpercyc = calculate_geomean(all_brpercyc)
    geomean_numbr = calculate_geomean(all_numbr)
    geomean_cycles = calculate_geomean(all_cycles)
    geomean_mpki = calculate_geomean(all_mpki)
    
    print(f"\n{'Metric':<20} {'Geometric Mean':>20} {'Sample Size':>15}")
    print("-"*120)
    print(f"{'BrPerCyc':<20} {geomean_brpercyc:>20.6f} {len(all_brpercyc):>15.0f}")
    print(f"{'NumBr':<20} {geomean_numbr:>20.2f} {len(all_numbr):>15.0f}")
    print(f"{'Cycles':<20} {geomean_cycles:>20.2f} {len(all_cycles):>15.0f}")
    print(f"{'MPKI':<20} {geomean_mpki:>20.6f} {len(all_mpki):>15.0f}")
    print("="*120)
    
    # Also show per-category geometric means for these metrics
    print(f"\n{'PER-CATEGORY GEOMETRIC MEANS':^120}")
    print("="*120)
    print(f"{'Category':<12} {'BrPerCyc':>15} {'NumBr':>15} {'Cycles':>15} {'MPKI':>15}")
    print("-"*120)
    
    for category in stats_df['Category'].str.lower():
        cat_data = data[data['Category'] == category]
        
        if len(cat_data) == 0:
            continue
        
        cat_brpercyc = calculate_geomean(cat_data['BrPerCyc'].values)
        cat_numbr = calculate_geomean(cat_data['NumBr'].values)
        cat_cycles = calculate_geomean(cat_data['Cycles'].values)
        cat_mpki = calculate_geomean(cat_data['MPKI'].values)
        
        print(f"{category.upper():<12} {cat_brpercyc:>15.6f} {cat_numbr:>15.2f} {cat_cycles:>15.2f} {cat_mpki:>15.6f}")
    
    print("="*120)
    
    # Detailed breakdown by benchmark
    print(f"\n{'DETAILED BREAKDOWN BY BENCHMARK':^120}")
    print("="*120)
    
    data_sorted = data.sort_values(['Category', 'BrPerCyc'], ascending=[True, False])
    
    for category in stats_df['Category'].str.lower():
        cat_data = data_sorted[data_sorted['Category'] == category]
        
        if len(cat_data) == 0:
            continue
        
        print(f"\n{category.upper()} Benchmarks (n={len(cat_data)}):")
        print("-"*120)
        print(f"{'Rank':<6} {'Benchmark':<30} {'BrPerCyc':>12} {'NumBr':>15} {'Cycles':>15} {'MPKI':>10}")
        print("-"*120)
        
        for i, (_, row) in enumerate(cat_data.iterrows(), 1):
            print(f"{i:<6} {row['Run'].replace('_trace', ''):<30} {row['BrPerCyc']:>12.4f} "
                  f"{row['NumBr']:>15.0f} {row['Cycles']:>15.0f} {row['MPKI']:>10.4f}")
    
    print("\n" + "="*120)


def plot_branch_density(stats_df, data):
    """Create visualizations for branch density analysis."""
    
    # Sort by median for consistent ordering
    stats_df = stats_df.sort_values('Median', ascending=False)
    categories = stats_df['Category'].tolist()
    
    # Calculate geometric mean and std deviation for each category
    def calculate_geomean(values):
        """Calculate geometric mean using logarithms, filtering out zeros and negatives"""
        filtered = values[values > 0]
        if len(filtered) == 0:
            return 0.0
        return np.exp(np.mean(np.log(filtered)))
    
    geomeans = []
    stds = []
    
    for category in categories:
        cat_data = data[data['Category'].str.upper() == category]['BrPerCyc']
        geomeans.append(calculate_geomean(cat_data.values))
        stds.append(cat_data.std())
    
    # Create single plot
    fig, ax = plt.subplots(figsize=(12, 7))
    
    x_pos = range(len(categories))
    colors = plt.cm.Set3(range(len(categories)))
    
    bars = ax.bar(x_pos, geomeans, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    ax.errorbar(x_pos, geomeans, yerr=stds, fmt='none', ecolor='black', 
                capsize=6, linewidth=2.5, capthick=2, label='±1 Std Dev')
    
    ax.set_ylabel('BrPerCyc Geometric Mean', fontsize=13, fontweight='bold')
    ax.set_xlabel('Benchmark Category', fontsize=13, fontweight='bold')
    ax.set_title('Branch Density (BrPerCyc): Geometric Mean by Category', 
                 fontsize=15, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(categories, rotation=0, fontsize=11, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='upper right')
    
    # Add value labels on bars
    for bar, val in zip(bars, geomeans):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height, f'{val:.4f}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    
    # Save to Reports directory
    import os
    output_dir = '../Reports/02_branch_density/graphs'
    os.makedirs(output_dir, exist_ok=True)
    filename = f'{output_dir}/branch_density_geomean.png'
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close figure instead of showing
    
    print(f"\nVisualization saved as: {filename}")


def main():
    """Main function to analyze branch density."""
    
    # Load results
    file_path = 'results.csv'
    data = load_file(file_path)
    
    # Analyze
    stats_df, data_with_category = analyze_branch_density(data)
    
    # Print analysis
    print_analysis(stats_df, data_with_category)
    
    # Create visualization
    plot_branch_density(stats_df, data_with_category)
    
    return stats_df


if __name__ == "__main__":
    stats = main()