# Find benchmarks within same category with similar metric but very different another metric
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - no popups
import matplotlib.pyplot as plt


def load_file(file_path):
    data = pd.read_csv(file_path)
    return data


def get_category(run_name):
    """Extract category from run name (e.g., 'int_11_trace' -> 'int')"""
    return run_name.split('_')[0]


def compare_benchmarks(data, category, similar_metric, different_metric, 
                       similarity_margin=None, difference_threshold=0.20):
    """
    Find benchmarks within a category with similar values for one metric 
    but very different values for another metric.
    
    Uses relative difference rule: |value₁ - value₂| / max(value₁, value₂) ≤ threshold
    
    Args:
        data: DataFrame with benchmark results
        category: Benchmark category to filter (e.g., 'int', 'fp', 'compress', 'infra', 'media', 'web')
        similar_metric: Metric that should be similar (e.g., 'MPKI')
        different_metric: Metric that should be different (e.g., 'Cycles')
        similarity_margin: Maximum % difference for similar_metric 
                          Default: 0.10-0.15 for INT, 0.05-0.10 for FP
        difference_threshold: Minimum % difference for different_metric (default: 0.20 = 20%)
    
    Returns:
        List of dictionaries containing comparison results
    """
    comparison_results = []
    
    # Set default similarity margins based on category
    if similarity_margin is None:
        if category == 'int':
            similarity_margin = 0.125  # 12.5% (midpoint of 10-15%)
        elif category == 'fp':
            similarity_margin = 0.075  # 7.5% (midpoint of 5-10%)
        else:
            similarity_margin = 0.10   # 10% for other categories
    
    # Filter data by category
    data['Category'] = data['Run'].apply(get_category)
    filtered_data = data[data['Category'] == category].copy()
    
    if len(filtered_data) < 2:
        print(f"Warning: Only {len(filtered_data)} benchmark(s) found in category '{category}'")
        return comparison_results
    
    # Get the run names and metrics
    runs = filtered_data['Run'].tolist()
    similar_values = filtered_data[similar_metric].tolist()
    different_values = filtered_data[different_metric].tolist()
    
    # Compare each benchmark with others in same category
    for i in range(len(runs)):
        for j in range(i + 1, len(runs)):
            # Relative difference rule for similar_metric: |val₁ - val₂| / max(val₁, val₂)
            similar_diff = abs(similar_values[i] - similar_values[j])
            similar_max = max(similar_values[i], similar_values[j])
            similar_percent = (similar_diff / similar_max) if similar_max != 0 else 0
            
            # Relative difference rule for different_metric
            different_diff = abs(different_values[i] - different_values[j])
            different_max = max(different_values[i], different_values[j])
            different_percent = (different_diff / different_max) if different_max != 0 else 0
            
            # Find pairs with similar metric AND very different metric
            if similar_percent <= similarity_margin and different_percent >= difference_threshold:
                comparison_results.append({
                    'category': category,
                    'run1': runs[i],
                    'run2': runs[j],
                    f'{similar_metric}_1': similar_values[i],
                    f'{similar_metric}_2': similar_values[j],
                    f'{similar_metric}_diff_percent': similar_percent * 100,
                    f'{different_metric}_1': different_values[i],
                    f'{different_metric}_2': different_values[j],
                    f'{different_metric}_diff': different_values[i] - different_values[j],
                    f'{different_metric}_diff_percent': different_percent * 100
                })
    
    return comparison_results


def plot_comparison(comparison_results, category, similar_metric, different_metric):
    """
    Plot comparison results showing discrepancies in different_metric 
    for benchmarks with similar similar_metric values.
    
    Args:
        comparison_results: List of comparison dictionaries
        category: Benchmark category
        similar_metric: Name of the similar metric
        different_metric: Name of the different metric
    """
    if not comparison_results:
        print(f"No discrepancies found in '{category}' benchmarks with similar {similar_metric}.")
        return
    
    # Sort by different_metric difference percentage (descending)
    comparison_results.sort(key=lambda x: abs(x[f'{different_metric}_diff_percent']), reverse=True)
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot 1: Bar chart comparing different_metric for benchmark pairs
    labels = [f"{r['run1'].replace('_trace', '')}\nvs\n{r['run2'].replace('_trace', '')}" 
              for r in comparison_results]
    x_pos = range(len(labels))
    
    diff_1 = [r[f'{different_metric}_1'] for r in comparison_results]
    diff_2 = [r[f'{different_metric}_2'] for r in comparison_results]
    
    bar_width = 0.35
    ax1.bar([x - bar_width/2 for x in x_pos], diff_1, bar_width, label='Benchmark 1', alpha=0.8, color='steelblue')
    ax1.bar([x + bar_width/2 for x in x_pos], diff_2, bar_width, label='Benchmark 2', alpha=0.8, color='coral')
    
    ax1.set_xlabel('Benchmark Pairs', fontsize=12)
    ax1.set_ylabel(different_metric, fontsize=12)
    ax1.set_title(f'{category.upper()} Benchmarks: {different_metric} Discrepancies with Similar {similar_metric}', 
                  fontsize=14, fontweight='bold')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Show both metrics side by side
    similar_1 = [r[f'{similar_metric}_1'] for r in comparison_results]
    similar_2 = [r[f'{similar_metric}_2'] for r in comparison_results]
    diff_percents = [r[f'{different_metric}_diff_percent'] for r in comparison_results]
    
    # Create secondary axis
    ax2_twin = ax2.twinx()
    
    # Plot similar metric (should be close)
    x_pairs = [[x - 0.2, x + 0.2] for x in x_pos]
    for i, (x_pair, s1, s2) in enumerate(zip(x_pairs, similar_1, similar_2)):
        ax2.plot(x_pair, [s1, s2], 'o-', color='green', alpha=0.6, linewidth=2, markersize=8)
    
    ax2.set_ylabel(f'{similar_metric} (similar)', fontsize=12, color='green')
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.set_xlabel('Benchmark Pairs', fontsize=12)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax2.grid(axis='y', alpha=0.3)
    
    # Plot difference percentage on right axis
    colors = ['red' if d > 50 else 'orange' for d in diff_percents]
    ax2_twin.bar(x_pos, diff_percents, alpha=0.3, color=colors, width=0.6)
    ax2_twin.set_ylabel(f'{different_metric} Difference %', fontsize=12, color='red')
    ax2_twin.tick_params(axis='y', labelcolor='red')
    
    ax2.set_title(f'Verification: {similar_metric} is Similar (green lines) but {different_metric} is Very Different (bars)', 
                  fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    # Save to Reports directory
    import os
    output_dir = '../Reports/05_comparative_analysis/divergent_behavior'
    os.makedirs(output_dir, exist_ok=True)
    filename = f'{output_dir}/discrepancy_{category}_{different_metric}_vs_{similar_metric}.png'
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close figure instead of showing
    
    print(f"\nPlot saved as: {filename}")


def main():
    """
    Main function to find discrepancies in benchmarks.
    """
    # Load results
    file_path = 'results.csv'
    data = load_file(file_path)
    
    print("Finding benchmark discrepancies...")
    print("=" * 80)
    
    # Configuration - CHANGE THESE VALUES AS NEEDED
    category = 'fp'              # Category: 'int', 'fp', 'compress', 'infra', 'media', 'web'
    similar_metric = 'MPKI'       # Metric that should be similar
    different_metric = 'CycWPPKI'   # Metric that should be very different
    similarity_margin = None      # Auto: 10-15% for INT, 5-10% for FP (or set manually)
    difference_threshold = 0.20   # 20% minimum difference for different metric
    
    comparison_results = compare_benchmarks(
        data, 
        category, 
        similar_metric, 
        different_metric,
        similarity_margin,
        difference_threshold
    )
    
    # Print results
    print(f"\nCategory: {category.upper()}")
    print(f"Finding benchmarks with:")
    if similarity_margin is None:
        margin_text = "7.5% for FP (auto)"
    else:
        margin_text = f"{similarity_margin*100}%"
    print(f"  - Similar {similar_metric} (±{margin_text}, using |val₁-val₂|/max formula)")
    print(f"  - Very different {different_metric} (at least {difference_threshold*100}% difference)")
    print("-" * 80)
    print(f"\nFound {len(comparison_results)} discrepancy pairs:")
    
    for i, result in enumerate(comparison_results, 1):
        print(f"\n{i}. {result['run1']} vs {result['run2']}:")
        print(f"   {similar_metric}: {result[f'{similar_metric}_1']:.4f} vs {result[f'{similar_metric}_2']:.4f} " +
              f"(diff: {result[f'{similar_metric}_diff_percent']:.2f}%)")
        print(f"   {different_metric}: {result[f'{different_metric}_1']:.0f} vs {result[f'{different_metric}_2']:.0f} " +
              f"(diff: {result[f'{different_metric}_diff_percent']:.2f}%)")
    
    # Plot results
    if comparison_results:
        plot_comparison(comparison_results, category, similar_metric, different_metric)
    
    return comparison_results


if __name__ == "__main__":
    results = main()
        