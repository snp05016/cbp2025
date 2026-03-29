#!/usr/bin/env python3
"""
Predictor Complementarity Analysis
Analyzes which traces are improved by which predictors and identifies
complementary predictor combinations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter
import sys

# Try to import upsetplot, fall back if not available
try:
    from upsetplot import UpSet, from_indicators
    UPSETPLOT_AVAILABLE = True
except ImportError:
    UPSETPLOT_AVAILABLE = False
    print("Warning: upsetplot not available. Will use alternative visualization.")

# Configuration
RESULTS_DIR = Path("results")
OUTPUT_DIR = Path("scripts/analysis/complementarity_analysis")
IMPROVEMENT_THRESHOLD = 0.1  # MPKI threshold for meaningful improvement

# Predictor file mappings
PREDICTOR_FILES = {
    'baseline': 'baseline.csv',
    'MPP': 'multi-perspective-predictor.csv',
    'Bullseye': 'bullseye-predictor.csv',
    'BALL': 'ball-fun.csv',
    'PIP': 'programming-idiom-predictor-mose.csv'
}

# Predictor names in order (excluding baseline)
PREDICTOR_NAMES = ['MPP', 'Bullseye', 'BALL', 'PIP']

def load_predictor_data():
    """Load all predictor CSV files and merge on trace name."""
    print("Loading predictor data...")

    # Load baseline
    baseline_df = pd.read_csv(RESULTS_DIR / PREDICTOR_FILES['baseline'])
    baseline_df = baseline_df.rename(columns={'MPKI': 'baseline_MPKI'})

    # Start with baseline Run column
    merged_df = baseline_df[['Run', 'Workload', 'baseline_MPKI']].copy()

    # Load and merge each predictor
    for pred_name in PREDICTOR_NAMES:
        pred_df = pd.read_csv(RESULTS_DIR / PREDICTOR_FILES[pred_name])
        pred_df = pred_df[['Run', 'MPKI']].rename(columns={'MPKI': f'{pred_name}_MPKI'})
        merged_df = merged_df.merge(pred_df, on='Run', how='inner')

    print(f"Loaded {len(merged_df)} traces")
    return merged_df

def compute_improvements(df):
    """Compute MPKI improvement over baseline for each predictor."""
    print("\nComputing improvements over baseline...")

    for pred_name in PREDICTOR_NAMES:
        # Improvement = baseline - predictor (positive means predictor is better)
        df[f'{pred_name}_improvement'] = df['baseline_MPKI'] - df[f'{pred_name}_MPKI']
        # Mark meaningful improvements
        df[f'{pred_name}_improved'] = df[f'{pred_name}_improvement'] >= IMPROVEMENT_THRESHOLD

    return df

def categorize_traces(df):
    """Categorize traces by which predictors improved them."""
    print("\nCategorizing traces by predictor combinations...")

    def get_category(row):
        improved = []
        for pred_name in PREDICTOR_NAMES:
            if row[f'{pred_name}_improved']:
                improved.append(pred_name)

        if len(improved) == 0:
            return "None"
        elif len(improved) == 1:
            return f"{improved[0]} Only"
        elif len(improved) == 2:
            return f"{improved[0]} + {improved[1]}"
        elif len(improved) == 3:
            return f"{improved[0]} + {improved[1]} + {improved[2]}"
        else:  # All 4
            return "All Four"

    df['category'] = df.apply(get_category, axis=1)

    # Extract workload category from trace name
    df['workload_category'] = df['Run'].str.split('_').str[0]

    return df

def print_summary(df):
    """Print comprehensive summary statistics."""
    print("\n" + "="*80)
    print("PREDICTOR COMPLEMENTARITY ANALYSIS SUMMARY")
    print("="*80)

    # Category distribution
    print("\nTRACE DISTRIBUTION BY PREDICTOR COMBINATION:")
    print("-" * 80)
    category_counts = df['category'].value_counts().sort_index()

    for category, count in category_counts.items():
        pct = (count / len(df)) * 100
        print(f"\n{category}: {count} traces ({pct:.1f}%)")
        traces = df[df['category'] == category]['Run'].tolist()
        # Print traces in columns for readability
        for i in range(0, len(traces), 3):
            print("  ", ", ".join(traces[i:i+3]))

    # Overall statistics
    print("\n" + "="*80)
    print("KEY INSIGHTS:")
    print("-" * 80)

    # How many traces no predictor improved
    no_improvement = len(df[df['category'] == 'None'])
    print(f"• Traces no predictor improved: {no_improvement} ({no_improvement/len(df)*100:.1f}%)")

    # Which predictor has most unique improvements
    unique_counts = {}
    for pred_name in PREDICTOR_NAMES:
        unique_counts[pred_name] = len(df[df['category'] == f'{pred_name} Only'])

    best_unique = max(unique_counts.items(), key=lambda x: x[1])
    print(f"• Predictor with most unique improvements: {best_unique[0]} ({best_unique[1]} traces)")

    # How many traces all four improved
    all_four = len(df[df['category'] == 'All Four'])
    print(f"• Traces all four predictors improved: {all_four} ({all_four/len(df)*100:.1f}%)")

    # Individual predictor coverage
    print("\n• Individual predictor improvement coverage:")
    for pred_name in PREDICTOR_NAMES:
        count = df[f'{pred_name}_improved'].sum()
        pct = (count / len(df)) * 100
        avg_improvement = df[df[f'{pred_name}_improved']][f'{pred_name}_improvement'].mean()
        print(f"  - {pred_name}: {count} traces ({pct:.1f}%), avg improvement: {avg_improvement:.3f} MPKI")

    # Pairwise overlaps
    print("\n• Pairwise predictor overlaps:")
    for i, pred1 in enumerate(PREDICTOR_NAMES):
        for pred2 in PREDICTOR_NAMES[i+1:]:
            overlap = (df[f'{pred1}_improved'] & df[f'{pred2}_improved']).sum()
            only_pred1 = (df[f'{pred1}_improved'] & ~df[f'{pred2}_improved']).sum()
            only_pred2 = (~df[f'{pred1}_improved'] & df[f'{pred2}_improved']).sum()
            print(f"  - {pred1} & {pred2}: {overlap} shared, {only_pred1} {pred1}-only, {only_pred2} {pred2}-only")

    print("="*80 + "\n")

def create_upset_plot(df, output_dir):
    """Create UpSet plot showing predictor intersections."""
    print("Creating UpSet plot...")

    if UPSETPLOT_AVAILABLE:
        # Prepare data for UpSet plot
        upset_data = df.set_index('Run')[[f'{pred}_improved' for pred in PREDICTOR_NAMES]]
        upset_data.columns = PREDICTOR_NAMES

        # Create UpSet plot
        fig = plt.figure(figsize=(14, 8))
        upset = UpSet(from_indicators(PREDICTOR_NAMES, data=upset_data),
                      subset_size='count',
                      show_counts=True,
                      sort_by='cardinality',
                      sort_categories_by='cardinality')
        upset.plot(fig=fig)
        plt.suptitle('Predictor Complementarity: UpSet Plot\n' +
                    f'Meaningful Improvement Threshold: {IMPROVEMENT_THRESHOLD} MPKI',
                    fontsize=14, y=0.98)

        plt.tight_layout()
        plt.savefig(output_dir / '01_upset_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_dir / '01_upset_plot.png'}")
    else:
        # Fallback: grouped bar chart
        category_counts = df['category'].value_counts().sort_values(ascending=False)

        fig, ax = plt.subplots(figsize=(14, 8))
        bars = ax.bar(range(len(category_counts)), category_counts.values)

        # Color bars by number of predictors involved
        colors = []
        for cat in category_counts.index:
            if cat == 'None':
                colors.append('#cccccc')
            elif 'Only' in cat:
                colors.append('#4CAF50')
            elif cat == 'All Four':
                colors.append('#F44336')
            elif '+' in cat:
                num_predictors = cat.count('+') + 1
                if num_predictors == 2:
                    colors.append('#2196F3')
                else:  # 3 predictors
                    colors.append('#FF9800')
            else:
                colors.append('#9C27B0')

        for bar, color in zip(bars, colors):
            bar.set_color(color)

        ax.set_xticks(range(len(category_counts)))
        ax.set_xticklabels(category_counts.index, rotation=45, ha='right')
        ax.set_ylabel('Number of Traces', fontsize=12)
        ax.set_title(f'Predictor Combination Distribution\nThreshold: {IMPROVEMENT_THRESHOLD} MPKI',
                    fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)

        # Add count labels on bars
        for i, (bar, count) in enumerate(zip(bars, category_counts.values)):
            ax.text(i, count + 0.5, str(count), ha='center', va='bottom', fontsize=10)

        plt.tight_layout()
        plt.savefig(output_dir / '01_category_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_dir / '01_category_distribution.png'}")

def create_heatmap(df, output_dir):
    """Create heatmap of improvements sorted by baseline MPKI."""
    print("Creating improvement heatmap...")

    # Sort by baseline MPKI descending
    df_sorted = df.sort_values('baseline_MPKI', ascending=False).copy()

    # Prepare improvement matrix
    improvement_matrix = df_sorted[[f'{pred}_improvement' for pred in PREDICTOR_NAMES]].values
    trace_names = df_sorted['Run'].values

    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(12, len(df_sorted) * 0.15)))

    # Create heatmap with diverging colormap
    im = ax.imshow(improvement_matrix, aspect='auto', cmap='RdBu_r',
                   vmin=-2, vmax=max(2, improvement_matrix.max()))

    # Set ticks
    ax.set_xticks(range(len(PREDICTOR_NAMES)))
    ax.set_xticklabels(PREDICTOR_NAMES, fontsize=11, fontweight='bold')
    ax.set_yticks(range(len(trace_names)))
    ax.set_yticklabels(trace_names, fontsize=7)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('MPKI Improvement over Baseline', rotation=270, labelpad=20, fontsize=11)

    # Add title
    ax.set_title('Predictor Improvement Heatmap\n(Traces sorted by baseline MPKI, descending)',
                fontsize=13, fontweight='bold', pad=15)

    # Add grid
    ax.set_xticks(np.arange(len(PREDICTOR_NAMES)) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(trace_names)) - 0.5, minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / '02_improvement_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_dir / '02_improvement_heatmap.png'}")

def create_workload_category_plot(df, output_dir):
    """Create stacked bar chart by workload category."""
    print("Creating workload category analysis plot...")

    # Get unique workload categories
    workload_categories = sorted(df['workload_category'].unique())

    # For each workload category, count predictor combinations
    category_data = {}
    for wl in workload_categories:
        wl_df = df[df['workload_category'] == wl]
        category_counts = wl_df['category'].value_counts()
        category_data[wl] = category_counts

    # Get all unique predictor combinations across all workloads
    all_combinations = set()
    for counts in category_data.values():
        all_combinations.update(counts.index)
    all_combinations = sorted(list(all_combinations))

    # Prepare data for stacked bar chart (proportions)
    proportions = []
    for combo in all_combinations:
        combo_props = []
        for wl in workload_categories:
            wl_df = df[df['workload_category'] == wl]
            count = category_data[wl].get(combo, 0)
            prop = count / len(wl_df) if len(wl_df) > 0 else 0
            combo_props.append(prop)
        proportions.append(combo_props)

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(14, 8))

    # Color scheme
    color_map = {
        'None': '#cccccc',
        'All Four': '#F44336',
    }
    # Assign colors to single predictor
    single_colors = ['#4CAF50', '#8BC34A', '#CDDC39', '#FFEB3B']
    for i, pred in enumerate(PREDICTOR_NAMES):
        color_map[f'{pred} Only'] = single_colors[i]

    # Pairs get blue shades
    pair_colors = ['#2196F3', '#03A9F4', '#00BCD4', '#009688', '#00796B', '#004D40']
    pair_idx = 0

    # Triples get orange shades
    triple_colors = ['#FF9800', '#FF5722', '#F44336', '#E91E63']
    triple_idx = 0

    for combo in all_combinations:
        if combo not in color_map:
            if 'Only' in combo:
                continue  # Already handled
            elif combo.count('+') == 1:  # Pair
                color_map[combo] = pair_colors[pair_idx % len(pair_colors)]
                pair_idx += 1
            elif combo.count('+') == 2:  # Triple
                color_map[combo] = triple_colors[triple_idx % len(triple_colors)]
                triple_idx += 1
            else:
                color_map[combo] = '#9C27B0'

    # Plot stacked bars
    bottom = np.zeros(len(workload_categories))
    for i, combo in enumerate(all_combinations):
        color = color_map.get(combo, '#9C27B0')
        ax.bar(workload_categories, proportions[i], bottom=bottom,
               label=combo, color=color, edgecolor='white', linewidth=0.5)
        bottom += proportions[i]

    ax.set_ylabel('Proportion of Traces', fontsize=12, fontweight='bold')
    ax.set_xlabel('Workload Category', fontsize=12, fontweight='bold')
    ax.set_title('Predictor Complementarity by Workload Category\n' +
                f'(Threshold: {IMPROVEMENT_THRESHOLD} MPKI)',
                fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.grid(axis='y', alpha=0.3)

    # Legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9,
             framealpha=0.9, title='Predictor Combination')

    # Add trace counts on x-axis labels
    new_labels = [f"{wl}\n(n={len(df[df['workload_category']==wl])})"
                  for wl in workload_categories]
    ax.set_xticklabels(new_labels, fontsize=10)

    plt.tight_layout()
    plt.savefig(output_dir / '03_workload_category_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_dir / '03_workload_category_analysis.png'}")

def save_results(df, output_dir):
    """Save categorized results to CSV."""
    print("\nSaving results...")

    # Prepare output dataframe
    output_cols = ['Run', 'workload_category', 'baseline_MPKI']

    # Add predictor MPKIs
    for pred in PREDICTOR_NAMES:
        output_cols.append(f'{pred}_MPKI')

    # Add improvements
    for pred in PREDICTOR_NAMES:
        output_cols.append(f'{pred}_improvement')

    # Add improved flags
    for pred in PREDICTOR_NAMES:
        output_cols.append(f'{pred}_improved')

    # Add category
    output_cols.append('category')

    output_df = df[output_cols].copy()

    # Sort by baseline MPKI descending
    output_df = output_df.sort_values('baseline_MPKI', ascending=False)

    output_path = output_dir / 'complementarity_analysis.csv'
    output_df.to_csv(output_path, index=False)
    print(f"  Saved: {output_path}")

    return output_df

def main():
    """Main analysis pipeline."""
    print("\n" + "="*80)
    print("PREDICTOR COMPLEMENTARITY ANALYSIS")
    print("="*80)
    print(f"Improvement threshold: {IMPROVEMENT_THRESHOLD} MPKI")
    print(f"Predictors: {', '.join(PREDICTOR_NAMES)}")
    print("="*80 + "\n")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load and process data
    df = load_predictor_data()
    df = compute_improvements(df)
    df = categorize_traces(df)

    # Print summary
    print_summary(df)

    # Create visualizations
    create_upset_plot(df, OUTPUT_DIR)
    create_heatmap(df, OUTPUT_DIR)
    create_workload_category_plot(df, OUTPUT_DIR)

    # Save results
    save_results(df, OUTPUT_DIR)

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"Output directory: {OUTPUT_DIR}")
    print("\nGenerated files:")
    print("  • complementarity_analysis.csv - Full categorized data")
    if UPSETPLOT_AVAILABLE:
        print("  • 01_upset_plot.png - UpSet plot of predictor intersections")
    else:
        print("  • 01_category_distribution.png - Category distribution bar chart")
    print("  • 02_improvement_heatmap.png - Improvement heatmap sorted by baseline")
    print("  • 03_workload_category_analysis.png - Workload category analysis")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
