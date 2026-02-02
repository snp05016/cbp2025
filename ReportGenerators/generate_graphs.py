import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def load_results(csv_path):
    df = pd.read_csv(csv_path)
    if 'MR' in df.columns:
        df['MR_numeric'] = df['MR'].str.rstrip('%').astype(float)
    if '50PercMR' in df.columns:
        df['50PercMR_numeric'] = df['50PercMR'].str.rstrip('%').astype(float)
    return df

def create_misprediction_rate_graphs(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    categories = df['Workload'].unique()
    category_mr = df.groupby('Workload')['MR_numeric'].mean().sort_values(ascending=False)
    
    plt.figure(figsize=(10, 6))
    plt.bar(category_mr.index, category_mr.values, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'])
    plt.xlabel('Benchmark Category')
    plt.ylabel('Average Misprediction Rate (%)')
    plt.title('Average Misprediction Rate by Benchmark Category')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/category_avg_misprediction_rate.png', dpi=150)
    plt.close()
    
    for category in categories:
        cat_data = df[df['Workload'] == category].sort_values('MR_numeric', ascending=False)
        plt.figure(figsize=(14, 6))
        plt.bar(range(len(cat_data)), cat_data['MR_numeric'], color='steelblue')
        plt.xlabel('Trace')
        plt.ylabel('Misprediction Rate (%)')
        plt.title(f'Misprediction Rate for {category} Traces')
        plt.xticks(range(len(cat_data)), cat_data['Run'], rotation=90, fontsize=8)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{category}_misprediction_rates.png', dpi=150)
        plt.close()

def create_performance_metric_graphs(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    categories = df['Workload'].unique()
    
    for category in categories:
        cat_data = df[df['Workload'] == category]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        axes[0, 0].bar(range(len(cat_data)), cat_data['IPC'], color='#2ca02c')
        axes[0, 0].set_xlabel('Trace Index')
        axes[0, 0].set_ylabel('IPC')
        axes[0, 0].set_title(f'{category}: Instructions Per Cycle')
        axes[0, 0].tick_params(axis='x', labelsize=8)
        
        axes[0, 1].bar(range(len(cat_data)), cat_data['MPKI'], color='#d62728')
        axes[0, 1].set_xlabel('Trace Index')
        axes[0, 1].set_ylabel('MPKI')
        axes[0, 1].set_title(f'{category}: Mispredictions Per Kilo Instructions')
        axes[0, 1].tick_params(axis='x', labelsize=8)
        
        axes[1, 0].bar(range(len(cat_data)), cat_data['MR_numeric'], color='#ff7f0e')
        axes[1, 0].set_xlabel('Trace Index')
        axes[1, 0].set_ylabel('Miss Rate (%)')
        axes[1, 0].set_title(f'{category}: Branch Misprediction Rate')
        axes[1, 0].tick_params(axis='x', labelsize=8)
        
        axes[1, 1].scatter(cat_data['IPC'], cat_data['MPKI'], alpha=0.6, s=100, color='#9467bd')
        axes[1, 1].set_xlabel('IPC')
        axes[1, 1].set_ylabel('MPKI')
        axes[1, 1].set_title(f'{category}: IPC vs MPKI')
        axes[1, 1].grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{category}_performance_metrics.png', dpi=150)
        plt.close()

def create_difficulty_analysis_graphs(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    categories = df['Workload'].unique()
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    for idx, category in enumerate(categories):
        cat_data = df[df['Workload'] == category]
        axes[idx].scatter(cat_data['MR_numeric'], cat_data['MPKI'], alpha=0.7, s=80)
        axes[idx].set_xlabel('Misprediction Rate (%)')
        axes[idx].set_ylabel('MPKI')
        axes[idx].set_title(f'{category}: MR vs MPKI')
        axes[idx].grid(alpha=0.3)
    
    for idx in range(len(categories), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/difficulty_mr_vs_mpki.png', dpi=150)
    plt.close()
    
    category_stats = df.groupby('Workload').agg({
        'MR_numeric': ['mean', 'std', 'max', 'min'],
        'MPKI': ['mean', 'std', 'max', 'min']
    }).round(3)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    mr_means = category_stats['MR_numeric']['mean'].sort_values(ascending=False)
    mr_stds = category_stats['MR_numeric']['std'].loc[mr_means.index]
    axes[0].bar(range(len(mr_means)), mr_means, yerr=mr_stds, capsize=5, color='steelblue', alpha=0.7)
    axes[0].set_xticks(range(len(mr_means)))
    axes[0].set_xticklabels(mr_means.index, rotation=45, ha='right')
    axes[0].set_ylabel('Misprediction Rate (%) Mean ± Std')
    axes[0].set_title('MR Variability by Category')
    axes[0].grid(axis='y', alpha=0.3)
    
    mpki_means = category_stats['MPKI']['mean'].sort_values(ascending=False)
    mpki_stds = category_stats['MPKI']['std'].loc[mpki_means.index]
    axes[1].bar(range(len(mpki_means)), mpki_means, yerr=mpki_stds, capsize=5, color='coral', alpha=0.7)
    axes[1].set_xticks(range(len(mpki_means)))
    axes[1].set_xticklabels(mpki_means.index, rotation=45, ha='right')
    axes[1].set_ylabel('MPKI Mean ± Std')
    axes[1].set_title('MPKI Variability by Category')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/difficulty_variability.png', dpi=150)
    plt.close()
    
    fig, ax = plt.subplots(figsize=(12, 8))
    for category in categories:
        cat_data = df[df['Workload'] == category]
        ax.scatter(cat_data['IPC'], cat_data['MR_numeric'], label=category, alpha=0.6, s=60)
    ax.set_xlabel('IPC')
    ax.set_ylabel('Misprediction Rate (%)')
    ax.set_title('IPC vs MR Across All Categories')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/difficulty_ipc_vs_mr_all.png', dpi=150)
    plt.close()

def create_branch_prediction_analysis(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    categories = df['Workload'].unique()
    
    for category in categories:
        cat_data = df[df['Workload'] == category]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        axes[0, 0].bar(range(len(cat_data)), cat_data['NumBr'], color='#1f77b4')
        axes[0, 0].set_xlabel('Trace Index')
        axes[0, 0].set_ylabel('Number of Branches')
        axes[0, 0].set_title(f'{category}: Total Branches')
        axes[0, 0].tick_params(axis='x', labelsize=8)
        
        axes[0, 1].bar(range(len(cat_data)), cat_data['MispBr'], color='#ff7f0e')
        axes[0, 1].set_xlabel('Trace Index')
        axes[0, 1].set_ylabel('Mispredicted Branches')
        axes[0, 1].set_title(f'{category}: Mispredicted Branches')
        axes[0, 1].tick_params(axis='x', labelsize=8)
        
        axes[1, 0].bar(range(len(cat_data)), cat_data['BrPerCyc'], color='#2ca02c')
        axes[1, 0].set_xlabel('Trace Index')
        axes[1, 0].set_ylabel('Branches Per Cycle')
        axes[1, 0].set_title(f'{category}: Branch Density')
        axes[1, 0].tick_params(axis='x', labelsize=8)
        
        axes[1, 1].scatter(cat_data['NumBr'], cat_data['MispBr'], alpha=0.6, s=80, color='#d62728')
        axes[1, 1].set_xlabel('Total Branches')
        axes[1, 1].set_ylabel('Mispredicted Branches')
        axes[1, 1].set_title(f'{category}: Branches vs Mispredictions')
        axes[1, 1].grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{category}_branch_analysis.png', dpi=150)
        plt.close()

def create_cycles_analysis(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    categories = df['Workload'].unique()
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    for idx, category in enumerate(categories):
        cat_data = df[df['Workload'] == category]
        axes[idx].scatter(cat_data['Cycles'], cat_data['Instr'], alpha=0.6, s=60)
        axes[idx].set_xlabel('Cycles')
        axes[idx].set_ylabel('Instructions')
        axes[idx].set_title(f'{category}: Cycles vs Instructions')
        axes[idx].grid(alpha=0.3)
    
    for idx in range(len(categories), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/cycles_vs_instructions_all.png', dpi=150)
    plt.close()
    
    for category in categories:
        cat_data = df[df['Workload'] == category]
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        axes[0].bar(range(len(cat_data)), cat_data['CycWPAvg'], color='purple', alpha=0.7)
        axes[0].set_xlabel('Trace Index')
        axes[0].set_ylabel('Avg Cycles Wasted Per Misprediction')
        axes[0].set_title(f'{category}: CycWPAvg')
        axes[0].tick_params(axis='x', labelsize=8)
        
        axes[1].scatter(cat_data['MR_numeric'], cat_data['CycWPAvg'], alpha=0.6, s=80, color='purple')
        axes[1].set_xlabel('Misprediction Rate (%)')
        axes[1].set_ylabel('CycWPAvg')
        axes[1].set_title(f'{category}: MR vs CycWPAvg')
        axes[1].grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{category}_cycle_waste.png', dpi=150)
        plt.close()
def create_benchmark_analysis(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract benchmark category from Run column (e.g., 'int_1' -> 'int')
    df['BenchCategory'] = df['Run'].str.extract(r'^([a-z]+)_')[0]
    df['BenchCategory'] = df['BenchCategory'].fillna('other')
    
    # Focus on specific categories
    categories_of_interest = ['int', 'infra', 'compress', 'web', 'fp']
    colors_map = {
        'int': '#1f77b4',
        'infra': '#ff7f0e', 
        'compress': '#2ca02c',
        'web': '#d62728',
        'fp': '#9467bd'
    }
    
    # Filter to only include categories of interest
    df_filtered = df[df['BenchCategory'].isin(categories_of_interest)]
    
    # 1. MPKI Comparison - Main Bar Chart with error bars
    mpki_stats = df_filtered.groupby('BenchCategory')['MPKI'].agg(['mean', 'std', 'min', 'max']).reset_index()
    mpki_stats = mpki_stats.sort_values('mean', ascending=False)
    
    plt.figure(figsize=(12, 7))
    x = np.arange(len(mpki_stats))
    bars = plt.bar(x, mpki_stats['mean'], 
                   yerr=mpki_stats['std'],
                   capsize=8,
                   color=[colors_map[cat] for cat in mpki_stats['BenchCategory']],
                   alpha=0.8,
                   edgecolor='black',
                   linewidth=1.5)
    
    plt.xlabel('Benchmark Category', fontsize=12, weight='bold')
    plt.ylabel('MPKI (Mispredictions Per Kilo Instructions)', fontsize=12, weight='bold')
    plt.title('MPKI Comparison Across Benchmark Categories', fontsize=14, weight='bold', pad=20)
    plt.xticks(x, [cat.upper() for cat in mpki_stats['BenchCategory']], fontsize=11)
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add value labels on bars
    for i, (bar, val, std) in enumerate(zip(bars, mpki_stats['mean'], mpki_stats['std'])):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + std + 0.1,
                f'{val:.2f}±{std:.2f}',
                ha='center', va='bottom', fontsize=9, weight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/mpki_comparison_benchmarks.png', dpi=150)
    plt.close()
    
    # 2. MPKI Box Plot
    plt.figure(figsize=(12, 7))
    data_to_plot = [df_filtered[df_filtered['BenchCategory'] == cat]['MPKI'].dropna() 
                    for cat in categories_of_interest if cat in df_filtered['BenchCategory'].values]
    labels = [cat.upper() for cat in categories_of_interest if cat in df_filtered['BenchCategory'].values]
    
    bp = plt.boxplot(data_to_plot, labels=labels, patch_artist=True, showmeans=True)
    
    for patch, cat in zip(bp['boxes'], [c for c in categories_of_interest if c in df_filtered['BenchCategory'].values]):
        patch.set_facecolor(colors_map[cat])
        patch.set_alpha(0.7)
    
    plt.ylabel('MPKI', fontsize=12, weight='bold')
    plt.title('MPKI Distribution Across Benchmark Categories', fontsize=14, weight='bold', pad=20)
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/mpki_boxplot_benchmarks.png', dpi=150)
    plt.close()
    
    # 3. MPKI per individual run within each category - Representative Sample
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    axes = axes.flatten()
    
    for idx, category in enumerate(categories_of_interest):
        if category not in df_filtered['BenchCategory'].values:
            axes[idx].axis('off')
            continue
            
        cat_data = df_filtered[df_filtered['BenchCategory'] == category].sort_values('MPKI', ascending=False)
        
        # Select representative sample: hardest 2, median 2, easiest 2
        n = len(cat_data)
        if n <= 6:
            # If 6 or fewer benchmarks, show all
            sample_data = cat_data
        else:
            # Select representative benchmarks
            hardest_2 = cat_data.head(2)
            easiest_2 = cat_data.tail(2)
            # Get 2 from middle
            mid_start = (n - 2) // 2
            middle_2 = cat_data.iloc[mid_start:mid_start+2]
            sample_data = pd.concat([hardest_2, middle_2, easiest_2])
        
        # Calculate geomean for ALL benchmarks in category
        def calculate_geomean(values):
            return np.exp(np.mean(np.log(values)))
        
        geomean_mpki = calculate_geomean(cat_data['MPKI'])
        
        axes[idx].bar(range(len(sample_data)), sample_data['MPKI'], 
                     color=colors_map[category], alpha=0.7, edgecolor='black')
        axes[idx].set_xlabel('Representative Traces', fontsize=10)
        axes[idx].set_ylabel('MPKI', fontsize=10)
        axes[idx].set_title(f'{category.upper()} - MPKI (n={n}, showing {len(sample_data)} representative)', 
                           fontsize=10, weight='bold')
        axes[idx].tick_params(axis='x', labelsize=7)
        axes[idx].grid(axis='y', alpha=0.3)
        
        # Add mean and geomean lines
        mean_val = cat_data['MPKI'].mean()
        axes[idx].axhline(y=mean_val, color='red', linestyle='--', linewidth=1.5, 
                         label=f'Mean: {mean_val:.2f}')
        axes[idx].axhline(y=geomean_mpki, color='darkgreen', linestyle=':', linewidth=2, 
                         label=f'Geomean: {geomean_mpki:.2f}')
        axes[idx].legend(fontsize=8)
        
        # Add labels to show which samples these are
        labels = []
        for i, (_, row) in enumerate(sample_data.iterrows()):
            if i < 2:
                labels.append(f"H{i+1}")  # Hardest
            elif i < 4:
                labels.append(f"M{i-1}")  # Middle
            else:
                labels.append(f"E{i-3}")  # Easiest
        axes[idx].set_xticks(range(len(sample_data)))
        axes[idx].set_xticklabels(labels, fontsize=8)
    
    axes[-1].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/mpki_per_category_detailed.png', dpi=150)
    plt.close()
    
    # 4. Summary Statistics Table with Geometric Mean
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.axis('tight')
    ax.axis('off')
    
    table_data = [['Category', 'Mean MPKI', 'Geomean MPKI', 'Std MPKI', 'Min MPKI', 'Max MPKI', 'Count']]
    
    def calculate_geomean(values):
        return np.exp(np.mean(np.log(values)))
    
    for category in categories_of_interest:
        if category in mpki_stats['BenchCategory'].values:
            row_data = mpki_stats[mpki_stats['BenchCategory'] == category].iloc[0]
            count = len(df_filtered[df_filtered['BenchCategory'] == category])
            
            # Calculate geomean from all data
            cat_mpki_values = df_filtered[df_filtered['BenchCategory'] == category]['MPKI']
            geomean = calculate_geomean(cat_mpki_values)
            
            table_data.append([
                category.upper(),
                f"{row_data['mean']:.3f}",
                f"{geomean:.3f}",
                f"{row_data['std']:.3f}",
                f"{row_data['min']:.3f}",
                f"{row_data['max']:.3f}",
                str(count)
            ])
    
    table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                    colWidths=[0.12, 0.15, 0.15, 0.15, 0.14, 0.14, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2.5)
    
    # Color header row
    for i in range(7):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Color data rows
    for i in range(1, len(table_data)):
        for j in range(7):
            table[(i, j)].set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
    
    plt.title('MPKI Statistics by Benchmark Category (Full Dataset)', pad=20, fontsize=14, weight='bold')
    plt.savefig(f'{output_dir}/mpki_summary_table.png', dpi=150, bbox_inches='tight')
    plt.close()
        
def create_first_vs_second_half_analysis(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    categories = df['Workload'].unique()
    
    for category in categories:
        cat_data = df[df['Workload'] == category]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        x = np.arange(len(cat_data))
        width = 0.35
        
        axes[0, 0].bar(x - width/2, cat_data['IPC'], width, label='Full', alpha=0.8)
        axes[0, 0].bar(x + width/2, cat_data['50PercIPC'], width, label='First 50%', alpha=0.8)
        axes[0, 0].set_xlabel('Trace Index')
        axes[0, 0].set_ylabel('IPC')
        axes[0, 0].set_title(f'{category}: IPC Comparison')
        axes[0, 0].legend()
        axes[0, 0].tick_params(axis='x', labelsize=8)
        
        axes[0, 1].bar(x - width/2, cat_data['MR_numeric'], width, label='Full', alpha=0.8)
        axes[0, 1].bar(x + width/2, cat_data['50PercMR_numeric'], width, label='First 50%', alpha=0.8)
        axes[0, 1].set_xlabel('Trace Index')
        axes[0, 1].set_ylabel('Misprediction Rate (%)')
        axes[0, 1].set_title(f'{category}: MR Comparison')
        axes[0, 1].legend()
        axes[0, 1].tick_params(axis='x', labelsize=8)
        
        axes[1, 0].bar(x - width/2, cat_data['MPKI'], width, label='Full', alpha=0.8)
        axes[1, 0].bar(x + width/2, cat_data['50PercMPKI'], width, label='First 50%', alpha=0.8)
        axes[1, 0].set_xlabel('Trace Index')
        axes[1, 0].set_ylabel('MPKI')
        axes[1, 0].set_title(f'{category}: MPKI Comparison')
        axes[1, 0].legend()
        axes[1, 0].tick_params(axis='x', labelsize=8)
        
        mr_diff = cat_data['50PercMR_numeric'] - cat_data['MR_numeric']
        axes[1, 1].bar(range(len(mr_diff)), mr_diff, color='coral', alpha=0.7)
        axes[1, 1].axhline(y=0, color='black', linestyle='--', linewidth=0.8)
        axes[1, 1].set_xlabel('Trace Index')
        axes[1, 1].set_ylabel('MR Difference (First50% - Full)')
        axes[1, 1].set_title(f'{category}: Warmup Effect on MR')
        axes[1, 1].tick_params(axis='x', labelsize=8)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{category}_first_vs_second_half.png', dpi=150)
        plt.close()

def create_aggregate_comparison_graphs(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    category_stats = df.groupby('Workload').agg({
        'IPC': 'mean',
        'MR_numeric': 'mean',
        'MPKI': 'mean',
        'BrPerCyc': 'mean',
        'CycWPAvg': 'mean'
    }).round(3)
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    metrics = ['IPC', 'MR_numeric', 'MPKI', 'BrPerCyc', 'CycWPAvg']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for idx, metric in enumerate(metrics):
        row, col = idx // 3, idx % 3
        data = category_stats[metric].sort_values(ascending=False)
        axes[row, col].bar(range(len(data)), data.values, color=colors[idx], alpha=0.8)
        axes[row, col].set_xticks(range(len(data)))
        axes[row, col].set_xticklabels(data.index, rotation=45, ha='right')
        axes[row, col].set_ylabel(metric)
        axes[row, col].set_title(f'Average {metric} by Category')
        axes[row, col].grid(axis='y', alpha=0.3)
    
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/aggregate_category_comparison.png', dpi=150)
    plt.close()

def main():
    csv_path = 'results.csv'
    output_dir = 'analysis_graphs'
    
    print(f'Loading {csv_path}...')
    df = load_results(csv_path)
    
    print('Creating misprediction rate graphs...')
    create_misprediction_rate_graphs(df, output_dir)
    
    print('Creating performance metric graphs...')
    create_performance_metric_graphs(df, output_dir)
    
    print('Creating difficulty analysis graphs...')
    create_difficulty_analysis_graphs(df, output_dir)
    
    print('Creating branch prediction analysis...')
    create_branch_prediction_analysis(df, output_dir)
    
    print('Creating benchmark type comparison (INT vs FP vs OTHER)...')
    create_benchmark_analysis(df, output_dir)
    
    print('Creating cycles analysis...')
    create_cycles_analysis(df, output_dir)
    
    print('Creating first vs second half analysis...')
    create_first_vs_second_half_analysis(df, output_dir)
    
    print('Creating aggregate comparison graphs...')
    create_aggregate_comparison_graphs(df, output_dir)
    
    print(f'\nAll graphs saved to {output_dir}/')

if __name__ == '__main__':
    main()
