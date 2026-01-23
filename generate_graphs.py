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
    
    print('Creating cycles analysis...')
    create_cycles_analysis(df, output_dir)
    
    print('Creating first vs second half analysis...')
    create_first_vs_second_half_analysis(df, output_dir)
    
    print('Creating aggregate comparison graphs...')
    create_aggregate_comparison_graphs(df, output_dir)
    
    print(f'\nAll graphs saved to {output_dir}/')

if __name__ == '__main__':
    main()
