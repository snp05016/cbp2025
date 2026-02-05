#!/usr/bin/env python3
"""
Comprehensive Predictor Comparison Script
Analyzes and compares 5 different branch predictors across multiple metrics
Generates text reports and visualizations without popups
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

# Define predictor files and names
PREDICTORS = {
    'Baseline': 'baseline.csv',
    'Programming Idiom Predictor-Mose': 'programming-idiom-predictor-mose.csv',
    'TAGE-SC-L-Alberto-Ros': 'tage-sc-l-alberto-ros-results.csv',
    'TAGE-SCL-Andrez-Seznec': 'tage-scl-andrez-seznec-results.csv',
    'Load-Value-Correlator-Man': 'load-value-correlator-man-results.csv',
    'Register-Value-Aware-Toru': 'register-value-aware-toru-results.csv'
}

# Key metrics to analyze
METRICS = {
    'MPKI': 'Mispredictions Per Kilo Instructions',
    'MR': 'Miss Rate (%)',
    'IPC': 'Instructions Per Cycle',
    'Cycles': 'Total Cycles',
    'MispBr': 'Mispredicted Branches',
    'NumBr': 'Total Branches',
    'BrPerCyc': 'Branches Per Cycle',
    'MispBrPerCyc': 'Mispredictions Per Cycle',
    'CycWPPKI': 'Cycles Wasted Per Kilo Instructions'
}

class PredictorComparator:
    def __init__(self, base_dir='.', output_dir='reports/comparison-predictors'):
        self.base_dir = Path(base_dir)
        self.output_dir = Path(output_dir)
        self.data = {}
        self.load_data()
        self.create_output_structure()
        
    def load_data(self):
        """Load all predictor CSV files"""
        print("Loading predictor data...")
        for name, filename in PREDICTORS.items():
            filepath = self.base_dir / filename
            if filepath.exists():
                df = pd.read_csv(filepath)
                # Clean percentage values
                if 'MR' in df.columns:
                    df['MR'] = df['MR'].str.rstrip('%').astype('float')
                if '50PercMR' in df.columns:
                    df['50PercMR'] = df['50PercMR'].str.rstrip('%').astype('float')
                self.data[name] = df
                print(f"  ✓ Loaded {name}: {len(df)} traces")
            else:
                print(f"  ✗ File not found: {filename}")
                
    def create_output_structure(self):
        """Create directory structure for outputs"""
        # Main output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Overall comparisons
        (self.output_dir / 'overall').mkdir(exist_ok=True)
        (self.output_dir / 'overall' / 'figures').mkdir(exist_ok=True)
        (self.output_dir / 'overall' / 'text_reports').mkdir(exist_ok=True)
        
        # Pairwise comparisons
        predictor_names = list(self.data.keys())
        for i, pred1 in enumerate(predictor_names):
            for pred2 in predictor_names[i+1:]:
                pair_dir = self.output_dir / f'{pred1}_vs_{pred2}'
                pair_dir.mkdir(exist_ok=True)
                (pair_dir / 'figures').mkdir(exist_ok=True)
                (pair_dir / 'text_reports').mkdir(exist_ok=True)
                
    def generate_overall_comparison(self):
        """Generate overall comparison across all predictors"""
        print("\nGenerating overall comparison reports...")
        
        # 1. Summary Statistics Report
        self._generate_summary_statistics()
        
        # 2. Misprediction Analysis
        self._generate_misprediction_analysis()
        
        # 3. Performance Metrics Comparison
        self._generate_performance_comparison()
        
        # 4. Trace Category Analysis
        self._generate_category_analysis()
        
        # 5. Overall Rankings
        self._generate_rankings()
        
    def _generate_summary_statistics(self):
        """Generate summary statistics for all predictors"""
        output_file = self.output_dir / 'overall' / 'text_reports' / '01_summary_statistics.txt'
        
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("OVERALL SUMMARY STATISTICS - ALL PREDICTORS\n")
            f.write("=" * 80 + "\n\n")
            
            for pred_name, df in self.data.items():
                f.write(f"\n{pred_name}\n")
                f.write("-" * 80 + "\n")
                
                # Key statistics
                stats = {
                    'Total Traces': len(df),
                    'Avg MPKI': df['MPKI'].mean(),
                    'Median MPKI': df['MPKI'].median(),
                    'Std MPKI': df['MPKI'].std(),
                    'Avg Miss Rate (%)': df['MR'].mean(),
                    'Avg IPC': df['IPC'].mean(),
                    'Total Branches': df['NumBr'].sum(),
                    'Total Mispredictions': df['MispBr'].sum(),
                    'Avg Cycles': df['Cycles'].mean(),
                }
                
                for key, value in stats.items():
                    if 'Total' in key:
                        f.write(f"  {key:30s}: {value:>15,.0f}\n")
                    else:
                        f.write(f"  {key:30s}: {value:>15.4f}\n")
                        
        print(f"  ✓ Summary statistics saved to {output_file}")
        
        # Generate comparison plot
        self._plot_summary_comparison()
        
    def _plot_summary_comparison(self):
        """Plot summary comparison of key metrics"""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Overall Predictor Comparison - Key Metrics', fontsize=16, fontweight='bold')
        
        metrics_to_plot = ['MPKI', 'MR', 'IPC', 'Cycles', 'BrPerCyc', 'CycWPPKI']
        
        for idx, metric in enumerate(metrics_to_plot):
            ax = axes[idx // 3, idx % 3]
            
            # Prepare data
            data_to_plot = []
            labels = []
            for pred_name, df in self.data.items():
                if metric in df.columns:
                    data_to_plot.append(df[metric].values)
                    labels.append(pred_name)
            
            # Create box plot
            bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
            
            # Color the boxes
            colors = plt.cm.Set3(np.linspace(0, 1, len(data_to_plot)))
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
            
            ax.set_title(f'{METRICS.get(metric, metric)}', fontweight='bold')
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
            ax.grid(True, alpha=0.3)
            
        plt.tight_layout()
        output_file = self.output_dir / 'overall' / 'figures' / '01_summary_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Summary comparison plot saved to {output_file}")
        
    def _generate_misprediction_analysis(self):
        """Detailed misprediction analysis"""
        output_file = self.output_dir / 'overall' / 'text_reports' / '02_misprediction_analysis.txt'
        
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("MISPREDICTION ANALYSIS - ALL PREDICTORS\n")
            f.write("=" * 80 + "\n\n")
            
            # Compare misprediction metrics
            f.write("Average Misprediction Metrics:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Predictor':<35} {'MPKI':<12} {'Miss Rate %':<12} {'Mispred/Cyc':<12}\n")
            f.write("-" * 80 + "\n")
            
            for pred_name, df in self.data.items():
                mpki = df['MPKI'].mean()
                mr = df['MR'].mean()
                mpc = df['MispBrPerCyc'].mean()
                f.write(f"{pred_name:<35} {mpki:<12.4f} {mr:<12.4f} {mpc:<12.6f}\n")
            
            f.write("\n\nBest and Worst Cases:\n")
            f.write("-" * 80 + "\n")
            
            for pred_name, df in self.data.items():
                f.write(f"\n{pred_name}:\n")
                
                # Best case (lowest MPKI)
                best_idx = df['MPKI'].idxmin()
                best_trace = df.loc[best_idx, 'Run']
                best_mpki = df.loc[best_idx, 'MPKI']
                
                # Worst case (highest MPKI)
                worst_idx = df['MPKI'].idxmax()
                worst_trace = df.loc[worst_idx, 'Run']
                worst_mpki = df.loc[worst_idx, 'MPKI']
                
                f.write(f"  Best:  {best_trace:<30} MPKI = {best_mpki:.4f}\n")
                f.write(f"  Worst: {worst_trace:<30} MPKI = {worst_mpki:.4f}\n")
                
        print(f"  ✓ Misprediction analysis saved to {output_file}")
        
        # Generate visualization
        self._plot_misprediction_comparison()
        
    def _plot_misprediction_comparison(self):
        """Plot misprediction comparison"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Misprediction Analysis - All Predictors', fontsize=16, fontweight='bold')
        
        # 1. MPKI Distribution
        ax = axes[0, 0]
        for pred_name, df in self.data.items():
            ax.hist(df['MPKI'], alpha=0.5, label=pred_name, bins=30)
        ax.set_xlabel('MPKI')
        ax.set_ylabel('Frequency')
        ax.set_title('MPKI Distribution')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # 2. Miss Rate Distribution
        ax = axes[0, 1]
        for pred_name, df in self.data.items():
            ax.hist(df['MR'], alpha=0.5, label=pred_name, bins=30)
        ax.set_xlabel('Miss Rate (%)')
        ax.set_ylabel('Frequency')
        ax.set_title('Miss Rate Distribution')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # 3. Average MPKI Comparison
        ax = axes[1, 0]
        pred_names = list(self.data.keys())
        avg_mpki = [self.data[name]['MPKI'].mean() for name in pred_names]
        colors = plt.cm.Set3(np.linspace(0, 1, len(pred_names)))
        bars = ax.bar(range(len(pred_names)), avg_mpki, color=colors)
        ax.set_xticks(range(len(pred_names)))
        ax.set_xticklabels(pred_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Average MPKI')
        ax.set_title('Average MPKI by Predictor')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for bar, val in zip(bars, avg_mpki):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.3f}', ha='center', va='bottom', fontsize=9)
        
        # 4. Mispredictions per Cycle
        ax = axes[1, 1]
        avg_mpc = [self.data[name]['MispBrPerCyc'].mean() for name in pred_names]
        bars = ax.bar(range(len(pred_names)), avg_mpc, color=colors)
        ax.set_xticks(range(len(pred_names)))
        ax.set_xticklabels(pred_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Mispredictions Per Cycle')
        ax.set_title('Average Mispredictions Per Cycle')
        ax.grid(True, alpha=0.3, axis='y')
        
        for bar, val in zip(bars, avg_mpc):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.5f}', ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        output_file = self.output_dir / 'overall' / 'figures' / '02_misprediction_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Misprediction comparison plot saved to {output_file}")
        
    def _generate_performance_comparison(self):
        """Generate performance metrics comparison"""
        output_file = self.output_dir / 'overall' / 'text_reports' / '03_performance_metrics.txt'
        
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("PERFORMANCE METRICS COMPARISON\n")
            f.write("=" * 80 + "\n\n")
            
            # IPC Analysis
            f.write("Instructions Per Cycle (IPC) Analysis:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Predictor':<35} {'Avg IPC':<12} {'Median':<12} {'Std Dev':<12}\n")
            f.write("-" * 80 + "\n")
            
            for pred_name, df in self.data.items():
                avg_ipc = df['IPC'].mean()
                med_ipc = df['IPC'].median()
                std_ipc = df['IPC'].std()
                f.write(f"{pred_name:<35} {avg_ipc:<12.4f} {med_ipc:<12.4f} {std_ipc:<12.4f}\n")
            
            # Cycle Analysis
            f.write("\n\nCycle Count Analysis:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Predictor':<35} {'Avg Cycles':<15} {'Total Cycles':<15}\n")
            f.write("-" * 80 + "\n")
            
            for pred_name, df in self.data.items():
                avg_cyc = df['Cycles'].mean()
                tot_cyc = df['Cycles'].sum()
                f.write(f"{pred_name:<35} {avg_cyc:<15,.0f} {tot_cyc:<15,.0f}\n")
            
            # Branches Per Cycle
            f.write("\n\nBranches Per Cycle Analysis:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Predictor':<35} {'Avg BrPerCyc':<15}\n")
            f.write("-" * 80 + "\n")
            
            for pred_name, df in self.data.items():
                avg_bpc = df['BrPerCyc'].mean()
                f.write(f"{pred_name:<35} {avg_bpc:<15.6f}\n")
            
            # Wasted Cycles
            f.write("\n\nWasted Cycles Analysis (CycWPPKI):\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Predictor':<35} {'Avg CycWPPKI':<15}\n")
            f.write("-" * 80 + "\n")
            
            for pred_name, df in self.data.items():
                avg_cwp = df['CycWPPKI'].mean()
                f.write(f"{pred_name:<35} {avg_cwp:<15.4f}\n")
                
        print(f"  ✓ Performance metrics saved to {output_file}")
        
        # Generate visualization
        self._plot_performance_comparison()
        
    def _plot_performance_comparison(self):
        """Plot performance metrics comparison"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Performance Metrics Comparison', fontsize=16, fontweight='bold')
        
        pred_names = list(self.data.keys())
        colors = plt.cm.Set3(np.linspace(0, 1, len(pred_names)))
        
        # 1. IPC Comparison
        ax = axes[0, 0]
        avg_ipc = [self.data[name]['IPC'].mean() for name in pred_names]
        bars = ax.bar(range(len(pred_names)), avg_ipc, color=colors)
        ax.set_xticks(range(len(pred_names)))
        ax.set_xticklabels(pred_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Average IPC')
        ax.set_title('Instructions Per Cycle (Higher is Better)')
        ax.grid(True, alpha=0.3, axis='y')
        
        for bar, val in zip(bars, avg_ipc):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.3f}', ha='center', va='bottom', fontsize=9)
        
        # 2. Average Cycles
        ax = axes[0, 1]
        avg_cycles = [self.data[name]['Cycles'].mean() for name in pred_names]
        bars = ax.bar(range(len(pred_names)), avg_cycles, color=colors)
        ax.set_xticks(range(len(pred_names)))
        ax.set_xticklabels(pred_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Average Cycles')
        ax.set_title('Average Cycle Count (Lower is Better)')
        ax.grid(True, alpha=0.3, axis='y')
        
        for bar, val in zip(bars, avg_cycles):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val/1e6:.2f}M', ha='center', va='bottom', fontsize=9)
        
        # 3. Branches Per Cycle
        ax = axes[1, 0]
        avg_bpc = [self.data[name]['BrPerCyc'].mean() for name in pred_names]
        bars = ax.bar(range(len(pred_names)), avg_bpc, color=colors)
        ax.set_xticks(range(len(pred_names)))
        ax.set_xticklabels(pred_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Branches Per Cycle')
        ax.set_title('Average Branches Per Cycle')
        ax.grid(True, alpha=0.3, axis='y')
        
        for bar, val in zip(bars, avg_bpc):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.4f}', ha='center', va='bottom', fontsize=9)
        
        # 4. Wasted Cycles
        ax = axes[1, 1]
        avg_cwp = [self.data[name]['CycWPPKI'].mean() for name in pred_names]
        bars = ax.bar(range(len(pred_names)), avg_cwp, color=colors)
        ax.set_xticks(range(len(pred_names)))
        ax.set_xticklabels(pred_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('CycWPPKI')
        ax.set_title('Cycles Wasted Per Kilo Instructions (Lower is Better)')
        ax.grid(True, alpha=0.3, axis='y')
        
        for bar, val in zip(bars, avg_cwp):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.2f}', ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        output_file = self.output_dir / 'overall' / 'figures' / '03_performance_metrics.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Performance metrics plot saved to {output_file}")
        
    def _generate_category_analysis(self):
        """Analyze performance by trace category"""
        output_file = self.output_dir / 'overall' / 'text_reports' / '04_category_analysis.txt'
        
        # Extract category from trace name
        def get_category(trace_name):
            if 'int_' in trace_name:
                return 'Integer'
            elif 'fp_' in trace_name:
                return 'Floating Point'
            elif 'web_' in trace_name:
                return 'Web'
            elif 'compress_' in trace_name:
                return 'Compress'
            elif 'infra_' in trace_name:
                return 'Infrastructure'
            elif 'media_' in trace_name:
                return 'Media'
            else:
                return 'Other'
        
        # Add category to each dataframe
        for pred_name, df in self.data.items():
            df['Category'] = df['Run'].apply(get_category)
        
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("TRACE CATEGORY ANALYSIS\n")
            f.write("=" * 80 + "\n\n")
            
            categories = set()
            for df in self.data.values():
                categories.update(df['Category'].unique())
            categories = sorted(categories)
            
            for category in categories:
                f.write(f"\n{category} Traces:\n")
                f.write("-" * 80 + "\n")
                f.write(f"{'Predictor':<35} {'Count':<8} {'Avg MPKI':<12} {'Avg IPC':<12}\n")
                f.write("-" * 80 + "\n")
                
                for pred_name, df in self.data.items():
                    cat_df = df[df['Category'] == category]
                    if len(cat_df) > 0:
                        count = len(cat_df)
                        avg_mpki = cat_df['MPKI'].mean()
                        avg_ipc = cat_df['IPC'].mean()
                        f.write(f"{pred_name:<35} {count:<8} {avg_mpki:<12.4f} {avg_ipc:<12.4f}\n")
                        
        print(f"  ✓ Category analysis saved to {output_file}")
        
        # Generate visualization
        self._plot_category_analysis()
        
    def _plot_category_analysis(self):
        """Plot category-wise comparison"""
        categories = set()
        for df in self.data.values():
            if 'Category' in df.columns:
                categories.update(df['Category'].unique())
        categories = sorted(categories)
        
        if len(categories) == 0:
            return
        
        fig, axes = plt.subplots(2, 1, figsize=(14, 12))
        fig.suptitle('Performance by Trace Category', fontsize=16, fontweight='bold')
        
        # 1. MPKI by Category
        ax = axes[0]
        x = np.arange(len(categories))
        width = 0.15
        
        for idx, (pred_name, df) in enumerate(self.data.items()):
            mpki_by_cat = []
            for cat in categories:
                cat_df = df[df['Category'] == cat]
                mpki = cat_df['MPKI'].mean() if len(cat_df) > 0 else 0
                mpki_by_cat.append(mpki)
            
            offset = (idx - len(self.data) / 2) * width
            ax.bar(x + offset, mpki_by_cat, width, label=pred_name)
        
        ax.set_xlabel('Category')
        ax.set_ylabel('Average MPKI')
        ax.set_title('MPKI by Trace Category')
        ax.set_xticks(x)
        ax.set_xticklabels(categories, rotation=45, ha='right')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')
        
        # 2. IPC by Category
        ax = axes[1]
        for idx, (pred_name, df) in enumerate(self.data.items()):
            ipc_by_cat = []
            for cat in categories:
                cat_df = df[df['Category'] == cat]
                ipc = cat_df['IPC'].mean() if len(cat_df) > 0 else 0
                ipc_by_cat.append(ipc)
            
            offset = (idx - len(self.data) / 2) * width
            ax.bar(x + offset, ipc_by_cat, width, label=pred_name)
        
        ax.set_xlabel('Category')
        ax.set_ylabel('Average IPC')
        ax.set_title('IPC by Trace Category')
        ax.set_xticks(x)
        ax.set_xticklabels(categories, rotation=45, ha='right')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        output_file = self.output_dir / 'overall' / 'figures' / '04_category_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Category analysis plot saved to {output_file}")
        
    def _generate_rankings(self):
        """Generate rankings of predictors across metrics"""
        output_file = self.output_dir / 'overall' / 'text_reports' / '05_predictor_rankings.txt'
        
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("PREDICTOR RANKINGS\n")
            f.write("=" * 80 + "\n\n")
            
            # Metrics for ranking (lower is better for some, higher for others)
            ranking_metrics = {
                'MPKI': ('lower', 'Mispredictions Per Kilo Instructions'),
                'MR': ('lower', 'Miss Rate'),
                'IPC': ('higher', 'Instructions Per Cycle'),
                'Cycles': ('lower', 'Average Cycles'),
                'CycWPPKI': ('lower', 'Cycles Wasted Per Kilo Instructions'),
            }
            
            for metric, (direction, description) in ranking_metrics.items():
                f.write(f"\n{description} ({direction} is better):\n")
                f.write("-" * 80 + "\n")
                
                # Calculate average for each predictor
                results = []
                for pred_name, df in self.data.items():
                    if metric in df.columns:
                        avg_val = df[metric].mean()
                        results.append((pred_name, avg_val))
                
                # Sort based on direction
                results.sort(key=lambda x: x[1], reverse=(direction == 'higher'))
                
                # Print rankings
                for rank, (pred_name, value) in enumerate(results, 1):
                    medal = "🥇" if rank == 1 else "🥈" if rank == 2 else "🥉" if rank == 3 else "  "
                    f.write(f"  {rank}. {medal} {pred_name:<40} {value:>12.4f}\n")
            
            # Overall Score (composite ranking)
            f.write("\n\n" + "=" * 80 + "\n")
            f.write("COMPOSITE RANKING (Lower score is better)\n")
            f.write("=" * 80 + "\n\n")
            
            composite_scores = {}
            for pred_name in self.data.keys():
                composite_scores[pred_name] = 0
            
            # Calculate composite score from all rankings
            for metric, (direction, _) in ranking_metrics.items():
                results = []
                for pred_name, df in self.data.items():
                    if metric in df.columns:
                        avg_val = df[metric].mean()
                        results.append((pred_name, avg_val))
                
                results.sort(key=lambda x: x[1], reverse=(direction == 'higher'))
                
                for rank, (pred_name, _) in enumerate(results, 1):
                    composite_scores[pred_name] += rank
            
            # Sort by composite score
            final_ranking = sorted(composite_scores.items(), key=lambda x: x[1])
            
            for rank, (pred_name, score) in enumerate(final_ranking, 1):
                medal = "🥇" if rank == 1 else "🥈" if rank == 2 else "🥉" if rank == 3 else "  "
                f.write(f"  {rank}. {medal} {pred_name:<40} Score: {score}\n")
                
        print(f"  ✓ Rankings saved to {output_file}")
        
        # Generate visualization
        self._plot_rankings(final_ranking)
        
    def _plot_rankings(self, final_ranking):
        """Plot overall rankings"""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        pred_names = [x[0] for x in final_ranking]
        scores = [x[1] for x in final_ranking]
        
        colors = plt.cm.RdYlGn_r(np.linspace(0.3, 0.9, len(pred_names)))
        bars = ax.barh(range(len(pred_names)), scores, color=colors)
        
        ax.set_yticks(range(len(pred_names)))
        ax.set_yticklabels(pred_names)
        ax.set_xlabel('Composite Score (Lower is Better)')
        ax.set_title('Overall Predictor Rankings - Composite Score', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='x')
        
        # Add score labels
        for i, (bar, score) in enumerate(zip(bars, scores)):
            width = bar.get_width()
            ax.text(width, bar.get_y() + bar.get_height()/2.,
                   f' {score}', ha='left', va='center', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        output_file = self.output_dir / 'overall' / 'figures' / '05_overall_rankings.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Rankings plot saved to {output_file}")
        
    def generate_pairwise_comparisons(self):
        """Generate pairwise comparisons between predictors"""
        print("\nGenerating pairwise comparisons...")
        
        predictor_names = list(self.data.keys())
        
        for i, pred1 in enumerate(predictor_names):
            for pred2 in predictor_names[i+1:]:
                print(f"\n  Comparing {pred1} vs {pred2}...")
                self._compare_pair(pred1, pred2)
                
    def _compare_pair(self, pred1, pred2):
        """Compare two predictors in detail"""
        df1 = self.data[pred1]
        df2 = self.data[pred2]
        
        pair_dir = self.output_dir / f'{pred1}_vs_{pred2}'
        
        # 1. Generate comparison report
        self._generate_pair_report(pred1, pred2, df1, df2, pair_dir)
        
        # 2. Generate comparison plots
        self._generate_pair_plots(pred1, pred2, df1, df2, pair_dir)
        
        # 3. Generate trace-by-trace comparison
        self._generate_trace_comparison(pred1, pred2, df1, df2, pair_dir)
        
    def _generate_pair_report(self, pred1, pred2, df1, df2, pair_dir):
        """Generate text report for predictor pair"""
        output_file = pair_dir / 'text_reports' / '01_comparison_summary.txt'
        
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"COMPARISON: {pred1} vs {pred2}\n")
            f.write("=" * 80 + "\n\n")
            
            # Summary statistics comparison
            f.write("Summary Statistics:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Metric':<30} {pred1:<20} {pred2:<20} {'Difference':<15}\n")
            f.write("-" * 80 + "\n")
            
            metrics = ['MPKI', 'MR', 'IPC', 'Cycles', 'BrPerCyc', 'MispBrPerCyc', 'CycWPPKI']
            
            for metric in metrics:
                if metric in df1.columns and metric in df2.columns:
                    val1 = df1[metric].mean()
                    val2 = df2[metric].mean()
                    diff = val2 - val1
                    pct_diff = (diff / val1 * 100) if val1 != 0 else 0
                    
                    f.write(f"{metric:<30} {val1:<20.4f} {val2:<20.4f} ")
                    if metric in ['MPKI', 'MR', 'Cycles', 'CycWPPKI']:
                        # Lower is better
                        if diff < 0:
                            f.write(f"↓ {abs(pct_diff):.2f}% better\n")
                        else:
                            f.write(f"↑ {pct_diff:.2f}% worse\n")
                    else:
                        # Higher is better
                        if diff > 0:
                            f.write(f"↑ {pct_diff:.2f}% better\n")
                        else:
                            f.write(f"↓ {abs(pct_diff):.2f}% worse\n")
            
            # Wins/Losses
            f.write("\n\nTrace-by-Trace Comparison (MPKI):\n")
            f.write("-" * 80 + "\n")
            
            # Merge dataframes on Run
            merged = df1.merge(df2, on='Run', suffixes=('_1', '_2'))
            
            wins_pred1 = (merged['MPKI_1'] < merged['MPKI_2']).sum()
            wins_pred2 = (merged['MPKI_1'] > merged['MPKI_2']).sum()
            ties = (merged['MPKI_1'] == merged['MPKI_2']).sum()
            
            f.write(f"{pred1} wins: {wins_pred1} traces\n")
            f.write(f"{pred2} wins: {wins_pred2} traces\n")
            f.write(f"Ties: {ties} traces\n")
            
            # Biggest improvements and regressions
            merged['MPKI_diff'] = merged['MPKI_2'] - merged['MPKI_1']
            
            f.write(f"\n\nBiggest Improvements ({pred2} better):\n")
            f.write("-" * 80 + "\n")
            top_improvements = merged.nsmallest(5, 'MPKI_diff')
            for _, row in top_improvements.iterrows():
                f.write(f"  {row['Run']:<40} {row['MPKI_1']:>8.4f} → {row['MPKI_2']:>8.4f} "
                       f"(↓ {abs(row['MPKI_diff']):.4f})\n")
            
            f.write(f"\n\nBiggest Regressions ({pred2} worse):\n")
            f.write("-" * 80 + "\n")
            top_regressions = merged.nlargest(5, 'MPKI_diff')
            for _, row in top_regressions.iterrows():
                f.write(f"  {row['Run']:<40} {row['MPKI_1']:>8.4f} → {row['MPKI_2']:>8.4f} "
                       f"(↑ {row['MPKI_diff']:.4f})\n")
                       
        print(f"    ✓ Comparison summary saved")
        
    def _generate_pair_plots(self, pred1, pred2, df1, df2, pair_dir):
        """Generate comparison plots for predictor pair"""
        
        # Merge dataframes
        merged = df1.merge(df2, on='Run', suffixes=('_1', '_2'))
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f'{pred1} vs {pred2}', fontsize=16, fontweight='bold')
        
        # 1. MPKI Scatter Plot
        ax = axes[0, 0]
        ax.scatter(merged['MPKI_1'], merged['MPKI_2'], alpha=0.6, s=50)
        
        # Add diagonal line
        max_val = max(merged['MPKI_1'].max(), merged['MPKI_2'].max())
        ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='Equal performance')
        
        ax.set_xlabel(f'{pred1} MPKI')
        ax.set_ylabel(f'{pred2} MPKI')
        ax.set_title('MPKI Comparison (Points below line = Pred2 better)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. IPC Scatter Plot
        ax = axes[0, 1]
        ax.scatter(merged['IPC_1'], merged['IPC_2'], alpha=0.6, s=50, color='green')
        
        max_val = max(merged['IPC_1'].max(), merged['IPC_2'].max())
        ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='Equal performance')
        
        ax.set_xlabel(f'{pred1} IPC')
        ax.set_ylabel(f'{pred2} IPC')
        ax.set_title('IPC Comparison (Points above line = Pred2 better)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 3. MPKI Difference Distribution
        ax = axes[1, 0]
        mpki_diff = merged['MPKI_2'] - merged['MPKI_1']
        ax.hist(mpki_diff, bins=30, alpha=0.7, edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', linewidth=2, label='No difference')
        ax.axvline(mpki_diff.mean(), color='green', linestyle='-', linewidth=2, label=f'Mean diff: {mpki_diff.mean():.4f}')
        ax.set_xlabel('MPKI Difference (Pred2 - Pred1)')
        ax.set_ylabel('Frequency')
        ax.set_title('MPKI Difference Distribution (Negative = Pred2 better)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 4. Metric Comparison Bar Chart
        ax = axes[1, 1]
        metrics = ['MPKI', 'MR', 'CycWPPKI']
        x = np.arange(len(metrics))
        width = 0.35
        
        vals1 = [df1[m].mean() for m in metrics]
        vals2 = [df2[m].mean() for m in metrics]
        
        # Normalize for visualization
        max_vals = [max(v1, v2) for v1, v2 in zip(vals1, vals2)]
        norm_vals1 = [v1/mv for v1, mv in zip(vals1, max_vals)]
        norm_vals2 = [v2/mv for v2, mv in zip(vals2, max_vals)]
        
        bars1 = ax.bar(x - width/2, norm_vals1, width, label=pred1, alpha=0.8)
        bars2 = ax.bar(x + width/2, norm_vals2, width, label=pred2, alpha=0.8)
        
        ax.set_ylabel('Normalized Value')
        ax.set_title('Key Metrics Comparison (Lower is Better)')
        ax.set_xticks(x)
        ax.set_xticklabels(metrics)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        output_file = pair_dir / 'figures' / '01_comparison_plots.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    ✓ Comparison plots saved")
        
    def _generate_trace_comparison(self, pred1, pred2, df1, df2, pair_dir):
        """Generate detailed trace-by-trace comparison"""
        output_file = pair_dir / 'text_reports' / '02_trace_by_trace.txt'
        
        merged = df1.merge(df2, on='Run', suffixes=('_1', '_2'))
        merged['MPKI_diff'] = merged['MPKI_2'] - merged['MPKI_1']
        merged['IPC_diff'] = merged['IPC_2'] - merged['IPC_1']
        
        # Sort by MPKI difference
        merged_sorted = merged.sort_values('MPKI_diff')
        
        with open(output_file, 'w') as f:
            f.write("=" * 100 + "\n")
            f.write(f"TRACE-BY-TRACE COMPARISON: {pred1} vs {pred2}\n")
            f.write("=" * 100 + "\n\n")
            
            f.write(f"{'Trace':<40} {pred1+' MPKI':<12} {pred2+' MPKI':<12} {'Diff':<12} {'Status':<20}\n")
            f.write("-" * 100 + "\n")
            
            for _, row in merged_sorted.iterrows():
                trace = row['Run']
                mpki1 = row['MPKI_1']
                mpki2 = row['MPKI_2']
                diff = row['MPKI_diff']
                
                if diff < -0.1:
                    status = f"{pred2} better"
                    marker = "✓"
                elif diff > 0.1:
                    status = f"{pred1} better"
                    marker = "✗"
                else:
                    status = "Similar"
                    marker = "≈"
                
                f.write(f"{trace:<40} {mpki1:<12.4f} {mpki2:<12.4f} {diff:<12.4f} {marker} {status:<18}\n")
                
        print(f"    ✓ Trace-by-trace comparison saved")
        
    def run_all_analyses(self):
        """Run all analyses and generate all reports"""
        print("\n" + "=" * 80)
        print("COMPREHENSIVE PREDICTOR COMPARISON")
        print("=" * 80)
        
        if len(self.data) == 0:
            print("ERROR: No predictor data loaded!")
            return
        
        # Overall comparison
        self.generate_overall_comparison()
        
        # Pairwise comparisons
        self.generate_pairwise_comparisons()
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"\nAll reports and figures saved to: {self.output_dir}")
        print("\nDirectory structure:")
        print("  - overall/")
        print("    - text_reports/")
        print("    - figures/")
        print("  - <Predictor1>_vs_<Predictor2>/")
        print("    - text_reports/")
        print("    - figures/")

if __name__ == '__main__':
    comparator = PredictorComparator()
    comparator.run_all_analyses()
