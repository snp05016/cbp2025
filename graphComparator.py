import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
df = pd.read_csv('results.csv')
if df['MR'].dtype == object:
    df['MR'] = df['MR'].str.replace('%', '').astype(float)
def generate_comparison_plot(data, workload_type):
    data = data.sort_values(by='MPKI', ascending=False)
    runs = data['Run']
    total_mpki = data['MPKI']
    geometric_mean_mpki = (np.prod(total_mpki))**(1/len(total_mpki))
    print(f"GeoMean MPKI for {workload_type}: {geometric_mean_mpki:.2f}")
    steady_mpki = data['50PercMPKI']
    geometric_mean_50perc_mpki = (np.prod(steady_mpki))**(1/len(steady_mpki))
    print(f"GeoMean 50PercMPKI for {workload_type}: {geometric_mean_50perc_mpki:.2f}")
    x = np.arange(len(runs)) * 1.5  # Label locations with spacing
    width = 0.50              # Width of the individual bars
    fig, ax = plt.subplots(figsize=(max(20, len(runs) * 0.4), 8))
    ax.bar(x - width/2, total_mpki, width, label='Total MPKI', color='#2c3e50')
    ax.bar(x + width/2, steady_mpki, width, label='Steady-State (50%) MPKI', color='#e67e22')
    ax.axhline(y=geometric_mean_mpki, color='#2980b9', linestyle='-', label=f'Geometric Mean Total MPKI: {geometric_mean_mpki:.2f}')
    ax.axhline(y=geometric_mean_50perc_mpki, color='#c0392b', linestyle='-', label=f'Geometric Mean 50PercMPKI: {geometric_mean_50perc_mpki:.2f}')
    ax.set_ylabel('MPKI (Mispredictions Per Kilo-Instruction)')
    ax.set_title(f'MPKI vs. 50PercMPKI {workload_type.upper()}')
    ax.set_xticks(x)
    ax.set_xticklabels(runs, rotation=90, ha='right', fontsize=8)
    ax.legend()
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    plt.subplots_adjust(bottom=0.25)
    plt.savefig(f'mpki_comparison_{workload_type}.png')

for workload in df['Workload'].unique():
    workload_data = df[df['Workload'] == workload]
    generate_comparison_plot(workload_data, workload)
