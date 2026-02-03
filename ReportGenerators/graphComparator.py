import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - no popups
import matplotlib.pyplot as plt
import numpy as np

# Load Data
df = pd.read_csv('results.csv')

# Clean Data
if df['MR'].dtype == object:
    df['MR'] = df['MR'].str.replace('%', '').astype(float)

def generate_comparison_plot(data, workload_type, full_data=None):
    if data.empty:
        print(f"No data found for {workload_type}")
        return

    data = data.sort_values(by='MPKI', ascending=False)
    runs = data['Run']
    total_mpki = data['MPKI']
    steady_mpki = data['50PercMPKI']
    
    # Calculate geometric means
    if full_data is not None:
        mean_src_mpki = full_data['MPKI']
        mean_src_steady = full_data['50PercMPKI']
    else:
        mean_src_mpki = total_mpki
        mean_src_steady = steady_mpki
        
    geometric_mean_mpki = (np.prod(mean_src_mpki))**(1/len(mean_src_mpki))
    print(f"GeoMean MPKI for {workload_type}: {geometric_mean_mpki:.2f}")
    
    geometric_mean_50perc_mpki = (np.prod(mean_src_steady))**(1/len(mean_src_steady))
    print(f"GeoMean 50PercMPKI for {workload_type}: {geometric_mean_50perc_mpki:.2f}")
    
    # Chart setup
    x = np.arange(len(runs))
    width = 0.35
    
    fig_width = max(12, len(runs) * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 8))
    
    # Plot bars
    ax.bar(x - width/2, total_mpki, width, label='Total MPKI', color='#2c3e50')
    ax.bar(x + width/2, steady_mpki, width, label='Steady-State (50%) MPKI', color='#e67e22')
    
    ax.axhline(y=geometric_mean_mpki, color='#2980b9', linestyle='-', label=f'GeoMean Total: {geometric_mean_mpki:.2f}')
    ax.axhline(y=geometric_mean_50perc_mpki, color='#c0392b', linestyle='-', label=f'GeoMean 50%: {geometric_mean_50perc_mpki:.2f}')
    
    ax.set_ylabel('MPKI')
    ax.set_title(f'MPKI Comparison: {workload_type}')
    ax.set_xticks(x)
    ax.set_xticklabels(runs, rotation=45, ha='right')
    ax.legend()
    ax.grid(axis='y', linestyle='--', alpha=0.4)
    plt.subplots_adjust(bottom=0.25)
    
    # Save to Reports directory
    import os
    output_dir = '../Reports/01_misprediction_analysis/graphs'
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f'{output_dir}/mpki_comparison_{workload_type}.png')
    plt.close()  # Close figure instead of showing

# --- FIX STARTS HERE ---
selected_runs = []

def get_spaced_items(items, n):
    if len(items) == 0: return []
    if len(items) <= n: return list(items)
    step = len(items) / n
    return [items[int(i * step)] for i in range(n)]

# We iterate through the desired prefixes ('int', 'fp', etc.)
# And check if the 'Run' column string starts with that prefix
categories = [('int', 6), ('fp', 6), ('web', 1), ('infra', 1), ('media', 1)]

for workload, count in categories:
    # Filter by checking the string content of 'Run' rather than the 'Workload' column
    # We add '_' to ensure 'int' doesn't match 'integer_math' if that existed, just 'int_...'
    runs = df[df['Run'].str.startswith(workload + '_')]['Run'].unique()
    
    # If your CSV names are inconsistent (e.g. sometimes 'int-11'), remove the + '_'
    if len(runs) == 0:
         runs = df[df['Run'].str.contains(workload)]['Run'].unique()

    selected_runs.extend(get_spaced_items(runs, count))

plot_data = df[df['Run'].isin(selected_runs)]
# --- FIX ENDS HERE ---

generate_comparison_plot(plot_data, "Representative_Mix", full_data=df)