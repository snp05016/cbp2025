import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# ==========================================
# CONFIGURATION & CONSTANTS
# ==========================================

# Input Files (must be in the same directory as script)
FILES = {
    'Baseline': 'baseline.csv',
    'LVCP': 'load-value-correlator-man-results.csv',
    'RVA-Toru': 'register-value-aware-toru-results.csv',
    'TAGE-SCL (Andrez)': 'tage-scl-andrez-seznec-results.csv',
    'TAGE-SC-L (Alberto)': 'tage-sc-l-alberto-ros-results.csv'
}

# Output Directory
OUTPUT_DIR = 'subsumption_analysis/figures/'

# Plot Styling for High Readability
# Removed 'grid.axis' to fix KeyError
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 14,
    'axes.labelsize': 15,
    'axes.titlesize': 18,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'axes.grid': False,       # We will turn this on manually for Y-axis only
    'grid.alpha': 0.4,
    'grid.linestyle': '--'
})

# ==========================================
# DATA PROCESSING
# ==========================================

def load_and_merge_data():
    """
    Loads CSVs, extracts 'Run' and 'MPKI', merges on 'Run'.
    """
    dfs = []
    
    try:
        base_df = pd.read_csv(FILES['Baseline'])
        base_df = base_df[['Run', 'MPKI']].rename(columns={'MPKI': 'Baseline'})
        dfs.append(base_df)
    except FileNotFoundError:
        print(f"Error: mandatory file {FILES['Baseline']} not found.")
        return None

    for name, filename in FILES.items():
        if name == 'Baseline': continue
        if not os.path.exists(filename):
            print(f"Error: file {filename} not found.")
            continue
        df = pd.read_csv(filename)
        df = df[['Run', 'MPKI']].rename(columns={'MPKI': name})
        dfs.append(df)

    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='Run', how='inner')

    return merged_df

def calculate_deltas(df):
    """Calculating Delta = MPKI_A - MPKI_B"""
    df['LVCP_to_RVAToru'] = df['LVCP'] - df['RVA-Toru']
    df['Andrez_to_RVAToru'] = df['TAGE-SCL (Andrez)'] - df['RVA-Toru']
    df['Andrez_to_Alberto'] = df['TAGE-SCL (Andrez)'] - df['TAGE-SC-L (Alberto)']
    return df

# ==========================================
# PLOTTING FUNCTION (Violin Plot)
# ==========================================

def plot_readable_violin_summary(df):
    """
    Generates a highly readable violin plot showing delta distributions.
    """
    comparisons = [
        ('LVCP_to_RVAToru', 'LVCP\n↓\nRVA-Toru'),
        ('Andrez_to_RVAToru', 'TAGE-SCL (Andrez)\n↓\nRVA-Toru'),
        ('Andrez_to_Alberto', 'TAGE-SCL (Andrez)\n↓\nTAGE-SC-L (Alberto)')
    ]

    # Prepare data: filter zeros and create list of arrays
    data_to_plot = []
    labels = []
    medians = []
    
    for col, label in comparisons:
        clean_data = df[col][df[col] != 0].values
        data_to_plot.append(clean_data)
        labels.append(label)
        medians.append(np.median(clean_data))

    # Setup Figure
    fig, ax = plt.subplots(figsize=(11, 7))

    # Manually enable grid for Y-axis only (Replaces the failed rcParams setting)
    ax.grid(True, axis='y', linestyle='--', alpha=0.4)

    # Create Violin Plot
    parts = ax.violinplot(data_to_plot, showmeans=False, showmedians=False, vert=True)

    # --- Custom Styling for Readability ---
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c'] # Distinct color for each comparison

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)

    # Add standard lines (min/max/bars) back in, styled black
    for partname in ('cbars', 'cmins', 'cmaxes'):
        v = parts[partname]
        v.set_edgecolor('black')
        v.set_linewidth(1)

    # Add clear white dots for the medians
    inds = np.arange(1, len(medians) + 1)
    ax.scatter(inds, medians, marker='o', color='white', s=60, zorder=3, edgecolors='black', linewidth=1.5, label='Median Δ')

    # Add heavy zero line
    ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.8)

    # Labels and Ticks
    ax.set_xticks(inds)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r'$\Delta$ MPKI (Higher is better for 2nd predictor)')
    ax.set_title('Delta MPKI Distribution Summary (Violin View)', pad=20)
    
    # Add an annotation explaining the zero line
    ax.text(len(labels) + 0.6, 0, 'Zero Gain\n(No Difference)', 
            verticalalignment='center', fontsize=11, color='black', style='italic')
    
    # Adjust layout to prevent cutting off labels
    plt.tight_layout()

    # Save
    filepath = os.path.join(OUTPUT_DIR, 'delta_violin_summary_readable.png')
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated readable graph: {filepath}")

# ==========================================
# MAIN
# ==========================================

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print("Processing data...")
    df = load_and_merge_data()
    if df is None: return
    df = calculate_deltas(df)
    
    print("Generating highly readable violin plot...")
    plot_readable_violin_summary(df)
    print("Done.")

if __name__ == "__main__":
    main()