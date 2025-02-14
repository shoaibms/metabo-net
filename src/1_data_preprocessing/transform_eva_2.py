"""
Transformation Evaluation Script for Metabolomics Data
---------------------------------------------------

This script evaluates the effectiveness of data transformations in metabolomics analysis
using MA plots and Relative Standard Deviation (RSD) visualisation approaches.

Key Visualisations:
1. MA Plots: Assess intensity-dependent bias in the data
   - M (log ratio) vs A (average log-intensity)
   - LOWESS smoothing curve shows systematic trends
   
2. Violin Plots: Compare RSD distributions before/after transformation
   - Shows full distribution shape with density estimation
   - Includes swarm overlay for individual data points
   
Methods implemented:
- MA transformation for paired sample comparison
- RSD calculation for measurement precision assessment
- LOWESS smoothing for trend visualisation
- Violin plots with swarm overlay for distribution comparison

Metrics saved:
- Transformation method
- RSD before/after transformation
- rMAD percentage of median before/after
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import os

def load_data(file_path):
    """Load dataset from CSV file."""
    return pd.read_csv(file_path)

def extract_n_cluster_data(df):
    """Extract columns containing N_Cluster pattern."""
    return df.filter(regex="N_Cluster")

def plot_ma_transform(before_df, after_df, title, cmap='Greens', title_size=16, 
                     label_size=14, tick_size=12):
    """
    Create MA plots comparing data before and after transformation.
    
    MA plots are used to visualize intensity-dependent ratio of changes:
    - M (y-axis): log ratio between conditions
    - A (x-axis): average log intensity
    
    Args:
        before_df: DataFrame before transformation
        after_df: DataFrame after transformation
        title: Plot title
        cmap: Colormap for scatter points
        title_size: Font size for title
        label_size: Font size for axis labels
        tick_size: Font size for tick labels
    """
    # Add 1 to avoid log(0)
    before = before_df + 1
    after = after_df + 1

    # Calculate M and A values
    M_before = np.log2(before.iloc[:, 0]) - np.log2(before.iloc[:, 1])
    A_before = 0.5 * (np.log2(before.iloc[:, 0]) + np.log2(before.iloc[:, 1]))
    M_after = np.log2(after.iloc[:, 0]) - np.log2(after.iloc[:, 1])
    A_after = 0.5 * (np.log2(after.iloc[:, 0]) + np.log2(after.iloc[:, 1]))

    # Create subplot figure
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))

    # Plot before transformation
    scatter1 = axes[0].scatter(A_before, M_before, alpha=0.5, c=M_before, cmap=cmap)
    lowess_before = sm.nonparametric.lowess(M_before, A_before, frac=0.3)
    axes[0].plot(lowess_before[:, 0], lowess_before[:, 1], color='red')
    axes[0].set_title(f'MA-transform Before Transformation - {title}', fontsize=title_size)
    axes[0].set_xlabel('A (average log-intensity)', fontsize=label_size)
    axes[0].set_ylabel('M (log ratio)', fontsize=label_size)
    axes[0].tick_params(axis='both', which='major', labelsize=tick_size)
    fig.colorbar(scatter1, ax=axes[0], label='M Before')

    # Plot after transformation
    scatter2 = axes[1].scatter(A_after, M_after, alpha=0.5, c=M_after, cmap=cmap)
    lowess_after = sm.nonparametric.lowess(M_after, A_after, frac=0.3)
    axes[1].plot(lowess_after[:, 0], lowess_after[:, 1], color='red')
    axes[1].set_title(f'MA-transform After Transformation - {title}', fontsize=title_size)
    axes[1].set_xlabel('A (average log-intensity)', fontsize=label_size)
    axes[1].set_ylabel('M (log ratio)', fontsize=label_size)
    axes[1].tick_params(axis='both', which='major', labelsize=tick_size)
    fig.colorbar(scatter2, ax=axes[1], label='M After')

    plt.tight_layout()
    plt.show()

def calculate_rsd(data):
    """
    Calculate Relative Standard Deviation (RSD) for each column.
    
    RSD = (Standard Deviation / Mean) * 100
    """
    mean = np.mean(data, axis=0)
    std_dev = np.std(data, axis=0)
    rsd = (std_dev / mean) * 100
    return rsd

def plot_rsd_violin(before_df, after_df, title, before_color='#649e68', 
                   after_color='#b3e6b5', title_size=16, label_size=14, tick_size=12):
    """
    Create violin plots comparing RSD distributions before and after transformation.
    
    Args:
        before_df, after_df: DataFrames before and after transformation
        title: Plot title
        before_color, after_color: Colors for before/after distributions
        title_size, label_size, tick_size: Font sizes for different elements
        
    Returns:
        Tuple of (rsd_before, rsd_after) arrays
    """
    # Calculate RSD values
    rsd_before = calculate_rsd(before_df)
    rsd_after = calculate_rsd(after_df)

    # Prepare data for plotting
    results_df = pd.DataFrame({
        'RSD': np.concatenate([rsd_before, rsd_after]),
        'Condition': ['Before'] * len(rsd_before) + ['After'] * len(rsd_after)
    })

    # Remove infinite values and NaNs
    results_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    results_df.dropna(inplace=True)

    # Create violin plot
    plt.figure(figsize=(10, 6))
    sns.violinplot(x='Condition', y='RSD', data=results_df, inner=None, 
                  palette=[before_color, after_color])
    sns.swarmplot(x='Condition', y='RSD', data=results_df, color='k', alpha=0.5)
    
    # Set plot attributes
    plt.title(f'Violin Plot of Relative Standard Deviation (RSD) - {title}', 
             fontsize=title_size)
    plt.xlabel('Condition', fontsize=label_size)
    plt.ylabel('RSD (%)', fontsize=label_size)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.show()

    return rsd_before, rsd_after

def save_metrics_to_csv(metrics, file_path):
    """Save transformation metrics to CSV file."""
    metrics_df = pd.DataFrame.from_dict(metrics, orient='index').transpose()
    metrics_df.to_csv(file_path, index=False)
    print(f"Metrics saved to {file_path}")

def main():
    # Define file paths for all transformations
    file_paths = {
        'original': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if.csv",
        'anscombe': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_anscombe.csv",
        'asinh': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_asinh.csv",
        'boxcox': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_boxcox.csv",
        'glog': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_glog.csv",
        'log': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_log.csv",
        'sqrt': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_sqrt.csv",
        'yeojohnson': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_yeojohnson.csv"
    }

    # Load original data
    original_df = load_data(file_paths['original'])
    original_n_cluster = extract_n_cluster_data(original_df)
    original_median = original_n_cluster.median().mean()

    # Initialise metrics dictionary
    metrics = {
        'Transformation': [],
        'RSD_Before': [],
        'RSD_After': [],
        'rMAD_Percentage_of_Median_Before': [],
        'rMAD_Percentage_of_Median_After': []
    }

    # Process each transformation
    for key, path in file_paths.items():
        if key == 'original':
            continue
        
        transformed_df = load_data(path)
        transformed_n_cluster = extract_n_cluster_data(transformed_df)
        title = key.replace('_', ' ').title()

        # Generate plots and calculate metrics
        plot_ma_transform(original_n_cluster, transformed_n_cluster, title)
        rsd_before, rsd_after = plot_rsd_violin(original_n_cluster, transformed_n_cluster, title)

        # Store metrics
        metrics['Transformation'].append(title)
        metrics['RSD_Before'].append(rsd_before.mean())
        metrics['RSD_After'].append(rsd_after.mean())
        metrics['rMAD_Percentage_of_Median_Before'].append(rsd_before.mean())
        metrics['rMAD_Percentage_of_Median_After'].append(rsd_after.mean())

    # Save results
    save_metrics_to_csv(metrics, os.path.join(os.path.dirname(file_paths['original']), 
                                            'transformation_metrics.csv'))

if __name__ == "__main__":
    main()