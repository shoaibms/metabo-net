"""
Transformation Evaluation Script for Metabolomics Data Analysis
------------------------------------------------------------

This script evaluates the effectiveness of data transformations (specifically asinh) 
by generating two key visualisations:

1. MA Plot (Mean-Difference Plot):
   - Visualises intensity-dependent ratio of two datasets
   - X-axis: Average of log intensities (A)
   - Y-axis: Difference of log intensities (M)
   - Helps identify intensity-dependent bias

2. Relative Standard Deviation (RSD) Plot:
   - Compares data variability before and after transformation
   - Uses violin plots to show full distribution
   - Includes swarm plots for individual data points

The analysis focuses on comparing original vs asinh-transformed data to assess
transformation effectiveness in improving data properties.

Required Libraries:
- pandas: Data manipulation
- numpy: Numerical operations
- matplotlib: Basic plotting
- seaborn: Statistical visualisation
- statsmodels: LOWESS smoothing for MA plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import os

def load_data(file_path):
    """
    Load data from CSV file.
    
    Args:
        file_path (str): Path to CSV file
        
    Returns:
        pandas.DataFrame: Loaded data
    """
    return pd.read_csv(file_path)

def extract_n_cluster_data(df):
    """
    Extract columns containing N_Cluster data.
    
    Args:
        df (pandas.DataFrame): Input dataframe
        
    Returns:
        pandas.DataFrame: Subset containing only N_Cluster columns
    """
    return df.filter(regex="N_Cluster")

def plot_ma_transform(before_df, after_df, title, cmap='Greens', 
                     title_size=16, label_size=15, tick_size=14, 
                     legend_size=14, cbar_label_size=14):
    """
    Generate MA plots comparing before and after transformation.
    
    MA plots (Bland-Altman plots) are used to visualise intensity-dependent
    ratio of two datasets, helpful for identifying systematic bias.
    
    Args:
        before_df (pandas.DataFrame): Original data
        after_df (pandas.DataFrame): Transformed data
        title (str): Plot title
        cmap (str): Colormap for scatter plots
        Various size parameters for plot elements
    """
    # Add 1 to avoid log(0)
    before = before_df + 1
    after = after_df + 1

    # Calculate M (log ratio) and A (average log intensity) values
    M_before = np.log2(before.iloc[:, 0]) - np.log2(before.iloc[:, 1])
    A_before = 0.5 * (np.log2(before.iloc[:, 0]) + np.log2(before.iloc[:, 1]))

    M_after = np.log2(after.iloc[:, 0]) - np.log2(after.iloc[:, 1])
    A_after = 0.5 * (np.log2(after.iloc[:, 0]) + np.log2(after.iloc[:, 1]))

    # Create comparison plots
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))

    # Plot before transformation
    scatter1 = axes[0].scatter(A_before, M_before, alpha=0.5, c=M_before, cmap=cmap)
    lowess_before = sm.nonparametric.lowess(M_before, A_before, frac=0.3)
    axes[0].plot(lowess_before[:, 0], lowess_before[:, 1], color='red')
    axes[0].set_title(f'MA-transform Before Transformation - {title}', fontsize=title_size)
    axes[0].set_xlabel('A (average log-intensity)', fontsize=label_size)
    axes[0].set_ylabel('M (log ratio)', fontsize=label_size)
    axes[0].tick_params(axis='both', which='major', labelsize=tick_size)
    cbar1 = fig.colorbar(scatter1, ax=axes[0], label='M Before')
    cbar1.ax.tick_params(labelsize=legend_size)
    cbar1.set_label('M Before', fontsize=cbar_label_size)

    # Plot after transformation
    scatter2 = axes[1].scatter(A_after, M_after, alpha=0.5, c=M_after, cmap=cmap)
    lowess_after = sm.nonparametric.lowess(M_after, A_after, frac=0.3)
    axes[1].plot(lowess_after[:, 0], lowess_after[:, 1], color='red')
    axes[1].set_title(f'MA-transform After Transformation - {title}', fontsize=title_size)
    axes[1].set_xlabel('A (average log-intensity)', fontsize=label_size)
    axes[1].set_ylabel('M (log ratio)', fontsize=label_size)
    axes[1].tick_params(axis='both', which='major', labelsize=tick_size)
    cbar2 = fig.colorbar(scatter2, ax=axes[1], label='M After')
    cbar2.ax.tick_params(labelsize=legend_size)
    cbar2.set_label('M After', fontsize=cbar_label_size)

    plt.tight_layout()
    plt.show()

def calculate_rsd(data):
    """
    Calculate Relative Standard Deviation.
    
    RSD = (Standard Deviation / Mean) * 100
    
    Args:
        data (pandas.DataFrame): Input data
        
    Returns:
        numpy.array: RSD values
    """
    mean = np.mean(data, axis=0)
    std_dev = np.std(data, axis=0)
    rsd = (std_dev / mean) * 100
    return rsd

def plot_rsd_violin(before_df, after_df, title, 
                   before_color='#f1facf', after_color='#68f7d8', 
                   title_size=18, label_size=17, tick_size=16):
    """
    Generate violin plots comparing RSD distributions before and after transformation.
    
    Args:
        before_df (pandas.DataFrame): Original data
        after_df (pandas.DataFrame): Transformed data
        title (str): Plot title
        Various color and size parameters for plot elements
    
    Returns:
        tuple: RSD values before and after transformation
    """
    # Calculate RSD for both datasets
    rsd_before = calculate_rsd(before_df)
    rsd_after = calculate_rsd(after_df)

    # Prepare data for plotting
    results_df = pd.DataFrame({
        'RSD': np.concatenate([rsd_before, rsd_after]),
        'Condition': ['Before'] * len(rsd_before) + ['After'] * len(rsd_after)
    })

    # Remove infinities and NaNs
    results_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    results_df.dropna(inplace=True)

    # Create violin plot with swarm overlay
    plt.figure(figsize=(7, 4))
    sns.violinplot(x='Condition', y='RSD', data=results_df, 
                  inner=None, palette=[before_color, after_color])
    sns.swarmplot(x='Condition', y='RSD', data=results_df, 
                 color='k', alpha=0.5)
    
    plt.title(f'Relative Standard Deviation (RSD) - {title}', fontsize=title_size)
    plt.xlabel('Data type', fontsize=label_size)
    plt.ylabel('RSD (%)', fontsize=label_size)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.tight_layout()
    plt.show()

    return rsd_before, rsd_after

def save_metrics_to_csv(metrics, file_path):
    """
    Save calculated metrics to CSV file.
    
    Args:
        metrics (dict): Dictionary containing metric values
        file_path (str): Path for output CSV
    """
    metrics_df = pd.DataFrame.from_dict(metrics, orient='index').transpose()
    metrics_df.to_csv(file_path, index=False)
    print(f"Metrics saved to {file_path}")

def main():
    """Main execution function."""
    # Define file paths
    file_paths = {
        'original': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if.csv",
        'asinh': r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_asinh.csv"
    }

    # Load original data and calculate median
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
        rsd_before, rsd_after = plot_rsd_violin(original_n_cluster, 
                                              transformed_n_cluster, title)

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