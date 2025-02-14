"""
Metabolomics Data Transformation Evaluation Script
-----------------------------------------------

This script evaluates the effectiveness of various data transformations by comparing
statistical metrics and visualizations before and after transformation. It generates
multiple plots and calculates key metrics:

1. Coefficient of Variation (CV)
2. MA plots (intensity-dependent ratio plots)
3. Relative Standard Deviation (RSD)
4. Relative Median Absolute Deviation (rMAD)

Note: While this script produces CV and rMAD plots, violin facet plots (implemented elsewhere)
provide better visualization for these metrics.

Dependencies:
    - pandas: Data manipulation
    - numpy: Numerical computations
    - matplotlib: Basic plotting
    - seaborn: Enhanced plotting
    - statsmodels: LOWESS smoothing for MA plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import os

def load_data(file_path):
    """Load data from CSV file."""
    return pd.read_csv(file_path)

def extract_n_cluster_data(df):
    """Extract N_Cluster columns using regex pattern matching."""
    return df.filter(regex="N_Cluster")

def plot_cv(before_df, after_df, title, cmap='BuGn', title_size=16, label_size=14, tick_size=12):
    """
    Plot Coefficient of Variation comparison before and after transformation.
    
    CV = (standard deviation / mean) represents relative variability.
    Lower values indicate better precision.
    """
    cv_before = before_df.std() / before_df.mean()
    cv_after = after_df.std() / after_df.mean()

    plt.figure(figsize=(7, 6))
    plt.scatter(cv_before, cv_after, alpha=0.5, c=cv_after, cmap=cmap)
    
    # Add diagonal reference line
    max_val = max(cv_before.max(), cv_after.max())
    plt.plot([0, max_val], [0, max_val], 'r--')
    
    plt.title(f'CV Before vs After Transformation - {title}', fontsize=title_size)
    plt.xlabel('CV Before Transformation', fontsize=label_size)
    plt.ylabel('CV After Transformation', fontsize=label_size)
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.grid(True)
    plt.colorbar(label='CV After')
    plt.show()

    return cv_before, cv_after

def plot_ma_transform(before_df, after_df, title, cmap='BuGn', title_size=16, label_size=14, tick_size=12):
    """
    Generate MA plots (intensity-dependent ratio plots) before and after transformation.
    
    MA plots show:
    - A: average of log intensities ((log2(x1) + log2(x2))/2)
    - M: log ratio (log2(x1) - log2(x2))
    
    Used to visualize intensity-dependent bias in the data.
    """
    # Add 1 to avoid log(0)
    before = before_df + 1
    after = after_df + 1

    # Calculate M and A values
    M_before = np.log2(before.iloc[:, 0]) - np.log2(before.iloc[:, 1])
    A_before = 0.5 * (np.log2(before.iloc[:, 0]) + np.log2(before.iloc[:, 1]))
    M_after = np.log2(after.iloc[:, 0]) - np.log2(after.iloc[:, 1])
    A_after = 0.5 * (np.log2(after.iloc[:, 0]) + np.log2(after.iloc[:, 1]))

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
    Calculate Relative Standard Deviation (RSD).
    RSD = (standard deviation / mean) * 100
    """
    mean = np.mean(data, axis=0)
    std_dev = np.std(data, axis=0)
    rsd = (std_dev / mean) * 100
    return rsd

def calculate_rmad(data):
    """
    Calculate Relative Median Absolute Deviation (rMAD).
    rMAD = (median(|x - median(x)|) / median(x)) * 100
    
    More robust to outliers than RSD.
    """
    median = np.median(data, axis=0)
    mad = np.median(np.abs(data - median), axis=0)
    rmad = (mad / median) * 100
    return rmad

def plot_rsd(before_df, after_df, title, before_color='#8ae391', after_color='#0c3b10', 
             title_size=16, label_size=14, tick_size=12):
    """
    Generate boxplot and histogram of Relative Standard Deviation (RSD).
    """
    rsd_before = calculate_rsd(before_df)
    rsd_after = calculate_rsd(after_df)

    results_df = pd.DataFrame({
        'RSD_Before': rsd_before,
        'RSD_After': rsd_after
    })

    # Remove infinite values and NaNs
    results_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    results_df.dropna(inplace=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Boxplot
    axes[0].boxplot([results_df['RSD_Before'], results_df['RSD_After']], 
                    labels=['Before', 'After'], patch_artist=True,
                    boxprops=dict(facecolor=before_color, color=before_color))
    axes[0].set_title(f'Boxplot of Relative Standard Deviation (RSD) - {title}', 
                      fontsize=title_size)
    axes[0].set_ylabel('RSD (%)', fontsize=label_size)
    axes[0].tick_params(axis='both', which='major', labelsize=tick_size)

    # Histogram
    axes[1].hist([results_df['RSD_Before'], results_df['RSD_After']], 
                 bins=30, alpha=0.5, label=['Before', 'After'], 
                 color=[before_color, after_color])
    axes[1].set_title(f'Histogram of Relative Standard Deviation (RSD) - {title}', 
                      fontsize=title_size)
    axes[1].set_ylabel('Frequency', fontsize=label_size)
    axes[1].legend(fontsize=label_size)
    axes[1].tick_params(axis='both', which='major', labelsize=tick_size)

    plt.tight_layout()
    plt.show()

    return rsd_before, rsd_after

def plot_rmad(before_df, after_df, title, before_color='#8ae391', after_color='#0c3b10', 
              title_size=16, label_size=14, tick_size=12):
    """
    Generate boxplot and histogram of Relative Median Absolute Deviation (rMAD).
    """
    rmad_before = calculate_rmad(before_df)
    rmad_after = calculate_rmad(after_df)

    results_df = pd.DataFrame({
        'rMAD_Before': rmad_before,
        'rMAD_After': rmad_after
    })

    # Remove infinite values and NaNs
    results_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    results_df.dropna(inplace=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Boxplot
    axes[0].boxplot([results_df['rMAD_Before'], results_df['rMAD_After']], 
                    labels=['Before', 'After'], patch_artist=True,
                    boxprops=dict(facecolor=before_color, color=before_color))
    axes[0].set_title(f'Boxplot of Relative Median Absolute Deviation (rMAD) - {title}', 
                      fontsize=title_size)
    axes[0].set_ylabel('rMAD (%)', fontsize=label_size)
    axes[0].tick_params(axis='both', which='major', labelsize=tick_size)

    # Histogram
    axes[1].hist([results_df['rMAD_Before'], results_df['rMAD_After']], 
                 bins=30, alpha=0.5, label=['Before', 'After'], 
                 color=[before_color, after_color])
    axes[1].set_title(f'Histogram of Relative Median Absolute Deviation (rMAD) - {title}', 
                      fontsize=title_size)
    axes[1].set_ylabel('Frequency', fontsize=label_size)
    axes[1].legend(fontsize=label_size)
    axes[1].tick_params(axis='both', which='major', labelsize=tick_size)

    plt.tight_layout()
    plt.show()

    return rmad_before, rmad_after

def save_metrics_to_csv(metrics, file_path):
    """Save computed metrics to CSV file."""
    metrics_df = pd.DataFrame.from_dict(metrics, orient='index').transpose()
    metrics_df.to_csv(file_path, index=False)
    print(f"Metrics saved to {file_path}")

def main():
    """
    Main execution function that:
    1. Loads original and transformed datasets
    2. Generates visualization plots
    3. Calculates metrics
    4. Saves results to CSV
    """
    # Define file paths for original and transformed data
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

    # Load and process original data
    original_df = load_data(file_paths['original'])
    original_n_cluster = extract_n_cluster_data(original_df)
    original_median = original_n_cluster.median().mean()

    # Initialize metrics dictionary
    metrics = {
        'Transformation': [],
        'CV_Before': [],
        'CV_After': [],
        'RSD_Before': [],
        'RSD_After': [],
        'rMAD_Before': [],
        'rMAD_After': [],
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
        cv_before, cv_after = plot_cv(original_n_cluster, transformed_n_cluster, title)
        plot_ma_transform(original_n_cluster, transformed_n_cluster, title)
        rsd_before, rsd_after = plot_rsd(original_n_cluster, transformed_n_cluster, title)
        rmad_before, rmad_after = plot_rmad(original_n_cluster, transformed_n_cluster, title)

        # Store metrics
        metrics['Transformation'].append(title)
        metrics['CV_Before'].append(cv_before.mean())
        metrics['CV_After'].append(cv_after.mean())
        metrics['RSD_Before'].append(rsd_before.mean())
        metrics['rMAD_Before'].append(rmad_before.mean())
        metrics['RSD_After'].append(rsd_after.mean())
        metrics['rMAD_After'].append(rmad_after.mean())
        metrics['rMAD_Percentage_of_Median_Before'].append(rmad_before.mean())
        metrics['rMAD_Percentage_of_Median_After'].append(rmad_after.mean())

    # Save results
    save_metrics_to_csv(metrics, os.path.join(os.path.dirname(file_paths['original']), 
                                             'transformation_metrics.csv'))

if __name__ == "__main__":
    main()


