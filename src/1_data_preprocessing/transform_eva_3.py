"""
Data Transformation Evaluation Script for Metabolomics Analysis
------------------------------------------------------------

This script creates facet grid visualisations to compare original and transformed
metabolomics data using two key metrics:
1. Coefficient of Variation (CV): Measures relative variability
2. Relative Median Absolute Deviation (rMAD): Robust measure of dispersion

The script generates violin plots arranged in a facet grid to visualise the 
distribution of these metrics before and after transformation, enabling direct
comparison of transformation effects on data variability.

Input:
- Original data CSV file
- Transformed data CSV file
- Both files should contain N_Cluster columns for analysis

Output:
- Facet grid plots showing:
  * CV distribution comparison
  * rMAD distribution comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def extract_n_cluster_data(df):
    """
    Extract columns containing N_Cluster in their names.
    
    Args:
        df (pd.DataFrame): Input DataFrame
        
    Returns:
        pd.DataFrame: DataFrame containing only N_Cluster columns
    """
    return df.filter(regex="N_Cluster")

def calculate_cv(series):
    """
    Calculate Coefficient of Variation (CV) for a data series.
    
    CV = (standard deviation / mean) × 100%
    
    Args:
        series (pd.Series): Input data series
        
    Returns:
        float: CV value as percentage
    """
    return np.std(series) / np.mean(series)

def calculate_rmad(series):
    """
    Calculate relative Median Absolute Deviation (rMAD) for a data series.
    
    rMAD = median(|x - median(x)|) / median(x)
    
    Args:
        series (pd.Series): Input data series
        
    Returns:
        float: rMAD value as ratio
    """
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    return mad / median

def calculate_metrics(df):
    """
    Calculate CV and rMAD for all columns in DataFrame.
    
    Args:
        df (pd.DataFrame): Input DataFrame
        
    Returns:
        tuple: (CV values, rMAD values) for all columns
    """
    cv_values = df.apply(calculate_cv)
    rmad_values = df.apply(calculate_rmad)
    return cv_values, rmad_values

def load_data(original_path, transformed_path):
    """
    Load original and transformed datasets.
    
    Args:
        original_path (str): Path to original data CSV
        transformed_path (str): Path to transformed data CSV
        
    Returns:
        tuple: (Original DataFrame, Transformed DataFrame)
    """
    original_df = pd.read_csv(original_path)
    transformed_df = pd.read_csv(transformed_path)
    return original_df, transformed_df

def prepare_data_for_plotting(original_df, transformed_df):
    """
    Prepare data for visualisation by calculating metrics and formatting for seaborn.
    
    Args:
        original_df (pd.DataFrame): Original data
        transformed_df (pd.DataFrame): Transformed data
        
    Returns:
        tuple: (CV DataFrame, rMAD DataFrame) ready for plotting
    """
    orig_n_cluster = extract_n_cluster_data(original_df)
    orig_cv, orig_rmad = calculate_metrics(orig_n_cluster)
    
    trans_n_cluster = extract_n_cluster_data(transformed_df)
    trans_cv, trans_rmad = calculate_metrics(trans_n_cluster)
    
    # Format CV data
    cv_df = pd.DataFrame({
        'Original': orig_cv,
        'Transformed': trans_cv
    }).reset_index().melt(id_vars='index', var_name='Type', value_name='CV')
    cv_df['Transformation'] = 'asinh'  # Set transformation type
    
    # Format rMAD data
    rmad_df = pd.DataFrame({
        'Original': orig_rmad,
        'Transformed': trans_rmad
    }).reset_index().melt(id_vars='index', var_name='Type', value_name='rMAD')
    rmad_df['Transformation'] = 'asinh'  # Set transformation type
    
    return cv_df, rmad_df

def plot_facet_grid(df, metric, font_size):
    """
    Create facet grid plot comparing original and transformed distributions.
    
    Args:
        df (pd.DataFrame): Data prepared for plotting
        metric (str): Metric name ('CV' or 'rMAD')
        font_size (dict): Font sizes for different plot elements
    """
    g = sns.FacetGrid(df, col="Transformation", sharex=False, sharey=False, height=4)
    g.map_dataframe(sns.violinplot, x='Type', y=metric, palette='BuGn', split=True)
    g.set_axis_labels("Data Type", metric)
    g.set_titles(col_template="{col_name}", size=font_size['title'])
    
    # Customise plot appearance
    for ax in g.axes.flat:
        ax.set_xlabel('Data Type', fontsize=font_size['xlabel'])
        ax.set_ylabel(metric, fontsize=font_size['ylabel'])
        ax.tick_params(axis='both', which='major', labelsize=font_size['ticks'])
    
    g.fig.subplots_adjust(top=0.85, hspace=0.4)
    g.fig.suptitle(f'Violin Plot of {metric}', fontsize=font_size['title'])
    plt.tight_layout()
    plt.show()

def main():
    # File paths
    original_path = r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if.csv"
    transformed_path = r"C:\Users\ms\Desktop\data_chem\data\transformation\n_l_if_asinh.csv"
    
    # Define font sizes for consistent plotting
    font_size = {
        'title': 18,
        'xlabel': 16,
        'ylabel': 16,
        'ticks': 14,
        'legend': 14
    }
    
    # Load and process data
    original_df, transformed_df = load_data(original_path, transformed_path)
    cv_df, rmad_df = prepare_data_for_plotting(original_df, transformed_df)
    
    # Generate visualisation
    plot_facet_grid(cv_df, 'CV', font_size)
    plot_facet_grid(rmad_df, 'rMAD', font_size)

if __name__ == "__main__":
    main()