"""
Metabolomics Data Quality Assessment Script
-----------------------------------------

This script performs initial quality assessment of metabolomics data by calculating two key metrics:
1. Coefficient of Variation (CV): Measures relative variability as (standard deviation / mean) × 100%
2. Relative Median Absolute Deviation (rMAD): A robust measure of variability that's less sensitive 
   to outliers than CV, calculated as (median absolute deviation / median) × 100%

These metrics help identify metabolites with high variability that might need additional
preprocessing or filtering before downstream analysis.

Thresholds:
- CV > 30%: Indicates high relative variability
- rMAD > 30%: Indicates high robust relative variability

Input:
- CSV files containing metabolomics data with N_Cluster or P_Cluster columns
- Non-cluster columns are ignored in the analysis

Output:
- Printed summary of the percentage of variables exceeding quality thresholds
"""

import pandas as pd
import numpy as np
import os
from typing import Tuple, List

def calculate_cv(data: pd.DataFrame) -> pd.Series:
    """
    Calculate Coefficient of Variation (CV) for each column.
    
    CV = (standard deviation / mean) × 100%
    
    Args:
        data: DataFrame containing metabolite measurements
        
    Returns:
        Series containing CV values for each column
    """
    cv = (data.std() / data.mean()) * 100
    return cv

def calculate_rmad(data: pd.DataFrame) -> pd.Series:
    """
    Calculate Relative Median Absolute Deviation (rMAD) for each column.
    
    rMAD = (median(|x - median(x)|) / median(x)) × 100%
    
    This is a more robust measure of variability than CV as it's less sensitive to outliers.
    
    Args:
        data: DataFrame containing metabolite measurements
        
    Returns:
        Series containing rMAD values for each column
    """
    median = data.median()
    mad = (np.abs(data - median)).median()
    rmad = (mad / median) * 100
    return rmad

def get_cluster_columns(data: pd.DataFrame) -> List[str]:
    """
    Identify metabolite cluster columns in the dataset.
    
    Args:
        data: DataFrame containing metabolomics data
        
    Returns:
        List of column names that start with either N_Cluster or P_Cluster
    """
    return [col for col in data.columns if col.startswith(('N_Cluster', 'P_Cluster'))]

def analyse_variability(data: pd.DataFrame, 
                       cluster_columns: List[str], 
                       cv_threshold: float, 
                       rmad_threshold: float) -> Tuple[float, float]:
    """
    Analyse variability in metabolite data using CV and rMAD metrics.
    
    Args:
        data: DataFrame containing metabolomics data
        cluster_columns: List of metabolite columns to analyse
        cv_threshold: Threshold for CV values (in percentage)
        rmad_threshold: Threshold for rMAD values (in percentage)
        
    Returns:
        Tuple containing:
        - Percentage of variables exceeding CV threshold
        - Percentage of variables exceeding rMAD threshold
    """
    cluster_data = data[cluster_columns]
    
    # Calculate metrics
    cv_values = calculate_cv(cluster_data)
    rmad_values = calculate_rmad(cluster_data)
    
    # Calculate percentage exceeding thresholds
    cv_exceeding = (cv_values > cv_threshold).mean() * 100
    rmad_exceeding = (rmad_values > rmad_threshold).mean() * 100
    
    return cv_exceeding, rmad_exceeding

def main():
    # Define quality thresholds
    CV_THRESHOLD = 30  # 30% threshold for Coefficient of Variation
    RMAD_THRESHOLD = 30  # 30% threshold for relative Median Absolute Deviation
    
    # List of files to analyse
    file_paths = [
        r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\n_l_if_asinh.csv',
        r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\n_r_if_asinh.csv',
        r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\p_l_if_asinh.csv',
        r'C:\Users\ms\Desktop\data_chem\data\CV_rMAD\p_r_if_asinh.csv'
    ]
    
    # Process each file
    for file_path in file_paths:
        # Load data
        data = pd.read_csv(file_path)
        file_name = os.path.basename(file_path)
        
        # Get metabolite columns
        cluster_columns = get_cluster_columns(data)
        
        # Analyse variability
        cv_exceeding, rmad_exceeding = analyse_variability(
            data, cluster_columns, CV_THRESHOLD, RMAD_THRESHOLD
        )
        
        # Print results
        print(f"\nResults for {file_name}:")
        print(f"Percentage of metabolites with CV > {CV_THRESHOLD}%: "
              f"{cv_exceeding:.2f}%")
        print(f"Percentage of metabolites with rMAD > {RMAD_THRESHOLD}%: "
              f"{rmad_exceeding:.2f}%")

if __name__ == "__main__":
    main()