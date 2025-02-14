"""
Missing At Random (MAR) Analysis Script
-------------------------------------

This script tests whether missing data follows a Missing At Random (MAR) pattern
using logistic regression. MAR occurs when the probability of missingness
depends on observed variables but not on the missing values themselves.

Methodology:
1. For each metabolite, creates a binary indicator (0/1) for missingness
2. Uses logistic regression to model missingness based on experimental factors:
   - Genotype
   - Treatment (TMT)
   - Day
   - Replication
3. Significant associations between these factors and missingness suggest MAR pattern

Input:
- CSV file containing metabolomics data with N_Cluster columns and experimental factors
- Factors must include: Genotype, TMT, Day, Rep

Output:
- CSV file with regression coefficients for each metabolite-factor combination
- Coefficients indicate strength and direction of association with missingness
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

def load_data(file_path):
    """
    Load and validate the input data file.
    
    Args:
        file_path (str): Path to the input CSV file
        
    Returns:
        pd.DataFrame: Loaded data
    """
    return pd.read_csv(file_path)

def create_missingness_indicator(data, metabolite_column):
    """
    Create binary indicator for missing values in specified metabolite.
    
    Args:
        data (pd.DataFrame): Input dataframe
        metabolite_column (str): Name of metabolite column
        
    Returns:
        pd.Series: Binary indicator (1 for missing, 0 for present)
    """
    return data[metabolite_column].isnull().astype(int)

def prepare_predictors(data, predictor_columns):
    """
    Prepare predictor variables with one-hot encoding.
    
    Args:
        data (pd.DataFrame): Input dataframe
        predictor_columns (list): List of predictor column names
        
    Returns:
        pd.DataFrame: One-hot encoded predictor variables
    """
    return pd.get_dummies(data[predictor_columns], drop_first=True)

def fit_logistic_model(X_train, X_test, y_train, y_test):
    """
    Fit logistic regression model with standardized features.
    
    Args:
        X_train, X_test: Training and test predictor variables
        y_train, y_test: Training and test response variables
        
    Returns:
        tuple: Fitted model and scaler object
    """
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    model = LogisticRegression()
    model.fit(X_train_scaled, y_train)
    
    return model, scaler

def extract_model_results(model, column):
    """
    Extract regression coefficients and create results DataFrame.
    
    Args:
        model: Fitted logistic regression model
        column (str): Name of metabolite column
        
    Returns:
        pd.DataFrame: Results with coefficients for each predictor
    """
    return pd.DataFrame({
        'Metabolite': [column],
        'Intercept': [model.intercept_[0]],
        'Genotype_G2': [model.coef_[0][0]] if len(model.coef_[0]) > 0 else [None],
        'TMT_1': [model.coef_[0][1]] if len(model.coef_[0]) > 1 else [None],
        'Day_2': [model.coef_[0][2]] if len(model.coef_[0]) > 2 else [None],
        'Rep': [model.coef_[0][3]] if len(model.coef_[0]) > 3 else [None],
    })

def create_null_results(column):
    """
    Create DataFrame with null results for cases with insufficient variation.
    
    Args:
        column (str): Name of metabolite column
        
    Returns:
        pd.DataFrame: Results DataFrame with null values
    """
    return pd.DataFrame({
        'Metabolite': [column],
        'Intercept': [None],
        'Genotype_G2': [None],
        'TMT_1': [None],
        'Day_2': [None],
        'Rep': [None],
    })

def main():
    # Input and output file paths
    input_path = r'C:\Users\ms\Desktop\data_chem\data\old_1\n_column_data_r_clean.csv'
    output_path = r'C:\Users\ms\Desktop\data_chem\result\n_column_data_r_clean_results.csv'
    
    # Load data
    data = load_data(input_path)
    
    # Initialize results DataFrame
    results = pd.DataFrame()
    
    # Analyze each metabolite column
    for column in data.columns:
        if not column.startswith('N_Cluster_'):
            continue
            
        # Create missingness indicator
        data['missing_' + column] = create_missingness_indicator(data, column)
        
        # Prepare predictor variables
        predictors = ['Genotype', 'TMT', 'Day', 'Rep']
        X = prepare_predictors(data, predictors)
        y = data['missing_' + column]
        
        # Skip if no variation in missingness
        if len(y.unique()) <= 1:
            result = create_null_results(column)
        else:
            try:
                # Split data and fit model
                X_train, X_test, y_train, y_test = train_test_split(
                    X, y, test_size=0.3, random_state=42
                )
                model, _ = fit_logistic_model(X_train, X_test, y_train, y_test)
                result = extract_model_results(model, column)
            except Exception as e:
                print(f"Error fitting model for {column}: {e}")
                continue
        
        # Add non-null results to output
        if not result.isnull().all().all():
            results = pd.concat([results, result], ignore_index=True)
    
    # Save results
    results.to_csv(output_path, index=False)
    print(f"Analysis complete. Results saved to: {output_path}")

if __name__ == "__main__":
    main()