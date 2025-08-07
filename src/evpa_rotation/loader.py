"""
Data loading and preprocessing module for EVPA rotation analysis.
"""

import pandas as pd
import numpy as np


def load_data(file_path):
    """
    Load and preprocess monitoring data from CSV file.
    
    Parameters:
    -----------
    file_path : str
        Path to the CSV file containing monitoring data
        
    Returns:
    --------
    pd.DataFrame
        Preprocessed data with MJD and time differences calculated,
        filtered for PD/error_PD >= 3
    """
    data = pd.read_csv(file_path)
    
    # Convert Julian date to Modified Julian Date
    data['MJD'] = data['Julian_date'] - 2400000.5
    
    # Calculate time differences between consecutive observations
    data['Time_Diff'] = data['MJD'].diff()
    
    # Filter data where PD/error_PD ratio is at least 3 (significant polarization detection)
    data = data[data['PD[%]'] / data['err_PD[%]'] >= 3]
    
    print("Data loaded and preprocessed.")
    print(f"Total observations: {len(data)}")
    print(f"Sources: {data['J2000_name'].nunique()}")
    print(f"Date range: {data['MJD'].min():.1f} - {data['MJD'].max():.1f} MJD")
    
    return data


def segment_data_by_gaps(df_source, gap_threshold=30):
    """
    Segment source data into continuous observation periods based on time gaps.
    
    Parameters:
    -----------
    df_source : pd.DataFrame
        Data for a single source, sorted by MJD
    gap_threshold : float, default=30
        Maximum gap in days to consider observations as continuous
        
    Returns:
    --------
    pd.DataFrame
        Input dataframe with added 'segment' column indicating observation periods
    """
    df_source = df_source.copy()
    df_source['MJD_diff'] = df_source['MJD'].diff()
    df_source['segment'] = (df_source['MJD_diff'] > gap_threshold).cumsum()
    
    return df_source


def prepare_source_data(data, source_name):
    """
    Prepare data for a specific source with segmentation and EVPA adjustment.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Full dataset
    source_name : str
        Name of the source to analyze
        
    Returns:
    --------
    pd.DataFrame
        Prepared source data with segments and adjusted EVPA
    """
    from .angles import adjust_angles
    
    df_source = data[data["J2000_name"] == source_name].sort_values("MJD").copy()
    
    # Segment data based on observation gaps
    df_source = segment_data_by_gaps(df_source)
    
    # Adjust EVPA angles to resolve 180° ambiguity
    df_source['adjusted_EVPA'] = adjust_angles(
        df_source['EVPA[deg]'].tolist(),
        df_source['err_EVPA[deg]'].tolist()
    )
    
    return df_source
