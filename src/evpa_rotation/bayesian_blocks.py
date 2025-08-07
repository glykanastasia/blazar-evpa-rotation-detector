"""
Bayesian Blocks analysis for EVPA time series.
"""

import numpy as np
import pandas as pd
import astropy.stats
from scipy.signal import argrelextrema


def compute_bayesian_blocks(subset, p0=0.001):
    """
    Computes the Bayesian Blocks for a given source, after preprocessing the data.
    
    Parameters:
    -----------
    subset : pd.DataFrame
        Data subset containing MJD, adjusted_EVPA, and err_EVPA[deg] columns
    p0 : float, default=0.001
        Prior probability for change points
        
    Returns:
    --------
    tuple
        (edges, bin_means, bin_errors) where:
        - edges: array of bin edges
        - bin_means: list of mean EVPA values in each bin
        - bin_errors: list of mean errors in each bin
    """
    # Remove duplicate MJDs: average EVPA and errors
    subset = subset.groupby('MJD', as_index=False).agg({
        'adjusted_EVPA': 'mean',
        'err_EVPA[deg]': 'mean'
    }).sort_values('MJD').copy()

    # Extract arrays
    mjd = subset['MJD'].to_numpy()
    adjusted_evpa = subset['adjusted_EVPA'].to_numpy()
    error = subset['err_EVPA[deg]'].to_numpy()

    # Compute Bayesian Blocks
    edges = astropy.stats.bayesian_blocks(
        mjd, adjusted_evpa, sigma=error, fitness='measures', p0=p0
    )

    # Adjust first and last edges to match data limits
    if len(edges) > 1:
        edges[0] = mjd.min()
        edges[-1] = mjd.max()

    # Compute bin means and errors
    bin_means = []
    bin_errors = []
    
    for i in range(len(edges) - 1):
        if i < len(edges) - 2:
            mask = (mjd >= edges[i]) & (mjd < edges[i+1])
        else:
            mask = (mjd >= edges[i]) & (mjd <= edges[i+1])
        
        if np.any(mask):
            bin_means.append(np.median(adjusted_evpa[mask]))
            bin_errors.append(np.mean(error[mask]))
        else:
            bin_means.append(np.nan)
            bin_errors.append(np.nan)

    return edges, bin_means, bin_errors


def compute_local_extrema(bin_means, bin_centers):
    """
    Finds local maxima and minima in Bayesian Block bins and ensures 
    first/last points are included if they are extrema.
    
    Parameters:
    -----------
    bin_means : array-like
        Mean EVPA values in each bin
    bin_centers : array-like
        Center time points of each bin
        
    Returns:
    --------
    tuple
        (maxima_dict, minima_dict) containing MJD and EVPA arrays for extrema
    """
    bin_means = np.array(bin_means)
    bin_centers = np.array(bin_centers)
    
    # Identify local maxima and minima
    maxima_indices = argrelextrema(bin_means, np.greater)[0].tolist()
    minima_indices = argrelextrema(bin_means, np.less)[0].tolist()

    # Check first and last bins
    if len(bin_means) > 1:
        if bin_means[0] > bin_means[1]:  # First point is a maximum
            maxima_indices.insert(0, 0)
        elif bin_means[0] < bin_means[1]:  # First point is a minimum
            minima_indices.insert(0, 0)

        if bin_means[-1] > bin_means[-2]:  # Last point is a maximum
            maxima_indices.append(len(bin_means) - 1)
        elif bin_means[-1] < bin_means[-2]:  # Last point is a minimum
            minima_indices.append(len(bin_means) - 1)

    # Create dictionaries
    maxima_dict = {
        'MJD': bin_centers[maxima_indices], 
        'EVPA': bin_means[maxima_indices]
    }
    minima_dict = {
        'MJD': bin_centers[minima_indices], 
        'EVPA': bin_means[minima_indices]
    }
    
    return maxima_dict, minima_dict


def analyze_extrema_pairs(maxima_dict, minima_dict, min_amplitude=0):
    """
    Analyze pairs of adjacent extrema to identify potential rotation events.
    
    Parameters:
    -----------
    maxima_dict : dict
        Dictionary with 'MJD' and 'EVPA' keys for maxima
    minima_dict : dict
        Dictionary with 'MJD' and 'EVPA' keys for minima
    min_amplitude : float, default=0
        Minimum EVPA amplitude to consider as rotation
        
    Returns:
    --------
    list
        List of tuples (t_start, evpa_start, t_end, evpa_end, amplitude, type)
        for potential rotation events
    """
    # Combine extrema
    extrema = []
    for mjd, evpa in zip(maxima_dict['MJD'], maxima_dict['EVPA']):
        extrema.append((mjd, evpa, 'max'))
    for mjd, evpa in zip(minima_dict['MJD'], minima_dict['EVPA']):
        extrema.append((mjd, evpa, 'min'))
    
    # Sort by time
    extrema.sort(key=lambda x: x[0])
    
    rotation_candidates = []
    
    # Find adjacent extrema pairs
    for i in range(len(extrema) - 1):
        t_start, evpa_start, type_start = extrema[i]
        t_end, evpa_end, type_end = extrema[i + 1]
        
        # Only consider max→min or min→max transitions
        if not ((type_start == 'max' and type_end == 'min') or
                (type_start == 'min' and type_end == 'max')):
            continue
        
        amplitude = abs(evpa_end - evpa_start)
        if amplitude >= min_amplitude:
            rotation_type = f"{type_start}→{type_end}"
            rotation_candidates.append((
                t_start, evpa_start, t_end, evpa_end, amplitude, rotation_type
            ))
    
    return rotation_candidates


def get_bin_boundaries_for_interval(edges, t_start, t_end):
    """
    Find the appropriate bin boundaries that contain the given time interval.
    
    Parameters:
    -----------
    edges : array
        Bin edges from Bayesian Blocks
    t_start : float
        Start time of interval
    t_end : float
        End time of interval
        
    Returns:
    --------
    tuple
        (start_edge, end_edge) - actual bin boundaries
    """
    # Find the bin edges that contain the extrema
    start_edge = edges[0]  # default
    end_edge = edges[-1]   # default
    
    # Find start edge
    for k in range(len(edges) - 1):
        if edges[k] <= t_start < edges[k + 1]:
            start_edge = edges[k]
            break
    
    # Find end edge
    for k in range(len(edges) - 1):
        if edges[k] <= t_end < edges[k + 1]:
            end_edge = edges[k + 1]
            break
    
    return start_edge, end_edge
