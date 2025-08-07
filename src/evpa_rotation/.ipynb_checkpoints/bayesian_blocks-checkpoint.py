import numpy as np
import pandas as pd
import astropy.stats
from scipy.signal import argrelextrema
from typing import Tuple, List, Dict

def compute_bayesian_blocks(
    subset: pd.DataFrame,
    p0: float = 0.001
) -> Tuple[np.ndarray, List[float], List[float]]:
    """
    Compute Bayesian Blocks on adjusted EVPA vs. MJD.
    Returns (edges, bin_means, bin_errors).
    """
    df = subset.groupby('MJD', as_index=False).agg({
        'adjusted_EVPA': 'mean',
        'err_EVPA[deg]': 'mean'
    }).sort_values('MJD')
    mjd = df['MJD'].to_numpy()
    evpa = df['adjusted_EVPA'].to_numpy()
    err = df['err_EVPA[deg]'].to_numpy()

    edges = astropy.stats.bayesian_blocks(
        mjd, evpa, sigma=err, fitness='measures', p0=p0
    )
    if len(edges) > 1:
        edges[0], edges[-1] = mjd.min(), mjd.max()

    means, errs = [], []
    for i in range(len(edges) - 1):
        if i < len(edges) - 2:
            mask = (mjd >= edges[i]) & (mjd < edges[i+1])
        else:
            mask = (mjd >= edges[i]) & (mjd <= edges[i+1])
        if np.any(mask):
            means.append(np.median(evpa[mask]))
            errs.append(np.mean(err[mask]))
        else:
            means.append(np.nan)
            errs.append(np.nan)

    return edges, means, errs

def compute_local_extrema(
    bin_means: List[float],
    bin_centers: List[float]
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Identify local maxima and minima in bin_means.
    Returns two dicts: maxima and minima with keys 'MJD' and 'EVPA'.
    """
    arr = np.array(bin_means)
    centers = np.array(bin_centers)

    max_idx = argrelextrema(arr, np.greater)[0].tolist()
    min_idx = argrelextrema(arr, np.less)[0].tolist()

    if len(arr) > 1:
        if arr[0] > arr[1]:
            max_idx.insert(0, 0)
        elif arr[0] < arr[1]:
            min_idx.insert(0, 0)
        if arr[-1] > arr[-2]:
            max_idx.append(len(arr)-1)
        elif arr[-1] < arr[-2]:
            min_idx.append(len(arr)-1)

    maxima = {'MJD': centers[max_idx], 'EVPA': arr[max_idx]}
    minima = {'MJD': centers[min_idx], 'EVPA': arr[min_idx]}
    return maxima, minima
