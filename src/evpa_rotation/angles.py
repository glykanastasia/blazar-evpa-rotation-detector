"""
Angle processing utilities for EVPA rotation analysis.
"""

import numpy as np


def adjust_angles(angles, errors):
    """
    Adjust the EVPA angles to resolve the 180° ambiguity based on minimal changes
    between consecutive measurements and accounting for measurement errors.
    
    Parameters:
    -----------
    angles : list or array
        EVPA measurements in degrees
    errors : list or array
        Corresponding errors for the EVPA measurements
    
    Returns:
    --------
    list
        The adjusted EVPA measurements
        
    Procedure:
    ----------
    - Compute the error-weighted EVPA variation between successive angles
    - If the variation Δθ = |θₙ₊₁ − θₙ| - sqrt(σ(θₙ₊₁)² + σ(θₙ)²) ≤ 90°,
      no shift is applied
    - If Δθ > 90°, try candidate shifts of multiples of 180° (n = -5,...,5)
      and choose the candidate that minimizes Δθ
    """
    if len(angles) == 0:
        return []
    
    adjusted = [angles[0]]  # Use the first measurement as is

    for i in range(1, len(angles)):
        prev = adjusted[-1]
        current = angles[i]
        
        # Calculate the error term for the current and previous measurement
        error_term = np.sqrt(errors[i]**2 + errors[i-1]**2)
        
        # Compute the error-weighted difference for the unshifted value
        raw_diff = abs(current - prev)
        
        # If the raw difference is smooth, leave the value unchanged
        if raw_diff <= 90:
            adjusted.append(current)
        else:
            # Search for a candidate that minimizes the weighted difference
            best_candidate = current  # default to no shift
            best_diff = float('inf')
            
            # Try shifts over a reasonable range (from -5 to +5 multiples of 180°)
            for n in range(-5, 6):
                candidate = current + 180 * n
                candidate_diff = abs(candidate - prev)
                
                # Look at the absolute value of the weighted difference
                if abs(candidate_diff) < best_diff:
                    best_diff = abs(candidate_diff)
                    best_candidate = candidate
            
            adjusted.append(best_candidate)
    
    return adjusted


def compute_angle_statistics(angles, errors):
    """
    Compute basic statistics for angle measurements.
    
    Parameters:
    -----------
    angles : array-like
        Angle measurements in degrees
    errors : array-like
        Corresponding measurement errors
        
    Returns:
    --------
    dict
        Dictionary containing mean, std, min, max, and amplitude
    """
    angles = np.array(angles)
    errors = np.array(errors)
    
    stats = {
        'mean': np.mean(angles),
        'std': np.std(angles),
        'min': np.min(angles),
        'max': np.max(angles),
        'amplitude': np.max(angles) - np.min(angles),
        'mean_error': np.mean(errors),
        'n_points': len(angles)
    }
    
    return stats


def circular_mean(angles, degrees=True):
    """
    Compute circular mean of angles.
    
    Parameters:
    -----------
    angles : array-like
        Angles in degrees (default) or radians
    degrees : bool, default=True
        Whether input angles are in degrees
        
    Returns:
    --------
    float
        Circular mean in same units as input
    """
    if degrees:
        angles = np.deg2rad(angles)
    
    # Compute circular mean
    mean_angle = np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles)))
    
    if degrees:
        mean_angle = np.rad2deg(mean_angle)
    
    return mean_angle


def angular_distance(angle1, angle2, degrees=True):
    """
    Compute the minimum angular distance between two angles.
    
    Parameters:
    -----------
    angle1, angle2 : float
        Angles in degrees (default) or radians
    degrees : bool, default=True
        Whether input angles are in degrees
        
    Returns:
    --------
    float
        Minimum angular distance in same units as input
    """
    if degrees:
        angle1 = np.deg2rad(angle1)
        angle2 = np.deg2rad(angle2)
    
    # Compute angular distance
    diff = np.abs(angle1 - angle2)
    distance = np.minimum(diff, 2*np.pi - diff)
    
    if degrees:
        distance = np.rad2deg(distance)
    
    return distance
