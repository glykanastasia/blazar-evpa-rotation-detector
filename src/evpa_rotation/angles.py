"""
Angle processing utilities for EVPA rotation analysis.
"""

import numpy as np


def adjust_angles(angles, errors):
    """
    Improved EVPA angle adjustment accounting for measurement errors.
    """
    if len(angles) == 0:
        return []
    
    adjusted = [angles[0]]
    
    for i in range(1, len(angles)):
        prev = adjusted[-1]
        current = angles[i]
        
        # Calculate error-weighted difference
        error_term = np.sqrt(errors[i]**2 + errors[i-1]**2)
        raw_diff = abs(current - prev)
        
        # Implement the error-weighted condition mentioned in comments
        weighted_diff = raw_diff - error_term
        
        # Only apply shifts if the change is significant beyond measurement errors
        if weighted_diff <= 90:
            adjusted.append(current)
        else:
            # Search for best candidate considering error-weighted differences
            best_candidate = current
            best_weighted_diff = float('inf')
            
            for n in range(-5, 6):
                candidate = current + 180 * n
                candidate_diff = abs(candidate - prev)
                candidate_weighted_diff = candidate_diff - error_term
                
                if abs(candidate_weighted_diff) < best_weighted_diff:
                    best_weighted_diff = abs(candidate_weighted_diff)
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
