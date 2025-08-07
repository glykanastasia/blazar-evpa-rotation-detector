def adjust_angles(angles, errors):
    """
    Adjust the EVPA angles to resolve the 180° ambiguity based on minimal changes
    between consecutive measurements and accounting for measurement errors.
    
    Parameters:
        angles (list or array): EVPA measurements in degrees.
        errors (list or array): Corresponding errors for the EVPA measurements.
    
    Returns:
        adjusted (list): The adjusted EVPA measurements.
        
    Procedure:
      - Compute the error-weighted EVPA variation between successive angles.
      - If the variation Δθ = |θₙ₊₁ − θₙ|  ≤ 90°,
        no shift is applied.
      - If Δθ > 90°, try candidate shifts of multiples of 180° (n = -5,...,5)
        and choose the candidate that minimizes Δθ.
    """
    if len(angles) == 0:
        return []
    
    adjusted = [angles[0]]  # Use the first measurement as is

    for i in range(1, len(angles)):
        prev = adjusted[-1]
        current = angles[i]
        # Compute the error-weighted difference for the unshifted value
        raw_diff = abs(current - prev) 
        
        # If the raw difference is smooth, leave the value unchanged
        if raw_diff <= 90:
            adjusted.append(current)
        else:
            # Otherwise, search for a candidate that minimizes the weighted difference.
            best_candidate = current  # default to no shift
            best_diff = float('inf')
            # Try shifts over a reasonable range (here from -5 to +5 multiples of 180°)
            for n in range(-5, 6):
                candidate = current + 180 * n
                candidate_diff = abs(candidate - prev)
                # We look at the absolute value of the weighted difference
                if abs(candidate_diff) < best_diff:
                    best_diff = abs(candidate_diff)
                    best_candidate = candidate
            adjusted.append(best_candidate)
    
    return adjusted