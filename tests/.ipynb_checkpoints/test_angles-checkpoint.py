import numpy as np
from evpa_rotation.angles import adjust_angles

def test_empty():
    assert adjust_angles([], []) == []

def test_simple_shift():
    angles = [10, 200]
    errs = [1, 1]
    adj = adjust_angles(angles, errs)
    # Expect 200→20 to minimize delta
    assert abs(adj[1] - 20) < 1e-6
