import numpy as np  # type: ignore[import]
from optproj import l2_norm

def test_l2_norm_basic():
    a = np.array([3.0, 4.0], dtype=np.float64)
    assert abs(l2_norm(a) - 5.0) < 1e-12

def test_l2_norm_zero():
    a = np.zeros(10, dtype=np.float64)
    assert abs(l2_norm(a)) < 1e-12
