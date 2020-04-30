import kea.hist
import numpy as np
import pytest


# Test the creation of the histogram
def test_histogram_creation():
    x = np.linspace(0,100,101)
    h = kea.hist.histogram(x, np.full(100,1))
    assert (h.values == np.full(100,1)).all()
    assert (h.edges == x).all()

# Test the manipulation of the histogram
def test_histogram_manipulation():
    x = np.linspace(0,100,101)
    h = kea.hist.histogram(x, np.full(100,1))
    assert (h == np.full(100,1)).all()
    assert h[0] == 1
    h[1] = 5
    assert h[1] == 5
    h[1:10] = 1
    assert (h == np.full(100,1)).all()

# Test the multiplication
def test_histogram_operators():
    x = np.linspace(0,100,101)
    h = kea.hist.histogram(x, np.full(100,1))
    assert (h * 10 == np.full(100,10)).all()
    h *= 10
    assert (h == np.full(100,10)).all()
    h /= 10
    assert (h == np.full(100,1)).all()

# Test get_bin
def test_histogram_bin():
    x = np.linspace(0,100,1001)
    h = kea.hist.histogram(x, np.full(1000,1))
    assert h.get_bin(0.0) == 0             # lowest bound
    assert h.get_bin(1.121e-03) == 0       # check bins
    assert h.get_bin(0.1) == 1       #
    assert h.get_bin(0.15) == 1    # check if edge in correct bin
    assert h.get_bin(100) == 999 # last bin moved outside
    for num, val in enumerate(x[:-1]):
        assert h.get_bin(val) == num
    with pytest.raises(Exception, match="x outside of range"):
        h.get_bin(-1)
        h.get_bin(101)



# Test integration
def test_histogram_integration():
    x = np.linspace(0,100,101)
    h = kea.hist.histogram(x, np.full(100,1))
    assert h.integral(0,1) == 1
    assert h.integral(0,0.5) == 0.5
    assert h.integral(0,1.5) == 1.5
    assert h.integral(0,10) == 10
    assert h.integral(0,100) == 100
    h = kea.hist.histogram(np.linspace(0,100,101), np.full(100,2))
    assert h.integral(0,1) == 2

def test_hitogram_bin_center():
    x = np.linspace(0,100,101)
    h = kea.hist.histogram(x, np.full(100,1))
    assert h.get_bin_center(0) == 0.5
    for i, _ in enumerate(h):
        assert h.get_bin_center(i) == i+0.5

def test_histogram_copy():
    x = np.linspace(0,100,101)
    h = kea.hist.histogram(x, np.full(100,1))
    h2 = h.copy()
    assert (h.values == h2.values).all()
    h[10] = 10
    assert ((h.values == h2.values).all() == False)
    h2[10] = 10
    assert (h.values == h2.values).all()


def test_histogram_BPASS():
    h = kea.hist.histogram_BPASS(np.full(51,1))
    assert h.integral(0,10)
