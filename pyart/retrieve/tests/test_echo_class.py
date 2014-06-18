
import pyart
import numpy as np

def test_steiner_conv_strat():
    grid = pyart.testing.make_target_grid()
    eclass = pyart.retrieve.steiner_conv_strat(grid, refl_field='reflectivity')
    assert np.all(eclass['data'] == 1)
