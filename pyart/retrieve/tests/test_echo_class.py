
import pyart
import numpy as np
from numpy.testing.decorators import skipif


# TODO replace with a faster and more appropiate unit test.
@skipif(not pyart.retrieve._F90_EXTENSIONS_AVAILABLE)
def test_steiner_conv_strat():
    grid = pyart.testing.make_target_grid()
    eclass = pyart.retrieve.steiner_conv_strat(grid, refl_field='reflectivity')
    assert np.all(eclass['data'] == 1)
