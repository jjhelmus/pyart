""" Unit tests for Py-ART's retrieve/echo_classification.py module. """

import numpy as np

import pyart


def test_echo_classification_steiner():
    grid = pyart.testing.make_target_grid()
    edic = pyart.retrieve.echo_classification_steiner(
        grid, 'reflectivity', 100.0, intense=30.0)
    edata = edic['data']
    assert edata.shape == (400, 320)
    assert edata[0, 0] == -1.
    assert edata[200, 160] == 1.


def test_echo_classification_steiner_masked():
    grid = pyart.testing.make_target_grid()
    grid.fields['reflectivity']['data'] = np.ma.masked_equal(
        grid.fields['reflectivity']['data'], 0)
    edic = pyart.retrieve.echo_classification_steiner(
        grid, 'reflectivity', 100.0, intense=30.0)
    edata = edic['data']
    assert edata.shape == (400, 320)
    assert np.ma.is_masked(edata[0, 0])
    assert not np.ma.is_masked(edata[200, 160])
    assert edata[100, 160] == -1.
    assert edata[200, 160] == 1.
