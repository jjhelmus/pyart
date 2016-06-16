""" Unit Tests for Py-ART's core/sweep.py module. """

import pickle
import sys
# we need a class which excepts str for writing in Python 2 and 3
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import inspect

import numpy as np
from numpy.testing import assert_raises, assert_allclose, assert_almost_equal
import pyart
from pyart.lazydict import LazyLoadDict


def test_sweeps():
    # verify that Radar instances are picklable
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep1 = pyart.core.Sweep(radar, 0)
    sweep2 = pyart.core.Sweep(radar, 1)


def test_sweep_sweep_number():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)

    assert sweep.sweep_number == 0
    assert radar.sweep_number['data'][0] == 0
    sweep.sweep_number = 42
    assert sweep.sweep_number == 42
    assert radar.sweep_number['data'][0] == 42


def test_sweep_sweep_mode():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)

    assert sweep.sweep_mode == 'azimuth_surveillance'
    assert radar.sweep_mode['data'][0] == 'azimuth_surveillance'
    sweep.sweep_mode = 'rhi'
    assert sweep.sweep_mode == 'rhi'
    assert radar.sweep_mode['data'][0] == 'rhi'


def test_sweep_fixed_angle():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)

    assert_almost_equal(sweep.fixed_angle, 0.75, 2)
    assert_almost_equal(radar.fixed_angle['data'][0], 0.75, 2)
    sweep.fixed_angle = 42.12
    assert_almost_equal(sweep.fixed_angle, 42.12, 2)
    assert_almost_equal(radar.fixed_angle['data'][0], 42.12, 2)


def test_sweep_target_scan_rate():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)

    assert radar.target_scan_rate is None
    assert sweep.target_scan_rate is None
    with assert_raises(AttributeError):
        sweep.target_scan_rate = 42

    radar.target_scan_rate = {'data': np.ones(2)}
    assert radar.target_scan_rate['data'][0] == 1
    assert sweep.target_scan_rate == 1
    sweep.target_scan_rate = 42
    assert radar.target_scan_rate['data'][0] == 42
    assert sweep.target_scan_rate == 42


def test_sweep_rays_are_indexed():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)

    assert radar.rays_are_indexed is None
    assert sweep.rays_are_indexed is None
    with assert_raises(AttributeError):
        sweep.rays_are_indexed = 'true'

    radar.rays_are_indexed = {'data': np.array(['true', 'false'])}
    assert radar.rays_are_indexed['data'][0] == 'true'
    assert sweep.rays_are_indexed == 'true'
    sweep.rays_are_indexed = 'false'
    assert radar.rays_are_indexed['data'][0] == 'false'
    assert sweep.rays_are_indexed == 'false'


def test_sweep_ray_angle_res():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)

    assert radar.ray_angle_res is None
    assert sweep.ray_angle_res is None
    with assert_raises(AttributeError):
        sweep.ray_angle_res = 42

    radar.ray_angle_res = {'data': np.ones(2)}
    assert radar.ray_angle_res['data'][0] == 1
    assert sweep.ray_angle_res == 1
    sweep.ray_angle_res = 42
    assert radar.ray_angle_res['data'][0] == 42
    assert sweep.ray_angle_res == 42


def test_sweep_ngates():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)
    assert sweep.ngates == 5
    with assert_raises(AttributeError):
        sweep.ngates = 42


def test_sweep_nrays():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)
    assert sweep.nrays == 4
    with assert_raises(AttributeError):
        sweep.nrays = 42


def test_sweep_picklable():
    # verify that Sweep instances are picklable
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    sweep = pyart.core.Sweep(radar, 0)
    picklestring = pickle.dumps(sweep)
    sweep_new = pickle.loads(picklestring)
    assert 'data' in sweep.gate_x
    assert 'data' in sweep_new.gate_x
