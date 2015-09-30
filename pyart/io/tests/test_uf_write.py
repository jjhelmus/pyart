""" Unit Tests for Py-ART's uf_write module. """

from __future__ import unicode_literals

try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO

import datetime

import numpy as np
from numpy.testing import assert_raises, assert_almost_equal
import netCDF4

import pyart
from pyart.io.uffile import UFFile
from pyart.io.uf_write import UFRayCreator, write_uf


TEMPLATES_EXTRA = {
    'mandatory_header': {
        'radar_name': b'xsapr-sg',
        'site_name': b'xsapr-sg',
        'generation_year': 15,
        'generation_month': 8,
        'generation_day': 19,
        'generation_facility_name': b'RSLv1.48',
    },
    'optional_header': {
        'project_name': b'TRMMGVUF',
        'tape_name': b'RADAR_UF',
    }
}

FIELD_MAPPING = {
    'DZ': 'DZ',
    'VR': 'VR',
    'SW': 'SW',
    'CZ': 'CZ',
    'ZT': 'ZT',
    'DR': 'DR',
    'ZD': 'ZD',
    'RH': 'RH',
    'PH': 'PH',
    'KD': 'KD',
    'SQ': 'SQ',
    'HC': 'HC',
}


def test_ray_section_by_section():

    ufile = UFFile(pyart.testing.UF_FILE)
    uray = ufile.rays[0]
    ref_ray_buf = uray._buf
    ufile.close()

    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)
    volume_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    volume_start -= datetime.timedelta(seconds=8)
    nfields = len(radar.fields)
    field_write_order = [d['data_type'] for d in uray.field_positions]
    ufraycreator = UFRayCreator(
        radar, FIELD_MAPPING, field_write_order, volume_start=volume_start,
        templates_extra=TEMPLATES_EXTRA)

    # mandatory header
    ref_man_header = ref_ray_buf[:90]
    tst_man_header = ufraycreator.make_mandatory_header(0)
    assert tst_man_header == ref_man_header

    # optional header
    ref_opt_header = uray._buf[90:118]
    tst_opt_header = ufraycreator.make_optional_header()
    assert tst_opt_header == ref_opt_header

    # data headers
    ref_data_header = uray._buf[118:124]
    tst_data_header = ufraycreator.make_data_header()
    assert tst_data_header == ref_data_header

    # field position info
    ref_field_position = uray._buf[124:124 + 4 * nfields]
    tst_field_position = ufraycreator.make_field_position()
    assert tst_field_position == ref_field_position

    # DZ field header
    ref_field_header = uray._buf[172:172+38]
    ufraycreator.field_header_template['edit_code'] = b'\x00\x00'
    tst_field_header = ufraycreator.make_field_header(106, 0, 100)
    assert tst_field_header == ref_field_header

    # DZ data
    ref_dz_data = uray.field_raw_data[0]
    tst_dz_data = ufraycreator.make_data_array('DZ', 0)
    assert np.array_equal(ref_dz_data, tst_dz_data)

    # DZ data buffer
    ref_dz_data_buf = uray._buf[210:1544]
    assert ref_dz_data_buf == ref_dz_data.tostring()
    assert ref_dz_data_buf == tst_dz_data.tostring()

    # VR field header
    ref_field_header = uray._buf[1544:1544+42]
    ufraycreator.field_header_template['edit_code'] = b'  '
    tst_field_header = ufraycreator.make_field_header(794, 0, 100)
    vel_header = ufraycreator.make_fsi_vel(0, 100)
    assert tst_field_header + vel_header == ref_field_header

    # VR data
    ref_vr_data = uray.field_raw_data[1]
    tst_vr_data = ufraycreator.make_data_array('VR', 0)
    assert np.array_equal(ref_vr_data, tst_vr_data)

    # VR data buffer
    ref_vr_data_buf = uray._buf[1586:2920]
    assert ref_vr_data_buf == ref_vr_data.tostring()
    assert ref_vr_data_buf == tst_vr_data.tostring()

    # SW field header
    ref_field_header = uray._buf[2920:2920+38]
    ufraycreator.field_header_template['edit_code'] = b'  '
    tst_field_header = ufraycreator.make_field_header(1480, 0, 100)
    assert tst_field_header == ref_field_header

    # SW data
    ref_sw_data = uray.field_raw_data[2]
    tst_sw_data = ufraycreator.make_data_array('SW', 0)
    assert np.array_equal(ref_sw_data, tst_sw_data)

    # SW data buffer
    ref_sw_data_buf = uray._buf[2958:4292]
    assert ref_sw_data_buf == ref_sw_data.tostring()
    assert ref_sw_data_buf == tst_sw_data.tostring()

    # ZT field header
    ref_field_header = uray._buf[5664:5664+38]
    ufraycreator.field_header_template['edit_code'] = b'\x00\x00'
    tst_field_header = ufraycreator.make_field_header(2852, 0, 100)
    assert tst_field_header == ref_field_header

    # PH field header
    ref_field_header = uray._buf[11152:11152+38]
    ufraycreator.field_header_template['edit_code'] = b'  '
    tst_field_header = ufraycreator.make_field_header(5596, 0, 10)
    assert tst_field_header == ref_field_header

    # PH data
    ref_ph_data = uray.field_raw_data[8]
    tst_ph_data = ufraycreator.make_data_array('PH', 0, 10.)
    assert np.array_equal(ref_ph_data, tst_ph_data)

    # PH data buffer
    ref_ph_data_buf = uray._buf[11152+38:12524]
    assert ref_ph_data_buf == ref_ph_data.tostring()
    assert ref_ph_data_buf == tst_ph_data.tostring()


def test_ray_full():

    ufile = UFFile(pyart.testing.UF_FILE)
    ref_ray = ufile.rays[0]._buf[:]
    ufile.close()

    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)
    radar.fields['PH']['_UF_scale_factor'] = 10
    field_write_order = ['DZ', 'VR', 'SW', 'CZ', 'ZT', 'DR', 'ZD', 'RH', 'PH',
                         'KD', 'SQ', 'HC']
    volume_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    volume_start -= datetime.timedelta(seconds=8)
    ufraycreator = UFRayCreator(
        radar, FIELD_MAPPING, field_write_order, volume_start=volume_start,
        templates_extra=TEMPLATES_EXTRA)
    tst_ray = ufraycreator.make_ray(0)
    tst_ray = tst_ray[:204] + b'\x00\x00' + tst_ray[206:]   # DZ edit_code
    tst_ray = tst_ray[:5696] + b'\x00\x00' + tst_ray[5698:]     # ZT edit_code
    assert ref_ray == tst_ray


def test_complete_file():

    with open(pyart.testing.UF_FILE, 'rb') as fh:
        ref_file = fh.read()

    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)
    radar.fields['PH']['_UF_scale_factor'] = 10
    field_write_order = ['DZ', 'VR', 'SW', 'CZ', 'ZT', 'DR', 'ZD', 'RH', 'PH',
                         'KD', 'SQ', 'HC']

    in_mem = StringIO()
    volume_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    volume_start -= datetime.timedelta(seconds=8)
    write_uf(in_mem, radar, uf_field_names=FIELD_MAPPING,
             field_write_order=field_write_order, volume_start=volume_start,
             templates_extra=TEMPLATES_EXTRA)
    in_mem.seek(0)
    tst_file = in_mem.read()
    tst_file = tst_file[:208] + b'\x00\x00' + tst_file[210:]    # DZ edit_code
    tst_file = tst_file[:5700] + b'\x00\x00' + tst_file[5702:]  # ZT edit_code

    assert ref_file == tst_file


def test_complete_file_standard_names():

    with open(pyart.testing.UF_FILE, 'rb') as fh:
        ref_file = fh.read()

    radar = pyart.io.read_uf(pyart.testing.UF_FILE)
    radar.fields['differential_phase']['_UF_scale_factor'] = 10

    field_write_order = [
        'reflectivity',
        'velocity',
        'spectrum_width',
        'corrected_reflectivity',
        'total_power',
        'corrected_differential_reflectivity',
        'differential_reflectivity',
        'cross_correlation_ratio',
        'differential_phase',
        'specific_differential_phase',
        'normalized_coherent_power',
        'radar_echo_classification',
    ]
    in_mem = StringIO()
    volume_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    volume_start -= datetime.timedelta(seconds=8)
    write_uf(in_mem, radar, field_write_order=field_write_order,
             volume_start=volume_start, templates_extra=TEMPLATES_EXTRA)
    in_mem.seek(0)
    tst_file = in_mem.read()
    tst_file = tst_file[:208] + b'\x00\x00' + tst_file[210:]    # DZ edit_code
    tst_file = tst_file[:5700] + b'\x00\x00' + tst_file[5702:]  # ZT edit_code

    assert ref_file == tst_file


def test_write_radar_field_names():
    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)
    in_mem = StringIO()
    write_uf(in_mem, radar, radar_field_names=True)
    assert in_mem.tell() == 16648


def test_write_defaults():
    radar = pyart.io.read_uf(pyart.testing.UF_FILE)
    in_mem = StringIO()
    write_uf(in_mem, radar)
    assert in_mem.tell() == 16648
