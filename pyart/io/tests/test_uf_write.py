""" Unit Tests for Py-ART's io/uf.py and io/uffile.py modules. """
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO

import numpy as np
from numpy.testing import assert_raises, assert_almost_equal

import pyart
from pyart.io.uffile import UFFile, UFRay
from pyart.io.uf_write import UFRayCreator, UFFileCreator

import struct

def test_ray_section_by_section():

    ufile = UFFile(pyart.testing.UF_FILE)
    uray = ufile.rays[0]
    ref_ray_buf = uray._buf
    ufile.close()

    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)
    nfields = len(radar.fields)
    field_order = [d['data_type'] for d in uray.field_positions]
    ufraycreator = UFRayCreator(radar, field_order)

    # mandatory header
    ref_man_header = ref_ray_buf[:90]
    tst_man_header = ufraycreator.make_uf_ray_mandatory_header()
    assert tst_man_header == ref_man_header

    # optional header
    ref_opt_header = uray._buf[90:118]
    tst_opt_header = ufraycreator.make_uf_ray_optional_header()
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
    ufraycreator._field_header_template['edit_code'] = '\x00\x00'
    tst_field_header = ufraycreator.make_field_header(106)
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
    ufraycreator._field_header_template['edit_code'] = '  '
    tst_field_header = ufraycreator.make_field_header(794)
    vel_header = ufraycreator.make_fsi_vel()
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
    ufraycreator._field_header_template['edit_code'] = '  '
    tst_field_header = ufraycreator.make_field_header(1480)
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
    ufraycreator._field_header_template['edit_code'] = '\x00\x00'
    tst_field_header = ufraycreator.make_field_header(2852)
    assert tst_field_header == ref_field_header

    # PH field header
    ref_field_header = uray._buf[11152:11152+38]
    ufraycreator._field_header_template['edit_code'] = '  '
    tst_field_header = ufraycreator.make_field_header(5596, 10)
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
    field_order = ['DZ', 'VR', 'SW', 'CZ', 'ZT', 'DR', 'ZD', 'RH', 'PH',
                   'KD', 'SQ', 'HC']
    ufraycreator = UFRayCreator(radar, field_order)

    tst_ray = ufraycreator.make_ray(0)
    tst_ray = tst_ray[:204] + b'\x00\x00' + tst_ray[206:]   # DZ edit_code
    tst_ray = tst_ray[:5696] + b'\x00\x00' + tst_ray[5698:]     # ZT edit_code
    assert ref_ray == tst_ray


def test_complete_file():

    with open(pyart.testing.UF_FILE, 'rb') as fh:
        ref_file = fh.read()

    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)
    field_order = ['DZ', 'VR', 'SW', 'CZ', 'ZT', 'DR', 'ZD', 'RH', 'PH',
                   'KD', 'SQ', 'HC']

    in_mem = StringIO()
    UFFileCreator(in_mem, radar, field_order)
    in_mem.seek(0)
    tst_file = in_mem.read()
    tst_file = tst_file[:208] + b'\x00\x00' + tst_file[210:]    # DZ edit_code
    tst_file = tst_file[:5700] + b'\x00\x00' + tst_file[5702:]  # ZT edit_code

    assert ref_file == tst_file
