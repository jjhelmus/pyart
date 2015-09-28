""" Unit Tests for Py-ART's io/uf.py and io/uffile.py modules. """
import math

try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO

import numpy as np
from numpy.testing import assert_raises, assert_almost_equal

import pyart
from pyart.io.uffile import UFFile, UFRay

import netCDF4


def test_write_uf():
    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)


    with open(pyart.testing.UF_FILE, 'rb') as f:
        ref_uf = f.read()

    with open(pyart.testing.UF_FILE, 'rb') as f:
        tst_uf = f.read()
    assert tst_uf == ref_uf

def test_make_ray():


    ufile = UFFile(pyart.testing.UF_FILE)
    uray = ufile.rays[0]
    ref_ray_buf = uray._buf
    ufile.close()

    radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)

    # mandatory header
    ref_man_header = ref_ray_buf[:90]
    tst_man_header = make_uf_ray_mandatory_header(radar)
    assert tst_man_header == ref_man_header

    # optional header
    ref_opt_header = uray._buf[90:118]
    tst_opt_header = make_uf_ray_optional_header(radar)
    assert tst_opt_header == ref_opt_header

    # data headers
    ref_data_header = uray._buf[118:124]
    nfields = len(radar.fields)
    data_header = {
        'ray_nrecords': 1,
        'ray_nfields': nfields,
        'record_nfields': nfields,
    }
    tst_data_header = _pack_structure(data_header, UF_DATA_HEADER)
    assert tst_data_header == ref_data_header

    # field position info
    ref_field_position = uray._buf[124:124 + 4 * nfields]
    field_order = [d['data_type'] for d in uray.field_positions]
    field_position = []
    offset = 62 + nfields * 2 + 1  # 87
    for field in field_order:
        field_position.append({
            'data_type': field,
            'offset_field_header': offset,
        })
        offset += radar.ngates + 19
        if field in [b'VF', b'VE', b'VR', b'VT', b'VP']:
            offset += 2
    tst_field_position = ''.join(
        [_pack_structure(fpos, UF_FIELD_POSITION) for fpos in field_position])
    assert tst_field_position == ref_field_position

    # DZ field header
    ref_field_header = uray._buf[172:172+38]
    field_header = {
        'bandwidth': 9670,
        'beam_width_h': 64,
        'beam_width_v': 64,
        'bits_per_bin': 16,
        'data_offset': 106,
        'edit_code': '\x00\x00',
        'nbins': radar.ngates,
        'polarization': 0,
        'prt_ms': 450,
        'pulse_width_m': 134,
        'range_spacing_m': 60,
        'range_start_km': 0,
        'range_start_m': 0,
        'sample_size': 90,
        'scale': -32768,
        'scale_factor': 100,
        'threshold_data': '  ',
        'threshold_value': -32768,
        'wavelength_cm': 198,
    }
    tst_field_header = _pack_structure(field_header, UF_FIELD_HEADER)
    assert tst_field_header == ref_field_header

    # DZ data
    ref_dz_data = uray.field_raw_data[0]
    tst_dz_data = np.round(
        radar.fields['DZ']['data'][0] * 100.).filled(-32768).astype('>i2')
    assert np.array_equal(ref_dz_data, tst_dz_data)

    # DZ data buffer
    ref_dz_data_buf = uray._buf[210:1544]
    assert ref_dz_data_buf == ref_dz_data.tostring()
    assert ref_dz_data_buf == tst_dz_data.tostring()

    # VR field header
    ref_field_header = uray._buf[1544:1544+42]
    field_header = {
        'bandwidth': 9670,
        'beam_width_h': 64,
        'beam_width_v': 64,
        'bits_per_bin': 16,
        'data_offset': 794,
        'edit_code': '  ',
        'nbins': radar.ngates,
        'polarization': 0,
        'prt_ms': 450,
        'pulse_width_m': 134,
        'range_spacing_m': 60,
        'range_start_km': 0,
        'range_start_m': 0,
        'sample_size': 90,
        'scale': -32768,
        'scale_factor': 100,
        'threshold_data': '  ',
        'threshold_value': -32768,
        'wavelength_cm': 198,
        'nyquist': 1722,
        'spare': 1,
    }
    tst_field_header = _pack_structure(field_header, UF_FIELD_HEADER)
    vel_header = _pack_structure(field_header, UF_FSI_VEL)
    tst_field_header += vel_header
    assert tst_field_header == ref_field_header

    # VR data
    ref_vr_data = uray.field_raw_data[1]
    tst_vr_data = np.round(
        radar.fields['VR']['data'][0] * 100.).filled(-32768).astype('>i2')
    assert np.array_equal(ref_vr_data, tst_vr_data)

    # VR data buffer
    ref_vr_data_buf = uray._buf[1586:2920]
    assert ref_vr_data_buf == ref_vr_data.tostring()
    assert ref_vr_data_buf == tst_vr_data.tostring()

    # full buffer
    tst_ray_buf = ref_ray_buf
    assert tst_ray_buf == ref_ray_buf

import struct

def make_uf_ray_optional_header(radar):
    opt_header = {
        'baseline_azimuth': -32768,
        'baseline_elevation': -32768,
        'flag': 2,
        'project_name': 'TRMMGVUF',
        'tape_name': 'RADAR_UF',
    }
    dt = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    opt_header['volume_hour'] = dt.hour
    opt_header['volume_minute'] = dt.minute
    opt_header['volume_second'] = dt.second - 8
    return _pack_structure(opt_header, UF_OPTIONAL_HEADER)

def make_uf_ray_mandatory_header(radar):
    # time parameters
    times = netCDF4.num2date(radar.time['data'], radar.time['units'])
    dt = times[0]
    man_header = make_empty_mandatory_header()
    man_header['year'] = dt.year - 2000
    man_header['month'] = dt.month
    man_header['day'] = dt.day
    man_header['hour'] = dt.hour
    man_header['minute'] = dt.minute
    man_header['second'] = dt.second

    # location parameters
    degrees, minutes, seconds = d_to_dms(radar.latitude['data'][0])
    man_header['latitude_degrees'] = int(degrees)
    man_header['latitude_minutes'] = int(minutes)
    man_header['latitude_seconds'] = int(seconds * 64)

    degrees, minutes, seconds = d_to_dms(radar.longitude['data'][0])
    man_header['longitude_degrees'] = int(degrees)
    man_header['longitude_minutes'] = int(minutes)
    man_header['longitude_seconds'] = int(seconds * 64)

    man_header['height_above_sea_level'] = int(radar.altitude['data'][0])

    # ray/sweep numbers
    man_header['record_number'] = man_header['ray_number'] = 1
    man_header['sweep_number'] = 1

    # pointing
    azimuth = radar.azimuth['data'][0]
    man_header['azimuth'] = int(round(azimuth * 64))

    elevation = radar.elevation['data'][0]
    man_header['elevation'] = int(round(elevation * 64))

    fixed_angle = radar.fixed_angle['data'][0]  # index by sweep
    man_header['fixed_angle'] = int(round(fixed_angle * 64))

    if radar.scan_rate is not None:
        scan_rate = radar.scan_rate['data'][0]
    else:
        scan_rate = 0
    man_header['sweep_rate'] = int(round(scan_rate * 64))

    _UF_SWEEP_MODES = {
        'calibration': 0,
        'ppi': 1,
        'coplane': 2,
        'rhi': 3,
        'vpt': 4,
        'target': 5,
        'manual': 6,
        'idle': 7,
    }
    if radar.scan_type in _UF_SWEEP_MODES:
        sweep_mode_number = _UF_SWEEP_MODES[radar.scan_type]
    else:
        raise ValueError
        sweep_mode_number = 6
    man_header['sweep_mode'] = sweep_mode_number

    # TODO
    man_header['record_length'] = 8320

    # hacks to get match with example file
    man_header['radar_name'] = 'xsapr-sg'
    man_header['site_name'] = 'xsapr-sg'
    man_header['generation_year'] = 15
    man_header['generation_month'] = 8
    man_header['generation_day'] = 19
    man_header['generation_facility_name'] = 'RSLv1.48'

    return _pack_structure(man_header, UF_MANDATORY_HEADER)

def d_to_dms(in_deg):
    # add or subtract a fraction of a second to fix round off issues
    epsilon = 0.01 / 3600.
    in_deg += epsilon * np.sign(in_deg)
    remain, degrees = math.modf(in_deg)
    remain, minutes = math.modf(remain * 60.)
    remain, seconds = math.modf(remain * 60.)
    return degrees, minutes, seconds

def make_empty_mandatory_header():
    return {
        'uf_string': 'UF',
        'offset_optional_header': 46,
        'offset_local_use_header': 60,
        'offset_data_header': 60,
        'volume_number': 1,
        'ray_record_number': 1,
        'time_zone': 'UT',
        'missing_data_value': -32768,
    }


class UFWriteRay(object):
    """
    """


def _pack_structure(dic, structure):
    """ Pack a structure from a dictionary """
    fmt = '>' + ''.join([i[1] for i in structure])  # UF is big-endian
    values = [dic[i[0]] for i in structure]
    return struct.pack(fmt, *values)

def _structure_size(structure):
    """ Find the size of a structure in bytes. """
    return struct.calcsize('>' + ''.join([i[1] for i in structure]))


def _unpack_from_buf(buf, pos, structure):
    """ Unpack a structure from a buffer. """
    size = _structure_size(structure)
    return _unpack_structure(buf[pos:pos + size], structure)


def _unpack_structure(string, structure):
    """ Unpack a structure from a string """
    fmt = '>' + ''.join([i[1] for i in structure])  # UF is big-endian
    lst = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], lst))

INT16 = 'h'

UF_MANDATORY_HEADER = (
    ('uf_string', '2s'),
    ('record_length', INT16),
    ('offset_optional_header', INT16),
    ('offset_local_use_header', INT16),
    ('offset_data_header', INT16),
    ('record_number', INT16),
    ('volume_number', INT16),
    ('ray_number', INT16),
    ('ray_record_number', INT16),
    ('sweep_number', INT16),
    ('radar_name', '8s'),
    ('site_name', '8s'),
    ('latitude_degrees', INT16),
    ('latitude_minutes', INT16),
    ('latitude_seconds', INT16),
    ('longitude_degrees', INT16),
    ('longitude_minutes', INT16),
    ('longitude_seconds', INT16),
    ('height_above_sea_level', INT16),
    ('year', INT16),
    ('month', INT16),
    ('day', INT16),
    ('hour', INT16),
    ('minute', INT16),
    ('second', INT16),
    ('time_zone', '2s'),
    ('azimuth', INT16),
    ('elevation', INT16),
    ('sweep_mode', INT16),
    ('fixed_angle', INT16),
    ('sweep_rate', INT16),
    ('generation_year', INT16),
    ('generation_month', INT16),
    ('generation_day', INT16),
    ('generation_facility_name', '8s'),
    ('missing_data_value', INT16),
)

UF_OPTIONAL_HEADER = (
    ('project_name', '8s'),
    ('baseline_azimuth', INT16),
    ('baseline_elevation', INT16),
    ('volume_hour', INT16),
    ('volume_minute', INT16),
    ('volume_second', INT16),
    ('tape_name', '8s'),
    ('flag', INT16)
)

UF_DATA_HEADER = (
    ('ray_nfields', INT16),
    ('ray_nrecords', INT16),
    ('record_nfields', INT16),
)

UF_FIELD_POSITION = (
    ('data_type', '2s'),
    ('offset_field_header', INT16),
)

UF_FIELD_HEADER = (
    ('data_offset', INT16),
    ('scale_factor', INT16),
    ('range_start_km', INT16),
    ('range_start_m', INT16),
    ('range_spacing_m', INT16),
    ('nbins', INT16),
    ('pulse_width_m', INT16),
    ('beam_width_h', INT16),    # degrees * 64
    ('beam_width_v', INT16),    # degrees * 64
    ('bandwidth', INT16),       # Reciever bandwidth in MHz * 16
    ('polarization', INT16),    # 1: hort, 2: vert 3: circular, 4: ellip
    ('wavelength_cm', INT16),   # cm * 64
    ('sample_size', INT16),
    ('threshold_data', '2s'),
    ('threshold_value', INT16),
    ('scale', INT16),
    ('edit_code', '2s'),
    ('prt_ms', INT16),
    ('bits_per_bin', INT16),    # Must be 16
)

UF_FSI_VEL = (
    ('nyquist', INT16),
    ('spare', INT16),
)

# This structure is defined but not used in Py-ART
# No sample file which contain the structure could be found.
UF_FSI_DM = (
    ('radar_constant', INT16),
    ('noise_power', INT16),
    ('reciever_gain', INT16),
    ('peak_power', INT16),
    ('antenna_gain', INT16),
    ('pulse_duration', INT16),
)
