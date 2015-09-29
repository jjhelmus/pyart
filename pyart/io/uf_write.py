import math
import struct
import warnings

import numpy as np
import netCDF4


class UFFileCreator(object):
    """
    A class for creating UF files

    """

    def __init__(self, filename, radar, field_order, volume_start=None):
        if hasattr(filename, 'write'):
            self._fh = filename
        else:
            self._fh = open(filename, 'wb')
        self.filename = filename
        self.radar = radar

        raycreator = UFRayCreator(
            radar, field_order, volume_start=volume_start)

        # XXX
        # hacks to get match with example file
        header = raycreator._mandatory_header_template
        header['radar_name'] = 'xsapr-sg'
        header['site_name'] = 'xsapr-sg'
        header['generation_year'] = 15
        header['generation_month'] = 8
        header['generation_day'] = 19
        header['generation_facility_name'] = 'RSLv1.48'

        header = raycreator._optional_header_template
        header['project_name'] = 'TRMMGVUF'
        header['tape_name'] = 'RADAR_UF'

        for ray_num in range(radar.nrays):

            pad = struct.pack('>i', 16640)
            ray_bytes = raycreator.make_ray(ray_num)

            self._fh.write(pad)
            self._fh.write(ray_bytes)
            self._fh.write(pad)

    def close():
        self._fh.close()


class UFRayCreator(object):
    """
    A class for generating UF file rays for writing to disk

    Parameters
    ----------

    """

    def __init__(self, radar, field_order, volume_start=None):
        """
        """
        self.radar = radar
        self.nfields = len(radar.fields)
        self.field_order = field_order

        self._mandatory_header_template = UF_MANDATORY_HEADER_TEMPLATE.copy()
        self._field_header_template = UF_FIELD_HEADER_TEMPLATE.copy()
        self._optional_header_template = UF_OPTIONAL_HEADER_TEMPLATE.copy()


        self._set_mandatory_header_location()
        self._set_optional_header_volume_time(volume_start)

        return

    def _set_optional_header_volume_time(self, volume_start):
        header = self._optional_header_template
        if volume_start is None:
            volume_start = netCDF4.num2date(
                self.radar.time['data'][0], self.radar.time['units'])
        header['volume_hour'] = volume_start.hour
        header['volume_minute'] = volume_start.minute
        header['volume_second'] = volume_start.second

    def _set_mandatory_header_location(self):
        """
        """
        header = self._mandatory_header_template

        degrees, minutes, seconds = d_to_dms(self.radar.latitude['data'][0])
        header['latitude_degrees'] = int(degrees)
        header['latitude_minutes'] = int(minutes)
        header['latitude_seconds'] = int(seconds * 64)

        degrees, minutes, seconds = d_to_dms(self.radar.longitude['data'][0])
        header['longitude_degrees'] = int(degrees)
        header['longitude_minutes'] = int(minutes)
        header['longitude_seconds'] = int(seconds * 64)

        header['height_above_sea_level'] = int(self.radar.altitude['data'][0])
        return

    def make_ray(self, ray_num):

        ray = self.make_mandatory_header(ray_num)
        ray += self.make_uf_ray_optional_header()
        ray += self.make_data_header()
        field_positions = self.make_field_position_list()
        ray += self.make_field_position()

        for field_info in field_positions:

            data_type = field_info['data_type']
            offset = field_info['offset_field_header'] + 19
            if data_type in [b'VF', b'VE', b'VR', b'VT', b'VP']:
                offset += 2
                vel_header = self.make_fsi_vel()
            else:
                vel_header = b''

            if data_type in [b'PH']:    # XXX hack
                scale = 10
            else:
                scale = 100
            field_header = self.make_field_header(offset, scale)
            data_array = self.make_data_array(data_type, ray_num, scale)

            ray += field_header
            ray += vel_header
            ray += data_array.tostring()

        return ray

    def make_mandatory_header(self, ray_num):

        # time parameters
        dt = netCDF4.num2date(
            self.radar.time['data'][ray_num], self.radar.time['units'])
        header = self._mandatory_header_template
        header['year'] = dt.year - 2000
        header['month'] = dt.month
        header['day'] = dt.day
        header['hour'] = dt.hour
        header['minute'] = dt.minute
        header['second'] = dt.second

        # ray/sweep numbers
        header['record_number'] = header['ray_number'] = ray_num + 1
        header['sweep_number'] = 1  # XXX sweep

        # pointing
        azimuth = self.radar.azimuth['data'][ray_num]
        header['azimuth'] = int(round(azimuth * 64))

        elevation = self.radar.elevation['data'][ray_num]
        header['elevation'] = int(round(elevation * 64))

        fixed_angle = self.radar.fixed_angle['data'][0]  # XXX sweep
        header['fixed_angle'] = int(round(fixed_angle * 64))

        if self.radar.scan_rate is not None:
            scan_rate = self.radar.scan_rate['data'][ray_num]
        else:
            scan_rate = UF_MISSING_VALUE
        header['sweep_rate'] = int(round(scan_rate * 64))

        if self.radar.scan_type in UF_SWEEP_MODES:
            sweep_mode_number = UF_SWEEP_MODES[self.radar.scan_type]
        else:
            warnings.warn(
                'Unknown scan_type: %s, defaulting to PPI' %
                (self.radar.scan_type))
            sweep_mode_number = UF_SWEEP_MODES['ppi']
        header['sweep_mode'] = sweep_mode_number

        header['record_length'] = 8320  # XXX

        return _pack_structure(header, UF_MANDATORY_HEADER)

    def make_uf_ray_optional_header(self):
        header = self._optional_header_template
        return _pack_structure(header, UF_OPTIONAL_HEADER)

    def make_data_header(self):
        data_header = {
            'ray_nrecords': 1,
            'ray_nfields': self.nfields,
            'record_nfields': self.nfields,
        }
        return _pack_structure(data_header, UF_DATA_HEADER)

    def make_field_position(self):

        field_pos = self.make_field_position_list()
        return ''.join(
            [_pack_structure(fp, UF_FIELD_POSITION) for fp in field_pos])

    def make_field_position_list(self):
        field_pos = []
        offset = 62 + self.nfields * 2 + 1  # 87
        for field in self.field_order:
            field_pos.append({
                'data_type': field,
                'offset_field_header': offset,
            })
            offset += self.radar.ngates + 19
            if field in [b'VF', b'VE', b'VR', b'VT', b'VP']:
                offset += 2
        return field_pos

    def make_field_header(self, data_offset, scale_factor=100):
        field_header = self._field_header_template
        field_header['nbins'] = self.radar.ngates
        field_header['data_offset'] = data_offset
        field_header['scale_factor'] = scale_factor
        return _pack_structure(field_header, UF_FIELD_HEADER)

    def make_fsi_vel(self):
        field_header = {
            'nyquist': 1722,
            'spare': 1,
        }
        return _pack_structure(field_header, UF_FSI_VEL)

    def make_data_array(self, field, ray_num, scale=100.):
        field_data = np.round(self.radar.fields[field]['data'][ray_num]*scale)
        return field_data.filled(-32768).astype('>i2')


UF_MISSING_VALUE = -32768

UF_SWEEP_MODES = {
    'calibration': 0,
    'ppi': 1,
    'coplane': 2,
    'rhi': 3,
    'vpt': 4,
    'target': 5,
    'manual': 6,
    'idle': 7,
}

UF_MANDATORY_HEADER_TEMPLATE = {
    'uf_string': 'UF',
    'record_length': 999,
    'offset_optional_header': 46,   # Always include the optional read
    'offset_local_use_header': 60,  # Never include a local use header
    'offset_data_header': 60,       # Data header follows optional header
    'record_number': 999,
    'volume_number': 1,             # always a single volume
    'ray_number': 999,
    'ray_record_number': 1,         # always one recorded per ray
    'sweep_number': 999,
    'radar_name': 'XXXXXXXX',
    'site_name': 'XXXXXXXX',
    'latitude_degrees': 999,
    'latitude_minutes': 999,
    'latitude_seconds': 999,
    'longitude_degrees': 999,
    'longitude_minutes': 999,
    'longitude_seconds': 999,
    'height_above_sea_level': 999,
    'year': 999,
    'month': 999,
    'day': 999,
    'hour': 999,
    'minute': 999,
    'second': 999,
    'time_zone': 'UT',
    'azimuth': 999,
    'elevation': 999,
    'sweep_mode': 999,
    'fixed_angle': 999,
    'sweep_rate': 999,
    'generation_year': 999,
    'generation_month': 999,
    'generation_day': 999,
    'generation_facility_name': 'XXXXXXXX',
    'missing_data_value': UF_MISSING_VALUE,     # standard missing value
}

UF_OPTIONAL_HEADER_TEMPLATE = {
    'project_name': 'XXXXXXXX',
    'baseline_azimuth': UF_MISSING_VALUE,
    'baseline_elevation': UF_MISSING_VALUE,
    'volume_hour': 999,
    'volume_minute': 999,
    'volume_second': 999,
    'tape_name': 'XXXXXXXX',
    'flag': 2,   # default used by RSL
}

UF_FIELD_HEADER_TEMPLATE = {
    'data_offset': 999,
    'scale_factor': 100,
    'range_start_km': 0,
    'range_start_m': 0,
    'range_spacing_m': 60,
    'nbins': 999,
    'pulse_width_m': 134,
    'beam_width_h': 64,
    'beam_width_v': 64,
    'bandwidth': 9670,
    'polarization': 0,
    'wavelength_cm': 198,
    'sample_size': 90,
    'threshold_data': '  ',
    'threshold_value': -32768,
    'scale': -32768,
    'edit_code': '  ',
    'prt_ms': 450,
    'bits_per_bin': 16,
}


def d_to_dms(in_deg):
    # add or subtract a fraction of a second to fix round off issues
    epsilon = 0.01 / 3600.
    in_deg += epsilon * np.sign(in_deg)
    remain, degrees = math.modf(in_deg)
    remain, minutes = math.modf(remain * 60.)
    remain, seconds = math.modf(remain * 60.)
    return degrees, minutes, seconds


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

from pyart.io.uffile import UF_MANDATORY_HEADER
from pyart.io.uffile import UF_OPTIONAL_HEADER



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
