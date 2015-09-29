import math

import struct

import numpy as np
import netCDF4


class UFRayCreator(object):
    """
    A class for generating UF file rays for writing to disk

    Parameters
    ----------

    """

    def __init__(self, radar, field_order):
        """
        """
        self.radar = radar
        self.nfields = len(radar.fields)
        self.field_order = field_order
        self._field_header_template = UF_FIELD_HEADER_TEMPLATE.copy()
        return

    def make_data_header(self):
        data_header = {
            'ray_nrecords': 1,
            'ray_nfields': self.nfields,
            'record_nfields': self.nfields,
        }
        return _pack_structure(data_header, UF_DATA_HEADER)

    def make_uf_ray_mandatory_header(self):
        # time parameters
        times = netCDF4.num2date(
            self.radar.time['data'], self.radar.time['units'])
        dt = times[0]
        man_header = make_empty_mandatory_header()
        man_header['year'] = dt.year - 2000
        man_header['month'] = dt.month
        man_header['day'] = dt.day
        man_header['hour'] = dt.hour
        man_header['minute'] = dt.minute
        man_header['second'] = dt.second

        # location parameters
        degrees, minutes, seconds = d_to_dms(self.radar.latitude['data'][0])
        man_header['latitude_degrees'] = int(degrees)
        man_header['latitude_minutes'] = int(minutes)
        man_header['latitude_seconds'] = int(seconds * 64)

        degrees, minutes, seconds = d_to_dms(self.radar.longitude['data'][0])
        man_header['longitude_degrees'] = int(degrees)
        man_header['longitude_minutes'] = int(minutes)
        man_header['longitude_seconds'] = int(seconds * 64)

        man_header['height_above_sea_level'] = int(
            self.radar.altitude['data'][0])

        # ray/sweep numbers
        man_header['record_number'] = man_header['ray_number'] = 1
        man_header['sweep_number'] = 1

        # pointing
        azimuth = self.radar.azimuth['data'][0]
        man_header['azimuth'] = int(round(azimuth * 64))

        elevation = self.radar.elevation['data'][0]
        man_header['elevation'] = int(round(elevation * 64))

        fixed_angle = self.radar.fixed_angle['data'][0]  # index by sweep
        man_header['fixed_angle'] = int(round(fixed_angle * 64))

        if self.radar.scan_rate is not None:
            scan_rate = self.radar.scan_rate['data'][0]
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
        if self.radar.scan_type in _UF_SWEEP_MODES:
            sweep_mode_number = _UF_SWEEP_MODES[self.radar.scan_type]
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

    def make_uf_ray_optional_header(self):
        opt_header = {
            'baseline_azimuth': -32768,
            'baseline_elevation': -32768,
            'flag': 2,
            'project_name': 'TRMMGVUF',
            'tape_name': 'RADAR_UF',
        }
        dt = netCDF4.num2date(
            self.radar.time['data'][0], self.radar.time['units'])
        opt_header['volume_hour'] = dt.hour
        opt_header['volume_minute'] = dt.minute
        opt_header['volume_second'] = dt.second - 8
        return _pack_structure(opt_header, UF_OPTIONAL_HEADER)

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

    def make_ray(self, ray_num):

        ray = self.make_uf_ray_mandatory_header()
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
