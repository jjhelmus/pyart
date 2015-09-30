"""
pyart.io.uf_write
=================

Functions for writing UF files.


.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    UFRayCreator

.. autosummary::
    :toctree: generated/

    _structure_size
    _unpack_from_buf
    _unpack_structure

"""

from __future__ import division, unicode_literals

import math
import struct
import warnings

import numpy as np
import netCDF4

from pyart.io.uffile import UF_MANDATORY_HEADER
from pyart.io.uffile import UF_OPTIONAL_HEADER
from pyart.io.uffile import UF_DATA_HEADER
from pyart.io.uffile import UF_FIELD_POSITION
from pyart.io.uffile import UF_FIELD_HEADER
from pyart.io.uffile import UF_FSI_VEL
from pyart.io.uf import _LIGHT_SPEED
from pyart.io.uffile import POLARIZATION_STR


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

        # TODO
        # hacks to get match with example file
        header = raycreator.mandatory_header_template
        header['radar_name'] = b'xsapr-sg'
        header['site_name'] = b'xsapr-sg'
        header['generation_year'] = 15
        header['generation_month'] = 8
        header['generation_day'] = 19
        header['generation_facility_name'] = b'RSLv1.48'

        header = raycreator.optional_header_template
        header['project_name'] = b'TRMMGVUF'
        header['tape_name'] = b'RADAR_UF'

        for ray_num in range(radar.nrays):

            ray_bytes = raycreator.make_ray(ray_num)
            pad = struct.pack('>i', raycreator.record_length * 2)

            self._fh.write(pad)
            self._fh.write(ray_bytes)
            self._fh.write(pad)

    def close(self):
        self._fh.close()


class UFRayCreator(object):
    """
    A class for generating UF file rays for writing to disk

    Parameters
    ----------
    radar : Radar
        Radar used to create rays.
    field_order : list of bytes

    volume_start :

    """

    def __init__(self, radar, field_order, volume_start=None):
        """ Initialize the object. """
        self.radar = radar
        self.field_order = field_order
        self.record_length = self._calc_record_length(radar, field_order)
        self.ray_num_to_sweep_num = self._calc_ray_num_to_sweep_num(radar)

        self.mandatory_header_template = UF_MANDATORY_HEADER_TEMPLATE.copy()
        self.optional_header_template = UF_OPTIONAL_HEADER_TEMPLATE.copy()
        self.field_header_template = UF_FIELD_HEADER_TEMPLATE.copy()

        self._set_mandatory_header_location()
        self._set_optional_header_time(volume_start)
        self._set_field_header()

        return

    @staticmethod
    def _calc_ray_num_to_sweep_num(radar):
        """ Return an array mapping ray number to sweep numbers. """
        ray_num_to_sweep_num = np.zeros((radar.nrays, ), dtype='int32')
        for isweep, sweep_slice in enumerate(radar.iter_slice()):
            ray_num_to_sweep_num[sweep_slice] = isweep
        return ray_num_to_sweep_num

    @staticmethod
    def _calc_record_length(radar, field_order):
        """ Return the record length in 2-byte words. """
        # record length is given by the sum of
        # 45 word (90 byte) mandatory header
        # 14 word (28 byte) optional header which is always written
        # 3 word (6 byte) data header
        # For each field:
        #   * nbin words for the data (stored as int16)
        #   * 2 words for the data type name and offset in the data header
        #   * 19 words for the field header
        # For each velocity-like field:
        #   * 2 words (4-bytes) for the FSI velocity structure
        nfields = len(field_order)
        nvel = sum([field in UF_VEL_DATA_TYPES for field in field_order])
        return 45+14+3 + (radar.ngates+2+19)*nfields + 2*nvel

    def _set_optional_header_time(self, volume_start):
        """ Populate the optional header template with the volume start. """
        header = self.optional_header_template
        if volume_start is None:
            volume_start = netCDF4.num2date(
                self.radar.time['data'][0], self.radar.time['units'])
        header['volume_hour'] = volume_start.hour
        header['volume_minute'] = volume_start.minute
        header['volume_second'] = volume_start.second

    def _set_mandatory_header_location(self):
        """ Populate the mandatory header template with the location. """
        header = self.mandatory_header_template

        degrees, minutes, seconds = _d_to_dms(self.radar.latitude['data'][0])
        header['latitude_degrees'] = int(degrees)
        header['latitude_minutes'] = int(minutes)
        header['latitude_seconds'] = int(seconds * 64)

        degrees, minutes, seconds = _d_to_dms(self.radar.longitude['data'][0])
        header['longitude_degrees'] = int(degrees)
        header['longitude_minutes'] = int(minutes)
        header['longitude_seconds'] = int(seconds * 64)

        header['height_above_sea_level'] = int(self.radar.altitude['data'][0])
        return

    def _set_field_header(self):
        """ Populate the field header template with radar parameters. """
        # these parameter are assumed to be constant for all rays in
        # the volume, UF can store volumes where these change between
        # rays but Py-ART does not support writing such volumes
        header = self.field_header_template

        range_step = self.radar.range['meters_between_gates']
        range_start = self.radar.range['meters_to_center_of_first_gate']
        range_start -= (range_step / 2.)  # range bin center to edge
        header['range_start_km'] = int(round(range_start)) * 1000
        header['range_start_m'] = int(round(range_start))
        header['range_spacing_m'] = int(round(range_step))

        iparams = self.radar.instrument_parameters

        if iparams is not None and 'radar_beam_width_h' in iparams:
            beam_width_h = iparams['radar_beam_width_h']['data'][0]
            header['beam_width_h'] = int(round(beam_width_h * 64))
        else:
            header['beam_width_h'] = UF_MISSING_VALUE

        if iparams is not None and 'radar_beam_width_v' in iparams:
            beam_width_v = iparams['radar_beam_width_v']['data'][0]
            header['beam_width_v'] = int(round(beam_width_v * 64))
        else:
            header['beam_width_v'] = UF_MISSING_VALUE

        if iparams is not None and 'radar_receiver_bandwidth' in iparams:
            bandwidth = iparams['radar_receiver_bandwidth']['data'][0]
            header['bandwidth'] = int(round(bandwidth / 1.e6 * 16))
        else:
            header['bandwidth'] = UF_MISSING_VALUE

        if iparams is not None and 'frequency' in iparams:
            frequency = iparams['frequency']['data'][0]
            header['wavelength_cm'] = int(round(
                _LIGHT_SPEED / frequency * 100 * 64))
        else:
            header['wavelength_cm'] = UF_MISSING_VALUE

        return

    def make_ray(self, ray_num):
        """ Return a byte string representing a complete UF ray. """
        ray = self.make_mandatory_header(ray_num)
        ray += self.make_optional_header()
        ray += self.make_data_header()
        field_positions = self.make_field_position_list()
        ray += self.make_field_position()

        for field_info in field_positions:

            data_type = field_info['data_type']
            offset = field_info['offset_field_header'] + 19
            if '_UF_scale_factor' in self.radar.fields[data_type]:
                scale = self.radar.fields[data_type]['_UF_scale_factor']
            else:
                scale = UF_DEFAULT_SCALE_FACTOR

            if data_type in UF_VEL_DATA_TYPES:
                offset += 2
                vel_header = self.make_fsi_vel(ray_num, scale)
            else:
                vel_header = b''

            field_header = self.make_field_header(offset, ray_num, scale)
            data_array = self.make_data_array(data_type, ray_num, scale)

            ray += field_header
            ray += vel_header
            ray += data_array.tostring()

        return ray

    def make_mandatory_header(self, ray_num):
        """ Return a byte string representing a UF mandatory header. """

        # time parameters
        ray_time = netCDF4.num2date(
            self.radar.time['data'][ray_num], self.radar.time['units'])
        header = self.mandatory_header_template
        header['year'] = ray_time.year - 2000
        header['month'] = ray_time.month
        header['day'] = ray_time.day
        header['hour'] = ray_time.hour
        header['minute'] = ray_time.minute
        header['second'] = ray_time.second

        # ray/sweep numbers
        sweep_num = self.ray_num_to_sweep_num[ray_num]
        header['record_number'] = header['ray_number'] = ray_num + 1
        header['sweep_number'] = sweep_num + 1

        # pointing
        azimuth = self.radar.azimuth['data'][ray_num]
        header['azimuth'] = int(round(azimuth * 64))

        elevation = self.radar.elevation['data'][ray_num]
        header['elevation'] = int(round(elevation * 64))

        fixed_angle = self.radar.fixed_angle['data'][sweep_num]
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

        header['record_length'] = self.record_length

        return _pack_structure(header, UF_MANDATORY_HEADER)

    def make_optional_header(self):
        """ Return a byte string representing a UF optional header. """
        header = self.optional_header_template
        return _pack_structure(header, UF_OPTIONAL_HEADER)

    def make_data_header(self):
        """ Return a byte string representing a UF data header. """
        header = UF_DATA_HEADER_TEMPLATE.copy()
        nfields = len(self.field_order)
        header['ray_nfields'] = nfields
        header['record_nfields'] = nfields
        return _pack_structure(header, UF_DATA_HEADER)

    def make_field_position_list(self):
        """ Return a list of field position dictionaries. """
        # 62 words (124 bytes) are occuplied by the:
        # * mandatory header (90 bytes)
        # * optional header (28)
        # * data header (6)
        # This is followed by 2 words (4 bytes) for each field name/offset
        # Finally one must be added since the offset has origin of 1
        offset = 62 + len(self.field_order) * 2 + 1
        field_positions = []
        for data_type in self.field_order:
            field_position = UF_FIELD_POSITION_TEMPLATE.copy()
            field_position['data_type'] = bytes(data_type)
            field_position['offset_field_header'] = offset
            field_positions.append(field_position)

            # account for field header and data
            offset += self.radar.ngates + 19
            if data_type in UF_VEL_DATA_TYPES:
                offset += 2  # account for the FSI_VEL structure
        return field_positions

    def make_field_position(self):
        """ Return a byte string representing the UF field positions. """
        fps = self.make_field_position_list()
        return b''.join([_pack_structure(fp, UF_FIELD_POSITION) for fp in fps])

    def make_field_header(self, data_offset, ray_num, scale_factor):
        """ Return a byte string representing a field header. """
        field_header = self.field_header_template
        field_header['nbins'] = self.radar.ngates
        field_header['data_offset'] = data_offset
        field_header['scale_factor'] = scale_factor

        iparams = self.radar.instrument_parameters
        if iparams is not None and 'pulse_width' in iparams:
            pulse_width = iparams['pulse_width']['data'][ray_num]
            field_header['pulse_width_m'] = int(round(
                pulse_width * _LIGHT_SPEED))
        else:
            field_header['pulse_width_m'] = UF_MISSING_VALUE

        if iparams is not None and 'prt' in iparams:
            prt = iparams['prt']['data'][ray_num]
            field_header['prt_ms'] = int(round(prt * 1.e6))
        else:
            field_header['prt_ms'] = UF_MISSING_VALUE

        field_header['polarization'] = 1  # default to horizontal polarization
        sweep_num = self.ray_num_to_sweep_num[ray_num]
        if iparams is not None and 'polarization_mode' in iparams:
            mode = str(iparams['polarization_mode']['data'][sweep_num])
            if mode in POLARIZATION_STR:
                field_header['polarization'] = POLARIZATION_STR.index(mode)
        return _pack_structure(field_header, UF_FIELD_HEADER)

    def make_fsi_vel(self, ray_num, scale):
        """ Return a byte string representing a UF FSI velocity structure. """
        fsi_vel = UF_FSI_VEL_TEMPLATE.copy()
        iparams = self.radar.instrument_parameters
        if iparams is not None and 'nyquist_velocity' in iparams:
            nyquist = iparams['nyquist_velocity']['data'][ray_num]
            fsi_vel['nyquist'] = int(round(nyquist * scale))
        else:
            fsi_vel['nyquist'] = UF_MISSING_VALUE
        return _pack_structure(fsi_vel, UF_FSI_VEL)

    def make_data_array(self, field, ray_num, scale=100.):
        """ Return an array of UF field data. """
        field_data = np.round(self.radar.fields[field]['data'][ray_num]*scale)
        return field_data.filled(-32768).astype('>i2')


def _d_to_dms(in_deg):
    """ Degrees to degree, minutes, seconds. """
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


# Constants
UF_MISSING_VALUE = -32768
UF_DEFAULT_SCALE_FACTOR = 100   # default field scale factor
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
UF_VEL_DATA_TYPES = [b'VF', b'VE', b'VR', b'VT', b'VP']

# Structure dictionary templates
UF_MANDATORY_HEADER_TEMPLATE = {
    'uf_string': b'UF',
    'record_length': 999,
    'offset_optional_header': 46,   # Always include the optional read
    'offset_local_use_header': 60,  # Never include a local use header
    'offset_data_header': 60,       # Data header follows optional header
    'record_number': 999,
    'volume_number': 1,             # always a single volume
    'ray_number': 999,
    'ray_record_number': 1,         # always one recorded per ray
    'sweep_number': 999,
    'radar_name': b'XXXXXXXX',
    'site_name': b'XXXXXXXX',
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
    'time_zone': b'UT',
    'azimuth': 999,
    'elevation': 999,
    'sweep_mode': 999,
    'fixed_angle': 999,
    'sweep_rate': 999,
    'generation_year': 999,
    'generation_month': 999,
    'generation_day': 999,
    'generation_facility_name': b'XXXXXXXX',
    'missing_data_value': UF_MISSING_VALUE,     # standard missing value
}

UF_OPTIONAL_HEADER_TEMPLATE = {
    'project_name': b'XXXXXXXX',
    'baseline_azimuth': UF_MISSING_VALUE,
    'baseline_elevation': UF_MISSING_VALUE,
    'volume_hour': 999,
    'volume_minute': 999,
    'volume_second': 999,
    'tape_name': 'XXXXXXXX',
    'flag': 2,   # default used by RSL
}

UF_DATA_HEADER_TEMPLATE = {
    'ray_nfields': 999,
    'ray_nrecords': 1,  # 1 record per ray
    'record_nfields': 999,
}

UF_FIELD_POSITION_TEMPLATE = {
    'data_type': b'XX',
    'offset_field_header': 999,
}

UF_FIELD_HEADER_TEMPLATE = {
    'data_offset': 999,
    'scale_factor': 999,
    'range_start_km': 999,
    'range_start_m': 999,
    'range_spacing_m': 999,
    'nbins': 999,
    'pulse_width_m': 999,
    'beam_width_h': 999,
    'beam_width_v': 999,
    'bandwidth': 999,
    'polarization': 999,
    'wavelength_cm': 999,
    'sample_size': 90,          # Apparently defaults to 90?
    'threshold_data': b'  ',    # No thresholding field
    'threshold_value': UF_MISSING_VALUE,
    'scale': UF_MISSING_VALUE,
    'edit_code': b'  ',         # unknown use, typically blank or null
    'prt_ms': 999,
    'bits_per_bin': 16,         # 16 bits (2 bytes, 1 word) per bin or gate
}

UF_FSI_VEL_TEMPLATE = {
    'nyquist': 999,
    'spare': 1,     # 1 if commonly used for this unused element
}