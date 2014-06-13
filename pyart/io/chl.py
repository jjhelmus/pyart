"""
pyart.io.chl
============

Utilities for reading CSU-CHILL CHL files.

"""

import struct
import numpy as np

from ..config import FileMetadata, get_fillvalue
from .radar import Radar
from .common import make_time_unit_str
from .common import radar_coords_to_cart


def read_chl(filename):
    """
    Read a CHL file.

    Parameters
    ----------
    filename : str
        Name of CHL file.

    Returns
    -------
    radar : Radar
        Radar object containing data from CHL file.

    Notes
    -----
    This is still an alpha-level function so use with caution.

    """

    # create metadata retrival object
    filemetadata = FileMetadata('chill')    # XXX additional parameters

    # read data
    chl_file = CHLfile(filename)

    # time
    time = filemetadata('time')
    time['data'] = chl_file.time
    time['units'] = 'seconds since 1970-01-01 00:00 UTC'    # XXX

    # range
    _range = filemetadata('range')
    _range['data'] = chl_file._range
    _range['meters_between_gates'] = np.array(chl_file.dr)
    _range['meters_to_center_of_first_gate'] = 0.0

    chl_file._sweep_number = {
        'data': range(chl_file.sweep_num), 'long_name': 'Sweep_number',
        'standard_name': 'sweep_number', 'units': 'count'}

    # scan_type
    scan_type = {'data': chl_file.scan_mode}    # XXX

    # fields XXX
    fields = {}
    for i in range(0, len(chl_file.field_scale_list)):

        chl_field_name = chl_file.field_scale_list[i]['name']
        if(chl_field_name not in chl_file.fields.keys()):
            continue

        field_name = CHILL_FIELD_MAPPING[chl_field_name]
        field_scale = chl_file.field_scale_list[i]
        field_dic = {
            'coordinates': 'elevation azimuth range',
            'data': np.ma.masked_array(chl_file.fields[chl_field_name]),
            'long_name': CHILL_FIELD_LONG_NAME[chl_field_name],
            'standard_name': field_name,
            'units': field_scale['units'],
            'valid_max': field_scale['max_val'],
            'valid_min': field_scale['min_val']
        }
        fields[field_name] = field_dic

    # metadata
    metadata = filemetadata('metadata')
    metadata.update(chl_file.metadata)

    # longitude, latitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = chl_file._radar_info['latitude']
    longitude['data'] = chl_file._radar_info['longitude']
    altitude['data'] = chl_file._radar_info['altitude']

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_number['data'] = chl_file._sweep_number
    sweep_mode['data'] = chl_file.scan_mode
    fixed_angle['data'] = chl_file.fixed_angle
    sweep_start_ray_index['data'] = chl_file.sweep_start
    sweep_end_ray_index['data'] = chl_file.sweep_end

    # azimuth, elevation
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    azimuth['data'] = chl_file.azimuth
    elevation['data'] = chl_file.elevation

    # instrument parameters
    instrument_parameters = None

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


class CHLfile(object):
    """
    A file object for CHL data.

    A `CHLFile` object stores metadata and data from a CHL file.  Metadata is
    stored in dictionaries as attributes of the object, field data is
    stored as NumPy ndarrays as attributes with the field name.

    Parameters
    ----------
    filename : str
        Name of CHL file to read.
    """

    def __init__(self, filename):
        self.filename = filename

        self.field_scale_list = {}
        self.num_sweeps = 0
        self.radar_info = []
        self.processor_info = []
        self.azimuth = []
        self.elevation = []
        self.fixed_angle = []
        self.sweep_end = []
        self.data_array = []
        self.time = []
        self.fields = {}
        self._current_ray_num = 0
        self.metadata = {}

        self._chl_arch_open_archive()
        self._chl_close_archive()

        self.sweep_end.append(self._current_ray_num)
        del(self.sweep_end[0])
        self._process_data_blocks()
        self._range = np.array(range(self.num_gates)) * self.dr

        self.sweep_start = [0, ]
        self.sweep_start.extend([idx + 1 for idx in self.sweep_end[0:-1]])

    def _chl_arch_open_archive(self):
        self.f = open(self.filename, "rb")
        packet = 1
        while packet is not None:
            packet = self._chl_arch_read_block()

    def _chl_close_archive(self):
        self.f.close()

    # I need to extract this to several shorter funcs
    def _chl_arch_read_block(self):

        pld = self.f.read(8)
        if pld == '':
            return None
        id, length = struct.unpack("<2i", pld)
        payload = self.f.read(length - 8)

        if id is None:  # redundant sanity check
            return None

        packet = {}
        if hex(id) == '0x5aa80004':  # arch_file_hdr_t
            packet = _unpack_structure(payload, ARCH_FILE_HDR_T)

        elif hex(id) == '0x5aa80002':  # field_scale_t
            packet = _unpack_structure(payload, FIELD_SCALE_T)
            self._parse_field_struct_packet(packet)

        elif hex(id) == '0x5aa80003':  # arch_ray_header
            packet = _unpack_structure(payload, ARCH_RAY_HEADER)

            self._current_ray_num = packet['ray_number']
            fmat_string = self._format_string_from_bitmask(packet['bit_mask'])
            data_packet = self.f.read(
                struct.calcsize(fmat_string) * packet['gates'])
            self.data_array.append(
                struct.unpack(fmat_string * packet['gates'], data_packet))
            # We should only do this once, but this will work for now.
            self.field_count = self._num_fields_from_bitmask(
                packet['bit_mask'])
            self.bit_mask = packet['bit_mask']
            self.time.append(packet['time'])
            self.azimuth.append(packet['azimuth'])
            self.elevation.append(packet['elevation'])
            # This assumes number of gates is constant. Not the best assumption
            # however.
            self.num_gates = packet['gates']

        elif hex(id) == '0x5aa50001':  # radar_info_t
            packet = _unpack_structure(payload, RADAR_INFO_T)
            self._radar_info = packet.copy()
            self.metadata['instrument_name'] = packet[
                'radar_name'].rstrip('\x00')
            self.metadata['original_container'] = 'CHL'

        elif hex(id) == '0x5aa50003':  # processor_info
            packet = _unpack_structure(payload, PROCESSOR_INFO)
            self.dr = packet['gate_spacing']

        elif hex(id) == '0x5aa50002':  # scan_seg
            packet = _unpack_structure(payload, SCAN_SEG)
            self.sweep_num = packet['sweep_num']
            self.sweep_end.append(self._current_ray_num)
            self.fixed_angle.append(packet['current_fixed_angle'])
            self.scan_mode = self._scan_mode_names[packet['scan_type']]

        elif hex(id) == '0x5aa80005':
            packet['num_sweeps'] = struct.unpack('I', payload[0:4])[0]
            packet['swp_offsets'] = struct.unpack(
                str(packet['num_sweeps']) + 'Q', payload[4:])
            self.num_sweeps = packet['num_sweeps']

        packet['id'] = hex(id)
        packet['length'] = length
        return packet

    def _process_data_blocks(self):
        darray = np.array(self.data_array)
        cvi = 0
        for cv in range(0, len(self.field_scale_list)):
            if(bool(self.bit_mask & 2 ** cv)):
                fsl = self.field_scale_list[cv]

                if(self.field_scale_list[cv]['data_size'] == 'f'):
                    dat = darray[:, cvi::self.field_count].copy()
                    dat[dat == 0] = np.nan
                else:
                    dat = darray[:, cvi::self.field_count].copy()
                    dat[dat == 0] = np.nan
                    # I need to come back and do whole np.ma thing
                    dat = (dat * fsl['dat_factor'] + fsl['dat_bias']) / \
                        float(fsl['fld_factor'])
                self.fields[fsl['name']] = dat
                cvi += 1

    def _format_string_from_bitmask(self, bitmask):
        format_str = ''
        for b in range(0, 38):
            if(bool(bitmask & 2 ** b)):
                format_str += self.field_scale_list[b]['data_size']
        return format_str

    def _num_fields_from_bitmask(self, bitmask):
        count = 0
        for b in range(0, 38):
            if(bool(bitmask & 2 ** b)):
                count += 1
        return count

    def _parse_field_struct_packet(self, packet):
        packet['name'] = packet['name'].rstrip('\x00')
        packet['units'] = packet['units'].rstrip('\x00')
        packet['descr'] = packet['descr'].rstrip('\x00')
        packet['data_size'] = self.format_length_lookup_list[packet['format']]
        self.field_scale_list[packet['bit_mask_pos']] = packet

    format_length_lookup_list = [
        'c',
        'Q',
        'f',
        'H'
    ]

    arch_housekeeping_t = (     # XXX not used XXX
        "id",
        "length",
        "rest_of_packet",
        # Need to finish this
    )

    packet_type = {
        '0x00010000': 'ARCH_FORMAT_VERSION',
        '0x5aa80001': 'ARCH_ID_CONTROL',
        '0x5aa80002': 'ARCH_ID_FIELD_SCALE',
        '0x5aa80003': 'ARCH_ID_RAY_HDR',
        '0x5aa80004': 'ARCH_ID_FILE_HDR',
        '0x5aa80005': 'ARCH_ID_SWEEP_BLOCK',
        '0x5aa50003': 'HSK_ID_PROCESSOR_INFO',
        '0x5aa50001': 'HSK_ID_RADAR_INFO',
        '0x5aa50002': 'HSK_ID_SCAN_SEG'
    }

    _scan_mode_names = ['ppi', 'rhi', 'fixed',
                        'manual ppi', 'manual rhi', 'idle']

##############
# Structures #
##############


def _unpack_structure(string, structure):
    """ Unpack a structure """
    fmt = ''.join([i[1] for i in structure])
    l = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], l))

ARCH_FILE_HDR_T = (
    ('version', 'I'),
    ('creation_version', 'I'),
    ('creator_id', '32s'),
    ('sweep_table_offset', 'Q'),
)

ARCH_RAY_HEADER = (
    ('azimuth', 'f'),
    ('elevation', 'f'),
    ('azimuth_width', 'f'),
    ('elevation_width', 'f'),
    ('gates', 'H'),
    ('beam_index', 'H'),
    ('ns_time', 'I'),
    ('time', 'Q'),
    ('bit_mask', 'Q'),
    ('ray_number', 'I'),
    ('num_pulses', 'I'),
)

FIELD_SCALE_T = (
    ('format', 'i'),
    ('min_val', 'f'),
    ('max_val', 'f'),
    ('bit_mask_pos', 'i'),
    ('type_hint', 'i'),
    ('fld_factor', 'i'),
    ('dat_factor', 'i'),
    ('dat_bias', 'i'),
    ('name', '32s'),
    ('units', '32s'),
    ('descr', '128s')
)

RADAR_INFO_T = (
    ('radar_name', '32s'),
    ('latitude', 'f'),
    ('longitude', 'f'),
    ('altitude', 'f'),
    ('beamwidth', 'f'),
    ('wavelength_cm', 'f'),
    ('gain_ant_h', 'f'),
    ('gain_ant_v', 'f'),
    ('zdr_cal_base', 'f'),
    ('phidp_rot', 'f'),
    ('base_radar_constant', 'f'),
    ('power_measurement_loss_h', 'f'),
    ('power_measurement_loss_v', 'f'),
    ('zdr_cal_base_vhs', 'f'),
    ('test_power_h', 'f'),
    ('test_power_v', 'f'),
    ('dc_loss_h', 'f'),
    ('dc_loss_v', 'f'),
    ('unknown_0', 'f'),
    ('unknown_1', 'f'),
    ('unknown_2', 'f'),
    ('unknown_3', 'f'),
    ('unknown_4', 'f'),
)

PROCESSOR_INFO = (
    ('polarization_mode', 'i'),
    ('processing_mode', 'i'),
    ('pulse_type', 'i'),
    ('test_type', 'i'),
    ('integration_cycle_pulses', 'I'),
    ('clutter_filter_number', 'I'),
    ('range_gate_averaging', 'I'),
    ('indexed_beam_width', 'f'),
    ('gate_spacing', 'f'),
    ('prt_usec', 'f'),
    ('range_start', 'f'),
    ('range_stop', 'f'),
    ('max_gate', 'I'),
    ('test_power', 'f'),
    ('test_pulse_range', 'f'),
    ('test_pulse_length', 'f'),
    ('prt2', 'f'),
    ('range_offset', 'f'),
    ('unknown_0', 'f'),
    ('unknown_1', 'f'),
)

SCAN_SEG = (
    ('az_manual', 'f'),
    ('el_manual', 'f'),
    ('az_start', 'f'),
    ('el_start', 'f'),
    ('scan_rate', 'f'),
    ('segname', '24s'),
    ('opt', 'i'),
    ('follow_mode', 'i'),
    ('scan_type', 'i'),
    ('scan_flags', 'I'),
    ('volume_num', 'I'),
    ('sweep_num', 'I'),
    ('time_limit', 'I'),
    ('webtilt', 'I'),
    ('left_limit', 'f'),
    ('right_limit', 'f'),
    ('up_limit', 'f'),
    ('down_limit', 'f'),
    ('step', 'f'),
    ('max_sweeps', 'I'),
    ('filter_break_sweep', 'I'),
    ('clutter_filter1', 'I'),
    ('clutter_filter2', 'I'),
    ('project', '16s'),
    ('current_fixed_angle', 'f'),
)

# XXX

_FIELD_TABLE = {
    # Chill field name : (Py-ART field name, field long_name attribute)
    'Z': ('DBZ', 'equivalent_reflectivity_factor'),
    'V': ('VEL', 'radial_velocity_of_scatterers_away_from_instrument'),
    'W': ('WIDTH', 'doppler_spectrum_width'),
    'ZDR': ('ZDR', 'log_differential_reflectivity_hv'),
    'LDRH': ('LDRH', 'log_linear_depolarization_ratio_h'),
    'LDRV': ('LDRV', 'log_linear_depolarization_artio_v'),
    '\xce\xa8 DP': ('PHIDP', 'differential_phase_hv'),
    'KDP': ('KDP', 'specific_differential_phase_hv'),
    '\xcf\x81 HV': ('RHOHV', 'cross_correlation_ratio_hv'),
    'NCP': ('NCP', 'normalized_coherent_power'),
    'H Re(lag 1)': ('H Re(lag 1)', 'real_part_of_lag_1_correlation_h'),
    'V Re(lag 2)': ('V Re(lag 2)', 'real_part_of_lag_2_correlation_v'),
    'VAvgQ': ('VAvgQ', 'v_average_quadrature'),
    'V Im(lag 1)': ('V Im(lag 1)', 'imaginary_part_of_v_at_lag_1'),
    'HAvgQ': ('HAvgQ', 'h_average_quadrature'),
    'H Im(lag 2)': ('H Im(lag 2)', 'imaginary_part_lag_2_correlation_h'),
    'V lag 0': ('V lag 0', 'absolute_value_of_lag_0_correlation_v'),
    'H lag 0': ('H lag 0', 'absolute_value_of_lag_0_correlation_h'),
    'H lag 0 cx':
    ('H lag 0 cx', 'absolute_value_of_lag_0_cross_correlation_h'),
    'H Im(lag 1)':
    ('H Im(lag 1)', 'imaginary_part_of_lag_1_correlation_h'),
    'H Re(lag 2)':
    ('H Re(lag 2)', 'real_part_of_lag_2_correlation_h'),
    'V lag 0 cx':
    ('V lag 0 cx', 'absolute_value_of_lag_0_cross_correlation_v'),
    'V Re(lag 1)': ('V Re(lag 1)', 'real_part_of_lag_1_correlation_v'),
    'V Im(lag 2)':
    ('V Im(lag 2)', 'imaginary_part_of_lag_2_correlation_v'),
    'HV lag 0 I':
    ('HV lag 0 I', 'real_part_of_cross_channel_correlation_at_lag_0'),
    'HV lag 0 Q':
    ('HV lag 0 Q',
     'imaginary_part_of_cross_channel_correlation_at_lag_0'),
    'VAvgI': ('VAvgI', 'v_average_inphase'),
    'HAvgI': ('HAvgI', 'h_average_inphase'),
    '\xcf\x81 HCX': ('RHOHCX', 'lag_0_h_co_to_cross_correlation'),
    '\xcf\x81 VCX': ('RHOVCX', 'lag_0_v_co_to_cross_correlation'),
}
CHILL_FIELD_MAPPING = dict((k, v[0]) for k, v in _FIELD_TABLE.items())
CHILL_FIELD_LONG_NAME = dict((k, v[1]) for k, v in _FIELD_TABLE.items())
