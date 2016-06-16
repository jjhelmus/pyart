
from ..config import get_metadata
from ..lazydict import LazyLoadDict
from .transforms import antenna_vectors_to_cartesian, cartesian_to_geographic


_TIME_DIMENTIONED_INSTRUMENT_PARAMETERS = [
    'pulse_width', 'prt', 'prt_ratio', 'nyquist_velocity',
    'unambiguous_range', 'n_samples', 'sampling_ratio']

# TODO:
# * detach sweep from radar
# * make iterable
# * make Radar produce sweep object during iteration
# * gate_edges?
# * refactor display classes to use this object?


class Sweep(object):
    """
    A class storing data from a single radar sweep in antenna coordinates.

    Parameters
    ----------
    radar : Radar
        Radar instance from which sweep is derived.
    sweep_num : int
        Number of sweep in radar to extract.
    instr_param_keys : list, optional
        List of entries in the radar's instrument_parameters attribute to
        include in the sweep's instrument_parameters attribute. These entries
        should have a leading 'sweep' dimensions. A value of None, the default,
        will use a predefined list of instrument parameter which are typically
        indexed by sweep.

    Attributes
    ----------
    sweep_number : int
        The number of the sweep, 0-based.
    sweep_mode : str
        Sweep mode for the sweep.
    fixed_angle : float
        Target angle for the sweep. Azimuth angle in RHI sweeps, elevation
        angle in all other modes.
    target_scan_rate : float or None
        Intended scan rate for each sweep.  If not provided this attribute is
        set to None, indicating this parameter is not available.
    rays_are_indexed : str or None
        Indication of whether ray angles are indexed to a regular grid in
        each sweep.  If not provided this attribute is set to None, indicating
        ray angle spacing is not determined.
    ray_angle_res : float or None
        If rays_are_indexed is not None, this provides the angular resolution
        of the grid.  If not provided or available this attribute is set to
        None.
    ngates : int
        Number of gates (bins) in the sweep, cannot be changed.
    nrays : int
        Number of rays in the volume, cannot be changed.

    time : dict
        Time at the center of each ray.
    range : dict
        Range to the center of each gate (bin).
    fields : dict of dicts
        Moment fields.
    metadata : dict
        Metadata describing the instrument and data.
    scan_type : str
        Type of scan, one of 'ppi', 'rhi', 'sector' or 'other'.  If the scan
        volume contains multiple sweep modes this should be 'other'.
    latitude : dict
        Latitude of the instrument.
    longitude : dict
        Longitude of the instrument.
    altitude : dict
        Altitude of the instrument, above sea level.
    altitude_agl : dict or None
        Altitude of the instrument above ground level.  If not provided this
        attribute is set to None, indicating this parameter not available.
    azimuth : dict
        Azimuth of antenna, relative to true North.
    elevation : dict
        Elevation of antenna, relative to the horizontal plane.
    gate_x, gate_y, gate_z : LazyLoadDict
        Location of each gate in a Cartesian coordinate system assuming a
        standard atmosphere with a 4/3 Earth's radius model. The data keys of
        these attributes are create upon first access from the data in the
        range, azimuth and elevation attributes. If these attributes are
        changed use :py:func:`init_gate_x_y_z` to reset.
    gate_longitude, gate_latitude : LazyLoadDict
        Geographic location of each gate.  The projection parameter(s) defined
        in the `projection` attribute are used to perform an inverse map
        projection from the Cartesian gate locations relative to the radar
        location to longitudes and latitudes. If these attributes are changed
        use :py:func:`init_gate_longitude_latitude` to reset the attributes.
    projection : dic or str
        Projection parameters defining the map projection used to transform
        from Cartesian to geographic coordinates.  The default dictionary sets
        the 'proj' key to 'pyart_aeqd' indicating that the native Py-ART
        azimuthal equidistant projection is used. This can be modified to
        specify a valid pyproj.Proj projparams dictionary or string.
        The special key '_include_lon_0_lat_0' is removed when interpreting
        this dictionary. If this key is present and set to True, which is
        required when proj='pyart_aeqd', then the radar longitude and
        latitude will be added to the dictionary as 'lon_0' and 'lat_0'.
    gate_altitude : LazyLoadDict
        The altitude of each radar gate as calculated from the altitude of the
        radar and the Cartesian z location of each gate.  If this attribute
        is changed use :py:func:`init_gate_altitude` to reset the attribute.
    scan_rate : dict or None
        Actual antenna scan rate.  If not provided this attribute is set to
        None, indicating this parameter is not available.
    antenna_transition : dict or None
        Flag indicating if the antenna is in transition, 1 = yes, 0 = no.
        If not provided this attribute is set to None, indicating this
        parameter is not available.
    instrument_parameters : dict of dicts or None
        Instrument parameters, if not provided this attribute is set to None,
        indicating these parameters are not avaiable.  This dictionary also
        includes variables in the radar_parameters CF/Radial subconvention.
    radar_calibration : dict of dicts or None
        Instrument calibration parameters.  If not provided this attribute is
        set to None, indicating these parameters are not available

    """

    def __init__(self, radar, sweep_num, instr_param_keys=None):
        """
        """

        self._radar = radar
        self._sweep_num = sweep_num

        # no change from complete volume
        self.range = radar.range
        self.metadata = radar.metadata
        self.radar_calibration = radar.radar_calibration

        # locations are assumed fixed (not a moving platform)
        self.latitude = radar.latitude
        self.longitude = radar.longitude
        self.altitude = radar.altitude
        self.altitude_agl = radar.altitude_agl
        self.projection = radar.projection

        # slice data per ray
        # TODO make these linked to the radar dictionaries
        sweep_slice = radar.get_slice(sweep_num)

        self.time = _dict_from_dict(radar.time, sweep_slice)
        self.azimuth = _dict_from_dict(radar.azimuth, sweep_slice)
        self.elevation = _dict_from_dict(radar.elevation, sweep_slice)
        self.scan_rate = _dict_from_dict(radar.scan_rate, sweep_slice)
        self.antenna_transition = _dict_from_dict(
            radar.antenna_transition, sweep_slice)

        self.rotation = _dict_from_dict(radar.rotation, sweep_slice)
        self.tilt = _dict_from_dict(radar.tilt, sweep_slice)
        self.roll = _dict_from_dict(radar.roll, sweep_slice)
        self.drift = _dict_from_dict(radar. drift, sweep_slice)
        self.heading = _dict_from_dict(radar.heading, sweep_slice)
        self.pitch = _dict_from_dict(radar.pitch, sweep_slice)
        self.georefs_applied = _dict_from_dict(
            radar.georefs_applied, sweep_slice)

        self.fields = {}
        for field_name, field_dict in radar.fields.items():
            self.fields[field_name] = _dict_from_dict(field_dict, sweep_slice)

        if radar.instrument_parameters is None:
            self.instrument_parameters = None
        else:
            instr_params = {}
            if instr_param_keys is None:
                instr_param_keys = _TIME_DIMENTIONED_INSTRUMENT_PARAMETERS
            for key in instr_param_keys:
                if key not in radar.instrument_parameters:
                    continue
                ip_dict = radar.instrument_parameters[key]
                instr_params[key] = _dict_from_dict(ip_dict, sweep_slice)
            self.instrument_parameters = instr_params

        # initalize attributes with lazy load dictionaries
        self.init_gate_x_y_z()
        self.init_gate_longitude_latitude()
        self.init_gate_altitude()

        return

    # read only properties

    @property
    def ngates(self):
        return self._radar.ngates

    @property
    def nrays(self):
        return self._radar.rays_per_sweep['data'][self._sweep_num]

    # single values extracted/set from the underlying radar instance

    @property
    def sweep_number(self):
        return int(self._radar.sweep_number['data'][self._sweep_num])

    @sweep_number.setter
    def sweep_number(self, value):
        self._radar.sweep_number['data'][self._sweep_num] = value

    @property
    def sweep_mode(self):
        return str(self._radar.sweep_mode['data'][self._sweep_num])

    @sweep_mode.setter
    def sweep_mode(self, value):
        self._radar.sweep_mode['data'][self._sweep_num] = value

    @property
    def fixed_angle(self):
        return float(self._radar.fixed_angle['data'][self._sweep_num])

    @fixed_angle.setter
    def fixed_angle(self, value):
        self._radar.fixed_angle['data'][self._sweep_num] = value

    @property
    def target_scan_rate(self):
        if self._radar.target_scan_rate is None:
            return None
        return float(self._radar.target_scan_rate['data'][self._sweep_num])

    @target_scan_rate.setter
    def target_scan_rate(self, value):
        if self._radar.target_scan_rate is None:
            raise AttributeError("Radar attribute is None")
        self._radar.target_scan_rate['data'][self._sweep_num] = value

    @property
    def rays_are_indexed(self):
        if self._radar.rays_are_indexed is None:
            return None
        return str(self._radar.rays_are_indexed['data'][self._sweep_num])

    @rays_are_indexed.setter
    def rays_are_indexed(self, value):
        if self._radar.rays_are_indexed is None:
            raise AttributeError("Radar attribute is None")
        self._radar.rays_are_indexed['data'][self._sweep_num] = value

    @property
    def ray_angle_res(self):
        if self._radar.ray_angle_res is None:
            return None
        return float(self._radar.ray_angle_res['data'][self._sweep_num])

    @ray_angle_res.setter
    def ray_angle_res(self, value):
        if self._radar.ray_angle_res is None:
            raise AttributeError("Radar attribute is None")
        self._radar.ray_angle_res['data'][self._sweep_num] = value

    def __getstate__(self):
        """ Return object's state which can be pickled. """
        state = self.__dict__.copy()  # copy the objects state
        # Remove unpicklable entries (those which are lazily loaded)
        del state['gate_x']
        del state['gate_y']
        del state['gate_z']
        del state['gate_longitude']
        del state['gate_latitude']
        del state['gate_altitude']
        return state

    def __setstate__(self, state):
        """ Restore unpicklable entries from pickled object. """
        self.__dict__.update(state)
        self.init_gate_x_y_z()
        self.init_gate_longitude_latitude()
        self.init_gate_altitude()

    # Attribute init/reset method
    def init_gate_x_y_z(self):
        """ Initialize or reset the gate_{x, y, z} attributes. """
        gate_x = LazyLoadDict(get_metadata('gate_x'))
        gate_x.set_lazy('data', _gate_data_factory(self, 0))
        self.gate_x = gate_x

        gate_y = LazyLoadDict(get_metadata('gate_y'))
        gate_y.set_lazy('data', _gate_data_factory(self, 1))
        self.gate_y = gate_y

        gate_z = LazyLoadDict(get_metadata('gate_z'))
        gate_z.set_lazy('data', _gate_data_factory(self, 2))
        self.gate_z = gate_z

    def init_gate_longitude_latitude(self):
        """
        Initialize or reset the gate_longitude and gate_latitude attributes.
        """
        gate_longitude = LazyLoadDict(get_metadata('gate_longitude'))
        gate_longitude.set_lazy('data', _gate_lon_lat_data_factory(self, 0))
        self.gate_longitude = gate_longitude

        gate_latitude = LazyLoadDict(get_metadata('gate_latitude'))
        gate_latitude.set_lazy('data', _gate_lon_lat_data_factory(self, 1))
        self.gate_latitude = gate_latitude

    def init_gate_altitude(self):
        """ Initialize the gate_altitude attribute. """
        gate_altitude = LazyLoadDict(get_metadata('gate_altitude'))
        gate_altitude.set_lazy('data', _gate_altitude_data_factory(self))
        self.gate_altitude = gate_altitude


def _dict_from_dict(radar_dict, sweep_slice):
    if radar_dict is None:
        return None
    sweep_dict = radar_dict.copy()
    sweep_dict['data'] = radar_dict['data'][sweep_slice]
    return sweep_dict


def _gate_data_factory(radar, coordinate):
    """ Return a function which returns the Cartesian locations of gates. """
    def _gate_data():
        """ The function which returns the Cartesian locations of gates. """
        ranges = radar.range['data']
        azimuths = radar.azimuth['data']
        elevations = radar.elevation['data']
        cartesian_coords = antenna_vectors_to_cartesian(
            ranges, azimuths, elevations, edges=False)
        # load x, y, and z data except for the coordinate in question
        if coordinate != 0:
            radar.gate_x['data'] = cartesian_coords[0]
        if coordinate != 1:
            radar.gate_y['data'] = cartesian_coords[1]
        if coordinate != 2:
            radar.gate_z['data'] = cartesian_coords[2]
        return cartesian_coords[coordinate]
    return _gate_data


def _gate_lon_lat_data_factory(radar, coordinate):
    """ Return a function which returns the geographic locations of gates. """
    def _gate_lon_lat_data():
        """ The function which returns the geographic locations gates. """
        x = radar.gate_x['data']
        y = radar.gate_y['data']
        projparams = radar.projection.copy()
        if projparams.pop('_include_lon_0_lat_0', False):
            projparams['lon_0'] = radar.longitude['data'][0]
            projparams['lat_0'] = radar.latitude['data'][0]
        geographic_coords = cartesian_to_geographic(x, y, projparams)
        # set the other geographic coordinate
        if coordinate == 0:
            radar.gate_latitude['data'] = geographic_coords[1]
        else:
            radar.gate_longitude['data'] = geographic_coords[0]
        return geographic_coords[coordinate]
    return _gate_lon_lat_data


def _gate_altitude_data_factory(radar):
    """ Return a function which returns the gate altitudes. """
    def _gate_altitude_data():
        """ The function which returns the gate altitudes. """
        return radar.altitude['data'] + radar.gate_z['data']
    return _gate_altitude_data
