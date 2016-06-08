
import netCDF4
import numpy as np

from ..config import get_metadata
from ..lazydict import LazyLoadDict
from .transforms import antenna_vectors_to_cartesian, cartesian_to_geographic


class Ray:
    """
    A class storing data from a single radar ray in antenna coordinates.

    Attributes
    ----------
    time : datetime
        Time at the center of the ray.
    range : dict
        Range to the center of each gate (bin).
    fields : dict of dicts
        Moment fields.
    metadata : dict
        Metadata describing the instrument and data.
    latitude : dict
        Latitude of the instrument.
    longitude : dict
        Longitude of the instrument.
    altitude : dict
        Altitude of the instrument, above sea level.
    altitude_agl : dict or None
        Altitude of the instrument above ground level.  If not provided this
        attribute is set to None, indicating this parameter not available.
    azimuth : float
        Azimuth of antenna, relative to true North.
    elevation : float
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
    scan_rate : float or None
        Actual antenna scan rate.  If not provided this attribute is set to
        None, indicating this parameter is not available.
    antenna_transition : int or None
        Flag indicating if the antenna is in transition, 1 = yes, 0 = no.
        If not provided this attribute is set to None, indicating this
        parameter is not available.
    instrument_parameters : dict of numerical values or None
        Instrument parameters, if not provided this attribute is set to None,
        indicating these parameters are not avaiable.  This dictionary also
        includes variables in the radar_parameters CF/Radial subconvention.
    radar_calibration : dict of dicts or None
        Instrument calibration parameters.  If not provided this attribute is
        set to None, indicating these parameters are not available
    ngates : int
        Number of gates (bins) in the ray.

    """

    def __init__(self, sweep, ray_number, copy=False):
        """
        """
        # TODO copy

        # sizes
        self.ngates = sweep.ngates

        # no change from sweep
        self.range = sweep.range
        self.metadata = sweep.metadata
        self.radar_calibration = sweep.radar_calibration

        # locations are assumed fixed (not a moving platform)
        self.latitude = sweep.latitude
        self.longitude = sweep.longitude
        self.altitude = sweep.altitude
        self.altitude_agl = sweep.altitude_agl
        self.projection = sweep.projection

        # single values
        # TODO make these properties?
        self.time = netCDF4.num2date(
            sweep.time['data'][ray_number], sweep.time['units'])

        self.azimuth = _none_or_float(sweep.azimuth, ray_number)
        self.elevation = _none_or_float(sweep.elevation, ray_number)
        self.scan_rate = _none_or_float(sweep.scan_rate, ray_number)

        self.antenna_transition = _none_or_int(
            sweep.antenna_transition, ray_number)

        self.rotation = _none_or_float(sweep.rotation, ray_number)
        self.tilt = _none_or_float(sweep.tilt, ray_number)
        self.roll = _none_or_float(sweep.roll, ray_number)
        self.drift = _none_or_float(sweep. drift, ray_number)
        self.heading = _none_or_float(sweep.heading, ray_number)
        self.pitch = _none_or_float(sweep.pitch, ray_number)
        self.georefs_applied = _none_or_int(sweep.georefs_applied, ray_number)

        sweep_ip = sweep.instrument_parameters
        if sweep_ip is None:
            self.instrument_parameters = None
        else:
            ray_ip = {k: v['data'][0] for k, v in sweep_ip.items()}
            self.instrument_parameters = ray_ip

        self.fields = {}
        for field_name, field_dict in sweep.fields.items():
            # TODO link these dictionaries to the sweep dictionary
            ray_field_dict = field_dict.copy()
            ray_field_dict['data'] = field_dict['data'][ray_number]
            self.fields[field_name] = ray_field_dict

        # initalize attributes with lazy load dictionaries
        self.init_gate_x_y_z()
        self.init_gate_longitude_latitude()
        self.init_gate_altitude()

        return

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


def _none_or_float(sweep_dict, ray_number):
    if sweep_dict is None:
        return None
    else:
        return float(sweep_dict['data'][ray_number])


def _none_or_int(sweep_dict, ray_number):
    if sweep_dict is None:
        return None
    else:
        return int(sweep_dict['data'][ray_number])


def _gate_data_factory(ray, coordinate):
    """ Return a function which returns the Cartesian locations of gates. """
    def _gate_data():
        """ The function which returns the Cartesian locations of gates. """
        ranges = ray.range['data']
        azimuths = ray.azimuth
        elevations = ray.elevation
        cartesian_coords = antenna_vectors_to_cartesian(
            ranges, azimuths, elevations, edges=False)
        # load x, y, and z data except for the coordinate in question
        if coordinate != 0:
            ray.gate_x['data'] = np.squeeze(cartesian_coords[0])
        if coordinate != 1:
            ray.gate_y['data'] = np.squeeze(cartesian_coords[1])
        if coordinate != 2:
            ray.gate_z['data'] = np.squeeze(cartesian_coords[2])
        return np.squeeze(cartesian_coords[coordinate])
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
