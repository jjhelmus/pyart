"""
pyart.retrieve.winds
====================

"""

import time
import numpy as np
import collections

from datetime import datetime
from copy import deepcopy
from scipy.optimize import minimize
from mpl_toolkits.basemap import pyproj
from netCDF4 import num2date

from ..io import Grid
from ..util.datetime_utils import datetime_from_grid, datetimes_from_dataset
from ..config import get_fillvalue, get_field_name, get_metadata
from ..retrieve import cga, divergence, continuity, gradient

class Cost(object):

    __init__(self)

class Jacobian(object):

    __init__(self)


def _radar_coverage(grids, fill_value=None, refl_field=None, vel_field=None):
    """
    Parameters
    ----------
    grids : list
        List of all available radar grids used to determine the
        coverage within the analysis domain.
        
    Optional parameters
    -------------------
    refl_field, vel_field : str
        Name of fields which will be used to determine coverage. A value of
        None will use the default field name as defined in the Py-ART
        configuration file.
    
    Returns
    -------
    cover : dict
        Radar coverage data dictionary.
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
    
    # Initialize coverage array
    cover = np.zeros_like(grids[0].fields[vel_field]['data'])
    
    # Loop over all grids to get the total number of observations at
    # each grid point in the analysis domain
    for grid in grids:
        # Get data
        ze = np.ma.filled(grid.fields[refl_field]['data'], fill_value)
        vr = np.ma.filled(grid.fields[vel_field]['data'], fill_value)
        
        # Create appropriate boolean arrays
        is_bad_refl = ze == fill_value
        is_bad_vel = vr == fill_value
        
        has_obs = ~np.logical_or(is_bad_refl, is_bad_vel)
        cover = np.where(has_obs, cover + 1, cover)
        
    # Return dictionary of results
    return {'data': cover.astype(np.int32),
            'standard_name': 'radar_coverage',
            'long_name': 'Number of radar observations',
            'valid_min': 0,
            'valid_max': len(grids)}
    

def _radar_components(grids, proj='lcc', datum='NAD83', ellps='GRS80'):
    """
    Add a Cartesian components field to Grid objects.
    
    Parameters
    ----------
    grids : list
        List of all available radar grids which will have a Cartesian
        components field added to their existing fields.
        
    Optional parameters
    -------------------
    proj : str
    
    datum : str
    
    ellps : str
        
    Returns
    -------
    grids : list
        List of grids with added Cartesian components fields.
    """
    
    # Get analysis domain axes
    x = grids[0].axes['x_disp']['data']
    y = grids[0].axes['y_disp']['data']
    z = grids[0].axes['z_disp']['data']
    
    # Loop over all grids
    for i, grid in enumerate(grids):
        # Get latitude and longitude of analysis domain origin and the
        # latitude and longitude of the current grid (radar)
        lat_0 = grid.axes['lat']['data'][0]
        lon_0 = grid.axes['lon']['data'][0]
        lat_r = grid.metadata['radar_0_lat']
        lon_r = grid.metadata['radar_0_lon']
        
        # Create map projection centered at the analysis domain origin
        pj = pyproj.Proj(proj=proj, datum=datum, ellps=ellps, lat_0=lat_0,
                           lon_0=lon_0, x_0=0.0, y_0=0.0)
        
        # Get the (x, y) location of the radar in the analysis domain from
        # the projection
        x_r, y_r = pj(lon_r, lat_r)
        
        # Create the grid mesh which has an origin at the radar location
        Z, Y, X = np.meshgrid(z, y-y_r, x-x_r, indexing='ij')
        R = np.sqrt(X**2 + Y**2 + Z**2)
        
        # Compute Cartesian components for the current grid
        ic = X / R
        jc = Y / R
        kc = Z / R
        
        # Create dictionaries of results
        ic = {'data': ic.astype(np.float64),
              'standard_name': 'x_component',
              'long_name': 'Eastward component',
              'valid_min': -1.0,
              'valid_max': 1.0}
        jc = {'data': jc.astype(np.float64),
              'standard_name': 'y_component',
              'long_name': 'Northward component',
              'valid_min': -1.0,
              'valid_max': 1.0}
        kc = {'data': kc.astype(np.float64),
              'standard_name': 'z_component',
              'long_name': 'Vertical component',
              'valid_min': -1.0,
              'valid_max': 1.0}
              
        # Add new fields to current grid and update grids list
        grid.fields.update({'x_component': ic})
        grid.fields.update({'y_component': jc})
        grid.fields.update({'z_component': kc})
        
    return


def _echo_bounds(grids, mds=0.0, min_layer=1500.0, top_offset=500.0,
                  proc=1, fill_value=None, refl_field=None):
    """
    Determine the echo base and top heights from the large-scale
    coverage of the radar network.
    
    Parameters
    ----------
    grids : list
    
    
    Optional Parameters
    -------------------
    mds : float
        Minimum detectable signal in dBZ used to define the noise cut-off.
    min_layer : float
        Minimum cloud layer allowed in meters.
    top_offset : float
        Value in meters added to the echo top height. See Protat and Zawadzki
        (1999) for more information.
    proc : int
        Number of processors requested.
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.
    refl_field : str
        Name of reflectivity field which will be used to estimate the fall
        speed. A value of None will use the default field name as defined in
        the Py-ART configuration file.

    References
    ----------
    Protat, A. and I. Zawadzki, 1999:

    Returns
    -------
    base, top : dict
        Echo base and top data dictionaries.
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    
    # Get analysis domain axes
    z = grids[0].axes['z_disp']['data']
    
    # Compute the maximum reflectivity observed at each grid point
    # using all grids (radars)
    ze = np.ma.max([grid.fields[refl_field]['data'] for
                    grid in grids], axis=0)
    ze = np.ma.filled(ze, fill_value).astype(np.float64)
        
    # Estimate echo base and top heights, and add offset to echo top heights
    base, top = continuity.boundary_conditions(ze, z, mds=mds,
                                min_layer=min_layer, proc=proc,
                                fill_value=fill_value)
    
    top = np.where(top != fill_value, top + top_offset, top)
    
    base = np.ma.masked_equal(base, fill_value)
    top = np.ma.masked_equal(top, fill_value)
        
    # Create dictionaries of results
    base = {'data': base,
            'standard_name': 'echo_base_height',
            'long_name': 'Height of echo base',
            'valid_min': z.min(),
            'valid_max': z.max(),
            '_FillValue': base.fill_value,
            'units': 'meters'}
    top = {'data': top,
           'standard_name': 'echo_top_height',
           'long_name': 'Height of echo top',
           'valid_min': z.min(),
           'valid_max': z.max(),
           '_FillValue': top.fill_value,
           'units': 'meters'}
        
    return base, top

            
def _observation_weight(grids, wgt_o=1.0, fill_value=None,
                        refl_field=None, vel_field=None):
    """
    Add an observation weight field to Grid objects. Grid points
    with valid observations should be given a scalar weight greater
    than 0, and grid points with missing data should be given a scalar
    weight of 0.
    
    Parameters
    ----------
    grids : list
        List of all available grids which will have an observation
        weight field added to their existing fields.
    
    Optional parameters
    -------------------
    wgt_o : float
        Observation weight used at each grid point with valid observations.
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.
    refl_field, vel_field : str
        Name of fields which will be used to determine grid points with
        valid observations. A value of None will use the default field name
        as defined in the Py-ART configuration file.
        
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
        
    # Loop over all grids
    for grid in grids:
        # Get appropriate data
        ze = np.ma.filled(grid.fields[refl_field]['data'], fill_value)
        vr = np.ma.filled(grid.fields[vel_field]['data'], fill_value)
        
        # Initialize observation weight array for each grid (radar)
        grid_weight = np.zeros_like(ze)
        
        # Create appropriate boolean arrays
        is_bad_refl = ze == fill_value
        is_bad_vel = vr == fill_value
        
        is_good = ~np.logical_or(is_bad_refl, is_bad_vel)
        
        grid_weight[is_good] = wgt_o
        
        # Create dictionary of results
        grid_weight = {'data': grid_weight,
                       'standard_name': 'observation_weight',
                       'long_name': 'Radar observation weight',
                       'valid_min': 0.0,
                       'valid_max': wgt_o}
                 
        # Add new field to grid object
        grid.fields.update({'observation_weight': grid_weight})

    return
    
            
def _radar_qc(grids, mds=0.0, vel_max=55.0, vel_grad_max=10.0, ncp_min=0.3,
              rhv_min=0.7, window_size=5, noise_ratio=50.0, fill_value=None,
              refl_field=None, vel_field=None, ncp_field=None,
              rhv_field=None):
    """
    Parameters
    ----------
    grids : list
        
    
    Optional parameters
    -------------------
    mds : float
        Minimum detectable signal in dBZ.
    vel_max : float
        Maximum absolute radial velocity allowed in m/s.
    vel_grad_max : float
        Maximum velocity gradient magnitude allowed in m/s.
    ncp_min, rhv_min : float
        Minimum values allowed for normalized coherent power and correlation
        coefficient.
    window_size : int
    noise_ratio : float
    fill_value : float
    
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
        
    # Get analysis domain dimensions
    nz, ny, nx = grids[0].fields[refl_field]['data'].shape
    
    # Get analysis domain heights and create its mesh
    z = grids[0].axes['z_disp']['data']
    Z = np.repeat(z, ny*nx, axis=0).reshape(nz, ny, nx)
    
    # Loop over all grids
    for grid in grids:
        # Compute the magnitude of the radial velocity gradient and determine
        # where it is bad
        vr = grid.fields[vel_field]['data']
        dvrz, dvry, dvrx = np.gradient(vr)
        grad_mag = np.ma.sqrt(dvrx**2 + dvry**2 + dvrz**2)
        grad_mag = np.ma.filled(grad_mag, fill_value)
        is_bad_vel_grad = np.logical_or(grad_mag > vel_grad_max,
                                        grad_mag == fill_value)
        
        # Fill data
        ze = np.ma.filled(grid.fields[refl_field]['data'], fill_value)
        vr = np.ma.filled(grid.fields[vel_field]['data'], fill_value)
        ncp = np.ma.filled(grid.fields[ncp_field]['data'], fill_value)
        rhv = np.ma.filled(grid.fields[rhv_field]['data'], fill_value)
        
        # Create appropriate boolean arrays
        is_noise = np.logical_or(ze < mds, ze == fill_value)
        is_high_vel = np.logical_or(np.abs(vr) > vel_max, vr == fill_value)
        is_vel_artifact = np.logical_or(is_high_vel, is_bad_vel_grad)
        is_bad_ncp = np.logical_or(ncp < ncp_min, ncp == fill_value)
        is_bad_rhv = np.logical_or(rhv < rhv_min, rhv == fill_value)
        is_non_meteo = np.logical_or(is_bad_ncp, is_bad_rhv)
        
        is_bad_refl = np.logical_or(is_noise, is_non_meteo)
        is_bad_vel = np.logical_or(is_vel_artifact, is_non_meteo)
        
        # Special attention to heights below 2000 m where velocity artifacts
        # can especially be a problem
        is_bad_vel[np.logical_and(np.abs(vr) > 30.0, Z < 2000.0)] = True
        
        # Update masks
        grid.fields[refl_field]['data'] = np.ma.masked_where(is_bad_refl, ze)
        grid.fields[vel_field]['data'] = np.ma.masked_where(is_bad_vel, vr)
        
        # Despeckle reflectivity and Doppler velocity fields
        grid.despeckle_field(refl_field, window_size=window_size,
                             noise_ratio=noise_ratio,
                             fill_value=fill_value)
        grid.despeckle_field(vel_field, window_size=window_size,
                             noise_ratio=noise_ratio,
                             fill_value=fill_value)
        
    return

    
def _column_types(cover, base, top, fill_value=None):
    """
    Parameters
    ----------
    cover : dict
    base, top : dict
    
    Optional parameters
    -------------------
    
    Returns
    -------
    column : dict
        Column classification dictionary.
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
        
    # Get data
    cover = cover['data'].astype(np.int32)
    base = base['data'].astype(np.float64)
    top = top['data'].astype(np.float64)
    
    # Get the column types using the radar coverage and echo base height
    column = continuity.column_type(cover, base, fill_value=fill_value)
     
    return {'data': column.astype(np.int32),
            'standard_name': 'column_type',
            'long_name': 'Column classifications',
            'valid_min': 0,
            'valid_max': 5,
            'comment': ('0 = Undefined, 1 = Well-defined, '
                        '2 = Top-defined, 3 = Anvil-like, '
                        '4 = Transition-like, '
                        '5 = Discontinuous')}


def _arm_interp_sonde(grid, sonde, target=None, rho0=1.2, H=10000.0,
                      standard_density=False, finite_scheme='basic',
                      fill_value=None, debug=False, verbose=False):
    """
    Parameters
    ----------
    grid : Grid
    
    Optional parameters
    -------------------
    target : datetime
    rho0, H : float
        Reference density in kg m^-3 and scale height in meters. Only
        applicable if 'standard_density' is True.
    finite_scheme : 'basic' or 'high-order'
        Finite difference scheme to use when calculating density gradient.
        Only applicable if 'standard' is False.
    standard_density : bool
        If True, the returned density profile is from a standard atmosphere.
        False uses the sounding data.
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.
    debug : bool
        Print debugging information
    verbose : bool
       If True print as much information as possible.
    
    Returns
    -------
    T : np.ndarray
        Temperature profile in degrees Celsius.
    P : np.ndarray
        Pressure profile in kPa.
    rho : np.ndarray
        Density profile in kg m^-3.
    drho : np.ndarray
        Profile of the rate of change of density with respect to height
        in kg m^-4.
    u, v : np.ndarray
        Eastward and northward wind component profiles in meters per second.
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Get axes. The heights given in the interpolated or merged sounding 
    # product are given in kilometers above mean sea level
    z_grid = grid.axes['z_disp']['data'] + grid.axes['alt']['data'] # in (m)
    z_sonde = 1000.0 * sonde.variables['height'][:] # in (m)
    
    # Get closest time index in sounding to the target time
    if target is None:
        target = datetime_from_grid(grid)
    dt_sonde = datetimes_from_dataset(sonde)
    t = np.abs(dt_sonde - target).argmin()
    if verbose:
        print 'Target time is %s' %target
        print 'Closest merged sounding time to target is %s' %dt_sonde[t]
    
    # Get data from sounding
    temp = sonde.variables['temp'][t,:] # (C)
    pres = sonde.variables['bar_pres'][t,:] # (kPa)
    u = sonde.variables['u_wind'][t,:] # (m/s)
    v = sonde.variables['v_wind'][t,:] # (m/s)
    
    # Now compute the density of dry air and its first derivative with
    # respect to height
    if standard_density:
        rho = rho0 * np.exp(-z_sonde / H)
        drho = -(rho0 / H) * np.exp(-z_sonde / H)
        
    else:
        R = 287.058 # (J K^-1 mol^1)
        rho = (pres * 1000.0) / (R * (temp + 273.15)) # (kg m^-3)
        drho = gradient.density1d(rho, z_sonde, fill_value=fill_value,
                                  finite_scheme=finite_scheme)
    
    # Interpolate sounding data to vertical dimension of grid
    temp = np.interp(z_grid, z_sonde, temp, left=None, right=None)
    pres = np.interp(z_grid, z_sonde, pres, left=None, right=None)
    rho = np.interp(z_grid, z_sonde, rho, left=None, right=None)
    drho = np.interp(z_grid, z_sonde, drho, left=None, right=None)
    u = np.interp(z_grid, z_sonde, u, left=None, right=None)
    v = np.interp(z_grid, z_sonde, v, left=None, right=None)
    
    if debug:
        print 'Minimum air temperature is %.2f C' %temp.min()
        print 'Maximum air temperature is %.2f C' %temp.max()
        print 'Minimum air pressure is %.2f kPa' %pres.min()
        print 'Maximum air pressure is %.2f kPa' %pres.max()
        print 'Minimum air density is %.3f kg/m^3' %rho.min()
        print 'Maximum air density is %.3f kg/m^3' %rho.max()
        print 'Minimum eastward wind component is %.2f m/s' %u.min()
        print 'Maximum eastward wind component is %.2f m/s' %u.max()
        print 'Minimum northward wind component is %.2f m/s' %v.min()
        print 'Maximum northward wind component is %.2f m/s' %v.max()
    
    return temp, pres, rho, drho, u, v


def _standard_atmosphere(grid, rho0=1.2, H=10000.0):
    """
    Parameters
    ----------
    grid : Grid

    Optional parameters
    -------------------
    rho0 : float
        Reference density in kg m^-3.
    H : float
        Scale height in meters.

    Returns
    -------
    rho : np.ndarray
        Density profile in kg m^-3.
    drho : np.ndarray
        Profile of the rate of change of density with respect to height
        in kg m^-4.
    """

    # Get analysis domain heights
    z = grid.axes['z_disp']['data']
    
    # Compute the air density and its first derivative with respect to height
    rho = rho0 * np.exp(-z / H)
    drho = -(rho0 / H) * np.exp(-z / H)

    return rho, drho

        
def _fall_speed_caya(grids, temp, fill_value=None, refl_field=None):
    """
    Parameters
    ----------
    grids : list
        List of all available grids which will have a hydrometeor fall
        speed field added to their existing fields.
    temp : np.ndarray
        Temperature profile in degrees Celsius.
    
    Optional parameters
    -------------------
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.
    refl_field : str
        Name of reflectivity field which will be used to estimate the fall
        speed. A value of None will use the default field name as defined in
        the Py-ART configuration file.
        
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    
    # Get analysis domain dimensions
    nz, ny, nx = grids[0].fields[refl_field]['data'].shape
    
    # Get analysis domain heights and create its mesh
    z = grids[0].axes['z_disp']['data']
    Z = np.repeat(z, ny*nx, axis=0).reshape(nz, ny, nx)
    
    # Create the temperature mesh
    T = np.repeat(temp, ny*nx, axis=0).reshape(nz, ny, nx)
    
    # Compute the maximum reflectivity observed at each grid point
    # using all grids (radars)
    ze = np.ma.max([grid.fields[refl_field]['data'] for
                    grid in grids], axis=0)
    
    # Compute the precipitation concentration
    M = np.ma.exp((ze - 43.1) / 7.6)
    
    # Define liquid and ice relations
    liquid = lambda M: -5.94 * M**(1.0 / 8.0) * np.exp(Z / 20000.0)
    ice = lambda M: -1.15 * M**(1.0 / 12.0) * np.exp(Z / 20000.0)
    
    # Compute the fall speed of hydrometeors
    vt = np.ma.where(T >= 0.0, liquid(M), ice(M))
    vt.set_fill_value(fill_value)
    
    vt = {'data': vt,
          'standard_name': 'hydrometeor_fall_velocity',
          'long_name': 'Hydrometeor fall velocity',
          'valid_min': vt.min(),
          'valid_max': vt.max(),
          '_FillValue': vt.fill_value,
          'units': 'meters_per_second',
          'comment': 'Fall speed relations from Caya (2001)'}
    
    # Loop over grids and add the hydrometeor fall velocity field
    for grid in grids:
        grid.fields.update({'hydrometeor_fall_velocity': vt})
    
    return


def _horizontal_divergence(grid, finite_scheme='basic', proc=1,
                           fill_value=None, u_field=None, v_field=None):
    """
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
        
    # Parse field parameters
    if u_field is None:
        u_field = get_field_name('eastward_wind_component')
    if v_field is None:
        v_field = get_field_name('northward_wind_component')
        
    # Get wind data
    u = grid.fields[u_field]['data']
    v = grid.fields[v_field]['data']
    u = np.ma.filled(u, fill_value).astype(np.float64)
    v = np.ma.filled(v, fill_value).astype(np.float64)
    
    # Compute horizontal wind divergence
    div, du, dv = divergence.horizontal_wind(u, v, dx=dx, dy=dy, proc=proc,
                                             finite_scheme=finite_scheme,
                                             fill_value=fill_value)
    
    div = np.ma.masked_equal(div, fill_value)
    
    div = {'data': div,
           'standard_name': 'horizontal_wind_divergence',
           'long_name': 'Divergence of horizontal wind field',
           'valid_min': div.min(),
           'valid_max': div.max(),
           '_FillValue': div.fill_value,
           'units': 'per_second'}
    
    grid.fields.update({'horizontal_divergence': div})
    
    return
        
    
def _check_analysis(grids, conv, sonde, target=None, fall_speed='Caya',
                    finite_scheme='basic', proc=1, standard_density=False,
                    fill_value=None, refl_field=None, vel_field=None,
                    u_field=None, v_field=None, w_field=None, verbose=False):
    """
    Parameters
    ----------
    grids : list
    conv : Grid
    sonde : netCDF4.Dataset
    
    Optional parameters
    -------------------
    
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
    if u_field is None:
        u_field = get_field_name('eastward_wind_component')
    if v_field is None:
        v_field = get_field_name('northward_wind_component')
    if w_field is None:
        w_field = get_field_name('vertical_wind_component')
        
    # Get dimensions
    nz, ny, nx = grids[0].fields[vel_field]['data'].shape
    
    # Get target time. If no target time is provided, use the earliest time
    # out of the list of grids (radars)
    if target is None:
        target = min([datetime_from_grid(grid) for grid in grids])
    
    # Use the ARM interpolated or merged sounding product to get the
    # atmospheric thermodynamic and horizontal wind profiles
    res = _arm_interp_sonde(grids[0], sonde, target=target,
                            standard_density=standard_density,
                            fill_value=fill_value, debug=False,
                            verbose=verbose)
    temp, pres, rho, drho, us, vs = res
    
    # Get wind data
    u = np.ma.filled(conv.fields[u_field]['data'], fill_value)
    v = np.ma.filled(conv.fields[v_field]['data'], fill_value)
    w = np.ma.filled(conv.fields[w_field]['data'], fill_value)
    
    # Initialize metrics dictionary
    metrics = collections.defaultdict(list)
    
    # Compute the RMSD of the radial velocity field for each grid (radar)
    for grid in grids:
        # Get appropriate grid (radar) data
        vr_obs = grid.fields[vel_field]['data']
        vt = grid.fields['hydrometeor_fall_velocity']['data']
        ic = grid.fields['x_component']['data']
        jc = grid.fields['y_component']['data']
        kc = grid.fields['z_component']['data']
        radar_name = grid.metadata['radar_0_instrument_name']
        
        # Compute the projected radial velocity from the wind analysis
        vr = u * ic + v * jc + (w + vt) * kc
        
        # Compute radial velocity RMSD
        rmse_vr = np.ma.sqrt(((vr - vr_obs)**2).mean())
        metrics['RMSD %s' %radar_name] = rmse_vr
        
        if verbose:
            print ('The radial velocity RMSD for radar %s is %.3f m/s'
                   %(radar_name, rmse_vr))
            
    # Compute the normalized anelastic continuity residual
    res = divergence.full_wind(u, v, w, dx=dx, dy=dy, dz=dz, proc=proc,
                               finite_scheme=finite_scheme,
                               fill_value=fill_value)
    div, du, dv, dw = res
    
    for k in xrange(nz):
        # At each height compute anelastic air mass continuity (D) as well as
        # the square of each of its individual terms
        D = div[k,:,:] + w[k,:,:] * drho[k] / rho[k]
        du2 = du[k,:,:]**2
        dv2 = dv[k,:,:]**2
        dw2 = dw[k,:,:]**2
        compress2 = (w[k,:,:] * drho[k] / rho[k])**2

        # Now compute the normalized anelastic continuity residual for the
        # current height
        numer = np.ma.sqrt((D**2).mean())
        denom = np.ma.sqrt((du2 + dv2 + dw2 + compress2).mean())
        normalized_residual = 100.0 * numer / denom
        metrics['normalized continuity residual'].append(normalized_residual)
    
    if verbose:
        min_nd = min(metrics['normalized continuity residual'])
        max_nd = max(metrics['normalized continuity residual'])
        print 'Minimum normalized continuity residual = %.3f%%' %min_nd
        print 'Maximum normalized continuity residual = %.3f%%' %max_nd
              
    # Compute the RMSD of the impermeability condition at the surface to the
    # analysis vertical velocity at the surface
    rmsd_w0 = np.sqrt((w[0,:,:]**2).mean())
    metrics['RMSD impermeability'] = rmsd_w0
    
    if verbose:
        print 'The impermeability condition RMSD is %.3f m/s' %rmsd_w0
    
    return metrics
    

def solve_wind_field(grids, sonde=None, target=None, technique='3d-var',
                     solver='scipy.fmin_cg', first_guess='zero',
                     background='sounding', sounding='ARM interpolated',
                     fall_speed='Caya', finite_scheme='basic',
                     continuity='potvin', smooth='potvin',
                     impermeability='strong', first_pass=True,
                     sub_beam=False, gtol=1.0e-5, ftol=1.0e7, maxiter=100,
                     maxcor=10, disp=False, length_scale=None, use_qc=True,
                     mds=0.0, ncp_min=0.4, rhv_min=0.8, vel_max=40.0,
                     vel_grad_max=10.0, noise_ratio=50.0,
                     window_size=10, min_layer=1500.0, top_offset=500.0,
                     standard_density=False, save_refl=True, proc=1,
                     proj='lcc', datum='NAD83', ellps='GRS80',
                     fill_value=None, refl_field=None, vel_field=None,
                     ncp_field=None, rhv_field=None, u_field=None,
                     v_field=None, w_field=None, debug=False,
                     verbose=False,):
    """
    Parameters
    ----------
    grids : list
        All available radar grids to use in the wind retrieval.
        
    Optional parameters
    -------------------
    sonde : netCDF4.Dataset
        Radiosonde dataset.
    target : datetime
        Target date and time. If target is not provided, the earliest time
        out of all the grids will be used.
    technique : '3d-var'
        Technique used to derive the 3-D wind field.
    solver : 'scipy.fmin_cg', 'scipy.fmin_bfgs'
        Algorithm used to solve the posed problem.
    first_guess : 'zero', 'sounding'
        Define what to use as a first guess field.
    background : 'zero', 'sounding'
        Define what to use as a background field.
    sounding : 'ARM interpolated'
        Define the source for the radiosonde data.
    continuity : 'potvin' or 'iterative'
        Define what method to use as the anelastic air mass continuity
        constraint. The 'iterative' method is a so-called non-simultaneaous
        method, while the 'potvin' method is a simultaneaous method.
    smooth : 'potvin' or 'collis'
    impermeability : 'strong' or 'weak'
        Determines whether surface impermeability is considered as a strong
        or weak constraint.
    wgt_o : float
        Observation weight used at each grid point with valid observations.
        Only applicable when 'technique' is '3d-var'.
    wgt_c : float
        Weight given to the anelastic air mass continuity constraint. Only
        applicable when 'technique' is '3d-var'.
    wgt_s : list
        Weights given to the smoothness constraint. It should be a list of
        five floats, the first only applicable when smooth is 'collis', and
        the other four only applicable when smooth is 'potvin'. When smooth
        is 'potvin', the four floats correspond to S1, S2, S3, and S4 from
        Equation (6) in Potvin et al. (2012). Only applicable when 'technique'
        is '3d-var'.
    wgt_b : list of 3 floats
        Weight given to the background field. Only applicable when
        'technique' is '3d-var'.
    wgt_w0 : float
        Weight given to satisfying the impermeability condition at the
        surface. Only applicable when 'technique' is '3d-var' and
        impermeability is 'weak'.
    length_scale : float
        Only applicable when 'technique' is '3d-var'.
    first_pass : bool
        True to perform a heavily-smoothed first pass which is designed
        to retrieve the large-scale horizontal flow. Only applicable when
        'technique' is '3d-var'.
    mds : float
        Minimum detectable signal in dBZ used to define the minimum
        reflectivity value.
    ncp_min, rhv_min : float
        Minimum values allowed in the normalized coherent power
        and correlation coefficient fields, respectively.
    min_layer : float
        Minimum cloud thickness in meters allowed
    save_refl : bool
        If true, then the reflectivity of the radar network is returned as
        a field in the grid object. The reflectivity of the radar network is
        defined as the maximum reflectivity value observed from all input
        grids (radars).
    gtol : float
    
    ftol: float
    
    maxiter : int
    
    maxcor : int
        Only applicable when solver is 'scipy'.
    disp : bool
        Only applicable when 'solver' is 'scipy'.
    retall : bool
        Only applicable when 'solver' is 'scipy'.
    proc : int
        Number of processors requested.
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.

    Returns
    -------
    conv : Grid
        A grid containing the 3-D Cartesian wind components.

    References
    ----------
    Potvin, C., 2012:
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if u_field is None:
        u_field = get_field_name('eastward_wind_component')
    if v_field is None:
        v_field = get_field_name('northward_wind_component')
    if w_field is None:
        w_field = get_field_name('vertical_wind_component')
        
    if verbose:
        print 'Observations from %i radar(s) will be used' %len(grids)
    
    # Get dimensions of posed problem. We will eventually have to permute the
    # problem from the initial grid space which is in (nz, ny, nx) to a vector
    # space which is 1-D
    #
    # If the total number of grid points in the grid space is N, then we
    # will need a vector of length 3N since we are solving for the 3
    # Cartesian wind components
    nz, ny, nx = grids[0].fields[refl_field]['data'].shape
    N = nz * ny * nx
    if verbose:
        print 'We have to minimize a function of %i variables' %(3 * N)
        
    # Get target time. If no target time is provided, use the earliest time
    # out of the list of grids (radars)
    if target is None:
        target = min([datetime_from_grid(grid) for grid in grids])

    # If no radiosonde data is provided, we will use a standard atmosphere
    # profile to define the density and its derivative with respect to height
    if sonde is None:
        if first_guess == 'sounding':
            raise ValueError('Radiosonde data is required for the initial '
                             '(first) guess field but none is provided')
        if background == 'sounding':
            raise ValueError('Radiosonde data is required for the background '
                             'field but none is provided')
        temp, pres, rho, drho = _standard_atmosphere(grids[0])
    
    # Use the ARM interpolated or merged sounding product to get the
    # atmospheric thermodynamic and horizontal wind profiles
    if sounding == 'ARM interpolated':
        res = _arm_interp_sonde(grids[0], sonde, target=target,
                                standard_density=standard_density,
                                fill_value=fill_value, debug=debug,
                                verbose=verbose)
        temp, pres, rho, drho, us, vs = res
    else:
        raise ValueError('Unsupported radiosonde data source')
              
    # Get the first guess field. Here we will put the variables into their
    # vector space,
    #
    # u0 = (u1, u2, ... , uN)
    # v0 = (v1, v2, ... , vN)
    # w0 = (w1, w2, ... , wN)
    if first_guess == 'zero':
        u0 = np.zeros(N, dtype=np.float64)
        v0 = np.zeros(N, dtype=np.float64)
        w0 = np.zeros(N, dtype=np.float64)
    elif first_guess == 'sounding':
        u0 = np.ravel(np.repeat(us, ny*nx, axis=0).reshape(nz, ny, nx))
        v0 = np.ravel(np.repeat(vs, ny*nx, axis=0).reshape(nz, ny, nx))
        w0 = np.zeros(N, dtype=np.float64)
    else:
        raise ValueError('Unsupported initial (first) guess field')
        
    # Get the background field. Here we will put the variables into their
    # grid space,
    #
    # ub = ub(z, y, x)
    # vb = vb(z, y, x)
    # wb = wb(z, y, x)
    if background == 'zero':
        ub = np.zeros((nz,ny,nx), dtype=np.float64)
        vb = np.zeros((nz,ny,nx), dtype=np.float64)
        wb = np.zeros((nz,ny,nx), dtype=np.float64)
    elif background == 'sounding':
        ub = np.repeat(us, ny*nx, axis=0).reshape(nz, ny, nx)
        vb = np.repeat(vs, ny*nx, axis=0).reshape(nz, ny, nx)
        wb = np.zeros((nz,ny,nx), dtype=np.float64)
    else:
        raise ValueError('Unsupported background field')
    
    # Quality control procedures. This attempts to remove noise and
    # non-meteorological returns
    if use_qc:
        _radar_qc(grids, mds=mds, vel_max=vel_max, vel_grad_max=vel_grad_max,
                  ncp_min=ncp_min, rhv_min=rhv_min, window_size=window_size,
                  noise_ratio=noise_ratio, refl_field=refl_field,
                  vel_field=vel_field, ncp_field=ncp_field,
                  rhv_field=rhv_field)
        
    # Add the Cartesian components field to all the grids
    _radar_components(grids, proj=proj, datum=datum, ellps=ellps)
    
    # Get fall speed of hydrometeors
    if fall_speed == 'Caya':
        _fall_speed_caya(grids, temp, fill_value=fill_value,
                         refl_field=refl_field)
    else:
        raise ValueError('Unsupported fall speed relation')
        
    # Get radar coverage from all grids
    cover = _radar_coverage(grids, refl_field=refl_field, vel_field=vel_field)
    
    # Get echo base and top heights. This is only applicable if we are using
    # an iterative technique since this technique has to explicitly integrate
    # the anelastic continuity equation and therefore we require boundary
    # conditions
    if continuity == 'iterative':
        base, top = _echo_bounds(grids, mds=mds, min_layer=min_layer,
                                 top_offset=top_offset, fill_value=fill_value,
                                 refl_field=refl_field)
        base['data'] = np.ma.filled(base['data'], fill_value)
        top['data'] = np.ma.filled(top['data'], fill_value)

        # Get column types
        column = _column_types(cover, base, top, fill_value=fill_value)
    
    # This is an important step. We make sure that every field for every grid
    # (radar) is a NumPy array and not a masked array. First make a copy of
    # the original grids so as not to disturb their data array types
    grids_nd = deepcopy(grids)
    
    for grid in grids_nd:
        for field, field_dic in grid.fields.iteritems():
            field_dic['data'] = np.ma.filled(field_dic['data'], fill_value)
            grid.fields[field]['data'] = field_dic['data'] 
    
    # Now the very important step of concatenating the arrays of the
    # initial guess field so that we now have,
    #
    # xo = (u1, ... , uN, v1, ... , vN, w1, ... , wN)
    #
    # which is the space in which we solve the wind retrieval problem
    x0 = np.concatenate((u0, v0, w0), axis=0)
    
    # The first block is for when a 3-D variational algorithm with conjugate
    # gradient minimization from the SciPy optimization toolkit is used
    if technique == '3d-var' and 'scipy' in solver:
        
        # Multiply each weighting coefficient (tuning parameter) by the length
        # scale if necessary. The length scale is designed to make the
        # dimensionality of each cost uniform, as well as bring each cost
        # within 1-3 orders of magnitude of each other
        if length_scale is not None:
        
            if continuity == 'potvin':
                wgt_c = wgt_c * length_scale**2
            if smooth == 'potvin':
                wgt_s1 = wgt_s1 * length_scale**4
                wgt_s2 = wgt_s2 * length_scale**4
                wgt_s3 = wgt_s3 * length_scale**4
                wgt_s4 = wgt_s4 * length_scale**4
                
        # Add observation weight field to all grids
        _observation_weight(grids_nd, wgt_o=wgt_o, fill_value=fill_value,
                            refl_field=refl_field, vel_field=vel_field)
        
        # SciPy nonlinear conjugate gradient method
        if solver == 'scipy.fmin_cg':
            method = 'CG'
            opts = {'maxiter': maxiter, 'gtol': gtol, 'disp': disp}
        
        # SciPy Broyden-Fletcher-Goldfarb-Shanno method. Using this algorithm
        # requires storage of a 3N x 3N matrix (the Hessian matrix), and in
        # most wind retrieval applications will throw a MemoryError
        elif solver == 'scipy.fmin_bfgs':
            method = 'BFGS'
            opts = {'maxiter': maxiter, 'gtol': gtol, 'disp': disp}
            
        # SciPy limited memory Broyden-Fletcher-Goldfarb-Shanno method
        elif solver == 'scipy.fmin_l_bfgs_b':
            method = 'L-BFGS-B'
            opts = {'maxiter': maxiter, 'ftol': ftol, 'gtol': gtol,
                    'maxcor': maxcor, 'disp': disp}
            
        # SciPy truncated Newton method
        elif solver == 'scipy.fmin_tnc':
            method = 'TNC'
            opts = {'maxiter': maxiter, 'gtol': gtol, 'ftol': ftol,
                    'disp': disp}
            
        else:
            raise ValueError('Unsupported SciPy solver')
        
        # This is a very important step. Group the required arguments for the
        # cost function and gradient together. The order at which the
        # arguments are entered is very important, and must be adhered to by
        # the functions that use it as arguments
        args = (nx, ny, nz, N, grids_nd, ub, vb, wb, rho, drho, base['data'],
                top['data'], column['data'], wgt_o, wgt_c, wgt_s, wgt_b,
                wgt_w0, continuity_cost, smooth_cost, dx, dy, dz, sub_beam,
                finite_scheme, fill_value, proc, vel_field, debug, verbose)
        
        # Get the appropriate cost function and the function that computes
        # the gradient of the cost function
        f = cga._cost_wind
        jac = cga._grad_wind
        
        # Check if the user wants to perform a first pass which retrieves a
        # heavily-smoothed horizontal wind field and a w field of 0 everywhere
        if first_pass:
            if verbose:
                print 'Performing heavily-smoothed first pass'
            # Here we set the appropriate weights for the individual cost
            # functions. Recall that the first pass is meant to retrieve
            # the large-scale horizontal wind flow. This means that we
            # want to ignore the continuity cost, as well as set the
            # smoothness weight "very high" to produce a heavily-smoothed
            # horizontal wind field
            continuity_cost0 = None
            wgt_s0 = [10.0 * wgt for wgt in wgt_s]
            
            # This is an important step. We have to change the arguments to
            # account for different costs and weighting coefficients for the
            # first pass
            args0 = (nx, ny, nz, N, grids_nd, ub, vb, wb, rho, drho,
                     base['data'], top['data'], column['data'], wgt_o,
                     wgt_c, wgt_s0, wgt_b, wgt_w0, continuity_cost0,
                     smooth_cost, dx, dy, dz, sub_beam, finite_scheme,
                     fill_value, proc, vel_field, debug, verbose)
            
            # Call the SciPy solver
            res = minimize(f, x0, args=args0, method=method, jac=jac,
                           hess=None, hessp=None, bounds=None,
                           constraints=None, options=opts)
            
            # Unpack the results from the first pass
            x0 = res.x
            x0 = x0.astype(np.float64)
            
            # Set the vertical velocity field to 0 m/s everywhere after the
            # first pass
            x0[2*N:3*N] = 0.0
        
        # Call the SciPy solver to perform the retrieval of the full 3-D
        # wind field     
        res = minimize(f, x0, args=args, method=method, jac=jac, hess=None,
                       hessp=None, bounds=None, constraints=None,
                       options=opts)
        
        # Unpack the full 3-D wind field results
        xopt = res.x
        xopt = xopt.astype(np.float64)
       
    else:
        raise ValueError('Unsupported technique and solver combination')
    
    # Now the last step: prepare the necessary data for output. First get the
    # wind retrieval metadata from configuration file
    u = get_metadata(u_field)
    v = get_metadata(v_field)
    w = get_metadata(w_field)
    
    u['_FillValue'] = fill_value
    v['_FillValue'] = fill_value
    w['_FillValue'] = fill_value
    
    # This is an important step. Get the control variables from the analysis
    # vector. This requires us to keep track of how the analysis vector is
    # ordered, since the way in which we slice it requires this knowledge. We
    # assume the analysis vector is of the form,
    #
    # x = x(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
    #
    # so u is packed first, then v, and finally w. At the same time we will
    # permute the control variables back to the grid space (3-D)
    u['data'] = np.reshape(xopt[0:N], (nz,ny,nx))
    v['data'] = np.reshape(xopt[N:2*N], (nz,ny,nx))
    w['data'] = np.reshape(xopt[2*N:3*N], (nz,ny,nx))
    
    u['data'] = np.ma.masked_equal(u['data'], fill_value)
    v['data'] = np.ma.masked_equal(v['data'], fill_value)
    w['data'] = np.ma.masked_equal(w['data'], fill_value)
    
    # Get hydrometeor fall velocity since we will add this to the fields
    # of the output grid object
    vt = grids[0].fields['hydrometeor_fall_velocity']
    
    # Define the fields
    fields = {
            u_field: u,
            v_field: v,
            w_field: w,
            'radar_coverage': cover,
            'hydrometeor_fall_velocity': vt
            }
    
    # Calculate the maximum reflectivity from all grids, which we will save
    # as a field in the output grid object
    if save_refl:
        # Get metadata from configuration file
        ze = get_metadata(refl_field)
        
        # Get data
        ze['data'] = np.ma.max([grid.fields[refl_field]['data'] for
                                grid in grids], axis=0)
        ze['data'].set_fill_value(fill_value)
        
        # Define more metadata
        ze['long_name'] = 'Reflectivity of radar network'
        ze['_FillValue'] = fill_value
        ze['comment'] = ('Reflectivity values are maximum values, '
                         'not mean values')
        
        fields[refl_field] = ze
    
    # Define the axes
    # We will populate the axes with that of the first grid and then
    # update this accordingly
    axes = deepcopy(grids[0].axes)
    
    # Remove only the time axes since we will populate these ourself
    [axes.pop(key) for key in axes.keys() if 'time' in key]
    
    # Find the grid with the earliest volume start time and the grid with
    # the latest volume end time. This will become the start and end times
    volume_starts = [num2date(grid.axes['time_start']['data'][0],
                              grid.axes['time_start']['units'])
                     for grid in grids]
    volume_ends = [num2date(grid.axes['time_end']['data'][0],
                            grid.axes['time_end']['units'])
                   for grid in grids]
    grid_start = grids[np.argmin(volume_starts)]
    grid_end = grids[np.argmax(volume_ends)]
    
    # Populate the start and end time axes
    axes['time_start'] = grid_start.axes['time_start']
    axes['time_end'] = grid_end.axes['time_end']
    
    # Populate the time axis 
    # This time will correspond to the time half-way between
    # the start and end times
    td = max(volume_ends) - min(volume_starts)
    seconds_since_start = (td.seconds + td.days * 24 * 3600) / 2 
    axes['time'] = {
            'data': np.array(seconds_since_start, np.float64),
            'long_name': 'Time in seconds since volume start',
            'calendar': 'gregorian',
            'units': axes['time_start']['units']
            }
    
    # ARM time variables
    dt = num2date(axes['time']['data'], axes['time']['units'])
    td = dt - datetime.utcfromtimestamp(0)
    td = td.seconds + td.days * 24 * 3600
    
    axes['base_time'] = {
            'data': np.array(td, np.int32),
            'long_name': 'Base time in Epoch',
            'string': dt.strftime('%d-%b-%Y,%H:%M:%S GMT'),
            'units': 'seconds since 1970-1-1 0:00:00 0:00',
            'ancillary_variables': 'time_offset'
            }
    
    axes['time_offset'] = {
            'data': np.array(axes['time']['data'], np.float64),
            'long_name': 'Time offset from base_time',
            'units': axes['time']['units'].replace('T',' ').replace('Z',''),
            'ancillary_variables': 'base_time',
            'calendar': 'gregorian'
            }
    
    # Define the metadata
    datastreams_description = ('A string consisting of the datastream(s), '
                               'datastream version(s), and datastream '
                               'date (range).')
    
    metadata = {
        'title': 'Convective Vertical Velocities',
        'dod_version': '',
        'process_version': '',
        'command_line': '',
        'site_id': '',
        'facility_id': '',
        'source': '',
        'Conventions': '',
        'references': '',
        'input_datastreams_description': datastreams_description,
        'input_datastreams_num': '',
        'input_datastreams': '',
        'state': '',
        'history': ''
        }
    
    return Grid(fields, axes, metadata)
