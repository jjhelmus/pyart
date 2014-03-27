"""
pyart.retrieve.winds
====================

"""

import time
import numpy as np

from warnings import warn
from copy import deepcopy
from scipy.optimize import minimize
from mpl_toolkits.basemap import pyproj

from ..io import Grid
from ..util import datetime_utils
from ..config import get_fillvalue, get_field_name, get_metadata
from ..retrieve import cga, divergence, gradient
                   

def _radar_coverage(grids, refl_field=None, vel_field=None):
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
        
        # Use reflectivity and Doppler velocity masks to find where
        # each radar has coverage, i.e. where each radar has valid
        # observations
        ze = grid.fields[refl_field]['data']
        vr = grid.fields[vel_field]['data']
        has_obs = np.logical_or(~ze.mask, ~vr.mask)
        cover = np.where(has_obs, cover + 1, cover)
        
    # Return dictionary of results
    return {'data': cover.astype(np.int32),
            'standard_name': 'radar_coverage',
            'long_name': 'Radar coverage flags',
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
    
    # Get axes
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
        proj = pyproj.Proj(proj=proj, datum=datum, ellps=ellps, lat_0=lat_0,
                           lon_0=lon_0, x_0=0.0, y_0=0.0)
        
        # Get the (x, y) location of the radar in the analysis domain from
        # the projection
        x_r, y_r = proj(lon_r, lat_r)
        
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
        grid.add_field('x_component', ic)
        grid.add_field('y_component', jc)
        grid.add_field('z_component', kc)
        
    return


def _echo_bounds(network, mds=0.0, min_layer=1500.0, top_offset=500.0,
                 fill_value=None, proc=1, refl_field=None):
    """
    Determine the echo base and top heights from the large-scale
    coverage of the radar network.
    
    Parameters
    ----------
    network : Grid object
        This grid should represent the large-scale coverage within
        the analysis domain. Note that this is only applicable to
        fields like reflectivity.
    
    Optional Parameters
    -------------------
    mds : float
        Minimum detectable signal in dBZ used to define the noise cut-off.
    min_layer : float
        Minimum cloud layer allowed in meters.
    top_offset : float
        Value in meters added to the echo top height.
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.
    proc : int
        Number of processors requested.
    refl_field : str
        Name of reflectivity field which will be used to estimate the fall
        speed. A value of None will use the default field name as defined in
        the Py-ART configuration file.
    
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
    
    # Get axes
    z = network.axes['z_disp']['data']
    
    # Get reflectivity data
    ze = deepcopy(network.fields[refl_field]['data'])
    ze = np.ma.filled(ze, fill_value).astype(np.float64)
        
    # Estimate echo base and top heights. Add offset to echo top heights
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

            
def _observation_weight(grids, wgt_o=1.0, refl_field=None, vel_field=None):
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
    refl_field, vel_field : str
        Name of fields which will be used to determine grid points with
        valid observations. A value of None will use the default field name
        as defined in the Py-ART configuration file.
    
    Returns
    -------
    grids : list
        List of grids with updated observation weight fields.
    """
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
        
    # Initialize observation weight array
    lam_o = np.zeros_like(grids[0].fields[vel_field]['data'])
        
    # Loop over all grids
    for grid in grids:
        
        # Get reflectivity and Doppler velocity data. Use their masks to
        # compute the observation weight for each grid
        ze = grid.fields[refl_field]['data']
        vr = grid.fields[vel_field]['data']
        
        is_good = np.logical_or(~ze.mask, ~vr.mask)
        lam_o[is_good] = wgt_o
        
        # Create dictionary of results
        lam_o = {'data': lam_o.astype(np.float64),
                 'standard_name': 'observation_weight',
                 'long_name': 'Radar observation weight',
                 'valid_min': 0.0,
                 'valid_max': wgt_o}
                 
        # Add new field to grid object
        grid.add_field('observation_weight', lam_o)

    return
    
            
def _radar_qc(grids, mds=0.0, vel_max=55.0, ncp_min=0.3, rhv_min=0.7, 
              window_size=5, noise_ratio=85.0, fill_value=None, 
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
    
    ncp_min, rhv_min : float
    
    window_size : int
    
    noise_ratio : float
    
    fill_value : float
    
    Returns
    -------
    
    
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
    
    # Loop over all grids
    for grid in grids:
        
        # Get data
        ze = grid.fields[refl_field]['data']
        vr = grid.fields[vel_field]['data']
        ncp = grid.fields[ncp_field]['data']
        rhv = grid.fields[rhv_field]['data']
        
        # Create appropriate masks
        is_noise = ze < mds
        is_bad_vel = np.abs(vr) > vel_max
        is_non_meteo = np.logical_or(ncp < ncp_min, rhv < rhv_min)
        
        # Update masks
        ze.mask = np.logical_or(is_noise, is_non_meteo)
        vr.mask = np.logical_or(is_bad_vel, is_non_meteo)
        
        # Despeckle reflectivity and Doppler velocity fields
        grid.despeckle_field(refl_field, window_size=window_size,
                             noise_ratio=noise_ratio,
                             fill_value=fill_value)
        grid.despeckle_field(vel_field, window_size=window_size,
                             noise_ratio=noise_ratio,
                             fill_value=fill_value)
        
    return

    
def _column_types(cover, base, fill_value=None):
    """
    Parameters
    ----------
    
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


def _arm_interp_sonde(grid, sonde, target, fill_value=None,
                      rho0=1.2, H=10000.0, standard_density=False,
                      finite_scheme='basic', debug=False, verbose=False):
    """
    Parameters
    ----------
    
    Optional parameters
    -------------------
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.
    finite_scheme : 'basic' or 'high-order'
        Finite difference scheme to use when calculating density gradient.
        Only applicable if 'standard' is False.
    standard_density : bool
        If True, the returned density profile is from a standard atmosphere.
        False uses the sounding data.
    rho0, H : float
        Reference density in kg m^-3 and scale height in meters. Only
        applicable if 'standard' is True.
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
    z_grid = grid.axes['z_disp']['data'] + grid.axes['alt']['data'] # (m)
    z_sonde = 1000.0 * sonde.variables['height'][:] # (m)
    
    # Get closest time index in sounding to the target time
    dt_sonde = datetime_utils.datetimes_from_dataset(sonde)
    t = np.abs(dt_sonde - target).argmin()
    
    if verbose:
        print 'Closest merged sounding time to target is %s' %dt_sonde[t]
    
    # Get data from sounding
    T = sonde.variables['temp'][t,:] # (C)
    P = sonde.variables['bar_pres'][t,:] # (kPa)
    u = sonde.variables['u_wind'][t,:] # (m/s)
    v = sonde.variables['v_wind'][t,:] # (m/s)
    
    # Now compute the density of dry air and its first derivative with
    # respect to height
    if standard_density:
        rho = rho0 * np.exp(-z_sonde / H)
        drho = -(rho0 / H) * np.exp(-z_sonde / H)
        
    else:
        R = 287.058 # (J K^-1 mol^1)
        rho = (P * 1000.0) / (R * (T + 273.15)) # (kg m^-3)
        drho = gradient.density1d(rho, z_sonde, fill_value=fill_value,
                                  finite_scheme=finite_scheme)
    
    # Interpolate sounding data to vertical dimension of grid
    T = np.interp(z_grid, z_sonde, T, left=None, right=None)
    P = np.interp(z_grid, z_sonde, P, left=None, right=None)
    rho = np.interp(z_grid, z_sonde, rho, left=None, right=None)
    drho = np.interp(z_grid, z_sonde, drho, left=None, right=None)
    u = np.interp(z_grid, z_sonde, u, left=None, right=None)
    v = np.interp(z_grid, z_sonde, v, left=None, right=None)
    
    if debug:
        print 'Minimum air temperature is %.2f C' %T.min()
        print 'Maximum air temperature is %.2f C' %T.max()
        print 'Minimum air pressure is %.2f kPa' %P.min()
        print 'Maximum air pressure is %.2f kPa' %P.max()
        print 'Minimum air density is %.3f kg/m^3' %rho.min()
        print 'Maximum air density is %.3f kg/m^3' %rho.max()
        print 'Minimum eastward wind component is %.2f m/s' %u.min()
        print 'Maximum eastward wind component is %.2f m/s' %u.max()
        print 'Minimum northward wind component is %.2f m/s' %v.min()
        print 'Maximum northward wind component is %.2f m/s' %v.max()
    
    return T, P, rho, drho, u, v

        
def _fall_speed_caya(grid, temp, fill_value=None, refl_field=None):
    """
    Parameters
    ----------
    grid : Grid
        This grid should represent the large-scale coverage within
        the analysis domain. Note that this is only applicable to
        fields like reflectivity.
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
    
    Returns
    -------
    vt : dict
        Hydrometeor fall speed data dictionary.
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
    
    # Parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
    
    # Get dimensions
    nz, ny, nx = grid.fields[refl_field]['data'].shape
    
    # Get height axis and determine its mesh
    z = grid.axes['z_disp']['data']
    Z = np.repeat(z, ny*nx, axis=0).reshape(nz, ny, nx)
    
    # Determine the temperature mesh
    T = np.repeat(temp, ny*nx, axis=0).reshape(nz, ny, nx)
    
    # Get reflectivity data and compute precipitation concentration
    ze = grid.fields[refl_field]['data']
    M = np.exp((ze - 43.1) / 7.6)
    
    # Define liquid and ice relations
    liquid = lambda M: -5.94 * M**(1.0 / 8.0) * np.exp(Z / 20000.0) 
    ice = lambda M: -1.15 * M**(1.0 / 12.0) * np.exp(Z / 20000.0)
    
    # Compute the fall speed of hydrometeors
    vt = np.ma.where(T >= 0.0, liquid(M), ice(M))
        
    # Return dictionary of results
    return {'data': vt,
            'standard_name': 'hydrometeor_fall_velocity',
            'long_name': 'Hydrometeor fall velocity',
            'valid_min': vt.min(),
            'valid_max': vt.max(),
            '_FillValue': vt.fill_value,
            'units': 'meters_per_second',
            'comment': 'Fall speed relations from Caya (2001)'}


def _hor_divergence(grid, dx=500.0, dy=500.0, finite_scheme='basic',
                    fill_value=None, u_field=None, v_field=None):
    """
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
        
    # Parse field parameters
    if u_field is None:
        u_field = get_field_name('u_wind')
    if v_field is None:
        v_field = get_field_name('v_wind')
        
    # Get axes
    z = grid.axes['z_disp']['data']
        
    # Get wind data
    u = grid.fields[u_field]['data']
    v = grid.fields[v_field]['data']
    u = np.ma.filled(u, fill_value).astype(np.float64)
    v = np.ma.filled(v, fill_value).astype(np.float64)
    
    # Compute horizontal wind divergence
    div, du, dv = divergence.horiz_wind(u, v, dx=dx, dy=dy, proc=proc,
                                        finite_scheme=finite_scheme,
                                        fill_value=fill_value)
    
    div = np.ma.masked_equal(div, fill_value)
    
    return {'data': div,
           'standard_name': 'horizontal_wind_divergence',
           'long_name': 'Divergence of horizontal wind field',
           'valid_min': div.min(),
           'valid_max': div.max(),
           '_FillValue': div.fill_value,
           'units': 'per_second'}
   
    
def _cost_magnitudes(grids, dx=500.0, dy=500.0, dz=500.0, fill_value=None,
                     continuity_cost='original', smooth_cost='original',
                     finite_scheme='basic', verbose=False, vel_field=None):
    """
    Parameters
    ----------
    
    Optional parameters
    -------------------
    
    Returns
    -------
    
    References
    ----------
    Shapiro, A., C. K. Potvin, and J. Gao, 2009: Use of a Vertical Vorticity
    Equation in Variational Dual-Doppler Wind Analysis. J. Atmos. Oceanic
    Technol., 26, 2089-2106
    """
    
    # Get fill value
    
    if fill_value is None: fill_value = get_fillvalue()
    
    # Parse the field parameters
    
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
        
    
    # Get the size of the analysis domain
    
    nz, ny, nx = grids[0].fields[vel_field]['data'].shape
    
    
    # Estimate the magnitude of the observation cost function
    #
    # We do this by using the actual radial velocity observations for each
    # radar (grid). Recall that the observation cost is given by,
    #
    # Jo = 0.5 * sum( wgt_o * [ vr - vr_obs ]**2 )
    #
    # so we can get an idea of magnitude of the observation cost by
    # neglecting the analysis radial velocity and summing the squares of the
    # observations for each radar
    
    N = 0.0
    Jo = 0.0
    
    for grid in grids:
        
        vr = deepcopy(grid.fields[vel_field]['data'])
        
        N = N + np.logical_not(vr.mask).sum()
        
        Jo = Jo + np.sum(vr**2, dtype='float64')
        
    
    if verbose:
        
        print 'Total number of analysis grid points is %i' %(nz * ny * nx)
        print 'Total number of radar observations is   %i' %N
              
        print 'Expect the observation cost to be on the order %.0e' %Jo
        
        
def _check_analysis(grids, network, u, v, w, T, u0, v0, w0, dx=500.0,
                    dy=500.0, dz=500.0, finite_scheme='basic',
                    fall_speed='caya', fill_value=None, proc=1,
                    verbose=False, refl_field=None, vel_field=None):
    """
    Parameters
    ----------
    
    Optional parameters
    -------------------
    
    Returns
    -------
    """
    
    # Get fill value
    
    if fill_value is None: fill_value = get_fillvalue()
    
    # Parse the field parameters
    
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
        
    if vel_field is None:
        vel_field = get_field_name('corrected_velocity')
        
    
    # We copy the grids since they are in fact mutable objects and any
    # changes we make to them will be reflected in the outer scope
    # (assuming they are passed by reference)
    
    grids = deepcopy(grids)
    network = deepcopy(network)
        
    # Get grid dimensions
    
    nz, ny, nx = grids[0].fields[vel_field]['data'].shape
    
    
    # Calculate the fall speed of hydrometeors using the reflectivity
    # field from the network
    
    if fall_speed == 'caya':
        
        vt = _fall_speed_caya(network, T, fill_value=fill_value,
                              refl_field=refl_field)
    
    
    # Get coverage
    
    cover = _radar_coverage(grids, refl_field=refl_field, vel_field=vel_field)
    
    # Add the Cartesian components field to all the grids
    
    grids = _radar_components(grids, proc=1)
    
    # Compute the RMSE of the radial velocity field for each grid (radar)
    
    for grid in grids:
        
        vr_obs = grid.fields[vel_field]['data']
        
        ic = grid.fields['x_component']['data']
        jc = grid.fields['y_component']['data']
        kc = grid.fields['z_component']['data']
        
        vr = u['data'] * ic + v['data'] * jc + (w['data'] + vt['data']) * kc
        
        rmse_vr = np.sqrt(np.ma.mean((vr - vr_obs)**2))
        
        if verbose:
            
            print 'The radial velocity RMSE for radar %s is %.3f m/s' \
                    %(grid.metadata['radar_0_instrument_name'], rmse_vr)
            
        if rmse_vr > 2.0:
            
            warn('The radial velocity RMSE for radar %s is greater than 2 m/s'
                 %grid.metadata['radar_0_instrument_name'])
            
    # Compute the normalized divergence profile
    #
    # This is defined on each analysis height level as the ratio of RMS
    # velocity divergence to the root of the mean of the sum of the squares
    # of the three terms comprising the velocity divergence
    
    norm_div = np.zeros(nz, dtype='float64')
    
    div, du, dv, dw = divergence.full_wind(u['data'], v['data'], w['data'],
                                           dx=dx, dy=dy, dz=dz, proc=proc,
                                           finite_scheme=finite_scheme,
                                           fill_value=fill_value)
    
    div = np.ma.masked_where(cover['data'] < 1, div)
    du = np.ma.masked_where(cover['data'] < 1, du)
    dv = np.ma.masked_where(cover['data'] < 1, dv)
    dw = np.ma.masked_where(cover['data'] < 1, dw)
    
    for k in xrange(nz):
        
        num = np.sqrt(np.ma.mean(div[k,:,:]**2))
        den = np.sqrt(np.ma.mean(du[k,:,:]**2 + dv[k,:,:]**2 + dw[k,:,:]**2))
        
        norm_div[k] = 100.0 * num / den
        
    if verbose:
        
        print 'Minimum normalized divergence = %.3f%%' %norm_div.min()
        print 'Maximum normalized divergence = %.3f%%' %norm_div.max()
        
    if norm_div.max() > 5.0:
        
        warn('The maximum normalized divergence is > 5%')
        
    
    # Compute the RMSD between the initial guess field to the analysis field
    
    if u0.ndim == v0.ndim == w0.ndim == 1:
        
        u0 = np.repeat(u0, ny * nx, axis=0).reshape(nz, ny, nx)
        v0 = np.repeat(v0, ny * nx, axis=0).reshape(nz, ny, nx)
        w0 = np.repeat(w0, ny * nx, axis=0).reshape(nz, ny, nx)
        
    
    rmsd_u = np.sqrt(np.ma.mean((u['data'] - u0)**2))
    rmsd_v = np.sqrt(np.ma.mean((v['data'] - v0)**2))
    rmsd_w = np.sqrt(np.ma.mean((w['data'] - w0)**2))
    
    if verbose:
        
        print 'The RMSD between initial guess for the eastward, northward, and ' + \
              'vertical wind components are %.3f, %.3f, and %.3f m/s' \
              %(rmsd_u, rmsd_v, rmsd_w)
              
    # Compute the RMSE of the impermeability condition at the surface to the
    # analysis vertical velocity at the surface
    
    rmsd_w0 = np.sqrt(np.ma.mean(w['data'][0,:,:]**2))
    
    if verbose:
        
        print 'The RMSD between the impermeability condition at the ' + \
              'surface and the analyzed vertical velocity at the ' + \
              'surface is %.3f m/s' %rmsd_w0
              
    
    return norm_div
    

def solve_wind_field(grids, network, sonde, target, dx=500.0, dy=500.0,
                     dz=500.0, mds=0.0, vel_max=50.0, ncp_min=0.3,
                     rhv_min=0.7, technique='3d-var', solver='scipy.fmin_cg',
                     gtol=1.0e-5, ftol=1.0e7, maxcor=10, maxiter=200,
                     disp=True, retall=False, finite_scheme='basic',
                     first_guess='zero', background='sounding', wgt_o=3.0,
                     wgt_c=3.0, wgt_s=[0.05,1.0,1.0,1.0,0.1],
                     wgt_b=[0.01,0.01,0.0], wgt_w0=0.0, length_scale=None,
                     smooth_cost='original', continuity_cost='integrate',
                     fall_speed='caya', min_layer=1500.0, top_offset=500.0,
                     first_pass=True, use_qc=True, sub_beam=True,
                     standard_density=False, window_size=5, window_max=0.85,
                     fill_value=None, proc=1, verbose=False, debug=False,
                     refl_field=None, vel_field=None, ncp_field=None,
                     rhv_field=None):
    """
    Parameters
    ----------
    grids : list
        All available radar grids to use in the wind retrieval.
    network : Grid
        This grid should represent the large-scale coverage within
        the analysis domain. Note that this is only applicable to
        fields like reflectivity. It will be used to estimate the echo
        base and top heights.
    sonde : netCDF4.Dataset
        Sounding dataset.
    target : datetime
        
    Optional parameters
    -------------------
    dx, dy, dz : float
        Grid resolution in x-, y-, and z-dimension, respectively.
    wgt_o : float
        Observation weight used at each grid point with valid observations.
        Only applicable when 'technique' is '3d-var'.
    wgt_c : float
        Weight given to the anelastic air mass continuity constraint. Only
        applicable when 'technique' is '3d-var'.
    wgt_s : list of 5 floats
        Weight given to the smoothness constraint.
    wgt_b : list of 3 floats
        Weight given to the background field.
    wgt_w0 : float
        Weight given to satisfying the impermeability condition at the surface.
    mds : float
        Minimum detectable signal in dBZ used to define the minimum
        reflectivity value.
    ncp_min, rhv_min : float
        Minimum values allowed in the normalized coherent power
        and correlation coefficient fields, respectively
    technique : '3d-var'
        Technique used to derive wind field
    solver : 'scipy.fmin_cg', 'scipy.fmin_bfgs'
        Algorithm used to solve the posed problem
    first_guess : 'zero', 'sounding'
        Define what to use as a first guess field
    background : sounding'
        Define what to use as a background field
    length_scale : float
    
    first_pass : bool
        True to perform a heavily-smoothed first pass which is designed
        to retrieve the large-scale horizontal flow
    min_layer : float
        Minimum cloud thickness in meters allowed
    gtol : float
    
    ftol: float
    
    maxiter : int
    
    maxcor : int
    
    disp : bool
        Only applicable when 'solver' is 'scipy'
    
    retall : bool
        Only applicable when 'solver' is 'scipy'
    fill_value : float
        Missing value used to signify bad data points. A value of None
        will use the default fill value as defined in the Py-ART
        configuration file.
    proc : int
        Number of processors available.
    
    Returns
    -------
    conv : Grid
        A grid containing the 3-D Cartesian wind components.
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
        
    if verbose:
        print 'Observations from %i radar(s) will be used' %len(grids)
    
    # We copy the grids since they are mutable objects when they are passed
    # by reference and therefore any changes we make to them will be reflected
    # in the outer scope
    grids = deepcopy(grids)
    network = deepcopy(network)
    
    # Get dimensions of problem. We will eventually have to permute the
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
        
    # Multiply each weighting coefficient (tuning parameter) by the length
    # scale, if necessary. The length scale is designed to make the
    # dimensionality of each cost uniform, as well as bring each cost within
    # 1-3 orders of magnitude of each other
    if length_scale is not None:
        
        if continuity_cost == 'potvin':
            wgt_c = wgt_c * length_scale**2
            
        if smooth_cost == 'potvin':
            wgt_s = [wgt * length_scale**4 for wgt in wgt_s]
    
    # Use the ARM interpolated or merged sounding product to get the
    # atmospheric thermodynamic and horizontal wind profiles
    T, P, rho, drho, us, vs = _arm_interp_sonde(grids[0], sonde, target,
                                    standard_density=standard_density,
                                    fill_value=fill_value, debug=debug,
                                    verbose=verbose)
              
    # Get the first guess field. Here we will put the variables into their
    # vector space,
    #
    # u0 = (u1,u2,...,uN)
    # v0 = (v1,v2,...,vN)
    # w0 = (w1,w2,...,wN)
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
    # ub = ub(z,y,x)
    # vb = vb(z,y,x)
    # wb = wb(z,y,x)
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
    
    # Quality control procedures. This attempts to remove noise from gridded
    # reflectivity and non-meteorological returns using polarization data
    if use_qc:
        _radar_qc(grids, mds=mds, vel_max=vel_max, ncp_min=ncp_min,
                  rhv_min=rhv_min, refl_field=refl_field, vel_field=vel_field,
                  ncp_field=ncp_field, rhv_field=rhv_field)
 
    # Add observation weight field to all grids
    _observation_weight(grids, wgt_o=wgt_o, refl_field=refl_field,
                        vel_field=vel_field)
        
    # Add the Cartesian components field to all the grids
    _radar_components(grids)
    
    # Get fall speed of hydrometeors
    if fall_speed == 'caya':
        vt = _fall_speed_caya(network, T, fill_value=fill_value,
                              proc=proc, refl_field=refl_field)
        
        vt['data'] = np.ma.filled(vt['data'], fill_value)
        
    else:
        raise ValueError('Unsupported fall speed relation')
        
    # Get radar coverage from all grids
    cover = _radar_coverage(grids, refl_field=refl_field, vel_field=vel_field)
    
        
    # Get echo base and top heights
    base, top = _echo_bounds(network, mds=mds, min_layer=min_layer,
                        top_offset=top_offset, fill_value=fill_value,
                        refl_field=refl_field)
    
    base['data'] = np.ma.filled(base['data'], fill_value)
    top['data'] = np.ma.filled(top['data'], fill_value)
    
        
    # Get column types
    column = _column_types(cover, base, fill_value=fill_value)
    
    # This is an important step. We turn the velocity fields for every grid
    # from NumPy masked arrays into NumPy arrays
    for grid in grids:
        vr = grid.fields[vel_field]['data']
        grid.fields[vel_field]['data'] = np.ma.filled(vr, fill_value) 
    
    # Now the very important step of concatenating the arrays of the
    # initial guess field so that we now have,
    #
    # xo = (u1,...,uN,v1,...,vN,w1,...,wN)
    #
    # which is the space in which we solve the wind retrieval problem
    x0 = np.concatenate((u0,v0,w0), axis=0)
    
    # The first block is for when a 3-D variational algorithm with conjugate
    # gradient minimization from the SciPy optimization toolkit is used
    if technique == '3d-var' and 'scipy' in solver:
        
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
        args = (nx, ny, nz, N, grids, ub, vb, wb, vt['data'], rho, drho,
                base['data'], top['data'], column['data'], wgt_o, wgt_c,
                wgt_s, wgt_b, wgt_w0, continuity_cost, smooth_cost, dx, dy,
                dz, sub_beam, finite_scheme, fill_value, proc, vel_field,
                debug, verbose)
        
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
            args0 = (nx, ny, nz, N, grids, ub, vb, wb, vt['data'], rho, drho,
                    base['data'], top['data'], column['data'], wgt_o, wgt_c,
                    wgt_s0, wgt_b, wgt_w0, continuity_cost0, smooth_cost, dx,
                    dy, dz, sub_beam, finite_scheme, fill_value, proc,
                    vel_field, debug, verbose)
            
            # Call the SciPy solver
            res = minimize(f, x0, args=args0, method=method, jac=jac,
                           hess=None, hessp=None, bounds=None,
                           constraints=None, options=opts)
            
            # Unpack the results from the first pass
            x0 = res.x
            x0 = x0.astype(np.float64)
            
            # Make sure to set the vertical velocity field to 0 m/s
            # everywhere after the first pass
            x0[2*N:3*N] = 0.0
        
        # The debugging flag set to True will time the minimization
        # process
        if debug:
            t0 = time.clock()
        
        # Call the SciPy solver to perform the retrieval of the full 3-D
        # wind field     
        res = minimize(f, x0, args=args, method=method, jac=jac, hess=None,
                       hessp=None, bounds=None, constraints=None,
                       options=opts)
        
        if debug:
            t1 = time.clock()
            print 'The minimization took %i seconds' %(t1 - t0)
        
        # Unpack the full 3-D wind field results
        xopt = res.x
        xopt = xopt.astype(np.float64)
       
    else:
        raise ValueError('Unsupported technique and solver combination')
    
    
    # Create dictionaries for wind field
    u = get_metadata('u_wind')
    v = get_metadata('v_wind')
    w = get_metadata('w_wind')
    
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
    
    # Create winds (grid) object
    fields = {'u_wind': u, 'v_wind': v, 'w_wind': w}
    axes = {}
    metadata = {}
    
    return Grid(fields, axes, metadata)
