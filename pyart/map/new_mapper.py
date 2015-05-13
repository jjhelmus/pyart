
import numpy as np
from ..io.common import radar_coords_to_cart
from ..config import get_fillvalue, get_field_name
from ..graph.common import corner_to_point
from ._map_gates_to_grid import GateMapper
from ._map_gates_to_grid import ROIFunction, ConstantROI, DistBeamROI, DistROI

from .grid_mapper import _gen_roi_func_constant
from .grid_mapper import _gen_roi_func_dist
from .grid_mapper import _gen_roi_func_dist_beam

#@profile
def new_mapper(radars, grid_shape, grid_limits, grid_origin=None,
                grid_origin_alt=None, fields=None, refl_filter_flag=True,
                refl_field=None, max_refl=None, map_roi=True,
                weighting_function='Barnes', toa=17000.0, copy_field_data=True,
                algorithm='kd_tree', leafsize=10., roi_func='dist_beam',
                constant_roi=500., z_factor=0.05, xy_factor=0.02,
                min_radius=500.0, h_factor=1.0, nb=1.5, bsp=1.0):
    # Supported
    # * grid_shape
    # * grid_limits
    # * grid_origin
    # * grid_origin_alt
    # * fields
    # * refl_field
    # * toa
    # * refl_filter_flag
    # * max_refl
    # * roi_func
    # * constant_roi, z_factor, xy_factor, min_radius, h_factor, nb, bsp

    # TODO Not yet implemented
    # * map_roi

    # Depreciate
    # * copy_field_data (don't want to implement)
    # * algorithm
    # * leafsize

    # parse the parameters
    if max_refl is None:
        max_refl = np.finfo('float32').max
    if weighting_function.upper() == 'CRESSMAN':
        cy_weighting_function = 1
    elif weighting_function.upper() == 'BARNES':
        cy_weighting_function = 0
    else:
        raise ValueError('unknown weighting_function')

    # find the grid origin if not given
    if grid_origin is None:
        lat = float(radars[0].latitude['data'])
        lon = float(radars[0].longitude['data'])
        grid_origin = (lat, lon)

    if grid_origin_alt is None:
        grid_origin_alt = float(radars[0].altitude['data'])

    # fields which should be mapped, None for fields which are in all radars
    if fields is None:
        fields = set(radars[0].fields.keys())
        for radar in radars[1:]:
            fields = fields.intersection(radar.fields.keys())
        fields = list(fields)
    nfields = len(fields)

    # find the reflectivity field, check that it is mapped and
    # move it to the front of the fields list
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if refl_field not in fields:
        raise ValueError('reflectivity field not mapped')
    fields.insert(0, fields.pop(fields.index(refl_field)))


    # loop over the radars finding offsets from the origin
    offsets = []    # offsets from the grid origin, in meters, for each radar
    for radar in radars:
        radar_lat = float(radar.latitude['data'])
        radar_lon = float(radar.longitude['data'])
        x_disp, y_disp = corner_to_point(grid_origin, (radar_lat, radar_lon))
        z_disp = float(radar.altitude['data']) - grid_origin_alt
        offsets.append((z_disp, y_disp, x_disp))

    # unpack the grid parameters
    nz, ny, nx = grid_shape
    zr, yr, xr = grid_limits
    z_start, z_stop = zr
    y_start, y_stop = yr
    x_start, x_stop = xr

    if nz == 1:
        z_step = 0.
    else:
        z_step = (z_stop - z_start) / (nz - 1.)
    if ny == 1:
        y_step = 0.
    else:
        y_step = (y_stop - y_start) / (ny - 1.)
    if nx == 1:
        x_step = 0.
    else:
        x_step = (x_stop - x_start) / (nx - 1.)

    # Radius of influence function
    if not isinstance(roi_func, ROIFunction):
        if roi_func == 'constant':
            roi_func = ConstantROI(constant_roi)
        elif roi_func == 'dist':
            roi_func = DistROI(z_factor, xy_factor, min_radius, offsets)
        elif roi_func == 'dist_beam':
            roi_func = DistBeamROI(h_factor, nb, bsp, min_radius, offsets)
        else:
            raise ValueError('unknown roi_func: %s' % roi_func)

    print offsets
    print fields

    # prepare grid storage arrays
    grid_sum = np.zeros((nz, ny, nx, nfields), dtype=np.float32)
    grid_wsum = np.zeros((nz, ny, nx, nfields), dtype=np.float32)
    gatemapper = GateMapper(
        x_step, y_step, z_step, x_start, y_start, z_start,
        nx, ny, nz, grid_sum, grid_wsum)

    # project gates from each radar onto the grid
    for radar, radar_offset in zip(radars, offsets):
        shape = (radar.nrays, radar.ngates, nfields)
        field_data = np.empty(shape, dtype='float32')
        field_mask = np.empty(shape, dtype='uint8')

        for i, field in enumerate(fields):
            fdata = radar.fields[field]['data']
            field_data[:, :, i] = np.ma.getdata(fdata)
            field_mask[:, :, i] = np.ma.getmaskarray(fdata)

        z_offset, y_offset, x_offset = radar_offset
        gatemapper.map_gates_to_grid(
            radar.elevation['data'].astype('float32'),
            radar.azimuth['data'].astype('float32'),
            radar.range['data'].astype('float32'),
            field_data, field_mask, radar_offset, toa, roi_func,
            refl_filter_flag, max_refl, cy_weighting_function)


    # create and return the grid dictionary
    mweight = np.ma.masked_equal(grid_wsum, 0)
    msum = np.ma.masked_array(grid_sum, mweight.mask)
    grids = dict(
        [(f, msum[..., i] / mweight[..., i]) for i, f in enumerate(fields)])
    #if map_roi:
        #grids['ROI'] = roi
    return grids
