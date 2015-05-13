
import numpy as np
from ..io.common import radar_coords_to_cart
from ..graph.common import corner_to_point
from ._map_gates_to_grid import GateMapper
from ._map_gates_to_grid import ConstantROI, DistBeamROI, DistROI

from .grid_mapper import _gen_roi_func_constant
from .grid_mapper import _gen_roi_func_dist
from .grid_mapper import _gen_roi_func_dist_beam

#@profile
def new_mapper(
        radars, grid_shape, grid_limits,

        grid_origin=None, grid_origin_alt=None,
        roi_func='dist_beam',
        constant_roi=500., z_factor=0.05, xy_factor=0.02,
        min_radius=500.0, h_factor=1.0, nb=1.5, bsp=1.0,
        ):
    # TODO multiple radars
    # TODO multiple fields
    # TODO radar/grid displacement

    radar = radars[0]
    if grid_origin is None:
        lat = float(radars[0].latitude['data'])
        lon = float(radars[0].longitude['data'])
        grid_origin = (lat, lon)

    if grid_origin_alt is None:
        grid_origin_alt = float(radars[0].altitude['data'])

    offsets = []
    radar_lat = float(radar.latitude['data'])
    radar_lon = float(radar.longitude['data'])
    x_disp, y_disp = corner_to_point(grid_origin, (radar_lat, radar_lon))
    z_disp = float(radar.altitude['data']) - grid_origin_alt
    offsets.append((z_disp, y_disp, x_disp))
    print offsets
    if not hasattr(roi_func, '__call__'):
        if roi_func == 'constant':
            roi_func = ConstantROI(constant_roi)
        elif roi_func == 'dist':
            roi_func = DistROI(z_factor, xy_factor, min_radius, offsets)
        elif roi_func == 'dist_beam':
            roi_func = DistBeamROI(h_factor, nb, bsp, min_radius, offsets)
        else:
            raise ValueError('unknown roi_func: %s' % roi_func)

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

    grid_sum = np.zeros((nz, ny, nx, 2), dtype=np.float32)
    grid_wsum = np.zeros((nz, ny, nx, 2), dtype=np.float32)

    field_data = np.empty((radar.nrays, radar.ngates, 2), dtype='float32')
    field_mask = np.empty((radar.nrays, radar.ngates, 2), dtype='uint8')
    field_data[:, :, 0] = np.ma.getdata(radar.fields['reflectivity']['data'])
    field_mask[:, :, 0] = np.ma.getmaskarray(
        radar.fields['reflectivity']['data'])
    field_data[:, :, 1] = np.ma.getdata(
        radar.fields['normalized_coherent_power']['data'])
    field_mask[:, :, 1] = np.ma.getmaskarray(
        radar.fields['normalized_coherent_power']['data'])

    z_offset, y_offset, x_offset = offsets[0]
    gatemapper = GateMapper(x_step, y_step, z_step, x_start, y_start, z_start,
                            nx, ny, nz, grid_sum, grid_wsum)
    gatemapper.map_gates_to_grid(
        radar.elevation['data'], radar.azimuth['data'], radar.range['data'],
        field_data, field_mask.astype(np.uint8),
        x_offset, y_offset, x_offset,
        roi_func)

    #import IPython; IPython.embed()
    mweight = np.ma.masked_equal(grid_wsum, 0)
    msum = np.ma.masked_array(grid_sum, mweight.mask)
    grid = msum / mweight
    return grid


def map_gates_to_grid_old(
        grid_sum, grid_wsum,
        nrays, ngates,
        elevations, azimuths, ranges,
        field,
        x_start, x_step, y_start, y_step, z_start, z_step, nx, ny, nz
        ):

    for nray in range(nrays):

        print "nray", nray
        elevation = elevations[nray]
        azimuth = azimuths[nray]

        for ngate in range(ngates):

            _range = ranges[ngate] / 1000.
            value = field[nray, ngate]
            if np.ma.is_masked(value):
                continue

            x, y, z = radar_coords_to_cart(_range, azimuth, elevation)
            # TODO apply radar displayment from grid origin

            #roi = roi_func(z, y, x)  # Region of interest for the radar gate
            # TODO true ROI calculation
            roi = 2000

            x_min = max(np.ceil((x - roi - x_start) / x_step), 0)
            x_max = min(np.floor((x + roi - x_start) / x_step), nx-1)
            if x_min > nx-1 or x_max < 0:
                continue

            y_min = max(np.ceil((y - roi - y_start) / y_step), 0)
            y_max = min(np.floor((y + roi - y_start) / y_step), ny-1)
            if y_min > ny-1 or y_max < 0:
                continue

            z_min = max(np.ceil((z - roi - z_start) / z_step), 0)
            z_max = min(np.floor((z + roi - z_start) / z_step), nz-1)
            if z_min > nz-1 or z_max < 0:
                continue

            grid_x_indices = range(int(x_min), int(x_max)+1)
            grid_y_indices = range(int(y_min), int(y_max)+1)
            grid_z_indices = range(int(z_min), int(z_max)+1)

            for xi in grid_x_indices:
                for yi in grid_y_indices:
                    for zi in grid_z_indices:

                        xg = x_start + x_step * xi
                        yg = y_start + y_step * yi
                        zg = z_start + z_step * zi

                        dist = np.sqrt((xg-x)**2 + (yg-y)**2 + (zg-z)**2)
                        weight = np.exp(-(dist*dist) / (2*roi*roi)) + 1e-5

                        grid_sum[zi, yi, xi] += weight * value
                        grid_wsum[zi, yi, xi] += weight

