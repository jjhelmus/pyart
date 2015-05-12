
import numpy as np
from ..io.common import radar_coords_to_cart


def new_mapper(radar, grid_shape, grid_limits):
    # TODO multiple radars
    # TODO multiple fields
    # TODO ROI function

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

    grid_sum = np.zeros((nz, ny, nx), dtype=np.float32)
    grid_wsum = np.zeros((nz, ny, nx), dtype=np.float32)

    for nray in range(radar.nrays):

        print "nray", nray
        elevation = radar.elevation['data'][nray]
        azimuth = radar.azimuth['data'][nray]

        for ngate in range(radar.ngates):

            _range = radar.range['data'][ngate] / 1000.
            value = radar.fields['reflectivity']['data'][nray, ngate]
            if np.ma.is_masked(value):
                continue

            x, y, z = radar_coords_to_cart(_range, azimuth, elevation)
            # TODO apply radar displayment from grid origin

            #roi = roi_func(z, y, x)  # Region of interest for the radar gate
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

    mweight = np.ma.masked_equal(grid_wsum, 0)
    msum = np.ma.masked_array(grid_sum, mweight.mask)
    grid = msum / mweight
    return grid
