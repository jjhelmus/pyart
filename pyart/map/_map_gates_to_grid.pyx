from libc.math cimport sqrt, exp, ceil, floor, sin, cos, asin

cimport cython
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def map_gates_to_grid(
        float [:, :, ::1] grid_sum, float[:, :, ::1] grid_wsum,
        int nrays, int ngates,
        float[::1] elevations, float [::1] azimuths, float[::1] ranges,
        float[:, ::1] field_data, char[:, ::1] field_mask,
        float x_start, float x_step, 
        float y_start, float y_step,
        float z_start, float z_step,
        int nx, int ny, int nz
        ):

    cdef float elevation, azimuth, _range
    cdef float weight, dist, value, roi
    cdef float xg, yg, zg
    cdef float x, y, z
    cdef int xi, yi, zi
    cdef int x_min, x_max
    cdef int y_min, y_max
    cdef int z_min, z_max
    cdef float pi=3.141592653589793
    for nray in range(nrays):

        elevation = elevations[nray]
        azimuth = azimuths[nray]

        for ngate in range(ngates):

            _range = ranges[ngate] / 1000.
            
            if field_mask[nray, ngate]:
                continue
            value = field_data[nray, ngate]

            # radar_coords_to_cart inline function???
            #x, y, z = radar_coords_to_cart(_range, azimuth, elevation)
            # elevation angle in radians.
            theta_e = elevation * pi / 180.0       
            # azimuth angle in radians.
            theta_a = azimuth * pi / 180.0        
            # effective radius of earth in meters.
            R = 6371.0 * 1000.0 * 4.0 / 3.0     
            # distances to gates in meters.
            r = _range * 1000.0                    
            z = (r ** 2 + R ** 2 + 2.0 * r * R * sin(theta_e)) ** 0.5 - R
            s = R * asin(r * cos(theta_e) / (R + z))  # arc length in m.
            x = s * sin(theta_a)
            y = s * cos(theta_a)

            # TODO apply radar displayment from grid origin

            #roi = roi_func(z, y, x)  # Region of interest for the radar gate
            # TODO true ROI calculation
            roi = 2000.

            x_min = <int>ceil((x - roi - x_start) / x_step)
            if x_min < 0:
                x_min = 0
            x_max = <int>floor((x + roi - x_start) / x_step)
            if x_max > nx-1:
                x_max = nx-1
            if x_min > nx-1 or x_max < 0:
                continue

            y_min = <int>ceil((y - roi - y_start) / y_step)
            if y_min < 0:
                y_min = 0
            y_max = <int>floor((y + roi - y_start) / y_step)
            if y_max > ny-1:
                y_max = ny-1
            if y_min > ny-1 or y_max < 0:
                continue

            if z_step == 0:
                z_min = 0
                z_max = 0
            else:
                z_min = <int>ceil((z - roi - z_start) / z_step)
                if z_min < 0:
                    z_min = 0
                z_max = <int>floor((z + roi - z_start) / z_step)
                if z_max > nz-1:
                    z_max = nz-1
                if z_min > nz-1 or z_max < 0:
                    continue
            
            for xi in range(x_min, x_max+1):
                for yi in range(y_min, y_max+1): 
                    for zi in range(z_min, z_max+1):

                        xg = x_start + x_step * xi
                        yg = y_start + y_step * yi
                        zg = z_start + z_step * zi

                        dist = sqrt((xg-x)**2 + (yg-y)**2 + (zg-z)**2)
                        if roi == 0:
                            weight = 1e-5
                        else:
                            weight = exp(-(dist*dist) / (2*roi*roi)) + 1e-5

                        grid_sum[zi, yi, xi] += weight * value
                        grid_wsum[zi, yi, xi] += weight

