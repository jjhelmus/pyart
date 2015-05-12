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
    cdef float value, roi
    cdef float x, y, z
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

            # shift positions so that grid starts at 0
            # TODO apply radar displayment from grid origin
            x = x - x_start
            y = y - y_start
            z = z - z_start

            #roi = roi_func(z, y, x)  # Region of interest for the radar gate
            # TODO true ROI calculation
            roi = 2000.
            map_gate(x, y, z, roi, value,
                     x_step, y_step, z_step, nx, ny, nz, grid_sum, grid_wsum)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int map_gate(float x, float y, float z, float roi, float value,
                  float x_step, float y_step, float z_step, 
                  int nx, int ny, int nz,
                  float [:, :, ::1] grid_sum, float[:, :, ::1] grid_wsum):
    """ Map a single point to the grid. """
    
    cdef float xg, yg, zg, dist
    cdef int x_min, x_max, y_min, y_max, z_min, z_max
    cdef int xi, yi, zi

    x_min = find_min(x, roi, x_step)
    x_max = find_max(x, roi, x_step, nx)
    if x_min > nx-1 or x_max < 0:
        return 0

    y_min = find_min(y, roi, y_step)
    y_max = find_max(y, roi, y_step, ny)
    if y_min > ny-1 or y_max < 0:
        return 0

    z_min = find_min(z, roi, z_step)
    z_max = find_max(z, roi, z_step, nz)
    if z_min > nz-1 or z_max < 0:
        return 0
    
    for xi in range(x_min, x_max+1):
        for yi in range(y_min, y_max+1): 
            for zi in range(z_min, z_max+1):
                xg = x_step * xi
                yg = y_step * yi
                zg = z_step * zi
                dist = sqrt((xg-x)**2 + (yg-y)**2 + (zg-z)**2)
                if roi == 0:
                    weight = 1e-5
                else:
                    weight = exp(-(dist*dist) / (2*roi*roi)) + 1e-5
                grid_sum[zi, yi, xi] += weight * value
                grid_wsum[zi, yi, xi] += weight
    return 1

@cython.cdivision(True)
cdef int find_min(float a, float roi, float step):
    cdef int a_min
    if step == 0:
        return 0
    a_min = <int>ceil((a - roi) / step)
    if a_min < 0:
        a_min = 0
    return a_min


@cython.cdivision(True)
cdef int find_max(float a, float roi, float step, int na):
    cdef int a_max
    if step == 0:
        return 0
    a_max = <int>floor((a + roi) / step)
    if a_max > na-1:
        a_max = na-1
    return a_max
