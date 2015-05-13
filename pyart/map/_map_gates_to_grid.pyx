# cython: profile=False
from libc.math cimport sqrt, exp, ceil, floor, sin, cos, tan, asin

cimport cython
from cython.view cimport array as cvarray


# This definition can be added to a .pxd file so others can defined fast 
# RoI functions
cdef class ROIFunction:

    cpdef float get_roi(self, float z, float y, float x):
        return 0


cdef class ConstantROI(ROIFunction):
    
    cdef float constant_roi

    def __init__(self, float constant_roi):
        self.constant_roi = constant_roi

    cpdef float get_roi(self, float z, float y, float x):
        return self.constant_roi


cdef class DistROI(ROIFunction):
    
    cdef float z_factor, xy_factor, min_radius
    cdef int num_offsets
    cdef float[:, :] offsets

    def __init__(self, z_factor, xy_factor, min_radius, offsets):
        self.z_factor = z_factor
        self.xy_factor = xy_factor
        self.min_radius = min_radius

        self.num_offsets = len(offsets)
        # does this array need to be explicitly de-allocated when the
        # class instance is removed?
        self.offsets = cvarray(
            shape=(self.num_offsets, 3), itemsize=sizeof(float), format='f')
        
        for i, (z_offset, y_offset, x_offset) in enumerate(offsets):
            self.offsets[i, 0] = z_offset
            self.offsets[i, 1] = y_offset
            self.offsets[i, 2] = x_offset

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef float get_roi(self, float z, float y, float x):
        
        cdef float min_roi, roi, z_offset, y_offset, x_offset
        cdef int i

        min_roi =  999999999.0
        for i in range(self.num_offsets):
            z_offset = self.offsets[i, 0]
            y_offset = self.offsets[i, 1]
            x_offset = self.offsets[i, 2]
            roi = (self.z_factor * (z - z_offset) + self.xy_factor *
                   sqrt((x - x_offset)**2 + (y - y_offset)**2) +
                   self.min_radius)
            if roi < min_roi:
                min_roi = roi

        return min_roi


cdef class DistBeamROI(ROIFunction):
    
    cdef float h_factor, min_radius, beam_factor
    cdef int num_offsets
    cdef float[:, :] offsets

    def __init__(self, h_factor, nb, bsp, min_radius, offsets):
        self.h_factor = h_factor
        self.min_radius = min_radius
        self.beam_factor = tan(nb * bsp * 3.141592653589793 / 180.)

        self.num_offsets = len(offsets)
        # does this array need to be explicitly de-allocated when the
        # class instance is removed?
        self.offsets = cvarray(
            shape=(self.num_offsets, 3), itemsize=sizeof(float), format='f')
        
        for i, (z_offset, y_offset, x_offset) in enumerate(offsets):
            self.offsets[i, 0] = z_offset
            self.offsets[i, 1] = y_offset
            self.offsets[i, 2] = x_offset

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef float get_roi(self, float z, float y, float x):
        
        cdef float min_roi, roi, z_offset, y_offset, x_offset
        cdef int i

        min_roi =  999999999.0
        for i in range(self.num_offsets):
            z_offset = self.offsets[i, 0]
            y_offset = self.offsets[i, 1]
            x_offset = self.offsets[i, 2]
            roi = (self.h_factor * ((z - z_offset) / 20.0) +
                   sqrt((y - y_offset)**2 + (x - x_offset)**2) *
                   self.beam_factor + self.min_radius)
            if roi < min_roi:
                min_roi = roi

        return min_roi



cdef class GateMapper:
    """
    """

    cdef float x_step, y_step, z_step
    cdef float x_start, y_start, z_start
    cdef int nx, ny, nz
    cdef float[:, :, ::1] grid_sum
    cdef float[:, :, ::1] grid_wsum

    def __init__(self, float x_step, float y_step, float z_step, 
                 float x_start, float y_start, float z_start,
                 int nx, int ny, int nz,
                 float [:, :, ::1] grid_sum, float[:, :, ::1] grid_wsum):
        """ initialize. """
        self.x_step = x_step
        self.y_step = y_step
        self.z_step = z_step
        self.x_start = x_start
        self.y_start = y_start
        self.z_start = z_start
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.grid_sum = grid_sum
        self.grid_wsum = grid_wsum
        return

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def map_gates_to_grid(self,
            float[::1] elevations, float[::1] azimuths, float[::1] ranges,
            float[:, ::1] field_data, char[:, ::1] field_mask,
            float x_offset, float y_offset, float z_offset,
            ROIFunction roi_func):

        cdef int nrays, ngates
        cdef float elevation, azimuth, r, s 
        cdef float pi = 3.141592653589793
        # 4/3 earths radius of 6371 km in meters
        cdef float R = 8494666.66666667
        cdef float value, roi
        cdef float x, y, z

        nrays = len(elevations)
        ngates = len(ranges)
        for nray in range(nrays):
            # elevation and azimuth angles in radians
            elevation = elevations[nray] * pi / 180.0
            azimuth = azimuths[nray] * pi / 180.0
            for ngate in range(ngates):

                if field_mask[nray, ngate]:
                    continue
                value = field_data[nray, ngate]

                # calculate cartesian coordinate assuming 4/3 earth radius 
                r = ranges[ngate] 
                z = (r**2 + R**2 + 2.0*r*R*sin(elevation))**0.5 - R
                s = R * asin(r * cos(elevation) / (R + z))  # arc length in m.
                x = s * sin(azimuth)
                y = s * cos(azimuth)

                # Region of influence 
                roi = roi_func.get_roi(z, y, x)

                # shift positions so that grid starts at 0
                # TODO apply radar displayment from grid origin
                x = x - self.x_start + x_offset
                y = y - self.y_start + y_offset
                z = z - self.z_start + z_offset

                # TODO dynamic ROI
                #roi = 2000.
                self.map_gate(x, y, z, roi, value)

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int map_gate(self, float x, float y, float z, float roi, float value):
        """ Map a single gate to the grid. """
        
        cdef float xg, yg, zg, dist
        cdef int x_min, x_max, y_min, y_max, z_min, z_max
        cdef int xi, yi, zi

        x_min = find_min(x, roi, self.x_step)
        x_max = find_max(x, roi, self.x_step, self.nx)
        if x_min > self.nx-1 or x_max < 0:
            return 0

        y_min = find_min(y, roi, self.y_step)
        y_max = find_max(y, roi, self.y_step, self.ny)
        if y_min > self.ny-1 or y_max < 0:
            return 0

        z_min = find_min(z, roi, self.z_step)
        z_max = find_max(z, roi, self.z_step, self.nz)
        if z_min > self.nz-1 or z_max < 0:
            return 0
        
        for xi in range(x_min, x_max+1):
            for yi in range(y_min, y_max+1): 
                for zi in range(z_min, z_max+1):
                    xg = self.x_step * xi
                    yg = self.y_step * yi
                    zg = self.z_step * zi
                    dist = sqrt((xg-x)**2 + (yg-y)**2 + (zg-z)**2)
                    if roi == 0:
                        weight = 1e-5
                    else:
                        weight = exp(-(dist*dist) / (2*roi*roi)) + 1e-5
                    self.grid_sum[zi, yi, xi] += weight * value
                    self.grid_wsum[zi, yi, xi] += weight
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
