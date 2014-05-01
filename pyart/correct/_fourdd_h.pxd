"""
Cython wrapper around University of Washington 4DD code.
"""

cimport pyart.io._rsl_h as _rsl_h

cdef extern from "fourdd_jjh.h":

    void unfoldVolume(_rsl_h.Volume* rvVolume, _rsl_h.Volume* soundVolume,
                      _rsl_h.Volume* lastVolume, float missingVal,
                      unsigned short rm, unsigned short* success)


cdef extern from "prepVolume.h":

    void prepVolume(_rsl_h.Volume* DBZVolume, _rsl_h.Volume* rvVolume,
                    float missingVal)

cdef extern from "sounding_to_volume.h":
    
    int sounding_to_volume(_rsl_h.Volume* soundVolume, float missingVal,
                           float *height_array, float *speed_array,
                           float *direction_array, int nlevels, 
                           float maxshear, int sign)


