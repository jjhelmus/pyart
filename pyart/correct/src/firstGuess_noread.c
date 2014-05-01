/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**      This routine creates a firstGuess radial velocity field given a
**      sounding or VAD. Assumes standard atmosphere refraction (4Rearth/3) 
**      and extrapolates sounding data to all radar bins (assuming the wind
**      is horizontally uniform).
**
**  DEVELOPER:
**  Curtis N. James     08 Dec 98
**  Jonathan J. Helmus Refactor May 2014
**
**
**  HISTORY:
**  An elaboration of the NSSL-Eilts algorithm.
**
**
**
*/
#include "helpers.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <rsl.h>


/******************************
 * Private data and functions *
 ******************************/

typedef struct Sounding {
    float data[1000][5];
    int numLevs;
    float meanShearU;
    float meanShearV;
} Sounding;

/* Load sounding data from arrays into Sounding structure.
 * Filters out data with shear greater than maxshear. */
static void load_sounding(
    float *height_array, float *speed_array, float *direction_array, 
    int nlevels, Sounding *sound)
{
    int i, idx;
    i=1;

    sound->data[0][0]=0.0;
    sound->data[0][1]=0.0;
    sound->data[0][2]=0.0;
    sound->data[0][3]=0.0;
    sound->data[0][4]=0.0;
    sound->meanShearU = 0;
    sound->meanShearV = 0;

    /* store sounding data which meets shear conditions */
    for (idx = 1; idx < nlevels; ++idx) {
        /* Height (m) */
        sound->data[i][0] = height_array[idx]; 
        /* U REMOVED - (negative)*/
        sound->data[i][1] = sin(direction_array[idx] * PI / 180.0) * 
                        speed_array[idx];
        /* V REMOVED - (negative)*/
        sound->data[i][2] = cos(direction_array[idx] * PI / 180.0) *
                        speed_array[idx];
        /* Shear U */
        sound->data[i][3] = (sound->data[i][1] - sound->data[i-1][1]) /
                        (sound->data[i][0] - sound->data[i-1][0]);
        /* Shear V */
        sound->data[i][4] = (sound->data[i][2] - sound->data[i-1][2]) /
                        (sound->data[i][0] - sound->data[i-1][0]);
        if (fabs(sound->data[i][3]) <= MAXSHEAR && 
            fabs(sound->data[i][4]) <= MAXSHEAR) {
            i++;
            sound->meanShearU = sound->meanShearU + sound->data[i][3];
            sound->meanShearV = sound->meanShearV + sound->data[i][4];
        }
    } 
    /* Force wind at ground to be same as that of first level: */
    sound->data[0][1] = sound->data[1][1];
    sound->data[0][2] = sound->data[1][2];
    sound->data[1][3] = 0.0;
    sound->data[1][4] = 0.0;
    
    /* Calculate mean shear: */
    if (i>2) {
        sound->meanShearU = sound->meanShearU / (sound->data[i][0] - sound->data[1][0]);
        sound->meanShearV = sound->meanShearV / (sound->data[i][0] - sound->data[1][0]);
    } else {
        sound->meanShearU=0.0;
        sound->meanShearV=0.0;
    }
    sound->numLevs = i-1;
    return;
}

/* Calculate direction and speed from wind parameters */
static void calc_wind_dir(
    float ua_data_i[5], float shearU, float shearV, float height,
    float *wind, float *dir)
{
    float U, V, offset;
    U = ua_data_i[1] + shearU * (height - ua_data_i[0]);
    V = ua_data_i[2] + shearV * (height - ua_data_i[0]);
    *wind = sqrt(pow(U,2) + pow(V,2));
    if (SIGN < 0) 
        offset=0.0;
    else 
        offset=PI;
    if (U>=0)
        *dir = (acos(V / *wind) + offset) * 180 / PI;
    else 
        *dir = (offset - acos(V / *wind)) * 180 / PI;
    return;
}

/* Calculate the wind speed and direction by interpolating from the 
 * nearest measurements in the sounding data. */
static void interpolate_wind(Sounding *sound,float height, float *wind, float *dir)
{
    static int index=0; 
    if (height >= sound->data[0][0] && height <= sound->data[sound->numLevs-1][0]) {
        /* height is between two soundings measurments, 
         * find and interpolate from lower measurement */
        while(1) {
            if (height >= sound->data[index][0] && height < sound->data[index+1][0]) {
                calc_wind_dir(sound->data[index], sound->data[index+1][3], 
                              sound->data[index+1][4], height, wind, dir);
                break;
            } else if (height < sound->data[index][0]) {
                index--;
            } else {
                index++;
            }
        }
        
    } else if (height > sound->data[sound->numLevs-1][0]) {
        /* Height is above the heightest sounding measurement, extrapolate
         * from the heighest measurments using mean shear. */
        calc_wind_dir(sound->data[sound->numLevs-1], sound->meanShearU, sound->meanShearV,
                      height, wind, dir);
    }
}



void firstGuessNoRead(
    Volume* soundVolume, float missingVal,
    float *height_array, float *speed_array, float *direction_array, 
    int nlevels, int VAD_time, unsigned short* sounding) 
{
    int alt, i, sweepIndex, currIndex;
    unsigned short flag = 0;
    float ke, dRdz, height, rnge, elev, az, start_range, h_range, gate_size,
          wind, wind_val_rv, dir, ang;
    Sweep *sweep;
    Sounding sound;

    load_sounding(height_array, speed_array, direction_array, nlevels, &sound);

    dRdz=-39.2464; /* Standard Atmosphere refractivity gradient in km^-1 */
    ke = 1/(1 + A * dRdz * pow(10,-6)); /* Doviak and Zrnic, 1993 */
    alt = soundVolume->sweep[0]->ray[0]->h.alt;
    wind = dir = 0.0;
    
    /* Create a first-guess velocity field using the sounding
     ** profile, assuming standard atmospheric refraction. */
    for(sweepIndex = 0; sweepIndex < soundVolume->h.nsweeps; sweepIndex++) {
        
        sweep = soundVolume->sweep[sweepIndex]; 
        start_range = sweep->ray[0]->h.range_bin1;
        gate_size = sweep->ray[0]->h.gate_size;
        elev = PI * (sweep->ray[0]->h.elev) / 180.0;
        
        for(i = 0; i < (sweep->ray[0]->h.nbins); i++) {
            
            /* To print out a range circle of radial velocity values: */
            rnge = start_range + i * gate_size + gate_size / 2.0;
            height = sqrt(pow(rnge,2) + pow(ke*A*1000, 2) + 
                          2*rnge*ke*A*1000*sin(elev)) - ke*A*1000 + alt;
            h_range = ke*A*1000*asin(rnge*cos(elev) / (ke*A*1000 + height - alt));
            ang = atan(cos(elev) * sin(elev+h_range/ke/A/1000) * 
                pow(cos(elev + h_range/ke/A/1000), -2) ); /* atan (dh/ds) */
            for(currIndex = 0; currIndex < sweep->h.nrays; currIndex++) {
                
                interpolate_wind(&sound, height, &wind, &dir);
                if (wind>=0.0 && dir>=0.0) {
                    az = PI * (sweep->ray[currIndex]->h.azimuth) / 180.0;
                    wind_val_rv = wind * cos(ang) * cos(PI*dir/180.0 - az); 
                    ray_set(sweep->ray[currIndex], i, wind_val_rv);
                    flag=1;
                } else {
                    ray_set(sweep->ray[currIndex], i, missingVal);
                }               
            }
        }
    }
    if (sound.numLevs==0) {
        flag=0;
    }
    if (flag) *sounding=1;
    return;
}





