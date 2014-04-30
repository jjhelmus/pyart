/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**      This routine averages the values in a range and azimuth window of a
**      sweep and computes the standard deviation.
**
**  DEVELOPER:
**	Curtis N. James    26 Jan 1999
**
**
**
*/
#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */

float window(Sweep *rv_sweep, int startray, 
	int endray, int firstbin, int lastbin, float missingVal,
    unsigned short* success) 
{
    /* Calculate the mean for gates within a window of a sweep */
    int num, currIndex, rangeIndex, numRays, numBins;
    float val, sum, sumsq, ave, NyqVelocity, std;
     
    *success=0;
    NyqVelocity = rv_sweep->ray[0]->h.nyq_vel;
    numRays = rv_sweep->h.nrays;
    numBins = rv_sweep->ray[0]->h.nbins;

    /* Now, sum the data in the window region between startray, 
     **  endray, firstbin, lastbin. */
    num=0;
    sum=0.0;
    sumsq=0.0;
       
    if (firstbin>=numBins || lastbin>=numBins || firstbin<0 || lastbin<0)
        return missingVal;
	for (rangeIndex=firstbin; rangeIndex<=lastbin; rangeIndex++) {
        if (startray>endray){
            for (currIndex=startray; currIndex<numRays; currIndex++) {
	            val = ray_val(rv_sweep->ray[currIndex], rangeIndex);
	            if (val!=missingVal) {
	                num++;
	                sum+=val;
	                sumsq+=(val*val);
	            }
	        }
            for (currIndex=0; currIndex<=endray; currIndex++) {
	            val = ray_val(rv_sweep->ray[currIndex], rangeIndex);
	            if (val!=missingVal) {
	                num++;
	                sum+=val;
	                sumsq+=(val*val);
	            }
	        }
        } else {
            for (currIndex=startray; currIndex<=endray; currIndex++) { 
	            val = ray_val(rv_sweep->ray[currIndex], rangeIndex);
	            if (val!=missingVal) {
	                num++;
	                sum+=val;
	                sumsq+=(val*val);
	            }
	        }
        }
    }
    if (num>=MINGOOD) {
        ave=sum/num;
        std=sqrt(fabs((sumsq-(sum*sum)/num)/(num-1)));
        if (std<=STDTHRESH*NyqVelocity) *success=1;
    } else {
        ave=missingVal;
        std=0.0;
        *success=1;
    }
    return ave; 
}
