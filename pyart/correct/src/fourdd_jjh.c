/*
 * Unfold the doppler velocity data for a single radar volume using 
 * sounding data and/or a priviously unfolded volume.
 *  
 * Adapted from routines from:
 * 
 * UW Radial Velocity Dealiasing Algorithm
 * Four-Dimensional Dealiasing (4DD)
 * 
 * Developer by Curtis N. James     Jan-Feb 1999
 *
 * See UWASHINGTON_4DD_README file for license
 *
 * Adapted for use in Py-ART by Jonathan J. Helmus May 2014
 */


/*
**
**
**  DESCRIPTION:
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and a higher elevation sweep to unfold some of the bins
**  in each sweep. Then, it completes the unfolding, assuming spatial
**  continuity around each bin.
**
**
*/

/* TODO
 * 
 * - Introduce Dealias structure to hold sweeps, etc
 * - Remove use of constants, add these to structure and make arguments to
 *   unfoldVolume function.
 * - Refactor and clean.
 * - change == and != missingVal to better floating point comparisons
 */

#include "helpers.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */ 


/******************************
 * Private data and functions *
 ******************************/
typedef struct DealiasParams {
    float missingVal;
    float compthresh;
    float ckval;
} DealiasParams;


/* Structure for storing data on sweeps being examined for unfolding */
typedef struct SweepCollection {
    Sweep *vals;
    Sweep *rv;
    Sweep *last;
    Sweep *sound;
    Sweep *above;
    float NyqVelocity;
    float NyqInterval;
    int nrays;
    int nbins;
} SweepCollection;


/*
 * Berger and Albers filter functions
 */

/* Count the number of neighbors which are not missing value */
static int count_nonmissing_neighbors(
    Sweep *sweep, int currIndex, int i, float missingVal)
{
    int countindex = 0;
    int left, right, next, prev;

    if (currIndex == 0) 
        left = sweep->h.nrays-1;
    else 
        left = currIndex - 1;
    if (currIndex == sweep->h.nrays-1) 
        right = 0;
    else 
        right = currIndex+1;
    next = i+1;
    prev = i-1;
    /* Look at all bins adjacent to current bin in question: */
    if (i != 0) {
    if (ray_val(sweep->ray[left], prev) != missingVal)
        countindex++;
    if (ray_val(sweep->ray[currIndex], prev) != missingVal)
        countindex++;
    if (ray_val(sweep->ray[right], prev) != missingVal)
        countindex++;
    }
    if (ray_val(sweep->ray[left], i) != missingVal)
        countindex++;
    if (ray_val(sweep->ray[right], i) != missingVal)
        countindex++;
    if (i < sweep->ray[0]->h.nbins-1) {  
        if (ray_val(sweep->ray[left], next) != missingVal)
            countindex++;
        if (ray_val(sweep->ray[currIndex], next)!= missingVal)
            countindex++;
        if (ray_val(sweep->ray[right], next) != missingVal)
            countindex++;
    }
    return countindex;
}

/* Perform a 3x3 missing value Berger and Albers filter */
/* TODO remove use of 5 and 3 magic values */
static void berger_albers_filter(
    Sweep *sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS])
{
    int currIndex, i, count;
    float val;
    for (currIndex=0; currIndex<sweep->h.nrays; currIndex++) {
        for (i=DELNUM; i<sweep->ray[0]->h.nbins; i++) {
            val = ray_val(sweep->ray[currIndex], i); 
            if (val==missingVal) {
                GOOD[i][currIndex]=-1;
            } else {
                count = count_nonmissing_neighbors(
                    sweep, currIndex, i, missingVal);
                if (count >= 5)
                    GOOD[i][currIndex] = 0;  
                else if (i == (sweep->ray[0]->h.nbins - 1) && count >= 3)
                    GOOD[i][currIndex] = 0;  
                else
                    GOOD[i][currIndex] = -1; 
            }
        }
    }
} 

/*
 * Zero-ing GOOD array function
 */

/* Fill the GOOD array with zeros at all location where the value in sweep 
 * is not missing */
static void zero_good(
    Sweep *sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS])
{
    int currIndex, i;
    float val;
    for (currIndex=0; currIndex<sweep->h.nrays; currIndex++) {
        for (i=DELNUM; i<sweep->ray[0]->h.nbins; i++) {
            val = ray_val(sweep->ray[currIndex], i); 
            if (val==missingVal)
                GOOD[i][currIndex]=-1;
            else
                GOOD[i][currIndex]=0;
        }
    }
}

/*
 * Sweep initialization
 */

/* Initalize (fill) a sweep with a given value */
static void init_sweep(Sweep *sweep, float val)
{
    int currIndex, i;
    for (currIndex=0; currIndex<sweep->h.nrays; currIndex++) {
        for (i=DELNUM; i<sweep->ray[0]->h.nbins; i++) {
            ray_set(sweep->ray[currIndex], i, val);
        }
    }
}

/*
 * Continuity Dealising (initial dealias) functions 
 */

/* Finds the rayindex of the nearest ray in sweep2 to
 * the currIndex ray in sweep1. */
static int findRay(Sweep *sweep1, Sweep *sweep2, int currIndex) 
{
    int numRays, rayIndex1;
    float az0, az1, diffaz;
    float spacing;
    short direction, lastdir;
    
    numRays = sweep2->h.nrays;

    if (currIndex<numRays) 
        rayIndex1=currIndex;
    else 
        rayIndex1=numRays-1;
    
    az0 = sweep1->ray[currIndex]->h.azimuth;
    az1 = sweep2->ray[rayIndex1]->h.azimuth;
    if (az0 == az1) {
        return rayIndex1;
    } else {
        
        /* Since the beamwidth is not necessarily the spacing between rays: */
        spacing = fabs(sweep2->ray[0]->h.azimuth - sweep2->ray[50]->h.azimuth);
        if (spacing>180)
            spacing = 360.0 - spacing;
        spacing = spacing / 50.0;

        /* Compute the difference in azimuth between the two rays: */
        diffaz=az0-az1;
        if (diffaz>=180.0) 
            diffaz=diffaz-360.0;
        else if (diffaz<-180.0) 
            diffaz=diffaz+360.0;
       
        /* Get close to the correct index: */
        rayIndex1=rayIndex1+(int) (diffaz/spacing);
        if (rayIndex1>=numRays) 
            rayIndex1=rayIndex1-numRays;
        if (rayIndex1<0) 
            rayIndex1=numRays+rayIndex1;
        az1=sweep2->ray[rayIndex1]->h.azimuth;
        diffaz=az0-az1;
        if (diffaz>=180.0)
            diffaz=diffaz-360.0;
        else if 
            (diffaz<-180.0) diffaz=diffaz+360.0;

        /* Now add or subtract indices until the nearest ray is found: */
        if (diffaz>=0) 
            lastdir=1;
        else 
            lastdir=-1;
        while (fabs(diffaz)>spacing/2.0) {
	        if (diffaz>=0) {
	            rayIndex1++;
	            direction=1;
	        } else {
	            rayIndex1--;
	            direction=-1;
	        }
	        if (rayIndex1>=numRays) 
                rayIndex1=rayIndex1-numRays;
	        if (rayIndex1<0) 
                rayIndex1=numRays+rayIndex1;
	        az1=sweep2->ray[rayIndex1]->h.azimuth;
	        diffaz=az0-az1;
	        if (diffaz>=180.0) 
                diffaz=diffaz-360.0;
	        else if 
                (diffaz<-180.0) diffaz=diffaz+360.0;
	        if (direction!=lastdir) 
                break;
	        else 
                lastdir=direction;
       }
       return rayIndex1;
     }
}

/* Minize the difference between val and cval by adding or substracting 
 * the Nyquist Interval, update val. */
static float min_diff(float *val, float cval, float NyqVelocity)
{
    int numtimes; 
    numtimes=0;
    while (fabs(cval-*val) > 0.99999 * NyqVelocity && numtimes <= MAXCOUNT) {
        numtimes++;
        if (*val > cval)
            *val-=2 * NyqVelocity;
        else
            *val+=2 * NyqVelocity;
    }
    return fabs(cval-*val);
}

/* Set a value in GOOD and the corresponding ray value in a sweep */
static void set_good_ray(
    short GOOD[][MAXRAYS], Sweep *sweep, 
    int ray_idx, int bin_idx, short good_val, float ray_val)
{
    ray_set(sweep->ray[ray_idx], bin_idx, ray_val);
    GOOD[bin_idx][ray_idx] = good_val;
    return;
}

/* Unfold bins where vertical and temporal continuity
 * produces the same value (i.e., bins where the radial velocity
 * agrees with both the previous sweep and previous volume within
 * a COMPTHRESH of the Nyquist velocity). This produces a number
 * of good data points from which to begin unfolding assuming spatial
 * continuity. */
static void continuity_dealias(
    SweepCollection *sweepc, DealiasParams *dp, 
    short GOOD[MAXBINS][MAXRAYS])
{
    int currIndex;
    int i, prevIndex = 0, abIndex = 0; 
    float val;
    float prevval, soundval, abval, diff, thresh;

    thresh = dp->compthresh * sweepc->NyqVelocity;

    for (currIndex=0; currIndex < sweepc->nrays; currIndex++) { 
        
        if (sweepc->last != NULL)
           prevIndex = findRay(sweepc->rv, sweepc->last, currIndex);
        if (sweepc->above != NULL)
            abIndex = findRay(sweepc->rv, sweepc->above, currIndex);
        
        for (i=DELNUM; i < sweepc->nbins; i++) { 
           
            if (GOOD[i][currIndex]!=0) continue;
            val = ray_val(sweepc->vals->ray[currIndex], i);
            if (val == dp->missingVal) continue;
            if (fabs(val) <= dp->ckval) continue;
            
            /* find velocities in related sweeps */
            prevval = soundval = abval = dp->missingVal;
            if (sweepc->last!=NULL)
                prevval=ray_val(sweepc->last->ray[prevIndex], i);
            if (sweepc->sound != NULL && sweepc->last == NULL)
                soundval=ray_val(sweepc->sound->ray[currIndex], i);
            if (sweepc->above != NULL)
                abval=ray_val(sweepc->above->ray[abIndex], i);
            
            if (sweepc->last == NULL) {
                if (soundval == dp->missingVal) continue;
                if (abval == dp->missingVal) {
                    diff = min_diff(&val, soundval, sweepc->NyqVelocity);
                    if (diff < thresh)
                        set_good_ray(GOOD, sweepc->rv, currIndex, i, 1, val);
                } else {
                    diff = min_diff(&val, abval, sweepc->NyqVelocity);
                    if (diff < thresh && fabs(soundval-val) < thresh)
                        set_good_ray(GOOD, sweepc->rv, currIndex, i, 1, val);
                }
            } else {
                if (prevval != dp->missingVal && abval != dp->missingVal) {
                    diff = min_diff(&val, prevval, sweepc->NyqVelocity);
                    if (diff < thresh && fabs(abval-val) < thresh)
                        set_good_ray(GOOD, sweepc->rv, currIndex, i, 1, val); 
                }
            }
        }
    }
}


/* Mark gates where all neighbors are GOOD == 0 with GOOD == -1 */
static void mark_neighborless(
    short GOOD[MAXBINS][MAXRAYS], int currIndex, int i, int nrays, int nbins)
{
    int next, prev, left, right;
    int countbins = 0;
    if (currIndex==0) left=nrays-1;
    else left=currIndex-1;
    if (currIndex==nrays-1) right=0;
    else right=currIndex+1;
    next = i+1;
    prev = i-1;
 
    if (i != 0) {
        if (GOOD[prev][left]==0)
            countbins++;
        if (GOOD[prev][currIndex]==0)
            countbins++;
        if (GOOD[prev][right]==0)
            countbins++;
    }
    if (GOOD[i][left]==0)
        countbins++;
    if (GOOD[i][right]==0)
        countbins++;
    if (i<nbins-1) {
        if (GOOD[next][left]==0)
            countbins++;
        if (GOOD[next][currIndex]==0)
            countbins++;
        if (GOOD[next][right]==0)
            countbins++;
    }
    if (countbins == 0)
        GOOD[i][currIndex]=-1;  
    return;
}

/* Record the location of a bin and ray index in the arrays */
static void record_gate(
    int countindex, int b_index, int r_index, 
    int binindex[8], int rayindex[8])
{
    binindex[countindex] = b_index;
    rayindex[countindex] = r_index;
    return;
}

/* Count the number of GOOD == 1 neighbors */
static int count_neighbors(
    short GOOD[MAXBINS][MAXRAYS],
    int currIndex, int i, int nrays, int nbins,
    int binindex[8], int rayindex[8]
    )
{
    int next, prev, left, right;
    int countindex = 0;
    
    if (currIndex==0) left = nrays-1;
    else left=currIndex-1;
    if (currIndex==nrays-1) right=0;
    else right=currIndex+1;
    next=i+1;
    prev=i-1;
 
    /* Look at all 8 bins adjacent to current bin in question: */
    if (i!=0) { 
        if (GOOD[prev][left]==1)
            record_gate(countindex++, prev, left, binindex, rayindex);
        if (GOOD[prev][currIndex]==1)
            record_gate(countindex++, prev, currIndex, binindex, rayindex);
        if (GOOD[prev][right]==1)
            record_gate(countindex++, prev, right, binindex, rayindex);
    }
    if (GOOD[i][left]==1)
        record_gate(countindex++, i, left, binindex, rayindex);
    if (GOOD[i][right]==1)
        record_gate(countindex++, i, right, binindex, rayindex);
    if (i<nbins-1) {  
        if (GOOD[next][left]==1)
            record_gate(countindex++, next, left, binindex, rayindex);
        if (GOOD[next][currIndex]==1)
            record_gate(countindex++, next, currIndex, binindex, rayindex);
        if (GOOD[next][right]==1)
            record_gate(countindex++, next, right, binindex, rayindex);
    }
    return countindex;
}


/* Unfold against all adjacent values where GOOD==1 */
static void unfold_adjacent(
    float val, Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval, 
    int i, int currIndex, int countindex, int loopcount,
    int binindex[8], int rayindex[8]
    )
{
    int in, out, numneg, numpos, l, numtimes;
    float diff;
    numtimes=0;
    while(val!=missingVal && GOOD[i][currIndex]==0) {
        numtimes++;
        in = out = numneg = numpos = 0;
        for (l=0; l<countindex; l++) {
            diff = ray_val(rv_sweep->ray[rayindex[l]], binindex[l]) - val;
            if (fabs(diff) < THRESH*NyqVelocity)
                in++;
            else
                out++;
            if (diff > NyqVelocity)
                numpos++;
            else if (diff < -NyqVelocity)
                numneg++; 
        }
        if (in>0 && out==0) {
            ray_set(rv_sweep->ray[currIndex], i, val);
            GOOD[i][currIndex]=1;
        } else {
            if (numtimes<=MAXCOUNT) {
                if ((numpos+numneg)<(in+out-(numpos+numneg))) {
                    if (loopcount>2) {
                        /* Keep the value after two passes through data. */
                        ray_set(rv_sweep->ray[currIndex], i, val);
                        GOOD[i][currIndex]=1;
                    }
                } else if (numpos>numneg) {
                    val=val+NyqInterval;
                } else if (numneg>numpos) {
                    val=val-NyqInterval;
                } else {
                    /* Save bin for windowing if unsuccessful after
                    ** four passes: */
                    if (loopcount>4) GOOD[i][currIndex]=-2;
                }
            } else {
                /* Remove bin: */
                GOOD[i][currIndex]=-2;
            }
        }
    }
}

/* XXX similar to unfold_adjacent but slighly different XXX */
static void unfold_adjacent2(
    float val, Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval, 
    int i, int currIndex, int countindex, int loopcount,
    int binindex[8], int rayindex[8]
    )
{
    int in, out, numneg, numpos, l, numtimes;
    float diff;
    numtimes=0;
    while(val!=missingVal && GOOD[i][currIndex]==0) {
        numtimes++;
        in = out = numneg = numpos = 0;
        for (l=0; l<countindex; l++) {
            diff = ray_val(rv_sweep->ray[rayindex[l]], binindex[l]) - val;
            if (fabs(diff) < THRESH*NyqVelocity) {
                in++;
            } else {
                out++;
            }
            if (diff>NyqVelocity) {
                numpos++;
            } else if (diff<-NyqVelocity) {
                numneg++;
            }
        }
        if (in>out) {
            ray_set(rv_sweep->ray[currIndex], i, val);
            GOOD[i][currIndex]=1;
        } else {
            if (numtimes<=MAXCOUNT) {
                if (numpos+numneg<(in+out-(numpos+numneg))) {
                    if (loopcount<=2) {
                        val=missingVal; /* Try later */
                    } else {
                        /* Keep the value after two passes through
                        ** data */
                        ray_set(rv_sweep->ray[currIndex], i, val);
                        GOOD[i][currIndex]=1;
                        
                    }
                } else if (numpos>numneg) {
                    val=val+NyqInterval;
                } else if (numneg>numpos) {
                    val=val-NyqInterval;
                } else {
                    /* Remove bin after four passes through data: */
                    if (loopcount>4) GOOD[i][currIndex]=-1;
                }
            } else {
                /* Remove bin: */
                GOOD[i][currIndex]=-1;
            }
        }
    }
}

/* Unfold GOOD=0 bins assuming spatial continuity: */
static void spatial_dealias(
    SweepCollection *sweepc, DealiasParams *dp, 
    short GOOD[MAXBINS][MAXRAYS],
    int *pstep, int pass)
{
    int loopcount, start, end, i, countindex, currIndex;
    int binindex[8], rayindex[8];
    float val;
    unsigned short flag;

    Sweep *rv_sweep;
    rv_sweep = sweepc->rv;

    loopcount=0;
    flag=1;
    while (flag==1) {
        loopcount++;
        flag=0;
        if (*pstep==1) {
            *pstep=-1;
            start=sweepc->nrays-1;
            end=-1;
        } else {
            *pstep=1;
            start=0;
            end=sweepc->nrays;
        }
        if (pass == 1) {
            for (i=DELNUM; i < sweepc->nbins; i++) {
                for (currIndex=start;currIndex!=end;currIndex=currIndex+*pstep) {
                    if (GOOD[i][currIndex] != 0)
                        continue;
                
                    countindex = count_neighbors(
                        GOOD, currIndex, i, sweepc->nrays, sweepc->nbins,
                        binindex, rayindex);
                
                    if (countindex>=1) {
                        flag=1;
                        val=ray_val(sweepc->vals->ray[currIndex], i); 
                        unfold_adjacent(
                            val, sweepc->rv, dp->missingVal, GOOD,
                            sweepc->NyqVelocity, sweepc->NyqInterval,
                            i, currIndex, countindex, loopcount,
                            binindex, rayindex);
                    }
                
                    /* This only needs to be performed on the first pass of the spatial
                    dealiasing */
                    if (loopcount==1 && countindex<1) {
                    mark_neighborless(GOOD, currIndex, i, 
                                      sweepc->nrays, sweepc->nbins);
                    }
                }
            }
        } else {  /* 2nd pass, reverse ray/bin loop order, different unfold 
                     and no marking of neighborless gates */
            for (currIndex=start;currIndex!=end;currIndex=currIndex+*pstep) {
                for (i=DELNUM; i < sweepc->nbins; i++) {
                    if (GOOD[i][currIndex] != 0)
                        continue;
                
                    countindex = count_neighbors(
                        GOOD, currIndex, i, sweepc->nrays, sweepc->nbins,
                        binindex, rayindex);
                
                    if (countindex>=1) {
                        flag=1;
                        val=ray_val(sweepc->vals->ray[currIndex], i); 
                        unfold_adjacent2(
                            val, sweepc->rv, dp->missingVal, GOOD,
                            sweepc->NyqVelocity, sweepc->NyqInterval,
                            i, currIndex, countindex, loopcount,
                            binindex, rayindex);
                    }
                }
            }
        }    
    }
}

/* Calculate the mean for gates within a window of a sweep */
static float window(Sweep *rv_sweep, int startray, 
	int endray, int firstbin, int lastbin, float missingVal,
    unsigned short* success) 
{
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

/* Detemine the window mean and if windowing is succeeded */
static float find_window(
    Sweep *rv_sweep, float missingVal, int currIndex, int i, unsigned short *wsuccess)
{
    int startray, endray, firstbin, lastbin;
    float winval;
    int numRays = rv_sweep->h.nrays;
    int numBins = rv_sweep->ray[0]->h.nbins;

    startray = currIndex - PROXIMITY;
    endray = currIndex + PROXIMITY;
    firstbin = i - PROXIMITY;
    lastbin = i + PROXIMITY;
    if (startray < 0) startray = numRays + startray;
    if (endray > numRays - 1) endray = endray - numRays;
    if (firstbin < 0) firstbin = 0;
    if (lastbin > numBins-1) lastbin = numBins - 1;
    winval=window(rv_sweep, startray, endray, 
                   firstbin, lastbin, missingVal, wsuccess);
    if (winval == missingVal && *wsuccess == 1) {     
        /* Expand the window: */
        startray = currIndex - 2 * PROXIMITY;
        endray = currIndex + 2 * PROXIMITY;
        firstbin = i - 2 * PROXIMITY;
        lastbin = i + 2 * PROXIMITY;
        if (startray < 0) startray = numRays + startray;
        if (endray > numRays - 1) endray = endray - numRays;
        if (firstbin < 0) firstbin = 0;
        if (lastbin > numBins - 1) lastbin = numBins - 1;
        winval=window(rv_sweep, startray, endray, 
                       firstbin, lastbin, missingVal, wsuccess);
    }
    return winval;
}

/* Unfold remote bins or those that were previously unsuccessful
 * using a window with dimensions 2(PROXIMITY)+1 x 2(PROXIMITY)+1:
 * if still no luck delete data (or unfold against VAD if PASS2). */
static void unfold_remote(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval)
{
    
    int i, j, currIndex; 
    unsigned short wsuccess;
    float val, diff, winval;
    int numRays = rv_sweep->h.nrays;
    int numBins = rv_sweep->ray[0]->h.nbins;

    for (i=DELNUM; i<numBins; i++) {
        for (currIndex=0; currIndex<numRays; currIndex++) { 
            
            if (GOOD[i][currIndex]!=0 && GOOD[i][currIndex]!=-2)
                continue;
            val = ray_val(vals_sweep->ray[currIndex], i); 
            winval = find_window(rv_sweep, missingVal, currIndex, i, &wsuccess);
            if (winval!=missingVal) {
                
                /* Add or subtract Nyquist intervals until within range */
                for (j=0; j<MAXCOUNT; j++) {
                    if (fabs(winval - val) <= 0.99999 * NyqVelocity)
                        break;
                    if (val > winval)
                        val -= NyqInterval;
                    else
                        val += NyqInterval;
                }
                diff = fabs(winval - val);
                
                if (diff < THRESH * NyqVelocity) {
                    /* Return the value. */
                    ray_set(rv_sweep->ray[currIndex], i, val);
                    GOOD[i][currIndex]=1;
                } else if (diff<(1.0 - (1.0 - THRESH)/2.0)*NyqVelocity) {
                    /* If within relaxed threshold, then return value, but
                    **   do not use to dealias other bins. */
                    ray_set(rv_sweep->ray[currIndex], i, val);
                    GOOD[i][currIndex]=-1;
                } else {
                    /* Remove bin */
                    GOOD[i][currIndex]=-1;
                }
            } else {
                if (wsuccess==0) {
                    /* Remove bin */
                    GOOD[i][currIndex]=-1; 
                    /* I don't think this ever happens -jjh */
                } else if (sound_sweep==NULL || last_sweep==NULL) {
                    if (GOOD[i][currIndex]==0 && RM!=1) {
                        /* Leave bin untouched. */
                        val = ray_val(vals_sweep->ray[currIndex], i);
                        ray_set(rv_sweep->ray[currIndex], i, val);
                        GOOD[i][currIndex]=-1; /* Don't use to unfold other bins*/
                    } else {
                        /* Remove bin */
                        GOOD[i][currIndex]=-1;
                    }
                } else if (GOOD[i][currIndex]==0 && PASS2 &&
                           sound_sweep!=NULL && last_sweep!=NULL) {
                    /* Leave GOOD[i][currIndex]=0 bins for a second pass.
                    ** In the second pass, we repeat unfolding, except this
                    ** time we use soundVolume for comparison instead of
                    ** lastVolume. */
                } else {
                    /* Remove bin */
                    GOOD[i][currIndex]=-1;
                }
            }
        }
    }
}

/* Second pass dealiasing using only sounding data */
static void second_pass(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval)
{

    int i, direction, currIndex; 
    unsigned short numtimes;
    float val, diff, valcheck, soundval;
    int numRays = rv_sweep->h.nrays;
    int numBins = rv_sweep->ray[0]->h.nbins;
    
     /* Beginning second pass, this time using sounding only: */
    for (currIndex=0;currIndex<numRays;currIndex++) {
        for (i=DELNUM;i<numBins;i++) {
            if (GOOD[i][currIndex]==0) {
                val = ray_val(vals_sweep->ray[currIndex], i);
                valcheck=val;
                soundval = ray_val(sound_sweep->ray[currIndex], i);
         
                if (soundval!=missingVal && val!=missingVal) {
                    diff=soundval-val;        
                    if (diff<0.0) {
                        diff=-diff;
                        direction=-1;
                    } else {
                        direction=1;
                    }
                    numtimes=0;
                    while (diff>0.99999*NyqVelocity && numtimes<=MAXCOUNT) {
                        val=val+NyqInterval*direction;
                        numtimes=numtimes+1;
                        diff=soundval-val;
                        if (diff<0.0) {
                            diff=-diff;
                            direction=-1;
                        } else {
                            direction=1;
                        }
                    }
                    if (diff<COMPTHRESH2*NyqVelocity&&fabs(valcheck)>CKVAL) {
                        /* Return the good value. */
                        ray_set(rv_sweep->ray[currIndex], i, val); 
		                GOOD[i][currIndex]=1;
                    }
                }
            }
        }
    }
}



/********************
 * Public functions *
*********************/

/*
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and the previous elevation sweep to unfold some of gates 
**  in each sweep. Then, it spreads outward from the 'good' gates, completing
**  the unfolding using gate-to-gate continuity. Gates that still remain
**  unfolded are compared to an areal average of neighboring dealiased gates.
**  Isolated echoes that still remain uncorrected are dealiased against a VAD
**  (as a last resort).
**
** This routine performs a preliminary unfold on the data using the bin in the
** next highest sweep (already unfolded), and the last unfolded volume. If the 
** last unfolded volume is not available, the VAD is used (or sounding if VAD
** not available). If this is successful, the bin is considered GOOD. Every 
** other bin is set to GOOD=0 or GOOD=-1 (if the bin is missing or bad).
**
** Then, the algorithm scans azimuthally and radially. If the majority of the
** GOOD=1 bins (up to 8) adjacent to a particular GOOD=0 bin are within 
** THRESH*Vnyq it is saved and considered GOOD as well. Otherwise, Nyquist 
** intervals are added/subtracted until the condition is met. If unsuccessful, 
** then GOOD is set to -2 for that bin.
**
** When all GOOD=0 bins that are next to GOOD=1 bins have been examined, then
** GOOD=0 and GOOD=-2 bins are unfolded against a window of GOOD=1 values. If
** there are too few GOOD=1 values inside the window, then the data are
** unfolded against the VAD, if available.
** 
** PARAMETERS:
** rvVolume: The radial velocity field to be unfolded.
** soundVolume: A first guess radial velocity field using VAD (or sounding) 
** lastVolume: The last radial velocity field unfolded (time of lastVolume 
**    should be within 10 minutes prior to rvVolume)
** missingVal: The value for missing radial velocity data.
** filt: A flag that specifies whether or not to use Bergen/Albers filter.
** success: flag indicating whether or not unfolding is possible.
*/
void unfoldVolume(Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
      float missingVal, unsigned short filt, unsigned short* success) 
{

    int sweepIndex, numSweeps;
    int step = -1;
    short GOOD[MAXBINS][MAXRAYS];
    float NyqVelocity, NyqInterval;
    
    Volume* VALS;
    Sweep *rv_sweep, *vals_sweep, *last_sweep, *sound_sweep, *above_sweep;
    SweepCollection sweepc;
    DealiasParams dp;

    /* Fill in Dealiasing parameters from arguments */
    dp.missingVal = missingVal;
    dp.compthresh = COMPTHRESH;     /* XXX make this a function argument */
    dp.ckval = CKVAL;               /* XXX Make this a function argument  */

    // Either a sounding or last volume must be provided
    if (soundVolume==NULL && lastVolume==NULL) {
        printf("First guess not available.\n");
        *success=0;
        return;
    }
    
    numSweeps = rvVolume->h.nsweeps;
    VALS=RSL_copy_volume(rvVolume);
    
    for (sweepIndex=numSweeps-1;sweepIndex>=0;sweepIndex--) {
        
        last_sweep = sound_sweep = above_sweep = NULL;
        if (lastVolume!=NULL) 
            last_sweep = lastVolume->sweep[sweepIndex];
        if (soundVolume!=NULL) 
            sound_sweep = soundVolume->sweep[sweepIndex];
        if (sweepIndex<numSweeps-1) 
            above_sweep = rvVolume->sweep[sweepIndex+1];
        rv_sweep = rvVolume->sweep[sweepIndex];
        vals_sweep = VALS->sweep[sweepIndex];

        NyqVelocity =rvVolume->sweep[sweepIndex]->ray[0]->h.nyq_vel;
        if (NyqVelocity == 0.0) NyqVelocity=9.8;
        NyqInterval = 2.0 * NyqVelocity;

        sweepc.vals = vals_sweep;
        sweepc.rv = rv_sweep;
        sweepc.last = last_sweep;
        sweepc.sound = sound_sweep;
        sweepc.above = above_sweep;
        sweepc.NyqVelocity = NyqVelocity;
        sweepc.NyqInterval = NyqInterval;
        sweepc.nrays = rv_sweep->h.nrays;
        sweepc.nbins = rv_sweep->ray[0]->h.nbins;

        /* Fill GOOD with -1 for missing value, 0 where not missing.
         * If requested the Berger Albers 3x3 filter is applied. */
        if (filt == 1)
            berger_albers_filter(sweepc.vals, dp.missingVal, GOOD);
        else
            zero_good(sweepc.vals, dp.missingVal, GOOD);
    
        /* Initialize the output sweep with missingVal */
        init_sweep(sweepc.rv, missingVal);

        /* Find non-missing value gates where the velocity either agrees with
         * or can be unfolded to agree with the sounding, the last volume, 
         * or the sweep directly above the current sweep.  Record the unfolded
         * velocity in these cases and mark GOOD with 1. */
        continuity_dealias(&sweepc, &dp, GOOD);
        //continuity_dealias(
        //    vals_sweep, rv_sweep, last_sweep, sound_sweep, above_sweep,
        //    missingVal, GOOD, NyqVelocity, NyqInterval);

        spatial_dealias(&sweepc, &dp, GOOD, &step, 1);

        unfold_remote(
            vals_sweep, rv_sweep, last_sweep, sound_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval);   

        if (last_sweep!=NULL && sound_sweep!=NULL) {
	   
            second_pass(
                vals_sweep, rv_sweep, last_sweep, sound_sweep,
                missingVal, GOOD, NyqVelocity, NyqInterval);
            
            spatial_dealias( 
                &sweepc, &dp, GOOD,
                &step, 2);
        }
    } // end of loop over sweeps
    *success=1;
    return;
}
