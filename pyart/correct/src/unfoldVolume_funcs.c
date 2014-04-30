

#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h>


/* Private function prototypes */
float sweep_val(Sweep *sweep, int ray_index, int range_index);

void ray_set(Ray *ray, int index, float val); 

/* Gate value getters */

float sweep_val(Sweep *sweep, int ray_index, int range_index)
{
    return ray_val(sweep->ray[ray_index], range_index);
}

float ray_val(Ray *ray, int index) 
{
    /* Return the value of range bin at index in ray */
    return (float) ray->h.f(ray->range[index]);
}

void ray_set(Ray *ray, int index, float val) 
{
    /* Set the gate value in ray at index to val */
    ray->range[index]=(unsigned short) ray->h.invf(val);
    return;
}


void bergen_albers_filter(Sweep *vals_sweep, int currIndex, int i,
                          float missingVal, short GOOD[MAXBINS][MAXRAYS]);

void bergen_albers_filter(Sweep *sweep, int currIndex, int i,
                         float missingVal, short GOOD[MAXBINS][MAXRAYS])
{
    /* Perform a 3x3 filter, as proposed by Bergen & Albers 1988 */
    int countindex=0;
    int left, right, next, prev;

    if (currIndex==0) left=sweep->h.nrays-1;
    else left=currIndex-1;
    if (currIndex==sweep->h.nrays-1) right=0;
    else right=currIndex+1;
    next=i+1;
    prev=i-1;
    /* Look at all bins adjacent to current bin in question: */
    if (i>DELNUM) {
        if (ray_val(sweep->ray[left], prev) != missingVal){
            countindex++;
        } 
        if (ray_val(sweep->ray[currIndex], prev) != missingVal) {
            countindex++;
        }
        if (ray_val(sweep->ray[right], prev) != missingVal) {
            countindex++;
        }
    }
    if (ray_val(sweep->ray[left], i) != missingVal) {
        countindex++;
    }
    if (ray_val(sweep->ray[right], i) != missingVal) {
        countindex++;
    }
    if (i<sweep->ray[0]->h.nbins-1) {  
        if (ray_val(sweep->ray[left], next) != missingVal) {
            countindex++;
        }
        if (ray_val(sweep->ray[currIndex], next)!= missingVal) {
            countindex++;
        }
        if (ray_val(sweep->ray[right], next) != missingVal) {
            countindex++;
        }
    }
    if (((i==sweep->ray[0]->h.nbins-1 || i==DELNUM) && countindex>=3) || (countindex>=5)) {
      GOOD[i][currIndex] = 0;  /* Save the bin for dealiasing. */
    } else {
      GOOD[i][currIndex] = -1; /* Assign missing value to the current bin. */
    }   
    return;
}   


void continuity_dealias(
    Sweep *rv_sweep, Sweep *sound_sweep, Sweep *last_sweep, Sweep *above_sweep,
    int currIndex, int i, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float val, int prevIndex, int abIndex, float NyqVelocity,
    float NyqInterval, float valcheck);

void continuity_dealias(
    Sweep *rv_sweep, Sweep *sound_sweep, Sweep *last_sweep, Sweep *above_sweep,
    int currIndex, int i, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float val, int prevIndex, int abIndex, float NyqVelocity,
    float NyqInterval, float valcheck)
{
    /* Try to dealias the bin using vertical and temporal 
     **   continuity (this is initial dealiasing). */
    
    int direction;
    unsigned short numtimes, dcase; 
    float prevval, soundval, abval, cval, diff;

    if (val==missingVal) return;  /* return if val is missingVal */
     
    /* Find velocities in related volumes */
    prevval = soundval = abval = missingVal;
    if (last_sweep!=NULL)
        prevval=ray_val(last_sweep->ray[prevIndex], i);
    if (sound_sweep!=NULL && last_sweep==NULL)
        soundval=ray_val(sound_sweep->ray[currIndex], i);
    if (above_sweep!=NULL)
        abval=ray_val(above_sweep->ray[abIndex], i);
    
    if (last_sweep==NULL && soundval!=missingVal && abval==missingVal) {
        cval=soundval;
        dcase=1;
    } else if (last_sweep==NULL && soundval!= missingVal && abval!=missingVal) {
        cval=abval;
        dcase=2;
    } else if (prevval!=missingVal && abval!=missingVal) {
        cval=prevval;
        dcase=3;
    } else {
        cval=missingVal;
        dcase=0;
    }
    
    if (dcase>0) {
        diff=cval-val;      
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
            diff=cval-val;
            if (diff<0.0) {
                diff=-diff;
                direction=-1;
            } else {
                direction=1;
            }
        }
        if (dcase==1 && diff<COMPTHRESH*NyqVelocity && fabs(valcheck)>CKVAL) {
            ray_set(rv_sweep->ray[currIndex], i, val);
            GOOD[i][currIndex]=1;
        } else if (dcase==2 && diff<COMPTHRESH*NyqVelocity && 
                   fabs(soundval-val)<COMPTHRESH*NyqVelocity && 
                   fabs(valcheck)>CKVAL) {
            ray_set(rv_sweep->ray[currIndex], i, val);
            GOOD[i][currIndex]=1;
        } else if (dcase==3 && diff<COMPTHRESH*NyqVelocity && fabs(abval-val)<
                   COMPTHRESH*NyqVelocity&&fabs(valcheck)>CKVAL) {
            ray_set(rv_sweep->ray[currIndex], i, val);
            GOOD[i][currIndex]=1;
        }
    }
}



/* Public functions */

void foobar(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    Sweep *above_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval,
    unsigned short filt
    )
{
    /** Unfold bins where vertical and temporal continuity
     ** produces the same value (i.e., bins where the radial velocity
     ** agrees with both the previous sweep and previous volume within
     ** a COMPTHRESH of the Nyquist velocity). This produces a number
     ** of good data points from which to begin unfolding assuming spatial
     ** continuity. */

    int currIndex;
    int i, prevIndex = 0, abIndex = 0; 
    float val, valcheck;

    for (currIndex=0;currIndex<rv_sweep->h.nrays;currIndex++) {
       
        if (last_sweep!=NULL)
           prevIndex=findRay2(rv_sweep, last_sweep, currIndex, missingVal);
        if (above_sweep!=NULL)
            abIndex=findRay2(rv_sweep, rv_sweep, currIndex, missingVal);
        for (i=DELNUM;i<rv_sweep->ray[0]->h.nbins;i++) {
            /* Assign uncorrect velocity bins to the array VALS: */
            ray_set(rv_sweep->ray[currIndex], i, missingVal);
            valcheck = val = ray_val(vals_sweep->ray[currIndex], i);
            if (val==missingVal) {
                GOOD[i][currIndex]=-1;
            } else {
                if (filt==1)
                    bergen_albers_filter(
                        vals_sweep, currIndex, i, missingVal, GOOD);
                else
                    /* If no filter is being applied save bin for dealiasing: */
                    GOOD[i][currIndex]=0;
            }
            if (GOOD[i][currIndex]==0)      
                continuity_dealias(
                    rv_sweep, sound_sweep, last_sweep, above_sweep,
                    currIndex, i, 
                    missingVal, GOOD,
                    val, prevIndex, abIndex, NyqVelocity,
                    NyqInterval, valcheck);
        }
    }
}

// XXX XXX XXX XXX
int count_neighbors(    
    Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    int i, int currIndex, int loopcount,
    int binindex[8], int rayindex[8]
    );

int count_neighbors(
    Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    int i, int currIndex, int loopcount,
    int binindex[8], int rayindex[8]
    )
{
    int next, prev, left, right;
    int countindex, countbins;
    countindex = countbins = 0;
    
    if (currIndex==0) left=rv_sweep->h.nrays-1;
    else left=currIndex-1;
    if (currIndex==rv_sweep->h.nrays-1) right=0;
    else right=currIndex+1;
    next=i+1;
    prev=i-1;
 
    /* Look at all bins adjacent to current bin in question: */
    if (i>DELNUM) {
        if (GOOD[prev][left]==1) {
            binindex[countindex]=prev;
            rayindex[countindex]=left;
            countindex++;
        }
        if (GOOD[prev][currIndex]==1) {
            binindex[countindex]=prev;
            rayindex[countindex]=currIndex;
            countindex++;
        }
        if (GOOD[prev][right]==1) {
            binindex[countindex]=prev;
            rayindex[countindex]=right;
            countindex++;
        }
    }
    if (GOOD[i][left]==1) {
        binindex[countindex]=i;
        rayindex[countindex]=left;
        countindex++;
    }
    if (GOOD[i][right]==1) {
        binindex[countindex]=i;
        rayindex[countindex]=right;
        countindex++;
    }
    if (i<rv_sweep->ray[0]->h.nbins-1) {  
        if (GOOD[next][left]==1) {
            binindex[countindex]=next;
            rayindex[countindex]=left;
            countindex++;
        }
        if (GOOD[next][currIndex]==1) {
            binindex[countindex]=next;
            rayindex[countindex]=currIndex;
            countindex++;
        }
        if (GOOD[next][right]==1) {
            binindex[countindex]=next;
            rayindex[countindex]=right;
            countindex++;
        }
    }

    /* This only needs to be performed on the first pass of the spatial
       dealiasing */
    if (loopcount==1 && countindex<1) { 
        if (i>DELNUM) {
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
        if (i<rv_sweep->ray[0]->h.nbins-1) {
            if (GOOD[next][left]==0)
                countbins++;
            if (GOOD[next][currIndex]==0)
                countbins++;
            if (GOOD[next][right]==0)
                countbins++;
        }
        if (countbins == 0)
        GOOD[i][currIndex]=-1;  
    }
    return countindex;
}

void unfold_adjacent(
    float val, Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval, 
    int i, int currIndex, int countindex, int loopcount,
    int binindex[8], int rayindex[8]
    );

void unfold_adjacent(
    float val, Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval, 
    int i, int currIndex, int countindex, int loopcount,
    int binindex[8], int rayindex[8]
    )
{
    /* Unfold against all adjacent values where GOOD==1 */
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




void spatial_dealias(
    Sweep *vals_sweep, Sweep *rv_sweep,
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval, 
    int *pstep)
{

    int loopcount, start, end, i, countindex, currIndex;
    int binindex[8], rayindex[8];
    float val;
    unsigned short flag;

    /* Now, unfold GOOD=0 bins assuming spatial continuity: */
    loopcount=0;
    flag=1;
    while (flag==1) {
        loopcount++;
        flag=0;
        if (*pstep==1) {
            *pstep=-1;
            start=rv_sweep->h.nrays-1;
            end=-1;
        } else {
            *pstep=1;
            start=0;
            end=rv_sweep->h.nrays;
        }
        for (i=DELNUM; i<rv_sweep->ray[0]->h.nbins; i++) {
            for (currIndex=start;currIndex!=end;currIndex=currIndex+*pstep) {
                if (GOOD[i][currIndex] != 0)
                    continue;
                val=ray_val(vals_sweep->ray[currIndex], i); 
                if (val == missingVal)
                    continue; 
                countindex = count_neighbors(
                    rv_sweep, missingVal, GOOD,
                    i, currIndex, loopcount,
                    binindex, rayindex);
                if (countindex>=1) {
                    flag=1;
                    unfold_adjacent(
                        val, rv_sweep, missingVal, GOOD,
                        NyqVelocity, NyqInterval,
                        i, currIndex, countindex, loopcount,
                        binindex, rayindex);
                }
            }
        }
    }

}

void second_pass(
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
                    }
                }
            }
        }
    }
}


float find_window(Sweep *rv_sweep, float missingVal, int currIndex, int i, unsigned short *wsuccess);

float find_window(Sweep *rv_sweep, float missingVal, int currIndex, int i, unsigned short *wsuccess)
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

void unfold_remote(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval)
{
    /* Unfold remote bins or those that were previously unsuccessful
     **   using a window with dimensions 2(PROXIMITY)+1 x 2(PROXIMITY)+1:
     **   if still no luck delete data (or unfold against VAD if PASS2). */
    
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
                
                if (diff < THRESH  *NyqVelocity) {
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
