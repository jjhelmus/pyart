/*
**
**  UW Radial Velocity Unfolding Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and a higher elevation sweep to unfold some of the bins
**  in each sweep. Then, it completes the unfolding, assuming spatial
**  continuity around each bin.
**
**  DEVELOPER:
**     Curtis N. James            25 Jan 99
**
*/

#include "helpers.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */ 

/* Private function prototypes */
void foobar(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    Sweep *above_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval,
    unsigned short filt);

    int findRay(Sweep *sweep1, Sweep *sweep2, int currIndex1, float missingVal);
    
    void bergen_albers_filter(Sweep *vals_sweep, int currIndex, int i,
                          float missingVal, short GOOD[MAXBINS][MAXRAYS]);

    void continuity_dealias(
        Sweep *rv_sweep, Sweep *sound_sweep, Sweep *last_sweep, Sweep *above_sweep,
        int currIndex, int i, 
        float missingVal, short GOOD[MAXBINS][MAXRAYS],
        float val, int prevIndex, int abIndex, float NyqVelocity,
        float NyqInterval, float valcheck);

void spatial_dealias(
    Sweep *vals_sweep, Sweep *rv_sweep,
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval,
    int *step);

    int count_neighbors(    
        Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
        int i, int currIndex, int loopcount,
        int binindex[8], int rayindex[8]);

    void unfold_adjacent(
        float val, Sweep* rv_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
        float NyqVelocity, float NyqInterval, 
        int i, int currIndex, int countindex, int loopcount,
        int binindex[8], int rayindex[8]);


void unfold_remote(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval);

    float find_window(Sweep *rv_sweep, float missingVal, int currIndex, 
                      int i, unsigned short *wsuccess);

    float window(Sweep *rv_sweep, int startray, int endray,
                 int firstbin, int lastbin, float missingVal, unsigned 
                 short* success);

void second_pass(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval);

/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**     This algorithm unfolds a volume of single Doppler radial velocity data.
**  The algorithm uses a previously unfolded volume (or VAD if previous volume
**  is unavailable) and the previous elevation sweep to unfold some of gates 
**  in each sweep. Then, it spreads outward from the 'good' gates, completing
**  the unfolding using gate-to-gate continuity. Gates that still remain
**  unfolded are compared to an areal average of neighboring dealiased gates.
**  Isolated echoes that still remain uncorrected are dealiased against a VAD
**  (as a last resort).
**
**  DEVELOPER:
**  Curtis N. James     25 Jan 99
**
**  Refactor April, 2014 by Jonathan J. Helmus (jhelmus@anl.gov)
**
**
*/


/* This routine performs a preliminary unfold on the data using the bin in the
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

        foobar(
            vals_sweep, rv_sweep, last_sweep, sound_sweep, above_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval,
            filt);

        spatial_dealias( 
            vals_sweep, rv_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval, 
            &step);

        unfold_remote(
            vals_sweep, rv_sweep, last_sweep, sound_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval);   

        if (last_sweep!=NULL && sound_sweep!=NULL) {
            second_pass(
                vals_sweep, rv_sweep, last_sweep, sound_sweep,
                missingVal, GOOD, NyqVelocity, NyqInterval);
        
            spatial_dealias( 
                vals_sweep, rv_sweep,
                missingVal, GOOD, NyqVelocity, NyqInterval, 
                &step);
        }
    } // end of loop over sweeps
    *success=1;
    return;
}


/* Private functions */

/** Unfold bins where vertical and temporal continuity
 ** produces the same value (i.e., bins where the radial velocity
 ** agrees with both the previous sweep and previous volume within
 ** a COMPTHRESH of the Nyquist velocity). This produces a number
 ** of good data points from which to begin unfolding assuming spatial
 ** continuity. */
void foobar(
    Sweep *vals_sweep, Sweep *rv_sweep, Sweep *last_sweep, Sweep *sound_sweep,
    Sweep *above_sweep, float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval,
    unsigned short filt
    )
{

    int currIndex;
    int i, prevIndex = 0, abIndex = 0; 
    float val, valcheck;

    for (currIndex=0;currIndex<rv_sweep->h.nrays;currIndex++) {
       
        if (last_sweep!=NULL)
           prevIndex=findRay(rv_sweep, last_sweep, currIndex, missingVal);
        if (above_sweep!=NULL)
            abIndex=findRay(rv_sweep, rv_sweep, currIndex, missingVal);
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

/*
**
**  UW Radial Velocity Dealiasing Algorithm
**  Four-Dimensional Dealiasing (4DD)
**
**  DESCRIPTION:
**      This routine finds the rayindex of the nearest ray in sweepIndex2 of 
**  rvVolume2 to sweepIndex1 in rvVolume1.
**
**  DEVELOPER:
**	Curtis N. James    1 Feb 1999
**
**
**
*/

int findRay(Sweep *sweep1, Sweep *sweep2, int currIndex, float missingVal) 
{

     int numRays, rayIndex1;
     float az0, az1, diffaz;
     float spacing;
     short direction, lastdir;
    
     numRays = sweep2->h.nrays;

     az0 = sweep1->ray[currIndex]->h.azimuth;
     if (currIndex<numRays) rayIndex1=currIndex;
     else rayIndex1=numRays-1;
     az1 = sweep2->ray[rayIndex1]->h.azimuth;
     if (az0==az1) {
       return rayIndex1;
     } else {
       /* Since the beamwidth is not necessarily the spacing between rays: */
       spacing = fabs(sweep2->ray[0]->h.azimuth-
		              sweep2->ray[50]->h.azimuth); 
       if (spacing>180) spacing=360.0-spacing;
       spacing=spacing/50.0;

       /* Compute the difference in azimuth between the two rays: */
       diffaz=az0-az1;
       if (diffaz>=180.0) diffaz=diffaz-360.0;
       else if (diffaz<-180.0) diffaz=diffaz+360.0;
       
       /* Get close to the correct index: */
       rayIndex1=rayIndex1+(int) (diffaz/spacing);
       if (rayIndex1>=numRays) rayIndex1=rayIndex1-numRays;
       if (rayIndex1<0) rayIndex1=numRays+rayIndex1;
       az1=sweep2->ray[rayIndex1]->h.azimuth;
       diffaz=az0-az1;
       if (diffaz>=180.0) diffaz=diffaz-360.0;
       else if (diffaz<-180.0) diffaz=diffaz+360.0;

       /* Now add or subtract indices until the nearest ray is found: */
       if (diffaz>=0) lastdir=1;
       else lastdir=-1;
       while (fabs(diffaz)>spacing/2.0) {
	 if (diffaz>=0) {
	   rayIndex1++;
	   direction=1;
	 } else {
	   rayIndex1--;
	   direction=-1;
	 }
	 if (rayIndex1>=numRays) rayIndex1=rayIndex1-numRays;
	 if (rayIndex1<0) rayIndex1=numRays+rayIndex1;
	 az1=sweep2->ray[rayIndex1]->h.azimuth;
	 diffaz=az0-az1;
	 if (diffaz>=180.0) diffaz=diffaz-360.0;
	 else if (diffaz<-180.0) diffaz=diffaz+360.0;
	 if (direction!=lastdir) break;
	 else lastdir=direction;
       }
       return rayIndex1;
     }
}
/* Perform a 3x3 filter, as proposed by Bergen & Albers 1988 */
void bergen_albers_filter(Sweep *sweep, int currIndex, int i,
                         float missingVal, short GOOD[MAXBINS][MAXRAYS])
{
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

/* Try to dealias the bin using vertical and temporal 
 **   continuity (this is initial dealiasing). */
void continuity_dealias(
    Sweep *rv_sweep, Sweep *sound_sweep, Sweep *last_sweep, Sweep *above_sweep,
    int currIndex, int i, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float val, int prevIndex, int abIndex, float NyqVelocity,
    float NyqInterval, float valcheck)
{
    
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


/* Now, unfold GOOD=0 bins assuming spatial continuity: */
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

/* Unfold remote bins or those that were previously unsuccessful
 **   using a window with dimensions 2(PROXIMITY)+1 x 2(PROXIMITY)+1:
 **   if still no luck delete data (or unfold against VAD if PASS2). */
void unfold_remote(
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


float find_window(
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