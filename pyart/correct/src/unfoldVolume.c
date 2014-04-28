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
#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h> /* Sweep */

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

    int sweepIndex, currIndex, i, direction, numSweeps, numRays,
        numBins, rayindex[8], binindex[8],
        startray, endray, firstbin, 
        lastbin, step = -1, prevIndex, abIndex;
    
    unsigned short numtimes, flag=1, wsuccess;
    short GOOD[MAXBINS][MAXRAYS];
    float NyqVelocity, NyqInterval, val, diff, fraction, finalval, initval, 
        valcheck, winval, diffs[8], fraction2;
    float pfraction, std;
    
    Volume* VALS;

    abIndex = 0;        /* initialize to supress compiler warning */
    prevIndex = 0;      /* initialize to spresss compiler watning */

    // Either a sounding or last volume must be provided
    if (soundVolume==NULL && lastVolume==NULL) {
        printf("First guess not available.\n");
        *success=0;
        return;
    }
    
    numSweeps = rvVolume->h.nsweeps;
    VALS=RSL_copy_volume(rvVolume);
    if (COMPTHRESH>1.0 || COMPTHRESH<=0.0 ) fraction=0.25;
    else fraction=COMPTHRESH;
    if (COMPTHRESH2>1.0 || COMPTHRESH2<=0.0 ) fraction2=0.25;
    else fraction2=COMPTHRESH2;
    if (THRESH>1.0 || THRESH<=0.0 ) pfraction=0.5;
    else pfraction=THRESH;

    for (sweepIndex=numSweeps-1;sweepIndex>=0;sweepIndex--) {
        NyqVelocity =rvVolume->sweep[sweepIndex]->ray[0]->h.nyq_vel;
     if (NyqVelocity == 0.0) NyqVelocity=9.8;
         NyqInterval = 2.0 * NyqVelocity;
         printf("NyqVelocity %f \n", NyqVelocity);
     numRays = rvVolume->sweep[sweepIndex]->h.nrays;
     numBins = rvVolume->sweep[sweepIndex]->ray[0]->h.nbins;

     /* First, unfold bins where vertical and temporal continuity
     ** produces the same value (i.e., bins where the radial velocity
     ** agrees with both the previous sweep and previous volume within
     ** a COMPTHRESH of the Nyquist velocity). This produces a number
     ** of good data points from which to begin unfolding assuming spatial
         ** continuity. */
     printf("Sweep: %d\n", sweepIndex);
     if (VERBOSE) printf("NyqVelocity: %f, missingVal: %f\n",
                 NyqVelocity, missingVal);
     flag=1;

     for (currIndex=0;currIndex<numRays;currIndex++) {
       if (lastVolume!=NULL) prevIndex=findRay(rvVolume, lastVolume,
          sweepIndex, sweepIndex,currIndex, missingVal);
       if (sweepIndex<numSweeps-1) abIndex=findRay(rvVolume, rvVolume,
          sweepIndex, sweepIndex+1, currIndex, missingVal);
       for (i=DELNUM;i<numBins;i++) {
         /* Initialize Output Sweep with missing values: */
         initval=(float)rvVolume->sweep[sweepIndex]->
           ray[currIndex]->h.invf(missingVal);
         rvVolume->sweep[sweepIndex]->ray[currIndex]->
           range[i]=(unsigned short) (initval);

         /* Assign uncorrect velocity bins to the array VALS: */
         val=(float) VALS->sweep[sweepIndex]->
                 ray[currIndex]->h.f(VALS->sweep[sweepIndex]->ray
                 [currIndex]->range[i]);
         valcheck=val;
         if (val==missingVal) GOOD[i][currIndex]=-1;
         else {
           if (filt==1) {
               bergen_albers_filter(VALS, sweepIndex, currIndex, i, numRays,
                                    numBins, missingVal, GOOD);
           } else {
         /* If no filter is being applied save bin for dealiasing: */
         GOOD[i][currIndex]=0;
           }
         }
         
         if (GOOD[i][currIndex]==0) {
            
            continuity_dealias(rvVolume, soundVolume, lastVolume,
                sweepIndex, currIndex, i, numRays, numBins, 
                missingVal, GOOD,
                val, prevIndex, numSweeps, abIndex, NyqVelocity,
                NyqInterval, valcheck, fraction);
         }
       }
     }

     spatial_dealias( 
        VALS, rvVolume,
        sweepIndex, currIndex, numRays, numBins,
        missingVal, GOOD,
        NyqVelocity, NyqInterval, pfraction,
        &flag, &step, binindex, rayindex, diffs
             );

     /* Unfold remote bins or those that were previously unsuccessful
     **   using a window with dimensions 2(PROXIMITY)+1 x 2(PROXIMITY)+1:
     **   if still no luck delete data (or unfold against VAD if PASS2). */
     for (i=DELNUM;i<numBins;i++) {
       for (currIndex=0;currIndex<numRays;currIndex++) { 
         if (GOOD[i][currIndex]==0 || GOOD[i][currIndex]==-2) {
           val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
           h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range[i]);
           startray=currIndex-PROXIMITY;
           endray=currIndex+PROXIMITY;
           firstbin=i-PROXIMITY;
           lastbin=i+PROXIMITY;
           if (startray<0) startray=numRays+startray;
           if (endray>numRays-1) endray=endray-numRays;
           if (firstbin<0) firstbin=0;
           if (lastbin>numBins-1) lastbin=numBins-1;
           winval=window(rvVolume, sweepIndex, startray, endray, 
                 firstbin, lastbin, &std, missingVal, &wsuccess);
           if (winval==missingVal && wsuccess==1){ /* Expand the window: */
         startray=currIndex-2*PROXIMITY;
         endray=currIndex+2*PROXIMITY;
         firstbin=i-2*PROXIMITY;
         lastbin=i+2*PROXIMITY;
         if (startray<0) startray=numRays+startray;
         if (endray>numRays-1) endray=endray-numRays;
         if (firstbin<0) firstbin=0;
         if (lastbin>numBins-1) lastbin=numBins-1;
         winval=window(rvVolume, sweepIndex, startray, endray, 
                   firstbin, lastbin, &std, missingVal, &wsuccess);
           }
           if (winval!=missingVal) {
         diff=winval-val;
         if (diff<0.0) {
           diff=-diff;
           direction=-1;
         } else direction=1;
         numtimes=0;
         while (diff>0.99999*NyqVelocity && numtimes<=MAXCOUNT) {
           val=val+NyqInterval*direction;
           numtimes=numtimes+1;
           diff=winval-val;
           if (diff<0.0) {
             diff=-diff;
             direction=-1;
           } else direction=1;
         }
         if (diff<pfraction*NyqVelocity) {
           /* Return the value. */
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
           GOOD[i][currIndex]=1;
         } else if (diff<(1.0 - (1.0 - pfraction)/2.0)*NyqVelocity) {
           /* If within relaxed threshold, then return value, but
           **   do not use to dealias other bins. */
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
           GOOD[i][currIndex]=-1;          
         } else {
           /* Remove bin */
           GOOD[i][currIndex]=-1;
         }
           } else { 
         if (wsuccess==0) {
           /* Remove bin */
           GOOD[i][currIndex]=-1;
         } else if (soundVolume==NULL || lastVolume==NULL) {
           if (GOOD[i][currIndex]==0 && RM!=1) {
             /* Leave bin untouched. */
             val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
               h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range[i]);
             finalval=(float)rvVolume->sweep[sweepIndex]->
               ray[currIndex]->h.invf(val);
             rvVolume->sweep[sweepIndex]->ray[currIndex]->
               range[i]=(unsigned short) (finalval);
             GOOD[i][currIndex]=-1; /* Don't use to unfold other bins*/
           } else {
             /* Remove bin */
             GOOD[i][currIndex]=-1;
           }
         } else if (GOOD[i][currIndex]==0 && PASS2 &&
             soundVolume!=NULL && lastVolume!=NULL) {
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
     
     second_pass(
        VALS, rvVolume, soundVolume, lastVolume,
        sweepIndex, currIndex, numRays, numBins, 
        missingVal, GOOD,
        NyqVelocity, NyqInterval, pfraction, 
        numSweeps, fraction2, 
        &flag, &step,
        binindex, rayindex, diffs
             );
       } // end of loop over sweeps
    *success=1;
    return;
}
