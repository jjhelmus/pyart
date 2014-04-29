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

    int sweepIndex, numSweeps;
    int step = -1;
    
    unsigned short flag=1;
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

        flag=1;

        foobar(
            vals_sweep, rv_sweep, last_sweep, sound_sweep, above_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval,
            filt);

        spatial_dealias( 
            vals_sweep, rv_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval, 
            &flag, &step);

        unfold_remote(
            vals_sweep, rv_sweep, last_sweep, sound_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval);   

        if (last_sweep!=NULL && sound_sweep!=NULL) {
            flag=1;
            second_pass(
                vals_sweep, rv_sweep, last_sweep, sound_sweep,
                missingVal, GOOD, NyqVelocity, NyqInterval);
        }
 
        spatial_dealias( 
            vals_sweep, rv_sweep,
            missingVal, GOOD, NyqVelocity, NyqInterval, 
            &flag, &step);
    
    } // end of loop over sweeps
    *success=1;
    return;
}
