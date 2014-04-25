

#include "FourDD.h"
#include <rsl.h>
#include <stdio.h>

/* Perform a 3x3 filter, as proposed by Bergen & Albers 1988 
 * GOOD is updated 
 * */


void bergen_albers_filter(Volume* VALS, int sweepIndex, int currIndex, int i,
                         int numRays, int numBins, float missingVal, 
                         short GOOD[MAXBINS][MAXRAYS])
{
    int countindex=0;
    int left, right, next, prev;

    countindex=0;
    if (currIndex==0) left=numRays-1;
    else left=currIndex-1;
    if (currIndex==numRays-1) right=0;
    else right=currIndex+1;
    next=i+1;
    prev=i-1;
    /* Look at all bins adjacent to current bin in question: */
    if (i>DELNUM) {
      if (((float) VALS->sweep[sweepIndex]->
       ray[left]->h.f(VALS->sweep[sweepIndex]->ray
               [left]->range[prev]))!=missingVal) {
        countindex=countindex+1;
      }
      if (((float) VALS->sweep[sweepIndex]->
       ray[currIndex]->h.f(VALS->sweep[sweepIndex]->ray
           [currIndex]->range[prev]))!=missingVal) {
        countindex=countindex+1;
      }
      if (((float) VALS->sweep[sweepIndex]->
       ray[right]->h.f(VALS->sweep[sweepIndex]->ray
       [right]->range[prev]))!=missingVal) {
        countindex=countindex+1;
      }
    }
    if (((float) VALS->sweep[sweepIndex]->
         ray[left]->h.f(VALS->sweep[sweepIndex]->ray
         [left]->range[i]))!=missingVal) {
      countindex=countindex+1;
    }
    if (((float) VALS->sweep[sweepIndex]->
         ray[right]->h.f(VALS->sweep[sweepIndex]->ray
         [right]->range[i]))!=missingVal) {
      countindex=countindex+1;
    }
    if (i<numBins-1) {  
      if (((float) VALS->sweep[sweepIndex]->
       ray[left]->h.f(VALS->sweep[sweepIndex]->ray
                   [left]->range[next]))!=missingVal) {
        countindex=countindex+1;
      }
      if (((float) VALS->sweep[sweepIndex]->
       ray[currIndex]->h.f(VALS->sweep[sweepIndex]->ray
       [currIndex]->range[next]))!=missingVal) {
        countindex=countindex+1;
      }
      if (((float) VALS->sweep[sweepIndex]->
       ray[right]->h.f(VALS->sweep[sweepIndex]->ray
                   [right]->range[next]))!=missingVal) {
        countindex=countindex+1;
      }
    }
    if (((i==numBins-1 || i==DELNUM) && countindex>=3)||
        (countindex>=5)) {
      /* Save the bin for dealiasing: */
      GOOD[i][currIndex] = 0;
      return;
    } else {
      /* Assign missing value to the current bin. */
      GOOD[i][currIndex] = -1;
      return;
    }   
      /* End 3 x 3 filter */
}   

void continuity_dealias(
    Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
    int sweepIndex, int currIndex, int i, int numRays, int numBins, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float val, int prevIndex, int numSweeps, int abIndex, float NyqVelocity,
    float NyqInterval, float valcheck, float fraction)
{
    int direction;
    unsigned short numtimes, dcase; 
    float prevval, soundval, abval, cval, finalval, diff;


           /* Try to dealias the bin using vertical and temporal 
           **   continuity (this is initial dealiasing). */
           if (val!=missingVal && lastVolume!=NULL) {
         prevval=(float) lastVolume->sweep[sweepIndex]->ray[prevIndex]
           ->h.f(lastVolume->sweep[sweepIndex]->ray[prevIndex]->
           range[i]);
           } else {
         prevval=missingVal;
           }
           if (val!=missingVal && soundVolume!=NULL && lastVolume==NULL) {
         soundval=(float) soundVolume->
           sweep[sweepIndex]->ray[currIndex]->h.f(soundVolume->
           sweep[sweepIndex]->ray[currIndex]->range[i]);
           } else {
         soundval=missingVal;
           }
           if (val!=missingVal && sweepIndex<numSweeps-1) {
         abval=(float) rvVolume->sweep[sweepIndex+1]->ray[abIndex]->
           h.f(rvVolume->sweep[sweepIndex+1]->ray[abIndex]->range[i]);
           } else {
         abval=missingVal;
           }

           if (val!=missingVal&&lastVolume==NULL&&soundval!=
           missingVal&&abval==missingVal) {
         cval=soundval;
         dcase=1;
           }
           else if (val!=missingVal&&lastVolume==NULL&&soundval!=
            missingVal&&abval!=missingVal) {
         cval=abval;
         dcase=2;
           }
           else if (val!=missingVal&&prevval!=missingVal&&abval!=
            missingVal) {
         cval=prevval;
         dcase=3;
           } 
           else {
         cval=missingVal;
         dcase=0;
           }
           if (dcase>0) {
         diff=cval-val;      
         if (diff<0.0) {
           diff=-diff;
           direction=-1;
         } else direction=1;
         numtimes=0;
         while (diff>0.99999*NyqVelocity && numtimes<=MAXCOUNT) {
           val=val+NyqInterval*direction;
           numtimes=numtimes+1;
           diff=cval-val;
           if (diff<0.0) {
             diff=-diff;
             direction=-1;
           } else direction=1;
         }
         if (dcase==1&&diff<fraction*NyqVelocity&&
             fabs(valcheck)>CKVAL) {
           if (VERBOSE) printf("GOOD1: %f\n", val);
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
           GOOD[i][currIndex]=1;
         }
         else if (dcase==2&&diff<fraction*NyqVelocity&&fabs(soundval
           -val)<fraction*NyqVelocity&&fabs(valcheck)>CKVAL) {
           if (VERBOSE) printf("GOOD2: %f\n", val);
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
           GOOD[i][currIndex]=1;
         }
         else if (dcase==3&&diff<fraction*NyqVelocity&&fabs(abval-val)<
              fraction*NyqVelocity&&fabs(valcheck)>CKVAL) {
           if (VERBOSE) printf("GOOD3: %f\n",val);
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval); 
           GOOD[i][currIndex]=1;
         }
           }
}
