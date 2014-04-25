

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
