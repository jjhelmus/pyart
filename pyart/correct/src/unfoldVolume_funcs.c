

#include "FourDD.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rsl.h>

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

void foobar(
    Volume* VALS, Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
    int sweepIndex, int numRays, int numBins, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval,
    int numSweeps,
    unsigned short filt
    )
{
    /* First, unfold bins where vertical and temporal continuity
     ** produces the same value (i.e., bins where the radial velocity
     ** agrees with both the previous sweep and previous volume within
     ** a COMPTHRESH of the Nyquist velocity). This produces a number
     ** of good data points from which to begin unfolding assuming spatial
         ** continuity. */


    int currIndex;
    int i, prevIndex, abIndex; 
    float val, initval, valcheck;

    Sweep *rvsweep = rvVolume->sweep[sweepIndex];
    Sweep *valssweep = VALS->sweep[sweepIndex];
    Ray *rvray;
    Ray *valsray;

    abIndex = 0;        /* initialize to supress compiler warning */
    prevIndex = 0;      /* initialize to spresss compiler watning */

    for (currIndex=0;currIndex<numRays;currIndex++) {
       
        rvray = rvsweep->ray[currIndex];
        valsray = valssweep->ray[currIndex];
         
        if (lastVolume!=NULL) {
           prevIndex=findRay(rvVolume, lastVolume, sweepIndex, sweepIndex,
                             currIndex, missingVal);
        }
        if (sweepIndex<numSweeps-1) {
            abIndex=findRay(rvVolume, rvVolume, sweepIndex, sweepIndex+1,
                            currIndex, missingVal);
        }
        for (i=DELNUM;i<numBins;i++) {
         
            initval=(float)rvray->h.invf(missingVal);
            rvray->range[i]=(unsigned short) (initval);

            /* Assign uncorrect velocity bins to the array VALS: */
            val=(float) valsray->h.f(valsray->range[i]);
            valcheck=val;
            if (val==missingVal) {
                GOOD[i][currIndex]=-1;
            } else {
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
                    NyqInterval, valcheck);
            }
        }
    }
}


void continuity_dealias(
    Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
    int sweepIndex, int currIndex, int i, int numRays, int numBins, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float val, int prevIndex, int numSweeps, int abIndex, float NyqVelocity,
    float NyqInterval, float valcheck)
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
         if (dcase==1&&diff<COMPTHRESH*NyqVelocity&&
             fabs(valcheck)>CKVAL) {
           if (VERBOSE) printf("GOOD1: %f\n", val);
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
           GOOD[i][currIndex]=1;
         }
         else if (dcase==2&&diff<COMPTHRESH*NyqVelocity&&fabs(soundval
           -val)<COMPTHRESH*NyqVelocity&&fabs(valcheck)>CKVAL) {
           if (VERBOSE) printf("GOOD2: %f\n", val);
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
           GOOD[i][currIndex]=1;
         }
         else if (dcase==3&&diff<COMPTHRESH*NyqVelocity&&fabs(abval-val)<
              COMPTHRESH*NyqVelocity&&fabs(valcheck)>CKVAL) {
           if (VERBOSE) printf("GOOD3: %f\n",val);
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval); 
           GOOD[i][currIndex]=1;
         }
           }
}


void spatial_dealias(
    Volume* VALS, Volume* rvVolume,
    int sweepIndex, int numRays, int numBins, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval, 
    unsigned short *pflag, int *pstep,
    int binindex[8], int rayindex[8], float diffs[8]
    )
{

    int loopcount, startindex, endindex, i, countindex, countbins,
        next, prev, left, right, numneg, numpos, l, in, out, numtimes;
    int currIndex;
    //unsigned short flag = *pflag;
    //int step = *pstep;
    float val, finalval, goodval;

     /* Now, unfold GOOD=0 bins assuming spatial continuity: */
     loopcount=0;
     while (*pflag==1) {
       loopcount=loopcount+1;
       *pflag=0;
       if (*pstep==1) {
         *pstep=-1;
         startindex=numRays-1;
         endindex=-1;
       } else {
         *pstep=1;
         startindex=0;
         endindex=numRays;
       }
       for (i=DELNUM;i<numBins;i++) {
         for (currIndex=startindex;currIndex!=endindex;currIndex=currIndex
        +*pstep) {
           if (GOOD[i][currIndex]==0) {
         countindex=0;
         countbins=0;
         val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
           h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range[i]);
         if (VERBOSE) printf("\n");
             if (VERBOSE) printf("Startval: %f\n", val);
             if (currIndex==0) left=numRays-1;
         else left=currIndex-1;
         if (currIndex==numRays-1) right=0;
         else right=currIndex+1;
         next=i+1;
         prev=i-1;
         /* Look at all bins adjacent to current bin in question: */
         if (val != missingVal) {  // don't look for neighbors of missingVal
         if (i>DELNUM) {
           if (GOOD[prev][left]==1) {
             if (*pflag==0) *pflag=1;
             countindex=countindex+1;
             binindex[countindex-1]=prev;
             rayindex[countindex-1]=left;
             if (VERBOSE) printf("pl ");
           }
           if (GOOD[prev][left]==0) {
             countbins=countbins+1;
           }
           if (GOOD[prev][currIndex]==1) {
             if (*pflag==0) *pflag=1;
             countindex=countindex+1;
             binindex[countindex-1]=prev;
             rayindex[countindex-1]=currIndex;
             if (VERBOSE) printf("pc ");
           }
           if (GOOD[prev][currIndex]==0) {
             countbins=countbins+1;
           }
           if (GOOD[prev][right]==1) {
             if (*pflag==0) *pflag=1;
             countindex=countindex+1;
             binindex[countindex-1]=prev;
             rayindex[countindex-1]=right;
             if (VERBOSE) printf("pr ");
           }
               if (GOOD[prev][right]==0) {
             countbins=countbins+1;
           }
         }
         if (GOOD[i][left]==1) {
           if (*pflag==0) *pflag=1;
           countindex=countindex+1;
           binindex[countindex-1]=i;
           rayindex[countindex-1]=left;
           if (VERBOSE) printf("il ");
         }
         if (GOOD[i][left]==0) {
           countbins=countbins+1;
         }
         if (GOOD[i][right]==1) {
           if (*pflag==0) *pflag=1;
           countindex=countindex+1;
           binindex[countindex-1]=i;
           rayindex[countindex-1]=right;
           if (VERBOSE) printf("ir ");     
         }
         if (GOOD[i][right]==0) {
           countbins=countbins+1;
         }
         if (i<numBins-1) {  
           if (GOOD[next][left]==1) {
             if (*pflag==0) *pflag=1;
             countindex=countindex+1;
             binindex[countindex-1]=next;
             rayindex[countindex-1]=left;
             if (VERBOSE) printf("nl "); 
           }
           if (GOOD[next][left]==0) {
             countbins=countbins+1;
           }
           if (GOOD[next][currIndex]==1) {
             if (*pflag==0) *pflag=1;
             countindex=countindex+1;
             binindex[countindex-1]=next;
             rayindex[countindex-1]=currIndex;
             if (VERBOSE) printf("nc ");
           }
           if (GOOD[next][currIndex]==0) {
             countbins=countbins+1;
           }
           if (GOOD[next][right]==1) {
             if (*pflag==0) *pflag=1;
             countindex=countindex+1;
             binindex[countindex-1]=next;
             rayindex[countindex-1]=right;
             if (VERBOSE) printf("nr ");
           }
           if (GOOD[next][right]==0) {
             countbins=countbins+1;
           }
         }
         }
         /* Perform last step of Bergen and Albers filter: */
         if (loopcount==1 && countbins+countindex<1) 
               GOOD[i][currIndex]=-1;  

         if (countindex>=1) {
           /* Unfold against all adjacent values where GOOD==1 */
           numtimes=0;
           while(val!=missingVal&&GOOD[i][currIndex]==0) {
             numtimes=numtimes+1;
             in=0;
             out=0;
             numneg=0;
             numpos=0;
             if (VERBOSE) printf("countindex: %d\n", countindex); 
             for (l=0;l<countindex;l++) {
               if (VERBOSE) printf("%d %d %f\n",binindex[l],
              rayindex[l],goodval);
               goodval=(float)rvVolume->sweep[sweepIndex]->ray
             [rayindex[l]]->h.f(rvVolume->sweep[sweepIndex]->
             ray[rayindex[l]]->range[binindex[l]]);
                   diffs[l]=goodval-val;
               if (fabs(diffs[l])<THRESH*NyqVelocity) in=in+1;
               else {
             out=out+1;
             if (diffs[l]>NyqVelocity) {
               numpos=numpos+1;
             } else if (diffs[l]<-NyqVelocity){
               numneg=numneg+1;
             }
               }
             }
             if (in>0 && out==0) {
               finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
               rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
               GOOD[i][currIndex]=1;
               if (VERBOSE) printf("Finalval: %4.2f\n", val);
             } else {
               if (numtimes<=MAXCOUNT) {
             if ((numpos+numneg)<(in+out-(numpos+numneg))) {
               if (loopcount>2)
               { /* Keep the value after two passes through
                 ** data. */
                 finalval=(float)rvVolume->sweep[sweepIndex]->
                   ray[currIndex]->h.invf(val);
                 rvVolume->sweep[sweepIndex]->ray[currIndex]->
                   range[i]=(unsigned short) (finalval);
                 GOOD[i][currIndex]=1;
                 if (VERBOSE) printf(
                     "Keeping after two passes: %f\n",val);
               }
             } else if (numpos>numneg) {
               val=val+NyqInterval;
             } else if (numneg>numpos) {
               val=val-NyqInterval;
             } else {
               /* Save bin for windowing if unsuccessful after
               ** four passes: */
               if (loopcount>4) GOOD[i][currIndex]=-2;
               if (VERBOSE) printf("Saving for windowing\n");
             }
               } else {
             /* Remove bin: */
             GOOD[i][currIndex]=-2;
             if (VERBOSE) printf("Saving for windowing\n");
               }
             }
           }
         }
           }
         }
       }
     }

}

void second_pass(
    Volume* VALS, Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
    int sweepIndex, int numRays, int numBins, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval, 
    unsigned short *pflag, int *pstep,
    int binindex[8], int rayindex[8], float diffs[8] 
        )
{

    int i, l, direction, 
        left, right, next, prev,
        countindex, numneg, numpos, in, out;
    int step;
    int startindex;
    int endindex;
    /* int abIndex; */
    int currIndex;
    int loopcount;
    
    unsigned short numtimes, flag=1;
    float val, diff, finalval,
        valcheck, goodval;
    float soundval;

    flag = *pflag;
    step = *pstep;

     /* Beginning second pass through the data, this time using 
     **  soundVolume only: */

     if (lastVolume!=NULL && soundVolume!=NULL) {
       flag=1;
       for (currIndex=0;currIndex<numRays;currIndex++) {
         /*if (sweepIndex<numSweeps-1) abIndex=findRay(rvVolume, rvVolume,
             sweepIndex, sweepIndex+1, currIndex, missingVal); */
         for (i=DELNUM;i<numBins;i++) {
           if (GOOD[i][currIndex]==0) {
         val=(float) VALS->sweep[sweepIndex]->ray[currIndex]->
           h.f(VALS->sweep[sweepIndex]->ray[currIndex]->range
               [i]);
         valcheck=val;
         soundval=(float) soundVolume->
           sweep[sweepIndex]->ray[currIndex]->h.f(soundVolume->
            sweep[sweepIndex]->ray[currIndex]->range[i]);
             
         if (soundval!=missingVal&&val!=missingVal) {
           diff=soundval-val;        
           if (diff<0.0) {
             diff=-diff;
             direction=-1;
           } else direction=1;
           numtimes=0;
           while (diff>0.99999*NyqVelocity&&numtimes<=MAXCOUNT) {
             val=val+NyqInterval*direction;
             numtimes=numtimes+1;
             diff=soundval-val;
             if (diff<0.0) {
               diff=-diff;
               direction=-1;
             } else direction=1;
           }
           if (diff<COMPTHRESH2*NyqVelocity&&fabs(valcheck)>CKVAL) {
             /* Return the good value. */
             finalval=(float)rvVolume->sweep[sweepIndex]->
               ray[currIndex]->h.invf(val);
             rvVolume->sweep[sweepIndex]->ray[currIndex]->
               range[i]=(unsigned short) (finalval);
             GOOD[i][currIndex]=1;
           }
         }
           }
         }
       }
           
       /* Now, try to unfold the rest of the GOOD=0 bins, assuming spatial
          continuity: */
       loopcount=0;
       while (flag==1) {
         loopcount=loopcount+1;
         flag=0;
         if (step==1) {
           step=-1;
           startindex=numRays-1;
           endindex=-1;
         } else {
           step=1;
           startindex=0;
           endindex=numRays;
         }
         for (currIndex=startindex;currIndex!=endindex;currIndex=
              currIndex+step) {
           for (i=DELNUM;i<numBins;i++) {
         if (GOOD[i][currIndex]==0) {
           countindex=0;
           val=(float) VALS->sweep[sweepIndex]->ray[currIndex]
             ->h.f(VALS->sweep[sweepIndex]->ray[currIndex]->
               range[i]);
           valcheck=val;
           if (currIndex==0) left=numRays-1;
           else left=currIndex-1;
           if (currIndex==numRays-1) right=0;
           else right=currIndex+1;
           next=i+1;
           prev=i-1;
           /* Look at all bins adjacent to current bin in
              question: */
           if (i>DELNUM) {
             if (GOOD[prev][left]==1) {
               if (flag==0) flag=1;
               countindex=countindex+1;
               binindex[countindex-1]=prev;
               rayindex[countindex-1]=left;
               if (VERBOSE) printf("pl ");
             }
             if (GOOD[prev][currIndex]==1) {
               if (flag==0) flag=1;
               countindex=countindex+1;
               binindex[countindex-1]=prev;
               rayindex[countindex-1]=currIndex;
               if (VERBOSE) printf("pc ");
             }
             if (GOOD[prev][right]==1) {
               if (flag==0) flag=1;
               countindex=countindex+1;
               binindex[countindex-1]=prev;
               rayindex[countindex-1]=right;
               if (VERBOSE) printf("pr ");
             }
           }
           if (GOOD[i][left]==1) {
             if (flag==0) flag=1;
             countindex=countindex+1;
             binindex[countindex-1]=i;
             rayindex[countindex-1]=left;
             if (VERBOSE) printf("il ");
           }
           if (GOOD[i][right]==1) {
             if (flag==0) flag=1;
             countindex=countindex+1;
             binindex[countindex-1]=i;
             rayindex[countindex-1]=right;
             if (VERBOSE) printf("ir ");       
           }
           if (i<numBins-1) {  
             if (GOOD[next][left]==1) {
               if (flag==0) flag=1;
               countindex=countindex+1;
               binindex[countindex-1]=next;
               rayindex[countindex-1]=left;
               if (VERBOSE) printf("nl "); 
             }
             if (GOOD[next][currIndex]==1) {
               if (flag==0) flag=1;
               countindex=countindex+1;
               binindex[countindex-1]=next;
               rayindex[countindex-1]=currIndex;
               if (VERBOSE) printf("nc ");
             }
             if (GOOD[next][right]==1) {
               if (flag==0) flag=1;
               countindex=countindex+1;
               binindex[countindex-1]=next;
               rayindex[countindex-1]=right;
               if (VERBOSE) printf("nr ");
             }
           }
           /* Unfold against all adjacent values with GOOD==1*/
           if (countindex>=1) {
             numtimes=0;
             while(val!=missingVal&&GOOD[i][currIndex]==0) {
               numtimes=numtimes+1;
               in=0;
               out=0;
               numneg=0;
               numpos=0;
               if (VERBOSE) printf("%d: ", countindex); 
               for (l=0;l<countindex;l++) {
             goodval=(float)rvVolume->sweep[sweepIndex]->
               ray[rayindex[l]]->h.f(rvVolume->sweep
                [sweepIndex]->ray[rayindex[l]]->range
                [binindex[l]]);
             diffs[l]=goodval-val;
             if (fabs(diffs[l])<THRESH*NyqVelocity) in=in+1;
             else {
               out=out+1;
               if (diffs[l]>NyqVelocity) {
                 numpos=numpos+1;
               } else if (diffs[l]<-NyqVelocity){
                 numneg=numneg+1;
               }
             }
               }
               if (in>out) {
             finalval=(float)rvVolume->sweep[sweepIndex]->
               ray[currIndex]->h.invf(val);
             rvVolume->sweep[sweepIndex]->ray[currIndex]->
               range[i]=(unsigned short) (finalval);
             GOOD[i][currIndex]=1;
             if (VERBOSE) printf("Val: %4.2f\n", val);
               } else {
             if (numtimes<=MAXCOUNT) {
               if (numpos+numneg<(in+out-(numpos+numneg))) {
                 if (loopcount<=2) val=missingVal; /* Try later */
                 else {
                   /* Keep the value after two passes through
                   ** data */
                   finalval=(float)rvVolume->sweep[sweepIndex]->
                 ray[currIndex]->h.invf(val);
                   rvVolume->sweep[sweepIndex]->ray[currIndex]->
                 range[i]=(unsigned short) (finalval);
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
         }
           }
         }
       }
     }


    *pflag = flag;
    *pstep = step;

}

void unfold_remote(
    Volume* VALS, Volume* rvVolume, Volume* soundVolume, Volume* lastVolume,
    int sweepIndex, int numRays, int numBins, 
    float missingVal, short GOOD[MAXBINS][MAXRAYS],
    float NyqVelocity, float NyqInterval)
{
    int i, direction,
        startray, endray, firstbin, 
        lastbin;
    int currIndex; 
    unsigned short numtimes, wsuccess;
    float val, diff, finalval, winval;
    float std;
  
    
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
         if (diff<THRESH*NyqVelocity) {
           /* Return the value. */
           finalval=(float)rvVolume->sweep[sweepIndex]->
             ray[currIndex]->h.invf(val);
           rvVolume->sweep[sweepIndex]->ray[currIndex]->
             range[i]=(unsigned short) (finalval);
           GOOD[i][currIndex]=1;
         } else if (diff<(1.0 - (1.0 - THRESH)/2.0)*NyqVelocity) {
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
}
