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

#ifndef FDD_H
#define FDD_H
#define MAXRAYS   500      /* added by SRB 980310 */
#define MAXBINS 2048

#define PROXIMITY 5 /* For unfolding using windowing.*/
#define COMPTHRESH2 0.49 /* The threshold for performing initial dealiasing 
			** using sounding (or VAD). */
#define MINGOOD 5 /* Number of good values required within unfolding window
		  **  to unfold the current bin. */
#define STDTHRESH 0.8 /* Fraction of the Nyquist velocity to use as a standard
		      **  deviation threshold when windowing. */
#define DELNUM 0 /* The first DELNUM velocity bins will be deleted along each
		 **  ray (should be between 0 and 5). */
#define MAXCOUNT 10 /* Maximum number of folds. */
#include <rsl.h> /* Sweep */ 

float ray_val(Ray *ray, int index);
void ray_set(Ray *ray, int index, float val); 

#endif /* DEALIAS_H */
