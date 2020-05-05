/***********************************************************/
/** @file  frame.c
 * @author ksl
 * @date   May, 2020
 *
 * @brief
 * These routines have to do with switching from the observer 
 * to the local frame and back
 *
 * ###Notes###
 *
 * These routines are written so that if p and p_obs point 
 * to the same variable then everything occurs in place
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


double
observer_to_local_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double f;
  double v[3], vel;

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);

  vel = dot (p_in->lmn, v);
  f = p_out->freq = p_in->freq * (1. - vel / VLIGHT);

  return (f);
}



double
local_to_observer_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double f;
  double v[3], vel;

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);

  vel = dot (p_in->lmn, v);
  f = p_out->freq = p_in->freq / (1. - vel / VLIGHT);

  return (f);
}

double
local_to_observer_frame_disk (p_in, p_out)
     PhotPtr p_in, p_out;
{
  double f;
  double v[3], vel;

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);


  /* Calculate the local velocity of the disk at this position */
  vdisk (p_in->x, v);
  vel = dot (p_in->lmn, v);


  f = p_out->freq = p_in->freq / (1. - vel / VLIGHT);

  return (f);
}


/**********************************************************/
/**
 * @brief      calculate the  shift given the direction of the incoming
 * 	and outgoing photons and the local  velocity of the wind
 *
 * @param [in] PhotPtr  pin   The pre-scattered  photon (used for its direction)
 * @param [in,out] PhotPtr  pout   the scattered  photon (used for its direction)
 * @param [in] double  v[]   the velocity of the wind where the scatter occurred
 * @param [in] int  nres   either the number of the scatter
 * @return    Always returns 0
 *
 * pout->freq is updated
 *
 * @details
 * Given the incoming and the outgoing photon direction,
 * doppler calculates the outgoing frequency (to first order in beta).
 * There are two basic cases, depending on whether it was a resonant
 * or a nonresonant scatter.
 *
 * ### Notes ###
 * Called from scatter
 *
 **********************************************************/

int
doppler (p_in, p_out, nres)
     PhotPtr p_in, p_out;
     int nres;

{
  double dot ();
  WindPtr one;
  int ndom;
  double v[3];

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);




  if (nres == -1)               //Electron scattering (SS)
  {                             /*It was a non-resonant scatter */
    p_out->freq = p_in->freq * (1 - dot (v, p_in->lmn) / VLIGHT) / (1 - dot (v, p_out->lmn) / VLIGHT);


  }
  else if (nres > -1 && nres < nlines)
  {                             /* It was a resonant scatter. */
    p_out->freq = lin_ptr[nres]->freq / (1. - dot (v, p_out->lmn) / VLIGHT);
  }
  else if ((nres > NLINES && nres < NLINES + nphot_total + 1) || nres == -2)
    /* It was continuum emission - new comoving frequency has been chosen by
       the matom/kpkt routine, but now need to convert in the same way
       as for lines (SS) */
  {
    /*
       If a 2-level atom run, one should never arrive here.
       Just do a check that all is well - this can be removed eventually (SS)
     */
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      Error ("doppler: Not using macro atoms but trying to deexcite one? Abort.\n");
      Exit (0);
    }
    p_out->freq = p_out->freq / (1. - dot (v, p_out->lmn) / VLIGHT);
  }
/* Now do one final check that nothing is awry.  This is another
 * check added by SS that should probably be deleted or done before this point.
 * I have made it fatal so that we will pay attention to it if it occurs. ksl */

  else
  {
    Error ("doppler: nres %d > NLINES + nphot_total %d\n", nres, NLINES + nphot_total);
    Exit (0);
  }

  return (0);

}
