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
