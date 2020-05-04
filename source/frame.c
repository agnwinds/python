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
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


double
observer_to_local_frame (p)
     PhotPtr p;
{
  WindPtr one;
  int ndom;
  double f;
  double v[3],vel;

  one = &wmain[p->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p, v);
  vel = dot (p->lmn, v);
  f = p->freq * (1. - vel / VLIGHT);
  return (f);
}



double
local_to_observer_frame (p)
     PhotPtr p;
{
  WindPtr one;
  int ndom;
  double f;
  double v[3],vel;

  one = &wmain[p->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p, v);

  vel = dot (p->lmn, v);
  f = p->freq / (1. - vel / VLIGHT);
  return (f);
}
