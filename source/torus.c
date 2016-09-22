


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  These are routines coded specifically for including a scattering
  'torus' in a system.

  Description:	

  Arguments:  


  Returns:

  Notes:

 

  History:
11aug	ksl	Coding began     

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"


/* These are the routines to determine whether a photon hits the torus
as defined as a region defined as a kind of pillbox centered ont he
disk.

This is not the best place to put some of these routines.  In particular
ds_to_cylinder should probably be part of phot_util, and ds_to_torus
may belong in photon2d (where ds_to_wind) is located

ds_to_cylinder is modeled on ds_to_sphere

110930	ksl	Began coding
*/

double
ds_to_cylinder (rho, p)
     double rho;
     struct photon *p;
{
  double a, b, c, root[2];
  int i;

  a = 1. - p->lmn[2] * p->lmn[2];
  b = 2. * (p->lmn[0] * p->x[0] + p->lmn[1] * p->x[1]);
  c = (p->x[0] * p->x[0] + p->x[1] * p->x[1]) - rho * rho;

  i = quadratic (a, b, c, root);

/* root[i] is the smallest positive root unless i is
negative in which case either both roots were negative or 
both roots were imaginary */

  if (i >= 0)
    return (root[i]);
  return (VERY_BIG);
}
