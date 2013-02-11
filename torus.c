


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


double
torus_rho (x)
	     double x[];
{
	double z;
	double rho;

	/* !! ERROR - this is only approximate because of rho2nh */

	z=geo.compton_torus_rmax-geo.compton_torus_rmin;
	rho=geo.compton_torus_tau/(THOMPSON*rho2nh*z);
	return(rho);
}

