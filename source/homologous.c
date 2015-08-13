

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_homologous_params gets input data which is necessary for a homologous expansion law
   		v(r)= r / t0

Arguments:		

Returns:
 
Description:	
	The parameters, cl...,  obtained here are only used in the routines in homologous.c
	which calculate the velocity and density of the wind during the initialization process.
	Other portions of the structure, 

Notes:


History:
        13jul   sas     Created for SN test problem. Based on stellar_wind.c
	15aug	ksl	Modified to accept a domain number
**************************************************************/


int
get_homologous_params (ndom)
     int ndom;
{
  DomainPtr one;
  one = &zdom[ndom];

  Log ("Creating a homolgous wind model\n");





  one->stellar_wind_mdot = 100.;
  one->wind_rmin = geo.rstar;
  one->cl_v_zero = 200e5;
  one->cl_beta = 7.0;

  one->stellar_wind_mdot /= MSOL / YR;
  rddoub ("homologous_boundary_mdot(msol/yr)", &one->stellar_wind_mdot);
  one->stellar_wind_mdot *= MSOL / YR;

  rddoub ("homologous.radmin(cm)", &one->wind_rmin);	/*Radius where wind begins */
  if (geo.wind_rmin < geo.rstar)
    {
      Error
	("get_homologous_params: It is unreasonable to have the wind start inside the star!\n");
      Log ("Setting geo.wind_rmin to geo.rstar\n");
      one->wind_rmin = geo.rstar;
    }
  one->cl_rmin = one->wind_rmin;

  rddoub ("homologous.vbase(cm)", &one->cl_v_zero);	/* Velocity at base of the wind */
  rddoub ("homologous.density_exponent", &one->cl_beta);	/* Density law exponent */


/* Assign the generic parameters for the wind the generic parameters of the wind */
  one->wind_thetamin = 0.0;
  one->wind_thetamax = 90. / RADIAN;

/* define the the variables that determine the gridding */
  one->wind_rho_min = 0;
  one->wind_rho_max = geo.rmax;

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      one->xlog_scale = 0.3 * geo.rstar;
      one->zlog_scale = 0.3 * geo.rstar;
    }

  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double homologous_velocity(x,v) calulates the v the wind at a position 
	x (where both x and v in cartesian coordinates)
Arguments:		
	double x[]		the postion where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	
	The model is homolgous expansion

	v(r)=vmin r / R
	
	The values of the individiual constants should all be part of the structure geo.

	Vmin:  			cl_v_zero;		velocity at base of wind 
	R			cl_rmin	       	the inner radius of the wind

		
Notes:

History:
        13jul   sas     Coded for homolgous test
	15aug	ksl	Added domains support.  Note that this routine
			is only used to set up the grid and so we do not
			need to determin which domain we are in from x
 
**************************************************************/

double
homologous_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, speed;

  DomainPtr one;
  one = &zdom[ndom];

  if ((r = length (x)) == 0.0)
    {
      v[0] = 0.0;
      v[1] = 0.0;
      v[2] = 0.0;
      return (0.0);
    }


  if (r <= geo.rstar || r <= one->cl_rmin)
    speed = one->cl_v_zero;
  else
    {
      speed = r * one->cl_v_zero / one->cl_rmin;
    }
  v[0] = speed * x[0] / r;
  v[1] = speed * x[1] / r;
  v[2] = speed * x[2] / r;

  return (speed);

}

/***********************************************************
		Space Telescope Science Institute

 Synopsis:
	double homologous_rho(x) calculates the density of an stellar_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the density
Returns:
	The density at x is returned in gram/cm**3
	
Description:

	rho=rho_0;	
		
Notes:

History:
 	13Jul	sas	Modified from stellar_wind.c portions of the code.
	15aug 	ksl	Modified to accept domains
 
**************************************************************/

double
homologous_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rho, length ();
  DomainPtr one;
  one = &zdom[ndom];


  r = length (x);

  rho =
    one->stellar_wind_mdot / (4. * PI * one->cl_rmin * one->cl_rmin *
			      one->cl_v_zero) * pow ((r / one->cl_rmin),
						     -1. * one->cl_beta);

  return (rho);
}
