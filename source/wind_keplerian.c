/* The routines in this file define and summarize the properties of the wind.  The routines here are
   specific to the SV description of a wind. It is only useful in the 2-d version of the code.

   This file was created in 98apr in order to being to isolate the SV description from the more
   generic parts of the wind.  Major modifications, mostly moving new code here and
   creating the remaining subroutines to completely concentrate the sv dependent routines
   here were made in 98dec.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_sv_wind_params gets input data which is necessary for a Shlossman & Vitello 
	description of the wind
Arguments:		

Returns:
 
Description:	
	The parameters, geo.sv...,  obtained here are only used in the routines in stellar_winds.c
	which calculate the velocity and density of the wind during the initialization process.
	Other portions of the structure, geo defined here are more general purpose.		
Notes:


History:
 	98dec	ksl	Coded and debugged as part of major change in IO structure required when
 				adding a spherical wind
        080518  ksl     60a - geo should contain only cgs units
	11aug	ksl	70b - kluge to get better xscale with compton torus
**************************************************************/
double d_wind_keplerian_density;
double d_wind_keplerian_temperature;

int
get_sv_wind_params ()
{
  double windmin, windmax, theta_min, theta_max;
  double qromb (), sv_wind_mdot_integral ();

  Log ("Creating an SV wind model for a Cataclysmic Variable\n");

  geo.wind_mdot /= (MSOL / YR);	// Convert to MSOL/YR for easy of data entry
  rddoub ("wind.mdot(msol/yr)", &geo.wind_mdot);
  geo.wind_mdot *= MSOL / YR;


  geo.sv_rmin = 2.8e9;
  geo.sv_rmax = 8.4e9;
  geo.sv_thetamin = 20. / RADIAN;
  geo.sv_thetamax = 65. / RADIAN;
  geo.sv_gamma = 1.;
  geo.sv_v_zero = 6e5;		/* velocity at base of wind */
  geo.sv_r_scale = 7e10;	/*Accleration length scale for wind */
  geo.sv_alpha = 1.5;		/* Accleration scale exponent */
  geo.sv_v_infinity = 3;	/* Final speed of wind in units of escape velocity */
  geo.sv_lambda = 0.0;		/* Mass loss rate exponent */

  windmin = geo.sv_rmin / geo.rstar;
  windmax = geo.sv_rmax / geo.rstar;
  rddoub ("sv.diskmin(units_of_rstar)", &windmin);
  rddoub ("sv.diskmax(units_of_rstar)", &windmax);


  geo.sv_rmin = windmin * geo.rstar;
  geo.sv_rmax = windmax * geo.rstar;

  theta_min = geo.sv_thetamin * RADIAN;
  theta_max = geo.sv_thetamax * RADIAN;
  rddoub ("sv.thetamin(deg)", &theta_min);
  rddoub ("sv.thetamax(deg)", &theta_max);
  geo.sv_thetamin = theta_min / RADIAN;
  geo.sv_thetamax = theta_max / RADIAN;

  rddoub ("sv.mdot_r_exponent", &geo.sv_lambda);	/* Mass loss rate exponent */
  rddoub ("sv.v_infinity(in_units_of_vescape", &geo.sv_v_infinity);	/* Final speed of wind in units of escape velocity */

  rddoub ("sv.acceleration_length(cm)", &geo.sv_r_scale);	/*Accleration length scale for wind */
  rddoub ("sv.acceleration_exponent", &geo.sv_alpha);	/* Accleration scale exponent */
/* Assign the generic parameters for the wind the generic parameters of the wind */

  geo.wind_rmin = geo.rstar;
  geo.wind_rmax = geo.rmax;
  geo.wind_rho_min = geo.sv_rmin;
  geo.wind_rho_max = geo.sv_rmax;
  geo.wind_thetamin = geo.sv_thetamin;
  geo.wind_thetamax = geo.sv_thetamax;
  geo.xlog_scale = geo.sv_rmin;

  /* !! 70b - This change is to accomodate the torus, but it is not obvious this is the
   * best way to set the scales now. It might be better do do this in make_grid!!  */
  if (geo.compton_torus && geo.compton_torus_rmin < geo.xlog_scale)
    {
      geo.xlog_scale = geo.compton_torus_rmin;
    }

/*70d - ksl - This change made to give one a chance of being able to do an 
   agn and a CV with the sv model.  The underlying assumption is that the
   effective radius provides a good scale factor in the verticla direction.
   An alternative would be to use sv_rmin.
 */

//OLD70d  geo.zlog_scale = 1e7;
  geo.zlog_scale = geo.rstar;

/*Now calculate the normalization factor for the wind*/

  geo.mdot_norm =
    qromb (sv_wind_mdot_integral, geo.sv_rmin, geo.sv_rmax, 1e-6);
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double sv_velocity(x,v) calulates the v of a Schlossman Vitello wind from a position
	x
Arguments:		
	double x[]		the postion for which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	
		
Notes:

History:
 	98dec	ksl	Coded as part of effort to isolate SV portions of the code.
	01mar	ksl	Moved mdot_norm calculation from sv_rho to get_sv_params() 
			since it doesn't need to be recalulated every time
	04aug	ksl	52 -- Modified to allow for s thick disk, and to return
			the velocity in cartesian rather than cylindrical coords.
			These are the same if x is in xz plane.
	04dec	ksl	54a  -- Minor mod to avoid -O3 warning
	05jul	ksl	56d  -- Corrected error in the direction calculated for
			ptest
 
**************************************************************/

double wind_keplerian_velocity (double x[], double v[])
{
  double r, speed;
  struct photon ptest;
  double xtest[3];
  
  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  
  v[0] = 0;
  if (r > 0)
    v[1] = sqrt (G * geo.mstar * geo.sv_rmin) / r;
  else
    v[1] = 0;
  v[2] = 0;

  if (x[1] != 0.0)
    {
      project_from_cyl_xyz (x, v, xtest);
      stuff_v (xtest, v);
    }
  speed = (sqrt (v[0] * v[0] + v[1] * v[1]));
  return (speed);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double sv_rho(x) calculates the density of an sv_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the denisty
Returns:
	The density at x is returned in gram/cm**3
	
Description:	
		
Notes:
	52 -- If the disk is thick, our intention is that the stream lines
	trace back to the disk, and have the same properties at that radius
	as they would have had except for a vertical offset.

History:
 	98dec	ksl	Coded as part of effort to isolate SV portions of the code.
 
**************************************************************/
double wind_keplerian_rho (double x[])
{
  return (d_wind_keplerian_density);
}
