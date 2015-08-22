
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                     Space Telescope Science Institute

 Synopsis:
	These are the routines that are specific to an offiset wind model 
	intended to mimic a proposed model for AGN winds.  It is very similar
	to an SV wind except for an offset
Arguments:		

Returns:
 
Description:	

Notes:

	Although this is colloquially called an Elvis model; it is really just an offset
	wind model

	These routines were all developed from the original sv routines by SS and so 
	changes in one may well effect the other.  The velocity flow assumes that one
	calculates a "poloidal distance" starting with a vertical flow until you get
	to the offset height at which point the remaining distance is calculated according
	to that prescribed by the SV model 

	15aug The incorporation of domains is a little complicated for the offset model,
	because the integraal that must be carried out requires a number of parameters
	that need to be passed around the direct system call.  Instead of defining all
	of these individually, I have created a single external variable edom, which
	should only be used within the functios in this routine.  This is a safe only
	so long as no external routines call one of the routines that uses edom

History:
         Oct06 SS
	 15aug	ksl	Modififications to account for domains
**************************************************************/


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_elvis_wind_params gets input data which is necessary for a Shlossman & Vitello 
	description of the wind with a vertical offset as requested for the Elvis AGN model
	This is all copied almost directly from the equivalent SV routines.
Arguments:		

Returns:
 
Description:	
	The parameters, sv...,  needed to define the wind.
	The parameters are identical to that nneed for a SV wind with an additional offset.
Notes:


History:
         Oct06 SS
	 15uag	ksl	Modifications to allow for domains
**************************************************************/

int edom;			/* External variable which allows one to avoid passing the domain for calls within
				   the routines in this file  */

int
get_elvis_wind_params (ndom)
     int ndom;
{
  double windmin, windmax, theta_min, theta_max;

  edom = ndom;

  Log ("Creating an SV/Elvis wind model for an AGN\n");

  rddoub ("wind.mdot(msol/yr)", &geo.wind_mdot);
  geo.wind_mdot *= MSOL / YR;


  zdom[ndom].sv_rmin = 2.8e9;
  zdom[ndom].sv_rmax = 8.4e9;
  zdom[ndom].sv_thetamin = 20. / RADIAN;
  zdom[ndom].sv_thetamax = 65. / RADIAN;
  zdom[ndom].sv_gamma = 1.;
  zdom[ndom].sv_v_zero = 6e5;	/* velocity at base of wind */
  zdom[ndom].sv_r_scale = 7e10;	/*Accleration length scale for wind */
  zdom[ndom].sv_alpha = 1.5;	/* Accleration scale exponent */
  zdom[ndom].sv_v_infinity = 3;	/* Final speed of wind in units of escape velocity */
  zdom[ndom].sv_lambda = 0.0;	/* Mass loss rate exponent */
  zdom[ndom].elvis_offset = 1.e14;

  windmin = zdom[ndom].sv_rmin / geo.rstar;
  windmax = zdom[ndom].sv_rmax / geo.rstar;
  rddoub ("sv.diskmin(wd_rad)", &windmin);
  rddoub ("sv.diskmax(wd_rad)", &windmax);


  zdom[ndom].sv_rmin = windmin * geo.rstar;
  zdom[ndom].sv_rmax = windmax * geo.rstar;

  theta_min = zdom[ndom].sv_thetamin * RADIAN;
  theta_max = zdom[ndom].sv_thetamax * RADIAN;
  rddoub ("sv.thetamin(deg)", &theta_min);
  rddoub ("sv.thetamax(deg)", &theta_max);
  zdom[ndom].sv_thetamin = theta_min / RADIAN;
  zdom[ndom].sv_thetamax = theta_max / RADIAN;

  rddoub ("sv.mdot_r_exponent", &zdom[ndom].sv_lambda);	/* Mass loss rate exponent */
  rddoub ("sv.v_infinity(in_units_of_vescape", &zdom[ndom].sv_v_infinity);	/* Final speed of wind in units of escape velocity */

  rddoub ("sv.acceleration_length(cm)", &zdom[ndom].sv_r_scale);	/*Accleration length scale for wind */
  rddoub ("sv.acceleration_exponent", &zdom[ndom].sv_alpha);	/* Accleration scale exponent */

  rddoub ("elvis_offset(cm)", &zdom[ndom].elvis_offset);	/* This is the vertical offset - the length over which the wind rises vertically */

/* Assign the generic parameters for the wind the generic parameters of the wind */


  zdom[ndom].wind_rmin = geo.rstar;
  zdom[ndom].wind_rmax = geo.rmax;
  zdom[ndom].wind_rho_min =
    zdom[ndom].sv_rmin -
    (zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamin));
  zdom[ndom].wind_rho_max =
    zdom[ndom].sv_rmax -
    (zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamin));
  zdom[ndom].wind_thetamin = zdom[ndom].sv_thetamin;
  zdom[ndom].wind_thetamax = zdom[ndom].sv_thetamax;

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
    {
      zdom[ndom].xlog_scale = zdom[ndom].sv_rmin;
      zdom[ndom].zlog_scale = 1e15;	/* Big number - for AGN */
    }

/*Now calculate the normalization factor for the wind*/

  zdom[ndom].mdot_norm =
    qromb (elvis_wind_mdot_integral, zdom[ndom].sv_rmin, zdom[ndom].sv_rmax,
	   1e-6);
  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_velocity(ndom, x,v) calulates the v of an offset Schlossman Vitello wind from a position
	x
Arguments:		
	ndom			the relevant domain
	double x[]		the position where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	
		
Notes:

History:
 
**************************************************************/

double
elvis_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, rzero, theta, speed;
  double ldist, zzz, v_escape, vl;
  struct photon ptest;
  double xtest[3];
  double s;

  edom = ndom;			// External variable used only in the evlvis routines

  zzz = v_escape = -99.;

  rzero = elvis_find_wind_rzero (ndom, x);
  theta = elvis_theta_wind (ndom, rzero);

  /* Calculate the poloidal distance, assuming a stream line that is vertical 
   * up to the elvis offset and then along an SV stream line after that.  ksl 111123 
   */

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  if (fabs (x[2]) > zdom[ndom].elvis_offset)
    {
      ldist =
	zdom[ndom].elvis_offset + sqrt ((r - rzero) * (r - rzero) +
					(x[2] -
					 zdom[ndom].elvis_offset) * (x[2] -
								     zdom
								     [ndom].
								     elvis_offset));
    }
  else
    {
      ldist = fabs (x[2]);
    }

  /* If the disk is vertically extended ksl 111124 
   * ERROR? - This calculation of the poloical distance replaces the one above.  It does not take the offset into account*/
  if (geo.disk_type == 2)
    {
      xtest[0] = r;		// Define xtest in the +z plane
      xtest[1] = 0;
      xtest[2] = fabs (x[2]);
      ptest.x[0] = rzero;
      ptest.x[1] = 0.0;
      ptest.x[2] = EPSILON;
      ptest.lmn[0] = sin (theta);	// 56d -- ptest direction is along stream line
      ptest.lmn[1] = 0.0;
      ptest.lmn[2] = cos (theta);
      s = ds_to_disk (&ptest, 1);
      move_phot (&ptest, s);	// Now test photon is at disk surface
      vsub (ptest.x, xtest, xtest);
      ldist = length (x);
    }


  /* Having calculated the "poloidal distance" calculate the velocity ksl 111124 */
  vl = zdom[ndom].sv_v_zero;
  if (ldist > 0)
    {
      zzz = pow (ldist / zdom[ndom].sv_r_scale, zdom[ndom].sv_alpha);

      if (rzero < geo.rstar)
	v_escape = sqrt (2. * G * geo.mstar / geo.rstar);
      else
	v_escape = sqrt (2. * G * geo.mstar / rzero);

      vl =
	zdom[ndom].sv_v_zero + (zdom[ndom].sv_v_infinity * v_escape -
				zdom[ndom].sv_v_zero) * zzz / (1. + zzz);
    }

  /* Now calculate the direction of the velocity */
  if (fabs (x[2]) > zdom[ndom].elvis_offset)
    {
      v[0] = vl * sin (theta);
      v[2] = vl * cos (theta);
    }
  else
    {
      v[0] = 0;
      v[2] = vl;
    }


  if (r > 0)
    v[1] = sqrt (G * geo.mstar * rzero) / r;
  else
    v[1] = 0;

  if (x[2] < 0)
    v[2] *= (-1);
  /* 04aug -- ksl --52 At this point we have calculated the velocity in the xz plane, which
   * is identical to the statement that we have calculated it in
   * cylindrical coordinates.  The next bit projects back to xyz 
   * coordinates if x was not originally in the xz plane.
   */
  if (x[1] != 0.0)
    {
      project_from_cyl_xyz (x, v, xtest);
      stuff_v (xtest, v);
    }


  speed = (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
  if (sane_check (speed))
    {
      Error ("elvis_velocity:sane_check x %f %f %f v %f %f %f\n", x[0], x[1],
	     x[2], v[0], v[1], v[2]);
      Error
	("elvis_velocity: rzero %f theta %f ldist %f zzz %f v_escape %f vl %f\n",
	 rzero, theta, ldist, zzz, v_escape, vl);
    }


  return (speed);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_rho(x) calculates the density of an offset sv_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the denisty
Returns:
	The density at x is returned in gram/cm**3
	
Description:	
		
Notes:

History:
 
**************************************************************/

double
elvis_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rzero, theta;
  double ldist;
  double dmdot_da;
  double dtheta_drzero, dr_drzero;

  double v[3], rho;
  struct photon ptest;
  double xtest[3];
  double s;


  elvis_velocity (ndom, x, v);

  rzero = elvis_find_wind_rzero (ndom, x);

  /* ERROR - This is a mistake, and is a root cause of some of the problems we are seeing in the Elvis wind 111124 */
  if (rzero > zdom[ndom].sv_rmax)
    {
      rho = 1.e-20;		//This is a fudge which keeps the wind boundaries conical by filling in excess space with
      // "empty" space - a waste of memory and time but to do otherwise would require a 
      // more substantial change to the way wind boundaries are defined SSOct06
      return (rho);
    }
  if (rzero <
      (zdom[ndom].sv_rmin +
       (zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamin))))
    {
      rho = 1.e-20;		//fudge, as above
      return (rho);
    }




  theta = elvis_theta_wind (ndom, rzero);

  r = sqrt (x[0] * x[0] + x[1] * x[1]);
  if (fabs (x[2]) > zdom[ndom].elvis_offset)
    {
      ldist =
	zdom[ndom].elvis_offset + sqrt ((r - rzero) * (r - rzero) +
					(x[2] -
					 zdom[ndom].elvis_offset) * (x[2] -
								     zdom
								     [ndom].
								     elvis_offset));
    }
  else
    {
      ldist = fabs (x[2]);
    }

  if (geo.disk_type == 2)
    {
      xtest[0] = r;		// Define xtest in the +z plane
      xtest[1] = 0;
      xtest[2] = fabs (x[2]);
      ptest.x[0] = rzero;
      ptest.x[1] = 0.0;
      ptest.x[2] = EPSILON;
      ptest.lmn[0] = cos (theta);
      ptest.lmn[1] = 0.0;
      ptest.lmn[2] = sin (theta);
      s = ds_to_disk (&ptest, 1);
      move_phot (&ptest, s);	// Now test photon is at disk surface
      vsub (ptest.x, xtest, xtest);
      ldist = length (xtest);
      rzero = length (ptest.x);
    }

  if (fabs (x[2]) < zdom[ndom].elvis_offset)
    {
      dmdot_da =
	zdom[ndom].wind_mdot * pow (rzero,
				    zdom[ndom].sv_lambda) /
	zdom[ndom].mdot_norm / 2.;
      rho = dmdot_da / v[2];

      return (rho);
    }

  /* Reduced by a factor of 2 to account for radiation on both sides of the disk */
  dmdot_da =
    zdom[ndom].wind_mdot * pow (rzero,
				zdom[ndom].sv_lambda) * cos (theta) /
    zdom[ndom].mdot_norm / 2.;

  /* Although the current definition of sv_theta_wind is continuous, the derivative is not continuous accross the
     outer boundary of the wind and thus dtheta_drzero would appear to change discontinuously.   This created
     large apparent density jumps at the outside edge of the wind. We can't allow that and so we force the derivative to equal
     that at the edge of the wind if you try to calculate the density rho.  ksl 97apr23 */

  if (rzero > zdom[ndom].sv_rmax)
    rzero = zdom[ndom].sv_rmax;
  dtheta_drzero =
    (elvis_theta_wind (ndom, rzero) -
     elvis_theta_wind (ndom, (1. - EPSILON) * rzero)) / (EPSILON * rzero);

  dr_drzero = 1. + ldist * dtheta_drzero / cos (theta);
  /* Note VS93 eqn 8 is dr/drzero but equation  7 is drzero/dr   ksl 97 apr 19 */
  rho = rzero * dmdot_da / (dr_drzero * r * v[2]);

  return (rho);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	elvis_find_wind_rzero(p) locates the radial position in the disk from which the 
	stream line arises for an offset SV wind (the "elvis" model).  
 
Arguments:		
	double p[]	    A 3 vector position, and not a photon structure
Returns:
	elvis_find_wind_rzero returns the radius of the footpoint of the stream line
 
Description:	
		
Notes:

History:
 
**************************************************************/


double
elvis_find_wind_rzero (ndom, p)
     int ndom;
     double p[];		/* Note that p is a 3 vector and not a photon structure */
{
  double x, z;
  double rho_min, rho_max, rho;

  /* thetamin and theta max are defined w.r.t  z axis */


  z = fabs (p[2]);		/* This is necessary to get correct answer above
				   and below plane */

  /* If the position is underneath the elvis_offset, the stream line is assumed
   * to be vertical and so you just return rho 111124 ksl */
  if (z <= zdom[ndom].elvis_offset)
    {
      x = (sqrt (p[0] * p[0] + p[1] * p[1]));	// If p is in the xy plane, there is no need to think further
      return (x);
    }


  elvis_zero_init (p);		/* This initializes the routine sv_zero_r.  It is not
				   actually needed unless zbrent is called, but it
				   does allow you to check your answer otherwize

				   It does not require a domain number
				 */

  /* The next lines provide a value for the footpoint position when the positions
   * we want the footpoint for is outside of the wind 111124 ksl */
  rho_min = zdom[ndom].sv_rmin + z * tan (zdom[ndom].sv_thetamin);
  rho_max = zdom[ndom].sv_rmax + z * tan (zdom[ndom].sv_thetamax);
  rho = sqrt (p[0] * p[0] + p[1] * p[1]);

  /* Note that theta_min and theta_max are measured from the z axis */

  if (rho <= rho_min)
    {
      x = zdom[ndom].sv_rmin * rho / rho_min;
      return (x + (zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamin)));
    }
  if (rho >= rho_max)
    {
      x = zdom[ndom].sv_rmax + rho - rho_max;
      return (x + (zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamax)));
    }

  /* 100 here means that zbrent will end if one has a guess of rzero which is
     correct ot 100 cm */

  x =
    zbrent (elvis_zero_r,
	    zdom[ndom].sv_rmin +
	    (zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamin)),
	    zdom[ndom].sv_rmax +
	    (zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamax)), 100.);
  return (x);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	The next two subroutines are used to guess the position in the disk
	from which a streamline arises.  They work in the same manner as sv_zero_r 
	and sv_zero_init.
 
Arguments:		

Returns:
 
Description:	
	
		
Notes:

	This routine does not require a domain, indeed it is so short
	it is not entirely clear why it is required.  It simple sets 
	and internal vector zero_p

History:
 
**************************************************************/

double zero_p[3];

int
elvis_zero_init (p)
     double p[];
{
  stuff_v (p, zero_p);
  zero_p[2] = fabs (zero_p[2]);	/* Required to get correct 
				   answer below (in -z ) the disk */
  return (0);
}

/* This routine is used to test whether a guess of r_zero is correct.  If
   you have the answer exactly then elvis_zero_r will return 0 */

double
elvis_zero_r (r)
     double r;
{
  double theta;
  double rho, rho_guess;

  theta = elvis_theta_wind (edom, r);

  rho = sqrt (zero_p[0] * zero_p[0] + zero_p[1] * zero_p[1]);
  if (zero_p[2] < zdom[edom].elvis_offset)
    {
      rho_guess = r;
    }
  else
    {
      rho_guess = r + tan (theta) * (zero_p[2] - zdom[edom].elvis_offset);
    }
  return (rho_guess - rho);


}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_theta_wind(r) finds the angle at which the wind that rises 
	from a specific point on the disk bend to
 
Arguments:
	ndom   		the domain where the parameters are stored
	double r	a radial distance on the surface of the disk		

Returns:
 
Description:	
		
Notes:


History:
	15aug	ksl	Modifications made for domains.
 
**************************************************************/

double
elvis_theta_wind (ndom, r)
     int ndom;
     double r;
{
  double theta;
  if (r <=
      (zdom[ndom].sv_rmin +
       zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamin)))
    return (atan (tan (zdom[ndom].sv_thetamin * r / zdom[ndom].sv_rmin)));
  if (r >=
      (zdom[ndom].sv_rmax +
       zdom[ndom].elvis_offset * tan (zdom[ndom].sv_thetamax)))
    return (zdom[ndom].sv_thetamax);
  theta = zdom[ndom].sv_thetamin +
    (zdom[ndom].sv_thetamax -
     zdom[ndom].sv_thetamin) * pow ((r - zdom[ndom].sv_rmin -
				     zdom[ndom].elvis_offset *
				     tan (zdom[ndom].sv_thetamin)) /
				    (zdom[ndom].sv_rmax +
				     zdom[ndom].elvis_offset *
				     tan (zdom[ndom].sv_thetamax) -
				     zdom[ndom].sv_rmin -
				     zdom[ndom].elvis_offset *
				     tan (zdom[ndom].sv_thetamin)),
				    zdom[ndom].sv_gamma);
  return (theta);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double elvis_wind_mdot_integral(r) is the integrand of SV model for mdot as a function
	of radius - works in the same way as the SV routine but calls the elvis version of the theta routines. 
Arguments:		
	double r;
 
Returns:
 
Description:	
	
		
Notes:

	The relevant domain must stored in the external variable edom

	Note that this calls elvis_theta_wind()

History:

 
**************************************************************/

double
elvis_wind_mdot_integral (r)
     double r;
{
  double x;

  if (r <
      (zdom[edom].sv_rmin +
       zdom[edom].elvis_offset * tan (zdom[edom].sv_thetamin)))
    {
      x = 0.0;
      return (x);
    }
  if (r > zdom[edom].sv_rmax)
    {
      x = 0.0;
      return (x);
    }


  x =
    2 * PI * pow (r,
		  zdom[edom].sv_lambda + 1.) * cos (elvis_theta_wind (edom,
								      r));
  return (x);

}



/***********************************************************
	Space Telescope Science Institute

 Synopsis:
 	ds_to_pillbox calculates the distance to a pillbox, really
	an annular region with an inner radius and outer radius
	and a height.  The pillbox is symmetric with respect to the xy
Arguments:		
	pp	photon ptr containing a position and direction
	rmin	The minimum rho for the pillbox
	rmax	The maximum rho for the pillbox
	height 	Height of the pillbox
 
Returns:
 
Description:	
	
		
Notes:
	No allowance is made for a geometrically thick disk.  This is
	a purely geometrical construction.

	This is modeled on ds_to_torus, but here all the parameters
	are passed.  It's likely that this routine should replace
	that one

	This should almost surely be moved to phot_util

History:
	11nov	ksl	Began coding to deal with issues that
			currently exist with the elvis model

 
**************************************************************/

double
ds_to_pillbox (pp, rmin, rmax, height)
     PhotPtr pp;
     double rmin, rmax, height;
{

  struct photon ptest;
  double ds, ds_best, x;
  struct plane xplane;


  ds_best = VERY_BIG;

  /* Make sure we don't mess with pp */
  stuff_phot (pp, &ptest);
  ds = ds_to_cylinder (rmin, &ptest);

  /* Calculate the distance to the innner cylinder */
  if (ds < VERY_BIG)
    {
      /* Check whether we encounted the
       * part of the cylinder we are interested in
       */
      move_phot (&ptest, ds);
      if (fabs (ptest.x[2]) < height)
	{
	  ds_best = ds;
	}
      /* Now reinitialize ptest */
      stuff_phot (pp, &ptest);
    }

  /* Similarly calculate the distance to the outer
   * cylinder
   */
  ds = ds_to_cylinder (rmax, &ptest);
  if (ds < ds_best)
    {
      move_phot (&ptest, ds);
      if (fabs (ptest.x[2]) < height)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  /* At this point we know whether the photon has interecepted
   * the wall of the cylinder, but we do not know if it intercepted
   * the top or bottom of the cylinder earlier
   */

  xplane.x[0] = 0;
  xplane.x[1] = 0;
  xplane.x[2] = height;
  xplane.lmn[0] = 0.0;
  xplane.lmn[1] = 0.0;
  xplane.lmn[2] = 1.0;

  ds = ds_to_plane (&xplane, &ptest);
  // Note that ds to plane can return a negative number
  if (ds > 0 && ds < ds_best)
    {
      move_phot (&ptest, ds);
      x = fabs (ptest.x[0]);
      if (rmin < x && x < rmax)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  xplane.x[2] = (-height);

  ds = ds_to_plane (&xplane, &ptest);
  if (ds > 0 && ds < ds_best)
    {
      move_phot (&ptest, ds);
      x = fabs (ptest.x[0]);
      if (rmin < x && x < rmax)
	{
	  ds_best = ds;
	}
      stuff_phot (pp, &ptest);
    }

  return (ds_best);
}
