
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	get_corona_params gets input data which is necessary to describe a corona above
	the surface of the disk.        
Arguments:		

Returns:
 
Description:	
	The parameters obtained here are only used in the routines in corona.c
	Initially, we define it as a gaussian ring with and exponetial density distribution.   
Notes:
	The terminology here is a bit confusing because rmin and rmax as read in to the routine
	refer actually to rho, which is easy to get confused with other parameters like zdom[ndom].rmax
	which really are radii for the central source.


History:
 	00sep	ksl	Coding begun
	04jun	ksl	Moved from python.c to provide better uniformity of what files contain.
	15aug	ksl	Began mods to accomodate domains
	16apr	ksl	Added modifications for a uniform density grid
**************************************************************/


int
get_corona_params (ndom)
     int ndom;
{
  Log ("Creating a corona above a disk\n");

  // Start with reasonable values for everything which is important

  zdom[ndom].wind_thetamin = 0.0;
  zdom[ndom].wind_thetamax = 0.0;
  zdom[ndom].rmin = geo.rstar;

  zdom[ndom].corona_rmin = 1.e10;
  zdom[ndom].corona_rmax = 2.e10;
  zdom[ndom].corona_zmax = 2.e09;
  zdom[ndom].corona_base_density = 1.e13;
  zdom[ndom].corona_scale_height = 1.e9;

  rddoub ("corona.radmin(cm)", &zdom[ndom].corona_rmin);        /*Radius where corona begins */
  if (zdom[ndom].corona_rmin < geo.rstar)
  {
    Error ("get_corona_params: It is unreasonable to have the corona start inside the star!\n");
    Log ("Setting geo.corona_rmin to geo.rstar\n");
    zdom[ndom].corona_rmin = geo.rstar;
  }
  rddoub ("corona.radmax(cm)", &zdom[ndom].corona_rmax);        /*Radius where corona ends */
  rddoub ("corona.zmax(cm)", &zdom[ndom].corona_zmax);  /*Veritical heighe where corona ends */
  rddoub ("corona.base_den(cgs)", &zdom[ndom].corona_base_density);     /*Density at the base of the corona */
  rddoub ("corona.scale_height(cm)", &zdom[ndom].corona_scale_height);  /*Scale height of corona */
  rddoub ("corona.vel_frac", &zdom[ndom].corona_vel_frac);      /*fractional radial velocity of corona */

  zdom[ndom].rmin = zdom[ndom].corona_rmin;
  /* rmax here is the defines a radius beyond which this region does not exist, if the vertical height is large
   * compared to the horizontal size then one needs to include both */
  zdom[ndom].rmax = sqrt (zdom[ndom].corona_rmax * zdom[ndom].corona_rmax + zdom[ndom].corona_zmax * zdom[ndom].corona_zmax);
  zdom[ndom].zmax = zdom[ndom].corona_zmax;
  zdom[ndom].wind_rho_min = zdom[ndom].corona_rmin;
  zdom[ndom].wind_rho_max = zdom[ndom].corona_rmax;
  zdom[ndom].wind_thetamin = 0.0;
  zdom[ndom].wind_thetamax = 0.0;

  /* Set up wind planes for a layer with a specific height */

  zdom[ndom].windplane[0].x[0] = zdom[ndom].windplane[0].x[1] = zdom[ndom].windplane[0].x[2] = 0;
  zdom[ndom].windplane[0].lmn[0] = zdom[ndom].windplane[0].lmn[1] = 0;
  zdom[ndom].windplane[0].lmn[2] = 1;

  zdom[ndom].windplane[1].x[0] = zdom[ndom].windplane[0].x[1] = 0;
  zdom[ndom].windplane[1].x[2] = zdom[ndom].corona_zmax;
  zdom[ndom].windplane[1].lmn[0] = zdom[ndom].windplane[0].lmn[1] = 0;
  zdom[ndom].windplane[1].lmn[2] = 1;




  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = 0.3 * zdom[ndom].corona_rmin;
    if (zdom[ndom].corona_scale_height < zdom[ndom].corona_zmax)
      zdom[ndom].zlog_scale = 0.3 * zdom[ndom].corona_scale_height;
    else if (zdom[ndom].mdim > 0)
    {
      zdom[ndom].zlog_scale = zdom[ndom].corona_zmax / zdom[ndom].mdim;
    }
    else
    {
      Error ("corona: Cannot define z coordinates unless zdom[ndom].mdim is defined. Aborting\n");
      exit (0);
    }
  }


  zdom[ndom].wind_rho_min = zdom[ndom].corona_rmin;
  zdom[ndom].wind_rho_max = zdom[ndom].corona_rmax;
  zdom[ndom].wind_thetamin = 0;
  zdom[ndom].wind_thetamax = 0;

  return (0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double corona_velocity(x,v) calulates the v the wind at a position r
	x
Arguments:		
	ndom 			The domain number in which the corona is descibed
	double x[]		the position where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	

		
Notes:
	v is set to the Keplerian velocity at this radius

History:
 	00sept	ksl	Coded as part of effort to add a stellar wind option to python
	04aug	ksl	52 -- Modified to return xyz velocity in all
			cases.
	08aug	ksl	Modified for use with multiple domans
 
**************************************************************/

double
corona_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double rho, speed;
  double xtest[3];

  rho = sqrt (x[0] * x[0] + x[1] * x[1]);
  if (rho > 0.0)
    speed = sqrt (G * geo.mstar / rho);
  else
    speed = sqrt (G * geo.mstar / geo.rstar);

  v[0] = -zdom[ndom].corona_vel_frac * speed;
  v[2] = 0.0;
  v[1] = speed;

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


  return (speed);

}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double corona_rho(x) calculates the density of a corona at a position x
Arguments:		
	ndom 		The domain where the corona is descibed
	double x[]		the postion where for the which one desires the denisty
Returns:
	The density at x is returned in gram/cm**3
	
Description:

	At present the function is simply and exponential distribution with a constant
        scale height
		
Notes:

History:
 	00sept	ksl	Coded as part of effort to understand the FUSE data on U Gem.
 
**************************************************************/

double
corona_rho (ndom, x)
     double x[];
{
  double rho;
//OLD not used. 160705 ksl  double tref, t;
//OLD not used. 160705 ksl  double gref, g;
//OLD not used ksl double zscale;

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    Error ("corona_rho: Quitting. Need to think more about coronal model more with vertically extended disk\n");
    exit (0);
  }
//OLD not used 160705 ksl  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
//OLD not used 160705 ksl  t = teff (tref, x[0] / geo.rstar);

//OLD not used 160705 ksl  gref = gdisk (geo.mstar, geo.disk_mdot, geo.rstar);
//OLD not used 160705 ksl  g = geff (gref, x[0] / geo.rstar);

//OLD not used - ksl: zscale = BOLTZMANN * t / (MPROT * g);

  rho = zdom[ndom].corona_base_density * exp (-(x[2]) / zdom[ndom].corona_scale_height);

  if (rho < 1.e-10)
    rho = 1.e-10;               // A floor to the density appears to be needed for some of the
  // ionization calculations

  rho = rho / rho2nh;

  return (rho);
}
