
/***********************************************************/
/** @file  yso.c
 * @Author ksl
 * @date   January, 2018
 *
 * @brief  Simulate a young stellar object
 *
 * These are all of the routines necessary to define a YSO wind, which comprises
 * a KWD wind in some places and a stellar wind elsewhere.  The model is for
 * a stellar wind within the inner windcone.
 *
 * 	XXX The yso implemetentation works but could be reimplemented
 * 	with domains.
 * 
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** @name      get_yso_wind_params
 * @brief      gets input data that are necessary for a YSO
 * 	comprising a region that resembles a KWD wind and a region outside
 * 	this that contins a speherical wind
 *
 * @param [in] int  ndom   The domain number
 * @return     Always returns 0
 *
 * Basically this routine just reads the inputs from the two separate
 * 	wind models, since this is a combination of the two.
 *
 * ###Notes###
 *
 * The model has a spherical wind on the interior.  The windcones
 * 	are set to include the region from the polar axis to the outside
 * 	edge of the wind.  (Note that when filling out the grid regions
 * 	outside the wind may have non-zero values of quantities, such
 * 	as density.  This allows for interpolation to occur. But there
 * 	will not be (should not be) any scattering or whatever in this
 * 	region since that wind cones define the regions where photons
 * 	are scattered.  
 * 
 *
 **********************************************************/

int
get_yso_wind_params (ndom)
     int ndom;
{

/* The approach to get the input parameters is to call both input parameter routines
one after the other*/

  Log ("The yso model has a two component wind a spherical stellar wind, and\n a KWD wind\n");
  get_stellar_wind_params (ndom);
  get_knigge_wind_params (ndom);


  /* Next lines added by SS Sep 04. Changed the wind shape so that the boundary touches the outer 
     corner of the disk rather than the intersection of the disk edge with the xy-plane. */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    zdom[ndom].wind_thetamax = atan (geo.diskrad / (((zdom[ndom].kn_dratio * geo.rstar) + zdisk (geo.diskrad))));
  }


  /* The change in the boundary of the wind (as corner of disk -- see above) 
     means that wind_rho_max nees to be redefined so that it is used correctly
     to compute the boundary of the wind elsewhere. */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    zdom[ndom].rho_max = geo.diskrad - (zdisk (geo.diskrad) * tan (zdom[ndom].wind_thetamax));
  }

  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = geo.rstar;
    zdom[ndom].zlog_scale = 1e-4 * geo.rstar;
    /*
     * 1605- ksl - When Stuart wrote the yso paper, he hardwired zlong_xscale. When implementing domains, I noticed
     * this.  Generally, we would like variables like zlog_scale to depend on other input parameters, and so I fixed
     * this.  For the model calculated by Stuart, 1e7 below corresponds to a number of 2.6e-5, and but I arbitrarily
     * rounded up to 1e-4
     *
     zdom[ndom].zlog_scale = 1e7;
     */
  }

  return (0);
}



/**********************************************************/
/** @name      yso_velocity
 * @brief      double (x,v) calculates the v in cartesian coordinates
 * 	of a yso wind from a position x in cartesian coordinates.
 *
 * @param [in] int  ndom     The domain number 
 * @param [in] double  x[]   the position where for the which one desires the velocity
 * @param [out] double  v[]   The velocity in cartesian coordinates
 * @return     the amplitude of the velocity                                    
 * 	
 * 	
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
yso_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, rzero;
  double zzz;
  double dd;

/* First step is to determine whether we are inside the windcone of the kwd wind. 
This is done as it was done in kn_velocity by seeing if the kwd streamline hits
the disk.  If it does not, then we calculate the velocity for a stellar wind */

  dd = geo.rstar * zdom[ndom].kn_dratio;        //dd is the distance of the focus point in cm 
  r = sqrt (x[0] * x[0] + x[1] * x[1]); // rho of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));  // Nothing prevents this from being less than R_wd

  if (rzero < geo.rstar)
    zzz = stellar_velocity (ndom, x, v);
  else
    zzz = kn_velocity (ndom, x, v);


  return (zzz);

}



/**********************************************************/
/** @name      yso_rho
 * @brief      double (x) calculates the density of an kn_wind at a position x
 *
 * @param [in] int  ndom   The current domain number
 * @param [in] double  x[]   the position where for the which one desires the density
 * @return     The density at x is returned in gram/cm**3
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
yso_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rzero;
  double dd;
  double rho;

  dd = geo.rstar * zdom[ndom].kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]); //rho coordinate of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));  //rho at the base for this streamline

  if (rzero < geo.rstar)
  {
    rho = stellar_rho (ndom, x);
  }
  else
    rho = kn_rho (ndom, x);

  return (rho);
}
