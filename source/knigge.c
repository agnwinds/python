
/***********************************************************/
/** @file  knigge.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Routines needed to define a KWD wind
 *
 * The routines here are all that are needed to define a KWD-type
 * wind within Python.  Inputs are gathered, and the any values
 * such as the nomalization needed to convert a global to a local
 * mass loss rate.  Functions are provided to calculated the density
 * and velocity of the wind at any point.
 *
 * As written, the routines assume a standard Shakura-Sunyaeve disk
 * temperature distribution.
 *
 * Although a KWD wind is only really defined within wind cones,
 * the routines will return a non-zero velocity and density outside
 * these regions to allow for interpolating at the edges of the wind
 * cones.
 *
 ***********************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/// This is a variable required to do the integration of local mdot in KWD
double kn_lambda;



/**********************************************************/
/**
 * @brief      get the input data that is necessary to define KWD-type model
 * 	description of the wind
 *
 * @param [in] int  ndom   The domain number
 * @return     Alwasys returns 0
 *
 *
 * ###Notes###
 *
 * The routine first obtains wind_mdot and then the various parameters
 * required for the KWD model.   It then integrates the mass loss rate per unit area over
 * the region occupied by the wind to normalize the mass loss rate.
 *
 * @bug There may weel be errors associated with a vertically extended disk that need tobe addressed
 *
 **********************************************************/

int
get_knigge_wind_params (ndom)
     int ndom;
{
  double dmin;
  double disktheta, test;


  Log ("Creating KWD wind in domain %d\n", ndom);

  zdom[ndom].wind_mdot = 0.1 * geo.disk_mdot / (MSOL / YR);
  rddoub ("Wind.mdot(msol/yr)", &zdom[ndom].wind_mdot);
  zdom[ndom].wind_mdot *= MSOL / YR;


  /* Initialize the various parameters to something reasonable */

  zdom[ndom].kn_r_scale = 7e10; /*Accleration length scale for wind */
  zdom[ndom].kn_alpha = 1.5;    /* Accleration scale exponent for wind */
  zdom[ndom].kn_v_infinity = 3; /* Final speed of wind in units of escape velocity */
  zdom[ndom].kn_lambda = 0.0;   /* Mass loss rate exponent */
  zdom[ndom].kn_dratio = dmin = 0.5 * sqrt (geo.diskrad / geo.rstar);   /* Center of collimation in units of the WD
                                                                           radius. The value set here is for the minimum collimation, see KWD95.  The coefficient 0.5 is
                                                                           approximate */
  zdom[ndom].kn_v_zero = 1.0;   /* Parameter which
                                   can use to muliply the initial velocity of wind so that it
                                   is greater or less than the sound speed. This is not part of the standard KWD model */
  zdom[ndom].wind_rho_min = 1;  /* Innner and outer edges of the wind in stellar radii. These
                                   parameters were added to allow one to create models similar to those 
                                   used in the YSO paper (Sim+05)  */
  zdom[ndom].wind_rho_max = geo.diskrad / geo.rstar;


/* There is confusion in various papers concerning whether to use d or d/dmin.  In KWD95, d/dmin was
used but in later papers, e.g KD97 d in WD radii was used.   This is more natural and is used here.
but one should remember that this differs from KWD95.

To repeat, kn_dratio is the distance to the focus point in stellar radii!

*/

  rddoub ("KWD.d(in_units_of_rstar)", &zdom[ndom].kn_dratio);

  Log_silent ("dmin = %f so the ratio d/dmin here is %f  (%.2e %.2e) \n", dmin, zdom[ndom].kn_dratio / dmin, geo.diskrad, geo.rstar);


  rddoub ("KWD.mdot_r_exponent", &zdom[ndom].kn_lambda);        /* Mass loss rate exponent */
  rddoub ("KWD.v_infinity(in_units_of_vescape)", &zdom[ndom].kn_v_infinity);    /* Final speed of wind in units of escape velocity */
  if (zdom[ndom].kn_v_infinity < 0)
  {
    Log ("Since kn_v_infinity is less than zero, will use SV prescription for velocity law.\n Velocity at base remains the sound speed\n");
  }

  rddoub ("KWD.acceleration_length(cm)", &zdom[ndom].kn_r_scale);       /*Accleration length scale for wind */
  rddoub ("KWD.acceleration_exponent", &zdom[ndom].kn_alpha);   /* Accleration scale exponent */
  rddoub ("KWD.v_zero(multiple_of_sound_speed_at_base)", &zdom[ndom].kn_v_zero);

  rddoub ("KWD.rmin(in_units_of_rstar)", &zdom[ndom].wind_rho_min);
  rddoub ("KWD.rmax(in_units_of_rstar)", &zdom[ndom].wind_rho_max);

  zdom[ndom].wind_rho_min *= geo.rstar;
  zdom[ndom].wind_rho_max *= geo.rstar;
  zdom[ndom].wind_thetamin = atan (1. / zdom[ndom].kn_dratio);

/* Somewhat paradoxically diskrad is in cm, while dn_ratio which is really d in KWD95 is
in units of WD radii */

  zdom[ndom].wind_thetamax = atan (geo.diskrad / (zdom[ndom].kn_dratio * geo.rstar));

  /* Next lines added by SS Sep 04. Changed the wind shape so that the boundary touches the outer
     corner of the disk rather than the intersection of the disk edge with the xy-plane. */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    zdom[ndom].wind_thetamax = atan (geo.diskrad / (((zdom[ndom].kn_dratio * geo.rstar) + zdisk (geo.diskrad))));
  }


  zdom[ndom].rmin = zdom[ndom].wind_rho_min;

  /* The change in the boundary of the wind (as corner of disk -- see above)
     means that wind_rho_max nees to be redefined so that it is used correctly
     to compute the boundary of the wind elsewhere. */

  // XXX Next lines for a vertically extended disk supercede definitions above and look wrong
  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    zdom[ndom].wind_rho_max = geo.diskrad - (zdisk (geo.diskrad) * tan (zdom[ndom].wind_thetamax));
  }


  /* if modes.adjust_grid is 1 then we have already adjusted the grid manually */
  if (modes.adjust_grid == 0)
  {
    zdom[ndom].xlog_scale = geo.rstar;
    zdom[ndom].zlog_scale = geo.rstar;
  }


/*Now calculate the normalization factor for the wind*/

  test = geo.rstar;

  /* For non-flat disk some streamlines are missing (SS). */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    disktheta = atan (zdisk (geo.diskrad) / geo.diskrad);
    test =
      zdom[ndom].kn_dratio * geo.rstar * sin (zdom[ndom].wind_thetamin) *
      cos (disktheta) / sin ((PI / 2.) - zdom[ndom].wind_thetamin - disktheta);
  }

  kn_lambda = zdom[ndom].kn_lambda;
  zdom[ndom].mdot_norm = qromb (kn_wind_mdot_integral, test, geo.diskrad, 1e-6);

  return (0);
}

#define kfudge 1.01


/**********************************************************/
/**
 * @brief      Calculate the v at a postion in a KWD wind
 *
 * @param [in] int  ndom   the domain for the KN wind
 * @param [in] double  x[]   a position in cartesiatn coordiantes
 * @param [out] double  v[]   the calculated velocity at that position in cartesian coordinates
 * @return     The amplitude of the velocity is returned
 *
 * This implements equations 17 & 18 of KWD95.
 *
 *
 *
 * ###Notes###
 *
 * For a contant poloidal velocity along a stream line, kn_alpha should
 * be set to 0
 *
 * There is a largely hidden option that allows one to incoporated a SV velocity
 * law into a KWD model.  This is accomplished by settting 	v_inf<0. The code should
 * be examined for details.
 *
 * At the inner edge of the accretion disk the temperature of the disk
 * goes to zero, and so in principle the thermal velocity goes to zero.
 * The kfudge below keeps this from happening. What did Christian actually
 * do here??
 *
 *
 **********************************************************/

double
kn_velocity (ndom, x, v)
     int ndom;
     double x[], v[];
{
  double r, rzero, theta;
  double ldist, zzz, v_escape, vl;
  double dd;
  double vzero;
  struct photon ptest;
  double xtest[3];
  double s;

  DomainPtr one_dom;


  one_dom = &zdom[ndom];

/*dd is the distance of the focus point along the z axis in cm */

  dd = geo.rstar * one_dom->kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]); // rho of the point we have been given

/* rzero is the position where the stream line  hits the xy plane. This could less than R_wd. */

  rzero = r / (1. + fabs (x[2] / dd));
  theta = atan (rzero / dd);

/* ldist is the poloidal distance along the streamline */

  ldist = sqrt ((r - rzero) * (r - rzero) + x[2] * x[2]);

/* Take the thickness of the disk into account if that is necessary.  */

  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    xtest[0] = r;               // Define xtest in the +xz plane
    xtest[1] = 0;
    xtest[2] = fabs (x[2]);
    ptest.x[0] = rzero;         // Define ptest at the intersection of the streamline and x axis
    ptest.x[1] = 0.0;
    ptest.x[2] = EPSILON;
    ptest.lmn[0] = sin (theta); // lmn is along the streamline toward xtest
    ptest.lmn[1] = 0.0;
    ptest.lmn[2] = cos (theta);
    s = ds_to_disk (&ptest, 1);
    move_phot (&ptest, s);      // Now test photon is at disk surface
    vsub (ptest.x, xtest, xtest);
    ldist = length (xtest);
    rzero = length (ptest.x);
  }


/* The if statement, below  limits the poloidal component of the velocity to be no more than
the poloidal velocity at the inner edge of the wind. It is there for continuity reasons */

  if (rzero < geo.rstar)
    v_escape = sqrt (2. * G * geo.mstar / geo.rstar);
  else
    v_escape = sqrt (2. * G * geo.mstar / rzero);

/* Note that vzero for kwd the sound speed dos not depend on any of the kwd parameters */

  vzero = kn_vzero (rzero);

/* In the standard KWD model, the speed at the base of the wind is the
 sound speed.  The next line modifies this by the factor given by kn_v_zero, which
is an input to the code. */

  vzero *= one_dom->kn_v_zero;

/* 578 -- 06oct -- ksl -- The next lines are modified to allow one to create a SV style
velocity law if kn_v_infinity is less than 0 */

  if (one_dom->kn_v_infinity > 0.0)     // The normal case
  {
    zzz = ldist / (ldist + one_dom->kn_r_scale);
    zzz = pow (zzz, one_dom->kn_alpha); // In Knigge's notation this is actually beta
    vl = vzero + (one_dom->kn_v_infinity * v_escape - vzero) * zzz;
  }
  else                          // the SV options
  {
    zzz = pow (ldist / one_dom->kn_r_scale, one_dom->kn_alpha);
    vl = vzero + ((-one_dom->kn_v_infinity) * v_escape - vzero) * zzz / (1. + zzz);
  }


  v[0] = vl * sin (theta);


  if (rzero > geo.rstar)
    v[1] = sqrt (G * geo.mstar * rzero) / r;    // Eqn 8 KWD
  else if (r > 0)
    v[1] = sqrt (G * geo.mstar * geo.rstar) / r;
  else
    v[1] = 0;
  v[2] = vl * cos (theta);
  if (x[2] < 0)
    v[2] *= (-1);

/*

At this point we have calculated the velocity in the xz plane, which
is identical to the statement that we have calculated it in cylindrical coordinates.
We now project back to xyz coordinates if necessary.

In practice, the current version of python currently should never need to carry out
this projection.  However, it is needed, if we were to move to full 3d calculation and for certain
test programs.

*/
  if (x[1] != 0.0)
  {
    project_from_cyl_xyz (x, v, xtest);
    stuff_v (xtest, v);
  }

  return (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
}



/**********************************************************/
/**
 * @brief      double (x) calculates the density of an kn_wind at a position x
 *
 * @param [in] int  ndom   The domain where kn params are stored
 * @param [in] double  x[]   the position which one desires the density
 * @return     The density at x is returned in gram/cm**3
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
kn_rho (ndom, x)
     int ndom;
     double x[];
{
  double r, rzero;
  double dd;
  double vqr, vzero;
  double v[3], rho;
  struct photon ptest;
  double s, theta;

  DomainPtr one_dom;


  one_dom = &zdom[ndom];


  dd = geo.rstar * one_dom->kn_dratio;
  r = sqrt (x[0] * x[0] + x[1] * x[1]); //rho coordinate of the point we have been given
  rzero = r / (1. + fabs (x[2] / dd));  //rho at the base for this streamline

  /* If the disk is thick we need to modify the position of rzero */
  if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    theta = atan (rzero / dd);
    ptest.x[0] = rzero;
    ptest.x[1] = 0.0;
    ptest.x[2] = EPSILON;
    ptest.lmn[0] = cos (theta);
    ptest.lmn[1] = 0.0;
    ptest.lmn[2] = sin (theta);
    s = ds_to_disk (&ptest, 1);
    move_phot (&ptest, s);      // Now test photon is at disk surface
    rzero = sqrt (ptest.x[0] * ptest.x[0] + ptest.x[1] * ptest.x[1]);
  }

  kn_velocity (ndom, x, v);
  vqr = sqrt (v[0] * v[0] + v[2] * v[2]);       //poloidal velocity
  vzero = kn_vzero (rzero);     //polidal velocity at base
  if ((r * vqr) > 0)            // If statement to assure denominator is not zero
  {
    rho = kn_rho_zero (ndom, rzero) * (rzero * rzero * vzero) / (r * r * vqr);
  }
  else
    rho = 0.0;
  return (rho);
}






/**********************************************************/
/**
 * @brief      Calculate the sound speed (launch speed) of the wind at r
 *
 * @param [in] double  r   The position (radius) for caculating the sound speed
 * @return     The sound speed
 *
 * For KWD, the launch speed is the sound velocity
 *
 * ###Notes###
 *
 * Prior to calling this routine the properties of the disk must have been entered
 *
 * This routine assumes a standard Shakura-Sunyaev disk.
 *
 * In our implementation of the KWD model, the velocity at the footpoint
 * can be a mulitiple of the sound speed, but that correction is applied in kn_velocity
 **********************************************************/

double
kn_vzero (r)
     double r;
{
  double tref, tdisk ();
  double t, teff ();
  double ratio, v;
  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;
  if (ratio < kfudge)
    ratio = kfudge;
  t = teff (tref, ratio);
  v = 1.e6 * sqrt (t / 1e4);    // Frank, King & Raine 1985
  return (v);
}





/**********************************************************/
/**
 * @brief      the integrand of KWD model
 * 	for mdot as a function of radius
 *
 * @param [in] double  r   A position (radius) in the disk
 * @return     The value of the integrand
 *
 * The mass loss rate is proportional to T(r)**(4*kn_lambda)
 * as can be seen from KWD95 equation 10.  This functikon
 * has to be integrated over radius to normalize the local
 * mdot to the integrated value.
 *
 * ###Notes###
 *
 * The integration is carrid out in  get_knigge_wind_params
 *
 * There is an extra factor of 2 because the wind emerges
 *        from both sides of the disk
 *
 **********************************************************/

double
kn_wind_mdot_integral (r)
     double r;
{
  double t;
  double x, ratio;
  double tref;


  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;
  if (ratio < kfudge)
    ratio = kfudge;
  t = teff (tref, ratio);
  x = 4. * PI * pow (t, 4. * kn_lambda) * r;
  return (x);
}



/**********************************************************/
/**
 * @brief      Calculate rho at the base of the wind
 *
 * @param [in] int  ndom   The domain number
 * @param [in] double  r   A position in the disk
 * @return     The mass density at the base of the wind
 *
 *
 * ###Notes###
 *
 *
 **********************************************************/

double
kn_rho_zero (ndom, r)
     double r;
     int ndom;
{
  double tref;
  double t;
  double x, ratio;
  double vzero, dd, cosd;

  DomainPtr one_dom;


  one_dom = &zdom[ndom];


  tref = tdisk (geo.mstar, geo.disk_mdot, geo.rstar);
  ratio = r / geo.rstar;

  if (ratio < kfudge)
    ratio = kfudge;

  t = teff (tref, ratio);
  x = one_dom->wind_mdot * pow (t, 4. * one_dom->kn_lambda) / one_dom->mdot_norm;
  vzero = kn_vzero (r);
  dd = geo.rstar * one_dom->kn_dratio;
  cosd = dd / sqrt (r * r + dd * dd);
  x /= (vzero * cosd);
  return (x);
}
