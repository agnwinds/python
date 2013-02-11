

/* 
   This file was created in Feb 2011.    The purpose is to have a model where we have a single shell of material. This is different from the stellar wind because we do not want the inner surface of the wind to touch the star. This requires tight control of the grid and makes for a very prescriptive model. 
We also need a special grid, which is also stored in this file.

 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Routines to produce a thin shell of material at a given radius. 


Arguments:		

Returns:
 
Description:	
	

Notes:
	


History:
 	11feb	nsh	Coded as part of the effort to put power laws into python. It allows detailed testing.
 		
**************************************************************/


int
get_shell_wind_params ()
{
	double vtemp[3];
	double rhotemp[200];
	double rmin,dr,r[200];
	double postemp[3];
	double speedtemp;
	double cdensity;
	int	i;


  Log ("Creating wind with a single shell for testing purposes\n");


/* In order to make life simple, we are first going to check that we are in spherical coordinates, if not change!! */

  if (geo.coord_type != SPHERICAL)
      {
	Error("For the shell type wind, we should be in spherical coordinates - changing....\n");
	geo.coord_type = 0;
      }

/* This may not be important, but we sould make sure that NDIM is 4... */

  if (geo.ndim != 4)
     {
       Error("For the shell type wind, we take control of the grid, and need NDIM to be the minimum - 4 - changing\n");
       geo.ndim=4;
       geo.mdim=1;
     }
	

  geo.stellar_wind_mdot = 1.e-6;
  geo.wind_rmin = geo.rstar;
  geo.shell_rmin = 2.8e9;
  geo.shell_vmin = 200e5;
  geo.shell_vmax = 3000e5;

  geo.stellar_wind_mdot /= MSOL / YR;
  rddoub ("stellar_wind_mdot(msol/yr)", &geo.stellar_wind_mdot);
  geo.stellar_wind_mdot *= MSOL / YR;

  rddoub ("stellar.wind.radmin(cm)", &geo.wind_rmin);	/*Radius where wind begins */
  if (geo.wind_rmin < geo.rstar)
    {
      Error
	("get_stellar_wind_params: It is unreasonable to have the wind start inside the star!\n");
      Log ("Setting geo.wind_rmin to geo.rstar\n");
      geo.wind_rmin = geo.rstar;
    }
  geo.shell_rmin = geo.wind_rmin;

  rddoub ("shell.wind_v_at_rmin(cm)", &geo.shell_vmin);	/* Velocity at base of the wind */
  rddoub ("shell.wind.v_at_rmax(cm)", &geo.shell_vmax);	/* Final speed of wind in units of escape velocity */
  rddoub ("shell.wind.acceleration_exponent", &geo.shell_beta);	/* Accleration scale exponent */

/* Assign the generic parameters for the wind the generic parameters of the wind */
//  geo.wind_rmin = geo.rstar;
  geo.shell_rmax = geo.wind_rmax = geo.rmax;
  geo.wind_thetamin = 0.0;
  geo.wind_thetamax = 90. / RADIAN;

/* define the the variables that determine the gridding */
  geo.wind_rho_min = 0;
  geo.wind_rho_max = geo.rmax;
  geo.xlog_scale = 0.3 * geo.rstar;
  geo.zlog_scale = 0.3 * geo.rstar;


/* Since this is a diagnostic routine, we will write out some information to check it is doing what we think) */

	dr=(geo.shell_rmax-geo.shell_rmin)/100.0000;
	printf("dr= %e, root2= %10.30e\n",dr,pow(2.0,0.5));
	rmin=geo.shell_rmin-(dr);
	
	for (i=0;i<103;i++)
		{
		r[i]=rmin+i*dr;
	
		postemp[0]=postemp[2]=r[i]/pow(2.0,0.5);
		postemp[1]=0.0;
		speedtemp=shell_velocity(postemp,vtemp);
		rhotemp[i]=shell_rho(postemp)*rho2nh;
//		dvtemp[i]=shell_dv(postemp);
		printf("ring=%i,x=%e,r=%10.30e,speed=%10.20e,density=%10.20e\n",i,r[i]/pow(2.0,0.5),r[i],speedtemp,rhotemp[i]);
		}

	cdensity=0.0;
	for (i=1;i<100;i++)
		{
		cdensity+=((rhotemp[i]+rhotemp[i+1])/2.)*dr;
		}
	printf("Column density of hydrogen=%e\n",cdensity);


  return (0);
}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	double shell_velocity(x,v) calulates the v the wind at a position 
	x (where both x and v in cartesian coordinates)
Arguments:		
	double x[]		the postion where for the which one desires the velocity
Returns:
	double v[]		the calculated velocity
	
	The amplitude of the velocity is returned 
	
Description:	
	The model is a very simple velocity law of the form

	v(r)=V_min + ((V_max-V_min)/(rmax-rmin))*(r-rmin)**beta
	
	The values of the individiual constants should all be part of the structure geo.

	V_min:  			geo.shell_vmin;		velocity at the inner edge of the spherical shell 
	V_max:		geo.shell_vmax;	the velocity at the outer edge of the spherical shell
	Rmin			geo.shell_rmin	       	the inner radius of the shell
	rmax			geo.shell_rmax          the outer radius of the shell
	beta			geo.shell_beta		power law exponent for the velocity law

		
Notes:
	v is set to V_o inside of rstar or w_rscale, even though the should not really exist there!

	Since this is a radial wind, it is not affected by the disk at all. In particular there
	are no modifications to this routine if the disk is vertically extended.

History:
 	1112	ksl	Adapted by nsh from the same routine for a stellat wind, but the velocity law
			is not the same
 
**************************************************************/

double
shell_velocity (x, v)
     double x[], v[];
{
  double r, speed, zzz;
  double length ();

  if ((r = length (x)) == 0.0)
    {
      v[0] = 0.0;
      v[1] = 0.0;
      v[2] = 0.0;
      return (0.0);
    }

  if (r <= geo.rstar || r <= geo.shell_rmin)
    speed = geo.shell_vmin;
  else
    {
      zzz = pow ((r - geo.shell_rmin), geo.shell_beta);
      speed = geo.shell_vmin + ((geo.shell_vmax - geo.shell_vmin)/(pow((geo.shell_rmax - geo.shell_rmin),geo.shell_beta)) * zzz);
    }
  v[0] = speed * x[0] / r;
  v[1] = speed * x[1] / r;
  v[2] = speed * x[2] / r;

  return (speed);

}

/***********************************************************
		Space Telescope Science Institute

 Synopsis:
	double shell_rho(x) calculates the density of an test shell_wind at a position x
Arguments:		
	double x[]	the position where for the which one desires the density
Returns:
	The density at x is returned in gram/cm**3
	
Description:

	rho=mdot/(4PI*r*r*v);	
		
Notes:

History:
 	11Feb	nsh	Modified from stellar.c portions of the code.
 
**************************************************************/

double
shell_rho (x)
     double x[];
{
  double r, rho, v[3];
  double length (), shell_velocity ();
  r = length (x);
  if (r < geo.shell_rmin)
    {
      rho = geo.stellar_wind_mdot / (4. * PI * r * r * geo.shell_vmin);
    }
  else
    {
      rho =
	geo.stellar_wind_mdot / (4. * PI * r * r * shell_velocity (x, v));
    }

  return (rho);
}


/*
shell_vel_grad calculates the velocity gradient tensor at any point in
the flow

The velocity gradient is defined as a 3 x 3 tensor such that

	velgrad[i][j]= dv_i/dx_j

NB: in c the rightmost index changes changes most rapidly so that
	dv[i]= velgrad[i][j] * dx[j]
makes sense.

NB: Making ds too small can cause roundoff and/or precision errors.

        11feb   nsh     very simple modification of stellar_vel_grad to make sure it calls shell_velocity

*/
int
shell_vel_grad (x, velgrad)
     double x[], velgrad[][3];
{
  double v0[3], v1[3];
  double dx[3], dv[3];
  double ds;
  int i, j;
  int vsub (), stuff_v ();

  shell_velocity (x, v0);
//OLD    printf("HELLOOOO IM IN SELL_VEL_GRAD");
  ds = 1.e7;
  for (i = 0; i < 3; i++)
    {
      stuff_v (x, dx);
      dx[i] += ds;
      stellar_velocity (dx, v1);
      vsub (v1, v0, dv);
      for (j = 0; j < 3; j++)
	dv[j] /= ds;
      stuff_v (dv, &velgrad[i][0]);
    }

  return (0);
}
