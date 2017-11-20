#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
    Space Telescope Science Institute

Synopsis:
    These are general routines to read in a model that
    is in polar or rtheta coordinates
    grids

Arguments:		

Returns:
 
Description:	

Notes:

    There a various possibilities for how the
    velocities could be entered.  One possibility
    which is the way the zeus_python models work
    is for the velocity to be given in spherical
    polar coordinates.  

    However, internally, python uses xyz coordianes
    for velocites (as measured in the xz plane),
    and that is the model followed, here.  This also
    makes these routines similar to those used
    in imported cylindrical models.  

    This means that if the user provides a model
    where velocities are in spherical polar coordinates
    then one must translate them to the convention
    here before the model is read in



History:
	17nov   ksl Began coding.  Note that currently
                the routines here are mainly dummy
            routines     

**************************************************************/


# define LINELEN 512
# define NCELLS  512

struct
{
  int ndim, mdim, ncell;
  int i[NDIM_MAX * NDIM_MAX], j[NDIM_MAX * NDIM_MAX], inwind[NDIM_MAX * NDIM_MAX];
  double r[NDIM_MAX * NDIM_MAX], theta[NDIM_MAX * NDIM_MAX];
//OLD  double v_r[NDIM_MAX * NDIM_MAX], v_theta[NDIM_MAX * NDIM_MAX],
//OLD    v_phi[NDIM_MAX * NDIM_MAX];
  double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX],
    v_z[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];

  /* Note that the variables below look to be in xy, coordiantes
   * but they are really r, theta
   */

  double wind_x[NDIM_MAX], wind_z[NDIM_MAX], wind_midx[NDIM_MAX],
    wind_midz[NDIM_MAX];

} xx_polar;





/***********************************************************
    Space Telescope Science Institute

Synopsis:
    Read the an arbitray wind model in polar coordinates


Arguments:		

Returns:
 
Description:	

Notes:

    The basic data we need to read in are

//OLD    r theta v_r v_theta v_phi  rho (and optionally T)
    r theta v_x v_y v_z  rho (and optionally T)

    We assume that all of the variables are centered, that is
    we are not assuming that we are giving rho at the center of
    a cell, but that r and v_r are at the edges of a cell. 
    This is someghing that would presumable be easy to change


History:
	17nov   ksl Began coding                           
**************************************************************/

int
import_polar (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fopen (), *fptr;
  char line[LINELEN];
  int n, icell, jcell, ncell, inwind;
  int jz, jx;
  double delta;
  double q1, q2, q3, q4, q5, q6, q7;


  Log ("Reading a model %s in polar (r,theta) coordinates \n", filename);


  if ((fptr = fopen (filename, "r")) == NULL)
    {
      Error ("import_polar: No such file\n");
      exit (0);
    }


  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
    {
      n =
	sscanf (line, " %d %d %d %le %le %le %le %le %le %le", &icell, &jcell,
		&inwind, &q1, &q2, &q3, &q4, &q5, &q6, &q7);
      if (n < 4)
	{
	  continue;
	}
      else
	{
	  xx_polar.i[ncell] = icell;
	  xx_polar.j[ncell] = jcell;
	  xx_polar.r[ncell] = q1;
	  xx_polar.theta[ncell] = q2;
	  xx_polar.v_x[ncell] = q3;
	  xx_polar.v_y[ncell] = q4;
	  xx_polar.v_z[ncell] = q5;
	  xx_polar.rho[ncell] = q6;
	  if (n > 9)
	    {
	      xx_polar.t[ncell] = q7;
	    }
	  else
	    {
	      xx_polar.t[ncell] = 10000.;
	    }
	  ncell++;

	}
    }

  zdom[ndom].ndim = xx_polar.ndim = icell + 1;
  zdom[ndom].mdim = xx_polar.mdim = jcell + 1;
  xx_polar.ncell = ncell;
  zdom[ndom].ndim2 = zdom[ndom].ndim * zdom[ndom].mdim;
  jz = jx = 0;
  for (n = 0; n < xx_polar.ncell; n++)
    {
      if (xx_polar.i[n] == 0)
	{
	  xx_polar.wind_z[jz] = xx_polar.theta[n];
	  jz++;
	}
      if (xx_polar.j[n] == 0)
	{
	  xx_polar.wind_x[jx] = xx_polar.r[n];
	  jx++;
	}
    }

  for (n = 0; n < jz - 1; n++)
    {
      xx_polar.wind_midz[n] =
	0.5 * (xx_polar.wind_z[n] + xx_polar.wind_z[n + 1]);
    }

  delta = (xx_polar.wind_z[jz - 1] - xx_polar.wind_z[jz - 2]);
  xx_polar.wind_midz[jz] = xx_polar.wind_z[jz - 1] + 0.5 * delta;
  for (n = 0; n < jx - 1; n++)
    {
      xx_polar.wind_midx[n] =
	0.5 * (xx_polar.wind_x[n] + xx_polar.wind_x[n + 1]);
    }

  delta = (xx_polar.wind_x[n - 1] - xx_polar.wind_x[n - 2]);
  xx_polar.wind_midx[jx] = xx_polar.wind_x[jx - 1] + 0.5 * delta;
  return (0);
}


/* Create the grid f r a rtheta coordiante system..


 * */

int
polar_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{
  int n, nn;
  double theta;
  double rho_max, rho_min, r, rmin, rmax;
  double zmax;

  /* As in the case of other models we assume that the grid has been
   * read in correctly and so now that the WindPtrs have been generated
   * we can put a lot of information in directly
   */

  for (n = 0; n < xx_polar.ncell; n++)
    {
      wind_ij_to_n (ndom, xx_polar.i[n], xx_polar.j[n], &nn);
      w[nn].r = xx_polar.r[n];
      w[nn].theta = theta = xx_polar.theta[n];

      theta /= RADIAN;

      w[nn].x[0] = w[nn].r * sin (theta);
      w[nn].x[1] = 0;
      w[nn].x[2] = w[nn].r * cos (theta);
      w[nn].v[0] = xx_polar.v_x[n];
      w[nn].v[1] = xx_polar.v_y[n];
      w[nn].v[2] = xx_polar.v_z[n];
      w[nn].inwind = xx_polar.inwind[n];

      theta = xx_polar.wind_midz[xx_polar.j[n]] / RADIAN;


      w[nn].xcen[0] = xx_polar.wind_midx[xx_polar.i[n]] * sin (theta);
      w[nn].xcen[1] = 0;
      w[nn].xcen[2] = xx_polar.wind_midx[xx_polar.i[n]] * cos (theta);

      /* JM 1711 -- copy across the inwind variable to the wind pointer */
      w[nn].inwind = xx_polar.inwind[n];

      /* JM 1711 -- you're either in, or you're out. No part in wind cells allowed! 
       *          there is a question here about which choice (of not in wind or all in 
       *                   wind) is most appropriate */
      if (w[nn].inwind == W_PART_INWIND)
	w[nn].inwind = W_NOT_INWIND;
    }

  /* Now add information used in zdom */

  for (n = 0; n < zdom[ndom].ndim; n++)
    {
      zdom[ndom].wind_x[n] = xx_polar.wind_x[n];
    }



  for (n = 0; n < zdom[ndom].mdim; n++)
    {
      zdom[ndom].wind_z[n] = xx_polar.wind_z[n];
    }

  /* Now set up wind boundaries so they are harmless */


  rmax = rho_max = zmax = 0;
  rmin = rho_min = VERY_BIG;
    for (n = 0; n < xx_polar.ncell; n++)
            {
                      wind_ij_to_n (ndom, xx_polar.i[n], xx_polar.j[n], &nn);

      r=length(w[nn].x);


      if (w[nn].inwind >= 0)
	{
	  if (w[nn].x[0] > rho_max)
	    {
	      rho_max = w[nn].x[0];
	    }
	  if (w[nn].x[2] > zmax)
	    {
	      zmax = w[nn].x[2];
	    }
	  if (r > rmax)
	    {
	      rmax = r;
	    }
	}
      else
	{
	  if (rho_min > w[nn].x[0])
	    {
	      rho_min = w[nn].x[0];
	    }
	  if (rmin > r)
	    {
	      rmin = r;
	    }
	}
    }



  zdom[ndom].wind_rho_min = zdom[ndom].rho_min = rho_min;
  zdom[ndom].wind_rho_max = zdom[ndom].rho_max = rho_max;
  zdom[ndom].zmax = zmax;

  zdom[ndom].rmax = rmax;
  zdom[ndom].rmin = rmin;
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;

  return (0);
}


/* The next section calculates velocites.  We follow the hydro approach of
 * getting those velocities from the original grid.  This is really only
 * used for setting up the grid
 */

double
velocity_polar (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed = 0;
  Log ("Cannot make velocities for polar grid from model yet\n");
  return (speed);
}






/* Fill in plasma ptrs with densities.   
 *
 * For this we assume that the densities read in are 
 * given at the * midpoints of the grid
 *
 * 
 * */


double
rho_polar (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  Log ("Cannot make rho for polar grid from model yet\n");
  return (rho);
}
