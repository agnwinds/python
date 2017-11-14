#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
    Space Telescope Science Institute

Synopsis:
    These are general routines to read in a model that
    is in polar or rthet coordinates
    grids

Arguments:		

Returns:
 
Description:	

Notes:



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
  int ele_row[NDIM_MAX * NDIM_MAX], ele_col[NDIM_MAX * NDIM_MAX];
  double r[NDIM_MAX * NDIM_MAX], theta[NDIM_MAX * NDIM_MAX];
  double v_r[NDIM_MAX * NDIM_MAX], v_theta[NDIM_MAX * NDIM_MAX],
    v_phi[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];
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

    r theta v_r v_theta v_phi  rho (and optionally T)

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
  int n, icell, jcell, ncell;
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
	sscanf (line, " %d %d %le %le %le %le %le %le %le", &icell, &jcell,
		&q1, &q2, &q3, &q4, &q5, &q6, &q7);
      if (n < 4)
	{
	  continue;
	}
      else
	{
	  xx_polar.ele_row[ncell] = icell;
	  xx_polar.ele_col[ncell] = jcell;
	  xx_polar.r[ncell] = q1;
	  xx_polar.theta[ncell] = q2;
	  xx_polar.v_r[ncell] = q3;
	  xx_polar.v_theta[ncell] = q4;
	  xx_polar.v_phi[ncell] = q5;
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

  zdom[ndom].ndim = xx_polar.ndim = icell;
  zdom[ndom].mdim = xx_polar.mdim = jcell;

  zdom[ndom].ndim += 3;
  zdom[ndom].mdim += 3;



  return (0);
}


/* Create the grid for a rtheta coordiante system..


 * */

int
polar_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{

  Log ("Cannot make rtheta grid from model yet\n");
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
