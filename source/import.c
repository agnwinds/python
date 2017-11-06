#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

/***********************************************************
    Space Telescope Science Institute

Synopsis:
    These are general purpose routines for readin in model
    grids

Arguments:		

Returns:
 
Description:	

Notes:



History:
	17nov   ksl Began coding                           
**************************************************************/


# define LINELEN 512
# define NCELLS  512

/* The next variables have to be external because we need them to be available later on */

struct
{
  int element[NDIM_MAX];
  double r[NDIM_MAX], v[NDIM_MAX], rho[NDIM_MAX], t[NDIM_MAX];
} xx_1d;

struct
{
  int ele_row[NDIM_MAX], ele_col[NDIM_MAX];
  double r[NDIM_MAX * NDIM_MAX], z[NDIM_MAX * NDIM_MAX];
  double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX],
    v_z[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];
} xx_cyl;

struct
{
  int ele_row[NDIM_MAX], ele_col[NDIM_MAX];
  double r[NDIM_MAX * NDIM_MAX], theta[NDIM_MAX * NDIM_MAX];
  double v_r[NDIM_MAX * NDIM_MAX], v_theta[NDIM_MAX * NDIM_MAX],
    v_phi[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];
} xx_polar;


/* End of section defining where the raw input variables are stored */

int
import_wind (ndom)
     int ndom;


{
  char filename[LINELEN];

  rdstr ("File.with.model2read", filename);

  if (zdom[ndom].coord_type == SPHERICAL)
    {
      import_1d (ndom, filename);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      import_cylindrical (ndom, filename);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      import_polar (ndom, filename);
    }
  else
    {
      Error
	("import_wind: Do not know how to import a model of coor_type %d\n",
	 zdom[ndom].coord_type);
      exit (0);
    }



  exit (0);
}




/***********************************************************
    Space Telescope Science Institute

Synopsis:
    Read the an arbitray wind model intended to mimic a stellar
    wind or shell.


Arguments:		

Returns:
 
Description:	

Notes:

    The basic data we need to read in are

    i r v_r rho (and optionally T)

    We assume that all of the variables are centered, that is
    we are not assuming that we are giving rho at the center of
    a cell, but that r and v_r are at the edges of a cell. 
    This is someghing that would presumable be easy to change


History:
	17nov   ksl Began coding                           
**************************************************************/

int
import_1d (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fopen (), *fptr;
  char line[LINELEN];
  int n, icell, ncell;
  double q1, q2, q3, q4;


  Log ("Reading a 1d model %s\n", filename);


  if ((fptr = fopen (filename, "r")) == NULL)
    {
      printf ("Error: No such file\n");
      return (-1);
    }


  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
    {
      n = sscanf (line, " %d %le %le %le %le", &icell, &q1, &q2, &q3, &q4);
      printf ("ok %d\n", n);
      if (n < 4)
	{
	  printf ("Error. Ignore %s \n", line);
	  continue;
	}
      else
	{
	  xx_1d.element[ncell] = icell;
	  xx_1d.r[ncell] = q1;
	  xx_1d.v[ncell] = q2;
	  xx_1d.rho[ncell] = q3;
	  if (n > 4)
	    {
	      xx_1d.t[ncell] = q4;
	    }
	  else
	    {
	      xx_1d.t[ncell] = 10000.;
	    }
	  ncell++;

	}
    }

/* Having read in the data define some initial variables concerning the model. We cannot create
 * the wind grid or other things at this point, because we do not at this point know what
 * wind cells correspnd to whate elements of the grid */

  zdom[ndom].ndim = ncell + 3;




  return (0);
}


/***********************************************************
    Space Telescope Science Institute

Synopsis:
    Read the an arbitray wind model in cylindrical
    coordinates


Arguments:		

Returns:
 
Description:	

Notes:

    The basic data we need to read in are

    i, j, r z  v_x v_y v_z rho (and optionally T)

    where v_x,_v_y,v_z are the velocities in the x,z plane

    We assume that all of the variables are centered, that is
    we are not assuming that we are giving rho at the center of
    a cell, but that r and v_r are at the edges of a cell. 
    This is someghing that would presumable be easy to change


History:
	17nov   ksl Began coding                           
**************************************************************/

int
import_cylindrical (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fopen (), *fptr;
  char line[LINELEN];
  int n, icell, jcell, ncell;
  double q1, q2, q3, q4, q5, q6, q7;



  Log ("Reading a model in cylindrical coordinates %s\n", filename);



  if ((fptr = fopen (filename, "r")) == NULL)
    {
      printf ("Error: No such file\n");
      return (-1);
    }



  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
    {
      n =
	sscanf (line, " %d %d %le %le %le %le %le %le %le", &icell, &jcell,
		&q1, &q2, &q3, &q4, &q5, &q6, &q7);
      printf ("ok %d\n", n);
      if (n < 4)
	{
	  printf ("Error. Ignore %s \n", line);
	  continue;
	}
      else
	{
	  xx_cyl.ele_row[ncell] = icell;
	  xx_cyl.ele_col[ncell] = jcell;
	  xx_cyl.r[ncell] = q1;
	  xx_cyl.z[ncell] = q2;
	  xx_cyl.v_x[ncell] = q3;
	  xx_cyl.v_y[ncell] = q4;
	  xx_cyl.v_z[ncell] = q5;
	  xx_cyl.rho[ncell] = q6;
	  if (n > 9)
	    {
	      xx_cyl.t[ncell] = q7;
	    }
	  else
	    {
	      xx_cyl.t[ncell] = 10000.;
	    }
	  ncell++;

	}
    }


/* Having read in the data define some initial variables concerning the model. We cannot create
 * the wind grid or other things at this point, because we do not at this point know what
 * wind cells correspnd to whate elements of the grid */

  zdom[ndom].ndim = ncell + 3;

  return (0);
}


/***********************************************************
    Space Telescope Science Institute

Synopsis:
    Read the an arbitray wind model intended in polar coordinates


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
      printf ("Error: No such file\n");
      return (-1);
    }


  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
    {
      n =
	sscanf (line, " %d %d %le %le %le %le %le %le %le", &icell, &jcell,
		&q1, &q2, &q3, &q4, &q5, &q6, &q7);
      printf ("ok %d\n", n);
      if (n < 4)
	{
	  printf ("Error. Ignore %s \n", line);
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

/* Having read in the data define some initial variables concerning the model. We cannot create
 * the wind grid or other things at this point, because we do not at this point know what
 * wind cells correspnd to whate elements of the grid */

  zdom[ndom].ndim = ncell + 3;

  return (0);
}


/* The next section contains routines to make the grids for imported models */

int
spherical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{

  int j, n, ndim;

  ndim = zdom[ndom].ndim - 3;
  for (j = 0; j < ndim; j++)
    {
      n = j + zdom[ndom].nstart;
    w[n].r = xx_1d.r[j];}

  w[n].r = 1.01 * w[n - 1].r;
  w[n + 1].r = 1.02 * w[n - 1].r;
  w[n + 2].r = 1.03 * w[n - 1].r;

  ndim = zdom[ndom].ndim;	// Go now to the edge including buffer cells
  for (j = 0; j < ndim; j++)
    {
      n = j + zdom[ndom].nstart;
      w[n].x[1] = w[n].xcen[1] = 0.0;
      w[n].x[0] = w[n].x[2] = w[n].r * sin (PI / 4.);
      w[n].xcen[0] = w[n].xcen[2] = w[n].rcen * sin (PI / 4.);

    }




  return (0);
}
