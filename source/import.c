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
  int ndim;
  int element[NDIM_MAX];
  double r[NDIM_MAX], v[NDIM_MAX], rho[NDIM_MAX], t[NDIM_MAX];
} xx_1d;

struct
{
  int ndim, mdim;
  int ele_row[NDIM_MAX], ele_col[NDIM_MAX];
  double r[NDIM_MAX * NDIM_MAX], z[NDIM_MAX * NDIM_MAX];
  double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX],
    v_z[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];
} xx_cyl;

struct
{
  int ndim, mdim;
  int ele_row[NDIM_MAX], ele_col[NDIM_MAX];
  double r[NDIM_MAX * NDIM_MAX], theta[NDIM_MAX * NDIM_MAX];
  double v_r[NDIM_MAX * NDIM_MAX], v_theta[NDIM_MAX * NDIM_MAX],
    v_phi[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];
} xx_polar;


/* End of section defining where the raw input variables are stored */

/* The next routine just reads in the model into the right structure */

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



  return (0);
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
      Error ("import_1d: No such file\n");
      exit (0);
    }


  ncell = 0;
  while (fgets (line, 512, fptr) != NULL)
    {
      n = sscanf (line, " %d %le %le %le %le", &icell, &q1, &q2, &q3, &q4);
      if (n < 4)
	{
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


  xx_1d.ndim = ncell;
  zdom[ndom].ndim = ncell + 3;	// ADD Buffer
  zdom[ndom].mdim = 1;
  zdom[ndom].rho_min=zdom[ndom].rmin=xx_1d.r[0];
  zdom[ndom].zmax=zdom[ndom].rho_max=zdom[ndom].rmax=xx_1d.r[ncell-1];
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax=0.;







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
      Error ("import_cylindrical: No such file\n");
      exit (0);
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

  zdom[ndom].ndim = xx_cyl.ndim = icell;
  zdom[ndom].mdim = xx_cyl.mdim = jcell;

  zdom[ndom].ndim += 3;
  zdom[ndom].mdim += 3;



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


/* The next section contains routines to make the grids for imported models.

   We make some assumptions here.  We assume that every cell is in the wind.
   and we assume that r refers to the inside edge of the cell.

 * */

int
spherical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{

  int j, n;

  for (j = 0; j < xx_1d.ndim; j++)
    {
      n = j + zdom[ndom].nstart;
      w[n].r = xx_1d.r[j];

    }

  /* We have already added a buffer to zdom[ndom].ndim
   * so we can extend the grid
   */


  /* n has not been incremented past the end of input array
   * so we have to add 1 here 
   */

  w[n + 1].r = 1.01 * w[n - 1].r;
  w[n + 2].r = 1.02 * w[n - 1].r;
  w[n + 3].r = 1.03 * w[n - 1].r;


  for (j = 0; j < zdom[ndom].ndim; j++)
    {
      n = j + zdom[ndom].nstart;
      /* Need to define the midpoints of the grid */
      if (j < zdom[ndom].ndim - 1)
	{
	  w[n].rcen = 0.5*(w[n].r + w[n + 1].r);
	}
      else
	{
	  w[n].rcen = w[n].r * 1.005;
	}
      w[n].x[1] = w[n].xcen[1] = 0.0;
      w[n].x[0] = w[n].x[2] = w[n].r * sin (PI / 4.);
      w[n].xcen[0] = w[n].xcen[2] = w[n].rcen * sin (PI / 4.);
    }

  return (0);
}

int
cylindrical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{

  Log ("Cannot make cylindrical grid from model yet\n");

  return (0);
}

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
import_velocity (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed;


  if (zdom[ndom].coord_type == SPHERICAL)
    {
      speed = velocity_1d (ndom, x, v);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      speed = velocity_cylindrical (ndom, x, v);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      speed = velocity_polar (ndom, x, v);
    }
  else
    {
      Error
	("import_velocity: Do not know how to create velocities from model of coor_type %d\n",
	 zdom[ndom].coord_type);
      exit (0);
    }


  return (speed);
}


double
velocity_1d (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed;
  double r;

  int icell;

  r = length (x);

  icell = linterp (r, xx_1d.r, xx_1d.v, xx_1d.ndim, &speed, 0);

  v[0] = x[0] / r * speed;
  v[1] = x[1] / r * speed;
  v[2] = x[2] / r * speed;



  return (speed);
}


double
velocity_cylindrical (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed = 0;

  Log ("Cannot make velocities for cylindrical grid from model yet\n");

  return (speed);
}


double
velocity_polar (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed = 0;


  Log ("Cannot make velocities for polar grid from model yet\n");


  return (speed);
}


int
get_import_wind_params (ndom)
     int ndom;
{
  Log ("get_import_wind_params is currently a NOP\n");
  return (0);
}


/* Fill in the volumes */

//HOLD  int
//HOLD  import_volumes(ndom)

//HOLD       int ndom;
//HOLD  {


//HOLD    if (zdom[ndom].coord_type == SPHERICAL)
//HOLD      {
//HOLD        import_1d_volumes(ndom);
//HOLD      }
//HOLD    else if (zdom[ndom].coord_type == CYLIND)
//HOLD      {
//HOLD        import_cylindrical_volumes (ndom);
//HOLD      }
//HOLD    else if (zdom[ndom].coord_type == RTHETA)
//HOLD      {
//HOLD        import_polar_volumes(ndom);
//HOLD      }
//HOLD    else
//HOLD      {
//HOLD        Error
//HOLD  	("import_rho:  Do not know how to create velocities from model of coor_type %d\n",
//HOLD  	 zdom[ndom].coord_type);
//HOLD        exit (0);
//HOLD      }

//HOLD    return (0);
//HOLD  }

//HOLD  int import_1d_volumes(ndom)
//HOLD      int ndom;
//HOLD  {
//HOLD    Log ("Cannot make volumes for 1d  grid from model yet\n");
//HOLD    return(0);
//HOLD  }



//HOLD  int import_cylindrical_volumes(ndom)
//HOLD      int ndom;
//HOLD  {
//HOLD    Log ("Cannot make volumes for cylindrical  grid from model yet\n");
//HOLD    return(0);
//HOLD  }



//HOLD  int import_polar_volumes(ndom)
//HOLD      int ndom;
//HOLD  {
//HOLD    Log ("Cannot make volumes for polar grid from model yet\n");
//HOLD    return(0);
//HOLD  }




/* Fill in plasma ptrs with deesities.   
 *
 * For this we assume that the densities read in are 
 * given at the * midpoints of the grid
 *
 * 
 * */


double
import_rho(ndom, x)

     int ndom;
     double *x;
{
  double rho;


  if (zdom[ndom].coord_type == SPHERICAL)
    {
      rho = rho_1d (ndom, x);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      rho = rho_cylindrical (ndom, x);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      rho = rho_polar (ndom, x);
    }
  else
    {
      Error
	("import_rho:  Do not know how to create velocities from model of coor_type %d\n",
	 zdom[ndom].coord_type);
      exit (0);
    }


  return (rho);
}


double
rho_1d (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r;
  int n;

  r=length(x);

  n=0;
  while (r<xx_1d.r[n] && n<xx_1d.ndim)
  {
      n++;
  }

  if (n<xx_1d.ndim) 
  {
      rho=xx_1d.rho[n];
  }
  else
  {
      rho=xx_1d.rho[xx_1d.ndim-1];
  }
  

  Log ("rho %e \n",rho);


  return (rho);
}


double
rho_cylindrical (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;


  Log ("Cannot make rho for cylindrical  grid from model yet\n");


  return (rho);
}


double
rho_polar (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;


  Log ("Cannot make rho for polar grid from model yet\n");


  return (rho);
}

