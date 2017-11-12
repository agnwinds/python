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
  int ndim, mdim, ncell;
  int ele_row[NDIM_MAX * NDIM_MAX], ele_col[NDIM_MAX * NDIM_MAX],
    inwind[NDIM_MAX * NDIM_MAX];
  double x[NDIM_MAX * NDIM_MAX], z[NDIM_MAX * NDIM_MAX];
  double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX],
    v_z[NDIM_MAX * NDIM_MAX];
  double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];

  double wind_x[NDIM_MAX], wind_z[NDIM_MAX], wind_midx[NDIM_MAX],
    wind_midz[NDIM_MAX];
} xx_cyl;

struct
{
  int ndim, mdim, ncell;
  int ele_row[NDIM_MAX * NDIM_MAX], ele_col[NDIM_MAX * NDIM_MAX];
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
  zdom[ndom].wind_rho_min = zdom[ndom].rho_min = 0;
  zdom[ndom].rmin = xx_1d.r[0];
  zdom[ndom].wind_rho_max = zdom[ndom].zmax = zdom[ndom].rho_max =
    zdom[ndom].rmax = xx_1d.r[ncell - 1];
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;







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
    and where

    i is the column number  (Thus i corresponds to ndim)
    j is the row number     (and z coresponds to mdim)

    We assume that all of the variables are centered, that is
    we are not assuming that we are giving rho at the center of
    a cell, but that r and v_r are at the edges of a cell. 
    This is someghing that would presumable be easy to change

    Note that we assume that the data are being read in in the
    same order as printed out by windsave2table, that is that
    the first "column" is read in, and then the second "column".




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
  int n, icell, jcell, ncell, inwind;
  double q1, q2, q3, q4, q5, q6, q7;
  double rho_min, rho_max, zmax;
  double r, rmin, rmax, x[3];



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
	sscanf (line, " %d %d %d %le %le %le %le %le %le %le", &icell, &jcell,
		&inwind, &q1, &q2, &q3, &q4, &q5, &q6, &q7);
      if (n < 4)
	{
	  printf ("Error. Ignore %s \n", line);
	  continue;
	}
      else
	{
	  xx_cyl.ele_row[ncell] = icell;
	  xx_cyl.ele_col[ncell] = jcell;
	  xx_cyl.inwind[ncell] = inwind;
	  xx_cyl.x[ncell] = q1;
	  xx_cyl.z[ncell] = q2;
	  xx_cyl.v_x[ncell] = q3;
	  xx_cyl.v_y[ncell] = q4;
	  xx_cyl.v_z[ncell] = q5;
	  xx_cyl.rho[ncell] = q6;
	  if (n > 10)
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

  zdom[ndom].ndim = xx_cyl.ndim = icell + 1;
  zdom[ndom].mdim = xx_cyl.mdim = jcell + 1;
  xx_cyl.ncell = ncell;

  Log ("XX Dimensions of read in model: %d %d\n", zdom[ndom].ndim,
       zdom[ndom].mdim);

  zdom[ndom].ndim += 3;
  zdom[ndom].mdim += 3;


  rmax = rho_max = zmax = 0;
  rmin = rho_min = VERY_BIG;
  for (n = 0; n < ncell; n++)

    {
      x[0] = xx_cyl.x[n];
      x[1] = 0;
      x[2] = xx_cyl.z[n];

      r = length (x);

      if (xx_cyl.inwind[n] >= 0)
	{
	  if (xx_cyl.x[n] > rho_max)
	    {
	      rho_max = xx_cyl.x[n];
	    }
	  if (xx_cyl.z[n] > zmax)
	    {
	      zmax = xx_cyl.z[n];
	    }
	  if (r > rmax)
	    {
	      rmax = r;
	    }
	}
      else
	{
	  if (rho_min > xx_cyl.x[n])
	    {
	      rho_min = xx_cyl.x[n];
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
import_make_grid (w, ndom)
     WindPtr w;
     int ndom;
{
  if (zdom[ndom].coord_type == SPHERICAL)
    {
      spherical_make_grid_import (w, ndom);
    }
  else if (zdom[ndom].coord_type == CYLIND)
    {
      cylindrical_make_grid_import (w, ndom);
    }
  else if (zdom[ndom].coord_type == RTHETA)
    {
      polar_make_grid_import (w, ndom);
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


  w[n + 1].r = 1.01 * w[n].r;
  w[n + 2].r = 1.02 * w[n].r;
  w[n + 3].r = 1.03 * w[n].r;


  for (j = 0; j < zdom[ndom].ndim; j++)
    {
      n = j + zdom[ndom].nstart;
      /* Need to define the midpoints of the grid */
      if (j < zdom[ndom].ndim - 1)
	{
	  w[n].rcen = 0.5 * (w[n].r + w[n + 1].r);
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


// struct
// {
//   int ndim, mdim, ncell;
//   int ele_row[NDIM_MAX], ele_col[NDIM_MAX], inwind[NDIM_MAX * NDIM_MAX];
//   double x[NDIM_MAX * NDIM_MAX], z[NDIM_MAX * NDIM_MAX];
//   double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX],
//     v_z[NDIM_MAX * NDIM_MAX];
//   double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];
// } xx_cyl;

int
cylindrical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{
  int n;
  int jz, jx;
  int kz, kx;
  int nn;


/*  XXX This is an attempt to make the grid directly.  It's inconistent,
 *  somewhat with a separate attempt below.  The problem all
 *  has to do with what one does with the edge cells.
 *  */
  for (n=0;n<xx_cyl.ncell;n++){
      wind_ij_to_n(ndom,xx_cyl.ele_row[n],xx_cyl.ele_col[n],&nn);
      w[nn].x[0]=xx_cyl.wind_x[n];
      w[nn].x[1]=0;
      w[nn].x[2]=xx_cyl.wind_z[n];
      w[nn].v[0]=xx_cyl.v_x[n];
      w[nn].v[1]=xx_cyl.v_y[n];
      w[nn].v[2]=xx_cyl.v_z[n];
  }


  jz = jx = 0;
  for (n = 0; n < xx_cyl.ncell; n++)
    {
      if (xx_cyl.ele_row[n] == 0)
	{
	  xx_cyl.wind_z[jz] = xx_cyl.z[n];
	  jz++;
	}
      if (xx_cyl.ele_col[n] == 0)
	{
	  xx_cyl.wind_x[jx] = xx_cyl.x[n];
	  jx++;
	}
    }





  Log ("Gotcha %d %d %d\n", xx_cyl.ncell, jz, jx);

  /* Now we need to add some cells to wind_x and wind_y
   */

  xx_cyl.wind_x[jx] = 1.01 * xx_cyl.wind_x[jx - 1];
  xx_cyl.wind_x[jx + 1] = 1.02 * xx_cyl.wind_x[jx - 1];
  xx_cyl.wind_x[jx + 2] = 1.03 * xx_cyl.wind_x[jx - 1];

  xx_cyl.wind_z[jz] = 1.01 * xx_cyl.wind_z[jz - 1];
  xx_cyl.wind_z[jz + 1] = 1.02 * xx_cyl.wind_z[jz - 1];
  xx_cyl.wind_z[jz + 2] = 1.03 * xx_cyl.wind_z[jz - 1];

  jx += 3;
  jz += 3;


  /* At this point we should be able to fill the wind array with
   * positions 
   * XXX This is my original attempt.  It does deal with the
   * edge cells, does not put the velocities in
   * */

  n = zdom[ndom].nstart;
  for (kx = 0; kx < jx; kx++)
    {
      for (kz = 0; kz < jz; kz++)
	{
//	  w[n].x[0] = xx_cyl.wind_x[kx];
	  w[n].x[1] = 0;
//	  w[n].x[2] = xx_cyl.wind_z[kz];
	  n++;
	}

    }

  Log ("filled %d %d %d\n", n - zdom[ndom].nstart, jz, jx);


  /* Now add information used in zdom */

  for (n = 0; n < jx; n++)
    {
      zdom[ndom].wind_x[n] = xx_cyl.wind_x[n];
      if (n < jx - 1)
	{
	  zdom[ndom].wind_midx[n] =
	    0.5 * (xx_cyl.wind_x[n] + xx_cyl.wind_x[n + 1]);
	}
      else
	{
	  zdom[ndom].wind_midx[n] = 1.005 * xx_cyl.wind_x[n];
	}
    }



  for (n = 0; n < jz; n++)
    {
      zdom[ndom].wind_z[n] = xx_cyl.wind_z[n];
      if (n < jz - 1)
	{
	  zdom[ndom].wind_midz[n] =
	    0.5 * (xx_cyl.wind_z[n] + xx_cyl.wind_z[n + 1]);
	}
      else
	{
	  zdom[ndom].wind_midz[n] = 1.005 * xx_cyl.wind_z[n];
	}
    }

/* Now that we have the grid we can be clever and stuff the velocities in.  Then
 * We will be able to use interpolation routines */


  for (n = 0; n < xx_cyl.ncell; n++)
    {
      wind_ij_to_n (ndom, xx_cyl.ele_row[n], xx_cyl.ele_col[n], &nn);
      w[nn].x[0] = xx_cyl.v_x[n];
      w[nn].x[1] = xx_cyl.v_y[n];
      w[nn].x[2] = xx_cyl.v_z[n];
    }

  /* We have to do something about the velocities of the edge cells */




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

/* velocity for cylindrical coordinates only interpolates.  One has to
 * interpoalte for vgrad */


double
velocity_cylindrical (ndom, x, v)
     int ndom;
     double *x, *v;
{
  int j;
  int nn;
  int nnn[4], nelem;
  double frac[4];
  double speed;
  coord_fraction (ndom, 0, x, nnn, frac, &nelem);
  for (j = 0; j < 3; j++)
    {
      v[j] = 0;
      for (nn = 0; nn < nelem; nn++)
	v[j] += wmain[zdom[ndom].nstart + nnn[nn]].v[j] * frac[nn];
    }

  speed = length (v);
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


/* XXX This is a placeholder for Fill in the volumes.  Hopefully
 * it will not need to be used and Volumes use the standard 
 * routines for each coordinate system since the volumes just depend
 * on quantities if zdom, and can be calculated after these
 * are properly determined.  There may be a need to exclude
 * complicated regions, but this has not been addressed yet
 */





/* Fill in plasma ptrs with densities.   
 *
 * For this we assume that the densities read in are 
 * given at the * midpoints of the grid
 *
 * 
 * */


double
import_rho (ndom, x)
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
  r = length (x);
  n = 0;
  while (r > xx_1d.r[n] && n < xx_1d.ndim)
    {
      n++;
    }

  if (n < xx_1d.ndim)
    {
      rho = xx_1d.rho[n];
    }
  else
    {
      rho = xx_1d.rho[xx_1d.ndim - 1];
    }


  Log ("rho %e \n", rho);
  return (rho);
}

//OLD struct
//OLD {
//OLD   int ndim, mdim, ncell;
//OLD   int ele_row[NDIM_MAX * NDIM_MAX], ele_col[NDIM_MAX * NDIM_MAX],
//OLD     inwind[NDIM_MAX * NDIM_MAX];
//OLD   double x[NDIM_MAX * NDIM_MAX], z[NDIM_MAX * NDIM_MAX];
//OLD   double v_x[NDIM_MAX * NDIM_MAX], v_y[NDIM_MAX * NDIM_MAX],
//OLD     v_z[NDIM_MAX * NDIM_MAX];
//OLD   double rho[NDIM_MAX * NDIM_MAX], t[NDIM_MAX * NDIM_MAX];
//OLD 
//OLD   double wind_x[NDIM_MAX], wind_z[NDIM_MAX], wind_midx[NDIM_MAX],
//OLD     wind_midz[NDIM_MAX];
//OLD } xx_cyl;

/*  This routine should only be called to set up the plasma cells,
 *  and we assume that rho that was imported is the center of the 
 *  plasma cell, so there is no need to interpolate.
 */

double
rho_cylindrical (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r,z;
  int i,j,n;

  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  z=fabs(x[2]);

  i=0;
  while (z>xx_cyl.wind_z[i] && i<xx_cyl.mdim-1) {
      i++;
  }
  j=0;
  while (r>xx_cyl.wind_x[j] && j<xx_cyl.ndim-1) {
      j++;
  }

  n=j*xx_cyl.mdim+i;

  rho=xx_cyl.rho[n];


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
