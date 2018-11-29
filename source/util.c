
/***********************************************************/
/** @file  util.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  These are very simple utilities within Python that do not
 *   fall into a larger catergory like those in vector.
 *
 * ### Notes ###
 *
 * These routines should probably be refactored into other routines
 *
 * It may be useful to rethink them in addition, as they were all
 * written prior to the existence of domains, and domains are
 * just grafted on.
 *
 ***********************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "log.h"
#include "atomic.h"
#include "python.h"







/**********************************************************/
/**
 * @brief Perform linear/logarithmic interpolation of an array
 *
 * @param [in] double  value   The value used for interpoaltion
 * @param [in] double  array[]   An array containing a set of asceding values
 * @param [in] int  npts   The size of the array
 * @param [out] int *  ival  The lower index to use in the interpolation
 * of the array which need to be used to interpolate on
 * @param [out] double *  f   The fraction of the upper point to use in the interpolation
 * @param [in] int  mode  A switch to choose linear(0)  or lograrithmic(1) interpolation
 * @return     Usually returns 0, but returns -1 if the input value is less
 * than the first elememnt in the array, and 1 if it is greater than the
 * last element int he array.  In either of these cases, the fractions are
 * set up only to access one array element
 *
 * @details
 *
 *
 * Typically one has two parallel arrays, one containing a set of values,
 * in the original case frequencies, and another containing some function
 * of those values.  This routine finds the fraction of the function at
 * two point in the data array to interpolate.
 *
 * The values for the input array can be interpolated linear or logarithmically.
 * that is the fraction that is returned are based on the values in the
 * array or the logarithm of them.
 *
 *
 * The routine uses bisection
 *
 *
 *
 * ### Notes ###
 *
 * fraction is a utility written to speed up the search for the bracketing
 * topbase photionization x-sections, but in principle it should work in
 * other situations within python (which at one time at least contributed
 * significantly to the runtime for Python).
 *
 * Today, fraction is called directly in the routines that are used to set up coordinate
 * systems, and indirectly through linterp (below) which should be inspected
 * to see how the logarithmic interpolation is supposed to work.
 *
 * The routine is similar to the numerical * recipes routine locate.
 * There may be a gsl routine as well.  The routine
 * should probably be replaced.
 *
 *
 **********************************************************/

int
fraction (value, array, npts, ival, f, mode)
     double array[];            // The array in we want to search
     int npts, *ival;           // ival is the lower point
     double value;              // The value we want to index
     double *f;                 // The fractional "distance" to the next point in the array
     int mode;                  // 0 = compute in linear space, 1=compute in log space
{
  int imin, imax, ihalf;

  if (value < array[0])
  {
    *ival = 0;
    *f = 0.0;
    return (-1);
  }

  imax = npts - 1;
  if (value > array[imax])
  {
    *ival = npts - 2;
    *f = 1.0;
    return (1);
  }



  imin = 0;

/* In what follows, there is a specific issue as to how to treat
the situation where the value is exactly on an array element. Right
now this is set to try to identify the array element below the one
on which the value sits, and to set the fraction to 1.  This was
to reflect the behavior of the search routine in where_in_grid. */

  while (imax - imin > 1)
  {
    ihalf = (imin + imax) >> 1; // Compute a midpoint >> is a bitwise right shift
    if (value > array[ihalf])
    {
      imin = ihalf;
    }
    else
      imax = ihalf;
  }

// So array[imin] just <= value

  if (mode == 0)
    *f = (value - array[imin]) / (array[imax] - array[imin]);   //linear interpolation
  else if (mode == 1)
    *f = (log (value) - log (array[imin])) / (log (array[imax]) - log (array[imin]));   //log interpolation
  else
  {
    Error ("Fraction - unknown mode %i\n", mode);
    Exit (0);
    return (0);
  }

  *ival = imin;

  return (0);
}







/**********************************************************/
/**
 * @brief      Perform a linear interpolation on two parallel arrays, the fist
 * of which contains a set of values to be interpolated and the second of
 * which has the function at those values
 *
 * @param [in] double  x   A value
 * @param [in] double  xarray[]   The array that is interplated
 * @param [in] double  yarray[]   The array that contains a function of the values in xaray
 * @param [in] int  xdim   The length of the two arrays
 * @param [out] double *  y   The resulting intepolated value
 * @param [in] int  mode   A switch to choose linear(0) or "logarithmic" (1)
 * interpolation
 *
 * @return     The number of the array element that is used for the lower of the two
 * elements taht are intepolated on.
 *
 * @details
 * Given a number x, and an array of x's in xarray, and functional
 * values y = f(x) in yarray and the dimension of xarray and yarray,
 * linterp calculates y, the linearly interpolated value of y(x). It
 * also returns the element in the xarray so that one can consider
 * not doing a search to find the right array element in certain
 * circumstances

 *
 * ### Notes ###
 *
 * For mode 0, the value that is retuned is
 *
 * (1-f)*y[nelem]+f*[nelem+1)
 *
 * For mode 1, the value returned is
 * exp ((1. - f) * log (y[nelem]) + f * log (y[nelem + 1]))
 *
 *
 **********************************************************/

int
linterp (x, xarray, yarray, xdim, y, mode)
     double x;                  // The value that we wish to index i
     double xarray[], yarray[];
     int xdim;
     double *y;
     int mode;                  //0 = linear, 1 = log
{
  int nelem = 0;
  double frac;


  fraction (x, xarray, xdim, &nelem, &frac, mode);

  if (mode == 0)
    *y = (1. - frac) * yarray[nelem] + frac * yarray[nelem + 1];
  else if (mode == 1)
    *y = exp ((1. - frac) * log (yarray[nelem]) + frac * log (yarray[nelem + 1]));
  else
  {
    Error ("linterp - unknown mode %i\n", mode);
    Exit (0);
  }

  return (nelem);

}



int ierr_coord_fraction = 0;

/**********************************************************/
/**
 * @brief      calculate the fractional
 * 	contributions in interpolations of values of quantities,
 * 	such as electron density, by various cells in the
 * 	grid in a coordinate system independent way.
 *
 * @param [in] int  ndom   The domain in where the interpolation will take place
 * @param [in] int  ichoice  interpolate on vertices (0), interpolate on
 * centers (1)>
 * @param [in] double  x[]   the 3-vector position for which you want
 * the fractinal position
 * @param [out] int  ii[]    an array that contains the 1-d element
 * numbers that must be summed
 * @param [out] double  frac[]   the array that contains the fractional
 * contribution of the corresponding
 * element
 * @param [out] int *  nelem   the number of elements that must be
 * summed, nominally 2 for a spherical grid
 * and 4 for a two dimensional grid.  (For
 * a 3 d grid it would be 6, but we have not
 * implemented this.
 *
 * @return
 * 	 1 if  inside the grid
 * 	-2 if outside the grid
 * 	-1 if inside the grid
 *
 *
 * @details
 *
 * ### Notes ###
 * There are numerous times when one wants the value
 * 	of an interpoalted  variable in the wind.  There
 * 	is no easy way to interpolate the variable easily.
 * 	What this routine does is calculate the fractional
 * 	contributions of elements in the array to that
 * 	position.  Then one must sum up the actual variable
 * 	else where
 *
 * 	If positions are outside the grid, coord_fraction
 * 	attempts to give you the value at the edge of the
 * 	grid.
 *
 * 	It's possible that coord_fraction could be used
 * 	to interpolate beyond the edge of the grid where
 * 	a variable is defined, although this is not done
 * 	at present!
 *
 **********************************************************/

int
coord_fraction (ndom, ichoice, x, ii, frac, nelem)
     int ndom;
     int ichoice;
     double x[];
     int ii[];
     double frac[];
     int *nelem;
{
  double r, z;
  double *xx, *zz;
  int ix, iz;
  double dr, dz;
  int n;


  /* Jump to special routine if CYLVAR coords */

  if (zdom[ndom].coord_type == CYLVAR)
  {

    n = cylvar_coord_fraction (ndom, ichoice, x, ii, frac, nelem);
    if (n < 0 && ierr_coord_fraction < 1000)
    {
      Error ("coord_fraction: cylvar_coord fraction returning %d not in grid\n", n);
      ierr_coord_fraction++;
    }
    return (n);

  }

  /* Assign pointers to the xx and zz depending on whether
   * one wants to interpolate on vertex points (0) or
   * midpoints (1)
   */
  if (ichoice == 0)
  {
    xx = zdom[ndom].wind_x;
    zz = zdom[ndom].wind_z;
  }
  else
  {
    xx = zdom[ndom].wind_midx;
    zz = zdom[ndom].wind_midz;
  }

  /* Now convert x to the appropriate coordinate system */
  r = z = 0.0;
  if (zdom[ndom].coord_type == CYLIND)
  {
    r = sqrt (x[0] * x[0] + x[1] * x[1]);
    z = fabs (x[2]);
  }
  else if (zdom[ndom].coord_type == RTHETA)
  {
    r = length (x);
    z = acos (fabs (x[2]) / r) * RADIAN;
  }
  else if (zdom[ndom].coord_type == SPHERICAL)
  {
    r = length (x);
    z = 0;                      // To avoid -O3 warning
  }
  else
  {
    Error ("coord_fraction: Unknown coordinate type %d for doman\n", zdom[ndom].coord_type, ndom);
    Exit (0);
  }

  if (zdom[ndom].coord_type == SPHERICAL)
  {                             /* We are dealing with a 1d system */
    fraction (r, xx, zdom[ndom].ndim, &ix, &dr, 0);     //linear space
    ii[0] = ix;
    frac[0] = (1. - dr);
    ii[1] = ix + 1;
    frac[1] = dr;
    *nelem = 2;
    if (sane_check (dr))
    {
      Error ("coord_frac:sane_check dr=%f for spherical coords. \n", dr);
    }
  }
  else
  {                             /* We are dealing with a 2d system */
    fraction (r, xx, zdom[ndom].ndim, &ix, &dr, 0);
    fraction (z, zz, zdom[ndom].mdim, &iz, &dz, 0);

    ii[0] = ix * zdom[ndom].mdim + iz;
    frac[0] = (1. - dz) * (1. - dr);

    ii[1] = (ix + 1) * zdom[ndom].mdim + iz;
    frac[1] = (1. - dz) * dr;

    ii[2] = ix * zdom[ndom].mdim + iz + 1;
    frac[2] = (dz) * (1. - dr);

    ii[3] = (ix + 1) * zdom[ndom].mdim + iz + 1;
    frac[3] = (dz) * (dr);
    *nelem = 4;

    if (sane_check (dr) || sane_check (dz))
    {
      Error ("coord_frac:sane_check dr=%f dz=%f for 2d coords\n", dr, dz);
    }
  }
  /* At this point i,j are just outside the x position */
  /* Check to see if x is outside the region of the calculation */
  /* Note that this is a very incoplethe check in the sneste that
   * the posision could be out of the grid in other directions */

  if (r > xx[zdom[ndom].ndim - 1])
  {
    return (-2);                /* x is outside grid */
  }
  else if (r < xx[0])
  {
    return (-1);                /*x is inside grid */
  }

  return (1);

}



int ierr_where_in_2dcell = 0;

/**********************************************************/
/**
 * @brief      calculates the fractional
 * 	location of a position in any 2d grid.
 *
 * @param [in] int  ichoice   interpolate on vertices (0) or centers (1)
 * @param [in] double  x[]   the 3-vector position for which you want
 * the fractinal position
 * @param [inout] int  n   The cell number in wmain
 * @param [out] double *  fx   Fractional position in the x (1st) direction
 * @param [out] double *  fz   fracitionl postion in the z (2nd) direction
 * @return     0 if the positiion is within the cell as defined, non-zero otherwise
 *
 * @details
 *
 * The routine provides the fractions need to interpopolate
 * variables, such as ne, in the wind.  It is called when one
 * already knows what cell one is in.  Like many other routines
 * one can intepolate using the boundaries of the wind cell for
 * values which are defined there (like velocity), and values
 * defined at the cell centeres (like densities)
 *
 *
 * ### Notes ###
 * This routine is general.  It does not rely on the
 * predefined arrays that just contain the grid vertices
 * and centers. On the other hand, it does not calculate
 * which grid cell x lies in. It only tells you that x
 * is not in a particular cell if that is the case.
 *
 * A cell is organized as follows:
 *
 * * x01	x11
 * * x00	x10
 *
 * This routine was written to account for cylind_var coordiantes
 *
 **********************************************************/

int
where_in_2dcell (ichoice, x, n, fx, fz)
     int ichoice;
     double x[];
     int n;                     // A known wind cell
     double *fx, *fz;
{
  double *x00, *x01, *x10, *x11;
  double z[3];
  int i;
  int ndom, nstart, nstop, mdim;

  ndom = wmain[n].ndom;
  nstart = zdom[ndom].nstart;
  nstop = zdom[ndom].nstop;
  mdim = zdom[ndom].mdim;


  /* First check that n is reasonable. This semems almost impossible */
  if (n < nstart || n + mdim + 1 >= nstop)
  {
    if (ierr_where_in_2dcell < 100)
    {
      Error ("where_in_2dcell: Unreasonable n %d This should not happen. Please investigate \n", n);
      ierr_where_in_2dcell++;
    }
    return (n);                 // And hope that calling routine knows how to handle this.
  }

  /* Assign the corners ot the region we want to determine
   * the fractional position of
   */

  if (ichoice == 0)
  {
    x00 = wmain[n].x;
    x01 = wmain[n + 1].x;
    x10 = wmain[n + mdim].x;
    x11 = wmain[n + mdim + 1].x;
  }
  else
  {
    x00 = wmain[n].xcen;
    x01 = wmain[n + 1].xcen;
    x10 = wmain[n + mdim].xcen;
    x11 = wmain[n + mdim + 1].xcen;
  }

  /* Rotate the input vector onto the xz plane
   * which is the plane where the grid is defined
   */

  z[0] = sqrt (x[0] * x[0] + x[1] * x[1]);
  z[1] = 0.0;
  z[2] = fabs (x[2]);

  i = bilin (z, x00, x01, x10, x11, fx, fz);

  return (i);


}




/**********************************************************/
/**
 * @brief      Translate a cell number in the wind to the two-d grid
 * position in a specific domain
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] int  n   The 1-d position in the wind domain
 * @param [out] int *  i   The first index for a 2d position in a given domain
 * @param [out] int *  j   The second index for a 2d postion in a given domain
 * @return     i,j position in a given domain
 *
 * @details
 * wmain contains elements for the entire winds.  This routine translates
 * from a wmain index to the two-d index of a specific domain so that one
 * can interpolate.
 *
 * ### Notes ###
 * For 2d, the underlying 1-d grid is organized so that as n increases,
 * one goes up in z, and then steps out in rho., or in theta and then
 * z as the case may be.
 *
 * With domains, this routine has to be used carefully, and a better
 * statement would be that these routines map to and from the element
 * in wmain to the normally 2d element in an individual domain.
 *
 * So this means that if you want to get the ith and jth element of
 * domain 1, then one needs to give nstart+whatever to wind_n_to_ij
 *
 **********************************************************/

int
wind_n_to_ij (ndom, n, i, j)
     int n, *i, *j, ndom;
{
  int n_use;
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    *i = n - zdom[ndom].nstart;
    *j = 0;
  }
  n_use = n - zdom[ndom].nstart;
  *i = n_use / zdom[ndom].mdim;
  *j = n_use - (*i) * zdom[ndom].mdim;
  return (0);
}


/**********************************************************/
/**
 * @brief      Translate from the 2d element for an individual domain to the wind
 * element number
 *
 * @param [in, out] int  ndom   The domain of interest
 * @param [in, out] int  i   The element number for the first (x) dimension
 * @param [in, out] int  j   The element number fro the second (z or theta) dimension
 * @param [in, out] int *  n   The element number in the wind
 * @return     The element number is returned
 *
 * @details
 * This is just a simple translation to find the wind domain element
 * number
 *
 * ### Notes ###
 *
 **********************************************************/

int
wind_ij_to_n (ndom, i, j, n)
     int *n, i, j, ndom;
{
  if (zdom[ndom].coord_type == SPHERICAL)
  {
    Error ("Warning: wind_ij_to_n being called for spherical coordinates %d %d\n", i, j);
    *n = i;
    return (*n);
  }
  *n = zdom[ndom].nstart + i * zdom[ndom].mdim + j;     // MDIM because the array is in z order
  return (*n);
}


/**********************************************************/
/**
 * @brief      Determine the wind domain element number from a postion
 *
 * @param [in] double  x[]   A three vector describing a position
 * @param [out] int *  n   The element number in the wind domain
 * @return     Also retruns the element number in the wind domain
 *
 * @details
 *
 * The routine simply cycles through the various domains until
 * if finds one in where x is in that domain.
 *
 * ### Notes ###
 *
 * @bug  There are no checks for whether x is within the wind.
 * This seems very dangerous.
 *
 **********************************************************/

int
wind_x_to_n (double x[], int *n)
{
  int ndom, i, j;
  for (ndom = 0; ndom < geo.ndomain; ndom++)
  {
    for (i = 0; i < zdom[ndom].ndim - 1; i++)
    {
      if (x[0] > zdom[ndom].wind_x[i] && x[0] < zdom[ndom].wind_x[i + 1])
      {
        for (j = 0; j < zdom[ndom].mdim - 1; j++)
        {
          if (x[2] > zdom[ndom].wind_z[j] && x[2] < zdom[ndom].wind_z[j + 1])
          {
            wind_ij_to_n (ndom, i, j, n);
            return (*n);
          }
        }
      }
    }
  }
  return (*n);
}
