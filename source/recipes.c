
/***********************************************************/
/** @file  recipes.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  These are collection of Numberical Recipes routines
 * used in Python
 *
 *
 * The routines have been modified to make them ansi compatable
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>

#include "atomic.h"
#include "python.h"
//OLD #include "recipes.h"
#include "log.h"


/******************************
 * The next two routines were written by ksl.  They were not part of
   the recipes programs which I had but I think they are what was intended
   TODO EP: I think these to functions are redundant now - we don't seem to use them
********************************/

double *
vector (i, j)
     int i, j;
{
  double dummy, *d;
  d = calloc (sizeof (dummy), (j - i + 1) + 1);
  return (d);
}

void
free_vector (a, i, j)
     double *a;
     int i, j;
{
  free (a);
}



/**********************************************************/
/**
 * @brief      A wrapper function that carries out numerical integration on a supplied function between limits
 *
 *
 * @param [in] func - a function to integrate, needs to be of the form (double (double,void*))
 * @param [in] a - lower bound
 * @param [in] b - upper bound
 * @param [in] eps - the relative accuracy desired.
 * @return   The integral                          .
 *
 * @details
 * This routine carries out numerical integration of the function func from a to b. It currently
 * uses the gsl implementation of the ROMBERG integration scheme. Initial tests using QAG showed that
 * the type of integrals we do - often with exponential tails - is not well treated. There is an 
 * internal test to ensure that the user is not asking for an integral where the function is zero everywhere.
 * The gsl routine does not handle this well, and takes ages to run and (sometimes) return zero.
 *
 * ### Notes ###
 *
 **********************************************************/


double
num_int (func, a, b, eps)
     double (*func) (double, void *);
     double a, b;
     double eps;
{
  double result, error, result2;
  double alpha = 0.0;
  void *test = NULL;
  double delta;
  int zflag, i;
  int status = 0;
  int status2 = 0;

//OLD  int npoints;
  size_t neval;
  gsl_function F;
  F.function = func;
  F.params = &alpha;
  zflag = 1;
//OLD  npoints = 1000;
  if (func (a, test) == 0.0 && func (b, test) == 0.0)
  {
    zflag = 0;
    delta = (b - a) / 101;
    for (i = 0; i < 100; i++)
    {
      if (func (a + delta * i, test) != 0)
        zflag = 1;
    }
  }
  if (zflag == 1)
  {
    gsl_set_error_handler_off ();       //We need to be able to catch and handle gsl errors 

    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    status = gsl_integration_qags (&F, a, b, 0, eps, 1000, w, &result, &error);
    if (status)
    {
      if (status == GSL_EROUND) //The rounding error has been triggered - try a different integrator
      {
        gsl_integration_workspace_free (w);
        gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc (30);
        status2 = gsl_integration_romberg (&F, a, b, 0, eps, &result2, &neval, w);
        gsl_integration_romberg_free (w);
        if (status2)
        {
          Error ("num_init: some kind of error in romberg and qags integration\n");
        }
      }
    }
    else
    {
      gsl_integration_workspace_free (w);
    }
  }
  else
  {
    result = 0.0;
  }

  return (result);
}

/**********************************************************/
/**
* @brief      A wrapper function that finds the root of a function between two limits
*
*
* @param [in] *func - the function we want to find the root of
* @param [in] double x_lo - lower bound
* @param [in] double x_hi - upper bound
* @param [in] double bol - the relative accuracy desired.
* @param [out] int * ierr An error return, TRUE if an error
* @return   The location between a and b of the zero point                         .
*
* @details
* This routine finds the root  of the function func from x_lo to x_hi. It currently
* uses the gsl implementation of the BRENT root finding scheme. This replaced the 
* zbrentnumerical recipes, and the call is intended to be identical
*
* ### Notes ###
*
*  The function needs to be of the form 
*
* double (* function) (double x, void * params)
*
* See the documentation for gsl_function, but historically passed
* the information needed via external variables and so typically
* nothing is passed in the parameters.  This is something that
* ultimately should be changed.
*
**********************************************************/


double
zero_find (func, x1, x2, tol, ierr)
     double (*func) (double, void *);
     double x1, x2;
     double tol;
     int *ierr;
{
  double result;
  double x_below, x_above;
  double alpha = 0.0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  int iter = 0, max_iter = 100;
  int status;
  double x_lo, x_hi;

// If necessary, reorder the inputs to the gsl function

  if (x1 < x2)
  {
    x_lo = x1;
    x_hi = x2;
  }
  else
  {
    x_lo = x2;
    x_hi = x1;
  }



  *ierr = FALSE;
  // Check that the interval is bracketed

  x_below = func (x_lo, (void *) &alpha);
  x_above = func (x_hi, (void *) &alpha);

  if (x_below * x_above > 0.0)
  {
    Error ("zero_find: function not bracketed x_lo %e -> %e, x_hi %e -> %e\n", x_lo, x_below, x_hi, x_above);
    *ierr = TRUE;
  }


  gsl_function F;
  F.function = func;
  F.params = &alpha;



  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);



  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, tol, 0);



  }
  while (status == GSL_CONTINUE && iter < max_iter);

  if (status != GSL_SUCCESS)
  {
    x_below = func (x_lo, (void *) &alpha);
    x_above = func (x_hi, (void *) &alpha);
    Error ("zero_find failed: %d of %d brackets x_lo %e -> %e, x_hi %e -> %e\n", iter, max_iter, x_lo, x_below, x_hi, x_above);
    *ierr = TRUE;
  }

  result = (x_lo + x_hi) / 2.0;

  gsl_root_fsolver_free (s);

  return (result);
}


/**********************************************************/
/**
* @brief      Find the mimimum value of a function in an interval
*
*
* @param [in] a - lower bound
* @param [in] b - guess for the minumum
* @param [in] c - upper bound
* @param [in] f  the function to minimize
* @param [in] tol A tolerance
* @param [out] * xmin - The place where the minimum occors
* @return   The minimum of the function f                           .
*
* @details
* This is a wrapper routine which uses gsl routines to 
* find the minimum of a function
*
* ### Notes ###
*
* At of  212222, this routine is only use by the roche lobe 
* realated routines
* 
**********************************************************/



double
find_function_minimum (a, m, b, func, tol, xmin)
     double (*func) (double, void *);
     double a, m, b;
     double tol, *xmin;

{
  int status = 0;
  void *test = NULL;

  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;

  gsl_function F;
  F.function = func;
  F.params = 0;

  T = gsl_min_fminimizer_brent;

  s = gsl_min_fminimizer_alloc (T);
  gsl_set_error_handler_off ();
  status = gsl_min_fminimizer_set (s, &F, m, a, b);
  if (status)
  {
    if (status == GSL_EINVAL)   //THere is no minimum 
    {
      return fmin (func (a, test), func (b, test));     //keep old behaviour
    }
  }

  do
  {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
    m = gsl_min_fminimizer_x_minimum (s);
    a = gsl_min_fminimizer_x_lower (s);
    b = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (a, b, 0.0, tol);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  *xmin = m;

  gsl_min_fminimizer_free (s);

  return func (m, test);
}
