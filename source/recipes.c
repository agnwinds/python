
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
#include "atomic.h"
#include <math.h>
#include "recipes.h"
#include "log.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>







#define ITMAX 100
#define EPS 3.0e-8

/**********************************************************/
/**
 * @brief      A routine that finds the rooot of  a function f
 *
 *
 * @param [in] *funcd  The function which one wishes a root for
 * @param [in] x1  one end of the interval
 * @param [in] x2  the other end of the interval
 * @param [in] tol   A tolerence for the root
 * @return   The root of the function f                           .
 *
 * @details
 * A Numerical Recipes routine to find the root of a function known to
 * be between x1 and x2
 *
 * ### Notes ###
 *
 **********************************************************/


double
zbrent (func, x1, x2, tol)
     double x1, x2, tol;
     double (*func) (double);   /* ANSI: double (*func)(double); */
{
  int iter;
  double a = x1, b = x2, c, d, e, min1, min2;
  double fa = (*func) (a), fb = (*func) (b), fc, p, q, r, s, tol1, xm;

  c = d = e = 0;                // to avoid -03 warning


  if (fb * fa > 0.0)
  {
    Log ("ZBRENT: Min %e & Max %e must bracket zero, but got %e & %e\n", x1, x2, fa, fb);
  }
  fc = fb;
  for (iter = 1; iter <= ITMAX; iter++)
  {
    if (fb * fc > 0.0)
    {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs (fc) < fabs (fb))
    {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * EPS * fabs (b) + 0.5 * tol;
    xm = 0.5 * (c - b);
    if (fabs (xm) <= tol1 || fb == 0.0)
      return b;
    if (fabs (e) >= tol1 && fabs (fa) > fabs (fb))
    {
      s = fb / fa;
      if (a == c)
      {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0)
        q = -q;
      p = fabs (p);
      min1 = 3.0 * xm * q - fabs (tol1 * q);
      min2 = fabs (e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2))
      {
        e = d;
        d = p / q;
      }
      else
      {
        d = xm;
        e = d;
      }
    }
    else
    {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs (d) > tol1)
      b += d;
    else
      b += (xm > 0.0 ? fabs (tol1) : -fabs (tol1));
    fb = (*func) (b);
  }
  Error ("Maximum number of iterations exceeded in ZBRENT\n");
  return b;
}

#undef ITMAX
#undef EPS




#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER

/******************************
 * The next two routines were written by ksl.  They were not part of
   the recipes programs which I had but I think they are what was intended
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




#include <math.h>
#define MAXIT 100

/**********************************************************/
/**
 * @brief      A routine that mimimizes a function f
 *
 *
 * @param [in] *funcd
 * @param [in] x1
 * @param [in] x2
 * @param [in] xacc
 * @return   The minimum of the function f                           .
 *
 * @details
 * A Numerical Recipes routine to find the zero of an equation.  A function
 * which returns the value of the function and and the derivative at x is required.
 * The root must be bracketed by x1 and x2.
 *
 * ### Notes ###
 *
 * This function is a "fail safe" algorithm and is a hybrid between the bisection
 * and Newton-Raphson rooting finding algorithms. Whenever the N-R method would
 * taken the solution out of bounds, or whenever the N-R algorithm is not reducing
 * the size of the bracketed interval rapidly enough, the bisection method takes
 * over, which is, apparently, the fail-safe property of this algorithm.
 *
 * 15/02/19: EP added an Exit call as realistically if there is no change in
 *           sign between fl and fh, then there is no root to converge towards
 *           in the bracketed interval. Note that in Numerical Recipes in C, when
 *           this error is encounted the routine nerror is called, which leads
 *           to the program exiting.
 *
 **********************************************************/

double
rtsafe (void (*funcd) (double, double *, double *), double x1, double x2, double xacc)
{
  int j;
  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;

  (*funcd) (x1, &fl, &df);
  (*funcd) (x2, &fh, &df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
  {
    Error ("rtsafe: Root must be bracketed in RTSAFE\n");
    Error ("rtsafe: x1 %8.3e f1 %8.3e x2 %8.3e f2 %8.3e \n", x1, fl, x2, fh);
    Error ("rtsafe: expected f1 * f2 <= 0 but got %8.3e", fl * fh);
    Exit (1);                   // We should be exiting here
  }
  if (fl == 0.0)
    return x1;
  if (fh == 0.0)
    return x2;
  if (fl < 0.0)
  {
    xl = x1;
    xh = x2;
  }
  else
  {
    xh = x1;
    xl = x2;
  }
  rts = 0.5 * (x1 + x2);
  dxold = fabs (x2 - x1);
  dx = dxold;
  (*funcd) (rts, &f, &df);
  for (j = 1; j <= MAXIT; j++)
  {
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) || (fabs (2.0 * f) > fabs (dxold * df)))
    {
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts)
        return rts;
    }
    else
    {
      dxold = dx;
      dx = f / df;
      temp = rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }
    if (fabs (dx) < xacc)
      return rts;
    (*funcd) (rts, &f, &df);
    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  Error (" Maximum number of iterations exceeded in RTSAFE\n");
  return 0.0;
}



#undef MAXIT

#define R 0.61803399
#define CC (1.0-R)
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

/**********************************************************/
/**
 * @brief      A routine that mimimizes a function f
 *
 *
 * @param [in] ax
 * @param [in] bx
 * @param [in] cx
 * @param [in] f  the function to minimize
 * @param [in] tol A tolerance
 * @param [out] * xmin - The place where the minimum occors
 * @return   The minimum of the function f                           .
 *
 * @details
 * This is a Numerical Recipes routine.  Fiven a funYtion f,
 * and bracketing triplet of abscissas a,b,c, such that it is known
 * that f(b) < f(a) and f(b)<f(c), golden returns the minimum value
 * of the function, and the value xmin where the minimum occurs
 *
 * ### Notes ###
 *
 **********************************************************/

double
golden (ax, bx, cx, f, tol, xmin)
     double ax, bx, cx, tol, *xmin;
     double (*f) (double);      /* ANSI: double (*f)(double); */
{
  double f0, f1, f2, f3, x0, x1, x2, x3;

  x0 = ax;
  x3 = cx;
  if (fabs (cx - bx) > fabs (bx - ax))
  {
    x1 = bx;
    x2 = bx + CC * (cx - bx);
  }
  else
  {
    x2 = bx;
    x1 = bx - CC * (bx - ax);
  }
  f1 = (*f) (x1);
  f2 = (*f) (x2);
  while (fabs (x3 - x0) > tol * (fabs (x1) + fabs (x2)))
  {
    if (f2 < f1)
    {
    SHFT (x0, x1, x2, R * x1 + CC * x3) SHFT (f0, f1, f2, (*f) (x2))}
    else
    {
    SHFT (x3, x2, x1, R * x2 + CC * x0) SHFT (f3, f2, f1, (*f) (x1))}
  }
  if (f1 < f2)
  {
    *xmin = x1;
    return f1;
  }
  else
  {
    *xmin = x2;
    return f2;
  }
}

#undef CC
#undef R


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
 * internal test to ensure that the user isnt asking for an integral where the function is zero everywhere.
 * The gsl routine doesnt handle this well, and takes ages to (sometimes) return zero.
 *
 * ### Notes ###
 *
 **********************************************************/


double num_int(func, a, b, eps)
    double (*func) (double,void*);	
    double a, b;
    double eps;
	{
	double result;
	double alpha=0.0;
	void *test=NULL;
	double delta;
	int zflag,i;
	size_t  neval;
	
    gsl_function F;
    F.function = func;
    F.params = &alpha;
	zflag=1;
	if (func(a,test)==0.0 && func(b,test)==0.0)
	{
		zflag=0;
		delta=(b-a)/101;
		for (i=0;i<100;i++)
		{
			if (func(a+delta*i,test)!=0) zflag=1;			
		}
	}
	if (zflag==1)
		{	
    	gsl_integration_romberg_workspace * w  = gsl_integration_romberg_alloc (30);
    	gsl_integration_romberg (&F, a, b, 0, eps, &result, &neval,w);
    	gsl_integration_romberg_free (w);
	}
	else
	{
		result=0.0;
	}
	

	
	return(result);
	}

		
	
double zero_find(func, x_lo, x_hi, tol)
    double (*func) (double,void*);	
    double x_lo, x_hi;
    double tol;
	{
	double result;
	double alpha=0.0;
	double r = 0;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	int iter = 0, max_iter = 100;
	int status;
	

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
	      r = gsl_root_fsolver_root (s);
	      x_lo = gsl_root_fsolver_x_lower (s);
	      x_hi = gsl_root_fsolver_x_upper (s);
	      status = gsl_root_test_interval (x_lo, x_hi,
	                                        tol,0);



	    }
	  while (status == GSL_CONTINUE && iter < max_iter);

	  result=(x_lo+x_hi)/2.0;

	return(result);
	}	
	





