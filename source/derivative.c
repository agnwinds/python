

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


double f (double x, void * params)
{
  (void)(params); /* avoid unused parameter warning */
  return pow (x, 1.5);
}

int
get_derivative (void)
{
  gsl_function F;
  double result, abserr;

  F.function = &f;
  F.params = 0;

  printf ("f(x) = x^(3/2)\n");

  gsl_deriv_central (&F, 2.0, 1e-8, &result, &abserr);
  printf ("x = 2.0\n");
  printf ("f'(x) = %.10f +/- %.10f\n", result, abserr);
  printf ("exact = %.10f\n\n", 1.5 * sqrt(2.0));

  gsl_deriv_forward (&F, 0.0, 1e-8, &result, &abserr);
  printf ("x = 0.0\n");
  printf ("f'(x) = %.10f +/- %.10f\n", result, abserr);
  printf ("exact = %.10f\n", 0.0);

  return 0;
}





/**********************************************************/
/** 
 * @brief      calculate  the velocity gradient at positions in the flow based on
 * the analytic wind models and imported models.  
 *
 * @param [in] int  ndom   The domain of interest
 * @param [in] double  x[]   A position in the domain
 * @param [out] double  v_grad[][3]   The velocity gradient tensor at the position
 * @return   Always returns 0  
 *
 * @details
 * 
 * The routine calls model_velocity multiple times to calculate the velocity
 * gradient tensor at a particular position.
 *
 * ### Notes ###
 *
 * This routine is normally used only during the initialization
 * of the wind
 *
 * Since the routine calls model_velocity directly, and not vwind_xyz
 * it can be called before wmain.v has been populated.

 *
 **********************************************************/

int
xmodel_vgrad (ds_fraction,ndom, x, v_grad)
    double ds_fraction;
     double x[], v_grad[][3];
     int ndom;
{

  double v0[3], v_forward[3],v_reverse[3];
  double dx_forward[3], dx_reverse[3];
  double dv_forward[3], dv_reverse[3];;
  double dv[3];
  double ds;
  int i, j;
  int vsub (), stuff_v ();

  model_velocity (ndom, x, v0);

  /* get a small distance, 1/1000th of the cell distance from origin */
  ds = ds_fraction * length (x);

  for (i = 0; i < 3; i++)
  {
    /* first create a vector dx, which is the position x but moved ds in direction i */
    stuff_v (x, dx_forward);
    dx_forward[i] += ds;
    stuff_v (x, dx_reverse);
    dx_reverse[i] -= ds;

    /* calculate the velocity at position dx */
    model_velocity (ndom, dx_forward, v_forward);
    model_velocity (ndom, dx_reverse, v_reverse);



    observer_to_local_frame_velocity (v_forward, v0, dv_forward);
    observer_to_local_frame_velocity (v_reverse, v0, dv_reverse);


    vsub(dv_forward,dv_reverse,dv);



    for (j = 0; j < 3; j++)
      dv[j] /= (2*ds);


    stuff_v (dv, v_grad[i]);
  }


  return (0);



}

