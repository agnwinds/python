#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



/* A basis is defined such that if x is a 3 vector as expressed an unprimed cartesian coordinate
   frame, and if y is the same vector in some rotated frame, then
   x[i] = a[i][j] y[j]
   Basis is defined in python.h
 */

/* randvec returns the 3 vector "a" whose direction will be random and 
   whose length will be r.


   Notes: Routine was incorrect for a long time.  Fixed by ck  June 1998

   The current routine, as corrected by CK, creates a random vector
   in spherical coordinates, noting that an element of solid angle i
   s sin(theta) d(theta) d(phi).  But sin(theta)d(theta) is simply 
   -d(cos(theta)), and so a random vector is created by generating
   a random number between 0 and 2 PI representing phi and between
   -1 and 1 representing cos(phi).   

 */

int
randvec (a, r)
     double a[], r;
{

  double costheta, sintheta, phi, sinphi, cosphi;

//  phi = 2. * PI * (rand () / MAXRAND); //DONE
  phi = 2. * PI *gsl_rng_get(rng)/randmax;
  sinphi = sin (phi);
  cosphi = cos (phi);
//  costheta = 2. * (rand () / MAXRAND) - 1.; //DONE
  costheta = 2. * (gsl_rng_get(rng)/randmax) -1.;
  sintheta = sqrt (1. - costheta * costheta);
  a[0] = r * cosphi * sintheta;
  a[1] = r * sinphi * sintheta;
  a[2] = r * costheta;

  return (0);

}

/* Create a photon direction "lmn" in the hemisphere with the vector "north" pointing to the "north
   pole" pole of the hemispere in the case where the photon is originating in a photosphere.
Another way of saying this is that north is the normal to surface at the point
at which the photon originates.  

   The photon directions will be distributed according to the Eddinton approximation

History:
	02jan	ksl	Add jumps array and modified call to pdf_gen_from_func so to taylor the
			pdf_array to include more points near 90 degrees.
*/

double zzz[] = { 0.0, 0.0, 1.0 };


int init_vcos = 0;

int
randvcos (lmn, north)
     double lmn[], north[];
{
  double x[3];                  /* the photon direction in the rotated frame */
  double l, m, n;               /* the individual direction cosines in the rotated frame */
  double q, jumps[5];
  struct basis nbasis;
  int echeck, cdf_gen_from_func ();
  int create_basis (), project_from ();
  double vcos (), cdf_get_rand ();
  double phi;

  /* XXX - It seems unlikely that jumps are necessary with the current verisonof 
   * pdf_gen_gen_from_func, which does not resample the orginal grid.  But this
   * should be checked.  171012 ksl 
   */

  if (init_vcos == 0)
  {
    jumps[0] = 0.01745;
    jumps[1] = 0.03490;
    jumps[2] = 0.05230;
    jumps[3] = 0.06976;
    jumps[4] = 0.08716;

    if ((echeck = cdf_gen_from_func (&cdf_vcos, &vcos, 0., 1., 5, jumps)) != 0)
    {
      Error ("Randvcos: return from cdf_gen_from_func %d\n", echeck);;
    }
    init_vcos = 1;
  }


  n = cdf_get_rand (&cdf_vcos);
  q = sqrt (1. - n * n);

// The is the correct approach to generating a uniform azimuthal distribution

//  phi = 2. * PI * (rand () / MAXRAND); //DONE
  phi = 2. * PI * gsl_rng_get(rng)/randmax;
  l = q * cos (phi);
  m = q * sin (phi);

/* So at this point we have the direction cosines in a frame in which
the z axis is along the normal to the surface.  We must now put the
direction in the cartesian frame.  If north is in the +-z direction
this is simple. Otherwise one must do a coordinate rotation. */

  if (north[0] == 0 && north[1] == 0)
  {                             /* Deal with it as a special case */
    lmn[0] = l;
    lmn[1] = m;
    if (north[2] > 0)
      lmn[2] = n;
    else
      lmn[2] = -n;
  }
  else
  {
    create_basis (north, zzz, &nbasis); /* Create a basis with the first axis in 
                                           direction of "north" and the second axis in the 
                                           yz plane */
    x[0] = n;
    x[1] = l;
    x[2] = m;

    project_from (&nbasis, x, lmn);     /* Project the vector back to the standard
                                           frame */
  }
  return (0);

}


/* This is the appropriate distribution for Eddington limb darkining I believe ksl 97july06 */

double
vcos (x)
     double x;
{
  double a, b;
  double z;

  a = 0.5;
  b = 1.5;
  z = x * (a * (1. + b * x));
  return (z);
}
