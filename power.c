

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	Compute the abundances baset on the power law approximation
	use by Stuart in one of his early AGN papers.
	
 Arguments:		
	WindPtr w;

Returns:
 
Description:	
	
Notes:

History:
	111126	ksl	This was originally coded by nsh, but it
			was inserted directly into what was intended
			to be the steering routine for all ionization
			calculations.  I have moved it to its own
			routine in python_71
	111227	ksl	Eliminated n as an external variable
			because it really is not one
	111229	ksl	Small modifications made to reflect moving
			nxfreq and xfreq into the geo structure

**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

// external variables set up so zbrent can solve for alpha using sim_appha as integrand.
double sim_numin, sim_numax, sim_meanfreq;	



int
power_abundances (xplasma,mode)
     PlasmaPtr xplasma;
     int mode;
{
  int ireturn;
  double alphamin, alphamax, alphatemp, sim_w_temp, j;
  int n;

// One shot at updating t_e before calculating densities using the SIM correction

/* This call is after a photon flight, so we *should* have access to j and ave freq, and so we can calculate proper values for W and alpha
To avoid problems with solving, we need to find a reasonable range of values within which to search for a solution to eq18. A reasonable guess is that it is around the current value....
*/



/* We loop over all of the bands, the first band is band number 0, and the last is band nxfreq-1 */
/* 71 - 111229 - ksl - Small modification to reflect moving nxfreq, etc in the geo structure */
  for (n = 0; n < geo.nxfreq; n++)
    {

      if (xplasma->nxtot[n] == 0)
	{
	  Error
	    ("power_abundances: no photons in band %d \n",n);
	  sim_numin = xband.f1[0];	/*NSH 1108 Use the lower bound of the lowest band for sumnumin */
	  sim_numax = xband.f2[xband.nbands - 1];	/*NSH 1108 and the upper bound of the upperband for max */
	  sim_meanfreq = xplasma->ave_freq;
	  j = 0;		// If there are no photons then our guess is that there is no flux in this band
	  xplasma->sim_w[n] = 0;	//We also want to make sure that the weight will be zero, this way we make sure there is no contribution to the ionization balance from this frequency.
	}
      else
	{
	  sim_numin = geo.xfreq[n];	/*1108 NSH n is defined in python.c, and says which band of radiation estimators we are interested in using the for power law ionisation calculation */
	  sim_numax = geo.xfreq[n + 1];
	  sim_meanfreq = xplasma->xave_freq[n];
	  j = xplasma->xj[n];
	}

      Log
	("NSH We are about to calculate w and alpha, j=%10.2e, mean_freq=%10.2e, numin=%10.2e, numax=%10.2e, number of photons in band=%i\n",
	 j, sim_meanfreq, sim_numin, sim_numax, xplasma->nxtot[n]);


      /*1108 NSH ?? this could be a problem. At the moment, it relies on sim_alpha being defined at this point. */
      /*We should initialise it somewhere so that it *always* has a value, not just when the power law */

      alphamin = xplasma->sim_alpha[n] - 0.1;
      alphamax = xplasma->sim_alpha[n] + 0.1;

      while (sim_alpha_func (alphamin) * sim_alpha_func (alphamax) > 0.0)
	{
	  alphamin = alphamin - 1.0;
	  alphamax = alphamax + 1.0;
	}


/* We compute temporary values for sim alpha and sim weight. This will allow us to 
 * check that they are sensible before reassigning them */

      alphatemp = zbrent (sim_alpha_func, alphamin, alphamax, 0.00001);

      if (alphatemp > 3.0)
	alphatemp = 3.0;	//110818 nsh check to stop crazy values for alpha causing problems
      if (alphatemp < -3.0)
	alphatemp = -3.0;

/*This next line computes the sim weight using an external function. Note that xplasma->j already 
 * contains the volume of the cell and a factor of 4pi, so the volume sent to sim_w is set to 1 
 * and j has a factor of 4PI reapplied to it. This means that the equation still works in balance. 
 * It may be better to just implement the factor here, rather than bother with an external call.... */

      sim_w_temp = sim_w (j * 4 * PI, 1, 1, alphatemp, sim_numin, sim_numax);

      if (sane_check (sim_w_temp))
	{
	  Error
	    ("New sim parameters unreasonable, using existing parameters. Check number of photons in this cell\n");
	}
      else
	{
	  xplasma->sim_alpha[n] = alphatemp;
	  xplasma->sim_w[n] = sim_w_temp;
	}

    }

  xplasma->dt_e_old = xplasma->dt_e;
  xplasma->dt_e = xplasma->t_e - xplasma->t_e_old;	//Must store this before others
  xplasma->t_e_old = xplasma->t_e;
  xplasma->t_r_old = xplasma->t_r;
  xplasma->lum_rad_old = xplasma->lum_rad;


  Log ("NSH in this cell, we have %e AGN photons and %e disk photons\n",
       xplasma->ntot_agn, xplasma->ntot_disk);

  ireturn = one_shot (xplasma, mode);


/* Convergence check */
  convergence (xplasma);

  return (ireturn);

}



/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 	integrand used with zbrent to calulate alpha from
	a mean frequency and a minimum and maximum frequecy
	for a band
	
 Arguments:		

Returns:
 
Description:	
	
Notes:

History:

**************************************************************/

double
sim_alpha_func (alpha)
     double alpha;
{
  double answer;
  answer =
    ((alpha + 1.) / (alpha + 2.)) * ((pow (sim_numax, (alpha + 2.)) -
				      pow (sim_numin,
					   (alpha + 2.))) / (pow (sim_numax,
								  (alpha +
								   1.)) -
							     pow (sim_numin,
								  (alpha +
								   1.))));
  answer = answer - sim_meanfreq;
//      printf("NSH alpha=%.3f,f1=%10.2e,f2=%10.2e,meanfreq=%10.2e,ans=%.3f\n",alpha,sim_numin,sim_numax,sim_meanfreq,answer);
  return (answer);
}
