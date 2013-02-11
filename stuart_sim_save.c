#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"
#include "balance_templates.h"

struct photoionization *sim_xver;	//Verner & Ferland description of a photoionization x-section
struct topbase_phot *sim_xtop;	//Topbase description of a photoionization x-section
double sim_te;


#define MIN_FUDGE  1.e-10
#define MAX_FUDGE  10.

int
sim_driver (xplasma)
     PlasmaPtr xplasma;
{
  int nelem, nion, niterate;
  double xne, xnew;
  double newden[NIONS];
  double t_r, nh;
  double t_e, www;

  t_r = xplasma->t_r;
  t_e = xplasma->t_e;
  www = xplasma->w;


  /* Initally assume electron density from the LTE densities */

  xne = xplasma->ne;
  if (xne < DENSITY_MIN)
    {
      Error
	("nebular_concentrations: Very low ionization: ne initially %8.2e\n",
	 xne);
      xne = DENSITY_MIN;
    }

  nh = xplasma->rho * rho2nh;	//LTE -- Not clear needed at this level

  /* Begin iteration loop to find ne */
  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
      for (nelem = 0; nelem < nelements; nelem++)
	{
	  sim_pl (nh, t_r, t_e, www, nelem, xplasma->ne,
			 xplasma->density, xne, newden);

	  /* Re solve for the macro atom populations with the current guess for ne */
	  if (geo.macro_ioniz_mode == 1)
	    {
	      macro_pops (xplasma, xne);
	    }
	}
      for (nion = 0; nion < nions; nion++)
	{
	  /* if the ion is being treated by macro_pops then use the populations just computed */
	  if ((ion[nion].macro_info == 1) && (geo.macro_simple == 0)
	      && (geo.macro_ioniz_mode == 1))
	    {
	      newden[nion] = xplasma->density[nion];
	    }

	  /*Set some floor so future divisions are sensible */
	  if (newden[nion] < DENSITY_MIN)
	    newden[nion] = DENSITY_MIN;
	}
      xnew = get_ne (newden);
      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;

      /* Check to see whether the search for xne has converged and if so exit loop */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      /* else star another iteration of the main loop */
      xne = (xnew + xne) / 2.;	/* Make a new estimate of xne */
      niterate++;
    }
  /* End of main iteration loop */

  if (niterate == MAXITERATIONS)
    {
      Error
	("nebular_concentrations: failed to converge:nh %8.2e www %8.2e t_e %8.2e  t_r %8.2e \n",
	 nh, www, t_e, t_r);
      return (-1);
    }

/* Finally transfer the calculated densities to the real density array */

  xplasma->ne = xnew;
  for (nion = 0; nion < nions; nion++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nion].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nion] = newden[nion];
	}
    }
  return (0);
}








int
sim_pl (xplasma)
     PlasmaPtr xplasma; 
			
{
	double freq1,freq2,freq;
        double weight,t_e;
	int    nion;
	double simfactor[nions];
	double tbtest[10000];
        double ntmin,ntmax,fmin,fmax,dfreq;
	int j,n;
	int    first,last,nelem;
	

  	t_e = xplasma->t_e;
  	weight = xplasma->w;



	freq1=1e14;
	freq2=1e19;		
	freq=1e15;	


	printf ("weight = %e\n",weight);
	printf ("power law exponent = %e\n",geo.alpha_agn);
	printf ("t_e=%e\n",t_e);	

		
  

	for (nion = 0; nion < nions; nion++)	
		{
    		printf ("we know about ion %i which is atom %i in state %i with a density of %e\n",nion,ion[nion].z,ion[nion].istate,xplasma->density[nion]);
		if (ion[nion].phot_info == -1)
			{
			if (ion[nion].istate == (ion[nion].z)+1)
				printf ("Fully ionised so obviously no data\n");
			else
				printf ("No Topbase or VFKY information, so setting sim factor to 1.0\n");
			simfactor[nion]=1.0;
			}		
		if (ion[nion].phot_info == 0)
			{		

			printf ("And that ion has a VFKY ionisation record associated with it \n");
			simfactor[nion]= xinteg_sim(t_e, freq1,freq2, nion);
			printf ("Returned value for sim factor is %e\n",simfactor[nion]);


			}
		if (ion[nion].phot_info == 1)
			{			
			printf ("And that ion has a topbase ionisation record associated with it \n");
			printf ("It has a startingtopbase of %i, and a total of %i topbase levels\n",ion[nion].ntop_first,ion[nion].ntop); 
			simfactor[nion]= xinteg_sim(t_e, freq1,freq2, nion);

			printf ("Returned value for sim factor is %e\n",simfactor[nion]);
			}
		}

 /*	  for (j=0;j<=9999;j++)
	     tbtest[j]=0.0;
     ntmin = ion[6].ntop_first;
      ntmax = ntmin + ion[6].ntop;
	fmin=1e16;
	fmax=1e16;

	for (n = ntmin; n < ntmax; n++)

    	{
          sim_xtop = &phot_top[n];
			// loop over relevent Topbase photoionzation x-sections
	printf("fmin=%e,fmax=%e\n",sim_xtop->freq[0],sim_xtop->freq[sim_xtop->np - 1]);
	  if (sim_xtop->freq[0] < fmin)
		fmin=sim_xtop->freq[0];
	  if (sim_xtop->freq[sim_xtop->np - 1] > fmax)
		fmax=sim_xtop->freq[sim_xtop->np - 1];	// Argues that this should be part of structure

         }

	dfreq=(fmax-fmin)/10000;
	for (n = ntmin; n < ntmax; n++)
    	{
	  printf("%i\n",n);
          sim_xtop = &phot_top[n];
	  for (j=0;j<=9999;j++)
             {
             freq=fmin+(dfreq*j);
	    tbtest[j]+=sigma_phot_topbase(sim_xtop,freq);
             }
	}
	f = fopen ("topbase", "r");

	for (j=0;j<=9999;j++)
		{
		freq=fmin+(dfreq*j);
	//	printf ("Freq=%e crosssection=%e\n",freq,tbtest[j]);
  		}
*/
//    Ready to apply the modifications to densities

	printf ("Now we are going to apply the factors. We know about %i elements\n",nelements);

	for (nelem=0; nelem<nelements ; nelem++)
		{	
		first = ele[nelem].firstion;	/*first and last identify the postion in the array */
		last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */
		printf ("Element %i has %i ions, starting with %i\n",nelem,ele[nelem].nions,ele[nelem].firstion);		
		for (nion = first +1; nion < last; nion++)
			{
//			fudge2 =ground_frac[nion - 1].frac[ilow] +interpfrac * (ground_frac[nion - 1].frac[ihi] -ground_frac[nion - 1].frac[ilow]);
			printf ("Lower level is ion %i which is atom %i in state %i with a density of %e\n",nion,ion[nion-1].z,ion[nion-1].istate,xplasma->density[nion-1]);
			printf ("Upper level is ion %i which is atom %i in state %i with a density of %e\n",nion,ion[nion].z,ion[nion].istate,xplasma->density[nion]);
			printf ("Current density ratio is %e\n",xplasma->density[nion]/xplasma->density[nion-1]);
			printf ("And our sim ratio is %e\n",simfactor[nion-1]);
//			printf ("And the lucy mazzali fudge factor is %e",fudge2);
			}
		}

/*
  while (density[first] < 1.1 * DENSITY_MIN)
    {
      newden[first] = DENSITY_MIN;
      first++;
    }
  while (density[last - 1] < 1.1 * DENSITY_MIN)
    {
      newden[last - 1] = DENSITY_MIN;
      last--;
    }

  sum = newden[first] = 1.;
  for (nion = first + 1; nion < last; nion++)
    {
      numerator = newden[nion - 1] * fudge * (ne) * density[nion];
      denominator = density[nion - 1] * xne;
      q = numerator / denominator;

      /* now apply the Mazzali and  Lucy fudge, i.e. zeta + w (1-zeta)
         first find fraction of recombinations going directly to
         the ground state for this temperature and ion */
//	fudge2 =
//	ground_frac[nion - 1].frac[ilow] +
//	interpfrac * (ground_frac[nion - 1].frac[ihi] -
//		      ground_frac[nion - 1].frac[ilow]);
      /*Since nion-1 points at state i-1 (saha: n_i/n_i-1) we want ground_frac[nion-1].
         Note that we NEVER access the last ion for any element in that way
         which is just as well since you can't have recombinations INTO
         a fully ionized state -- The corresponding lines in recombination.out 
         are just dummies to make reading the data easier */
//      fudge2 = fudge2 + www * (1.0 - fudge2);
//      newden[nion] = fudge2 * q;
//      sum += newden[nion];
// This is the equation being calculated
//              sum+=newden[nion]=newden[nion-1]*fudge*(*ne)*density[nion]/density[nion-1]/xne;
//    }
//  a = nh * ele[nelem].abun / sum;
//  for (nion = first; nion < last; nion++)
//    newden[nion] *= a;




return(0);
}



double
xinteg_sim (t, f1, f2, nion)
     double t;			// The temperature at which to integrate the cross sections
     double f1, f2;		// The frequencies overwhich to integrate the cross sections
     int nion;			// The ion for which the integrated cross section is calculated
{
  int n;
  double sim_num,sim_denom,sim_frac;    // The upper and lower part of the fraction used to recalculate the ion fractions
  double fthresh, fmax;
  double den_config ();
  double sigma_phot (), sigma_phot_topbase ();
  int ntmin, ntmax;		// These are the Topbase photo-ionization levels that are used
  int nvmin, nvmax;		// These are the limits on the Verland x-sections
  double qromb ();
  double fb_planck (), verner_planck ();


	printf ("Got here with t=%e,f1=%e,f2=%e,nion=%i\n",t, f1, f2, nion);
  if (-1 < nion && nion < nions)	//Get emissivity for this specific ion_number
   {
     ntmin = ion[nion].ntop_first;
      ntmax = ntmin + ion[nion].ntop;
      nvmin = nion;
      nvmax = nvmin + 1;
    }
  else				// Get the total emissivity
    {
      Error ("integ_fb: %d is unacceptable value of nion\n", nion);
      mytrap ();
      exit (0);
      return (0);
    }

// Put information where it can be used by the integrating function
	sim_te=t;

// Place over all limits on the integration interval if they are very large
/* We need to limit the frequency range to one that is reasonable
if we are going to integrate */
  if (f1 < 3e14)
    f1 = 3e14;			// 10000 Angstroms
  if (f2 > 3e19)
    f2 = 3e19;			// 10 Angstroms
  if (f2 < f1)
    return (0);			// Because either f2 corresponded to something redward of 1000 A or f1 
  // was blueward of 10 Angstroms

  sim_num = 0.0;
  sim_denom = 0.0;

  for (n = ntmin; n < ntmax; n++)
    {				// loop over relevent Topbase photoionzation x-sections
      sim_xtop = &phot_top[n];
      /* Adding an if statement here so that photoionization that's part of a macro atom is 
         not included here (these will be dealt with elsewhere). (SS, Apr04) */
      if (sim_xtop->macro_info == 0 || geo.macro_simple == 1 || geo.rt_mode == 1)	//Macro atom check. (SS)
	{
	  fthresh = sim_xtop->freq[0];
	  fmax = sim_xtop->freq[sim_xtop->np - 1];	// Argues that this should be part of structure
	  if (f1 > fthresh)
	    fthresh = f1;
	  if (f2 < fmax)
	    fmax = f2;

	  // Now calculate the emissivity as long as fmax exceeds xthreshold and there are ions to recombine
	  if (fmax > fthresh)
	printf ("We are going to topbase qromb with fthresh=%e and fmax=%e\n",fthresh,fmax);
	    sim_denom += qromb (tb_planck, fthresh, fmax, 1.e-4);
	    sim_num += qromb (tb_pow, fthresh, fmax, 1.e-4);
	printf ("And we are back in the room with PL=%e and Planck=%e\n",sim_num,sim_denom);
	}
  }
// This completes the calculation of those levels for which we have Topbase x-sections, now do Verner

  for (n = nvmin; n < nvmax; n++)
    {
      if (ion[n].phot_info == 0)
	{			// Only work on ions without Topbase and with Verner
	  sim_xver = &xphot[ion[n].nxphot];
	  fthresh = sim_xver->freq_t;
	  fmax = sim_xver->freq_max;	//So at this point these are the maximal allowable
	  if (f1 > fthresh)
	    fthresh = f1;	//So move fthresh up, if f1 was greater than this
	  if (f2 < fmax)	//Move fmax down if we wanted a narrower range
	    fmax = f2;
//	  // Now integrate only if its in allowable range  && there are ions to recombine
	  if (fmax > fthresh)
	printf ("We are going to verner qromb with fthresh=%e and fmax=%e\n",fthresh,fmax);
	    sim_denom += qromb (verner_planck, fthresh, fmax, 1.e-4);
	    sim_num += qromb (verner_pow, fthresh, fmax, 1.e-4);
	printf ("And we are back in the room with PL=%e and Planck=%e\n",sim_num,sim_denom);
	}
    }

	sim_frac=sim_num/sim_denom;


  return (sim_frac);
}

double
tb_planck(freq)
	double freq;
{
	double answer,bbe;	

	bbe=exp((H*freq)/(BOLTZMANN*sim_te));
	answer=(2.*H*pow(freq,3.))/(pow(C,2));
 	answer*=(1/(bbe-1));
	answer*=sigma_phot_topbase(sim_xtop,freq);

	return (answer);
}

double
verner_planck(freq)
	double freq;
{
	double answer,bbe;	
	bbe=exp((H*freq)/(BOLTZMANN*sim_te));
	answer=(2.*H*pow(freq,3.))/(pow(C,2));
 	answer*=(1/(bbe-1));
	answer*=sigma_phot(sim_xver,freq);

	return (answer);
}

double
tb_pow(freq)
	double freq;
{
	double answer;	
	double pl_weight;



	answer=geo.const_agn*(pow(freq,geo.alpha_agn)); 
	pl_weight = 1.33333*geo.d_agn*geo.d_agn;
	answer/=pl_weight;
	answer*=sigma_phot_topbase(sim_xtop,freq);


	return (answer);
}

double
verner_pow(freq)
	double freq;
{
	double answer;	
	double pl_weight;	

	answer=geo.const_agn*(pow(freq,geo.alpha_agn)); 
	pl_weight = 1.33333*geo.d_agn*geo.d_agn;
	answer/=pl_weight;
	answer*=sigma_phot(sim_xver,freq);


	return (answer);
}













