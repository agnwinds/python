#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"
#include "balance_templates.h"

struct photoionization *sim_xver;	//Verner & Ferland description of a photoionization x-section
struct topbase_phot *sim_xtop;	//Topbase description of a photoionization x-section
double sim_te;                  //This is a storage variable for the current electron temperature so it is available for qromb calls
double weight;                  //This is a storage variable for the current PL geometric weight so it can be used in qromb
double xsim_alpha;               //This is a storage variable so the value of alpha (PL spectral index) can be used in qromb
double xsim_w;                  //Storage variable for the sim W factor (either computed from photons going thruogh a cell, or before photon generation it is the power law constant x radiative weight / 4pi */

double fudge_store[300];           //this is a store so we can see how the sim fudge factor varies
double num_store[300];             // store for the numbnerator of the sim factor (PL)
double denom_store[300];           // store for the denominator of the sim factor (BB)
double ion_store[300];

#define SAHA 4.82907e15		/* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
#define MAXITERATIONS	200
#define FRACTIONAL_ERROR 0.03
#define THETAMAX	 1e4
#define MIN_TEMP         100.

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
  www = weight = xplasma->w;

  xsim_alpha = xplasma->sim_alpha;  //Set the shared variable to the alpha for the cell.
  xsim_w = xplasma->sim_w;  //Set the shared variable to the W for the cell.
	printf ("We are in sim, t_e=%f, t_r=%f, sim_alpha=%f, sim_w=%e,\n",t_e,t_r,xsim_alpha,xsim_w);


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


	printf("Starting off with ne=%e and nh=%e\n",xne,nh);

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
	printf ("current estimate of ne after loop %i=%e\n",niterate,xnew);
      if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;

      /* Check to see whether the search for xne has converged and if so exit loop */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	break;

      /* else start another iteration of the main loop */
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
//	printf("Density of ion %i following sim correction is %e\n",nion,newden[nion]);
	
	}
    }
  return (0);
}


int
sim_pl (nh, t_r, t_e, www, nelem, ne, density, xne, newden)
     double nh, t_r, t_e, www;
     int nelem;
     double ne, density[], xne, newden[];
{
  double fudge, dummy, interpfrac;
  double fudge2, q;
  double sum, a;
  int ilow, ihi;
  int first, last, nion;
  double numerator, denominator;
  double f1,f2,max_ratio,current_ratio;
    FILE *fopen(),*fudgefile,*numfile,*denomfile,*ionfile;

	f1=1.23e15;
	f2=1.21e19;

	weight=www;

  if (t_e > MIN_TEMP)
    {
	fudge=1.;

      /* now get the right place in the ground_frac tables  CK */
      dummy = t_e / TMIN - 1.;
      ilow = dummy;		/* have now truncated to integer below */
      ihi = ilow + 1;		/*these are the indeces bracketing the true value */
      interpfrac = (dummy - ilow);	/*this is the interpolation fraction */
      if (ilow < 0)
	{
	  ilow = 0;
	  ihi = 0;
	  interpfrac = 1.;
	}
      if (ihi > 19)
	{
	  ilow = 19;
	  ihi = 19;
	  interpfrac = 1.;
	}

    }

  else
    {
      Error_silent ("sim_pl: t_r too low www %8.2e t_e %8.2e  t_r %8.2e \n",
	     www, t_e, t_r);
      fudge = 0.0;
      interpfrac = 0.0;
      ihi = ilow = 0;
    }
  if (fudge < MIN_FUDGE || MAX_FUDGE < fudge)
    {
      Error_silent
	("sim_pl: fudge>10 www %8.2e t_e %8.2e  t_r %8.2e \n",
	 www, t_e, t_r);
    }

  /* Initialization of fudges complete open lots of files to log the sim factor*/
		fudgefile=fopen ("fudge_summary.out", "a");
		numfile  = fopen ("num_summary.out", "a");
		denomfile  = fopen ("denom_summary.out", "a");
	        ionfile = fopen ("ion_summary.out","a");
		fprintf (fudgefile,"Te= %e, Element %s",t_e,ele[nelem].name);
		fprintf (numfile,"Te= %e, Element %s",t_e,ele[nelem].name);
		fprintf (denomfile,"Te= %e, Element %s",t_e,ele[nelem].name);
		fprintf (ionfile,"Te= %e, Element %s",t_e,ele[nelem].name);
  first = ele[nelem].firstion;	/*first and last identify the postion in the array */
  last = first + (ele[nelem].nions);	/*  So for H which has 2 ions, H1 and H2, first will generally
					   be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */

//   printf ("WE are working on element %i (%s), which has %i ions, starting with %i\n",nelem,ele[nelem].name,ele[nelem].nions,ele[nelem].firstion);
//   printf ("Ion %i is of element %i, and Ion %i is of element %i\n",ele[nelem].firstion,ion[ele[nelem].firstion].z,ele[nelem].firstion+ele[nelem].nions-1,ion[ele[nelem].firstion+ele[nelem].nions-1].z);

	for (nion = ele[nelem].firstion; nion<(ele[nelem].firstion+ ele[nelem].nions); nion++)
		{
	ion_store[nion]=density[nion];
		}



  while (density[first] < 1.1 * DENSITY_MIN)
    {
     newden[first] = DENSITY_MIN;
      fudge_store[first]=-1.;
      num_store[first]=-1.;
      denom_store[first]=-1.;
      first++;

    }
  while (density[last-1] < 1.1 * DENSITY_MIN)
    {
      newden[last - 1] = DENSITY_MIN;
      fudge_store[last-1]=-1;
      num_store[last-1]=-1;
      denom_store[last-1]=-1;
      last--;
    }    

//   printf ("We have decided that we need to loop over ions %i (element %i) to %i (element %i)\n",
//	first,ion[first].z,last-1,ion[last-1].z);

   max_ratio=((ele[nelem].abun*nh)/DENSITY_MIN)/10.0;  //This is the maximum ratio between two ions of the same element
  sum = newden[first] = 1e-200;   //We need a very large dynamic range to make the sim correction factor work. 
	fudge_store[first]=0.0;
  for (nion = first + 1; nion < last; nion++)   
    {

//	printf("Going to xinteg_sim with nion=%i, d1=%e and d2=%e\n",nion,density[nion],density[nion-1]);
	current_ratio=(density[nion]/density[nion-1]);
      fudge = (xinteg_sim (t_e, f1, f2, nion-1,max_ratio,current_ratio));  //get the sim correction factor, we use the lower ionisartion cross section
	fudge_store[nion]=fudge;
	ion_store[nion]=density[nion];
//	printf("Sim factor for nion=%i,with current densities d1=%e and d2=%e is %e\n",nion,density[nion],density[nion-1],fudge);
    numerator = newden[nion - 1] * fudge * (ne) * density[nion];
      denominator = density[nion - 1] * xne;
      q = numerator / denominator;
//	printf ("q factor for ion %i of element %i is %e. Num=%e, 1=%e, 2=%e, 3=%e, 4=%e, Denom=%e, 1=%e, 2=%e\n",nion,ion[nion].z,q,numerator, newden[nion - 1] , fudge ,(ne) , density[nion]  ,          denominator,density[nion - 1] ,xne);
      /* find fraction of recombinations going directly to
         the ground state for this temperature and ion */
      fudge2 =
	ground_frac[nion - 1].frac[ilow] +
	interpfrac * (ground_frac[nion - 1].frac[ihi] -
		      ground_frac[nion - 1].frac[ilow]); 



      /*Since nion-1 points at state i-1 (saha: n_i/n_i-1) we want ground_frac[nion-1].
         Note that we NEVER access the last ion for any element in that way
         which is just as well since you can't have recombinations INTO
         a fully ionized state -- The corresponding lines in recombination.out 
         are just dummies to make reading the data easier */
 //     fudge2 = fudge2 + www * (1.0 - fudge2);
      newden[nion] = fudge2 * q;  //we multiply by the recombination fraction going directly to ground state
      sum += newden[nion];
// This is the equation being calculated
//              sum+=newden[nion]=newden[nion-1]*fudge*(*ne)*density[nion]/density[nion-1]/xne;
    }
	for (nion = ele[nelem].firstion; nion<(ele[nelem].firstion+ ele[nelem].nions); nion++)
		{
		fprintf (fudgefile,",%e",fudge_store[nion]);
		fprintf (numfile,",%e",num_store[nion]);
		fprintf (denomfile,",%e",denom_store[nion]);
		fprintf (ionfile,",%e",ion_store[nion]);
		}
		fprintf (fudgefile,"\n");
		fprintf (numfile,"\n");
		fprintf (denomfile,"\n");
		fprintf (ionfile,"\n");

		fclose (fudgefile);
		fclose (numfile);
		fclose (denomfile);
		fclose (ionfile);

  a = nh * ele[nelem].abun / sum;   //get the factor to ensure we dont end up with more than the abundance of the element
  for (nion = first; nion < last; nion++)
    newden[nion] *= a;



  return (0);
}









double
xinteg_sim (t, f1, f2, nion,max_ratio,current_ratio)
     double t;			// The temperature at which to integrate the cross sections
     double f1, f2,max_ratio,current_ratio;		// The frequencies overwhich to integrate the cross sections
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

//	printf ("Got here with t=%e,f1=%e,f2=%e,nion=%i,maxratio=%e,current_ratio=%e\n",t, f1, f2, nion,maxratio,current_ratio);


//	printf("We are in xinteg_sim with ion%i (element%i in state%i). This has a phot info of %i\n",nion,ion[nion].z,ion[nion].istate,ion[nion].phot_info);

  if (-1 < nion && nion < nions)	//Get cross section for this specific ion_number
   {
      ntmin = ion[nion].ntop_first;
      ntmax = ntmin + ion[nion].ntop;
      nvmin = nion;
      nvmax = nvmin + 1;
//;     printf ("For ion %i, the first topbase level is %i, and the first verner level is %i\n",nion,ion[nion].ntop_first,nvmin);
    }
  else				// Get the total emissivity
    {
      Error ("xinteg_sim: %d is unacceptable value of nion\n", nion);
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
    f1 = 3e14;			// 
  if (f2 > 3e19)
    f2 = 3e19;			// 
  if (f2 < f1)
    return (0);			// Because either f2 corresponded to something redward of 1000 A or f1 
  // was blueward of 10 Angstroms

  sim_num = 0.0;
  sim_denom = 0.0;

	if (ion[nion].phot_info == 1)
		{
  		n = ntmin;                    // we only want to use the ground state topbase ionisation cross section
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


	  		if (fmax > fthresh)
				{
//	printf ("We are going to topbase qromb for ion %i with fthresh=%e and fmax=%e\n",n,fthresh,fmax);
				sim_num = qromb (tb_pow, fthresh, fmax, 1.e-4);
				sim_denom = qromb (tb_planck, fthresh, fmax, 1.e-4);
				}


//	printf ("And we are back in the room with PL=%e and Planck=%e\n",sim_num,sim_denom);
//	printf ("for this ion, its current density is %e, and the minimum density is %e",density,DENSITY_MIN);
			}
   		}
// This completes the calculation of those levels for which we have Topbase x-sections, now do Verner

    
    
     else if (ion[nvmin].phot_info == 0)
	{			// Only work on ions without Topbase and with Verner
	  n = nvmin;  //just the ground state ioinzation fraction.
	  sim_xver = &xphot[ion[n].nxphot];
	  fthresh = sim_xver->freq_t;
//	printf("XXXXXXXXX fthresh=%e\n",fthresh);
	  fmax = sim_xver->freq_max;	//So at this point these are the maximal allowable
	  if (f1 > fthresh)
	    fthresh = f1;	//So move fthresh up, if f1 was greater than this
	  if (f2 < fmax)	//Move fmax down if we wanted a narrower range
	    fmax = f2;
//	printf("We want to do verner, with fmax=%e and fthresh=%e",fmax,fthresh);
//	  // Now integrate only if its in allowable range  && there are ions to recombine
	  if (fmax > fthresh)
//	printf ("We are going to verner qromb with ion %i, fthresh=%e and fmax=%e\n",n,fthresh,fmax);
		{
	    sim_denom = qromb (verner_planck, fthresh, fmax, 1.e-4);
	    sim_num = qromb (verner_pow, fthresh, fmax, 1.e-4);
		}
//	printf ("And we are back in the room with PL=%e and Planck=%e\n",sim_num,sim_denom);
	}
    else    
	{
//	printf ("No ionisation data for ion %i,setting sim_frac to 1\n",nion);
	sim_frac=1.;
	return (sim_frac);
	} 
//	printf ("Alert!!!!!!!!  test=%e, current ratio=%e, max_ratio=%e\n",((sim_num*current_ratio)/max_ratio),current_ratio,max_ratio);
	if (sim_denom < ((sim_num*current_ratio)/max_ratio))
		{ 
		printf ("Sim denom of %e is too low, resetting to %e\n",sim_denom,((sim_num*current_ratio)/max_ratio));
		sim_denom = ((sim_num*current_ratio)/max_ratio);
//		sim_denom = 1e100;
		}

	num_store[nion+1] = sim_num;
	denom_store[nion+1] = sim_denom;

	sim_frac=sim_num/sim_denom;

//	printf("Returning sim_frac=%e\n",sim_frac);

  return (sim_frac);
}

/* tb_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
    ionisation cross section is significant. This is the function for ions with a topbase cross section NSH 16/12/10 */

double
tb_planck(freq)
	double freq;
{
	double answer,bbe;	

	bbe=exp((H*freq)/(BOLTZMANN*sim_te));
	answer=(2.*H*pow(freq,3.))/(pow(C,2));
 	answer*=(1/(bbe-1));
//	answer*=weight;
	answer*=sigma_phot_topbase(sim_xtop,freq);
	answer/=freq;

	return (answer);
}


/* verner_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
    ionisation cross section is significant. This is the function for ions with a verner cross section NSH 16/12/10 */

double
verner_planck(freq)
	double freq;
{
	double answer,bbe;	
	bbe=exp((H*freq)/(BOLTZMANN*sim_te));
	answer=(2.*H*pow(freq,3.))/(pow(C,2));
 	answer*=(1/(bbe-1));
//	answer*=weight;
	answer*=sigma_phot(sim_xver,freq);
	answer/=freq;

	return (answer);
}


/* tb_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation cross section
   to allow the numerator of the sim correction factor to be calulated for ions with topbase cross sections.*/

double
tb_pow(freq)
	double freq;
{
	double answer;
//	printf("we have alpha=%f and w=%e\n",xsim_alpha,xsim_w);
//	answer=geo.const_agn*(pow(freq,(xsim_alpha-1.0))); 
//	answer*=weight;//divide by the weight, should normally just be the surface area of a sphere.
//	answer/=(4.0*PI); //divide by 4pi to get a mean intensity ber solid angle.
	answer=xsim_w*(pow(freq,(xsim_alpha-1.0)));
	answer*=sigma_phot_topbase(sim_xtop,freq); // and finally multiply by the cross section.


	return (answer);
}

/* verner_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation cross section
   to allow the numerator of the sim correction factor to be calulated for ions with verner cross sections.*/

double
verner_pow(freq)
	double freq;
{
	double answer;		
//	answer=geo.const_agn*(pow(freq,(xsim_alpha-1.0)));
//	answer*=weight;   //divide by the weight, should normally just be the surface area of a sphere.
//	answer/=(4.0*PI); //divide by 4pi to get a mean intensity ber solid angle.
	answer=xsim_w*(pow(freq,(xsim_alpha-1.0)));
	answer*=sigma_phot(sim_xver,freq); // and finally multiply by the cross section.
	return (answer);
}


// NSH 7/2/11
// This little routine solves equation 18 of Sim et at (2008) to get an estimate of the spectral index of a cell
// Equation 18 seems to be a monotonically increasing function in alpha, so this should work (famous last words)

double
sim_alphasolve(ratans,numin,numax)
	double ratans,numin,numax;
{
double alphamin,alphamax,alphamid;
double error,ratmid,alpha;


// Set crazy values of range of alpha - if this screws up I might need to put in a check here
alphamin=-10;
alphamax=10;
alpha=-0.2;
// Set a value of error that gets the loop running
error=10.;
	printf("Solving with ratans=%e, numin=%e, numax=%e\n",ratans,numin,numax);
	printf("Just out of interest, for alpha=-0.2, we get %e\n",((alpha+1)/(alpha+2))*((pow(numax,(alpha+2))-pow(numin,(alpha+2)))/(pow(numax,(alpha+1))-pow(numin,(alpha+1)))));
while (error>0.001)
	{
// Work out a mid point
	alphamid=(alphamax+alphamin)/2.0;
// What En(2)/En(1) is implied by this alpha
	ratmid=((alphamid+1)/(alphamid+2))*((pow(numax,(alphamid+2))-pow(numin,(alphamid+2)))/(pow(numax,(alphamid+1))-pow(numin,(alphamid+1))));
// If actual answer is below ratmid, our solution is in the lower half
	if (ratans<ratmid)
		alphamax=alphamid;
	else    
// If actual answer is above ratmin, our solution is in the upper half
		alphamin=alphamid;
	error=fabs(alphamax-alphamin);
	printf ("Searching for answer, alphamin=%f, alphamax=%f, error=%f, current answer=%e \n",alphamin,alphamax,error,ratmid);
	}	
// If we get here, we must have converged.
return (alphamid);
}


//NSH 7/2/11 - A little function to evaluate equation 19 in Sim et at (2008)
// This works out the effective radiative weight in a cell given a value of alpha and other bits.

double
sim_w(en1,v,dt,alpha,numin,numax)
	double en1; //Sum of energies of photons times frequency times path length in cell
	double v; //volume of cell
	double dt; //Effective time interval represented by the MC simulation - this gets us from an energy to a power, we already have powers, so this will be one in our formulation (in balance at least...)
	double alpha; //Computed spectral index for the cell
	double numin,numax; //Range of frequencies we are considering
	{
	double w; //the answer

	w=en1/(4.0*PI*v*dt);
	w*=(alpha+1.0);
	w/=(pow(numax,(alpha+1.0))-pow(numin,(alpha+1.0))); 

	return(w);
	}










