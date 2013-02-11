

/* This file was created to contain all of the routines associated
	with the attempt to calculate ionic abundances using a
	temperature fixed for each pair of ions to make the uncorrected
	abundance ration roughly equal to 1. There are then correction
	factors applied, either based upon a power law model of
	the true radiation field, or as a dilute blackbody.

	12Feb	nsh	Began coding 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"



struct photoionization *xver;	//Verner & Ferland description of a photoionization x-section
struct topbase_phot *xtop;	//Topbase description of a photoionization x-section
PlasmaPtr xxxplasma;
double qromb_temp;			//This is a storage variable for the current electron temperature so it is available for qromb calls




/***********************************************************
                                       Southampton University

  Synopsis:   

int
variable_temperature (xplasama, mode)  modifies the densities of ions, levels, and
	partition functions of ions within a cell of the wind based upon the mode, 
	and other data contained within the WindPtr itself, such as t_r, t_e, w, 
	based on the "mode".  

        Unlike nebular_concentrations, this code works on pairs of ions at a time, 
	and uses a temperature calculated to be suitable for that particular pair
	to have a ratio of abundances about equal to 1.
  
  Arguments:		
     PlasmaPtr ww;
     int mode;			// 6=correct using a dilute blackbody


  Returns:
 	variable temperature alters portions of the wind ptr.  Exactly how things 
	are changed depends on the mode.
 
 	variable_temperature returns 0 if it converged to an answer, -1 otherwise.  On an
 	abnormal return the density array and ne are left unchanged.
 	
  Description:	


 
  Notes:



  History:
	2012Feb	nsh	Coded and debugged as part of QSO effort. 

**************************************************************/


#define SAHA 4.82907e15		/* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
#define MAXITERATIONS	200
#define FRACTIONAL_ERROR 0.03
#define THETAMAX	 1e4
#define MIN_TEMP         100.

double xxxne,xip;

int
variable_temperature (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;			//   6=correct using dilute blackbody, 7=power law
{
  int nion, niterate;
  double xnew, xsaha, dne, dne_old;
  double theta, x;
  double get_ne();
  double t, t_e,t_r, xtemp, nh, xne, xxne, www;
  double a,b;
 // double xip; //ionzation potential of lower ion.
  int nelem,first,last;
  double *density;
  double *partition;
  double sum,big;
  double pi_fudge,recomb_fudge,gs_fudge,tot_fudge; /*The three correction factors for photoionization rate, recombination to ground state, and recombination rate */
  xxxplasma = xplasma;    /*Copy xplasma to local plasma varaible, used to communicate w and alpha to the power law correction routine. */

  dne_old=0.0;
  nh = xplasma->rho * rho2nh;	//LTE
  t_e = xplasma->t_e; 
  t_r = xplasma->t_r;
  www = xplasma->w;
  density = xplasma->density;  /* Set the local array density to the density held for the cell */
  partition = xplasma->partition; /* Set the partition function array to that held for the cell */


  /* make an initial estimate of ne based on H alone,  Our guess
     assumes ion[0] is H1.  Note that x below is the fractional
     ionization of H and should vary from 0 to 1

     Note that in a pure H plasma, the left hand side of the Saha 
     equation can be written as (x*x*nh*nh)/((1-x)*nh)  and can
     then be converted into a quadratic.  That is what is done
     below.  

     Since this is an initial estimate g factors have been ignored
     in what follows.
   */
  //  t_e=7.6753e3;
  if (t_e < MIN_TEMP)
    t_e = MIN_TEMP;		/* fudge to prevent divide by zeros */

  xsaha = SAHA * pow (t_e, 1.5);

  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t_e)) / nh;

  if (theta < THETAMAX)
    {
      x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
      xne = xxne = xxxne = x * nh;
    }
  else
    xne = xxne = xxxne = nh; /*xxne is just a store so the error can report the starting value of ne. xxxne is the shared variable so the temperature solver routine can access it*/

  if (xne < 1.e-6)
    xne = xxxne = 1.e-6;	   /* fudge to assure we can actually calculate
				   xne the first time through the loop */

  /* At this point we have an initial estimate of ne. */

  niterate = 0;
  while (niterate < MAXITERATIONS)
    {
       /* We loop over all elements */
    for (nelem = 0; nelem < nelements; nelem++)
    	{
      	first = ele[nelem].firstion;	/*first and last identify the position in the array */
      	last = first + ele[nelem].nions;	/*  So for H which has 2 ions, H1 and H2, first will generally be 0 and last will be 2 so the for loop below will just be done once for nion = 1 */
      	if (first < 0 || first >= nions || last < 0 || last > nions)
		{
	  	Error
	    	("variable_temperature: Confusion for element %d with first ion %d and last ion %d\n",
	     nelem, first, last);
	  	exit (0); 
		}
	/* and now we loop over all the ions of this element */
      	sum = density[first] = 1.0; /* set the density of the first ion of this element to 1.0 - this is (always??) the neutral ion */
      	big = pow (10., 250. / (last - first)); /*make sure that we will not overflow the last ion */

      	for (nion = first + 1; nion < last; nion++)  /*nion is the upper ion of the pair, last is one above the actual last ion, so the last ion to be considered is last-1 */
		{
		/* now we need to work out the correct temperature to use */
		xip=ion[nion-1].ip;  //the IP is that from the lower to the upper of the pair

		xtemp=zbrent(temp_func,MIN_TEMP,1e8,10);  //work out correct temperature
//		xtemp=1e8;
//		xtemp=(xip/H)/5.879e10; /* 71c use wiens law to calculate a temperature based on ioinzaation threshold of the lower ion which we will be using in a minute to calculate the denominator in the ionization rate correction factor */

		if (nion==1) printf ("HYDROGEN Best temp for ion %i is %e\n",nion,xtemp);

/* given this temperature, we need the pair of partition functions for these ions */		
		partition_functions_2 (xplasma, nion, xtemp);

/* and now we need the saha equation linking the two states at our chosen temp*/
		xsaha = SAHA * pow (xtemp, 1.5);
	  	b = xsaha * partition[nion]
	    	* exp (-ion[nion - 1].ip / (BOLTZMANN * xtemp)) / (xne *
							   partition[nion-1]);
		if (nion==1) printf ("HYDROGEN IP = %e\n",ion[nion - 1].ip/EV2ERGS);
		if (nion==1) printf ("HYDROGEN b = %e\n",b);

/* we now correct b to take account of the temperature and photon field 
		t_r and www give the actual radiation field in the cell, xtemp is the temp we used
  		t_e is the actual electron temperature of the cell*/


	      	recomb_fudge = sqrt (t_e / xtemp);
		if (nion==1) printf ("HYDROGEN t_e = %e xtemp= %e recomb_fudge= %e \n",t_e,xtemp,recomb_fudge);	
		gs_fudge = compute_zeta (t_e, nion, 1e14, 1e18, 1); /* Calculate the ground state recombination rate correction factor based on the cells true electron temperature. Zeta does the nion-1 bit */
//		printf ("recombination fudge = %e\n",recomb_fudge);
		if (nion==1) printf ("HYDROGEN gs_fudge= %e\n",gs_fudge);	


         	if (mode == 6) /* Correct the SAHA equation abundance pair using an actual radiation field modelled as a dilute black body */
			{
			pi_fudge = bb_correct_2 (xtemp, t_r, www, nion);
			tot_fudge=pi_fudge*recomb_fudge*(gs_fudge+www*(1-gs_fudge));
			}
		else if (mode == 7) /* correct the SAHA equation abundance pair using an actual radiation field modelled as a broken power law */
			{
			pi_fudge = pl_correct_2 (xtemp, nion);
			tot_fudge=pi_fudge*recomb_fudge*gs_fudge;
			}

 if (nion == 1)
	{
	printf ("nion=%2i, Z=%i, ion=%i, b=%e, pi=%e, recomb=%e, gs=%e, new=%e, big=%e, t_e=%f, xtemp=%f, ne=%e\n",nion,ion[nion].z,ion[nion].istate,b,pi_fudge,recomb_fudge,gs_fudge,b*tot_fudge,big,t_e,xtemp,xne);
}


		/* apply correction factors */
		
		b *= tot_fudge;



//		printf ("new ratio = %e\n",b);
          	if (b > big)	 
	   	b = big;		//limit step so there is no chance of overflow
	  	a = density[nion - 1] * b;
	  	sum += density[nion] = a;
	  	if (density[nion] < 0.0)
			{
	    		mytrap ();
			}
	  	sane_check (sum);
		}  //end of loop over ions - we now have a full set of ions for this element

      	a = nh * ele[nelem].abun / sum;   //the scaling factor to get the overall abundance right
      	for (nion = first; nion < last; nion++)
		{
	  	density[nion] *= a;  //apply scaling
	  	sane_check (density[nion]);  //check nothing has gone crazy
		}



      	if (geo.macro_ioniz_mode == 1)
		{
  		macro_pops (xplasma, xne);
		}

      /*Set some floor so future divisions are sensible */
      	for (nion = 0; nion < nions; nion++)
		{
	  	if (xplasma->density[nion] < DENSITY_MIN)
	    		xplasma->density[nion] = DENSITY_MIN;
		}

      	/* Now determine the new value of ne from the ion abundances */
	}   //end of loop over elements
    xnew = get_ne (xplasma->density);	/* determine the electron density for this density distribution */

    if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;	/* fudge to keep a floor on ne */
    if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	{
	break;
	}
    xne = xxxne = (xnew + xne)/2.;   /*New value of ne */

    niterate++;


    if (niterate == MAXITERATIONS)
    	{
      	Error ("variable_temperature: failed to converge t %.2g nh %.2g xnew %.2g\n",
	     t, nh, xnew);
      	Error ("variable_temperature: xxne %e theta %e\n", xxne, theta);
      	return (-1);
    	}
    }
  xplasma->ne = xnew;
  return (0);
}





/***********************************************************
                                       Southampton University

  Synopsis:   
	bb_correct_2(xtemp, t_r,www,nion) calculates the 
 	correction factor of the number density of two adjacent levels calculated
	at any temperature. LM is a special case of this where the temperatures
	are the same. 
 
  Arguments:		

	double xtemp		the temperature at which the saha abundances were calculated
	double t_r		temperature of the radiation field 
	double www		dilution factor for the radiation field
	int nion;		number of the upper ion in the pair for which the correction factor is 					being calculated.

  Returns:
 	j			the correction factor
 
 
 	
 Description:	


 
  Notes:



  History:
 
	12feb		NSH	Coded as part of the quasar project

**************************************************************/
int bb_correct_err=0;

double
bb_correct_2 (xtemp, t_r, www, nion)
     double t_r, xtemp, www;
     int nion;
{
  double q;
  int ion_lower,n;
  double numerator, denominator;
  int ntmin,nvmin;
  double fthresh,fmax;



  /* Initialization of fudges complete */

  ion_lower = nion-1;	/*the call uses nion, which is the upper ion in the pair */


  if (-1 < ion_lower && ion_lower < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[ion_lower].ntop_ground;
      nvmin = ion_lower;
    }
  else				
    {
      Error ("bb_correct_2: %d is unacceptable value of nion\n", ion_lower);
      mytrap ();
      exit (0);
      return (1.0);
    }


/*The integration is to correct for photoionization rates from the lower to the upper, so
 all integrals have to be over the lower ion cross section. The integrals use the code in
 power_sub, so the temperature has to be set to an external variable. Annoyingly this is 
 called sim_te....*/

 if (ion[ion_lower].phot_info == 1)  //topbase
    {
      n = ntmin;
      xtop = &phot_top[n];
      fthresh = xtop->freq[0];
      fmax = xtop->freq[xtop->np - 1];

      

      qromb_temp=t_r;  //The numerator is for the actual radiation temperature
      numerator=www*qromb (tb_planck1, fthresh, fmax, 1.e-4); //and is corrected for W
      if (numerator==0.0)
	return (0.0); //There is no need to waste time evaluating the denominator


      qromb_temp=xtemp; //The denominator is calculated for the LTE rate at our ideal temp

//	printf ("topbase n=%i,nion=%i,temp=%e,fthresh=%e,fmax=%e,hnu=%e,hnu/kt=%e,fmax=%e\n",n,ion_lower,temp,fthresh,fmax,fthresh*H,(fthresh*H)/(temp*BOLTZMANN),5.879e10*temp);



      denominator=qromb (tb_planck1, fthresh, fmax, 1.e-4);
    }
 else if (ion[ion_lower].phot_info == 0)  // verner
    {			
      n = nvmin;		//just the ground state ioinzation fraction.
      xver = &xphot[ion[n].nxphot];
      fthresh = xver->freq_t;
      fmax=xver->freq_max;


      qromb_temp=t_r;
      numerator = www*qromb (verner_planck1, fthresh, fmax, 1.e-4);
       if (numerator==0.0)
	return (0.0); //There is no need to waste time evaluating the denominator
      qromb_temp=xtemp;

//	printf ("verner n=%i,nion=%i,temp=%e,fthresh=%e,fmax=%e,hnu=%e,hnu/kt=%e,fmax=%e\n",n,ion_lower,temp,fthresh,fmax,fthresh*H,(fthresh*H)/(temp*BOLTZMANN),5.879e10*temp);



      denominator = qromb (verner_planck1, fthresh, fmax, 1.e-4);
    }

 else
    {
	    bb_correct_err++;
	    /* 71 - 111229 - ksl - Suppress this error after 100 instances so program does not bomb */
	    if (bb_correct_err<100){
      Error
	("bb_correct: No photoionization xsections for ion %d (element %d, ion state %d),setting sim_frac to 1\n",
	 nion,ion[nion].z,ion[nion].istate);
	    }
	    else if (bb_correct_err==100){
		    Error("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
	    }

	   
      q = 1.;
      return (q);
    }

  //    if (denominator < 1e-40)
//	denominator = 1e-40;


      q = numerator / denominator;


     


  return (q);
}

/* This next little bit of code was used to calculate an ideal temperature. This gets a temperature which gives a ratio of n to n+1 of close to 1, but at this temperature, the blackbody function used for the denominator of the correction factor is almost always zero. I've left the routine in, but it doesnt seem to work...*/



double
temp_func (solv_temp)
     double solv_temp;
{
  double answer;
//	printf ("xxxne=%e, xip=%e\n",xxxne,xip);
//	printf ("1=%e,2=%e,3=%e\n",log(4.83e15/xxxne),1.5*log(temp),(xip/(BOLTZMANN*temp)));
  answer = log(4.83e15/xxxne)+1.5*log(solv_temp)-(xip/(BOLTZMANN*solv_temp));

  return (answer);
}


/***********************************************************
                                       Southampton University

  Synopsis:   
	pl_correct_2(xtemp,nion) calculates the 
 	correction factor of the number density of two adjacent levels calculated
	at any temperature based upon a true radiation field modelled by a broken power law
 
  Arguments:		

	double xtemp		the temperature at which the saha abundances were calculated
	int nion;		the upper ion in an ion pair for which we want to correction factyor

  Returns:
 	q		        the correction factor
 	
 Description:	

	 

 
  Notes:
	The weights and alpha values needed for the integrals are communicated to this routine vias xxxplasma, an external varaible which is popullated in the variable_temperature routine.
	The integrals over power law require individual alphas and w, these are communicated via external variables xpl_alpha and xpl_w.	


  History:
 
	12feb	NSH	Coded as part of the quasar project

**************************************************************/
int pl_correct_err=0;
double xpl_alpha,xpl_w;



double
pl_correct_2 (xtemp, nion)
     double xtemp;
     int nion;
{
  double q;
  int ion_lower,n,j;
  double numerator, denominator;
  int ntmin,nvmin;
  double fthresh,fmax;



  /* Initialization of fudges complete */

  ion_lower = nion-1;	/*the call uses nion, which is the upper ion in the pair, we are interested in the PI from the lower ion */	

  if (-1 < ion_lower && ion_lower < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[ion_lower].ntop_ground;  /*We only ever use the ground state cross sections. This is for topbase */
      nvmin = ion_lower; /*and this is for verner cross sections*/
    }
  else				
    {
      Error ("bb_correct_2: %d is unacceptable value of nion\n", ion_lower);
      mytrap ();
      exit (0);
      return (1.0);
    }


/*The integration is to correct for photoionization rates from the lower to the upper, so
 all integrals have to be over the lower ion cross section. The integrals use external functions, so the temperature has to be set to an external variable temp. */

 if (ion[ion_lower].phot_info == 1)  //topbase
    {
      n = ntmin;

      xtop = &phot_top[n];
      fthresh = xtop->freq[0];
      fmax = xtop->freq[xtop->np - 1];
       if (nion==1) printf("HYDROGEN topbase threshold frequency= %e Hz %e ev\n",fthresh,fthresh*HEV);
	if (nion==1) printf("HYDROGEN fmax=%e \n",fmax);
	numerator=0;
	  for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
	    {
              if (nion==1) printf("HYDROGEN band %i runs from %e to %e\n",j,geo.xfreq[j],geo.xfreq[j+1]);
	      if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//Case 1- 
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator += qromb (tb_pow1, fthresh, fmax, 1.e-4);
		}
	      else if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j + 1] < fmax)	//case 2 
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator += qromb (tb_pow1, fthresh, geo.xfreq[j + 1], 1.e-4);
		}
	      else if (geo.xfreq[j] > fthresh && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//case 3
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator += qromb (tb_pow1, geo.xfreq[j], fmax, 1.e-4);
		}
	      else if (geo.xfreq[j] > fthresh && geo.xfreq[j + 1] < fmax)	// case 4
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator += qromb (tb_pow1, geo.xfreq[j], geo.xfreq[j + 1], 1.e-4);
		}
	      else		//case 5 - should only be the case where the band is outside the range for the integral.
		{
		  numerator += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		}
	    }
      if (numerator==0.0)
	return (0.0); //There is no need to waste time evaluating the denominator

      qromb_temp=xtemp; //The denominator is calculated for the LTE rate at our ideal temp. If we get this wrong, then divide by zeros abound!

//	printf ("topbase n=%i,nion=%i,temp=%e,fthresh=%e,fmax=%e,hnu=%e,hnu/kt=%e,fmax=%e\n",n,ion_lower,temp,fthresh,fmax,fthresh*H,(fthresh*H)/(temp*BOLTZMANN),5.879e10*temp);



      denominator=qromb (tb_planck1, fthresh, fmax, 1.e-4);
    }
 else if (ion[ion_lower].phot_info == 0)  // verner
    {			
      n = nvmin;		//just the ground state ioinzation fraction.
      xver = &xphot[ion[n].nxphot];
      fthresh = xver->freq_t;
      fmax=xver->freq_max;

	numerator=0;
	  for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
	    {
	      if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//Case 1
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator += qromb (verner_pow1, fthresh, fmax, 1.e-4);
		}
	      else if (geo.xfreq[j] < fthresh && fthresh < geo.xfreq[j + 1] && geo.xfreq[j + 1] < fmax)	//case 2 
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator +=
		    qromb (verner_pow1, fthresh, geo.xfreq[j + 1], 1.e-4);
		}
	      else if (geo.xfreq[j] > fthresh && geo.xfreq[j] < fmax && fmax < geo.xfreq[j + 1])	//case 3
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator += qromb (verner_pow1, geo.xfreq[j], fmax, 1.e-4);
		}
	      else if (geo.xfreq[j] > fthresh && geo.xfreq[j + 1] < fmax)	// case 4
		{
		  xpl_alpha = xxxplasma->sim_alpha[j];
		  xpl_w = xxxplasma->sim_w[j];
		  numerator +=
		    qromb (verner_pow1, geo.xfreq[j], geo.xfreq[j + 1], 1.e-4);
		}
	      else		//case 5 - should only be the case where the band is outside the range for the integral.
		{
		  numerator += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		}
	    }
       if (numerator==0.0)
	return (0.0); //There is no need to waste time evaluating the denominator
/* Denominator is integral at LTE of our chosen temperature. */

      qromb_temp=xtemp;

//	printf ("verner n=%i,nion=%i,temp=%e,fthresh=%e,fmax=%e,hnu=%e,hnu/kt=%e,fmax=%e\n",n,ion_lower,temp,fthresh,fmax,fthresh*H,(fthresh*H)/(temp*BOLTZMANN),5.879e10*temp);



      denominator = qromb (verner_planck1, fthresh, fmax, 1.e-4);
    }

 else
    {
	    pl_correct_err++; /* If we get here, there are no cross sections available */
	    if (pl_correct_err<100){
      Error
	("pl_correct: No photoionization xsections for ion %d (element %d, ion state %d),setting sim_frac to 1\n",
	 nion,ion[nion].z,ion[nion].istate);
	    }
	    else if (pl_correct_err==100){
		    Error("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
	    }

	   
      q = 1.; /*This is really bad actually, this will leave the abundances all wrong.... */
      return (q);
    }



 //     printf ("numerator = %e, denominator = %e\n",numerator,denominator);
      q = numerator / denominator; /*Just work out the correction factor - hope it isnt infinite! */

if (nion == 1)
	{
	printf ("HYDROGEN numerator = %e, denominator = %e, q=%e\n",numerator,denominator,q);
}
     


  return (q);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

 tb_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
   ionisation cross section is significant. This is the function for ions with a topbase cross section NSH 16/12/10 


  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.

 ************************************************************************/


double
tb_planck1 (freq)
     double freq;
{
  double answer, bbe;
  bbe = exp ((H * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot_topbase (xtop, freq);
  answer /= freq;

  return (answer);
}




/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	
	verner_planck is the function a\nuB\nu and is called by Qromb in order to integrate over the frequency range where the
	ionisation cross section is significant. This is the function for ions with a verner cross section NSH 16/12/10 

  Arguments:  


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.

 ************************************************************************/


double
verner_planck1 (freq)
     double freq;
{
  double answer, bbe;
  bbe = exp ((H * freq) / (BOLTZMANN * qromb_temp));
  answer = (2. * H * pow (freq, 3.)) / (pow (C, 2));
  answer *= (1 / (bbe - 1));
//      answer*=weight;
  answer *= sigma_phot (xver, freq);
  answer /= freq;

  return (answer);
}

/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	
	tb_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation 
	cross section
  	to allow the numerator of the sim correction factor to be calulated for ions with topbase cross sections.

  Arguments:  


  Returns:

 Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.
 ************************************************************************/



double
tb_pow1 (freq)
     double freq;
{
  double answer;

  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));
  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.

  return (answer);
}



/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  
	verner_pow is the function to allow integration of a power law photon distribution multiplied by the inoisation cross section
   	to allow the numerator of the sim correction factor to be calulated for ions with verner cross sections.

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Feb NSH - written as part of the varaible temperature effort.
 ************************************************************************/


double
verner_pow1 (freq)
     double freq;
{
  double answer;

  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));
  answer *= sigma_phot (xver, freq);	// and finally multiply by the cross section.

  return (answer);
}



