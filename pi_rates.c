#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"

struct topbase_phot *xtop;	//Topbase description of a photoionization x-section 

double qromb_temp;  //Temerature used in integrations.

double xpl_alpha, xpl_w, xpl_logw;
double xexp_temp, xexp_w;




double 
calc_pi_rate (ion_lower,xplasma,mode)
	PlasmaPtr xplasma;
	int ion_lower;
	int mode;
{
  double q;
  double xtemp;
  int  n, j;
  double pi_rate;
  int ntmin, nvmin;
  double fthresh, fmax, fmaxtemp;
  double f1, f2;
  double exp_qromb, pl_qromb;


  xtemp=xplasma->t_r;
  exp_qromb = 1e-4;
  pl_qromb = 1e-4;


  if (-1 < ion_lower && ion_lower < nions)	//Get cross section for this specific ion_number
    {
      ntmin = ion[ion_lower].ntop_ground;	/*We only ever use the ground state cross sections. This is for topbase */
      nvmin = ion_lower;	/*and this is for verner cross sections */
    }
  else
    {
      Error ("calc_pi_rate: %d is unacceptable value of nion\n", ion_lower);
      mytrap ();
      exit (0);
      return (1.0);
    }



  if (ion[ion_lower].phot_info == 1)	//topbase
    {
      n = ntmin;
      xtop = &phot_top[n];
    }
  else if (ion[ion_lower].phot_info == 0)	// verner
    {
      n = nvmin;		//just the ground state ionization fraction.
      xtop = &xphot_tab[ion[n].nxphot];
    }
  else
    {
//      pl_correct_err++;		/* If we get here, there are no cross sections available */
//      if (pl_correct_err < 100)
//	{
	  Error
	    ("calc_pi_rate: No photoionization xsections for ion %d (element %d, ion state %d), setting rate to 0\n",
	     ion_lower, ion[ion_lower].z, ion[ion_lower].istate);
	  pi_rate=0;
	  return (pi_rate);
//	}
//      else if (pl_correct_err == 100)
//	{
//	  Error
//	    ("xinteg_sim: Suppressing photoionization xsection error, but more photo xsections are needed in atomic data\n");
//	}
//    q = 1.;			/*This is really bad actually, this will leave the abundances all wrong.... */
//      return (q);
    }

  fthresh = xtop->freq[0];
  fmax = xtop->freq[xtop->np - 1];
  pi_rate = 0;





if (mode==1) //Modelled version of J
	{
      for (j = 0; j < geo.nxfreq; j++)	//We loop over all the bands
	{
	  xpl_alpha = xplasma->pl_alpha[j];
	  xpl_logw = xplasma->pl_log_w[j];
	  xexp_temp = xplasma->exp_temp[j];
	  xexp_w = xplasma->exp_w[j];
	  if (xplasma->spec_mod_type[j] > 0)	//Only bother doing the integrals if we have a model in this band
	    {
	      f1 = xplasma->fmin_mod[j]; //NSH 131114 - Set the low frequency limit to the lowest frequency that the model applies to
	      f2= xplasma->fmax_mod[j]; //NSH 131114 - Set the high frequency limit to the highest frequency that the model applies to
	      if (f1 < fthresh && fthresh < f2 && f1 < fmax && fmax < f2)	//Case 1- 
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, fthresh, fmax, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, fthresh, fmax, exp_qromb);
		    }
		}
	      else if (f1 < fthresh && fthresh < f2 && f2 < fmax)	//case 2 
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, fthresh, f2, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, fthresh, f2, exp_qromb);
		    }
		}
	      else if (f1 > fthresh && f1 < fmax && fmax < f2)	//case 3
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, f1, fmax, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, f1, fmax, exp_qromb);
		    }
		}
	      else if (f1 > fthresh && f2 < fmax)	// case 4
		{
		  if (xplasma->spec_mod_type[j] == SPEC_MOD_PL)
		    {
		      pi_rate += qromb (tb_logpow1, f1, f2, pl_qromb);
		    }
		  else
		    {
		      pi_rate += qromb (tb_exp1, f1, f2, exp_qromb);
		    }
		}
	      else		//case 5 - should only be the case where the band is outside the range for the integral.
		{
		  pi_rate += 0;	// Add nothing - bit of a null statement, but makes the code look nice.
		}
	    }			//End of loop to only integrate in this band if there is power
	}

}
else if (mode==2)  //planck
	{
  fmaxtemp = xtop->freq[xtop->np - 1];
  fmax = check_fmax (fthresh, fmaxtemp, xtemp);
  if (fthresh > fmax)
    {
      Error
	("pl_correct: After checking, fthresh has been set below fmin - we cannot compute denominator\n");
      q = 1.0;
      return (q);
    }
  qromb_temp = xtemp;		
  pi_rate = qromb (tb_planck1, fthresh, fmax, 1.e-4);
    }


return(pi_rate);
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







double
tb_logpow1 (freq)
     double freq;
{
  double answer;

//  answer = xpl_w * (pow (freq, (xpl_alpha - 1.0)));

  answer = pow(10,xpl_logw+(xpl_alpha-1.0)*log10(freq));


  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.

  return (answer);
}



/**************************************************************************
                    Southampton University


  Synopsis:  

  Description:	
	tb_exp is the function to allow integration of an exponential photon distribution multiplied by the inoisation 
	cross section
  	to allow the numerator of the correction factor to be calulated for ions with topbase cross sections.

  Arguments:  


  Returns:

 Notes:

This is almost identical to code written to compute the sim power law correction. It is copied here to make the coding easier, plus it is likely that it will supplant the earlier code so that can all be deleted.


 

  History:

12Aug Written by NSH as part of the effort to improve spectral modelling
 ************************************************************************/



double
tb_exp1 (freq)
     double freq;
{
  double answer;

  answer = xexp_w * exp ((-1.0 * H * freq) / (BOLTZMANN * xexp_temp));
  answer *= sigma_phot_topbase (xtop, freq);	// and finally multiply by the cross section.
  answer /= freq;
  return (answer);
}






