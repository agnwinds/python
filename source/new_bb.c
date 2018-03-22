
/***********************************************************/
/** @file  new_bb.c
 * @Author ksl
 * @date   January, 2018
 *
 * @brief  Routines concerning the generation of photons from blackbodies
 *
 * ???
 ***********************************************************/

//OLD /***********************************************************
//OLD                                        Space Telescope Science Institute
//OLD 
//OLD  Synopsis:
//OLD 	
//OLD double planck(t,freqmin,freqmax)	returns the frequency of a "pseudo-randomly generated photon. 
//OLD double t,freqmin,freqmax;		temperature of BB, and frequency limits
//OLD   
//OLD  					The frequency is pseudo random in the following limited sense.  
//OLD 					The photons are returned are weighted so that the energy distribution 
//OLD  					of a planck function is approximately reproduced. Within an energy 
//OLD  					bin the photon frequency is uniform.  
//OLD 
//OLD double
//OLD emittance_bb (freqmin, freqmax, t)
//OLD 
//OLD double freqmin, freqmax, t;
//OLD 					Calculate the emittance of a bb between freqmin and freqmax
//OLD 
//OLD double integ_planck_d(alphamin,alphamax)
//OLD double alphamin,alphamax;
//OLD 
//OLD 					returns the integral of the dimensionless blackbody function
//OLD 					between alphamin and alphamax.  
//OLD 
//OLD 					Because we cannot afford to keep reintegrating this integral, 
//OLD 					the following procedure is used.  The first time the function is
//OLD 					called the integ_planck[] is filled.  integ_plank[n]  contains the 
//OLD 					integral of the dimensionless planck function from  ALPHAMIN to 
//OLD 					ALPHAMAX.  Therefore if one wants to obtain the integral of the 
//OLD 					dimensionless bb function, one simply interpolates this array
//OLD 
//OLD 
//OLD init_integ_planck_d() 			performs the integration of the dimensionless bb function
//OLD 					in order to fill the arry integ_planck.  This routine needs only to be
//OLD 					 called onetime. The actuall integration is done with the numerical 
//OLD 					 recipes routine "qromb" .
//OLD 
//OLD check_fmax(fmin,fmax,temp)		This is a little helper routine written by NSH in August 2012. We
//OLD 					were having lots of problems with trying to integrate cross sections
//OLD 					x a plack function at temperatures where there is no flux at fmax. This
//OLD 					little routine just checks if fmax is going to be more than 100xt/(h/k)
//OLD 					which will return tiny numbers for planck, and make qromb return
//OLD 					nonsense.
//OLD 
//OLD Arguments:		
//OLD 
//OLD Returns:
//OLD  
//OLD Description:	
//OLD 		The way this works currently is a little obscure.  The primary routines are emittance_bb or planck.  
//OLD 		They call the other routines
//OLD 
//OLD 		The first call to either of these routines (surely emittance_bb) results in a call to integ_planck_d, 
//OLD 		which in turn calls integ_planck_init.  This initialization routine in turn populates the 
//OLD 		array integ_planck , which contains the integral of the bb function in an array.
//OLD 
//OLD 		bb_emittance continues to access the array integ_plank through integ_planck_d every 
//OLD 		future time is is called.
//OLD 
//OLD 		planck does the same thing albeit more indirectly. It sets up a cdf each time new frequency 
//OLD 		limits are placed on it.  planck therefore really uses the
//OLD 		cdf.
//OLD 
//OLD 		
//OLD Notes:
//OLD 		It seems possible that some of these routines could be combined in a clever way
//OLD 		
//OLD 
//OLD History:
//OLD  	97	ksl	The various portions of this file were written in 1997 as part of the python effort
//OLD  	98mar	ksl	Combined the routines here from bb_gen.c and planck.c to reduce the total number of
//OLD  			files.
//OLD  	98mar23	ksl	planck(t,freqmin,freqmax) modified to use new pdf routines
//OLD  	98apr25	ksl	integ_planck_d modified to fix an error when the ranges of alphamin and alphamax
//OLD  			were out of bounds
//OLD 	98jul27	ksl	planck moddified to force denser selection of gridpoints for small and high freq relative
//OLD 			to kt
//OLD 	01dec	ksl	Increased number of points at high end as part of exercise to force selection
//OLD 			of some high frequency photons (python_40)
//OLD 	04mar	ksl	Modified to provide a better treatment when the frequency range is
//OLD 			above or below the frequency limits of the cdf that was calculated 
//OLD 			initially.
//OLD 	05jul	ksl	56d -- Relegated writing the pdf file to a case where you were debugging
//OLD 	06aug	kls	57h -- Recommented and reorganized the file, so that could complete a partial
//OLD 			fix to a problem Sturat had identified.  The partial fix was generating 
//OLD 			photons outsde the min and maximum frequency
//OLD 	12aug	nsh	73e -- check_fmax written to check qromb calls will work properly.
//OLD 	12nov	ksl	74a_ksl -- Revised  rand_exp, which is used to allow us to reasonalbly 
//OLD 			approximate photons in the high alpha limit.  Eliminated some old notes about
//OLD 			this routine that seemed no longer very relevant.\
//OLD 	13mar	jm	74b5_JM -- fixed bug JM130303 in init_integ_planck. This bug was causing 
//OLD 			the emittance in the range wavemin-wavemax in the spectral cycles to be calculated 
//OLD 			wrongly. Was noticed when testing YSO models.
//OLD 
//OLD **************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "atomic.h"
#include "python.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//OLD /***********************************************************
//OLD                                        Space Telescope Science Institute
//OLD 
//OLD  Synopsis: planck returns the frequency of a "pseudo-randomly generated photon.  The
//OLD    frequency is pseudo random in the following limited sense.  The photons are returned
//OLD    are weighted so that the energy distribution of a planck function is approximately
//OLD    reproduced. Within an energy bin the photon frequency is uniform.  
//OLD 	
//OLD 
//OLD Arguments:		
//OLD 
//OLD      temperature, 
//OLD      freqmin, freqmax  - frequency interval 
//OLD 
//OLD Returns:
//OLD 
//OLD      The frequency of a photon  drawn randomly from a BB function with 
//OLD      temperature T
//OLD  
//OLD Description:	
//OLD 
//OLD 	The first time one enters the program, a cdf for a diminensionaless
//OLD 	BB function is created.  On subseqent entries, when the temperature
//OLD 	or frequency limits are changed, we use standard routines to limit
//OLD 	what portion of the dimensionless cdf to use.
//OLD 
//OLD 	Because of of the rapid fall off of the proabability density function
//OLD 	at very low and high frequency, special routines exist to deal with
//OLD 	photons in these regions of the spectrum.
//OLD 
//OLD Notes:
//OLD 
//OLD 	12nov - ksl - I have deleted some of the old notes associated with
//OLD 	this routine.  When this routine was written, computers are a lot
//OLD 	slower than they are today, and I was quite concerned about the
//OLD 	possibility that this routine would take a lot of computer time.
//OLD 	Computers are much faster today, and this routine does not 
//OLD 	consitute a significant fraction of the time of the entire 
//OLD 	program.  It is possible, that a more accurate routine could
//OLD 	be created more simply today.
//OLD 
//OLD History:
//OLD 	06aug	ksl	57h --Recommented, as moved to fix a problem
//OLD 			that was returning photons that were occassionally outside
//OLD 			of the frequency limits.  The problem had been introduced
//OLD 			to handle the YSO case, where Stuart had found the uniform
//OLD 			distribution used to interpolate between frequency points
//OLD 			annoyoing
//OLD 	06aug	ksl	Eliminated call to Stuart's attempt to better handle
//OLD 			the low frequency limit, in favor of a new routine
//OLD 			that I hope fixes the problem of the frequency
//OLD 			limits.  
//OLD 	06sep	ksl	Hopefully correct the remaining problems with this,
//OLD 			but there is a significant issue with all of our pdfs after 
//OLD 			we have created pdf arrays, we simply create a uniform 
//OLD 			distribution within an interval.
//OLD 	06sep	ksl	57i -- Have simplfied this to elimiante for now any lo or hi
//OLD 			frequency special cases
//OLD 	12nov	ksl	Removed some of the old Notes assocaited with this routine.
//OLD 			See version earlier than 74 for these old notes.
//OLD 	17jul	nsh - changed references to PDFs to CDFs
//OLD     17jul 	nsh - also changed how the BB cdf is generated - low alpha jumps were not really doing anything!
//OLD **************************************************************/

#define ALPHAMIN 0.4            // Region below which we will use a low frequency approximation
#define ALPHAMAX 30.            // Region above which we will use a high frequency approximation
#define ALPHABIG 100.           //  Region over which can maximmally integrate the Planck function
#define NJUMPS 30               //The number of 'jumps' - places in the CDF where we want to force points
#define NMAX 		1000        //The number of points at which the planck function is integrated between ALPHAMIN and ALPHAMAX for storage/

int ninit_planck = 0;           //A flag to say wether we have computed our stored blackbody integral

double old_t = 0;
double old_freqmin = 0;
double old_freqmax = 0;
double alphamin, alphamax;
double cdf_bb_lo, cdf_bb_hi, cdf_bb_tot;        // The precise boundaries in the the bb cdf 
double cdf_bb_ylo, cdf_bb_yhi;  // The places in the CDF defined by freqmin & freqmax
double lo_freq_alphamin, lo_freq_alphamax, hi_freq_alphamin, hi_freq_alphamax;  //  the limits to use for the low and high frequency values

// bb_set is thae array that cdf_gen_from_func uses to esablish the 
// specific points in the cdf of the dimensionless bb function.

/* These are what we call 'jumps' and are used by cdf_gen_from_func to 
ensure important parts of the CDF have points */
double bb_set[] = {

	  10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
	  19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29.
	};


int error_bb_hi = 0;
int error_bb_lo = 0;


/**********************************************************/
/** @name      planck
 * @brief      returns the frequency of a "pseudo-randomly generated photon.  The
 *    frequency is pseudo random in the following limited sense.  The photons are returned
 *    are weighted so that the energy distribution of a  function is approximately
 *    reproduced. Within an energy bin the photon frequency is uniform.
 *
 * @param [in out] double  t   ???
 * @param [in out] double  freqmin   ???
 * @param [in out] double  freqmax   - frequency interval
 * @return     The frequency of a photon  drawn randomly from a BB function with 
 *      temperature T
 *
 * @details
 * The first time one enters the program, a cdf for a diminensionaless
 * 	BB function is created.  On subseqent entries, when the temperature
 * 	or frequency limits are changed, we use standard routines to limit
 * 	what portion of the dimensionless cdf to use.
 * 
 * 	Because of of the rapid fall off of the proabability density function
 * 	at very low and high frequency, special routines exist to deal with
 * 	photons in these regions of the spectrum.
 *
 * ### Notes ###
 * 12nov - ksl - I have deleted some of the old notes associated with
 * 	this routine.  When this routine was written, computers are a lot
 * 	slower than they are today, and I was quite concerned about the
 * 	possibility that this routine would take a lot of computer time.
 * 	Computers are much faster today, and this routine does not 
 * 	consitute a significant fraction of the time of the entire 
 * 	program.  It is possible, that a more accurate routine could
 * 	be created more simply today.
 *
 **********************************************************/

double
planck (t, freqmin, freqmax)
     double t, freqmin, freqmax;
{
  FILE *fopen ();
  double freq, alpha, y;
  double planck_d (), cdf_get_rand_limit ();
  double get_rand_pow ();
  int cdf_gen_from_func (), cdf_to_file (), echeck;
  int cdf_limit ();
  
  


  /*First time through create the array containing the proper boundaries for the integral of the BB function,
     Note calling cdf_gen_from func also defines ylo and yhi */

  if (ninit_planck == 0)
  {                             /* First time through p_alpha must be initialized */
    if ((echeck = cdf_gen_from_func (&cdf_bb, &planck_d, ALPHAMIN, ALPHAMAX, 20, bb_set)) != 0)
    {
      Error ("Planck: on return from cdf_gen_from_func %d\n", echeck);
    }
    /* We need the integral of the bb function outside of the regions of interest as well */


    cdf_bb_tot = qromb (planck_d, 0, ALPHABIG, 1e-8);
    cdf_bb_lo = qromb (planck_d, 0, ALPHAMIN, 1e-8) / cdf_bb_tot;       //position in the full cdf of low frequcny boundary
    cdf_bb_hi = 1. - qromb (planck_d, ALPHAMAX, ALPHABIG, 1e-8) / cdf_bb_tot;   //postion in fhe full hi frequcny boundary

//      cdf_to_file (&cdf_bb, "cdf.out");
    ninit_planck++;

  }


/* If temperatures or frequencies have changed since the last call to planck
redefine various limitsi, including the region of the pcdfdf  to be used

Note - ksl - 1211 - It is not obvious why all of these parameters need to be
reset.  A careful review of them is warranted.
*/

  if (t != old_t || freqmin != old_freqmin || freqmax != old_freqmax)
  {
	  
    alphamin = H * freqmin / (BOLTZMANN * t);
    alphamax = H * freqmax / (BOLTZMANN * t);

    old_t = t;
    old_freqmin = freqmin;
    old_freqmax = freqmax;

    cdf_bb_ylo = cdf_bb_yhi = 1.0;
    if (alphamin < ALPHABIG)  //check to make sure we get a sensible number - planck_d(ALPHAMAX is too small to sensibly integrate)
    {
      cdf_bb_ylo = qromb (planck_d, 0, alphamin, 1e-8) / cdf_bb_tot;    //position in the full cdf of current low frequency boundary
      if (cdf_bb_ylo > 1.0)
        cdf_bb_ylo = 1.0;
    }
    if (alphamax < ALPHABIG)  //again, check to see that the integral will be sensible 
    {
      cdf_bb_yhi = qromb (planck_d, 0, alphamax, 1e-8) / cdf_bb_tot;    //position in the full cdf of currnet hi frequency boundary
      if (cdf_bb_yhi > 1.0)
        cdf_bb_yhi = 1.0;
    }

/* These variables are not always used */
    lo_freq_alphamin = alphamin;        //Set the minimum frequency to use the low frequency approximation to the lower band limit
    lo_freq_alphamax = alphamax;        //Set to a default value
    if (lo_freq_alphamax > ALPHAMIN)    //If the upper alpha for this band is above the loew frequency approximation lower limit
      lo_freq_alphamax = ALPHAMIN;      //Set the maximum alpha we will use the low frequency approximation to the default value

    hi_freq_alphamax = alphamax;         //Set the maximum frequency to use the high frequency approximation to to the upper band limit
    hi_freq_alphamin = alphamin;         //Set to a default value
    if (hi_freq_alphamin < ALPHAMAX)     //If the lower band limit is less than the high frequency limit
      hi_freq_alphamin = ALPHAMAX;       //Se the minimum alpha value to use the high frequency limit to the default value


    if (alphamin < ALPHAMAX && alphamax > ALPHAMIN) //Since alphamin is always below alphamax, this is saying we are within the 'normal' bb range.
    {
      cdf_limit (&cdf_bb, alphamin, alphamax);
    }

  }
  /* End of section redefining limits */







//  y = rand () / (MAXRAND);   //We get a random number between 0 and 1 - DONE
  y= random_number(0.0,1.0); //We get a random number between 0 and 1 (excl)

  y = cdf_bb_ylo * (1. - y) + cdf_bb_yhi * y;   // y is now in an allowed place in the cdf
  

	  

/* There are 3 cases to worry about
	The case where everything is in the low frequency limit
	The case where everything is in the normal limit
	The case where some photons are in the low regime and some are
	in the normal regime
*/

  if (y <= cdf_bb_lo || alphamax < ALPHAMIN) //we are in the low frequency limit
  {
    alpha = get_rand_pow (lo_freq_alphamin, lo_freq_alphamax, 2.);
  }
  else if (y >= cdf_bb_hi || alphamin > ALPHAMAX) //We are in the high frequency limit
  {
    alpha = get_rand_exp (hi_freq_alphamin, hi_freq_alphamax);
  }
  else 
  {
    alpha = cdf_get_rand_limit (&cdf_bb); //We are in the region where we use the BB function
  }

  freq = BOLTZMANN * t / H * alpha;
  if (freq < freqmin || freqmax < freq)
  {
    Error ("planck: freq %g out of range %g %g\n", freq, freqmin, freqmax);
  }
  return (freq);
}



//OLD /***********************************************************
//OLD 	Space Telescope Science Institute
//OLD 
//OLD  Synopsis: get_rand_pow obtains a random number between x1 and x2 
//OLD  	for a power law densiity distribution with index alpha
//OLD 	
//OLD 
//OLD Arguments:		
//OLD 
//OLD Returns:
//OLD  
//OLD Description:	
//OLD 
//OLD Notes:
//OLD 		
//OLD History:
//OLD 	06sep	ksl	Coded       
//OLD **************************************************************/

/**********************************************************/
/** @name      get_rand_pow
 * @brief      obtains a random number between x1 and x2 
 *  	for a power law densiity distribution with index alpha
 *
 * <NOTE: The [in out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @param [in out] double  x1   ???
 * @param [in out] double  x2   ???
 * @param [in out] double  alpha   ???
 * @return     ??? RETURNS ???
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
get_rand_pow (x1, x2, alpha)
     double x1, x2, alpha;
{
  double r;
  double a;

//  r = rand () / MAXRAND; DONE
  r = random_number(0.0,1.0); //This produces a random number between 0 and 1 excl

  if (alpha == -1)
  {
    x1 = log (x1);
    x2 = log (x2);
    a = (1. - r) * x1 + r * x2;
    a = exp (a);
  }
  else
  {
    x1 = pow (x1, alpha + 1.);
    x2 = pow (x2, alpha + 1.);

    a = (1. - r) * x1 + r * x2;

    a = pow (a, 1. / (alpha + 1.));
  }
  return (a);
}


//OLD /***********************************************************
//OLD 	Space Telescope Science Institute
//OLD 
//OLD  Synopsis: get_rand_exp obtains a random number between x1 and x2 
//OLD  	for an exp law density distribution with index alpha. 
//OLD 	
//OLD 
//OLD Arguments:		
//OLD 
//OLD 	alpha_min and alpha_max, the range of the exponential
//OLD 	function over which one wants the random number 
//OLD 	to be returned
//OLD 
//OLD Returns:
//OLD 
//OLD 	A double precision random number
//OLD  
//OLD Description:	
//OLD 
//OLD 	
//OLD Notes:
//OLD 
//OLD 	The cdf for an exponential distribution can be easily
//OLD 	shown to be given by solving this equation for alpha
//OLD 
//OLD 	r*(exp(-alpha_min)-exp(-alpha_max))=(exp(-alpha_min)-exp(alpha))
//OLD 
//OLD 	but you can recast this and solve for delta_alpha
//OLD 
//OLD 	exp(-delta_alpha)=(1-R)+R*exp(-(alpha_max-alpha_min))
//OLD 
//OLD 	This has the advantage that it depends only on the 
//OLD 	difference of alpha_max and alpha_min and not their
//OLD 	actual values, and as long as the exp of a very 
//OLD 	large number turns out to be zero and not 
//OLD 	not a number, it shuld not generate sane checks
//OLD 		
//OLD History:
//OLD 	06oct	ksl	Coded       
//OLD 	12nov	ksl	Added sane check to find a problem with
//OLD 			returning NaN.  The problem was due to
//OLD 			the fact that alpha_min and almax was so
//OLD 			large that x1 and x2 were both identically
//OLD 			zero and as a result the logarithm was
//OLD 			negative.  
//OLD 	12nov	ksl	Changed algorithm to minimimize the effects
//OLD 			of very large values of alpha_min and 
//OLD 			alpha_max.  Instead of calculating 
//OLD 			alpha directly, we calculate delta_alpha
//OLD 			the difference between alpha_min and the
//OLD 			value we want
//OLD **************************************************************/

/**********************************************************/
/** @name      get_rand_exp
 * @brief      obtains a random number between x1 and x2 
 *  	for an exp law density distribution with index alpha.
 *
 * <NOTE: The [in out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @param [in out] double  alpha_min   and alpha_max
 * @param [in out] double  alpha_max   ???
 * @return     A double precision random number
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * The cdf for an exponential distribution can be easily
 * 	shown to be given by solving this equation for alpha
 * 
 * 	r*(exp(-alpha_min)-exp(-alpha_max))=(exp(-alpha_min)-exp(alpha))
 * 
 * 	but you can recast this and solve for delta_alpha
 * 
 * 	exp(-delta_alpha)=(1-R)+R*exp(-(alpha_max-alpha_min))
 * 
 * 	This has the advantage that it depends only on the 
 * 	difference of alpha_max and alpha_min and not their
 * 	actual values, and as long as the exp of a very 
 * 	large number turns out to be zero and not 
 * 	not a number, it shuld not generate sane checks
 *
 **********************************************************/

double
get_rand_exp (alpha_min, alpha_max)
     double alpha_min, alpha_max;
{
  double r;
  //Old ksl 12nov double x1, x2;
  double x;
  double a, aa;
  double delta_alpha;

//  r = rand () / MAXRAND; //DONE - this used to produce a number betwen 0 and 1 incl
  r = random_number(0.0,1.0); //A random number between 0 and 1 excl

  x = exp (alpha_min - alpha_max);


  aa = (1. - r) + r * x;
  delta_alpha = -(log (aa));

  a = alpha_min + delta_alpha;

  if (sane_check (a))
  {
    Error ("get_rand_exp:sane_check %e %e %e %e %e\n", a, aa, delta_alpha, x, r);
  }
  return (a);
}

//OLD /***********************************************************
//OLD 	Space Telescope Science Institute
//OLD 
//OLD  Synopsis: integ_plank_d(alphamin,alphamax) returns the integral of the dimensionless blackbody function
//OLD    between alphamin and alphamax.  
//OLD 	
//OLD 
//OLD Arguments:		
//OLD 
//OLD Returns:
//OLD  
//OLD Description:	
//OLD 	To save computing time, the routine actually accesses an array that contains the
//OLD 	the integral of the bb function from 0 to x.
//OLD 
//OLD    	The first time the function is
//OLD 	called the integ_planck[] is filled.  integ_plank[n]  contains the integral of the
//OLD 	dimensionless planck function from  ALPHAMIN to ALPHAMAX.  Therefore if one
//OLD 	wants to obtain the integral of the dimensionless bb function, one simply interpolates
//OLD 	this array. 
//OLD Notes:
//OLD 		
//OLD History:
//OLD 	06aug	ksl	Recommented
//OLD **************************************************************/


double integ_planck[NMAX + 1];  //The dimensuionless planck function integrated from ALPHAMIN to a range of values of alpha
int i_integ_planck_d = 0;      //A flag to say wether we have initialised integ_planck

/**********************************************************/
/** @name      integ_planck_d
 * @brief      the integral of the dimensionless blackbody function
 *    between alphamin and alphamax.
 *
 * @param [in out] double  alphamin   dij
 * @param [in out] double  alphamax   ???
 * @return     ??? RETURNS ???
 *
 * @details
 * To save computing time, the routine actually accesses an array that contains the
 * 	the integral of the bb function from 0 to x.
 * 
 *    	The first time the function is
 * 	called the integ_planck[] is filled.  integ_plank[n]  contains the integral of the
 * 	dimensionless planck function from  ALPHAMIN to ALPHAMAX.  Therefore if one
 * 	wants to obtain the integral of the dimensionless bb function, one simply interpolates
 * 	this array.
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
integ_planck_d (alphamin, alphamax)
     double alphamin, alphamax;
{
  double x, z1, z2;  //z1 and z2 are the integrals of the dimensionless BB at alphamin and alphamax
  int n;
  int init_integ_planck_d ();
  if (i_integ_planck_d == 0)    //If the flag isnt set
  {                             /*First time through integ_planck must be defined */
    init_integ_planck_d ();      //Initialise the integ_planck array
    i_integ_planck_d++;
  }

  x = (alphamin - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;  //Find the required lower point in the integ_planck array
  if (x <= 0.0)   //The lowest frequency is off the bottom of the tabulated values
    z1 = 0.0;  //We treat the integral as zero off the bottom end
  else if (x >= (NMAX))  // are off the upper end
  {
    return (0.0);               /*The minimum frequency is off the top of the integ_planck array we cannot use the table*/
  }
  else
  {
    n = x;  //n is the array element
    x -= n;  //x is now the fractional distance between array elements
    z1 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;  //Interpolate between array elements to get the integral at alphamin
  }

  x = (alphamax - ALPHAMIN) / (ALPHAMAX - ALPHAMIN) * NMAX;  //Find the required upper point in the integ_planck array
  if (x < 0.0)    //The highest frequency is blow the bottom of the tabulated values
  {
    return (0.0);               /* Because the maximum frequency is too low */
  }
  else if (x >= (NMAX))  //We are off the top of the array
  {
    z2 = integ_planck[NMAX];  //We assume the integral has been maximiased - just set the value for alphamax to the highest element in integ_planck
  }
  else
  {
    n = x;
    x -= n;
    z2 = integ_planck[n] * (1. - x) + integ_planck[n + 1] * x;  //Interpolate
  }

  return (z2 - z1);  //The required integral is just the integral at alphamax - the ibntegral at alphamin.
}


//OLD /***********************************************************
//OLD        Space Telescope Science Institute
//OLD 
//OLD  Synopsis: init_integ_planck_d calulates integrals of the dimensionaless bb function and
//OLD  	stores it in the array integ_planck.  Each value of integ_plank contains
//OLD 	the integration from 0 to x
//OLD 	
//OLD 
//OLD Arguments:		
//OLD 
//OLD Returns:
//OLD  
//OLD Description:	
//OLD    This routine needs to be called once, since ALPHAMIN and ALPHAMAX are
//OLD    hardcoded. 
//OLD 
//OLD    The actual integration is done with the numerical recipes routine "qromb" 
//OLD 
//OLD Notes:
//OLD 		
//OLD History:
//OLD 	97+	ksl	Originally coded
//OLD 	06aug	ksl	Revised comments
//OLD 
//OLD **************************************************************/


/**********************************************************/
/** @name      init_integ_planck_d
 * @brief      calulates integrals of the dimensionaless bb function and
 *  	stores it in the array integ_planck.  Each value of integ_plank contains
 * 	the integration from 0 to x
 *
 * <NOTE: The [in out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @return     ??? RETURNS ???
 *
 * @details
 * This routine needs to be called once, since ALPHAMIN and ALPHAMAX are
 *    hardcoded. 
 * 
 *    The actual integration is done with the numerical recipes routine "qromb"
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

int
init_integ_planck_d ()
{
  double x;
  double planck_d (), qromb ();
  int n;
  for (n = 0; n < NMAX + 1; n++)
  {
    x = ALPHAMIN + n * (ALPHAMAX - ALPHAMIN) / NMAX;
// 1e-7 is the fractional accuracy in my modified version of qromb -- ksl
    integ_planck[n] = qromb (planck_d, 0.0, x, 1e-7);
  }

  return (0);
}


//OLD /***********************************************************
//OLD 	Space Telescope Science Institute
//OLD 
//OLD  Synopsis:
//OLD "planck_d" is the dimensionless bb function.  The total emittance
//OLD    is related to planck_d as follows:
//OLD 	
//OLD 
//OLD Arguments:		
//OLD 
//OLD Returns:
//OLD  
//OLD Description:	
//OLD    F_nu= 2*PI* (kT/h**3)*planck_d(h*freq/ k T)
//OLD Notes:
//OLD 		
//OLD History:
//OLD 	06aug	ksl	Recommented
//OLD **************************************************************/

#define EPSILON	1.e-6


/**********************************************************/
/** @name      planck_d
 * @brief      "" is the dimensionless bb function.  The total emittance
 *    is related to  as follows:
 *
 * <NOTE: The [in out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @param [in out] double  alpha   ???
 * @return     ??? RETURNS ???
 *
 * @details
 * F_nu= 2*PI* (kT/h**3)*planck_d(h*freq/ k T)
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
planck_d (alpha)
     double alpha;
{
  double x;
  if (alpha < EPSILON || alpha > ALPHABIG)
    return (0);
  x = (alpha * alpha * alpha) / (exp (alpha) - 1);
  return (x);
}

// Calculate the emittance of a bb between freqmin and freqmax
// Should integrate to sigma 

//NSH - 17Jul - made change to if/else statement so that qromb used whenever bounds are not completely within the tabulated bands.

/**********************************************************/
/** @name      emittance_bb
 * @brief      ??? SYNOPSIS ???
 *
 * <NOTE: The [in out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @param [in out] double  freqmin   ???
 * @param [in out] double  freqmax   ???
 * @param [in out] double  t   ???
 * @return     ??? RETURNS ???
 *
 * @details
 * ??? DESCRIPTION ???
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
emittance_bb (freqmin, freqmax, t)
     double freqmin, freqmax, t;
{
  double alphamin, alphamax, q1;
  double integ_planck_d ();
  q1 = 2. * PI * (BOLTZMANN * BOLTZMANN * BOLTZMANN * BOLTZMANN) / (H * H * H * C * C);




  alphamin = H * freqmin / (BOLTZMANN * t);
  alphamax = H * freqmax / (BOLTZMANN * t);
  
  
  
  
  if (alphamin > ALPHAMIN && alphamax < ALPHAMAX) //We are within the tabulated range
  {
    return (q1 * t * t * t * t * integ_planck_d (alphamin, alphamax));
  }
  else if (alphamax > ALPHABIG) 
  {
	  if (alphamin > ALPHABIG) //The whole band is above the point where we can sensibly integrate the BB function
	  	return(0);
	  else   //only the upper part of the band is above ALPHABIG
	      return (q1 * t * t * t * t * qromb (planck_d, alphamin, ALPHABIG, 1e-7));
  }
  else //We are outside the tabulated range and must integrate
  {
    return (q1 * t * t * t * t * qromb (planck_d, alphamin, alphamax, 1e-7));
  }
}


//OLD /***********************************************************
//OLD 	Southampton University
//OLD 
//OLD  Synopsis:
//OLD 
//OLD check_fmax decides whether a maximum frequency requested for an integral is sensible.
//OLD If it is too far off the end of the planck function, qromb will malfunction. We
//OLD just have to set it to a frequency where the BB function is tiny, say where hnu/kT =100.
//OLD At this point the bb function is
//OLD 
//OLD Arguments:		
//OLD 
//OLD Returns:
//OLD  
//OLD Description:	
//OLD We use alphabig to define the place in the BB spectrum where we want to give up
//OLD Notes:
//OLD 		
//OLD History:
//OLD 	12aug	nsh	written
//OLD **************************************************************/



/**********************************************************/
/** @name      check_fmax
 * @brief      decides whether a maximum frequency requested for an integral is sensible.
 * If it is too far off the end of the planck function, qromb will malfunction. We
 * just have to set it to a frequency where the BB function is tiny, say where hnu/kT =100.
 * At this point the bb function is
 *
 * <NOTE: The [in out] tag describes if the value of a parameter is used or altered. If it is used but not altered, delete 'OUT'. If the original value is not used and it is written to, delete 'IN'.>
 * @param [in out] double  fmin   ???
 * @param [in out] double  fmax   ???
 * @param [in out] double  temp   ???
 * @return     ??? RETURNS ???
 *
 * @details
 * We use alphabig to define the place in the BB spectrum where we want to give up
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

double
check_fmax (fmin, fmax, temp)
     double fmin, fmax, temp;
{
  double bblim;

  bblim = ALPHABIG * (temp / H_OVER_K); /*This is the frequency at which the exponent in the 
                                           planck law will be -100. This will give a *very* small b(nu). */
  if (bblim < fmax)
  {
    fmax = bblim;
  }

  return (fmax);

}




#undef NMAX
#undef ALPHAMIN
#undef ALPHAMAX
