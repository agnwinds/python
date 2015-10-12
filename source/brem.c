#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "atomic.h"
#include "python.h"

#include "log.h"


double constant;
double T_b;



double
emittance_brem (freqmin, freqmax, lum, t)
     double freqmin, freqmax, lum, t;
{
  double constant, emit;
  /* these are the frequencies over which the power law is defined - currently set to the
     equivalent of 2 to 10 keV */


#define   XFREQMIN  4.84e17
#define   XFREQMAX  2.42e18
  
  
  T_b=t;

  /* first we need to calculate the constant for the power law function */
  constant = lum / qromb(integ_brem,XFREQMIN,XFREQMAX,1e-4);
  /* now we need to work out the luminosity between our limited frequency range */
  /* we may need some checking routines to make sure that the requested frequency range is within the defined range,
     or it could default to zero outside the defined range */

  emit =  constant * qromb(integ_brem,freqmin,freqmax,geo.brem_temp);

  return (emit);
}




double
integ_brem (freq)
     double freq;
{
  double answer;
  answer = constant * pow(freq,-0.2) * exp ((-1.0 * H * freq) / (BOLTZMANN * geo.brem_temp));
  return (answer);
}