


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "atomic.h"
#include "python.h"


/***********************************************************
                                       Space Telescope Science Institute

Synopsis:
	partition_functions calculates the partition functions
	for a single ceell in the grid

Arguments:


Returns:
	
Description:

Mode here is identical to that used by nebular_concentrations, e.g

0 = LTE using t_r
1 = LTE using t_e
2 = Lucy and Mazzali



Notes:

	0800802 - This is a rewritten version of two routines, do_partitions
	and partition.   It is still true that levels is almost identical
	to some of the code here and therefore it is unclear why it is needed.
	levels is always called from partition_functions


History:
	080802	ksl	60b -- Brought levels routine into
			do_partitions since they are always
			called one after the other.  The routines
			are almost identical and it remains
			unclear why this was coded twice
	080802	ksl	Combined two routines do_partitions
			and partition into one in order to
			make this routine very similar to
			levels
	080804	ksl	Removed hubeny calculation of g from
			this code.  The advantage of the 
			Hubeny calculation was that it had
			a correction for density that we
			no longer had, but it was not being
			accessed by the code at all.

**************************************************************/



int
partition_functions (xplasma, mode)
     PlasmaPtr xplasma;
     int mode;
{
  int nion;
  double partition ();
  double t, weight;

  int n, m;
  int m_ground;			/* added by SS Jan 05 */
  double z, kt;


  if (mode == 0)
    {
      //LTE using t_r
      t = xplasma->t_r;
      weight = 1;
    }
  else if (mode == 1)
    {
      //LTE using t_e
      t = xplasma->t_e;
      weight = 1;
    }
  else if (mode == 2)
    {
      //Non LTE calculation with radiative weights
      t = xplasma->t_r;
      weight = xplasma->w;
    }
  else if (mode == 3)  
    {
     //Non LTE calculation with non BB radiation field. Use T_e to get partition functions, same as mode 1-
      t = xplasma->t_e;
      weight = 1;
    }
  else
    {
      Error ("partition_functions: Unknown mode %d\n", mode);
      exit (0);
    }

  /* Calculate the partition function for each ion in turn */
  kt = BOLTZMANN * t;

  for (nion = 0; nion < nions; nion++)
    {

      if (ion[nion].nlevels > 0)
	//Calculate data on levels using a weighed BB assumption
	{
	  m = ion[nion].firstlevel;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //Makes explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlevels; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else if (ion[nion].nlte > 0)
	//Calculate using "non-lte" levels
	{
	  m = ion[nion].first_nlte_level;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //This statement makes an explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlte; n++)
	    {
	      m++;
	      z +=
		weight * config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else
	{
	  z = ion[nion].g;
	}


      xplasma->partition[nion] = z;
    }

  levels (xplasma, mode);	//levels from reduced BB

  return (0);
}



/***********************************************************
                                      Southampton University 

Synopsis:
	cardona_part_func calculates partition functions using the data presetned in
		cardona et al 2010.

Arguments:
	xplasma - the cell under consideration

Returns:
	
Description:



Notes:

	This is experimental code - no mode is used. We may wish to extend it so that it uses 
	modes, calculates partition function using weghted BB where cardona data is not
	available etc.


History:
	130612	nsh	72b -- Written as part of the attempt to model AGN.
			problems with making python agree with cloudy, even
			in LTE prompted an investigation of wether better 
			partition function data would help. 


**************************************************************/



int
cardona_part_func (xplasma)
     PlasmaPtr xplasma;
{
  int nion;
  double partition ();
  double t;

  int m;
  int m_ground;			/* added by SS Jan 05 */
  double z, kt,x;
  double N;    /*Total number density */
  double q,nstar,Ehat,term1,term2,term3;    


      t = xplasma->t_e;
	kt = BOLTZMANN*t;
	printf ("t=%e,kt=%e\n",t,kt);
      N=xplasma->rho*rho2nh+xplasma->ne;  /*We will estimate the number density of partitcles in the cell as the sum of hydrogen density and electron density. This could be wrong by a factor of around 2.3/2 */
 // for (nion = 0; nion < nions; nion++)
//	{
 //	N += xplasma->density[nion];  //We need the number density of particles - excessive to do all ions, but accurate!
	printf ("N=%e\n",N);
//	}
 
  for (nion = 0; nion < nions; nion++)
    {
	printf ("test %i\n",ion[nion].cpartflag);
      if (ion[nion].cpartflag > 0) //If we have cardona partition functions try this
	{
	  m = ion[nion].firstlevel;
	  m_ground = m;
	  term1 = config[m].g; /*Assume the first term is just the ground state - this may slightly incorrect if the ground state is a multiplet */
	printf ("term1=%e\n",term1);
	  x=(ion[nion].istate/(2.0*PI*BOHR));  /*Our 'level' is actually the ionization state+1, so we dont need to add 1 to get effective Z */
	  q = sqrt(x)*pow(N,(-1.0/6.0));  /*q*/
          nstar = (q/2.0)*(1.0+sqrt(1.0+(4.0/q))); /*  N*  */
	  Ehat = (ion[nion].ip)-((ion[nion].istate*ion[nion].istate*RYD2ERGS)/(nstar*nstar));
	  term2=cpart[ion[nion].nxcpart].part_G*exp((-1.0*cpart[ion[nion].nxcpart].part_eps)/(t));
	  term3=(cpart[ion[nion].nxcpart].part_m/3.0)*(pow(nstar,3.0)-343.0)*exp((-1.0*Ehat)/(kt));
	  z=term1+term2+term3; /* Add all the terms together */
	printf ("term1=%e,term2=%e,term3=%e,z=%e\n",term1,term2,term3,z);
	}
 
      else
	{
	  z = ion[nion].g; /* This should be done better - we just default to the ground state if no cardona data */
	}

	printf ("Setting partition function for ion %i to %e\n",nion,z);
      xplasma->partition[nion] = z;
		printf ("Partition function for ion %i is %e\n",nion,xplasma->partition[nion]);
    }


  return (0);
}


/***********************************************************
                                       Southampton university

Synopsis:
	partition_functions_2 calculates the partition functions
	for a pair of states for a given temperature. This is needed
	to support the pairwise ioinzation state calculaions where 
	the saha equation is applied to a pair of states at a
	useful temperature, then corrected. Part of this is 
	a requirement to get the partition functions at that 
	temperature for just tow two states of interest. It is 
	wasteful to calculate all of the states at each temperature.

Arguments:

	xplasma - the cell we are interested in. Used to communicate back
		the new partition functions.
	xnion - the upper ion in the pair we are currently working on
	temp - the temperature we are using for the saha equation. In the
		original formulation of partition, this is got from xplasma,
		in this case we don't want to use the cell temperature.

Returns:

	changes the partition functions of just two ions we are currently working on
	
Description:



Notes:

	There is no need for the weight term in the calculation since this is only called when
		we are making the assumption that we are in LTE for the saha equation. 



History:
	2012Feb	nsh - began coding	

**************************************************************/



int
partition_functions_2 (xplasma, xnion, temp)
     PlasmaPtr xplasma;
     int xnion;	
     double temp;
{
  int nion;
  double partition ();

  int n, m;
  int m_ground;			
  double z, kt;

 
  /* Calculate the partition function for each ion in turn */
  kt = BOLTZMANN * temp;   

  for (nion = xnion-1; nion < xnion+1; nion++)
    {
	
      if (ion[nion].nlevels > 0)
	/*Calculate data on levels using blackbody (unweighted since the calculation is in LTE */
	{
	  m = ion[nion].firstlevel;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //Makes explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlevels; n++)
	    {
	      m++;
	      z +=
		config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else if (ion[nion].nlte > 0)
	//Calculate using "non-lte" levels
	{
	  m = ion[nion].first_nlte_level;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  z = config[m].g;
	  //This statement makes an explicit assumption that first level is ground

	  for (n = 1; n < ion[nion].nlte; n++)
	    {
	      m++;
	      z +=
		config[m].g *
		exp ((-config[m].ex + config[m_ground].ex) / kt);
	    }
	}
      else
	{
	  z = ion[nion].g;
	}


      xplasma->partition[nion] = z;
    }


  return (0);
}




/***********************************************************
                                      Southampton University 

Synopsis:
	cardona_part_func_2 calculates partition functions using the data presetned in
		cardona et al 2010 for a pair of ions. 

Arguments:
	xplasma - the cell under consideration
	xnion - the upper ion of the pair
	temp - the temperature we want to use

Returns:
	
Description:



Notes:

	This is experimental code - no mode is used. We may wish to extend it so that it uses 
	modes, calculates partition function using weghted BB where cardona data is not
	available etc.


History:
	130612	nsh	72b -- Written as part of the attempt to model AGN.
			problems with making python agree with cloudy, even
			in LTE prompted an investigation of wether better 
			partition function data would help. 


**************************************************************/



int
cardona_part_func_2 (xplasma,xnion,temp)
     PlasmaPtr xplasma;
	int xnion;
	float temp;
{
  int nion;
  double partition ();

  int m;
  int m_ground;			
  double z, kt,x;
  double N;    /*Total number density */
  double q,nstar,Ehat,term1,term2,term3;    


 
	kt = BOLTZMANN*temp;
      N=xplasma->rho*rho2nh+xplasma->ne;

	printf ("temp=%e, N=%e\n",temp,N);

 
  for (nion = xnion-1; nion < xnion+1; nion++)
    {

      if (ion[nion].cpartflag > 0) //If we have cardona partition functions try this
	{
	  m = ion[nion].firstlevel;
	  m_ground = m;
	  //store ground state - in case energy neq 0(SS)
	  term1 = config[m].g;
	  //Makes explicit assumption that first level is ground
	  x=(ion[nion].istate/(2.0*PI*BOHR));  //Our 'level' is actually the ionization state+1
	  q = sqrt(x)*pow(N,(-1.0/6.0)); 
          nstar = (q/2.0)*(1.0+sqrt(1.0+(4.0/q)));
	  Ehat = (ion[nion].ip)-((ion[nion].istate*ion[nion].istate*RYD2ERGS)/(nstar*nstar));
	  term2=cpart[ion[nion].nxcpart].part_G*exp((-1.0*cpart[ion[nion].nxcpart].part_eps)/(temp));
	  term3=(cpart[ion[nion].nxcpart].part_m/3.0)*(pow(nstar,3.0)-343.0)*exp((-1.0*Ehat)/(kt));
	printf ("term1=%e, term2=%e, term3=%e, temp=%e\n",term1,term2,term3,temp);
	  z=term1+term2+term3;
	}
 
      else
	{
	  z = ion[nion].g;
	}


      xplasma->partition[nion] = z;
    }


  return (0);
}

