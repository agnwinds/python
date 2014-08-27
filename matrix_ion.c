/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Saha_inv is a quick program to generate lots of saha abundances for a range of temperatures and electron densities.
 
Arguments:		


Returns:
 
Description:	
	
Notes:

History:
	 
	0212 - nsh wrote it
 	
 	Look in Readme.c for more text concerning the early history of the program.

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
//struct topbase_phot *xtop;
//double g1, g2, xne, xip;
#define SAHA 4.82907e15
//double xip, xxxne, qromb_temp;
#include "python.h"


#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>
#include "my_linalg.h"


int matrix_solv(xplasma,mode)
	PlasmaPtr xplasma;
	int mode;

{
int nn,mm,nrows;
double rate_matrix[nions][nions];
  double newden[NIONS];
double nh,t_e,t_r,www;
double xne,xxne,xxxne;
double xsaha,x,theta;
  int s;			//NEWKSL
  double b_temp[nions]; 
double *b_data,*a_data;
  gsl_permutation *p;		//NEWKSL
  gsl_matrix_view m;		//NEWKSL
  gsl_vector_view b;		//NEWKSL
  gsl_vector *populations, *test_vector;	//NEWKSL
  gsl_matrix *test_matrix;
  int ierr,niterate;
  double test_val,xnew;
  double xden[nelements];
  int xion[nions]; //This array keeps track of what ion is in each line
  int xelem[nions]; //This array keeps track of the element for each ion
  double pi_rates[nions];
  double rr_rates[nions];

  nh = xplasma->rho * rho2nh;	//LTE
  t_e = xplasma->t_e;
  t_r = xplasma->t_r;
  www = xplasma->w;
//printf ("nh=%e, t_e=%e t_r=%e\n",nh,t_e,t_r);

 compute_dr_coeffs (t_e);
 compute_di_coeffs (t_e);



 /* Copy the current densities into a temporary array */

  for (mm = 0; mm < nions; mm++)
    {
      newden[mm] = xplasma->density[mm]; 
      xion[mm] = mm;
      if (ion[mm].istate != 1) //We can recombine
		{
	  rr_rates[mm]=total_rrate (mm, xplasma->t_e);
//	printf ("RR_rate %i = %e\n",mm,rr_rates[mm]);
		}
      if (ion[mm].istate != ion[mm].z+1 ) //we can photoionize
	  {
	  if (mode == NEBULARMODE_MATRIX_BB)
		{
	  	pi_rates[mm]=calc_pi_rate (mm,xplasma,2);

		}
	  else if (mode == NEBULARMODE_MATRIX_SPECTRALMODEL)
		{
	  	pi_rates[mm]=calc_pi_rate (mm,xplasma,1);
//		printf ("PI_rate %i = %e\n",mm,pi_rates[mm]);
		}
  	  else
                {
// If reached this point the program does not understand what is desired.
                Error ("matrix_solv: Unknown mode %d\n", mode);
                exit (0);
                }
	   }
      for (nn=0;nn<nelements;nn++)
	{
  //      xden[nn]=0.0;
	if (ion[mm].z==ele[nn].z)
		{
		xelem[mm]=nn;
		}
	}
    }



xne = xxne = xxxne = get_ne (newden); /*Best guess for the starting n_e is from the last try! */

 //if (t_e < MIN_TEMP)
 //   t_e = MIN_TEMP;		/* fudge to prevent divide by zeros */

//  xsaha = SAHA * pow (t_e, 1.5);

//  theta = xsaha * exp (-ion[0].ip / (BOLTZMANN * t_e)) / nh;

//  if (theta < THETAMAX)
//    {
//      x = (-theta + sqrt (theta * theta + 4 * theta)) / 2.;
//      xne = xxne = xxxne = x * nh;
//    }
//  else
//    xne = xxne = xxxne = nh;	/*xxne is just a store so the error can report the starting value of ne. */
 


                             /*    xxxne is the shared variable so the temperature solver routine can access it */

 // if (xne < 1.e-6)
 //   xne = xxxne = 1.e-6;	/* Set a minimum ne to assure we can calculate
	//			   xne the first time through the loop */


//printf ("Initial guess at ne=%e\n",xne);




//printf ("Initial guess at ne=%e \n",xne);
  /* At this point we have an initial estimate of ne. */
//printf ("starting xne=%e\n",xne);

  niterate = 0;
  while (niterate < MAXITERATIONS)
{


/*First we initialise the matrix */
for (nn = 0; nn < nions; nn++)
	{
	  for (mm = 0; mm < nions; mm++)
	    {
	      rate_matrix[nn][mm] = 0.0;
	    }
	}


/*Now we populate the elements relating to PI depopulating a state*/


for (mm = 0; mm<nions; mm++)
	{
	if (ion[mm].istate != ion[mm].z+1) //we have electrons
		{
		rate_matrix[mm][mm]-=pi_rates[mm];
		}
	}





/*Now we populate the elements relating to PI populating a state*/

for (mm =0; mm<nions; mm++)
	{
	for (nn=0;nn<nions;nn++)
		{
		if (mm==nn+1 && ion[nn].istate != ion[nn].z+1 && ion[mm].z==ion[nn].z) 
			{
			rate_matrix[mm][nn]+=pi_rates[nn];
			}
		}
	}

/*Now we populate the elements relating to direct ionization depopulating a state*/


for (mm = 0; mm<nions; mm++)
	{
	if (ion[mm].istate != ion[mm].z+1 && ion[mm].dere_di_flag > 0) //we have electrons and a DI rate
		{
		rate_matrix[mm][mm]-=(xne*di_coeffs[mm]);
		}
	}


/*Now we populate the elements relating to direct ionization populating a state*/

for (mm =0; mm<nions; mm++)
	{
	for (nn=0;nn<nions;nn++)
		{
		if (mm==nn+1 && ion[nn].istate != ion[nn].z+1 && ion[mm].z==ion[nn].z && ion[nn].dere_di_flag > 0) 
			{
			rate_matrix[mm][mm]+=(xne*di_coeffs[nn]);
			}
		}
	}


/*Now we populate the elements relating to radiative recomb depopulating a state*/
		
for (mm = 0; mm<nions; mm++)
	{
	if (ion[mm].istate != 1) //we have space for electrons
		{
		rate_matrix[mm][mm]-=(xne*rr_rates[mm]);
		}
	}


/*Now we populate the elements relating to radiative recomb populating a state*/

for (mm =0; mm<nions; mm++)
	{
	for (nn=0;nn<nions;nn++)
		{
		if (mm==nn-1 && ion[nn].istate != 1 && ion[mm].z==ion[nn].z) 
			{
			rate_matrix[mm][nn]+=(xne*rr_rates[nn]);
			}
		}
	}

/*Now we populate the elements relating to dielectronic recombination depopulating a state*/
	

	
for (mm = 0; mm<nions; mm++)
	{
	if (ion[mm].istate != 1 && ion[mm].drflag > 0) //we have space for electrons
		{
		rate_matrix[mm][mm]-=(xne*dr_coeffs[mm]);
		}
	}


/*Now we populate the elements relating to dielectronic recombination populating a state*/



for (mm =0; mm<nions; mm++)
	{
	for (nn=0;nn<nions;nn++)
		{
		if (mm==nn-1 && ion[nn].istate != 1 && ion[mm].z==ion[nn].z && ion[mm].drflag > 0) 
			{
			rate_matrix[mm][nn]+=(xne*dr_coeffs[nn]);
			}
		}
	}






/*Now, we replace the first line for each element with 1's and 0's */

for (nn=0;nn<nions;nn++)
	{
	if (ion[nn].istate==1)
		{
		b_temp[nn]=nh * ele[xelem[nn]].abun;
		for (mm=0;mm<nions;mm++)
			{
			if (ion[mm].z==ion[nn].z)
				{
				rate_matrix[nn][mm]=1.0;
				}
			else
				{
				rate_matrix[nn][mm]=0.0;
				}
			}
		}
	else
		{
		b_temp[nn]=0.0;
		}
	}

//printf ("Populated - about to solve\n");

nrows=nions; //This is a placeholder, we may end up removing rows and columns

/*
printf ("About to solve the equation, M xne=%e\n",xne);
for (mm=0;mm<nrows;mm++)
	{
	for (nn=0;nn<nrows;nn++)
		{
		printf ("%e ",rate_matrix[mm][nn]);
		}
	printf ("\n");
	}
printf("\n");
*/

	  /********************************************************************************/
	  /* The block that follows (down to next line of ***s) is to do the
	     matrix inversion. It uses LU decomposition - the code for doing this is 
	     taken from the GSL manual with very few modifications. */
	  /* here we solve the matrix equation M x = b, where x is our vector containing
	     level populations as a fraction w.r.t the whole element */

	  /* Replaced inline array allocaation with calloc, which will work with older version of c compilers */


	  a_data =
	    (double *) calloc (sizeof (double), nrows * nrows);

	  for (mm = 0; mm < nrows; mm++)
	    {
	      for (nn = 0; nn < nrows; nn++)
		{
		  a_data[mm * nrows + nn] = rate_matrix[mm][nn];
		}
	    }



	  /* Replaced inline array allocaation with calloc, which will work with older version of c compilers 
	     calloc also sets the elements to zero, which is required */

	  b_data = (double *) calloc (sizeof (double), nrows);

/*We now put our b_data array into this array*/

for (nn=0;nn<nrows;nn++)
	{
	b_data[nn]=b_temp[nn];
//	printf ("b_data %i = %e \n",nn,b_data[nn]);
	}

//printf("populated solver matrices\n");


      /* create gsl matrix/vector views of the arrays of rates */
	  m = gsl_matrix_view_array (a_data, nrows, nrows);	//KSLNEW

	  /* these are used for testing the solution below */
	  test_matrix = gsl_matrix_alloc(nrows,nrows);
	  test_vector = gsl_vector_alloc (nrows);

	  gsl_matrix_memcpy(test_matrix, &m.matrix); 	// create copy for testing 

	  b = gsl_vector_view_array (b_data, nrows);	//KSLNEW

	  /* the populations vector will be a gsl vector which stores populations */
	  populations = gsl_vector_alloc (nrows);	//KSLNEW


	  p = gsl_permutation_alloc (nrows);	//NEWKSL

  
  int s;

  gsl_permutation * p = gsl_permutation_alloc (nrows);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, populations);

  gsl_permutation_free (p);

	  /**********************************************************/

      /* JM 140414 -- before we clean, we should check that the populations vector we 
         have just created really is a solution to the matrix equation */
      
      /* declaration contained in my_linalg.h, taken from gsl library */
      ierr = gsl_blas_dgemv(CblasNoTrans, 1.0, test_matrix, populations, 0.0, test_vector);

      if (ierr != 0)
        {
      	  Error_silent("matrix_solv: bad return when testing matrix solution to rate equations.\n");
        }

      /* now cycle through and check the solution to y = m * populations really is
         (1, 0, 0 ... 0) */
//printf ("\n");
      for (mm = 0; mm < nrows; mm++)
        {

          /* get the element of the vector we want to check */
          test_val = gsl_vector_get (test_vector, mm);
 //        printf ("test_val=%e\n",test_val);
  //       if ( (fabs((test_val - b_temp[mm]))/test_val) > EPSILON)
  //    	      Error_silent("matrix_solv: test solution fails for row %i %e != %e\n",
  //    	      	     mm,test_val,b_temp[mm]);    
          }     

        /* free memory */
        free (a_data);	
	    free (b_data);	
	    free (test_vector);
	    free (test_matrix);

//   for (nn=0;nn<nelements;nn++)
//	{
//        xden[nn]=0.0;
//	}

for (nn=0;nn<nions;nn++)
	{
	newden[nn]=0.0;
	for (mm=0;mm<nrows;mm++)
		{
		if (xion[mm]==nn)
			{
			newden[nn]=(gsl_vector_get (populations, mm));
			}
		}
	if (newden[nn] < DENSITY_MIN)
		newden[nn] = DENSITY_MIN;
	}


      xnew = get_ne (newden);	/* determine the electron density for this density distribution */
//printf ("new ne=%e\n",xne);
//printf ("new electron density is %e old is %e\n",xnew,xne);

     if (xnew < DENSITY_MIN)
	xnew = DENSITY_MIN;	/* fudge to keep a floor on ne */
      if (fabs ((xne - xnew) / (xnew)) < FRACTIONAL_ERROR || xnew < 1.e-6)
	{
	  break;
	}
      xne = xxxne = (xnew + xne) / 2.;	/*New value of ne */

      niterate++;


      if (niterate == MAXITERATIONS)
	{
	  Error
	    ("matrix_solv: failed to converge for cell %i t %e nh %e xnew %e\n",
	     xplasma->nplasma, t_e, nh, xnew);
 for (nn = 0; nn < geo.nxfreq; nn++)
    {
      Log ("numin= %e (%e) numax= %e (%e) Model= %2d PL_log_w= %e PL_alpha= %e Exp_w= %e EXP_temp= %e\n",xplasma-> fmin_mod[nn],geo.xfreq[nn],xplasma->fmax_mod[nn],geo.xfreq[nn+1],xplasma->spec_mod_type[nn],xplasma->pl_log_w[nn],xplasma->pl_alpha[nn],xplasma->exp_w[nn],xplasma->exp_temp[nn]);
}    




	  Error ("matrix_solv: xxne %e theta %e\n", xxne, theta);
	  return (-1);
	}
    }


  xplasma->ne = xnew;
  for (nn = 0; nn < nions; nn++)
    {
      /* If statement added here to suppress interference with macro populations (SS Apr 04) */
      if (ion[nn].macro_info == 0 || geo.macro_ioniz_mode == 0
	  || geo.macro_simple == 1)
	{
	  xplasma->density[nn] = newden[nn];
	}
    }







/*
for (nn=0;nn>nions;nn++)
	{
	nelem=

	  a = nh * ele[nelem].abun / sum;	//the scaling factor to get the overall abundance right
	  for (nion = first; nion < last; nion++)
	    {
	      newden[nion] *= a;	//apply scaling
	      if (sane_check (newden[nion]))	//check nothing has gone crazy
		Error
		  ("variable_temperature:sane check failed for density newden=%e, for ion=%i\n",
		   newden[nion], nion);
	    }

*/


return(0);
}



