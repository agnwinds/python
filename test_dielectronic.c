/***********************************************************
                                       Southampton university

 Synopsis:
	test_dielectronic is a stand alone code to produce data relating to trying to understand dielectronic recombination in all its forms!
 
Arguments:		


Returns:
 
Description:	
	
Notes:

History:
	 
	0512 - nsh wrote it
 	


**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
struct topbase_phot *fb_xtop;
double g1,g2,xne,xip;
#define SAHA 4.82907e15		
double xip,xxxne,qromb_temp;
double fbt;			// Temperature at which thee emissivity is calculated
int fbfr;
#include "python.h"
int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;
  PhotPtr p;
  int n,i,j,nion,ion_lower;
  int ntemp,nne,temp_start,tempdiv,ntmin;
  double ne[11],temp1,nh,xtemp,weight,xsaha,t_r,t_e,temp,zeta,alpha_all,alpha_gs,comp_zeta;
  double b,*partition,recomb_fudge,pi_fudge,gs_fudge,test,integral;
  double fthresh,fmax;
  double temp_func();
  int compute_dr_coeffs();
  FILE *fp_c_zeta;
  FILE *fp_c_die;
  FILE *fp_c_rec_all;
  FILE *fp_c_rec_gs;
  FILE *fp_c_comp_zeta;
  double qromb ();
  double fb_topbase_partial (); 
 
  PlasmaPtr xplasma;
  partition = xplasma->partition; /* Set the partition function array to that held for the cell */
  temp_start=299;
  ntemp=401;
  tempdiv=100.0;
  nne=1;
  nh=1e10;
  weight=0.01;
  xne=nh;

  fp_c_zeta=fopen("carbon_zeta.out","w");
  fp_c_die=fopen("carbon_die.out","w");
  fp_c_rec_all=fopen("carbon_rec_all.out","w");
  fp_c_rec_gs=fopen("carbon_rec_gs.out","w");
  fp_c_comp_zeta=fopen("carbon_comp_zeta.out","w");

  strcpy (geo.atomic_filename, "atomic/standard70");
  printf("atomic_filename=%s\n",geo.atomic_filename);
  get_atomic_data (geo.atomic_filename);


 //    init_freebound (1.e3, 1.e9, 0, 1e50);
	for (i=3;i<9;i++)
		{
		temp=pow(10,i/1.);
		compute_dr_coeffs (temp);

		fprintf(fp_c_zeta," %6.3e",temp);
		fprintf(fp_c_die," %6.3e",temp);
		fprintf(fp_c_rec_all," %6.3e",temp);
		fprintf(fp_c_rec_gs," %6.3e",temp);
		fprintf(fp_c_comp_zeta," %6.3e",temp);
 		fbt = temp;			/* Externally transmitted variable */
  		fbfr = 2;
		for (j=6;j<12;j++ ) 
			{
			printf ("Ion %i, is z=%i, istate=%i\n",j,ion[j].z,ion[j].istate);		
			zeta=compute_zeta (temp, j, 1e14, 1e19, 1);
     			alpha_all = xinteg_fb (temp, 1e14, 1e19, j-1, 2);
      			ntmin = ion[j-1].ntop_ground;
			printf ("For ion %i, ground topbase is %i\n",j-1,ntmin);
      			fb_xtop = &phot_top[ntmin];
			fthresh = fb_xtop->freq[0];
	  		fmax = fb_xtop->freq[fb_xtop->np - 1];	

			alpha_gs = qromb (fb_topbase_partial, fthresh,fmax, 1e-4);
			printf ("Integrating cross section %i from %e to %e and getting %e\n",ntmin,fthresh,fmax,alpha_gs);	
			comp_zeta=alpha_gs/alpha_all;
			fprintf(fp_c_zeta," %6.3e",zeta);
			fprintf(fp_c_die," %6.3e",dr_coeffs[j]);
			fprintf(fp_c_rec_all," %6.3e",alpha_all);
			fprintf(fp_c_rec_gs," %6.3e",alpha_gs);
			fprintf(fp_c_comp_zeta," %6.3e",comp_zeta);
			}
		fprintf(fp_c_rec_all,"\n");
		fprintf(fp_c_zeta,"\n");
		fprintf(fp_c_die,"\n");
		fprintf(fp_c_rec_gs,"\n");
		fprintf(fp_c_comp_zeta,"\n");
		}	





 return EXIT_SUCCESS;
}


