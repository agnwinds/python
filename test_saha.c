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

double g1,g2,xne,xip;


#include "python.h"
int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;
  PhotPtr p;
  int n,i,j,nion;
  int ntemp,nne,temp_start;
  double ne[11],temp,temp1,nh,xtemp,weight;
  FILE *fp_h;
  FILE *fp_he;
  FILE *fp_c;
  FILE *fp_n;
  FILE *fp_o;
  FILE *fp1_o;
//  PlasmaPtr xplasma;

  temp_start=399;
  ntemp=301;
  nne=1;
  nh=1e10;
  weight=1e-8;

  strcpy (geo.atomic_filename, "atomic/standard39");
   printf("atomic_filename=%s\n",geo.atomic_filename);
  get_atomic_data (geo.atomic_filename);

	printf("density_min=%e\n",1.1*DENSITY_MIN);


  p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);
  NDIM = ndim = 1;
  MDIM = mdim = 1;
  NDIM2 = NDIM * MDIM;

  calloc_wind (NDIM2);
  w = wmain;

  NPLASMA = 1;
  calloc_plasma (NPLASMA);
//  xplasma=plasmamain;


  plasmamain[0].rho=nh/rho2nh;
  
  fp_h=fopen("hydrogen_saha.out","w");
  fp_he=fopen("helium_saha.out","w");
  fp_c=fopen("carbon_saha.out","w");
  fp_n=fopen("nitrogen_saha.out","w");
  fp_o=fopen("oxygen_saha.out","w");
  fp1_o=fopen("oxygen_partition.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);
	fprintf(fp1_o,"%i %i\n",ntemp,nne);

		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",ion[j].ip);
			}
		fprintf(fp_o,"\n");

	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1++;
		temp=pow(10,temp1/100);
		plasmamain[0].t_e=temp;
  		plasmamain[0].ne=nh;
		plasmamain[0].t_r=temp;
		plasmamain[0].w=1;
		partition_functions (&plasmamain[0], 0);
  		saha (&plasmamain[0],  nh, temp);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);
		fprintf(fp1_o,"%e %e",temp,plasmamain[0].ne);		
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			fprintf(fp1_o," %6.3e",plasmamain[0].partition[j]);
			}
		fprintf(fp_o,"\n");
		fprintf(fp1_o,"\n");
		}

	
  fclose(fp_h);
  fclose(fp_he);
  fclose(fp_c);
  fclose(fp_n);
  fclose(fp_o);
  fclose(fp1_o);

/*
  fp_h=fopen("hydrogen_vt.out","w");
  fp_he=fopen("helium_vt.out","w");
  fp_c=fopen("carbon_vt.out","w");
  fp_n=fopen("nitrogen_vt.out","w");
  fp_o=fopen("oxygen_vt.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);

	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1=temp1+1;
		temp=pow(10,temp1/100);
		printf ("Temperature = %e\n",temp);
//		temp=1.380384e5;
		plasmamain[0].t_e=temp;
		plasmamain[0].t_r=temp;
		plasmamain[0].w=weight;
	//	printf("Going to the new code\n");
  		nebular_concentrations (&plasmamain[0], 6);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		}

  fclose(fp_h);
  fclose(fp_he);
  fclose(fp_c);
  fclose(fp_n);
  fclose(fp_o);*/


  fp_h=fopen("hydrogen_nc.out","w");
  fp_he=fopen("helium_nc.out","w");
  fp_c=fopen("carbon_nc.out","w");
  fp_n=fopen("nitrogen_nc.out","w");
  fp_o=fopen("oxygen_nc.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);
 


	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1++;
		temp=pow(10,temp1/100);
		plasmamain[0].t_e=temp;
		plasmamain[0].t_r=temp;
  		nebular_concentrations (&plasmamain[0], 0);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);		
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		}

  fclose(fp_h);
  fclose(fp_he);
  fclose(fp_c);
  fclose(fp_n);
  fclose(fp_o);


  fp_h=fopen("hydrogen_lm.out","w");
  fp_he=fopen("helium_lm.out","w");
  fp_c=fopen("carbon_lm.out","w");
  fp_n=fopen("nitrogen_lm.out","w");
  fp_o=fopen("oxygen_lm.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);

	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1=temp1+1;
		temp=pow(10,temp1/100);
//		temp=1.380384e5;
		plasmamain[0].t_e=temp;
		plasmamain[0].t_r=temp;
		plasmamain[0].w=weight;
  		nebular_concentrations (&plasmamain[0], 2);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		}

  fclose(fp_h);
  fclose(fp_he);
  fclose(fp_c);
  fclose(fp_n);
  fclose(fp_o);


  fp_h=fopen("hydrogen_pl.out","w");
  fp_he=fopen("helium_pl.out","w");
  fp_c=fopen("carbon_pl.out","w");
  fp_n=fopen("nitrogen_pl.out","w");
  fp_o=fopen("oxygen_pl.out","w");
	fprintf(fp_h,"%i %i\n",ntemp,nne);
	fprintf(fp_he,"%i %i\n",ntemp,nne);
	fprintf(fp_c,"%i %i\n",ntemp,nne);
	fprintf(fp_n,"%i %i\n",ntemp,nne);
	fprintf(fp_o,"%i %i\n",ntemp,nne);
 


	temp1=temp_start;
	for (i=0;i<ntemp;i++)
		{
		temp1=temp1+1;
		temp=pow(10,temp1/100);
		printf ("Temperature =%e\n",temp);
		geo.nxfreq=1;
		geo.xfreq[0]=1e14;
		geo.xfreq[1]=1e18;
		plasmamain[0].sim_alpha[0]=-1.5;
		plasmamain[0].sim_w[0]=1e15;
		xband.nbands=1;
    		xband.f1[0]=1e14;
		xband.f2[0]=1e18;
		plasmamain[0].t_e=temp;
		plasmamain[0].t_r=temp;
  		nebular_concentrations (&plasmamain[0], 7);
		fprintf(fp_h,"%6.3e %6.3e %6.3e %6.3e\n",temp,plasmamain[0].ne,plasmamain[0].density[0],plasmamain[0].density[1]);
		fprintf(fp_he,"%e %e",temp,plasmamain[0].ne);		
		for (j=2;j<5;j++ ) 
			{		
			fprintf(fp_he," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_he,"\n");
		fprintf(fp_c,"%e %e",temp,plasmamain[0].ne);		
		for (j=5;j<12;j++ ) 
			{		
			fprintf(fp_c," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_c,"\n");
		fprintf(fp_n,"%e %e",temp,plasmamain[0].ne);		
		for (j=12;j<20;j++ ) 
			{		
			fprintf(fp_n," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_n,"\n");
		fprintf(fp_o,"%e %e",temp,plasmamain[0].ne);		
		for (j=20;j<29;j++ ) 
			{		
			fprintf(fp_o," %6.3e",plasmamain[0].density[j]);
			}
		fprintf(fp_o,"\n");
		} 






 return EXIT_SUCCESS;
}


