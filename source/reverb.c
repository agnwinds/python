/* 
 * The subroutines in this file handle outputting spectrum array reverb data
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"
#include "python.h"

/***********************************************************
                                    Space Telescope Science Institute

Synopsis:
	delay_spectrum_summary(filename,mode,nspecmin,nspecmax,select_spectype,renorm,loglin)  
		writes out the spectrum to a file 

Arguments:		

	char filename[]		The name of the file to write
	char mode[];		The mode in which the file should be opened, usually 'w', but
					it could be 'a'
	int nspecmin,nspecmax	These two numbers define the spectra you want to write.    
	int select_spectype	The type of spectral file you want to create, 
					0 = raw, 1= flambda,2=fnu 
        double	renorm		This is renormalization which incrementally decreases to
				one as the detailed spectral calculation goes forward.  It
				was added to allow one to print out the spectrum at the
				end of each cycle, rather than the end of the entire
				calculation.
	char loglin[]		Are we outputting a log or a linear spectrum
				

Returns:
  
Description:	
	This simple routine simply writes the spectra to a file in an easily interpretable
	ascii format. Normally one would write all of the spectra in one go, but  one can use 
	spectrum summary to write various spectra to various files by using the variables
	nspecmin and nspecmax.. 
			
Notes:

History:
	14aug	-	Cloned from spectrum_summary by SWM
**************************************************************/
int 
delay_spectrum_summary (
    char filename[],
    char mode[],
    int nspecmin,
    int nspecmax,
    int select_spectype,
    double renorm,				/* parameter used to rescale spectrum as it is building up */
    int loglin					/* switch to tell the code if we are outputting a log or a lin */
	)
{
	FILE *fopen(), *fptr;
	int i, n;
	char string[LINELENGTH];
	double freq, freqmin, dfreq, freq1;
	double lfreqmin, lfreqmax, ldfreq;
	double x, dd;

	/* 
	 * Open or reopen a file for writing the spectrum 
	 */
	if ((fptr = fopen(filename, "w")) == NULL)
	{
		Error("delay_spectrum_summary: Unable to open %s for writing\n", filename);
		exit(0);
	}

	/* 
	 * Check that nspecmin and nspecmax are reasonable 
	 */
	if (nspecmin < 0 || nspecmax < 0 || nspecmin > nspecmax)
	{
		Error("delay_spectrum_summary: nspecmin %d or nspecmax %d not reasonable \n", nspecmin, nspecmax);
		exit(0);
	}

	/* 
	 * Construct and write a header string for the output file 
	 */
	fprintf(fptr, "# Python Version %s\n", VERSION);

	get_time(string);
	fprintf(fptr, "# Date	%s\n#  \n", string);

	/* 
	 * Write the rest of the header for the spectrum file 
	 */

	fprintf(fptr, "# \n# Freq.        Lambda");
	for (n = nspecmin; n <= nspecmax; n++)
	{
		fprintf(fptr, "   %8s", xxspec[n].name);
	}
	fprintf(fptr, " Weight\n");

	/* 
	 * Don't print out the end bins because they include all photons outside the frequency range and there may be some as a result of the fact that the bb function generate some IR photons 
	 */
	dd = 4. * PI * (100. * PC) * (100. * PC);

	if (loglin == 0)
	{
		freqmin = xxspec[nspecmin].freqmin;
		dfreq = (xxspec[nspecmin].freqmax - freqmin) / NWAVE;
		for (i = 1; i < NWAVE - 1; i++)
		{
			freq = freqmin + i * dfreq;
			fprintf(fptr, "%-8e %.3f ", freq, C * 1e8 / freq);
			for (n = nspecmin; n <= nspecmax; n++)
			{
				if(xxspec[n].delay[i] == 0.) x = 0;
				else x = xxspec[n].delay[i] / xxspec[n].delay_weight[i];
				fprintf(fptr, " %10.5g %10.5g", x * renorm, xxspec[n].f[i]);
			}
			fprintf(fptr, "\n");
		}
	}
	else if (loglin == 1)
	{
		lfreqmin = log10(xxspec[nspecmin].freqmin);
		freq1 = lfreqmin;
		lfreqmax = log10(xxspec[nspecmin].freqmax);
		ldfreq = (lfreqmax - lfreqmin) / NWAVE;

		for (i = 1; i < NWAVE - 1; i++)
		{
			freq = pow(10., (lfreqmin + i * ldfreq));
			dfreq = freq - freq1;
			fprintf(fptr, "%-8e %.3f ", freq, C * 1e8 / freq);
			for (n = nspecmin; n <= nspecmax; n++)
			{
				if(xxspec[n].delay[i] == 0.) x = 0;
				else x = xxspec[n].delay[i] / xxspec[n].delay_weight[i];
				fprintf(fptr, " %10.5g", x * renorm);	/* this really shouldn't get called if we are outputting log data */
			}


			fprintf(fptr, "\n");
			freq1 = freq;
		}
	}
	fclose(fptr);
	return (0);
}


/***********************************************************
Synopsis:
	delay_dump(filename, p, f1, f2, nspec)  
		Calculates the delay for the photons. Cloned from spectra.c

Arguments:		
	filename		File to output to
	PhotPtr p		Array of photons (pass global p)
	int nspec		Observer spectrum to dump for

Returns:
	To file only
  
Description:	
	Prep wipes clean the output file and writes the header, dump then
	builds up the output file in batches.

Notes:

History:
	14aug	-	Written by SWM
***********************************************************/
char delay_dump_file[LINELENGTH];
int delay_dump_bank_size=65535, delay_dump_bank_curr=0,delay_dump_spec;
PhotPtr delay_dump_bank;
int *delay_dump_bank_ex;

double
delay_to_observer(PhotPtr pp)
{
	plane_dummy observer;
	observer.x[0]   = pp->lmn[0]*geo.rmax;
	observer.x[1]   = pp->lmn[1]*geo.rmax;
	observer.x[2]   = pp->lmn[2]*geo.rmax;
	observer.lmn[0] = pp->lmn[0];
	observer.lmn[1] = pp->lmn[1];
	observer.lmn[2] = pp->lmn[2];
	return(pp->path + ds_to_plane(&observer,pp));
}

int 
delay_dump_prep (char filename[], int nspec, int restart_stat, int i_rank)
{
	FILE *fopen(), *fptr;
	char string[LINELENGTH], c_file[LINELENGTH], c_rank[LINELENGTH];
	int i;

	if(restart_stat) return(0);
	//Allocate and zero dump files and set extract status
	delay_dump_spec = nspec;
	delay_dump_bank = (PhotPtr) calloc(sizeof(p_dummy), delay_dump_bank_size);
	delay_dump_bank_ex = (int*) calloc(sizeof(int),		delay_dump_bank_size);
	for(i=0;i<delay_dump_bank_size;i++)	delay_dump_bank_ex[i] = 0;

	//Get output filename
	strcpy(c_file, filename);			//Copy filename to new string
	if(i_rank > 0)					
	{
		sprintf(c_rank,"%i",i_rank);	//Write thread # to string
		strcat(c_file,c_rank);			//Append thread # to filename
	}
	strcpy(delay_dump_file,c_file);		//Store modified filename for later

	/* If this isn't a continue run, prep the output file */
	if ((fptr = fopen(delay_dump_file, "w")) == NULL)
	{
		Error("delay_dump_prep: Unable to open %s for writing\n", delay_dump_file);
		exit(0);
	}

	if (nspec < MSPEC)					//Check that NSPEC is reasonable
	{
		Error("delay_dump_prep: nspec %d below MSPEC value \n", nspec);
		exit(0);
	}

	if(i_rank > 0) 
	{
		fprintf(fptr, "# Delay dump file for slave process %d\n",i_rank);
	}
	else 
	{									// Construct and write a header string for the output file 
		fprintf(fptr, "# Python Version %s\n", VERSION);	

		get_time(string);
		fprintf(fptr, "# Date	%s\n#  \n", string);
		fprintf(fptr, "# Spectrum: %s\n", xxspec[nspec].name);

		/* 
		 * Write the rest of the header for the spectrum file 
	 	*/
		fprintf(fptr, "# \n# Freq      Wavelength  Weight   "
				"  Last X       Last Y       Last Z     "
				"  Last L       Last M       Last N     "
				"Scatters RScatter Delay Extracted Spectrum\n");	
	}
	fclose(fptr);
	return(0);
}

/***********************************************************
Synopsis:
	delay_dump_finish()  
		Outputs any leftover photons.

Arguments:		
	None

Returns:
	To file only
  
Description:	
	Dumps last remaining extracted photons to file and frees arrays.

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int 
delay_dump_finish()
{
	if(delay_dump_bank_curr>0)
	{
		delay_dump(delay_dump_bank,delay_dump_bank_curr-1,delay_dump_spec, 1);
	}
	free(delay_dump_bank);
	free(delay_dump_bank_ex);
	return(0);
}

/***********************************************************
Synopsis:
	delay_dump_combine(iRanks)  
		Finishes off delay dump for multithreaded runs

Arguments:		
	int iRanks	Combines delay dump file from each thread

Returns:
	To file only
  
Description:	
	Uses system commands to achieve this quickly

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int 
delay_dump_combine(int iRanks)
{
	FILE *fopen(), *fptr;
	char string[LINELENGTH], cCall[LINELENGTH];
	
	//Yes this is done as a system call and won't work on Windows machines. Lazy solution!
	sprintf(cCall, "cat %s[0-9]* >> ", delay_dump_file);
	if(system(cCall)<0)
		Error("delay_dump_combine: Error calling system command '%s'",cCall);
	return(0);
}

/***********************************************************
Synopsis:
	delay_dump(PhotPtr p, int np, int nspec, int iExtracted)  
		Dumps array of photons provided to file

Arguments:		
	PhotPtr p 		Array of photons to dump
	int np 			Length of array
	int nspec 		Spectrum to dump for
	int iExtracted 	Whether or not this array is of extracted photons

Returns:
	To file only
  
Description:	
	Dumps to this thread's delay_dump_file the passed photon
	array.

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int 
delay_dump (PhotPtr p, int np, int nspec, int iExtracted)
{
	FILE *fopen(), *fptr;
	int nphot, mscat, mtopbot, i;
	double zangle;

	/* 
	 * Open a file for writing the spectrum 
	 */	
	if ((fptr = fopen(delay_dump_file, "a")) == NULL)
	{
		Error("delay_dump: Unable to reopen %s for writing\n", delay_dump_file);
		exit(0);
	}
	
	for (nphot = 0; nphot < np; nphot++)
	{
		if (p[nphot].istat == P_ESCAPE && p[nphot].nscat>0)
		{
			zangle = fabs(p[nphot].lmn[2]);
			/* 
			 * Complicated if statement to allow one to choose whether to construct the spectrum from all photons or just from photons which have scattered a specific number of times.  01apr13--ksl-Modified if statement to change behavior on negative numbers to say that a negative number for mscat implies that you accept any photon with |mscat| or more scatters 
			 */
			if (((mscat = xxspec[nspec].nscat) > 999 || p[nphot].nscat == mscat 
				|| (mscat < 0 && p[nphot].nscat >= (-mscat)))
				&& ((mtopbot = xxspec[nspec].top_bot) == 0 || (mtopbot * p[nphot].x[2]) > 0))
			{
				if (xxspec[nspec].mmin < zangle 
					&& zangle < xxspec[nspec].mmax)
				{	/* SWM 15/8/14 - Added path delay in comparison to photon heading straight from origin to rmax*/
					fprintf(fptr, "%10.5g %10.5g %10.5g %+10.5e %+10.5e %+10.5e %+10.5g %+10.5g %+10.5g %3d     %3d     %10.5g %5d %5d\n", 
						p[nphot].freq, C * 1e8 / p[nphot].freq, p[nphot].w, 
						p[nphot].x[0], p[nphot].x[1], p[nphot].x[2], 
						p[nphot].lmn[0], p[nphot].lmn[1], p[nphot].lmn[2], 
						p[nphot].nscat, p[nphot].nrscat, delay_to_observer(&p[nphot]),
						(iExtracted ? delay_dump_bank_ex[nphot] : 0), nspec);
				}
			}
		}
	}
	fclose(fptr);
	return (0);
}

/***********************************************************
Synopsis:
	delay_dump_single(PhotPtr pp, int extract_phot)  
		Finishes off delay dump for multithreaded runs

Arguments:		
	PhotPtr pp 			Photon to dump
	int extract_phot	Whether or not this was extracted

Returns:
	To internal arrays only
  
Description:	
	Pushes the passed photon into the delay_dump_bank. If
	this fills it, dumps current bank to file and restarts it.

Notes:

History:
	6/2/15	-	Written by SWM
***********************************************************/
int 
delay_dump_single (PhotPtr pp, int extract_phot)
{
	stuff_phot(pp, &delay_dump_bank[delay_dump_bank_curr]);			//Bank single photon in temp array
	delay_dump_bank_ex[delay_dump_bank_curr] = extract_phot;		//Record if it's extract photon
	if(delay_dump_bank_curr == delay_dump_bank_size-1)				//If temp array is full
	{																	//Dump to file, zero array position
		delay_dump(delay_dump_bank,delay_dump_bank_size,delay_dump_spec,1);
		delay_dump_bank_curr=0;
		int i;
		for(i=0;i<delay_dump_bank_size;i++) delay_dump_bank_ex[i]=0;	//Zero extract status of single array
	}
	else
	{
		delay_dump_bank_curr++;
	}
	return(0);
}

/***********************************************************
Synopsis:
	path_data_constructor(double r_rad_min, double r_rad_max,
						  int i_bins, int i_angles)  
		Returns global data container for wind path info

Arguments:		
	r_rad_min 		Minimum disk radius
	r_rad_max 		Maximum disk radius, used for making 
					the path distance arrays
	i_bins 			Number of path distance bins
	i_angles 		Number of observer angles

Returns:
	Pointer to prepped global path data info.
  
Description:	
	Sets up global data used by the individual cells for
	binning purposes.

Notes:

History:
	9/2/15	-	Written by SWM
***********************************************************/
Path_Data_Ptr
path_data_constructor (double r_rad_min, double r_rad_max, int i_bins, int i_angles)
{
	int i;
	Path_Data_Ptr data = (Path_Data_Ptr) calloc(sizeof(path_data_dummy),1);

	if(data == NULL)
	{
		Error("path_data_constructor: Could not allocate memory\n");
		exit(0);
	}

	data->i_obs = i_angles;
	data->i_path_bins=i_bins;
	data->ad_path_bin = (double*) calloc(sizeof(double),i_bins);
	for(i=0; i <= i_bins; i++){
		data->ad_path_bin[i] = r_rad_min + i*(r_rad_max*5.0-r_rad_min)/i_bins;
	}
	return(data);
}

/***********************************************************
Synopsis:
	wind_paths_constructor(WindPtr wind)  
		Constructs a wind path record and returns it.

Arguments:		
	Wind_Ptr wind	Wind cell the paths are being made for

Returns:
  	Pointer to finished wind paths

Description:	
	Constructs a wind cell's path data and returns a pointer
	to it. Also calculates the distance to each extract 
	viewpoint from the cell via a dummy photon.

Notes:

History:
	9/2/15	-	Written by SWM
***********************************************************/
Wind_Paths_Side_Ptr
wind_paths_side_constructor (WindPtr wind, int i_side)
{
	int i;
	PhotPtr p_test = (PhotPtr) calloc(sizeof(p_dummy), 1);
	Wind_Paths_Side_Ptr side = (Wind_Paths_Side_Ptr) calloc(sizeof(wind_paths_side_dummy), 1);

	if(side == NULL)
	{
		Error("wind_paths_side_constructor: Could not allocate memory for cell %d\n",wind->nwind);
		exit(0);
	}

	stuff_v(wind->xcen,p_test->x);
	p_test->x[0] *= i_side;
	side->ad_path_to_obs = (double*) calloc(sizeof(double), g_path_data->i_obs);

	for(i=0; i<g_path_data->i_obs;i++)
	{
		stuff_v(xxspec[MSPEC+i].lmn,p_test->lmn);
		side->ad_path_to_obs[i] = delay_to_observer(p_test);
	}	

	side->ad_freq_flux = (double*) calloc(sizeof(double), NWAVE-1);
	side->ad_freq_path_flux = (double*) calloc(sizeof(double), (NWAVE-1) * g_path_data->i_path_bins);

	if(side->ad_freq_path_flux == NULL)
	{
		Error("wind_paths_side_constructor: Could not allocate memory for cell %d bins\n",wind->nwind);
		exit(0);
	}
	free(p_test);
	return(side);
}

Wind_Paths_Ptr
wind_paths_constructor (WindPtr wind)
{
	Wind_Paths_Ptr paths = (Wind_Paths_Ptr) calloc(sizeof(wind_paths_dummy),1);

	if(paths == NULL)
	{
		Error("wind_paths_constructor: Could not allocate memory\n");
		exit(0);
	}
	paths->front = (Wind_Paths_Side_Ptr) wind_paths_side_constructor(wind,  1);
	paths->back	 = (Wind_Paths_Side_Ptr) wind_paths_side_constructor(wind, -1);
	return(paths);
}

/***********************************************************
Synopsis:
	wind_paths_add_phot(Wind_Paths_Ptr paths, PhotPtr pp)  
		Adds a photon's weight to a wind cell's delay array.

Arguments:		
	PhotPtr pp 				Photon to dump
	Wind_Paths_Ptr paths	Paths (from wind) to add to

Returns:
  
Description:	
	Adds photon to frequency and path bin as appropriate.
	As wind array is 2d, .5th dimension added via keeping
	two path sets for + and - x-axis positions). No full 3-d
	treatment required as Extract mode's viewpoint always
	lies along X-axis.

Notes:

History:
	9/2/15	-	Written by SWM
***********************************************************/
int
wind_paths_side_add_phot (Wind_Paths_Side_Ptr side, PhotPtr pp)
{
	int i,j;

	printf("wind_paths_side_add_phot: Photon at %g %g %g\n",
			pp->x[0], pp->x[1], pp->x[2]);

	for(i=0; i<NWAVE-1; i++)
	{
		if(pp->freq >= xxspec[0].f[i] && pp->freq <= xxspec[0].f[i+1])
		{
			for(j=0; j< g_path_data->i_path_bins; j++)
			{
				if(	pp->path >= g_path_data->ad_path_bin[i] && 
					pp->path <= g_path_data->ad_path_bin[i+1])
				{
					side->ad_freq_path_flux[i*(NWAVE-1)+j] += pp->w;
				}
			}
		}
	}
	return(0);
}
int
wind_paths_add_phot (WindPtr wind, PhotPtr pp)
{
	if(pp->x[0] >= 0.0)
		wind_paths_side_add_phot(wind->paths->front, pp);
	else
		wind_paths_side_add_phot(wind->paths->back,  pp);
	return(0);
}

/***********************************************************
Synopsis:
	wind_paths_add_phot(Wind_Paths_Ptr paths, PhotPtr pp)  
		Adds a photon's weight to a wind cell's delay array.

Arguments:		
	PhotPtr pp 				Photon to dump
	Wind_Paths_Ptr paths	Paths (from wind) to add to

Returns:
  
Description:	
	Adds photon to frequency and path bin as appropriate.
	As wind array is 2d, .5th dimension added via keeping
	two path sets for + and - x-axis positions). No full 3-d
	treatment required as Extract mode's viewpoint always
	lies along X-axis.

Notes:

History:
	10/2/15	-	Written by SWM
***********************************************************/
int
wind_paths_side_evaluate (Wind_Paths_Side_Ptr side)
{
	int i,j;

	side->d_flux  = 0.0;
	side->d_path  = 0.0;
	for(i=0; i<NWAVE-1; i++)
	{
		side->ad_freq_flux[i]=0.0;
		for(j=0; j< g_path_data->i_path_bins; j++)
		{
			side->d_path			+= side->ad_freq_path_flux[i*(NWAVE)+j]*
									(g_path_data->ad_path_bin[i] + g_path_data->ad_path_bin[i+1]) / 2.0;
			side->ad_freq_flux[i] 	+= side->ad_freq_path_flux[i*(NWAVE)+j];
		}
		side->d_flux += side->ad_freq_flux[i];
	}
	if(side->d_flux > 0.0) side->d_path /= side->d_flux;
	return(0);
}
int
wind_paths_evaluate(WindPtr wind)
{
	int i;
	for(i=0; i<geo.ndim*geo.mdim; i++)
	{
		if(wind[i].inwind) 
		{
			wind_paths_side_evaluate (wind[i].paths->front);
			wind_paths_side_evaluate (wind[i].paths->back);
		}
	}
	return(0);
}


int
wind_paths_init(WindPtr wind)
{
	int i;

	for(i=0; i<geo.ndim*geo.mdim; i++)
	{
		if(wind[i].inwind)wind_paths_constructor (&wind[i]);
	}
	return(0);
}