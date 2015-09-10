/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
 
Arguments:		

Returns:
 
Description:	
		
Notes:


History:

	15sep 	ksl	Moved calculating the detailed spectra to a separat routine
**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"


#include "python.h"

int make_spectra(restart_stat)
	int restart_stat;
{
  WindPtr w;
  PhotPtr p;

  double freqmin, freqmax;
  double renorm;
  long nphot_to_define;
  int iwind;

  /* Next three lines have variables that should be a structure, or possibly we
     should allocate the space for the spectra to avoid all this nonsense.  02feb ksl */


  int icheck;



/* XXXX - THE CALCULATION OF A DETAILED SPECTRUM IN A SPECIFIC REGION OF WAVELENGTH SPACE */

  p=photmain;
  w=wmain;

  freqmax = C / (geo.swavemin * 1.e-8);
  freqmin = C / (geo.swavemax * 1.e-8);


  /* Perform the initilizations required to handle macro-atoms during the detailed
     calculation of the spectrum.  

     Next lines turns off macro atom estimators and other portions of the code that are
     unnecessary during spectrum cycles.  */

  geo.ioniz_or_extract = 0;

/* 57h -- 07jul -- Next steps to speed up extraction stage */
  if (!modes.keep_photoabs)
    {
      DENSITY_PHOT_MIN = -1.0;	// Do not calculated photoabsorption in detailed spectrum 
    }

  /*Switch on k-packet/macro atom emissivities  SS June 04 */

  if (geo.rt_mode == 2)
    {
      geo.matom_radiation = 1;
    }

  /* Finished initializations required for macro-atom approach */

  /* Calculate and store which bf processess need to be considered in each cell
   * Note that this is not macro-specific but is just to speed the program up.
   */

  kbf_need (freqmin, freqmax);

  /* XXXX - BEGIN CYCLES TO CREATE THE DETAILED SPECTRUM */

  /* the next section initializes the spectrum array in two cases, for the
   * standard one where one is calulating the spectrum for the first time
   * and in the somewhat abnormal case where additional ionization cycles
   * were calculated for the wind
   */

  if (geo.pcycle == 0)
    {
      spectrum_init (freqmin, freqmax, geo.nangles, geo.angle, geo.phase,
		     geo.scat_select, geo.top_bot_select,
		     geo.select_extract, geo.rho_select, geo.z_select,
		     geo.az_select, geo.r_select);

      /* 68b - zero the portion of plasma main that records the numbers of scatters by
       * each ion in a cell
       */

      zero_scatters ();

    }

  /* the next condition should really when one has nothing more to do */

  else if (geo.pcycle >= geo.pcycles)
    xsignal (files.root,
	     "%-20s No spectrum   needed: pcycles(%d)==pcycles(%d)\n",
	     "COMMENT", geo.pcycle, geo.pcycles);

  else
    {
      /* Then we are restarting a run with more spectral cycles, but we 
         have already completed some. The memory for the spectral arrays
         should already have been allocated, and the spectrum was initialised
         on the original run, so we just need to renormalise the saved spectrum */
      /* See issue #134 (JM) */
      if (restart_stat == 0)
	Error
	  ("Not restarting, but geo.pcycle = %i and trying to renormalise!\n",
	   geo.pcycle);

      spectrum_restart_renormalise (geo.nangles);
    }


  while (geo.pcycle < geo.pcycles)
    {				/* This allows you to build up photons in bunches */

      xsignal (files.root, "%-20s Starting %d of %d spectral cycle \n",
	       "NOK", geo.pcycle, geo.pcycles);

      if (modes.ispy)
	ispy_init ("python", geo.pcycle + 1000);


      Log ("!!Cycle %d of %d to calculate a detailed spectrum\n",
	   geo.pcycle, geo.pcycles);
      Log_flush ();		/*NSH June 13 Added call to flush logfile */
      if (!geo.wind_radiation)
	iwind = -1;		/* Do not generate photons from wind */
      else if (geo.pcycle == 0)
	iwind = 1;		/* Create wind photons and force a reinitialization of wind parms */
      else
	iwind = 0;		/* Create wind photons but do not force reinitialization */

      /* Create the initial photon bundles which need to be trannsported through the wind 

         For the detailed spectra, NPHOT*pcycles is the number of photon bundles which will equal the luminosity, 
         1 implies that detailed spectra, as opposed to the ionization of the wind is being calculated

         JM 130306 must convert NPHOT and pcycles to double precision variable nphot_to_define

       */

      nphot_to_define = (long) NPHOT *(long) geo.pcycles;
      define_phot (p, freqmin, freqmax, nphot_to_define, 1, iwind, 0);

      for (icheck = 0; icheck < NPHOT; icheck++)
	{
	  if (sane_check (p[icheck].freq))
	    {
	      Error
		("python after define phot:sane_check unnatural frequency for photon %d\n",
		 icheck);
	    }
	}


      /* Tranport photons through the wind */

      trans_phot (w, p, geo.select_extract);

      if (modes.print_windrad_summary)
	wind_rad_summary (w, files.windrad, "a");


      spectrum_create (p, freqmin, freqmax, geo.nangles, geo.select_extract);

/* Write out the detailed spectrum each cycle so that one can see the statistics build up! */
      renorm = ((double) (geo.pcycles)) / (geo.pcycle + 1.0);

      /* Do an MPI reduce to get the spectra all gathered to the master thread */
#ifdef MPI_ON
      gather_spectra_para (spec_spec_helpers, nspectra);
#endif





#ifdef MPI_ON
      if (rank_global == 0)
	{
#endif
	  spectrum_summary (files.spec, "w", 0, nspectra - 1,
			    geo.select_spectype, renorm, 0);
#ifdef MPI_ON
	}
#endif
      Log ("Completed spectrum cycle %3d :  The elapsed TIME was %f\n",
	   geo.pcycle, timer ());

      /* SWM0215: Delay dump photons from this cycle */
      if (geo.reverb > REV_NONE)
	delay_dump (p, NPHOT, 0);	// SWM - Dump delay tracks from this iteration

      /* JM1304: moved geo.pcycle++ after xsignal to record cycles correctly. First cycle is cycle 0. */

      xsignal (files.root, "%-20s Finished %3d of %3d spectrum cycles \n",
	       "OK", geo.pcycle, geo.pcycles);

      geo.pcycle++;		// Increment the spectral cycles

#ifdef MPI_ON
      if (rank_global == 0)
	{
#endif
	  wind_save (files.windsave);	// This is only needed to update pcycle
	  spec_save (files.specsave);
#ifdef MPI_ON
	}
#endif
      check_time (files.root);
    }


/* XXXX -- END CYCLE TO CALCULATE DETAILED SPECTRUM */

  phot_gen_sum (files.phot, "a");

/* 57h - 07jul -- ksl -- Write out the freebound information */

#ifdef MPI_ON
  if (rank_global == 0)
    {
#endif
      fb_save ("recomb.save");
#ifdef MPI_ON
    }
#endif

  /* SWM0215: Dump the last photon path details to file */
  if (geo.reverb != REV_NONE)
    delay_dump_finish ();	// Each thread dumps to file
#ifdef MPI_ON
  MPI_Barrier (MPI_COMM_WORLD);	// Once all done
  if (my_rank == 0 && geo.reverb != REV_NONE)
    delay_dump_combine (np_mpi_global);	// Combine results if necessary
#endif


/* Finally done */

#ifdef MPI_ON
  sprintf (dummy, "End of program, Thread %d only", my_rank);	// added so we make clear these are just errors for thread ngit status    
  error_summary (dummy);	// Summarize the errors that were recorded by the program
  Log ("Run py_error.py for full error report.\n");
#else
  error_summary ("End of program");	// Summarize the errors that were recorded by the program
#endif


#ifdef MPI_ON
  MPI_Finalize ();
  Log_parallel ("Thread %d Finalized. All done\n", my_rank);
#endif


  xsignal (files.root, "%-20s %s\n", "COMPLETE", files.root);
  Log ("Completed entire program.  The elapsed TIME was %f\n", timer ());
  return EXIT_SUCCESS;
}
