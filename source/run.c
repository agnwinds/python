
/***********************************************************/
/** @file  run.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief  the driving routines to carry out calculation of 
 * the ionization of the plasma and also to extract detailed 
 * spectra after the inputs have been collected.
 *
 * ### Programming Comment ###
 * The name of this file is not really acurate.  The routines
 * here do drive the major portions of the calculation but they
 * are still run from python.c.  It might be better to move even
 * more of the running of the code to here.  Alternatively, one
 * might make python.c simpler, so that developers could see
 * the structure better, but moving the input section
 * into it's own file.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      run the ionization cycles for a 
 * python model
 *
 * @param [in] int  restart_stat   0 if the is run is beginning from
 * scratch,  non-zero if this was a restart
 * @return     Always returns 0 
 *
 * @details
 * This is the main routine for running the ionization
 * cycles in Python
 *
 * ### Notes ###
 *
 **********************************************************/

int
calculate_ionization (restart_stat)
     int restart_stat;
{
  int n, nn;
  double zz, z_abs_all, z_abs[N_ISTAT], z_else, ztot;
  int nphot_istat[N_ISTAT];
  WindPtr w;
  PhotPtr p;

  char dummy[LINELENGTH];

  double freqmin, freqmax, x;
  long nphot_to_define, nphot_min;
  int iwind;


#ifdef MPI_ON
  int ioniz_spec_helpers;
#endif

  /* Save the the windfile before the first ionization cycle in order to
   * allow investigation of issues that may have arisen at the very beginning
   */

#ifdef MPI_ON
  if (rank_global == 0)
  {
#endif
    wind_save (files.windsave); // This is only needed to update pcycle
#ifdef MPI_ON
  }
#endif



  p = photmain;
  w = wmain;

  freqmin = xband.f1[0];
  freqmax = xband.f2[xband.nbands - 1];

#ifdef MPI_ON
  /* the length of the big arrays to help with the MPI reductions of the spectra
     the variables for the estimator arrays are set up in the subroutines themselves */
  ioniz_spec_helpers = 2 * MSPEC * NWAVE;       //we need space for log and lin spectra for MSPEC XNWAVE
#endif

/* THE CALCULATION OF THE IONIZATION OF THE WIND */

  geo.ioniz_or_extract = 1;     //SS July 04 - want to compute MC estimators during ionization cycles
  //1 simply implies we are in the ionization section of the code
  //and allows routines to act accordinaly.


/* BEGINNING OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */

  if (geo.wcycle == geo.wcycles)
    xsignal (files.root, "%-20s No ionization needed: wcycles(%d)==wcyeles(%d)\n", "COMMENT", geo.wcycle, geo.wcycles);
  else
  {
    geo.pcycle = 0;             /* Set the spectrum cycles executed to 0, because 
                                   we are going to modify the wind and hence any
                                   previously calculated spectra must be recreated
                                 */
  }

  /* SWM - Setup for path tracking */
  if (geo.reverb > REV_NONE)
  {
    reverb_init (wmain);
    delay_dump_prep (restart_stat);
  }


  while (geo.wcycle < geo.wcycles)
  {                             /* This allows you to build up photons in bunches */

    xsignal (files.root, "%-20s Starting %3d of %3d ionization cycles \n", "NOK", geo.wcycle + 1, geo.wcycles);

    Log ("!!Python: Beginning cycle %d of %d for defining wind\n", geo.wcycle + 1, geo.wcycles);
    Log_flush ();               /* Flush the log file (so that we know where are if there are problems */

    /* Initialize all of the arrays, etc, that need initialization for each cycle
     */

    spectrum_init (freqmin, freqmax, geo.nangles, geo.angle, geo.phase,
                   geo.scat_select, geo.top_bot_select, geo.select_extract, geo.rho_select, geo.z_select, geo.az_select, geo.r_select);

    wind_rad_init ();           /*Zero the parameters pertaining to the radiation field */


    geo.n_ioniz = 0.0;
    geo.cool_tot_ioniz = 0.0;
    ztot = 0.0;                 /* ztot is the luminosity of the disk multipled by the number of cycles, which is used by save_disk_heating */

    if (!geo.wind_radiation || (geo.wcycle == 0 && geo.run_type != RUN_TYPE_PREVIOUS))
      iwind = -1;               /* Do not generate photons from wind */
    else
      iwind = 1;                /* Create wind photons and force a reinitialization of wind parms */


    /* If we are using photon speed up mode then the number of photons varies by cycle in the
     * ionization phase.  We set this up here
     */

    if (modes.photon_speedup)
    {
      nphot_min = NPHOT_MAX / pow (10., PHOT_RANGE);

      x = log10 (NPHOT_MAX / nphot_min) / (geo.wcycles - 1);
      NPHOT = nphot_min * pow (10., (x * geo.wcycle));
      if (NPHOT > NPHOT_MAX)
      {
        NPHOT = NPHOT_MAX;
      }
    }

    Log ("!!Python: %1.2e photons will be transported for cycle %i\n", (double) NPHOT, geo.wcycle);

    /* Create the photons that need to be transported through the wind
     *
     * NPHOT is the number of photon bundles which will equal the luminosity; 
     * 0 => for ionization calculation 
     */

    nphot_to_define = (long) NPHOT;

    define_phot (p, freqmin, freqmax, nphot_to_define, 0, iwind, 1);
    photon_checks (p, freqmin, freqmax, "Check before transport");

    /* Zero the arrays, and other variables that need to be zeroed after the photons are generated. */


    geo.lum_star_back = 0;
    geo.lum_disk_back = 0;


    for (n = 0; n < NRINGS; n++)
    {
      qdisk.heat[n] = qdisk.nphot[n] = qdisk.w[n] = qdisk.ave_freq[n] = 0;
    }


    zz = 0.0;
    for (nn = 0; nn < NPHOT; nn++)
    {
      zz += p[nn].w;
    }

    Log ("!!python: Total photon luminosity before transphot %18.12e\n", zz);
    Log_flush ();
    ztot += zz;                 /* Total luminosity in all cycles, used for calculating disk heating */

    /* kbf_need determines how many & which bf processes one needs to considere.  It was introduced
     * as a way to speed up the program.  It has to be recalculated evey time one changes
     * freqmin and freqmax
     */

    kbf_need (freqmin, freqmax);

    /* NSH 22/10/12  This next call populates the prefactor for free free heating for each cell in the plasma array */
    /* NSH 4/12/12  Changed so it is only called if we have read in gsqrd data */
    if (gaunt_n_gsqrd > 0)
      pop_kappa_ff_array ();

    /* Transport the photons through the wind */
    trans_phot (w, p, 0);

    /* Determine how much energy was absorbed in the wind. first zero counters. 
       There are counters for total energy absorbed and for each entry in the istat enum */
    z_abs_all = z_else = 0.0;
    for (nn = 0; nn < N_ISTAT; nn++)
    {
      z_abs[nn] = 0.0;
      nphot_istat[nn] = 0.0;
    }

    /* loop over the different photon istats to determine where the luminosity went */
    for (nn = 0; nn < NPHOT; nn++)
    {
      z_abs_all += p[nn].w;

      /* we want the istat to be >1 (not P_SCAT or P_INWIND) */
      if (p[nn].istat < N_ISTAT && p[nn].istat > 1)
      {
        z_abs[p[nn].istat] += p[nn].w;
        nphot_istat[p[nn].istat]++;
      }
      else
        z_else += p[nn].w;
    }

    Log
      ("!!python: Total photon luminosity after transphot  %18.12e (absorbed or lost  %18.12e). Radiated luminosity %18.12e\n",
       z_abs_all, z_abs_all - zz, z_abs[P_ESCAPE]);
    if (geo.rt_mode == RT_MODE_MACRO)
    {
      Log ("!!python: luminosity lost by adiabatic kpkt destruction %18.12e number of packets %d\n", z_abs[P_ADIABATIC],
           nphot_istat[P_ADIABATIC]);
      Log ("!!python: luminosity lost to low-frequency free-free    %18.12e number of packets %d\n", z_abs[P_LOFREQ_FF],
           nphot_istat[P_LOFREQ_FF]);
    }
    Log ("!!python: luminosity lost by being completely absorbed  %18.12e \n", z_abs[P_ABSORB]);
    Log ("!!python: luminosity lost by too many scatters          %18.12e \n", z_abs[P_TOO_MANY_SCATTERS]);
    Log ("!!python: luminosity lost by hitting the central object %18.12e \n", z_abs[P_HIT_STAR]);
    Log ("!!python: luminosity lost by hitting the disk           %18.12e \n", z_abs[P_HIT_DISK]);
    Log ("!!python: luminosity lost by errors                     %18.12e \n",
         z_abs[P_ERROR] + z_abs[P_ERROR_MATOM] + z_abs[P_REPOSITION_ERROR]);
    Log ("!!python: luminosity lost by the unknown                %18.12e \n", z_else);
    if (geo.binary == TRUE)
      Log ("!!python: luminosity lost by hitting the secondary %18.12e \n", z_abs[P_SEC]);


    photon_checks (p, freqmin, freqmax, "Check after transport");
    spectrum_create (p, geo.nangles, geo.select_extract);



    /* At this point we should communicate all the useful infomation 
       that has been accummulated on differenet MPI tasks */

#ifdef MPI_ON

    communicate_estimators_para ();

    communicate_matom_estimators_para ();       // this will return 0 if nlevels_macro == 0
#endif



    /* Calculate and store the amount of heating of the disk due to radiation impinging on the disk */
    /* We only want one process to write to the file, and we only do this if there is a disk */

#ifdef MPI_ON
    if (rank_global == 0)
    {
#endif
      if (geo.disk_type != DISK_NONE)
        qdisk_save (files.disk, ztot);
#ifdef MPI_ON
    }
    MPI_Barrier (MPI_COMM_WORLD);
#endif

/* Completed writing file describing disk heating */

    Log ("!!python: Number of ionizing photons %g lum of ionizing photons %g\n", geo.n_ioniz, geo.cool_tot_ioniz);

/* Note that this step is parallelized */

    wind_update (w);


    Log ("Completed ionization cycle %d :  The elapsed TIME was %f\n", geo.wcycle + 1, timer ());

    /* Do an MPI reduce to get the spectra all gathered to the master thread */

#ifdef MPI_ON

    gather_spectra_para (ioniz_spec_helpers, MSPEC);

#endif



#ifdef MPI_ON
    if (rank_global == 0)
    {
#endif

/* The variables for spectrum_sumamry are the filename, the attribute for the file write, the minimum and maximum spectra to write out, 
 * the type of spectrum (RAW meaning internal luminosity units, the amount by which to renormalize (1 means use the existing
 * values, loglin (0=linear, 1=log for the wavelength scale), all photons or just wind photons
 */

      spectrum_summary (files.wspec, 0, 6, SPECTYPE_RAW, 1., 0, 0);     /* .spec_tot */
      spectrum_summary (files.lwspec, 0, 6, SPECTYPE_RAW, 1., 1, 0);    /* .log_spec_tot */
      spectrum_summary (files.wspec_wind, 0, 6, SPECTYPE_RAW, 1., 0, 1);        /* .spec_tot_wind  */
      spectrum_summary (files.lwspec_wind, 0, 6, SPECTYPE_RAW, 1., 1, 1);       /* .log_spec_tot_wind */
      phot_gen_sum (files.phot, "w");   /* Save info about the way photons are created and absorbed
                                           by the disk */
#ifdef MPI_ON
    }
    MPI_Barrier (MPI_COMM_WORLD);
#endif

    /* Save everything after each cycle and prepare for the next cycle 
       JM1304: moved geo.wcycle++ after xsignal to record cycles correctly. First cycle is cycle 0. */
    /* NSH1306 - moved geo.wcycle++ back, but moved the log and xsignal statements */


    xsignal (files.root, "%-20s Finished %3d of %3d ionization cycles \n", "OK", geo.wcycle + 1, geo.wcycles);
    geo.wcycle++;               //Increment ionisation cycles


/* Save only the windsave file from thread 0, to prevent many processors from writing to the same
 * file. */

#ifdef MPI_ON
    if (rank_global == 0)
    {
#endif
      wind_save (files.windsave);
      Log_silent ("Saved wind structure in %s after cycle %d\n", files.windsave, geo.wcycle);

      /* In a diagnostic mode save the wind file for each cycle (from thread 0) */

      if (modes.keep_ioncycle_windsaves)
      {
        strcpy (dummy, "");
        sprintf (dummy, "python%02d.wind_save", geo.wcycle);
        wind_save (dummy);
        Log ("Saved wind structure in %s\n", dummy);
      }
      if (modes.make_tables)
      {
        strcpy (dummy, "");
        sprintf (dummy, "diag_%s/%s.%02d", files.root, files.root, geo.wcycle);
        do_windsave2table (dummy, 0);
      }

#ifdef MPI_ON
    }
    MPI_Barrier (MPI_COMM_WORLD);
#endif

    check_time (files.root);
    Log_flush ();               /*Flush the logfile */

  }                             // End of Cycle loop

/* END OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */


  Log (" Completed wind creation.  The elapsed TIME was %f\n", timer ());

  /* SWM - Evaluate wind paths for last iteration */
  if (geo.reverb == REV_WIND || geo.reverb == REV_MATOM)
  {
    wind_paths_evaluate (w, rank_global);
  }

  return (0);
}



/**********************************************************/
/** 
 * @brief      generates the detailed spectra
 *
 * @param [in, out] int  restart_stat   0 if the is run is beginning from
 * scratch, non-zero if this was a restart 
 * @return     Always returns EXIT_SUCCESS
 *
 * @details
 * This is the main routine for calculation detailed
 * spectra in Python.
 *
 * ### Notes ###
 *
 **********************************************************/

int
make_spectra (restart_stat)
     int restart_stat;
{
  WindPtr w;
  PhotPtr p;

  double freqmin, freqmax;
  double renorm;
  long nphot_to_define;
  int iwind;
  int n;

#ifdef MPI_ON
  char dummy[LINELENGTH];
  int spec_spec_helpers;
#endif

  int icheck;

  p = photmain;
  w = wmain;

  freqmax = VLIGHT / (geo.swavemin * 1.e-8);
  freqmin = VLIGHT / (geo.swavemax * 1.e-8);

#ifdef MPI_ON
  /* the length of the big arrays to help with the MPI reductions of the spectra
     the variables for the estimator arrays are set up in the subroutines themselves */
  spec_spec_helpers = (NWAVE * (MSPEC + geo.nangles));  //We need space for NWAVE wavelengths for nspectra, which will eventually equal nangles + MSPEC
#endif

  /* Perform the initilizations required to handle macro-atoms during the detailed
     calculation of the spectrum.  

     Next lines turns off macro atom estimators and other portions of the code that are
     unnecessary during spectrum cycles.  */

  geo.ioniz_or_extract = 0;

/* Next steps to speed up extraction stage */
  if (!modes.keep_photoabs)
  {
    DENSITY_PHOT_MIN = -1.0;    // Do not calculated photoabsorption in detailed spectrum 
  }

  /*Switch on k-packet/macro atom emissivities  SS June 04 */

  if (geo.rt_mode == RT_MODE_MACRO)
  {
    geo.matom_radiation = 1;
  }

  /* Finished initializations required for macro-atom approach */

  /* Calculate and store which bf processess need to be considered in each cell
   * Note that this is not macro-specific but is just to speed the program up.
   */

  kbf_need (freqmin, freqmax);

  /* force recalculation of kpacket rates */
  if (geo.rt_mode == RT_MODE_MACRO)
  {
    for (n = 0; n < NPLASMA; n++)
    {
      macromain[n].kpkt_rates_known = -1;
    }
  }

  /* BEGIN CYCLES TO CREATE THE DETAILED SPECTRUM */

  /* the next section initializes the spectrum array in two cases, for the
   * standard one where one is calulating the spectrum for the first time
   * and in the somewhat abnormal case where additional ionization cycles
   * were calculated for the wind
   */

  if (geo.pcycle == 0)
  {
    spectrum_init (freqmin, freqmax, geo.nangles, geo.angle, geo.phase,
                   geo.scat_select, geo.top_bot_select, geo.select_extract, geo.rho_select, geo.z_select, geo.az_select, geo.r_select);

    /* zero the portion of plasma main that records the numbers of scatters by
     * each ion in a cell
     */

    zero_scatters ();

  }

  /* the next condition should only occur when one has nothing more to do */

  else if (geo.pcycle >= geo.pcycles)
    xsignal (files.root, "%-20s No spectrum   needed: pcycles(%d)==pcycles(%d)\n", "COMMENT", geo.pcycle, geo.pcycles);

  else
  {
    /* Then we are restarting a run with more spectral cycles, but we 
       have already completed some. The memory for the spectral arrays
       should already have been allocated, and the spectrum was initialised
       on the original run, so we just need to renormalise the saved spectrum */
    /* See issue #134 and #503  */

    if (restart_stat == 0)
      Error ("Not restarting, but geo.pcycle = %i and trying to renormalise!\n", geo.pcycle);

    spectrum_restart_renormalise (geo.nangles);
  }


  while (geo.pcycle < geo.pcycles)
  {                             /* This allows you to build up photons in bunches */

    xsignal (files.root, "%-20s Starting %3d of %3d spectrum cycles \n", "NOK", geo.pcycle + 1, geo.pcycles);



    Log ("!!Cycle %d of %d to calculate a detailed spectrum\n", geo.pcycle + 1, geo.pcycles);
    Log_flush ();

    if (!geo.wind_radiation)
      iwind = -1;               /* Do not generate photons from wind */
    else if (geo.pcycle == 0)
      iwind = 1;                /* Create wind photons and force a reinitialization of wind parms */
    else
      iwind = 0;                /* Create wind photons but do not force reinitialization */

    /* Create the initial photon bundles which need to be trannsported through the wind 

       For the detailed spectra, NPHOT*pcycles is the number of photon bundles which will equal the luminosity, 
       1 implies that detailed spectra, as opposed to the ionization of the wind is being calculated

       JM 130306 must convert NPHOT and pcycles to double precision variable nphot_to_define

     */

    NPHOT = NPHOT_MAX;          // Assure that we really are creating as many photons as we expect.

    nphot_to_define = (long) NPHOT *(long) geo.pcycles;
    define_phot (p, freqmin, freqmax, nphot_to_define, 1, iwind, 0);

//OLD    if (modes.save_photons)
//OLD    {
//OLD      for (n = 0; n < NPHOT; n++)
//OLD      {
//OLD        save_photons (&p[n], "CREATE");
//OLD      }
//OLD    }

    for (icheck = 0; icheck < NPHOT; icheck++)
    {
      if (sane_check (p[icheck].freq))
      {
        Error ("python after define phot:sane_check unnatural frequency for photon %d\n", icheck);
      }
    }


    /* Tranport photons through the wind */

    trans_phot (w, p, geo.select_extract);

    spectrum_create (p, geo.nangles, geo.select_extract);

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

      spectrum_summary (files.spec, 0, nspectra - 1, geo.select_spectype, renorm, 0, 0);
      spectrum_summary (files.lspec, 0, nspectra - 1, geo.select_spectype, renorm, 1, 0);

      /* Next lines  produce spectra from photons in the wind only */
      spectrum_summary (files.spec_wind, 0, nspectra - 1, geo.select_spectype, renorm, 0, 1);
      spectrum_summary (files.lspec_wind, 0, nspectra - 1, geo.select_spectype, renorm, 1, 1);

#ifdef MPI_ON
    }
#endif
    Log ("Completed spectrum cycle %3d :  The elapsed TIME was %f\n", geo.pcycle + 1, timer ());

    /* JM1304: moved geo.pcycle++ after xsignal to record cycles correctly. First cycle is cycle 0. */

    xsignal (files.root, "%-20s Finished %3d of %3d spectrum cycles \n", "OK", geo.pcycle + 1, geo.pcycles);

    geo.pcycle++;               // Increment the spectral cycles

#ifdef MPI_ON
    if (rank_global == 0)
    {
#endif
      wind_save (files.windsave);       // This is only needed to update pcycle
      spec_save (files.specsave);
#ifdef MPI_ON
    }
#endif
    check_time (files.root);
  }


/* END CYCLE TO CALCULATE DETAILED SPECTRUM */

#ifdef MPI_ON
  if (rank_global == 0)
  {
#endif
    phot_gen_sum (files.phot, "a");
#ifdef MPI_ON
  }
#endif

  /* SWM0215: Dump the last photon path details to file */
  if (geo.reverb != REV_NONE)
    delay_dump_finish ();       // Each thread dumps to file
#ifdef MPI_ON
  MPI_Barrier (MPI_COMM_WORLD); // Once all done
  if (rank_global == 0 && geo.reverb != REV_NONE)
    delay_dump_combine (np_mpi_global); // Combine results if necessary
#endif


/* Finally done */

#ifdef MPI_ON
  sprintf (dummy, "End of program, Thread %d only", rank_global);       // added so we make clear these are just errors for thread ngit status    
  error_summary (dummy);        // Summarize the errors that were recorded by the program
  Log ("Run py_error.py for full error report.\n");
#else
  error_summary ("End of program");     // Summarize the errors that were recorded by the program
#endif


#ifdef MPI_ON
  MPI_Finalize ();
  Log_parallel ("Thread %d Finalized. All done\n", rank_global);
#endif


  xsignal (files.root, "%-20s %s\n", "COMPLETE", files.root);
  Log ("\nBrief Run Summary\nAt program completion, the elapsed TIME was %f\n", timer ());
  Log ("There were %d of %d ionization cycles and %d of %d spectral cycles run\n", geo.wcycle, geo.wcycles, geo.pcycle, geo.pcycles);
  if (geo.rt_mode == RT_MODE_MACRO)
  {
    if (nlevels_macro == 0)
    {
      Log ("THIS WAS A MACROATOM CALCULATION WITH NO MACROLEVELS. (Use for diagnostics only)\n");
    }
    else
    {
      Log ("This was a macro-atom calculation\n");
    }
  }
  else
  {
    Log ("This was a simple atom calculation\n");
  }

  Log ("Convergence statistics for the wind after the ionization calculation:\n");
  check_convergence ();
  Log ("Information about luminosities and apparent fluxes due to various portions of the system:\n");
  phot_status ();
  return EXIT_SUCCESS;
}
