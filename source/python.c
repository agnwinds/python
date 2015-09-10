/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Python is a program designed to simulate the transfer of radiation in a wind.  It uses the
	Sobolev approximation.  It models a wind as a biconical flow.     
	
	This is the "main" routine for Python.  It's basic purpose is to gather the input variables 
	and to control the flow of the program
 
Arguments:		

	Usage:  py [-h] [-r] [-d] [-f] [-t time_max] xxx  or simply py

	where xxx is the rootname or full name of a parameter file, e. g. test.pf

	and the switches have the following meanings

	-h 	to get this help message
	-r 	restart a run of the progarm reading the file xxx.windsave

	-t time_max	limit the total time to approximately time_max seconds.  Note that the program checks
		for this limit somewhat infrequently, usually at the ends of cycles, because it
		is attempting to save the program outputs so that the program can be restarted with
		-r if theat is desired.
	-v num  determines the amount of information that is printed out.  If num is small, then
		less information is printed out; if num is large more is printed out.  Setting
		v to 5 causes the routine to print out all the information which outputs have
		included previously.  The current default is set to 3 which suppresses Debug, Log_silent
		and Error_silent
	-d	Enters detailed or advanced mode. Allows one to access extra diagnositics and some
	    other advanced commands
    -f  Fixed temperature mode - does not attempt to chenge the temperature of cells.
	-e  Alter the maximum number of errors before the program quits
	-i  Diagnostic mode which quits after reading in inputs. Used for Travis test suite.


	
	if one simply types py or pyZZ where ZZ is the version number one is queried for a name
	of the parameter file.

	NOTE - If this is modified, please also modify the help message in help() below
Returns:
 
Description:	
	Python is far too complicated to describe.  Basically it simulates the radiative transfer
	of photons through the wind of a cataclysmic variable or a star.  The kinematic formulation for
	the CV wind is due to Schlossman and Vitello while that of the star is due to Castor & Larmors. 
	
	 Radiation from an optically thick disk, the WD star, a boundary layer and the wind itself
	 may be included
	
	There are 4 basic portions to the program which are easy to see in the main program.
	
	1. A data gathering stage
	
	2. A calculation of the state of ionization of the wind.
	
	3. A calculation of a detailed spectrum in which the ionization is held fixed.
	
		
Notes:
	The program has been designed and tested both on Suns and MacIntoshes (with Symentec's
	C compiler).  Some of its peculiarities are due to memory limitations and relatively small 
	stack sizes on the Mac.  When compiling, check to see that the global varible MAC is
	set properly in python.h.

History:
 	97jan   ksl	Coding on python began.
 	97jul	ksl	Added capability model the spectrum as a function of orbital phase.
 	97aug	ksl	Added boundary layer as additional source of radiation
 	97aug	ksl	Added capability to extract spectra with a specific number of scatters
 	97nov	ksl	Changed error and diagnostic logging
 	98july	ksl	Added radiation from the wind itself
 	98nov	ksl	Added spherical winds and greatly modified I/O  
	99jan	ksl	Added ability to enter parameter file (w/out) root on command line
	99jul	ksl	Begin eliminating MAC dependencies
	00sep	ksl	Began work on adding a new coronal possibility
	01mar	ksl	Added a knigge wind
        01sept	ksl	Added a thierry wind
	01dec	ksl	Incorporated improved atomic data reading codes, and
			made a variety of other major changes to the way in
			which many of the subroutines operate.
	01dec	ksl	Revised way in which old windfile option works.  It
			was not working at all originally, because of a
			change I made in April.  The new method is a bit dangerous
			because it does change a good bit of geo, since all of the
			input data is read in.  In principle, his can create an inconsistency
			between some of the info in geo and some in the wind array.
			But it keeps the wind the same, which is what was desired.  
	01dec	ksl	Added a general capability to python(40) so that during the
			the ionization cycle photons do not necessarily have identical
			weights.
	02jan	ksl	Moved specific inputs for KWD wind to knigge.c
	02feb	ksl	Added anisotropic wind scattering
	02feb	ksl	41.1 -- Allowed user to chose between anisotropic and
			isotropic scattering in the wind.
	02feb	ksl	41.2 -- Added aditional options for choosing what photons
			to count
	02apr	ksl	43.1 -- Added additional radiative transfer options.  Modified main
			so that a spectrum is printed out after each cycle
	02jul	ksl	43.7 -- Considerable changes associated with photoionization
			and recombination were made.  The initial attemps to allow
			for using a more accurate detailed balance approach were
			included in the ionizatio routines.
	03dec	ksl	Added back some timing ability
        04Mar   SS      Minor changes to set switch for macro atom method.
                        geo.line_mode=6 switches geo.rt_mode=2 (macro method)
                        and then geo.ine_mode back to =3. (SS)
        04Apr   SS      Minor change to set a switch for geo.line_mode=7 
                        which is the same as 6 but with anisotropic scattering
        04May   SS      geo.line_mode=8 added: this mode is the same as mode 6
                        BUT will treat all ions as simple - it is for 
                        test purposes only.
	04Jun	ksl	Added new disk + wind model for ysos.  In process, moved
			the variables describing the wind mdots to subroutines
	04Aug	ksl	52 -- Modified input variables to allow a variety of input models
			all ascii-based.
	04Aug	ksl	52 -- Modified input variables to allow for various disk illumination
			models
	04Aug	ksl	52 -- Modified to allow for vertically extended disks. Note that some of
			the input variables were rearranged in the process.
        04Aug   SS      Minor modifications to cope with finite disk thickness calculations.
	04dec	ksl	54a -- Minor mod to calculation of fraction of luminosity that hits
			disk so that correct values are printed out in situation where there
			are multiple subcycles.
	04dec	ksl	54e -- Minor mod to record master atomic file in geo, and, for the 
			case of fixed concentration to store the file where the concentrations
			are contaned.
	05apr	ksl	56 -- Added new variable geo.binary_system to avoid binary system
			questions in cases where one was modelling a single object.  Note
			that I did not remove seemingly superfluous quesions regarding
			phase because one can still extract photons in various azimuthal
			directions, and if we ever move to 3d situations one might want
			to deal with non-axially symmetric systems.
	05jul	ksl	56d -- Updated to include new windcone definitions.
	05jul	ksl	56d -- Completed basic work to allow vertically extended disk
	06may	ksl	57g -- Began work on modifying structures to lower memory requirments
			At empty wind cells have mostly been eliminated, and macroatoms
			have been split off into a separate structure so this is not 
			used.
	06jul	ksl	57h -- Have contined work to speed program up as much as possible
			by giving more control over the parameters that govern the speed
			such as SMAX_FRAC.
	06oct	ksl	58 -- This version is a cleaned up version of py57ib, which differs
			from 57h primarily in bb.c and pdf.c.   In addition, this version
			makes it possible to use the SV velocity law in what is otherwise
			a KWD description of the wind.  
	06nov	ksl	58c -- This version requires and updated version of kpar, that 
			keeps track of all the errors.  
	07aug	ksl	58f -- Modified the way in which the program handles a situation
			in which one attempts to run with a non-existent .pf file.
	08may	ksl	60a -- Significant modifications were made to the main routine
			to make sure we could restart models.  Unfortunaaely rhe point
			at which the user was asked for the atomic data needed to be changed,
			so earlier .pf files will need to be modified.  
	080808	ksl	62 -- Cleaned up wind ionization options
	081110	ksl	67 -- Revised a number of portions of the main python.c routine
			to enable restarts.  There are still inconsistencies in the
			way one reads information from the .pf file after having read
			in the windsave files that could be cleaned up.
	081218	ksl	67c -- Added a switch -h to provide information about usage.
	090202	ksl	68b -- Added switch -v  to control the level of verbosity.  Updated
			Plasma to allow routine to track scatters and absorption during
			generation of detailed spectrum
	090402	ksl	NSPEC has been moved to the main routine, as its only remaining
			purpose is to define some arrays in the main routine.  Note that
			MSPEC still has meaning as the number of spectra of various types
			that are construced without going through types.
	1108	ksl/nsh Adding code to keep track of gross spectra in a cell though xbands,
			xj and xave_freq.  Just set up the frequence limits here.	
	1112	ksl	Moved everything associated with frequency bands into bands_init
	1212	nsh	changed the way DFUDGE is defined.
        1302	jm	74b5 introduced double precision variable for photon number for calling 
			define_phot. Fixes Bug JM130302 in photon weights. It has been suggested
			that we should make NPHOT a long integer- at the moment I have not done 
			this and simply comverted to double before calling define_phot (otherwise 
			I believe rdint would have to be redone for long integers).
	1304	ksl	75 rewrote the fix above to use a long, instead of double.  
			This has plenty of range.  Notge that there is no need to make NPHOT
			a long, as suggested above.  We do not expect it to exceed 2e9,
			although some kind of error check might be reasonble.
	1306	ksl	Modified to allow a power law component to a stellar spectrum.  Made
			some changes use DEF variables instead of numbers to make choices
	1307	jm	SS Parallelized Python in June 2013, for release 76a. I have now introduced
			slightly altered reporting to allow more succinct reporting in parallel mode.
	1307	ksl	Removed the Thierry & Hubeny O-star models as an option from the code.
			This was never tested, and never really used.  Knox no longer even has the 
			models.  Note that Stuart is replacing this with a homologous expansion
			model
	1308	nsh	Added a call to generate rtheta wind cones - issue #41
	1309	nsh	Changed the loop around where disk parameters are read in - issue #44
	1309	nsh	Added commands to write out warning summary - relating to issue #47
  	1312	nsh	Added a new parameter file command to turn off heating and cooling mechanisms
			at the moment it only does adiabatc
			Also some modifications to the parallel communiactions to deal with some new
			plasma variabales, and the min and max frequency photons seen in bands.
	1409	ksl	Added new switch -d to enable use of a new Debug logging feature
	1411 	JM removed photons per cycle in favour of NPHOT. subcycles are now eliminated
	1501 	JM moved some parallelization stuff to subroutines in para_update.c
			functions are communicate_estimators_para, communicate_matom_estimators_para,
			and gather_spectra_para
	1502		Major reorganisation of input gathering and setup. See setup.c, #136 and #139
	1508	ksl	Instroduction of the concept of domains to handle the disk and wind as separate
			domains
 	
 	Look in Readme.c for more text concerning the early history of the program.

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"


#include "python.h"
#define NSPEC	20

int
main (argc, argv)
     int argc;
     char *argv[];
{
  WindPtr w;
  PhotPtr p;

  double freqmin, freqmax;
  double renorm;
  long nphot_to_define;
  int n;
  int iwind;
  int thermal_opt;		/*NSH 131213 - added to control options to turn on and off some heating and cooling mechanisms */

  /* Next three lines have variables that should be a structure, or possibly we
     should allocate the space for the spectra to avoid all this nonsense.  02feb ksl */

  char dummy[LINELENGTH];
  double x;

  int nn;
  double zz, zzz, zze, ztot, zz_adiab;
  int icheck, nn_adiab;
  FILE *fopen ();

  int disk_illum;
  int opar_stat, restart_stat;
  double time_max;		// The maximum time the program is allowed to run before halting
  double lstar;			// The luminosity of the star, iv it exists

  int my_rank;			// these two variables are used regardless of parallel mode
  int np_mpi;			// rank and number of processes, 0 and 1 in non-parallel
  int time_to_quit;
  int input_int;
  int ndomain = 0;		//Local variable for ndomain
  int ndom;


#ifdef MPI_ON
  int ioniz_spec_helpers, spec_spec_helpers;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &np_mpi);
#else
  my_rank = 0;
  np_mpi = 1;
#endif

  np_mpi_global = np_mpi;	// Global variable which holds the number of MPI processes
  rank_global = my_rank;	// Global variable which holds the rank of the active MPI process
  Log_set_mpi_rank (my_rank, np_mpi);	// communicates my_rank to kpar


  opar_stat = 0;		/* 59a - ksl - 08aug - Initialize opar_stat to indicate that if we do not open a rdpar file, 
				   the assumption is that we are reading from the command line */
  restart_stat = 0;		/* 67 -ksl - 08nov - Assume initially that these is a new run from scratch, and not 
				   a restart
				 */
  time_max = 13.8e9 * 3.2e7;	/* 67 - ksl - 08nov - The maximum time the program will run without stopping.  This
				   is initially set to the lifetime of the universe
				 */
  time_max = -1;
  time_to_quit = 100000;	// Initialise variable




  /* Set the verbosity level for logging.  To get more info raise the verbosity level to a higher number. To
     get less set the verbosity to a lower level. The verbosity can be reset from the comamnd line */
  verbosity = 3;
  Log_set_verbosity (verbosity);

  /* initialise options for advanced mode (all set to 0) */
  init_advanced_modes ();

  /* Parse the command line. Get the root. create files.diagfolder + diagfiles */

  restart_stat = parse_command_line (argc, argv);

  /* If the restart flag has been set, we check to see if a windsave file exists.  If it doues we will 
     we will restart from that point.  If the windsave file does not exist we will start from scratch */

  init_log_and_windsave (restart_stat);

  Log_parallel ("Thread %d starting.\n", my_rank);

  /* Start logging of errors and comments */

  Log ("!!Python Version %s \n", VERSION);	//54f -- ksl -- Now read from version.h
  Log ("!!Git commit hash %s\n", GIT_COMMIT_HASH);
  /* warn the user if there are uncommited changes */
  int git_diff_status = GIT_DIFF_STATUS;
  if (git_diff_status > 0)
    Log
      ("!!Git: This version was compiled with %i files with uncommitted changes.\n",
       git_diff_status);
  Log ("!!Python is running with %d processors\n", np_mpi_global);
  Log_parallel
    ("This is MPI task number %d (a total of %d tasks are running).\n",
     rank_global, np_mpi_global);
  Debug ("Debug statements are on. To turn off use lower verbosity (< 5).\n");

  /* Set the maximum time if it was defined */
  if (time_max > 0)
    {
      set_max_time (files.root, time_max);
    }

  xsignal (files.root, "%-20s Initializing variables for %s\n", "NOK",
	   files.root);

  opar_stat = setup_created_files ();

/* Set plausible values for everything in geo struct which basically defines the overall geometry
As of 1508,  init_geo() also allocates the memory for the domain structure */

  init_geo ();




  /* BEGIN GATHERING INPUT DATA */

  /* Describe the basic calculation in terms of the number of iterations which will
     be used to calculate the wind parameters and the number of iterations and wavelength
     range which will be used for the final spectrom.  Also describe the observer's views
     of the system */


  if (restart_stat == 0)	/* We are starting a new run from scratch, which is the normal
				   mode of operation */
    {

      /* First,  establish the overall system type . 
         Note 1509 - ksl - Exactly what we call a system type is a little bizarre. The original
         intent of this was to allow one to ignore a secondary star, but with addition of AGN it, really
         is a bit unclear what one would like to use here */

      geo.system_type = SYSTEM_TYPE_STAR;
      rdint ("System_type(0=star,1=binary,2=agn,3=previous)",
	     &geo.system_type);


// XXX it is not obious why run_type needs to be in geo.  It is used only in python and setup at present
      geo.run_type = 0;

      if (geo.system_type == SYSTEM_TYPE_PREVIOUS)
	{
	  /* This option is for the confusing case where we want to start with
	     a previous wind model, but we are going to write the result to a
	     new windfile. In other words it is not a restart where we would overwrite
	     the previous wind model.  */

	  strcpy (files.old_windsave, "earlier.run");
	  rdstr ("Old_windfile(root_only)", files.old_windsave);
	  strcat (files.old_windsave, ".wind_save");


	  Log
	    ("Starting a new run from scratch starting with previous windfile");

	  /* Note that wind_read also reads the atomic data file that was used to create the previous run of the data. */
	  if (wind_read (files.old_windsave) < 0)
	    {
	      Error ("python: Unable to open %s\n", files.old_windsave);	//program will exit if unable to read the file
	      exit (0);
	    }
	  geo.run_type = SYSTEM_TYPE_PREVIOUS;	// after wind_read one will have a different wind_type otherwise
	  w = wmain;
	  ndomain = geo.ndomain;	// XXX Needed because currently we set geo.ndomain=ndomain at the end of the inpusts


	}

      else
	{			/* Read the atomic datafile here, because for the cases where we have read
				   and old wind files, we also got the atomic data */

	  rdint
	    ("Wind_type(0=SV,1=Sphere,3=Hydro,4=Corona,5=knigge,6=homologous,7=yso,8=elvis,9=shell,10=None)",
	     &zdom[ndomain].wind_type);
	  if (zdom[ndomain].wind_type != 10)
	    {
	      strcat (zdom[ndomain].name, "Wind");
	      geo.wind_domain_number = ndomain;
	      ndomain++;
	    }
	  else
	    {
	      /* there's no wind, set wind domain_number to -1 */
	      geo.wind_domain_number = -1;
	    }

	  rdstr ("Atomic_data", geo.atomic_filename);

	  /* read a variable which controls whether to save a summary of atomic data
	     this is defined in atomic.h, rather than the modes structure */
	  // XXX Why is this an advanced option; it seems more like a debugging option to me.  
	  if (modes.iadvanced)
	    rdint ("write_atomicdata", &write_atomicdata);

	  if (write_atomicdata)
	    Log ("You have opted to save a summary of the atomic data\n");

	  get_atomic_data (geo.atomic_filename);

	}

      geo.wcycles = geo.pcycles = 1;
      geo.wcycle = geo.pcycle = 0;

    }

  else if (restart_stat == 1)	/* We want to continue a previous run. */
    {
      Log ("Continuing a previous run of %s \n", files.root);
      strcpy (files.old_windsave, files.root);
      strcat (files.old_windsave, ".wind_save");
      if (wind_read (files.old_windsave) < 0)
	{
	  Error ("python: Unable to open %s\n", files.old_windsave);	//program will exit if unable to read the file
	  exit (0);
	}
      w = wmain;
      ndomain = geo.ndomain;	// XXX Needed because currently we set geo.ndomain=ndomain at the end of the inpusts
      geo.run_type = SYSTEM_TYPE_PREVIOUS;	// We read the data from a file
      xsignal (files.root, "%-20s Read %s\n", "COMMENT", files.old_windsave);

      if (geo.pcycle > 0)
	{
	  spec_read (files.specsave);
	  xsignal (files.root, "%-20s Read %s\n", "COMMENT", files.specsave);
	}
    }




/* Get the remainder of the data.  Note that this is done whether or not the windsave file was read in */

  /* 140907 - ksl - Although photons_per_cycle is really an integer, 
     read in as a double so it is easier for input */

  x = 100000;
  rddoub ("photons_per_cycle", &x);
  NPHOT = x;			// NPHOT is photons/cycle

#ifdef MPI_ON
  Log ("Photons per cycle per MPI task will be %d\n", NPHOT / np_mpi_global);

  NPHOT /= np_mpi_global;
#endif

  rdint ("Ionization_cycles", &geo.wcycles);

  rdint ("spectrum_cycles", &geo.pcycles);


  if (geo.wcycles == 0 && geo.pcycles == 0)
    {
      Log
	("Both ionization and spectral cycles are set to 0; There is nothing to do so exiting\n");
      exit (0);			//There is really nothing to do!
    }

  /* Allocate the memory for the photon structure now that NPHOT is established */
  // XXX Not clear why we want to do this here; why not after all of the input data arre in hand

  p = (PhotPtr) calloc (sizeof (p_dummy), NPHOT);

  if (p == NULL)
    {
      Error
	("There is a problem in allocating memory for the photon structure\n");
      exit (0);
    }

  /* Define the coordinate system for the wind grid.  The wind array is allocated later */

  if (geo.wind_domain_number != -1 && geo.run_type != SYSTEM_TYPE_PREVIOUS)
    get_grid_params (ndomain - 1);	// JM PLACEHOLDER -- really we should change the input order here!


  /* ??? ksl - Acoording to line 110 of ionization. option 4 is LTE with SIM_correction.  It would be good to
   * know what this is actually.   Note that pairwise is the appraoch which cboses between pairwise_bb, and pairwise_pow.
   * Normally, any of the pairwise options should force use of a banding option with a broad set of bands
   *
   * NOte that in principle one can change the ioniztion mode between runs
   */

  // XXX  I is unclear to me why all of this dwon to the next XXX is not moved to a single subroutine.  It all
  // pertains to how the radiatiate tranfer is carreid out

  rdint
    ("Wind_ionization(0=on.the.spot,1=LTE,2=fixed,3=recalc_bb,6=pairwise_bb,7=pairwise_pow,8=matrix_bb,9=matrix_pow)",
     &geo.ioniz_mode);

  if (geo.ioniz_mode == IONMODE_FIXED)
    {
      rdstr ("Fixed.concentrations.filename", &geo.fixed_con_file[0]);
    }
  if (geo.ioniz_mode == 4 || geo.ioniz_mode == 5 || geo.ioniz_mode > 9)	/*NSH CLOUDY test - remove once done */
    {
      Log ("The allowed ionization modes are 0, 1, 2, 3, 6, 7\n");
      Error ("Unknown ionization mode %d\n", geo.ioniz_mode);
      exit (0);
    }



  /*Normally, geo.partition_mode is set to -1, which means that partition functions are calculated to take
     full advantage of the data file.  This means that in calculating the partition functions, the information
     on levels and their multiplicities is taken into account.   */

  geo.partition_mode = -1;	//?? Stuart, is there a reason not to move this earlier so it does not affect restart


  /* get_line_transfer_mode reads in the Line_transfer question from the user, 
     then alters the variables geo.line_mode, geo.scatter_mode, geo.rt_mode and geo.macro_simple */

  get_line_transfer_mode ();




  thermal_opt = 0;		/* NSH 131213 Set the option to zero - the default. The lines allow allow the
				   user to turn off mechanisms that affect the thermal balance. Adiabatic is the only one implemented
				   to start off with. */

  rdint
    ("Thermal_balance_options(0=everything.on,1=no.adiabatic)", &thermal_opt);

  if (thermal_opt == 1)
    {
      geo.adiabatic = 0;
    }

  else if (thermal_opt > 1 || thermal_opt < 0)
    {
      Error ("Unknown thermal balance mode %d\n", thermal_opt);
      exit (0);
    }


  /* 57h -- Next line prevents bf calculation of macro_estimaters when no macro atoms are present.   */

  if (nlevels_macro == 0)
    geo.macro_simple = 1;	// Make everything simple if no macro atoms -- 57h

  //SS - initalise the choice of handling for macro pops.
  if (geo.run_type == SYSTEM_TYPE_PREVIOUS)
    {
      geo.macro_ioniz_mode = 1;	// Now that macro atom properties are available for restarts
    }
  else
    {
      geo.macro_ioniz_mode = 0;
    }

  /* specify if there is a disk and what type */
  /* JM 1502 -- moved disk type question here- previously it was just before
     asking for disk radiation. See #8 and #44 */

  // XXX End here

  rdint
    ("disk.type(0=no.disk,1=standard.flat.disk,2=vertically.extended.disk)",
     &geo.disk_type);

  if (geo.disk_type == 0)
    {
      geo.disk_atmosphere = 0;
    }
  else
    {
      /* ksl 1508 Add parameters for a disk atmosphere XXX  */
      // XXX This looks like a problem for restarts
      zdom[ndomain].ndim = 30;
      zdom[ndomain].mdim = 10;

      rdint ("disk.atmosphere(0=no,1=yes)", &geo.disk_atmosphere);
    }
  if (geo.disk_atmosphere != 0)
    {
      /* specify the domain name and number */
      strcat (zdom[ndomain].name, "Disk Atmosphere");
      geo.atmos_domain_number = ndomain;
//XXX Whyi isn't get_grid_param used here.  We could pass a variable to it, to mofifiy the rdpar names if
// necessary
      input_int = 1;
      rdint ("atmos.coord.system(1=cylindrical,2=spherical_polar,3=cyl_var)",
	     &input_int);

      switch (input_int)
	{
	case 0:
	  zdom[ndomain].coord_type = SPHERICAL;
	  break;
	case 1:
	  zdom[ndomain].coord_type = CYLIND;
	  break;
	case 2:
	  zdom[ndomain].coord_type = RTHETA;
	  break;
	case 3:
	  zdom[ndomain].coord_type = CYLVAR;
	  break;
	default:
	  Error
	    ("Invalid parameter supplied for 'Coord_system'. Valid coordinate types are: \n\
                        0 = Spherical, 1 = Cylindrical, 2 = Spherical polar, 3 = Cylindrical (varying Z)");
	}


      rdint ("atmos.dim.in.x_or_r.direction", &zdom[ndomain].ndim);
      rdint ("atmos.dim.in.z_or_theta.direction", &zdom[ndomain].mdim);
      ndomain++;
    }




  /* Determine what radiation sources there are.  
     Note that most of these values are initilized in init_geo */

  get_radiation_sources ();


  if (geo.run_type == SYSTEM_TYPE_PREVIOUS)
    {
      disk_illum = geo.disk_illum;
    }


  if (geo.run_type != SYSTEM_TYPE_PREVIOUS)	// Start of block to define a model for the first time
    {

      /* get_stellar_params gets information like mstar, rstar, tstar etc.
         it returns the luminosity of the star */
      lstar = get_stellar_params ();


      /* Describe the disk */

      if (geo.disk_type)	/* Then a disk exists and it needs to be described */
	{
	  disk_illum = get_disk_params ();
	}

      else
	{
	  /* There is no disk so set variables accordingly */
	  geo.disk_radiation = 0;
	  geo.diskrad = 0;
	}

      /* describe the boundary layer / agn components to the spectrum if they exist. 
         reads in information specified by the user and sets variables in geo structure */
      get_bl_and_agn_params (lstar);



      /* Describe the wind. This routine reads in geo.rmax and geo.twind
         and then gets params by calling e.g. get_sv_wind_params() */
      /* PLACEHOLDER -- XXX call with wind domain number */
      get_wind_params (geo.wind_domain_number);

    }				// End of block to define a model for the first time

  else
    {
      if (geo.disk_type)	/* Then a disk exists and it needs to be described */
	{
	  if (geo.disk_radiation)
	    {
	      rdint
		("Disk.temperature.profile(0=standard;1=readin)",
		 &geo.disk_tprofile);
	      if (geo.disk_tprofile == 1)
		{
		  rdstr ("T_profile_file", files.tprofile);
		}
	    }
	}
    }


  /* Calculate additional parameters associated with the binary star system */

  // XXX This may be a problem for restarts, now that we have changed the meaing of 
  // SYSTEM type
  if (geo.system_type == SYSTEM_TYPE_BINARY)
    binary_basics ();

  /* Check that the parameters which have been supplied for the star, disk and boundary layer will
     allow generation of photons where that is appropriate */

  if (geo.tstar <= 0.0)
    geo.star_radiation = 0;
  if (geo.disk_mdot <= 0.0)
    geo.disk_radiation = 0;
  if (geo.t_bl <= 0.0 || geo.lum_bl <= 0.0)
    geo.bl_radiation = 0;

  /* Next block added by SS August 04 in case diskrad = 0 is used to mean no disk at any point. */

  if (geo.diskrad <= 0.0)
    {
      geo.disk_type = 0;
      geo.disk_radiation = 0;
    }

  if (geo.star_radiation)
    Log ("There is a star which radiates\n");
  else
    Log ("The star in the system does not radiate\n");

  if (!geo.disk_type)
    Log ("There is no disk in the system \n");
  else if (!geo.disk_radiation)
    Log ("The disk exists, but only absorbs photons\n");
  else
    Log ("There is a disk which radiates and absorbs\n");

  if (geo.bl_radiation)
    Log ("There is a boundary layer which radiates\n");
  else
    Log ("There is no boundary layer\n");

  if (geo.agn_radiation)
    Log ("There is a BH  which radiates\n");
  else
    Log ("There is no BH \n");

  /* Describe the spectra which will be extracted and the way it will be extracted */

  /* First initialise things to semi-reasonable values */
/* These two variables have to do with what types of spectra are created n the
 * spectrum files. They are not associated with the nature of the spectra that
 * are generated by say the boundary layer
 */

  geo.select_extract = 1;
  geo.select_spectype = 1;

/* Completed initialization of this section.  Note that get_spectype uses the source of the
 * ratiation and then value given to return a spectrum type. The output is not the same 
 * number as one inputs. It's not obvious that this is a good idea. */

  if (geo.pcycles > 0)
    {

      get_spectype (geo.star_radiation,
		    "Rad_type_for_star(0=bb,1=models,2=uniform)_in_final_spectrum",
		    &geo.star_spectype);


      get_spectype (geo.disk_radiation,
		    "Rad_type_for_disk(0=bb,1=models,2=uniform)_in_final_spectrum",
		    &geo.disk_spectype);


      get_spectype (geo.bl_radiation,
		    "Rad_type_for_bl(0=bb,1=models,2=uniform)_in_final_spectrum",
		    &geo.bl_spectype);

      geo.agn_spectype = 3;
      get_spectype (geo.agn_radiation,
		    "Rad_type_for_agn(3=power_law,4=cloudy_table)_in_final_spectrum",
		    &geo.agn_spectype);



      init_observers ();
    }

  geo.matom_radiation = 0;	//initialise for ionization cycles - don't use pre-computed emissivities for macro-atom levels/ k-packets.


  /* 57h -- New section of inputs to provide more control over how the program is
     run -- 07jul -- ksl
     1502 JM -- moved to subroutine
   */

  get_standard_care_factors ();

  /* 0415 SWM - Added metaparams */
  get_meta_params ();


/* 081221 - 67c - Establish limits on the frequency intervals to be used by the ionization cycles and 
 * the fraquency bands for stratified sampling.  Changes here were made to allow more control
 * over statified sampling, since we have expanded the temperature ranges of the types of systems
 * we would like to calculate.  This section of inputs might logically go earlier in the code, but
 * was put here so it would add on to existing .pf files.  It would be reasonble to consider moving
 * it to a more logical location
 */


/* Determine the frequency range which will be used to establish the ionization balance of the wind */

  // Note that bands_init asks .pf file or user what kind of banding is desired 

  bands_init (-1, &xband);

  /*if we have changed min and max in bands_init, we need to make sure this is reflected in the frequency bounds */
  freqmin = xband.f1[0];
  freqmax = xband.f2[xband.nbands - 1];

  /* 1112 - 71 - ksl Next routine sets up the frequencies that are used for charactizing the spectrum in a cell
   * These need to be coordinated with the bands that are set up for spectral gneration
   */
  freqs_init (freqmin, freqmax);


  if (modes.iadvanced)
    {
      /* Do we require extra diagnostics or not */
      rdint ("Extra.diagnostics(0=no,1=yes) ", &modes.diag_on_off);

      if (modes.diag_on_off)
	{
	  get_extra_diagnostics ();
	}
    }



  /* Wrap up and save all the inputs */

  if (strncmp (files.root, "mod", 3) == 0)
    cpar ("mod.pf");
  else if (strncmp (files.root, "dummy", 5) == 0)
    {
      cpar ("dummy.pf");
      exit (0);
    }
  else if (opar_stat == 1)
    {
      cpar (files.input);
    }
  else
    cpar ("python.pf");


  /* OK all inputs have been obtained at this point and the inputs have been copied to "mod.pf" or "python.pf" */
  /* JM 1502 -- if we have used the -i flag we want to quit after inputs as we were just testing readin */
  if (modes.quit_after_inputs)
    {
      Log ("Run with -i flag, so quitting now inputs have been gathered.\n");
      exit (0);
    }

/* INPUTS ARE FINALLY COMPLETE */

/* Print out some diagnositic infomration about the domains */

  // XXX This is clearly wroing for repeats
  geo.ndomain = ndomain;	// Store ndomain in geo so that it can be saved

  Log ("There are %d domains\n", geo.ndomain);
  for (n = 0; n < geo.ndomain; n++)
    {
      Log ("%20s %d %d %d %d\n", zdom[n].name, zdom[n].wind_type,
	   zdom[n].ndim, zdom[n].mdim, zdom[n].ndim2);
    }


  /* 121219 NSH Set up DFUDGE to be a value that makes some kind of sense
     given the scale of the wind. Up till py74b2 it was set to be fixed at
     1e5, so we ensure that this is a minimum, so any winds of CV type scale
     will keep the old dfudge, and hopefully look the same. We also need to
     set defudge slightly differently for the shell wind. */

  DFUDGE = setup_dfudge ();


  /* Now define the wind cones generically. modifies the global windcone structure */

  setup_windcone ();


  /*NSH 130821 broken out into a seperate routine added these lines to fix bug41, where
     the cones are never defined for an rtheta grid if the model is restarted. 

     XXX ksl is unclear why the wind cones ar being initilized here, rather than as part of
     routines located elsewhere, but I have follwoed previous practice and reinitialized them
     as part to the domain effort.  
   */

  for (ndom = 0; ndom < geo.ndomain; ndom++)
    {

      if (zdom[0].coord_type == RTHETA && geo.run_type == SYSTEM_TYPE_PREVIOUS)	//We need to generate an rtheta wind cone if we are restarting
	{
	  rtheta_make_cones (ndom, wmain);
	}
    }





  /* Next line finally defines the wind if this is the initial time this model is being run */
  if (geo.run_type != SYSTEM_TYPE_PREVIOUS)	// Define the wind and allocate the arrays the first time
    {
      define_wind ();
    }
  // Do not reinit if you want to use old windfile

  w = wmain;

  if (modes.save_cell_stats)
    {
      /* Open a diagnostic file or files.  These are all fixed files */
      open_diagfile ();
    }

  /* initialize the random number generator */
  //      srand( (n=(unsigned int) clock()));  
  srand (1084515760 + (13 * rank_global));

  /* 68b - 0902 - ksl - Start with photon history off */

  phot_hist_on = 0;

  /* If required, read in a non-standard disk temperature profile */

  if (geo.disk_tprofile == 1)
    {
      read_non_standard_disk_profile (files.tprofile);
    }



/* The next section sets up a structure qdisk to record the effects
 * of illumination on the disk.  disk_init is called primarily to get
 * a defined set of annular rings which are kept throughout the
 * ionization calculation.  A second structure qdisk is needed
 * because in the process of generating photons in various bands
 * the annular rings are changed
 *
 * disk_init calculates the flux from the disk in the energy range set by
 * freqmin and freqmax, and uses is this to identify the position of the
 * rings in the disk, so that each ring contributes the same amount to
 * the flux
 *
 * */


  disk_init (geo.rstar, geo.diskrad, geo.mstar, geo.disk_mdot, freqmin,
	     freqmax, 0, &geo.f_disk);

  qdisk_init ();		/* Initialize a disk qdisk to store the information about photons impinging on the disk */

/* 04aug -- ksl -- now that everything is initialized, we set geo.disk_illum
 *
 * 080518 - ksl - I believe that the reason for this somewhat weird logic is to
 * assure that models (e.g corona and knigge) where the base wind velocity
 * depends on teff are not altered by illumination, but since no photons
 * have been transported at this stage, it's hard to see why that would
 * matter.
 */

  geo.disk_illum = disk_illum;

  xsignal (files.root, "%-20s Finished initialization for %s\n", "NOK",
	   files.root);
  check_time (files.root);

#ifdef MPI_ON
  /* Since the wind is now set up can work out the length big arrays to help with the MPI reductions of the spectra
     the variables for the estimator arrays are set up in the subroutines themselves */
  ioniz_spec_helpers = 2 * MSPEC * NWAVE;	//we need space for log and lin spectra for MSPEC XNWAVE
  spec_spec_helpers = (NWAVE * (MSPEC + nangles));	//We need space for NWAVE wavelengths for nspectra, which will eventually equal nangles + MSPEC

#endif



/* XXXX - THE CALCULATION OF THE IONIZATION OF THE WIND */

  geo.ioniz_or_extract = 1;	//SS July 04 - want to compute MC estimators during ionization cycles
  //1 simply implies we are in the ionization section of the code
  //and allows routines to act accordinaly.

/* 67 -ksl- geo.wycle will start at zero unless we are completing an old run */

/* XXXX - BEGINNING OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */

  if (geo.wcycle == geo.wcycles)
    xsignal (files.root,
	     "%-20s No ionization needed: wcycles(%d)==wcyeles(%d)\n",
	     "COMMENT", geo.wcycle, geo.wcycles);
  else
    {
      geo.pcycle = 0;		/* Set the spectrum cycles executed to 0, because 
				   we are going to modify the wind and hence any
				   previously calculated spectra must be recreated
				 */
    }

  /* SWM - Setup for path tracking */
  if (geo.reverb > REV_NONE)
    {
      reverb_init (wmain, geo.nangles, freqmin, freqmax);
      delay_dump_prep (files.root, restart_stat, rank_global);
    }


  while (geo.wcycle < geo.wcycles)
    {				/* This allows you to build up photons in bunches */

      xsignal (files.root, "%-20s Starting %d of %d ionization cycle \n",
	       "NOK", geo.wcycle, geo.wcycles);


      Log ("!!Python: Begining cycle %d of %d for defining wind\n",
	   geo.wcycle, geo.wcycles);
      Log_flush ();		/*NH June 13 Added call to flush logfile */

      /* Initialize all of the arrays, etc, that need initialization for each cycle
       */

      spectrum_init (freqmin, freqmax, geo.nangles, geo.angle, geo.phase,
		     geo.scat_select, geo.top_bot_select,
		     geo.select_extract, geo.rho_select, geo.z_select,
		     geo.az_select, geo.r_select);


      wind_rad_init ();		/*Zero the parameters pertaining to the radiation field */



      if (modes.ispy)
	ispy_init ("python", geo.wcycle);


      geo.n_ioniz = 0.0;
      geo.lum_ioniz = 0.0;
      ztot = 0.0;		/* ztot is the luminosity of the disk multipled by the number of cycles, which is used by save_disk_heating */

      /* JM 1409 -- We used to execute subcycles here, but these have been removed */

      if (!geo.wind_radiation
	  || (geo.wcycle == 0 && geo.run_type != SYSTEM_TYPE_PREVIOUS))
	iwind = -1;		/* Do not generate photons from wind */
      else
	iwind = 1;		/* Create wind photons and force a reinitialization of wind parms */

      /* Create the photons that need to be transported through the wind
       *
       * NPHOT is the number of photon bundles which will equal the luminosity; 
       * 0 => for ionization calculation 
       */


      /* JM 130306 need to convert photons_per_cycle to double precision for define_phot */
      /* ksl 130410 - This is needed here not because we expect photons per cycle to 
       * exceed the size of an integer, but because of the call to define phot in the
       * spectrum cycle, which can exceed this
       */
      /* JM 1409 photons_per_cycle has been removed in favour of NPHOT */

      nphot_to_define = (long) NPHOT;

      define_phot (p, freqmin, freqmax, nphot_to_define, 0, iwind, 1);

      /* Zero the arrays that store the heating of the disk */

      /* 080520 - ksl - There is a conundrum here.  One should really zero out the 
       * quantities below each time the wind structure is updated.  But relatively
       * few photons hit the disk under normal situations, and therefore the statistcs
       * are not very good.  
       */

      /* 130213 JM -- previously this was done before define_phot, which meant that
         the ionization state was never computed with the heated disk */

      for (n = 0; n < NRINGS; n++)
	{
	  qdisk.heat[n] = qdisk.nphot[n] = qdisk.w[n] = qdisk.ave_freq[n] = 0;
	}



      photon_checks (p, freqmin, freqmax, "Check before transport");

      wind_ip ();


      zz = 0.0;
      for (nn = 0; nn < NPHOT; nn++)
	{
	  zz += p[nn].w;
	}

      Log
	("!!python: Total photon luminosity before transphot %18.12e\n", zz);
      Log_flush ();		/* NSH June 13 Added call to flush logfile */
      ztot += zz;		/* Total luminosity in all cycles, used for calculating disk heating */

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

      /*Determine how much energy was absorbed in the wind */
      zze = zzz = zz_adiab = 0.0;
      nn_adiab = 0;
      for (nn = 0; nn < NPHOT; nn++)
	{
	  zzz += p[nn].w;
	  if (p[nn].istat == P_ESCAPE)
	    zze += p[nn].w;
	  if (p[nn].istat == P_ADIABATIC)
	    {
	      zz_adiab += p[nn].w;
	      nn_adiab++;
	    }
	}

      Log
	("!!python: Total photon luminosity after transphot %18.12e (diff %18.12e). Radiated luminosity %18.12e\n",
	 zzz, zzz - zz, zze);
      if (geo.rt_mode == 2)
	Log
	  ("Luminosity taken up by adiabatic kpkt destruction %18.12e number of packets %d\n",
	   zz_adiab, nn_adiab);

      if (modes.print_windrad_summary)
	wind_rad_summary (w, files.windrad, "a");




      photon_checks (p, freqmin, freqmax, "Check after transport");

      spectrum_create (p, freqmin, freqmax, geo.nangles, geo.select_extract);



      /* At this point we should communicate all the useful infomation 
         that has been accummulated on differenet MPI tasks */

#ifdef MPI_ON

      communicate_estimators_para ();

      communicate_matom_estimators_para ();	// this will return 0 if nlevels_macro == 0
#endif




      if (modes.ispy)
	ispy_close ();


      /* Calculate and store the amount of heating of the disk due to radiation impinging on the disk */
      qdisk_save (files.disk, ztot);

/* Completed writing file describing disk heating */

      Log
	("!!python: Number of ionizing photons %g lum of ionizing photons %g\n",
	 geo.n_ioniz, geo.lum_ioniz);

/* This step shoudl be MPI_parallelised too */

      wind_update (w);

/* In a diagnostic mode save the wind file for each cycle (from thread 0) */

      if (modes.keep_ioncycle_windsaves)
	{
	  strcpy (dummy, "");
	  sprintf (dummy, "python%02d.wind_save", geo.wcycle);

#ifdef MPI_ON
	  if (rank_global == 0)
	    {
#endif
	      wind_save (dummy);
#ifdef MPI_ON
	    }
#endif
	  Log ("Saved wind structure in %s\n", dummy);
	}


      Log ("Completed ionization cycle %d :  The elapsed TIME was %f\n",
	   geo.wcycle, timer ());

      Log_silent ("Finished creating spectra\n");

      /* Do an MPI reduce to get the spectra all gathered to the master thread */

#ifdef MPI_ON

      gather_spectra_para (ioniz_spec_helpers, MSPEC);

#endif



#ifdef MPI_ON
      if (rank_global == 0)
	{
#endif

	  spectrum_summary (files.wspec, "w", 0, 6, 0, 1., 0);
	  spectrum_summary (files.lspec, "w", 0, 6, 0, 1., 1);	/* output the log spectrum */

#ifdef MPI_ON
	}
      MPI_Barrier (MPI_COMM_WORLD);
#endif
      phot_gen_sum (files.phot, "w");	/* Save info about the way photons are created and absorbed
					   by the disk */

      /* Save everything after each cycle and prepare for the next cycle 
         JM1304: moved geo.wcycle++ after xsignal to record cycles correctly. First cycle is cycle 0. */
      /* NSH1306 - moved geo.wcycle++ back, but moved the log and xsignal statements */


      xsignal (files.root, "%-20s Finished %d of %d ionization cycle \n",
	       "OK", geo.wcycle, geo.wcycles);
      geo.wcycle++;		//Increment ionisation cycles


/* NSH 1408 - Save only the windsave file from thread 0, to prevent many processors from writing to the same
 * file. */

#ifdef MPI_ON
      if (rank_global == 0)
	{
#endif
	  wind_save (files.windsave);
	  Log_silent ("Saved wind structure in %s after cycle %d\n",
		      files.windsave, geo.wcycle);
#ifdef MPI_ON
	}
      MPI_Barrier (MPI_COMM_WORLD);
#endif




      check_time (files.root);
      Log_flush ();		/*Flush the logfile */

    }				// End of Cycle loop

/* XXXX - END OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */


  Log (" Completed wind creation.  The elapsed TIME was %f\n", timer ());

  /* SWM - Evaluate wind paths for last iteration */
  if (geo.reverb == REV_WIND)
    {
      wind_paths_evaluate (w);
      wind_paths_output (w, files.root);
    }

/* XXXX - THE CALCULATION OF A DETAILED SPECTRUM IN A SPECIFIC REGION OF WAVELENGTH SPACE */

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
