/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Python is a program designed to simulate the transfer of radiation in a wind.  It uses the
	Sobolev approximation.  It models a wind as a biconical flow.     
	
	This is the "main" routine for Python.  It's basic purpose is to gather the input variables 
	and to control the flow of the program
 
Arguments:		

	Usage:  py [-h] [-r] [-d] [-t time_max] xxx  or simply py

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
	-e  Alter the maximum number of errors before the program quits


	
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
 	
 	Look in Readme.c for more text concerning the early history of the program.

**************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atomic.h"


#include "python.h"
#define NSPEC	20
#define SYSTEM_TYPE_STAR    0
#define SYSTEM_TYPE_BINARY  1
#define SYSTEM_TYPE_AGN     2

int main(argc, argv)
	 int argc;
	 char *argv[];
{
	WindPtr w;
	PhotPtr p;

	int i;
	double freqmin, freqmax;
	double swavemin, swavemax, renorm;
	long nphot_to_define;
	int n, nangles;
	int iwind;
	int thermal_opt;			/* NSH 131213 - added to control options to turn on and off some heating and cooling mechanisms */

	/* Next three lines have variables that should be a structure, or possibly we should allocate the space for the spectra to
	   avoid all this nonsense.  02feb ksl */

	double angle[NSPEC], phase[NSPEC];
	int scat_select[NSPEC], top_bot_select[NSPEC];
	double rho_select[NSPEC], z_select[NSPEC], az_select[NSPEC], r_select[NSPEC];

	char yesno[20];
	int select_extract, select_spectype;
	char root[LINELENGTH], input[LINELENGTH], wspecfile[LINELENGTH],
		lspecfile[LINELENGTH], specfile[LINELENGTH], diskfile[LINELENGTH];
	char windradfile[LINELENGTH], windsavefile[LINELENGTH];
	char specsavefile[LINELENGTH];
	char photfile[LINELENGTH], diagfile[LINELENGTH], old_windsavefile[LINELENGTH], diagfolder[LINELENGTH];
	char dummy[LINELENGTH];
	char tprofile[LINELENGTH];
	char delay_dump_file[LINELENGTH];	// SWM
	double x, xbl;

	int nn;
	double zz, zzz, zze, ztot, zz_adiab;
	int icheck, nn_adiab;
	FILE *fopen(), *qptr;

	int disk_illum;
	int istandard, keep_photoabs;
	int opar_stat, restart_stat;
	double time_max;			// The maximum time the program is allowed to run before halting
	double lstar;				// The luminosity of the star, iv it exists

	int my_rank;				// these two variables are used regardless of parallel mode
	int np_mpi;					// rank and number of processes, 0 and 1 in non-parallel
	int time_to_quit;

	int mkdir();

#ifdef MPI_ON
	int ioniz_spec_helpers, spec_spec_helpers;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np_mpi);
#else
	my_rank = 0;
	np_mpi = 1;
#endif

	np_mpi_global = np_mpi;		// / Glob al variable which holds the number of MPI processes
	rank_global = my_rank;		// / Global variable which holds the rank of the active MPI process

	Log_set_mpi_rank(my_rank, np_mpi);	// communicates my_rank to kpar


	opar_stat = 0;				/* 59a - ksl - 08aug - Initialize opar_stat to indicate that if we do not open a rdpar file, the
								   assumption is that we are reading from the command line */
	restart_stat = 0;			/* 67 -ksl - 08nov - Assume initially that these is a new run from scratch, and not a restart */
	time_max = 13.8e9 * 3.2e7;	/* 67 - ksl - 08nov - The maximum time the program will run without stopping.  This is initially
								   set to the lifetime of the universe */

	verbosity = 3;				/* Set the default verbosity to 3.  To get more info raise the verbosity level to a higher number.
								   To get less set the verbosity to a lower level. */

	time_to_quit = 100000;		// Initialise variable


	Log_set_verbosity(verbosity);


	Log_parallel("Thread %d starting.\n", my_rank);	// JM130723 moved this after verbosity switch



	restart_stat = 0;
	time_max = -1;

	/* initialise options for advanced mode (all set to 0) */
	init_advanced_modes();


	/* Parse the command line.  */

	if (argc == 1)
	{
		printf("Input file (interactive=stdin):");
		fgets(dummy, LINELENGTH, stdin);
		get_root(root, dummy);
		strcpy(diagfile, root);
		strcat(diagfile, ".diag");
	}
	else
	{

		for (i = 1; i < argc; i++)
		{

			if (strcmp(argv[i], "-h") == 0)
			{
				help();
			}
			else if (strcmp(argv[i], "-r") == 0)
			{
				Log("Restarting %s\n", root);
				restart_stat = 1;
			}
			else if (strcmp(argv[i], "-t") == 0)
			{
				if (sscanf(argv[i + 1], "%lf", &time_max) != 1)
				{
					Error("python: Expected time after -t switch\n");
					exit(0);
				}
				i++;

			}
			else if (strcmp(argv[i], "-v") == 0)
			{
				if (sscanf(argv[i + 1], "%d", &verbosity) != 1)
				{
					Error("python: Expected verbosity after -v switch\n");
					exit(0);
				}
				Log_set_verbosity(verbosity);
				i++;

			}
			else if (strcmp(argv[i], "-e") == 0)
			{
				if (sscanf(argv[i + 1], "%d", &time_to_quit) != 1)
				{
					Error("python: Expected max errors after -e switch\n");
					exit(0);
				}
				Log_quit_after_n_errors(time_to_quit);
				i++;

			}
			else if (strcmp(argv[i], "-d") == 0)
			{
				modes.iadvanced = 1;
				i++;
			}
			else if (strncmp(argv[i], "-", 1) == 0)
			{
				Error("python: Unknown switch %s\n", argv[i]);
				help();
			}
		}



		/* The last command line variable is always the .pf file */

		strcpy(dummy, argv[argc - 1]);
		get_root(root, dummy);

		/* This completes the parsing of the command line */

		/* JM130722 we now store diag files in a subdirectory if in parallel */
		/* ksl - I believe this is created in all cases, and that is what we want */

		sprintf(diagfolder, "diag_%s/", root);
		mkdir(diagfolder, 0777);
		strcpy(diagfile, diagfolder);
		sprintf(dummy, "_%d.diag", my_rank);
		strcat(diagfile, root);
		strcat(diagfile, dummy);


	}


	/* 0811 - ksl - If the restart flag has been set, we check to see if a windsave file exists.  If it doues we will we will
	   restart from that point.  If the windsave file does not exist we will start from scratch */

	if (restart_stat == 0)
	{							// Then we are simply running from a new model
		xsignal_rm(root);		// Any old signal file
		xsignal(root, "%-20s %s \n", "START", root);
		Log_init(diagfile);
	}
	else
	{
		/* Note that alghough we chekc that we dan open the windsave file, it is not read here.  */

		strcpy(windsavefile, root);
		strcat(windsavefile, ".wind_save");
		qptr = fopen(windsavefile, "r");

		if (qptr != NULL)
		{
			/* Then the file does exist and we can restart */
			fclose(qptr);
			xsignal(root, "%-20s %s\n", "RESTART", root);
			Log_append(diagfile);
		}
		else
		{
			/* It does not exist and so we start from scratch */
			restart_stat = 0;
			xsignal_rm(root);	// Any old signal file
			xsignal(root, "%-20s %s \n", "START", root);
			Log_init(diagfile);
		}
	}




	/* Start logging of errors and comments */

	Log("!!Python Version %s \n", VERSION);	// 54f -- ksl -- Now read from version.h
	Log("!!Git commit hash %s\n", GIT_COMMIT_HASH);
	Log("!!Python is running with %d processors\n", np_mpi_global);
	Log_parallel("This is MPI task number %d (a total of %d tasks are running).\n", rank_global, np_mpi_global);
	Debug("Debug statements are on. To turn off use lower verbosity (< 5).\n");
	/* Set the maximum time if it was defined */
	if (time_max > 0)
	{
		set_max_time(root, time_max);
	}


	xsignal(root, "%-20s Initializing variables for %s\n", "NOK", root);


	if (strncmp(root, "dummy", 5) == 0)
	{
		Log("Proceeding to create rdpar file in dummy.pf, but will not run prog\n");
	}
	else if (strncmp(root, "stdin", 5) == 0 || strncmp(root, "rdpar", 5) == 0 || root[0] == ' ' || strlen(root) == 0)
	{
		strcpy(root, "mod");
		Log("Proceeding in interactive mode\n Output files will have rootname mod\n");
	}
	else
	{
		strcpy(input, root);
		strcat(input, ".pf");

		if ((opar_stat = opar(input)) == 2)
		{
			Log("Reading data from file %s\n", input);
		}
		else
		{
			Log("Creating a new parameter file %s\n", input);
		}

	}

	/* Now create the names of all the files which will be written.  Note that some files have the same root as the input file,
	   while others have a generic name of python. This is intended so that files which you really want to keep have unique names,
	   while those which are for short-term diagnostics are overwritten.  ksl 97aug. */


	strcpy(basename, root);		// 56d -- ksl --Added so filenames could be created by routines as necessary

	strcpy(wspecfile, root);
	strcpy(lspecfile, root);

	strcpy(specfile, root);
	strcpy(windradfile, "python");
	strcpy(windsavefile, root);
	strcpy(specsavefile, root);
	strcpy(delay_dump_file, root); //SWM

	/* 130722 JM we now save python.phot and disk.diag files under diag_root folder */
	strcpy(photfile, diagfolder);
	strcpy(diskfile, diagfolder);
	strcat(photfile, "python");
	strcat(diskfile, root);

	strcat(wspecfile, ".spec_tot");
	strcat(lspecfile, ".log_spec_tot");
	strcat(specfile, ".spec");
	strcat(windradfile, ".wind_rad");
	strcat(windsavefile, ".wind_save");
	strcat(specsavefile, ".spec_save");
	strcat(photfile, ".phot");
	strcat(diskfile, ".disk.diag");
	strcat(delay_dump_file, ".delay_dump"); //SWM



	/* Provide plausible initial values for the sizes of the wind arrays.  This is desirable primarily for creating reasonable .pf
	   files */

	/* Set plausible values for everything in geo struct which basically defines the overall geometry */

	init_geo();

	/* Set the global variables that define the size of the grid as defined in geo.  These are used for convenience */

	NDIM = geo.ndim;
	MDIM = geo.mdim;
	NDIM2 = geo.ndim * geo.mdim;

	/* End of definition of wind arrays */


	/* Initialize variables which are used in the main routine */


	/* Initialize basis vectors for a cartesian coordinate system */

	x_axis[0] = 1.0;
	x_axis[1] = x_axis[2] = 0.0;
	y_axis[1] = 1.0;
	y_axis[0] = y_axis[2] = 0.0;
	z_axis[2] = 1.0;
	z_axis[1] = z_axis[0] = 0.0;





	/* BEGIN GATHERING INPUT DATA */

	/* Describe the basic calculation in terms of the number of iterations which will be used to calculate the wind parameters and
	   the number of iterations and wavelength range which will be used for the final spectrom.  Also describe the observer's views
	   of the system */


	if (restart_stat == 0)		/* We are starting a new run from scratch */
	{
		/* Note that these describe wind geometryies and not the type of object */


		rdint("Wind_type(0=SV,1=Sphere,2=Previous,3=Proga,4=Corona,5=knigge,6=homologous,7=yso,8=elvis,9=shell)", &geo.wind_type);


		if (geo.wind_type == 2)
		{
			/* This option is for the confusing case where we want to start with a previous wind model, but we are going to write
			   the result to a new windfile. In other words it is not a restart where we would overwrite the previous wind model.  */

			strcpy(old_windsavefile, "earlier.run");
			rdstr("Old_windfile(root_only)", old_windsavefile);
			strcat(old_windsavefile, ".wind_save");


			Log("Starting a new run from scratch starting with previous windfile");
			if (wind_read(old_windsavefile) < 0)
			{
				Error("python: Unable to open %s\n", old_windsavefile);	// program will exit if unable to read the file
				exit(0);
			}
			geo.wind_type = 2;	// after wind_read one will have a different wind_type otherwise
			w = wmain;


		}

		else
		{						/* Read the atomic datafile here, because for the cases where we have read and old wind files, we
								   also got the atomic data */

			rdstr("Atomic_data", geo.atomic_filename);

			/* read a variable which controls whether to save a summary of atomic data this is defined in atomic.h, rather than the 
			   modes structure */
			if (modes.iadvanced)
				rdint("write_atomicdata", &write_atomicdata);

			if (write_atomicdata)
				Log("You have opted to save a summary of the atomic data\n");

			get_atomic_data(geo.atomic_filename);

		}

		geo.wcycles = geo.pcycles = 1;
		geo.wcycle = geo.pcycle = 0;

	}

	else if (restart_stat == 1)	/* We want to continue a previous run */
	{
		Log("Continuing a previous run of %s \n", root);
		strcpy(old_windsavefile, root);
		strcat(old_windsavefile, ".wind_save");
		if (wind_read(old_windsavefile) < 0)
		{
			Error("python: Unable to open %s\n", old_windsavefile);	// program will exit if unable to read the file
			exit(0);
		}
		w = wmain;
		geo.wind_type = 2;		// We read the data from a file
		xsignal(root, "%-20s Read %s\n", "COMMENT", old_windsavefile);

		if (geo.pcycle > 0)
		{
			spec_read(specsavefile);
			xsignal(root, "%-20s Read %s\n", "COMMENT", specsavefile);
		}
	}




	/* Get the remainder of the data.  Note that this is done whether or not the windsave file was read in */

	/* 140907 - ksl - Although photons_per_cycle is really an integer, read in as a double so it is easier for input */

	x = 100000;
	rddoub("photons_per_cycle", &x);
	NPHOT = x;					// NPHOT is photons/cycle

#ifdef MPI_ON
	Log("Photons per cycle per MPI task will be %d\n", NPHOT / np_mpi_global);

	NPHOT /= np_mpi_global;
#endif

	rdint("Ionization_cycles", &geo.wcycles);

	rdint("spectrum_cycles", &geo.pcycles);


	Debug("Test %d %d \n", geo.wcycles, geo.pcycles);

	if (geo.wcycles == 0 && geo.pcycles == 0)
		exit(0);				// There is really nothing to do!

	/* Allocate the memory for the photon structure now that NPHOT is established */

	p = (PhotPtr) calloc(sizeof(p_dummy), NPHOT);

	if (p == NULL)
	{
		Error("There is a problem in allocating memory for the photon structure\n");
		exit(0);
	}


	if (geo.wind_type != 2)
	{
		/* Define the coordinate system for the grid and allocate memory for the wind structure */
		rdint("Coord.system(0=spherical,1=cylindrical,2=spherical_polar,3=cyl_var)", &geo.coord_type);

		rdint("Wind.dim.in.x_or_r.direction", &geo.ndim);
		if (geo.coord_type)
		{
			rdint("Wind.dim.in.z_or_theta.direction", &geo.mdim);
			if (geo.mdim < 4)
			{
				Error("python: geo.mdim must be at least 4 to allow for boundaries\n");
				exit(0);
			}
		}
		else
			geo.mdim = 1;

	}

	/* 130405 ksl - Check that NDIM_MAX is greater than NDIM and MDIM.  */

	if ((geo.ndim > NDIM_MAX) || (geo.mdim > NDIM_MAX))
	{
		Error("NDIM_MAX %d is less than NDIM %d or MDIM %d. Fix in python.h and recompile\n", NDIM_MAX, geo.ndim, geo.mdim);
		exit(0);
	}


	/* If we are in advanced then allow the user to modify scale lengths */
	if (modes.iadvanced)
	{
		rdint("adjust_grid(0=no,1=yes)", &modes.adjust_grid);

		if (modes.adjust_grid)
		{
			Log("You have opted to adjust the grid scale lengths\n");
			rddoub("geo.xlog_scale", &geo.xlog_scale);
			if (geo.coord_type)
				rddoub("geo.zlog_scale", &geo.zlog_scale);
		}
	}


	// 080808 - 62 - Ionization section has been cleaned up -- ksl
	/* ??? ksl - Acoording to line 110 of ioniztion. option 4 is LTE with SIM_correction.  It would be good to know what this is
	   actually.  Note that pairwise is the appraoch which cboses between pairwise_bb, and pairwise_pow. Normally, any of the
	   pairwise options should force use of a banding option with a broad set of bands */

	rdint("Wind_ionization(0=on.the.spot,1=LTE,2=fixed,3=recalc_bb,6=pairwise_bb,7=pairwise_pow)", &geo.ioniz_mode);

	if (geo.ioniz_mode == IONMODE_FIXED)
	{
		rdstr("Fixed.concentrations.filename", &geo.fixed_con_file[0]);
	}
	if (geo.ioniz_mode == 4 || geo.ioniz_mode == 5 || geo.ioniz_mode > 9)	/* NSH CLOUDY test - remove once done */
	{
		Log("The allowed ionization modes are 0, 1, 2, 3, 6, 7\n");
		Error("Unknown ionization mode %d\n", geo.ioniz_mode);
		exit(0);
	}

	/* Normally, geo.partition_mode is set to -1, which means that partition functions are calculated to take full advantage of the 
	   data file.  This means that in calculating the partition functions, the information on levels and their multiplicities is
	   taken into account.  */

	geo.partition_mode = -1;	// ?? Stuart, is there a reason not to move this earlier so it does not affect restart


	rdint
		("Line_transfer(0=pure.abs,1=pure.scat,2=sing.scat,3=escape.prob,6=macro_atoms,7=macro_atoms+aniso.scattering)",
		 &geo.line_mode);

	/* ?? ksl Next section seems rather a kluge.  Why don't we specifty the underlying variables explicitly It also seems likely
	   that we have mixed usage of some things, e.g geo.rt_mode and geo.macro_simple */

	/* JM 1406 -- geo.rt_mode and geo.macro_simple control different things. geo.rt_mode controls the radiative transfer and
	   whether or not you are going to use the indivisible packet constraint, so you can have all simple ions, all macro-atoms or
	   a mix of the two. geo.macro_simple just means one can turn off the full macro atom treatment and treat everything as
	   2-level simple ions inside the macro atom formalism */

	/* For now handle scattering as part of a hidden line transfermode ?? */
	if (geo.line_mode == 4)
	{
		geo.scatter_mode = 1;	// Turn on anisotropic scattering
		geo.line_mode = 3;		// Drop back to escape probabilities
		geo.rt_mode = 1;		// Not macro atom (SS)
	}
	else if (geo.line_mode == 5)
	{
		geo.scatter_mode = 2;	// Thermal trapping model
		geo.line_mode = 3;		// Single scattering model is best for this mode
		geo.rt_mode = 1;		// Not macro atom (SS) 
	}
	else if (geo.line_mode == 6)
	{
		geo.scatter_mode = 0;	// isotropic
		geo.line_mode = 3;		// Single scattering
		geo.rt_mode = 2;		// Identify macro atom treatment (SS)
		geo.macro_simple = 0;	// We don't want the all simple case (SS)
	}
	else if (geo.line_mode == 7)
	{
		geo.scatter_mode = 2;	// thermal trapping
		geo.line_mode = 3;		// Single scattering
		geo.rt_mode = 2;		// Identify macro atom treatment (SS)
		geo.macro_simple = 0;	// We don't want the all simple case (SS)
	}
	else if (geo.line_mode == 8)
	{
		geo.scatter_mode = 0;	// isotropic
		geo.line_mode = 3;		// Single scattering
		geo.rt_mode = 2;		// Identify macro atom treatment i.e. indivisible packets
		geo.macro_simple = 1;	// This is for test runs with all simple ions (SS)
	}
	else if (geo.line_mode == 9)	// JM 1406 -- new mode, as mode 7, but scatter mode is 1
	{
		geo.scatter_mode = 1;	// anisotropic scatter mode 1
		geo.line_mode = 3;		// Single scattering
		geo.rt_mode = 2;		// Identify macro atom treatment 
		geo.macro_simple = 0;	// We don't want the all simple case 
	}
	else
	{
		geo.scatter_mode = 0;	// isotropic
		geo.rt_mode = 1;		// Not macro atom (SS)
	}

	thermal_opt = 0;			/* NSH 131213 Set the option to zero - the default. The lines allow allow the user to turn off
								   mechanisms that affect the thermal balance. Adiabatic is the only one implemented to start off
								   with. */

	rdint("Thermal_balance_options(0=everything.on,1=no.adiabatic)", &thermal_opt);

	if (thermal_opt == 1)
	{
		geo.adiabatic = 0;
	}

	else if (thermal_opt > 1 || thermal_opt < 0)
	{
		Error("Unknown thermal balance mode %d\n", thermal_opt);
		exit(0);
	}


	/* 57h -- Next line prevents bf calculation of macro_estimaters when no macro atoms are present.  */

	if (nlevels_macro == 0)
		geo.macro_simple = 1;	// Make everything simple if no macro atoms -- 57h

	// SS - initalise the choice of handling for macro pops.
	if (geo.wind_type == 2)
	{
		geo.macro_ioniz_mode = 1;	// Now that macro atom properties are available for restarts
	}
	else
	{
		geo.macro_ioniz_mode = 0;
	}

	// Establish the overall system type - Added for python_69 to allow qso's have different inputs
	// Note - ksl - What happened to the possibility of a true single star with no disk - 110914

	rdint("System_type(0=star,1=binary,2=agn)", &geo.system_type);


	// Determine what radiation sources there are.  Note that most of these values are initilized in init_geo

	if (geo.system_type != SYSTEM_TYPE_AGN)
	{							/* If is a stellar system */
		rdint("Star_radiation(y=1)", &geo.star_radiation);
		rdint("Disk_radiation(y=1)", &geo.disk_radiation);
		rdint("Boundary_layer_radiation(y=1)", &geo.bl_radiation);
		rdint("Wind_radiation(y=1)", &geo.wind_radiation);
		geo.agn_radiation = 0;	// So far at least, our star systems don't have a BH
	}
	else						/* If it is an AGN */
	{
		geo.star_radiation = 0;	// 70b - AGN do not have a star at the center */
		rdint("Disk_radiation(y=1)", &geo.disk_radiation);
		geo.bl_radiation = 0;
		rdint("Wind_radiation(y=1)", &geo.wind_radiation);
		geo.agn_radiation = 1;
		rdint("QSO_BH_radiation(y=1)", &geo.agn_radiation);
	}

	if (!geo.star_radiation && !geo.disk_radiation && !geo.bl_radiation && !geo.bl_radiation && !geo.agn_radiation)
	{
		Error("python: No radiation sources so nothing to do but quit!\n");
		exit(0);
	}

	/* 
	   With the macro atom approach we won't want to generate photon bundles in the wind so switch it off here. (SS) */

	if (geo.rt_mode == 2)
	{
		Log("python: Using Macro Atom method so switching off wind radiation.\n");
		geo.wind_radiation = 0;
	}


	/* 080517 - ksl - Reassigning bb to -1, etc is to make room for reading in model grids, but complicates what happens if one
	   tries to restart a model.  This needs to be updated so one can re-read the geo file, proabbly by defining variaables BB
	   etc, and then by checking whether or not the type is assigned to BB or read in as 0.  Also need to store each of these model 
	   list names in geo structure. */

	get_spectype(geo.star_radiation, "Rad_type_for_star(0=bb,1=models)_to_make_wind", &geo.star_ion_spectype);

	get_spectype(geo.disk_radiation, "Rad_type_for_disk(0=bb,1=models)_to_make_wind", &geo.disk_ion_spectype);

	get_spectype(geo.bl_radiation, "Rad_type_for_bl(0=bb,1=models,3=pow)_to_make_wind", &geo.bl_ion_spectype);
	get_spectype(geo.agn_radiation,
				 "Rad_type_for_agn(0=bb,1=models,3=power_law,4=cloudy_table)_to_make_wind", &geo.agn_ion_spectype);


	/* 130621 - ksl - This is a kluge to add a power law to stellar systems.  What id done is to remove the bl emission, which we
	   always assume to some kind of temperature driven source, and replace it with a power law source

	   Note that the next 3 or 4 lines just tell you that there is supposed to be a power law source.  They don't teel you what the 
	   parameters are. */

	if (geo.bl_ion_spectype == SPECTYPE_POW)
	{
		geo.agn_radiation = 1;
		geo.agn_ion_spectype = SPECTYPE_POW;
		geo.bl_radiation = 0;
		Log("Trying to make a start with a power law boundary layer\n");
	}
	else
	{
		Log("NOt Trying to make a start with a power law boundary layer %d\n", geo.bl_ion_spectype);
	}



	if (geo.wind_type == 2)
	{
		disk_illum = geo.disk_illum;
	}


	if (geo.wind_type != 2)		// Start of block to define a model for the first time
	{

		/* Describe the basic binary star system */

		geo.mstar /= MSOL;		// Convert to MSOL for ease of data entry
		rddoub("mstar(msol)", &geo.mstar);
		geo.mstar *= MSOL;

		/* If a BH we want geo.rstar to be at least as large as the last stable orbit for a non-rotating BH */

		if (geo.system_type == SYSTEM_TYPE_AGN)
		{
			geo.rstar = 6. * G * geo.mstar / (C * C);	// correction - ISCO is 6x Rg NSH 121025
		}

		rddoub("rstar(cm)", &geo.rstar);


		geo.r_agn = geo.rstar;	/* At present just set geo.r_agn to geo.rstar */
		geo.rstar_sq = geo.rstar * geo.rstar;
		if (geo.star_radiation)
			rddoub("tstar", &geo.tstar);

		lstar = 4 * PI * geo.rstar * geo.rstar * STEFAN_BOLTZMANN * pow(geo.tstar, 4.);


		/* Describe the secondary if that is required */

		if (geo.system_type == SYSTEM_TYPE_BINARY)	/* It's a binary system */
		{

			geo.m_sec /= MSOL;	// Convert units for ease of data entry
			rddoub("msec(msol)", &geo.m_sec);
			geo.m_sec *= MSOL;

			geo.period /= 3600.;	// Convert units to hours for easy of data entry
			rddoub("period(hr)", &geo.period);
			geo.period *= 3600.;	// Put back to cgs immediately 
		}

		/* Describe the disk */

		rdint("disk.type(0=no.disk,1=standard.flat.disk,2=vertically.extended.disk)", &geo.disk_type);
		if (geo.disk_type)		/* Then a disk exists and it needs to be described */
		{
			// if (geo.disk_radiation) /*NSH 130906 - Commented out this if loop. It was causing problems with restart - bug #44
			// {
			geo.disk_mdot /= (MSOL / YR);	// Convert to msol/yr to simplify input
			rddoub("disk.mdot(msol/yr)", &geo.disk_mdot);
			geo.disk_mdot *= (MSOL / YR);
			disk_illum = 0;
			rdint("Disk.illumination.treatment(0=no.rerad,1=high.albedo,2=thermalized.rerad,3=analytic)", &disk_illum);
			rdint("Disk.temperature.profile(0=standard;1=readin)", &geo.disk_tprofile);
			if (geo.disk_tprofile == 1)
			{
				rdstr("T_profile_file", tprofile);
			}
			// }
			// else
			// {
			// geo.disk_mdot = 0;
			// disk_illum = 0;
			// }

			/* 04aug ksl ??? Until everything is initialized we need to stick to a simple disk, while teff is being set up..  This 
			   is because some of the models, e.g. knigge have wind structures that depend on teff. * 080518 - ksl - this is quite
			   confusing.  I understand that the KWD models have base velocities that are affected by t_eff, but we have not done
			   anything yet.  Possible this is a consideration for restart, but I would have guessed we recalculated all of that,
			   and in any event this is within the block to reset things.  this is likely a problem of some sort

			   However, the next line does force the illum to 0 while initialization is going on, unless ?? the illum is 3 */

			geo.disk_illum = 0;
			if (disk_illum == 3)	// 080518 - And why is it different in the analytic case?
			{
				geo.disk_illum = 3;
			}

			/* Set a default for diskrad for an AGN */
			if (geo.system_type == SYSTEM_TYPE_AGN)
			{
				geo.diskrad = 100. * geo.r_agn;
			}

			rddoub("disk.radmax(cm)", &geo.diskrad);
			Log("geo.diskrad  %e\n", geo.diskrad);

			/* If diskrad <= geo.rstar set geo.disk_type = 0 to make any disk transparent anyway. */

			if (geo.diskrad < geo.rstar)
			{
				Log("Disk radius is less than star radius, so assuming no disk)\n");
				geo.disk_type = 0;
			}

			if (geo.disk_type == 2)
			{					/* Get the additional variables need to describe a vertically extended disk */
				rddoub("disk.z0(fractional.height.at.diskrad)", &geo.disk_z0);
				rddoub("disk.z1(powerlaw.index)", &geo.disk_z1);
			}
		}

		else
		{						/* There is no disk so set variables accordingly */
			geo.disk_radiation = 0;
			geo.diskrad = 0;
		}


		/* Describe the boundary layer */

		// OLD 130622 if (geo.bl_radiation ) Change made to allow a power law boundary layer
		if (geo.bl_radiation && geo.bl_ion_spectype != SPECTYPE_POW)
		{
			xbl = geo.lum_bl = 0.5 * G * geo.mstar * geo.disk_mdot / geo.rstar;

			rddoub("lum_bl(ergs/s)", &geo.lum_bl);
			Log("OK, the bl lum will be about %.2e the disk lum\n", geo.lum_bl / xbl);
			rddoub("t_bl", &geo.t_bl);
		}
		else
		{
			geo.lum_bl = 0;
			geo.t_bl = 0;
		}

		/* Describe the agn */

		if (geo.agn_radiation && geo.system_type == SYSTEM_TYPE_AGN)	/* This peculiar line is to enamble us to add a star with a 
																		   power law component */
		{
			xbl = geo.lum_agn = 0.5 * G * geo.mstar * geo.disk_mdot / geo.r_agn;

			/* If there is no disk, initilize geo.lum to the luminosity of a star */
			if (geo.disk_type == 0)
			{
				geo.lum_agn = lstar;
			}

			// At present we have set geo.r_agn = geo.rstar, and encouraged the user
			// set the default for the radius of the BH to be 6 R_Schwartschild.
			// rddoub("R_agn(cm)",&geo.r_agn);

			rddoub("lum_agn(ergs/s)", &geo.lum_agn);
			Log("OK, the agn lum will be about %.2e the disk lum\n", geo.lum_agn / xbl);
			geo.alpha_agn = (-1.5);
			rddoub("agn_power_law_index", &geo.alpha_agn);

			/* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity. This is only used in the 
			   sim correction factor for the first time through. Afterwards, the photons are used to compute the sim parameters. */

			geo.const_agn =
				geo.lum_agn / (((pow(2.42e18, geo.alpha_agn + 1.)) - pow(4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
			Log("AGN Input parameters give a power law constant of %e\n", geo.const_agn);

			if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)	/* NSH 0412 - option added to allow direct comparison with cloudy power 
															   law table option */
			{
				geo.agn_cltab_low = 1.0;
				geo.agn_cltab_hi = 10000;
				rddoub("low_energy_break(ev)", &geo.agn_cltab_low);	/* lo frequency break - in ev */
				rddoub("high_energy_break(ev)", &geo.agn_cltab_hi);
				geo.agn_cltab_low_alpha = 2.5;	// this is the default value in cloudy
				geo.agn_cltab_hi_alpha = -2.0;	// this is the default value in cloudy
			}
		}
		else if (geo.agn_radiation)	/* We want to add a power law to something other than an AGN */
		{
			xbl = geo.lum_agn = 0.5 * G * geo.mstar * geo.disk_mdot / geo.r_agn;

			// At present we have set geo.r_agn = geo.rstar, and encouraged the user
			// set the default for the radius of the BH to be 6 R_Schwartschild.
			// rddoub("R_agn(cm)",&geo.r_agn);

			rddoub("lum_agn(ergs/s)", &geo.lum_agn);
			Log("OK, the agn lum will be about %.2e the disk lum\n", geo.lum_agn / xbl);
			geo.alpha_agn = (-1.5);
			rddoub("agn_power_law_index", &geo.alpha_agn);

			/* Computes the constant for the power law spectrum from the input alpha and 2-10 luminosity. This is only used in the 
			   sim correction factor for the first time through. Afterwards, the photons are used to compute the sim parameters. */

			geo.const_agn =
				geo.lum_agn / (((pow(2.42e18, geo.alpha_agn + 1.)) - pow(4.84e17, geo.alpha_agn + 1.0)) / (geo.alpha_agn + 1.0));
			Log("AGN Input parameters give a power law constant of %e\n", geo.const_agn);

			if (geo.agn_ion_spectype == SPECTYPE_CL_TAB)	/* NSH 0412 - option added to allow direct comparison with cloudy power 
															   law table option */
			{
				geo.agn_cltab_low = 1.0;
				geo.agn_cltab_hi = 10000;
				rddoub("low_energy_break(ev)", &geo.agn_cltab_low);	/* lo frequency break - in ev */
				rddoub("high_energy_break(ev)", &geo.agn_cltab_hi);
				geo.agn_cltab_low_alpha = 2.5;	// this is the default value in cloudy
				geo.agn_cltab_hi_alpha = -2.0;	// this is the default value in cloudy
			}
		}
		else
		{
			geo.r_agn = 0.0;
			geo.lum_agn = 0.0;
			geo.alpha_agn = 0.0;
			geo.const_agn = 0.0;
		}

		/* Describe the Compton torus */

		/* 70b - ksl - 1108067 - Here we add parameters for the compton torus or blocking region * * Note that the whole flow of
		   this may be a bit odd as it seems as if we have to keep checking for whether * we are modelling an agn * * Note that
		   these calls need to precede the calls below, because we want to keep the compton torus ??? * inside the actual wind, or
		   at least that's what ksl believes on 110809.  ??? */

		/* 140907: ksl - I have effectively commented out any consideration of the compton_torus because it is not debugged.  But
		   one can continue debugging it by setting the DEBUG variatble to true */
		geo.compton_torus = 0;

		/* JM 1411 -- we only ask about the Torus if we've used the -d flag */
		if (modes.iadvanced)
			rdint("Torus(0=no,1=yes)", &geo.compton_torus);

		if (geo.compton_torus)
		{
			rddoub("Torus.rmin(cm)", &geo.compton_torus_rmin);
			rddoub("Torus.rmax(cm)", &geo.compton_torus_rmax);
			rddoub("Torus.height(cm)", &geo.compton_torus_zheight);
			rddoub("Torus.optical_depth", &geo.compton_torus_tau);
			rddoub("Torus.tinit", &geo.compton_torus_te);
			if (geo.compton_torus_tau <= 0)
			{
				geo.compton_torus_tau = 0.001;
				Error("python: A torus with zero optical depth makes no sense. Setting to %f\n", geo.compton_torus_tau);
			}
		}


		/* Describe the wind */

		if (geo.system_type == SYSTEM_TYPE_AGN)
		{
			geo.rmax = 50. * geo.r_agn;
		}

		rddoub("wind.radmax(cm)", &geo.rmax);
		rddoub("wind.t.init", &geo.twind);

		geo.diskrad_sq = geo.diskrad * geo.diskrad;


		/* Now get parameters that are specific to a given wind model

		   Note: When one adds a new model, the only things that should be read in and modified are parameters in geo.  This is in
		   order to preserve the ability to continue a calculation with the same basic wind geometry, without reading in all of the 
		   input parameters. */

		if (geo.wind_type == 1)
		{
			get_stellar_wind_params();
		}
		else if (geo.wind_type == 0)
		{
			get_sv_wind_params();
		}
		else if (geo.wind_type == 3)
		{
			get_proga_wind_params();
		}
		else if (geo.wind_type == 4)
		{
			get_corona_params();
		}
		else if (geo.wind_type == 5)
		{
			get_knigge_wind_params();
		}
		else if (geo.wind_type == 6)
		{
			get_homologous_params();
		}
		else if (geo.wind_type == 7)
		{
			get_yso_wind_params();
		}
		else if (geo.wind_type == 8)
		{
			get_elvis_wind_params();
		}
		else if (geo.wind_type == 9)	// NSH 18/2/11 This is a new wind type to produce a thin shell.
		{
			get_shell_wind_params();
			/* NSH 121219 moved dfudge = (geo.wind_rmax - geo.wind_rmin) / 1000.0; Stop photons getting pushed out of the cell
			   Modified again in python 71b to take account of change in parametrisation of shell wind DFUDGE = dfudge; */
		}
		else if (geo.wind_type != 2)
		{
			Error("python: Unknown wind type %d\n", geo.wind_type);
			exit(0);
		}

		/* Get the filling factor of the wind */
		geo.fill = 1.;
		rddoub("wind.filling_factor(1=smooth,<1=clumped)", &geo.fill);

	}							// End of block to define a model for the first time
	else
	{
		if (geo.disk_type)		/* Then a disk exists and it needs to be described */
		{
			if (geo.disk_radiation)
			{
				rdint("Disk.temperature.profile(0=standard;1=readin)", &geo.disk_tprofile);
				if (geo.disk_tprofile == 1)
				{
					rdstr("T_profile_file", tprofile);
				}
			}
		}
	}

	/* 121219 NSH Set up DFUDGE to be a value that makes some kind of sense given the scale of the wind. Up till py74b2 it was set
	   to be fixed at 1e5, so we ensure that this is a minimum, so any winds of CV type scale will keep the old dfudge, and
	   hopefully look the same. We also need to set defudge slightly differently for the shell wind. */

	if (geo.wind_type == 9)
	{
		dfudge = (geo.wind_rmax - geo.wind_rmin) / 1000.0;
	}
	else
	{
		if (geo.rmax / 1.e10 < 1e5)
		{
			dfudge = 1e5;
			Log("DFUDGE set to minimum value of %e\n", dfudge);
		}
		else
		{
			dfudge = geo.rmax / 1.e10;
			Log("DFUDGE set to %e based on geo.rmax\n", dfudge);
		}
	}

	DFUDGE = dfudge;			// NSH Unsure why exactly this is done this way.





	/* Now define the wind cones generically.  The angles thetamin and thetamax are all defined from the z axis, so that an angle
	   of 0 is a flow that is perpeindicular to to the disk and one that is close to 90 degrees will be parallel to the plane of
	   the disk geo.wind_thetamin and max are defined in the routines that initialize the various wind models, e. g.
	   get_sv_wind_parameters. These have been called at this point.

	   z is the place where the windcone intercepts the z axis dzdr is the slope

	   111124 fixed notes on this - ksl */


	if (geo.wind_thetamin > 0.0)
	{
		windcone[0].dzdr = 1. / tan(geo.wind_thetamin);
		windcone[0].z = (-geo.wind_rho_min / tan(geo.wind_thetamin));
	}
	else
	{
		windcone[0].dzdr = VERY_BIG;
		windcone[0].z = -VERY_BIG;;
	}


	if (geo.wind_thetamax > 0.0)
	{
		windcone[1].dzdr = 1. / tan(geo.wind_thetamax);
		windcone[1].z = (-geo.wind_rho_max / tan(geo.wind_thetamax));
	}
	else
	{
		windcone[1].dzdr = VERY_BIG;
		windcone[1].z = -VERY_BIG;;
	}

	/* NSH 130821 broken out into a seperate routine added these lines to fix bug41, where the cones are never defined for an
	   rtheta grid if the model is restarted */

	if (geo.coord_type == RTHETA && geo.wind_type == 2)	// We need to generate an rtheta wind cone if we are restarting
	{
		rtheta_make_cones(wmain);
	}

	geo.rmax_sq = geo.rmax * geo.rmax;

	/* Calculate additional parameters associated with the binary star system */

	if (geo.system_type == SYSTEM_TYPE_BINARY)
		binary_basics();

	/* Check that the parameters which have been supplied for the star, disk and boundary layer will allow generation of photons
	   where that is appropriate */

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
		Log("There is a star which radiates\n");
	else
		Log("The star in the system does not radiate\n");

	if (!geo.disk_type)
		Log("There is no disk in the system \n");
	else if (!geo.disk_radiation)
		Log("The disk exists, but only absorbs photons\n");
	else
		Log("There is a disk which radiates and absorbs\n");

	if (geo.bl_radiation)
		Log("There is a boundary layer which radiates\n");
	else
		Log("There is no boundary layer\n");

	if (geo.agn_radiation)
		Log("There is a BH  which radiates\n");
	else
		Log("There is no BH \n");

	/* Describe the spectra which will be extracted and the way it will be extracted */

	/* First initialise things to semi-reasonable values */

	nangles = 4;
	angle[0] = 10;
	angle[1] = 30.;
	angle[2] = 60.;
	angle[3] = 80.;
	for (n = 4; n < NSPEC; n++)
		angle[n] = 45;
	for (n = 0; n < NSPEC; n++)
	{
		phase[n] = 0.5;
		scat_select[n] = 1000;
		top_bot_select[n] = 0;
	}
	swavemin = 1450;
	swavemax = 1650;

	/* These two variables have to do with what types of spectra are created n the spectrum files. They are not associated with the 
	   nature of the spectra that are generated by say the boundary layer */

	select_extract = 1;
	select_spectype = 1;

	/* Completed initialization of this section.  Note that get_spectype uses the source of the ratiation and then value given to
	   return a spectrum type. The output is not the same number as one inputs. It's not obvious that this is a good idea. */

	if (geo.pcycles > 0)
	{

		get_spectype(geo.star_radiation, "Rad_type_for_star(0=bb,1=models,2=uniform)_in_final_spectrum", &geo.star_spectype);


		get_spectype(geo.disk_radiation, "Rad_type_for_disk(0=bb,1=models,2=uniform)_in_final_spectrum", &geo.disk_spectype);


		get_spectype(geo.bl_radiation, "Rad_type_for_bl(0=bb,1=models,2=uniform)_in_final_spectrum", &geo.bl_spectype);

		geo.agn_spectype = 3;
		get_spectype(geo.agn_radiation, "Rad_type_for_agn(3=power_law,4=cloudy_table)_in_final_spectrum", &geo.agn_spectype);



		rddoub("spectrum_wavemin", &swavemin);
		rddoub("spectrum_wavemax", &swavemax);
		if (swavemin > swavemax)
		{
			swavemax = swavemin;
			swavemin = swavemax;
		}

		/* SS June 04: convert these to frequencies and store for use in computing macro atom and k-packet emissivities. */

		em_rnge.fmin = C / (swavemax * 1.e-8);
		em_rnge.fmax = C / (swavemin * 1.e-8);

		geo.matom_radiation = 0;	// initialise for ionization cycles - don't use pre-computed emissivities for macro-atom
									// levels/ k-packets.

		/* Note: Below here many of the variables which are read in are not currently part of geo stucture */

		rdint("no_observers", &nangles);

		if (nangles < 1 || nangles > NSPEC)
		{
			Error("no_observers %d should not be > %d or <0\n", nangles, NSPEC);
			exit(0);
		}

		for (n = 0; n < nangles; n++)
			rddoub("angle(0=pole)", &angle[n]);


    
      for (n = 0; n < nangles; n++)
	rddoub ("angle(0=pole)", &angle[n]);

      /* 05apr-ksl-56--For diagnositic reasons I have left questions regarding phase
       * even for systems which are not binaries.  Phase 0 in this case corresponds to
       * an extraction direction which is in the xz plane
       */
      /* JM 1502 -- change this so we only ask for phase if the system is a binary -- see #137 */

      if (geo.system_type == SYSTEM_TYPE_BINARY)
      {
        
        for (n = 0; n < nangles; n++)
	      rddoub ("phase(0=inferior_conjunction)", &phase[n]);
      }
      else
      	Log("No phase information needed as system type %i is not a binary\n",
      		 geo.system_type);


		if (modes.iadvanced)
		{
			strcpy(yesno, "n");
			rdstr("Select_specific_no_of_scatters_in_spectra(y/n)", yesno);
			if (yesno[0] == 'y')
			{
				Log("OK n>MAXSCAT->all; 0<=n<MAXSCAT -> n scatters; n<0 -> >= |n| scatters\n");
				for (n = 0; n < nangles; n++)
				{
					rdint("Select_scatters", &scat_select[n]);
				}
			}
			strcpy(yesno, "n");
			rdstr("Select_photons_by_position(y/n)", yesno);
			if (yesno[0] == 'y')
			{
				Log("OK 0->all; -1 -> below; 1 -> above the disk, 2 -> specific location in wind\n");
				for (n = 0; n < nangles; n++)
				{
					rdint("Select_location", &top_bot_select[n]);
					if (top_bot_select[n] == 2)
					{
						Log("Warning: Make sure that position will be in wind, or no joy will be obtained\n");
						rddoub("rho(cm)", &rho_select[n]);
						rddoub("z(cm)", &z_select[n]);
						rddoub("azimuth(deg)", &az_select[n]);
						rddoub("r(cm)", &r_select[n]);

					}
				}
			}
		}
	}

	/* Select the units of the output spectra.  This is always needed */

	rdint("spec.type(flambda(1),fnu(2),basic(other)", &select_spectype);
	if (select_spectype == 1)
	{
		Log("OK, generating flambda at 100pc\n");
	}
	else if (select_spectype == 2)
	{
		Log("OK, generating fnu at 100 pc\n");
	}
	else
		Log("OK, basic Monte Carlo spectrum\n");



	/* 57h -- New section of inputs to provide more control over how the program is run -- 07jul -- ksl */

	istandard = 1;
	SMAX_FRAC = 0.5;
	DENSITY_PHOT_MIN = 1.e-10;
	keep_photoabs = 1;

	/* 141116 - ksl - Made care factors and advanced command as this is clearly somethng that is diagnostic */

	if (modes.iadvanced)
	{
		rdint("Use.standard.care.factors(1=yes)", &istandard);

		if (!istandard)
		{
			rddoub("Fractional.distance.photon.may.travel", &SMAX_FRAC);
			rddoub("Lowest.ion.density.contributing.to.photoabsorption", &DENSITY_PHOT_MIN);
			rdint("Keep.photoabs.during.final.spectrum(1=yes)", &keep_photoabs);
		}
	}





	/* 081221 - 67c - Establish limits on the frequency intervals to be used by the ionization cycles and the fraquency bands for
	   stratified sampling.  Changes here were made to allow more control over statified sampling, since we have expanded the
	   temperature ranges of the types of systems we would like to calculate.  This section of inputs might logically go earlier in 
	   the code, but was put here so it would add on to existing .pf files.  It would be reasonble to consider moving it to a more
	   logical location */


	/* Determine the frequency range which will be used to establish the ionization balance of the wind */

	// Note that bands_init asks .pf file or user what kind of banding is desired 

	bands_init(-1, &xband);

	/* if we have changed min and max in bands_init, we need to make sure this is reflected in the frequency bounds */
	freqmin = xband.f1[0];
	freqmax = xband.f2[xband.nbands - 1];

	/* 1112 - 71 - ksl Next routine sets up the frequencies that are used for charactizing the spectrum in a cell These need to be
	   coordinated with the bands that are set up for spectral gneration */
	freqs_init(freqmin, freqmax);


	if (modes.iadvanced)
	{
		/* Do we require extra diagnostics or not */
		rdint("Extra.diagnostics(0=no,1=yes) ", &modes.diag_on_off);

		if (modes.diag_on_off)
		{
			get_extra_diagnostics();
		}
	}





	/* Wrap up and save all the inputs */

	if (strncmp(root, "mod", 3) == 0)
		cpar("mod.pf");
	else if (strncmp(root, "dummy", 5) == 0)
	{
		cpar("dummy.pf");
		exit(0);
	}
	else if (opar_stat == 1)
	{
		cpar(input);
	}
	else
		cpar("python.pf");

	/* OK all inputs have been obtained at this point and the inputs have been copied to "mod.pf" or "python.pf" */


	/* INPUTS ARE FINALLY COMPLETE */


	/* Next line finally defines the wind if this is the initial time this model is being run */
	if (geo.wind_type != 2)		// Define the wind and allocate the arrays the first time
		define_wind();
	// Do not reinit if you want to use old windfile

	w = wmain;

	if (modes.save_cell_stats)
	{
		/* Open a diagnostic file or files.  These are all fixed files */
		open_diagfile();
	}

	/* initialize the random number generator */
	// srand( (n=(unsigned int) clock())); 
	srand(1084515760 + (13 * rank_global));

	/* 68b - 0902 - ksl - Start with photon history off */

	phot_hist_on = 0;

	/* If required, read in a non-standard disk temperature profile */

	if (geo.disk_tprofile == 1)
	{
		read_non_standard_disk_profile(tprofile);
	}



	/* The next section sets up a structure qdisk to record the effects of illumination on the disk.  disk_init is called primarily 
	   to get a defined set of annular rings which are kept throughout the ionization calculation.  A second structure qdisk is
	   needed because in the process of generating photons in various bands the annular rings are changed disk_init calculates the 
	   flux from the disk in the energy range set by freqmin and freqmax, and uses is this to identify the position of the rings in 
	   the disk, so that each ring contributes the same amount to the flux */


	disk_init(geo.rstar, geo.diskrad, geo.mstar, geo.disk_mdot, freqmin, freqmax, 0, &geo.f_disk);

	qdisk_init();				/* Initialize a disk qdisk to store the information about photons impinging on the disk */

	/* 04aug -- ksl -- now that everything is initialized, we set geo.disk_illum 080518 - ksl - I believe that the reason for this 
	   somewhat weird logic is to assure that models (e.g corona and knigge) where the base wind velocity depends on teff are not
	   altered by illumination, but since no photons have been transported at this stage, it's hard to see why that would matter. */

	geo.disk_illum = disk_illum;

	xsignal(root, "%-20s Finished initialization for %s\n", "NOK", root);
	check_time(root);

#ifdef MPI_ON
	/* Since the wind is now set up can work out the length big arrays to help with the MPI reductions of the spectra the variables 
	   for the estimator arrays are set up in the subroutines themselves */
	ioniz_spec_helpers = 2 * MSPEC * NWAVE;	// we need space for log and lin spectra for MSPEC XNWAVE
	spec_spec_helpers = (NWAVE * (MSPEC + nangles));	// We need space for NWAVE wavelengths for nspectra, which will eventually
														// equal nangles + MSPEC

#endif



	/* XXXX - THE CALCULATION OF THE IONIZATION OF THE WIND */

	geo.ioniz_or_extract = 1;	// SS July 04 - want to compute MC estimators during ionization cycles
	// 1 simply implies we are in the ionization section of the code
	// and allows routines to act accordinaly.

	/* 67 -ksl- geo.wycle will start at zero unless we are completing an old run */

	/* XXXX - BEGINNING OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */

	if (geo.wcycle == geo.wcycles)
		xsignal(root, "%-20s No ionization needed: wcycles(%d)==wcyeles(%d)\n", "COMMENT", geo.wcycle, geo.wcycles);
	else
	{
		geo.pcycle = 0;			/* Set the spectrum cycles executed to 0, because we are going to modify the wind and hence any
								   previously calculated spectra must be recreated */
	}


#ifdef MPI_ON 					//SWM - If using MPI, dump each thread
	delay_dump_prep( delay_dump_file, 0, NSPEC -1, rank_global); //Currently hardcoded to last spectrum
#else 							//SWM - Dump master
	delay_dump_prep( delay_dump_file, NSPEC -1, restart_stat, -1); //Currently hardcoded to last spectrum
#endif

	while (geo.wcycle < geo.wcycles)
	{							/* This allows you to build up photons in bunches */

		xsignal(root, "%-20s Starting %d of %d ionization cycle \n", "NOK", geo.wcycle, geo.wcycles);


		Log("!!Python: Begining cycle %d of %d for defining wind\n", geo.wcycle, geo.wcycles);
		Log_flush();			/* NH June 13 Added call to flush logfile */

		/* Initialize all of the arrays, etc, that need initialization for each cycle */

		spectrum_init(freqmin, freqmax, nangles, angle, phase, scat_select,
					  top_bot_select, select_extract, rho_select, z_select, az_select, r_select);


		wind_rad_init();		/* Zero the parameters pertaining to the radiation field */



		if (modes.ispy)
			ispy_init("python", geo.wcycle);


		geo.n_ioniz = 0.0;
		geo.lum_ioniz = 0.0;
		ztot = 0.0;				/* ztot is the luminosity of the disk multipled by the number of cycles, which is used by
								   save_disk_heating */

		/* JM 1409 -- We used to execute subcycles here, but these have been removed */

		if (!geo.wind_radiation || (geo.wcycle == 0 && geo.wind_type != 2))
			iwind = -1;			/* Do not generate photons from wind */
		else
			iwind = 1;			/* Create wind photons and force a reinitialization of wind parms */

		/* Create the photons that need to be transported through the wind NPHOT is the number of photon bundles which will equal
		   the luminosity; 0 => for ionization calculation */


		/* JM 130306 need to convert photons_per_cycle to double precision for define_phot */
		/* ksl 130410 - This is needed here not because we expect photons per cycle to exceed the size of an integer, but because
		   of the call to define phot in the spectrum cycle, which can exceed this */
		/* JM 1409 photons_per_cycle has been removed in favour of NPHOT */

		nphot_to_define = (long)NPHOT;

		define_phot(p, freqmin, freqmax, nphot_to_define, 0, iwind, 1);

		/* Zero the arrays that store the heating of the disk */

		/* 080520 - ksl - There is a conundrum here.  One should really zero out the quantities below each time the wind structure 
		   is updated.  But relatively few photons hit the disk under normal situations, and therefore the statistcs are not very
		   good. */

		/* 130213 JM -- previously this was done before define_phot, which meant that the ionization state was never computed with
		   the heated disk */

		for (n = 0; n < NRINGS; n++)
		{
			qdisk.heat[n] = qdisk.nphot[n] = qdisk.w[n] = qdisk.ave_freq[n] = 0;
		}



		photon_checks(p, freqmin, freqmax, "Check before transport");

		wind_ip();


		zz = 0.0;
		for (nn = 0; nn < NPHOT; nn++)
		{
			zz += p[nn].w;
		}

		Log("!!python: Total photon luminosity before transphot %18.12e\n", zz);
		Log_flush();			/* NSH June 13 Added call to flush logfile */
		ztot += zz;				/* Total luminosity in all cycles, used for calculating disk heating */

		/* kbf_need determines how many & which bf processes one needs to considere.  It was introduced as a way to speed up the
		   program.  It has to be recalculated evey time one changes freqmin and freqmax */

		kbf_need(freqmin, freqmax);

		/* NSH 22/10/12 This next call populates the prefactor for free free heating for each cell in the plasma array */
		/* NSH 4/12/12 Changed so it is only called if we have read in gsqrd data */
		if (gaunt_n_gsqrd > 0)
			pop_kappa_ff_array();

		/* Transport the photons through the wind */
		trans_phot(w, p, 0);

		/* Determine how much energy was absorbed in the wind */
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

		Log("!!python: Total photon luminosity after transphot %18.12e (diff %18.12e). Radiated luminosity %18.12e\n",
			zzz, zzz - zz, zze);
		if (geo.rt_mode == 2)
			Log("Luminosity taken up by adiabatic kpkt destruction %18.12e number of packets %d\n", zz_adiab, nn_adiab);

		if (modes.print_windrad_summary)
			wind_rad_summary(w, windradfile, "a");




		photon_checks(p, freqmin, freqmax, "Check after transport");

		spectrum_create(p, freqmin, freqmax, nangles, select_extract);



		/* At this point we should communicate all the useful infomation that has been accummulated on differenet MPI tasks */

#ifdef MPI_ON

		communicate_estimators_para();

		communicate_matom_estimators_para();	// this will return 0 if nlevels_macro == 0
#endif




		if (modes.ispy)
			ispy_close();


		/* Calculate and store the amount of heating of the disk due to radiation impinging on the disk */
		qdisk_save(diskfile, ztot);

		/* Completed writing file describing disk heating */

		Log("!!python: Number of ionizing photons %g lum of ionizing photons %g\n", geo.n_ioniz, geo.lum_ioniz);

		/* This step shoudl be MPI_parallelised too */

		wind_update(w);

		/* In a diagnostic mode save the wind file for each cycle (from thread 0) */

		if (modes.keep_ioncycle_windsaves)
		{
			strcpy(dummy, "");
			sprintf(dummy, "python%02d.wind_save", geo.wcycle);

#ifdef MPI_ON
			if (rank_global == 0)
			{
#endif
				wind_save(dummy);
#ifdef MPI_ON
			}
#endif
			Log("Saved wind structure in %s\n", dummy);
		}


		Log("Completed ionization cycle %d :  The elapsed TIME was %f\n", geo.wcycle, timer());

		Log_silent("Finished creating spectra\n");

		/* Do an MPI reduce to get the spectra all gathered to the master thread */

#ifdef MPI_ON

		gather_spectra_para(ioniz_spec_helpers, MSPEC);

#endif



#ifdef MPI_ON
		if (rank_global == 0)
		{
#endif
			spectrum_summary(wspecfile, "w", 0, 5, 0, 1., 0);
			spectrum_summary(lspecfile, "w", 0, 5, 0, 1., 1);	/* output the log spectrum */

#ifdef MPI_ON
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		phot_gen_sum(photfile, "w");	/* Save info about the way photons are created and absorbed by the disk */

		/* Save everything after each cycle and prepare for the next cycle JM1304: moved geo.wcycle++ after xsignal to record
		   cycles correctly. First cycle is cycle 0. */
		/* NSH1306 - moved geo.wcycle++ back, but moved the log and xsignal statements */


		xsignal(root, "%-20s Finished %d of %d ionization cycle \n", "OK", geo.wcycle, geo.wcycles);
		geo.wcycle++;			// Increment ionisation cycles


		delay_dump(p, NPHOT, NSPEC -1, 0); //SWM - Dump delay tracks from this iteration

		/* NSH 1408 - Save only the windsave file from thread 0, to prevent many processors from writing to the same file. */

#ifdef MPI_ON
		if (rank_global == 0)
		{
#endif
			wind_save(windsavefile);
			Log_silent("Saved wind structure in %s after cycle %d\n", windsavefile, geo.wcycle);
#ifdef MPI_ON
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif




		check_time(root);
		Log_flush();			/* Flush the logfile */

	}							// End of Cycle loop

	/* XXXX - END OF CYCLE TO CALCULATE THE IONIZATION OF THE WIND */

#ifdef MPI_ON 						//SWM- If MPI is on
	delay_dump_finish();				//Each thread dumps to file
	MPI_Barrier(MPI_COMM_WORLD);		//Once all done
	delay_dump_combine(np_mpi_global);	//Combine results
#else 								//Otherwise
	delay_dump_finish(); 				//Just dump
#endif

	Log(" Completed wind creation.  The elapsed TIME was %f\n", timer());


	/* XXXX - THE CALCULATION OF A DETAILED SPECTRUM IN A SPECIFIC REGION OF WAVELENGTH SPACE */

	freqmax = C / (swavemin * 1.e-8);
	freqmin = C / (swavemax * 1.e-8);


	/* Perform the initilizations required to handle macro-atoms during the detailed calculation of the spectrum.

	   Next lines turns off macro atom estimators and other portions of the code that are unnecessary during spectrum cycles.  */

	geo.ioniz_or_extract = 0;

	/* 57h -- 07jul -- Next steps to speed up extraction stage */
	if (!keep_photoabs)
	{
		DENSITY_PHOT_MIN = -1.0;	// Do not calculated photoabsorption in detailed spectrum 
	}

	/* Switch on k-packet/macro atom emissivities SS June 04 */

	if (geo.rt_mode == 2)
	{
		geo.matom_radiation = 1;
	}

	/* Finished initializations required for macro-atom approach */

	/* Calculate and store which bf processess need to be considered in each cell Note that this is not macro-specific but is just
	   to speed the program up. */

	kbf_need(freqmin, freqmax);

	/* XXXX - BEGIN CYCLES TO CREATE THE DETAILED SPECTRUM */

	/* the next section initializes the spectrum array in two cases, for the standard one where one is calulating the spectrum for
	   the first time and in the somewhat abnormal case where additional ionization cycles were calculated for the wind */

	if (geo.pcycle == 0)
	{
		spectrum_init(freqmin, freqmax, nangles, angle, phase, scat_select,
					  top_bot_select, select_extract, rho_select, z_select, az_select, r_select);

		/* 68b - zero the portion of plasma main that records the numbers of scatters by each ion in a cell */

		zero_scatters();

	}

	/* the next condition should really when one has nothing more to do */

	else if (geo.pcycle >= geo.pcycles)
		xsignal(root, "%-20s No spectrum   needed: pcycles(%d)==pcycles(%d)\n", "COMMENT", geo.pcycle, geo.pcycles);

	else
	{
		/* Then we are restarting a run with more spectral cycles, but we have already completed some. The memory for the spectral 
		   arrays should already have been allocated, and the spectrum was initialised on the original run, so we just need to
		   renormalise the saved spectrum */
		/* See issue #134 (JM) */
		if (restart_stat == 0)
			Error("Not restarting, but geo.pcycle = %i and trying to renormalise!\n", geo.pcycle);

		spectrum_restart_renormalise(nangles);
	}


	while (geo.pcycle < geo.pcycles)
	{							/* This allows you to build up photons in bunches */

		xsignal(root, "%-20s Starting %d of %d spectral cycle \n", "NOK", geo.pcycle, geo.pcycles);

		if (modes.ispy)
			ispy_init("python", geo.pcycle + 1000);


		Log("!!Cycle %d of %d to calculate a detailed spectrum\n", geo.pcycle, geo.pcycles);
		Log_flush();			/* NSH June 13 Added call to flush logfile */
		if (!geo.wind_radiation)
			iwind = -1;			/* Do not generate photons from wind */
		else if (geo.pcycle == 0)
			iwind = 1;			/* Create wind photons and force a reinitialization of wind parms */
		else
			iwind = 0;			/* Create wind photons but do not force reinitialization */

		/* Create the initial photon bundles which need to be trannsported through the wind

		   For the detailed spectra, NPHOT*pcycles is the number of photon bundles which will equal the luminosity, 1 implies that 
		   detailed spectra, as opposed to the ionization of the wind is being calculated

		   JM 130306 must convert NPHOT and pcycles to double precision variable nphot_to_define

		 */

		nphot_to_define = (long)NPHOT *(long)geo.pcycles;
		define_phot(p, freqmin, freqmax, nphot_to_define, 1, iwind, 0);

		for (icheck = 0; icheck < NPHOT; icheck++)
		{
			if (sane_check(p[icheck].freq))
			{
				Error("python after define phot:sane_check unnatural frequency for photon %d\n", icheck);
			}
		}


		/* Tranport photons through the wind */

		trans_phot(w, p, select_extract);

		if (modes.print_windrad_summary)
			wind_rad_summary(w, windradfile, "a");


		spectrum_create(p, freqmin, freqmax, nangles, select_extract);

		/* Write out the detailed spectrum each cycle so that one can see the statistics build up! */
		renorm = ((double)(geo.pcycles)) / (geo.pcycle + 1.0);

		/* Do an MPI reduce to get the spectra all gathered to the master thread */
#ifdef MPI_ON
		gather_spectra_para(spec_spec_helpers, nspectra);
#endif





#ifdef MPI_ON
		if (rank_global == 0)
		{
#endif
			spectrum_summary(specfile, "w", 0, nspectra - 1, select_spectype, renorm, 0);
#ifdef MPI_ON
		}
#endif
		Log("Completed spectrum cycle %3d :  The elapsed TIME was %f\n", geo.pcycle, timer());



		/* JM1304: moved geo.pcycle++ after xsignal to record cycles correctly. First cycle is cycle 0. */

		xsignal(root, "%-20s Finished %3d of %3d spectrum cycles \n", "OK", geo.pcycle, geo.pcycles);

		geo.pcycle++;			// Increment the spectral cycles

#ifdef MPI_ON
		if (rank_global == 0)
		{
#endif
			wind_save(windsavefile);	// This is only needed to update pcycle
			spec_save(specsavefile);
#ifdef MPI_ON
		}
#endif
		check_time(root);
	}


	/* XXXX -- END CYCLE TO CALCULATE DETAILED SPECTRUM */

	phot_gen_sum(photfile, "a");

	/* 57h - 07jul -- ksl -- Write out the freebound information */

#ifdef MPI_ON
	if (rank_global == 0)
	{
#endif
		fb_save("recomb.save");
#ifdef MPI_ON
	}
#endif


	/* Finally done */
#ifdef MPI_ON
	sprintf(dummy, "End of program, Thread %d only", my_rank);	// added so we make clear these are just errors for thread ngit
																// status 
	error_summary(dummy);		// Summarize the errors that were recorded by the program
	Log("Run py_error.py for full error report.\n");
#else
	error_summary("End of program");	// Summarize the errors that were recorded by the program
#endif


#ifdef MPI_ON
	MPI_Finalize();
	Log_parallel("Thread %d Finalized. All done\n", my_rank);
#endif


	xsignal(root, "%-20s %s\n", "COMPLETE", root);
	Log("Completed entire program.  The elapsed TIME was %f\n", timer());
	return EXIT_SUCCESS;
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	print out some basic help on how to run the program
Arguments:		

Returns:
 
Description:	
		
Notes:

The easiest way to create the message, at least initially, is simply to to type
out what you want to appear on the screen and then as \n\ to all of the lines, including
the ones with nothing in them

History:
	081217	ksl	67c - Added so that ksl could remember all of the options
	09feb	ksl	68b - Added info on -v switch

**************************************************************/

int help()
{
	char *some_help;

	some_help = "\
\n\
This program simulates radiative transfer in a (biconical) CV, YSO, quasar or (spherical) stellar wind \n\
\n\
	Usage:  py [-h] [-r] [-t time_max] xxx  or simply py \n\
\n\
	where xxx is the rootname or full name of a parameter file, e. g. test.pf \n\
\n\
	and the switches have the following meanings \n\
\n\
	-h 	to ge this help message \n\
	-r 	restart a run of the progarm reading the file xxx.windsave \n\
	-e change the maximum number of errors before quit- don't do this unless you understand\
	the consequences! \n\
\n\
	-t time_max	limit the total time to approximately time_max seconds.  Note that the program checks \n\
		for this limit somewhat infrequently, usually at the ends of cycles, because it \n\
		is attempting to save the program outputs so that the program can be restarted with \n\
		-r if that is desired. \n\
\n\
	-v n	controls the amount of print out.  The default is 4.  Larger numbers increase  \n\
		the amount printed; smaller numbers decrease it.   \n\
	if one simply types py or pyZZ where ZZ is the version number one is queried for a name \n\
	of the parameter file. \n\
\n\
\n\
";								// End of string to provide one with help

	printf("%s\n", some_help);

	exit(0);
}

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	init_geo initializes the geo structure to something that is semi-reasonable
Arguments:		

Returns:
 
Description:	
	Initial values for all of the variables that are not part of the individual
	wind descriptions that are actully read(!) into the program should be 
	created here.  The derived values are not needed.

		
Notes:

	When initializing geo, be sure to initialize to cgs units, since this cgs units
	are the working units for the program.  This is necessary for consistncy when
	one tries to restart the program.

	Note that init_geo is set up for CVs and Stars and not AGN


History:
 	98dec	ksl	Coded and debugged.  Much of code was copied from old main routine for
			python
	04dec	ksl	This is probably still not completely up to date, but have
			added some initializations 
	080518	ksl	60a - modified to set all spectypes to SPECTYPE_BB, as part of 
			effort to get restarting models to work better
	080518	ksl	60a - modified all inputs to be in cgs units.  Added a small
      			amount of code to initialize model_list	
	081112	ksl	67 - moved more of initializaiton of geo into
			this routine as part of genearl cleanup of the main
			routine

**************************************************************/

int init_geo()
{
	geo.coord_type = 1;
	geo.ndim = 30;
	geo.mdim = 30;
	geo.disk_z0 = geo.disk_z1 = 0.0;	// 080518 - ksl - moved this up
	geo.adiabatic = 1;			// Default is now set so that adiabatic cooling is included in the wind
	geo.auger_ionization = 1;	// Default is on.


	geo.wind_type = 0;			// Schlossman and Vitello

	geo.star_ion_spectype = geo.star_spectype
		= geo.disk_ion_spectype = geo.disk_spectype = geo.bl_ion_spectype = geo.bl_spectype = SPECTYPE_BB;
	geo.agn_ion_spectype = SPECTYPE_POW;	// 130605 - nsh - moved from python.c

	geo.log_linear = 0;			/* Set intervals to be logarithmic */

	geo.rmax = 1e11;
	geo.rmax_sq = geo.rmax * geo.rmax;
	geo.rstar = 7e8;
	geo.rstar_sq = geo.rstar * geo.rstar;
	geo.mstar = 0.8 * MSOL;
	geo.m_sec = 0.4 * MSOL;
	geo.period = 3.2 * 3600;
	geo.tstar = 40000;
	geo.twind = 40000;
	geo.wind_mdot = 1.e-9 * MSOL / YR;

	geo.ioniz_mode = IONMODE_ML93;	/* default is on the spot and find the best t */
	geo.line_mode = 3;			/* default is escape probabilites */

	geo.star_radiation = 1;		/* 1 implies star will radiate */
	geo.disk_radiation = 1;		/* 1 implies disk will radiate */
	geo.bl_radiation = 0;		/* 1 implies boundary layer will radiate */
	geo.wind_radiation = 0;		/* 1 implies wind will radiate */

	geo.disk_type = 1;			/* 1 implies existence of a disk for purposes of absorption */
	geo.diskrad = 2.4e10;
	geo.disk_mdot = 1.e-8 * MSOL / YR;

	geo.t_bl = 100000.;


	strcpy(geo.atomic_filename, "data/standard77");
	strcpy(geo.fixed_con_file, "none");

	// Note that geo.model_list is initialized through get_spectype 


	return (0);
}

/* 
   Perform some simple checks on the photon distribution just produced.

   History: 01 ksl Removed from main routine 02jul ksl Loosened frequency limits to reflect the fact that in some cases, e.g. those 
   in which the photon distribution has been split into small energy segments, Doppler shifts move photons out of that region.
   08mar ksl Updated slightly, ane eliminated any frequency checks photons generated by macro atoms since these often get out of
   range. 090124 ksl Modified slightly to reduce output if all is OK and if not debugging

 */
int photon_checks(p, freqmin, freqmax, comment)
	 char *comment;
	 PhotPtr p;
	 double freqmin, freqmax;
{
	int nnn, nn;
	// double lum_ioniz; //NSH 16/2/2011 These are now declared externally to allow python to see them
	// int n_ioniz;
	int nlabel;

	geo.n_ioniz = 0;
	geo.lum_ioniz = 0.0;
	nnn = 0;
	nlabel = 0;


	/* Next two lines are to allow for fact that photons generated in a frequency range may be Doppler shifted out of that range,
	   especially if they are disk photons generated right up against one of the frequency limits 04aug--ksl-increased limit from
	   0.02 to 0.03, e.g from 6000 km/s to 9000 km/s 11apr--NSH-decreased freqmin to 0.4, to take account of double redshifted
	   photons. shift. */

	Debug("photon_checks: %s\n", comment);

	freqmax *= (1.8);
	freqmin *= (0.6);
	for (nn = 0; nn < NPHOT; nn++)
	{
		p[nn].np = nn;			/* NSH 13/4/11 This is a line to populate the new internal photon pointer */
		if (H * p[nn].freq > ion[0].ip)
		{
			geo.lum_ioniz += p[nn].w;
			geo.n_ioniz += p[nn].w / (H * p[nn].freq);
		}
		if (sane_check(p[nn].freq) != 0 || sane_check(p[nn].w))
		{
			if (nlabel == 0)
			{
				Error("photon_checks: nphot  origin  freq     freqmin    freqmax\n");
				nlabel++;
			}
			Error
				("photon_checks:sane_check %6d %5d %10.4e %10.4e %10.4e %5d w %10.4e \n",
				 nn, p[nn].origin, p[nn].freq, freqmin, freqmax, p[nn].w);
			p[nn].freq = freqmax;
			nnn++;
		}
		if (p[nn].origin < 10 && (p[nn].freq < freqmin || freqmax < p[nn].freq))
		{
			if (nlabel == 0)
			{
				Error("photon_checks: nphot  origin  freq     freqmin    freqmax\n");
				nlabel++;
			}
			Error
				("photon_checks: %6d %5d %10.4e %10.4e %10.4e freq out of range\n", nn, p[nn].origin, p[nn].freq, freqmin, freqmax);
			p[nn].freq = freqmax;
			nnn++;
		}
		if (nnn > 100)
		{
			Error("photon_checks: Exiting because too many bad photons generated\n");
			exit(0);
		}
	}
	Log("NSH Geo.n_ioniz=%e\n", geo.n_ioniz);

	if (nnn == 0)
		Debug("photon_checks: All photons passed checks successfully\n");

	return (0);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	Generalized routine to get the spectrum type and read the model files 
Arguments:		

Returns:
 
Description:	

		
Notes:

	The routine is slightly dangerous in the sense that if assumes that
	rdint always wants 0 for BB, models for 1, and 2 for uniform.  If
	we were to add another internally generated spectrum type one would
	have to carefull consider how to do this.  

	For models we have to handle two cases:
		A new model.  Here we want to start with a default value and
			to keep track of what was entered since it is likely
			we will want that together
		The continuation of an old model.  Here we need to expect the
			same choices as previously
081026 - Actual routine is still a mess.  
		* yesno is a statment of whether the component exists
		* question is really tightly associted with the way the program is 
		written
		* The progam returns the component type two ways, on through the call
		* and the other through the return.  It ssems like one would do.:




History:
	080518	ksl	Coded as part of effort to make models restart more
			easily, but also simplifies some of the code
        121025  nsh	added a mode for power law

**************************************************************/


char get_spectype_oldname[LINELENGTH] = "data/kurucz91.ls";	/* This is to assure that we read model lists in the same order
															   everytime */
int get_spectype_count = 0;
int get_spectype(yesno, question, spectype)
	 int yesno;
	 char *question;
	 int *spectype;
{
	char model_list[LINELENGTH];
	int stype;
	int get_models();			// Note: Needed because get_models cannot be included in templates.h
	if (yesno)
	{
		// First convert the spectype to the way the questionis supposed to be answered
		if (*spectype == SPECTYPE_BB || *spectype == SPECTYPE_NONE)
			stype = 0;
		else if (*spectype == SPECTYPE_UNIFORM)
			stype = 2;
		else if (*spectype == SPECTYPE_POW)
			stype = 3;
		else
			stype = 1;
		/* Now get the response */
		rdint(question, &stype);
		/* Now convert the response back to the values which python uses */
		if (stype == 0)
			*spectype = SPECTYPE_BB;	// bb
		else if (stype == 2)
			*spectype = SPECTYPE_UNIFORM;	// uniform
		else if (stype == 3)
			*spectype = SPECTYPE_POW;	// power law
		else if (stype == 4)
			*spectype = SPECTYPE_CL_TAB;
		else
		{
			if (geo.wind_type == 2)
			{					// Continuing an old model
				strcpy(model_list, geo.model_list[get_spectype_count]);
			}
			else
			{					// Starting a new model
				strcpy(model_list, get_spectype_oldname);
			}
			rdstr("Model_file", model_list);
			get_models(model_list, 2, spectype);
			strcpy(geo.model_list[get_spectype_count], model_list);	// Copy it to geo 
			strcpy(get_spectype_oldname, model_list);	// Also copy it back to the old name
			get_spectype_count++;
		}
	}
	else
	{
		*spectype = SPECTYPE_NONE;	// No radiation
	}

	return (*spectype);
}


/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:
	The next couple of routines are for recording information about photons/energy impinging
	on the disk, which is stored in a disk structure called qdisk.

	qdisk_init() just initializes the structure (once the disk structue has been initialized.

	qdisk_save records the results in a file

Arguments:		

Returns:
 
Description:	

		
Notes:



History:
	04mar	ksl	add section to dump heating of disk.  
 			Factor of 2 in area calculation reflects
 			fact that disk has two sides
	04dec	ksl	created variable ztot so fractional heating
			is correct if multiple subcyles
	080519	ksl	60a - Added code to calculate the irradiation of the disk
			in terms of t and w.  This should help to monitor the effect
			of irradiation on the disk

**************************************************************/


int qdisk_init()
{
	int n;
	for (n = 0; n < NRINGS; n++)
	{
		qdisk.r[n] = disk.r[n];
		qdisk.t[n] = disk.t[n];
		qdisk.g[n] = disk.g[n];
		qdisk.v[n] = disk.v[n];
		qdisk.heat[n] = 0.0;
		qdisk.nphot[n] = 0;
		qdisk.nhit[n] = 0;
		qdisk.w[n] = 0;
		qdisk.ave_freq[n] = 0;
		qdisk.t_hit[0] = 0;
	}
	return (0);
}

int qdisk_save(diskfile, ztot)
	 char *diskfile;
	 double ztot;
{
	FILE *qptr;
	int n;
	double area, theat;
	qptr = fopen(diskfile, "w");
	fprintf(qptr, "# r       zdisk     t_disk     heat      nhit nhit/nemit  t_heat    t_irrad  W_irrad\n");
	for (n = 0; n < NRINGS; n++)
	{
		area = (2. * PI * (qdisk.r[n + 1] * qdisk.r[n + 1] - qdisk.r[n] * qdisk.r[n]));
		theat = qdisk.heat[n] / area;
		theat = pow(theat / STEFAN_BOLTZMANN, 0.25);	// theat is temperature if no internal energy production
		if (qdisk.nhit[n] > 0)
		{

			qdisk.ave_freq[n] /= qdisk.heat[n];
			qdisk.t_hit[n] = H * qdisk.ave_freq[n] / (BOLTZMANN * 3.832);	// Basic conversion from freq to T
			qdisk.w[n] =
				qdisk.heat[n] / (4. * PI * STEFAN_BOLTZMANN * area *
								 qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n]);
		}

		fprintf(qptr,
				"%8.3e %8.3e %8.3e %8.3e %5d %8.3e %8.3e %8.3e %8.3e\n",
				qdisk.r[n], zdisk(qdisk.r[n]), qdisk.t[n],
				qdisk.heat[n], qdisk.nhit[n], qdisk.heat[n] * NRINGS / ztot, theat, qdisk.t_hit[n], qdisk.w[n]);
	}

	fclose(qptr);
	return (0);
}



/***********************************************************
	Space Telescope Science Institute

Synopsis:
	Stuart's routine to read a non-standard disk profile
	for YSO effort

Arguments:		

Returns:
 
Description:	

		
Notes:
	Originally part of main routine; moved to separate routine
	by ksl sometime in the fall of 08



History:

**************************************************************/

int read_non_standard_disk_profile(tprofile)
	 char *tprofile;
{

	FILE *fopen(), *fptr;
	int n;
	float dumflt1, dumflt2;
	int dumint;

	if ((fptr = fopen(tprofile, "r")) == NULL)
	{
		Error("Could not open filename %s\n", tprofile);
		exit(0);
	}

	fscanf(fptr, "%d\n", &dumint);
	blmod.n_blpts = dumint;
	for (n = 0; n < blmod.n_blpts; n++)
	{
		fscanf(fptr, "%g %g", &dumflt1, &dumflt2);
		blmod.r[n] = dumflt1 * 1.e11;
		blmod.t[n] = dumflt2 * 1.e3;
	}

	fclose(fptr);

	return (0);
}


/***********************************************************
				University of Southampton

Synopsis:
	init_advanced_modes simply initialises the set of 
	advanced modes stored in the modes structure to a 
	default value. For now, this is 0 (off).

Arguments:	
    none	

Returns:
    just initialises modes which is declared in python.h
 
Description:	
	
Notes:
    see #111 and #120

History:
    1410 -- JM -- Coded
**************************************************************/


int init_advanced_modes()
{
	modes.iadvanced = 0;		// this is controlled by the -d flag, global mode control.
	modes.save_cell_stats = 0;	// want to save photons statistics by cell
	modes.ispy = 0;				// want to use the ispy function
	modes.keep_ioncycle_windsaves = 0;	// want to save wind file each ionization cycle
	modes.track_resonant_scatters = 0;	// want to track resonant scatters
	modes.save_extract_photons = 0;	// we want to save details on extracted photons
	modes.print_windrad_summary = 0;	// we want to print the wind rad summary each cycle
	modes.adjust_grid = 0;		// the user wants to adjust the grid scale
	modes.diag_on_off = 0;		// extra diagnostics
	modes.use_debug = 0;
	modes.print_dvds_info = 0;	// print out information on velocity gradients
	write_atomicdata = 0;		// print out summary of atomic data 
	// note this is defined in atomic.h, rather than the modes structure 

	return (0);
}
