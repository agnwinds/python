#ifdef MPI_ON
#include "mpi.h"


#endif

#define UV_low 7.4e14           //The lower frequency bound of the UV band as defined in IOS 21348
#define UV_hi 3e16              //The lower frequency bound of the UV band as defined in IOS 21348

int q_test_count;

int np_mpi_global;              /// Global variable which holds the number of MPI processes

int rank_global;


int verbosity;                  /* verbosity level. 0 low, 10 is high */

/* the functions contained in log., rdpar.c and lineio.c are
   declare deparately from templates. This is because some functions
   only use log.h and don't use python.h due to repeated definitions */
#include "log.h"
#include "strict.h"

/* In python_43 the assignment of the WindPtr size has been moved from a fixed
value determined by values in python.h to a values which are adjustable from
within python */


#define REL_MODE_LINEAR 0      /*Only make v/c corrections when doing frame transfers*/
#define REL_MODE_FULL   1      /*Make full corrections for special relativity*/

  int rel_mode;                 /* How doppler effects are treated */



/* With domains NDIM and MDIM need to be removed but NDIM2 is the total number of cells in wmain, and there
are certain times we want to loop over everything.  The situation with NPLASMA is similar */

int NDIM2;                      //The total number of wind cells in wmain
int NPLASMA;                    //The number of cells with non-zero volume or the size of plasma structure

char basename[132];             // The root of the parameter file name being used by python

/* These are tunable parameters that control various aspects of python
 * and the assorted programs.  In some cases they affect the "care" with
 * which a calculation is made, and this can in principle affect the
 * speed of the program.  One day if may be desirable to change some
 * of these parameters with time in the program.  At present they are
 * simply collected here
 * 
 * */

double DFUDGE;
#define XFUDGE   1e-5           // The scale factor used in setting up cell x cell dfudge

#define VCHECK	1.e6            // The maximum allowable error in calculation of the velocity in calculate_ds


double SMAX_FRAC;               /* In translate_in_wind, a limit is placed on the maximum distance a
                                   photon can travel in one step.  It is a fraction SMAX_FRAC of the
                                   distance of the photon from the origin.  This had been hardwired to
                                   0.1 for up to 57h.  Changing it to 0.5 speeds up the current version
                                   of the code by as much as a factor of 2 for small sized grids.  This
                                   had been introduced as part of the attempt to assure ourselves that
                                   line shapes were calculated as accurately as possible.  The underlhying
                                   rational for having a maximum disstance is associated with the fact that
                                   we use linear interpolation along the line of sight to establish velocities
                                   and most of our grid cells in 2.5d are actually hoop shaped, which means
                                   one can travel a long distance within a hoop if the direction of the photon
                                   is not more or less radial, but if moving along the hoop. 
                                 */
double DENSITY_PHOT_MIN;        /* This constant is a minimum density for the purpose of calculating
                                   photoionization heating and recombination cooling.  It is important that heating and cooling
                                   be calculated self-consistently.  Program speed is somewhat sensitive 
                                   to this parameter, at the 10% level if raised from 1e-3 to 1.  There is a 
                                   trade-off since lower minima may give better results, especially for macro atoms. */

#define LDEN_MIN        1e-3    /* The minimum density required for a line to be conidered for scattering
                                   or emission in calculate_ds and lum_lines */


/* End of "care factor" definition */


#define RMIN   				1.e9
#define RMAX   				3.e10
#define VWIND  				2.e8
#define MDOT  				1.e-9*MSOL/YR
#define BETA  				1.0
#define KAPPA_CONT 			4.
#define EPSILON  			1.e-6   /* A general purpose fairly small number */
#define NSTAT 				10      // JM increased this to ten to allow for adiabatic
#define VMAX                		1.e9
#define TAU_MAX				20.     /* Sets an upper limit in extract on when
                                                   a photon can be assumed to be completely absorbed */

#define DELTA_V				1.      /*This is the accuracy in velocity space (cm/s) that we sample edges when producing freebound photons */

#define DANG_LIVE_OR_DIE   2.0  /* If constructing photons from a live or die run of the code, the
                                   angle over which photons will be accepted must be defined */

double PHOT_RANGE;              /* When a variable number of photons are called in different ionization
                                   cycles this is the log of the difference between NPHOT_MAX
                                   and the value in the first cycle
                                 */
int NPHOT_MAX;                  /* The maximum number of photon bundles created per cycle */
int NPHOT;                      /* The number of photon bundles created, defined in setup.c */

#define NWAVE  			  10000 //This is the number of wavelength bins in spectra that are produced
#define MAXSCAT 			2000

/* Define the structures */


/* Geometry is an actual structure which exists at the beginning of th program
   It carries the variables which define the geometry.  Reasonable values of each of
   these should be defined before it is altered with inputs fromt eh terminal. 
   The geometry structure is used to transfer all of the information about a wind


 */

/* Definitions of spectral types, which are all negative because when
 * one reads a spectrum from a list of models these are numbered beginning
 * with zero, see the discussion in get_models.c   080518 - ksl - 60a
 */
#define SPECTYPE_BB      -1
#define SPECTYPE_UNIFORM -2
#define SPECTYPE_POW     -4
#define SPECTYPE_CL_TAB  -5
#define SPECTYPE_BREM    -6
#define SPECTYPE_NONE	 -3
#define SPECTYPE_MODEL	 -99    // This is just used briefly, before a model number is assigned

/* Number of model_lists that one can have, should be the same as NCOMPS in models.h */
#define NCOMPS 	10
#define LINELENGTH 	256

/* This structure contains the information needed for each separate region of space, e.g the
 * wind and the disk
 */

// This is intialized in init_geo, but it my need to be in geo in order to be able to read
// everything back

enum coord_type_enum
{ SPHERICAL = 0,
  CYLIND = 1,
  RTHETA = 2,
  CYLVAR = 3
};


/* List of possible wind_types */

#define SV   			0
#define	STAR    		1
/* PREVIOUS is no longer an allowed type. Reading in an early model is now
 * handled as a system_type 
 */
// #define      PREVIOUS                2
#define	HYDRO 			3
#define	CORONA 			4
#define KNIGGE			5
#define	HOMOLOGOUS 		6
//OLD #define   YSO                     7
#define	SHELL 			9
#define IMPORT          10      // Model that is read in from a file
#define	DISK_ATMOS 		11


#define MaxDom			10

/* Next define structures that pertain to possilbe region geometries
 
   These definitions had to be moved up in python.h because they need to be defined 
   prior to defining the domains, which must contain these structures in the new
   schema  ksl 15aug
*/

typedef struct plane            /*SWM 10-10-14 - Switched to TypeDef */
{
  double x[3];                  /* A position included in the plane (usually the "center" */
  double lmn[3];                /* A unit vector perpendicular to the plane (usually in the "positive" direction */
} plane_dummy, *PlanePtr;
plane_dummy plane_l1, plane_sec, plane_m2_far;  /* these all define planes which are perpendicular to the line of sight from the 
                                                   primary to the seconday */


/* Note that since we are interested in biconical flows, our definition of a cone is not exactly
 * what one might guess.  The cone is defined in the positive z direction but reflected through 
 * the xy plane.  
 * 56d -- Beginning with 56d, ksl has switched to a new definition of cones, that is intended to
 * make it possible to use ds_to_cone easier as part of different coordinate systems.  The new definition
 * is based on the intersection of the cone with the z axis rather than the intersection with
 * the disk plane.  At present both definitions are used in the program and therefore both shold
 * be defined.  Once the new definition is promulgated through the entire program, and verified
 * the old definitions can be elimiated.  05jul -- ksl
 */

typedef struct cone
{
  double z;                     /* The place where the cone intersects the z axis (used after 56d) */
  double dzdr;                  /* the slope (used after 56d) */
}
cone_dummy, *ConePtr;


/* End of structures which are used to define boundaries to the emission regions */

#define NDIM_MAX 1000           // maximum size of the grid in each dimension

typedef struct domain
{
  char name[LINELENGTH];
  int wind_type;
  int ndim, mdim, ndim2;        //ndim is the size in the x direction, while mdim is the size in z or theta direction
  int nstart, nstop;            //the beginning and end (-1) location in wmain of this component
  enum coord_type_enum coord_type;
  int log_linear;               /*0 -> the grid spacing will be logarithmic in x and z, 1-> linear */
  double xlog_scale, zlog_scale;        /* Scale factors for setting up a logarithmic grid, the [1,1] cell
                                           will be located at xlog_scale,zlog_scale */

  /* The next few structures define the boundaries of an emission region */
  struct cone windcone[2];      /* The cones that define the boundary of winds like SV or kwd */
  struct plane windplane[2];    /* Planes which define the top and bottom of a layer */
  double rho_min, rho_max;      /* These are used for the inner and outer boundary of a pillbox */

  double wind_x[NDIM_MAX], wind_z[NDIM_MAX];    /* These define the edges of the cells in the x and z directions */
  double wind_midx[NDIM_MAX], wind_midz[NDIM_MAX];      /* These define the midpoints of the cells in the x and z directions */

  ConePtr cones_rtheta;         /*A ptr to the cones that define boundaries of cells in the theta direction 
                                   when rtheta coords  are being used */
/* Next two lines are for cyl_var coordinates.  They are used in locating the appropriate 
 * locating the appropriate cell, for example by cylvar_where_in_grid
 */

  double wind_z_var[NDIM_MAX][NDIM_MAX];
  double wind_midz_var[NDIM_MAX][NDIM_MAX];


/* Since in principle we can mix and match arbitrarily the next parameters now have to be part of the domain structure */

  /* Generic parameters for the wind */
  double wind_mdot, stellar_wind_mdot;  /* Mass loss rate in disk and stellar wind */
  double rmin, rmax;            /*Spherical extent of the wind */
  double zmin, zmax;            /* Vertical extent of the wind, often the same as rmax */
  double wind_rho_min, wind_rho_max;    /*Min/Max rho for wind in disk plane */
  double wind_thetamin, wind_thetamax;  /*Angles defining inner and outer cones of wind, measured from disk plane */
  double mdot_norm;             /*A normalization factor used in SV wind, and Knigge wind */

  double twind;                 // Initial temperature for a domain

  /* Parameters defining Shlossman & Vitello Wind */
  double sv_lambda;             /* power law exponent describing from  what portion of disk wind is radiated */
  double sv_rmin, sv_rmax, sv_thetamin, sv_thetamax, sv_gamma;  /* parameters defining the goemetry of the wind */
  double sv_v_zero;             /* velocity at base of wind */
  int sv_v_zero_mode;           /* use fixed initial velocity or multiple of sound speed */
#define FIXED 0
#define SOUND_SPEED 1
  double sv_r_scale, sv_alpha;  /* the scale length and power law exponent for the velocity law */
  double sv_v_infinity;         /* the factor by which the velocity at infinity exceeds the excape velocity */


  /* Parameters defining Knigge Wind */
  double kn_dratio;             /* parameter describing collimation of wind */
  double kn_lambda;             /* power law exponent describing from  what portion of disk wind is radiated */
  double kn_r_scale, kn_alpha;  /* the scale length and power law exponent for the velocity law */
  double kn_v_infinity;         /* the factor by which the velocity at infinity exceeds the excape velocity */
  double kn_v_zero;             /* NSH 19/04/11 - Added in as the multiple of the sound speed to use as the initial velocity */

  /* Parameters describing Castor and Larmors spherical wind */
  double cl_v_zero, cl_v_infinity, cl_beta;     /* Power law exponent */
  double cl_rmin, cl_rmax;

  /* Parameters describing a spherical shell test wind */
  double shell_vmin, shell_vmax, shell_beta;
  double shell_rmin, shell_rmax;

  /*Parameters defining a corona in a ring above a disk */
  double corona_rmin, corona_rmax;      /*the minimum and maximu radius of the corona */
  double corona_zmax;           /*The maximum vertical extent of the corona */
  double corona_base_density, corona_scale_height;      /*the density at the base of the corona and the scale height */
  double corona_vel_frac;       /* the radial velocity of the corona in units of the keplerian velocity */

  /* The filling factior for the wind or corona */
  /* JM 1601 -- Moved here from geo, see #212 */
  double fill;
}
domain_dummy, *DomainPtr;       // One structure for each domain

DomainPtr zdom;                 //This is the array pointer that contains the domains
int current_domain;             // This integer is used by py_wind only


/* the geometry structure contains information that applies to all domains, including
 * the basic system geometry, descriptions of the radition sources, and truly 
 * global information including how ionization calculations are caried out. 
 *
 * Information that is domain specific should be placed directly in the domain
 * structure.  ksl
 */

#define SYSTEM_TYPE_STAR   0
#define SYSTEM_TYPE_CV     1
#define SYSTEM_TYPE_BH     4
#define SYSTEM_TYPE_AGN    2
#define	SYSTEM_TYPE_PREVIOUS   	   3

/* RUN_TYPE differs from SYSTEM_TYPE in that
  it has implications on how the program is run
  wherease SYSTEM_TYPE refers (mainly) to the type
  of sytem, with the exception 
*/

#define RUN_TYPE_NEW       0
#define RUN_TYPE_RESTART   1
#define RUN_TYPE_PREVIOUS  3

#define TRUE  1
#define FALSE 0



struct geometry
{
  int system_type;              /* See allowed types above. system_type should only be used for setp */
  int binary;                   /* Indicates whether or not the system is a binary. TRUE or FALSE */

  int ndomain;                  /* The number of domains in a model */
  int ndim2;                    /* The total number of windcells in all domains */
  int nplasma, nmacro;          /* The total number of cells in the plasma and macro structures in all domains */

  /* variables which store the domain numbers of the wind, disk atmosphere.
     Other components should be added here.  Right now we need a wind_domain 
     number because the inputs for the disk and a putativel disk atmosphere are 
     interrsed.  The first step will be to put this information into alocal variale
     in python.c. We should not have to carry this forward */

  int wind_domain_number;
  /* Ultimately the next variable should not be needed but to avoid a bad
   * interaction with Nick's effort meld zeus calculations with Python
   * I have added a new variable.  It is likely that the domain number for
   * this will always be 0, but one could imagine cases where that might
   * not be the case  ksl -160927
   */

  int hydro_domain_number;      // Created for the special case of runs with Zeus


  /* This section allows for restarting the program, and adds parameters used
   * in the calculation */

  int wcycle, pcycle;           /* The number of completed ionization and spectrum cycles */
  int wcycles, pcycles, pcycles_renorm; /* The number of ionization and spectrum cycles desired, pcycles_renorm 
                                         * is only used on restarts.  See spectrum_restart_renormalize
                                         */


  /* This section stores information which specifies the spectra to be extracted.  Some of the parameters
   * are used only in advanced modes.  
   */

#define NSPEC   20
  int nangles;
  double angle[NSPEC], phase[NSPEC];
  int scat_select[NSPEC], top_bot_select[NSPEC];
  double rho_select[NSPEC], z_select[NSPEC], az_select[NSPEC], r_select[NSPEC];
  double swavemin, swavemax, sfmin, sfmax;      // The minimum and maximum wavelengths/freqs for detailed spectra
  int select_extract, select_spectype;

/* Begin description of the actual geometry */

/* The next variables refere to the entire space in which pbotons sill be tracked.  Photons
 * outside these regions are assumed to have hit something or be freely moving through space.
 */

  double rmax, rmax_sq;         /* The maximum distance to which a photon should be followed */

/* Basic paremeters of the system, as opposed to elements of the wind or winds */

  double mstar, rstar, rstar_sq, tstar, gstar;  /* Basic parameters for the star (often a WD) in the system */
  double tstar_init;            /* The temperature of the star, before backscattering is taken into account */
  double lum_star_init, lum_star_back;  /* The luminosity of the star as determined by tstar_init */

  double tmax;                  /*NSH 120817 the maximum temperature of any element of the model 
                                   - used to help estimate things for an exponential representation of the spectrum in a cell */


#define DISK_NONE   0
#define DISK_FLAT   1
#define DISK_VERTICALLY_EXTENDED   2

  int disk_type;

#define BACK_RAD_ABSORB_AND_DESTROY  0  /* Disk simply absorbs the radiation and it is lost */
#define BACK_RAD_SCATTER            1   /* Disk reradiates the radiation immediately via electron scattering */
#define BACK_RAD_ABSORB_AND_HEAT     2  /* Correct disk temperature for illumination by photons 
                                           which hit the dsik.  Disk radiation is absorbed and changes 
                                           the temperature of the disk for future ionization cycles
                                         */

  int absorb_reflect;           /*Controls what happens when a photon hits the disk or star
                                 */

#define DISK_TPROFILE_STANDARD          0       // This is a standard Shakura-Sunyaev disk. The profile depends on mstar and mdot_disk
#define DISK_TPROFILE_READIN            1       // Here the temperature profile for the disk is simply read in as a function of radius
//OLD #define DISK_TPROFILE_YSO               2       // The so-called YSO option was created for the YSO case
  int disk_tprofile;            /* This is an variable used to specify a standard accretion disk (0) or
                                   one that has been read in and stored. */
  double disk_mdot;             /* mdot of  DISK */
  double diskrad, diskrad_sq;
  double disk_z0, disk_z1;      /* For vertically extended disk, z=disk_z0*(r/diskrad)**disk_z1 *diskrad */
  double lum_disk_init, lum_disk_back;  /* The intrinsic luminosity of the disk, the back scattered luminosity */
  int run_type;                 /*1508 - New variable that describes whether this is a continuation of a previous run 
                                   Added in order to separate the question of whether we are continuing an old run fro
                                   the type of wind model.  Bascially if run_type is 0, this is a run from scratch,
                                   if SYSTEM_TYPE_PREVIOUS it is an old run     */
  int star_radiation, disk_radiation;   /* 1 means consider radiation from star, disk,  bl, and/or wind */
  int bl_radiation, wind_radiation, agn_radiation;
  int search_light_radiation;   /* 1605 - ksl - Added to implement 1d testing */
  int matom_radiation;          /* Added by SS Jun 2004: for use in macro atom computations of detailed spectra
                                   - 1 means use emissivities for BOTH macro atom levels and kpkts. 0 means don't
                                   (which is correct for the ionization cycles. */
  int ioniz_mode;               /* describes the type of ionization calculation which will
                                   be carried out.  The various ioniz_modes are defined by #defines IONMODE_MATRIX_BB
                                   etc.  See the documentation in this file for what these mean. */
  int macro_ioniz_mode;         /* Added by SS Apr04 to control the use of macro atom populations and
                                   ionization fractions. If it is set to 1 then macro atom populations
                                   computed from estimators are used. If set to 0 then the macro atom
                                   populations are computed as for minor ions. By default it is set to
                                   0 initially and then set to 1 the first time that
                                   Monte Carlo estimators are normalised. */
  int ioniz_or_extract;         /* Set to 1 (true) during ionization cycles, set to 0 (false) during calculation of
                                   detailed spectrum.  Originally introduced by SS in July04 as he added
                                   macro atoms.  Name changed by ksl (57h) since this variable can be used more
                                   generally to speed up the extract portion of the calculation.
                                 */
  int macro_simple;             /* Added by SS May04 for diagnostics. As default this is set to 0. A full
                                   Macro Atom calculation is performed in that case. If it is set to 1 it means
                                   that although Macro Atom data has been read in, all lines/continua are treated
                                   using the simplified two-level approximation. Such a calculation should reproduce
                                   the same results as pre-Macro Atom versions of the code. */
  int partition_mode;           /* Diagnostic to force the partition function to be calculated in
                                   a specific way. */

#define LINE_MODE_ABSORB      0
#define LINE_MODE_SCAT        1
#define LINE_MODE_SINGLE_SCAT 2
#define LINE_MODE_ESC_PROB    3



  int line_mode;                /*0, the atomosphere is a completely absorbing and no photons
                                   will be scattered.  In this mode, assuming the wind is a source
                                   of emission, the emissivity will be the Einstein A coefficient
                                   1, the atmosphere is completely scattering, there will be no
                                   interchange of energy between the photons and the electrons
                                   as a result of radiation transfer
                                   2, then a simple single scattering approximation is applied in which
                                   case the scattered flux is just  A21/(C21+A21). 
                                   3, then radiation trapping aka escape probabilities are included as well.
                                   6 - 9,  If set to 6-9 initially, this switches on the macro atom case
                                   and then behaves like LINE_MODE_ESC_PROB. 
                                 */
/* Note that the scatter_mode is actually a subsidiary variable of the line_mode.  Choosing a line_mode
 * results in the selection of a scatter_mode */

#define SCATTER_MODE_ISOTROPIC    0
#define SCATTER_MODE_THERMAL      2

  int scatter_mode;             /*The way in which scattering for resonance lines is treated 
                                   0  isotropic
                                   2  thermally broadened anisotropic
                                 */

#define RT_MODE_2LEVEL  1
#define RT_MODE_MACRO   2

  int rt_mode;                  /* radiative transfer mode. 2 for Macro Atom method,  1 for non-Macro Atom methods  */


  /* Define the choices for calculating the FB, see, e.g. integ_fb */

#define FB_FULL         0       /* Calculate fb emissivity including energy associated with the threshold */
#define FB_REDUCED      1       /* Calculqate the fb emissivity without the threshold energy */
#define FB_RATE         2       /* Calulate the fb recombinarion rate  */


  /* Define the modes for free bound integrals */
#define OUTER_SHELL  1
#define INNER_SHELL  2

  /* The frequency bands used when calculating parameters like a power law slope in limited regions. */

#define  NXBANDS 20             /* the maximum number of bands (frequency intervals that can be defined for
                                   storing coarse spectra for each plasma cell */

  int nxfreq;                   /* the number of frequency intervals actually used */
  double xfreq[NXBANDS + 1];    /* the frequency boundaries for the coarse spectra  */


  /* The next set pf variables assign a SPECTYPE (see above) for
     each possible source of radiation in a model.  The value assigned can be different for
     the ionization and detailed spectrum generation part of the code */

  int star_ion_spectype, star_spectype; /* The type of spectrum used to create the continuum
                                           for the star in the ionization and final spectrum calculation */
  int disk_ion_spectype, disk_spectype; /* Same as above but for the disk */
  int bl_ion_spectype, bl_spectype;     /* Same as above but for the boundary layer */
  int agn_ion_spectype, agn_spectype;   /* Same as above but for the AGN */
  int search_light_ion_spectype, search_light_spectype; /* Same as above but for the search_light. Created for 1d test */

  char model_list[NCOMPS][LINELENGTH];  /* The file which contains the model names and the associated values for the model */
  int model_count;              /*The number of distinct models that have been read in */

  double mdot_norm;             /*A normalization factor used in SV wind, and Knigge wind */
  int adiabatic;                /*0-> Do not include adiabatic heating in calculating the cooling of the wind
                                   1-> Use adiabatic heating in calculating the cooling of the wind
                                 */
  int nonthermal;               /* 0 --> No extra heating due to shocks
                                   1 --> Extra heating due to shocks (etc)  (Added for FU Ori)
                                 */

  double shock_factor;          /* A scaling factor used for including an extra heating term (for FU Ori stars
                                 */
  double frac_extra_kpkts;      /* in the case that we have extra heating and macro-atoms, the fraction of 
                                   photons to reserve for those generated directly by k-packets */

  int auger_ionization;         /*0 -> Do not include innershell photoionization /Auger effects; 1-> include them */

/* Initial values for defining wind structure for a planar geometry.  These are currently only used by balance and this
   may not be the best approach generally and depending on where this ends up. Some consolidation is desirable */
  double pl_vol, pl_vmax;
  double pl_t_r, pl_t_e, pl_w;
  double pl_nh;

  /* Variables having to do with heating and cooling */

  double lum_tot, lum_star, lum_disk, lum_bl, lum_wind; /* The total luminosities of the disk, star, bl, & wind 
                                                           are actually not used in a fundamental way in the program */
  double lum_agn;               /*The total luminosity of the AGN or point source at the center */
  double lum_ff, lum_rr, lum_lines;     /* The luminosity of the wind as a result of ff, fb, and line radiation */
  double cool_rr;               /*1706 NSH - the cooling rate due to radiative recombination - not the same as the luminosity */
  double cool_comp;             /*1108 NSH The luminosity of the wind as a result of compton cooling */
  double cool_di;               /* 1409 NSH The direct ionization luminosity */
  double cool_dr;               /*1109 NSH The luminosity of the wind due to dielectronic recombination */
  double cool_adiabatic;        /*1209 NSH The cooling of the wind due to adiabatic expansion */
  double heat_adiabatic;        /*1307 NSH The heating of the wind due to adiabatic heating - split out from cool_adiabatic to get an accurate idea of whether it is important */
  double heat_shock;            /*1806 - ksl - The amount of extra heating going into the wind due to shock heating. Added for FU Ori project */

  double f1, f2;                /* The freguency minimum and maximum for which the band limted luminosities have been calculated */
  double f_tot, f_star, f_disk, f_bl, f_agn, f_wind;    /* The integrated specific L between a freq min and max which are
                                                           used to establish the band limited luminosity  of photons of various types.
                                                           For detailed spectra cycles, this will the band limed luminosity between
                                                           the minimum and maximum wavelength of the detailed spectrum */

/* These variables are copies of the lum variables above, and are only calculated during ionization cycles
   This is a bugfix for JM130621, windsave bug */
  double lum_ff_ioniz, lum_rr_ioniz, lum_lines_ioniz;
  double cool_rr_ioniz;
  double cool_comp_ioniz;
  double cool_di_ioniz;         /* 1409 NSH The direct ionization luminosity */
  double cool_dr_ioniz;
  double cool_adiabatic_ioniz;
  double lum_wind_ioniz, lum_star_ioniz, lum_disk_ioniz, lum_bl_ioniz, lum_tot_ioniz;

  double f_matom, f_kpkt;       /*Added by SS Jun 2004 - to be used in computations of detailed spectra - the
                                   energy emitted in the band via k-packets and macro atoms respectively. */

//70i - nsh 111007 - put cool_tot_ioniz and n_ioniz into the geo structure. This will allow a simple estimate of ionisation parameter to be computed;

  double n_ioniz, cool_tot_ioniz;

/* The next set of parameters relate to the secondary
 */

  double m_sec, q;              /* Mass of the secondary, mass ratio of system */
  double period;                /* Period of the systems in seconds */
  double a, l1, l2, phi;        /* Separation of primary and secondary, distance of l1 from primary,phi at l1 */
  double l1_from_m2, r2_far;    /* Distance to l1 from m2, distance to far side of secondary from primary */
  double r2_width;              /* Maximum width of Roche surface of secondary in the plane of orbit */

  double t_bl;                  /*temperature of the boundary layer */
  double weight;                /*weight factor for photons/defined in define_phot */

/* The next set of parameters relate to the central source of an AGN
 */

  double brem_temp;             /*The temperature of a bremsstrahlung source */
  double brem_alpha;            /*The exponent of the nu term for a bremstrahlung source */

  double pl_low_cutoff;         /* accessible only in advanced mode- see #34. default to zero */

  double alpha_agn;             /*The power law index of a BH at the center of an AGN.  Note that the luminosity
                                   of the agn is elsewhere in the structure
                                 */
  double const_agn;             /*The constant for the Power law, there are lots of ways of defining the PL which is best? */
  double d_agn;                 /* the distance to the agn - only used in balance to calculate the ionization fraction */


  int pl_geometry;              /* geometry of X-ray point source */
#define PL_GEOMETRY_SPHERE 0
#define PL_GEOMETRY_LAMP_POST 1
  double lamp_post_height;      /* height of X-ray point source if lamp post */

/* The next four variables added by nsh Apr 2012 to allow broken power law to match the cloudy table command */
  double agn_cltab_low;         //break at which the low frequency power law ends
  double agn_cltab_hi;          //break at which the high frequency power law cuts in
  double agn_cltab_low_alpha;   //photon index for the low frequency end
  double agn_cltab_hi_alpha;    //photon index for the high frequency end       

// The next set of parameters describe the input datafiles that are read
  char atomic_filename[132];    /* The masterfile for the atomic data */
  char fixed_con_file[132];     /* For fixed concentrations, the file specifying concentrations */

  //Added by SWM for tracking C-IV/H-A hotspots
  int nres_halpha;

  /* Variables used for reverberation mapping */

  double fraction_converged, reverb_fraction_converged;
  int reverb_filter_lines, *reverb_filter_line;
  enum reverb_disk_enum
  { REV_DISK_CORRELATED = 0, REV_DISK_UNCORRELATED = 1, REV_DISK_IGNORE = 3 } reverb_disk;
  enum reverb_enum
  { REV_NONE = 0, REV_PHOTON = 1, REV_WIND = 2, REV_MATOM = 3 } reverb;
  enum reverb_vis_enum
  { REV_VIS_NONE = 0, REV_VIS_VTK = 1, REV_VIS_DUMP = 2, REV_VIS_BOTH = 3 } reverb_vis;
  int reverb_wind_cycles;
  int reverb_path_bins, reverb_angle_bins;      //SWM - Number of bins for path arrays, vtk output angular bins
  int reverb_dump_cells;        //SWM - Number of cells to dump
  double *reverb_dump_cell_x, *reverb_dump_cell_z;
  int *reverb_dump_cell;
  int reverb_lines, *reverb_line;       //SWM - Number of lines to track, and array of line 'nres' values

  int spec_mod;                 //A flag to say that we do hav spectral models  ??? What does this mean???
}
geo;

/*
 * EP: 27/09/19
 * Added enumerator to define different banding modes to make the banding
 * code more self-documenting
 */

enum band_definition_enum
{
  T_STAR_BAND = 0,
  MIN_MAX_FREQ_BAND = 1,
  CV_BAND = 2,
  YSO_BAND = 3,
  USER_DEF_BAND = 4,
  CLOUDY_TEST_BAND = 5,
  WIDE_BAND = 6,
  AGN_BAND = 7,
  LOG_USER_DEF_BAND = 8
};

/* xdisk is a structure that is used to store information about the disk in a system */
#define NRINGS	3001            /* The actual number of rings completely defined
                                   is NRINGS-1 ... or from 0 to NRINGS-2.  This is
                                   because you need an outer radius...but the rest
                                   of this element is not filled in. */

struct xdisk
{
  double r[NRINGS];             /* The inner radius of an annulus */
  double t[NRINGS];             /* The temperature at the middle of the annulus */
  double g[NRINGS];             /* The gravity at the middle of the annulus */
  double v[NRINGS];             /* The velocity at the middle of the annulus */
  double heat[NRINGS];          /* The total energy flux of photons hitting each annulus */
  double ave_freq[NRINGS];      /* The flux weighted average of frequency of photons hitting each annulus */
  double w[NRINGS];             /* The radiative weight of the photons that hit the disk */
  double t_hit[NRINGS];         /* The effective T of photons hitting the disk */
  int nphot[NRINGS];            /*The number of photons created in each annulus */
  int nhit[NRINGS];             /*The number of photons which hit each annulus */
}
disk, qdisk;                    /* disk defines zones in the disk which in a specified frequency band emit equal amounts
                                   of radiation. disk gets reinitialized whenever the frequency interval of interest
                                   is changed.  qdisk stores the amount of heating of the disk as a result of
                                   illumination by the star or wind. It's boundaries are fixed throughout a cycle */

/* the next structure is intended to store a non standard temperature
   profile for the disk
   */

#define NBLMODEL 5000

struct blmodel
{
  int n_blpts;
  double r[NBLMODEL];
  double t[NBLMODEL];
}
blmod;


/*
 * The next structure is associated with reverberation mappping.
    SWN 6-2-15
    Wind paths is defined per cell and contains a binned array holding the spectrum of paths. Layers are
    For each frequency:
      For each path bin:
        What's the total fluxback of all these photons entering the cell?
*/
typedef struct wind_paths
{
  double *ad_path_flux;         //Array[by frequency, then path] of total flux of photons with the given v&p
  double *ad_path_flux_disk;
  double *ad_path_flux_wind;
  double *ad_path_flux_cent;    // As above, by source
  int *ai_path_num;             //Array [by frequency, then path] of the number of photons in this bin
  int *ai_path_num_disk;
  int *ai_path_num_wind;
  int *ai_path_num_cent;        // As above, by source
  double d_flux, d_path;        //Total flux, average path
  int i_num;                    //Number of photons hitting this cell
} wind_paths_dummy, *Wind_Paths_Ptr;

/* 	This structure defines the wind.  The structure w is allocated in the main
	routine.  The total size of the structure will be NDIM x MDIM, and the two
	dimenssions do not need to be the same.  The order of the
	cells in the structure is such that the as you increse the cell number by one
	z increases the fastest.

57+ -- The wind structure was reduced to contain just information about the geometry.  
Variables for wind cells that actually have volume in the wind are now in plasmamain, 
or macromain.  The wind structure still contains a volume element, which is the volume
of that cell in the wind, This is used in some cases to determine whether the cell
has any portion in the wind.  

Note that most of the macro atom information is in a separate structure.  This was
done to make it easier to control the size of the entire structure   06jul-ksl

 */
#define NIONIZ	5               /*The number of ions (normally H and He) for which one separately tracks ionization 
                                   and recombinations */


/* 061104 -- 58b -- ksl -- Added definitions to characterize whether a cell is in the wind. */
/* 110810 -- ksl - these are assigned to w->inwind, and are used to help determine how photons 
that go through a cell are treated.  Most of the assignments are based on whether the
volume in the wind is zero, which is determined in cylind_volumes for the cylindrical wind
and in correpsonding reoutines elsewehere.  These routines are called from the routine define_wind.
W_IGNORE is currently set in define_wind itself.  The values of these variables are used
in translate_in_wind.

Note that where_in_wind has been modified to use some of the same returns.  In general, the idea
is that if a new component is to be added, it should be added with by with two varialles ALL in whatever
and PART in whatever, as n and n+1
*/

typedef struct wind
{
  int ndom;                     /*The domain associated with this element of the wind */
  int nwind;                    /*A self-reference to this cell in the wind structure */
  int nplasma;                  /*A cross refrence to the corresponding cell in the plasma structure */
  double x[3];                  /*position of inner vertex of cell */
  double xcen[3];               /*position of the "center" of a cell */
  double r, rcen;               /*radial location of cell (Used for spherical, spherical polar
                                   coordinates. */
  double theta, thetacen;       /*Angle of coordinate from z axis in degrees  */
  double dtheta, dr;            /* widths of bins, used in hydro import mode */
  struct cone wcone;            /* cone structure that defines the bottom edge of the cell in 
                                   CYLVAR coordinates */
  double v[3];                  /*velocity at inner vertex of cell.  For 2d coordinate systems this
                                   is defined in the xz plane */
  double v_grad[3][3];          /*velocity gradient tensor  at the inner vertex of the cell NEW */
  double div_v;                 /*Divergence of v at center of cell */
  double dvds_ave;              /* Average value of dvds */
  double dvds_max, lmn[3];      /*The maximum value of dvds, and the direction in a cell in cylindrical coords */
  double vol;                   /* valid volume of this cell (that is the volume of the cell that is considered
                                   to be in the wind.  This differs from the volume in the Plasma structure
                                   where the volume is the volume that is actually filled with material. */
  double dfudge;                /* A number which defines a push through distance for this cell, which replaces the
                                   global variable DFUDGE in many instances */
  enum inwind_enum
  { W_IN_DISK = -5, W_IN_STAR = -4, W_IGNORE = -2, W_NOT_INWIND = -1,
    W_ALL_INWIND = 0, W_PART_INWIND = 1, W_NOT_ASSIGNED = -999
  } inwind;
  Wind_Paths_Ptr paths, *line_paths;    // SWM 6-2-15 Path data struct for each cell
}
wind_dummy, *WindPtr;

WindPtr wmain;

/* Plasma is a structure that contains information about the properties of the
plasma in regions of the geometry that are actually included n the wind */

/* 70 - 1108 - Define wavelengths in which to record gross spectrum in a cell, see also xave_freq and xj in plasma structure */
/* ksl - It would probably make more sense to define these in the same ways that bands are done for the generation of photons, or to
 * do both the same way at least.  Choosing to do this in two different ways makes the program confusing. The structure that has
 * the photon generation is called xband */
/* 78 - 1407 - NSH - changed several elements (initially those of size nions) in the plasma array to by dynamically allocated.
They are now pointers in the array. */


typedef struct plasma
{
  int nwind;                    /*A cross reference to the corresponding cell in the  wind structure */
  int nplasma;                  /*A self reference to this  in the plasma structure */
  double ne;                    /* electron density in the shell */
  double rho;                   /*density at the center of the cell. For clumped models, this is rho of the clump */
  double vol;                   /* volume of this cell (more specifically the volume  that is filled with material
                                   which can differe from the valid volume of the cell due to clumping. */
  double *density;              /*The number density of a specific ion.  This needs to correspond
                                   to the ion order obtained by get_atomic_data. 78 - changed to dynamic allocation */
  double *partition;            /*The partition function for each  ion. 78 - changed to dynamic allocation */
  double *levden;               /*The number density (occupation number?) of a specific level */

  double kappa_ff_factor;       /* Multiplicative factor for calculating the FF heating for a photon. */


  double *recomb_simple;        /* "alpha_e - alpha" (in Leon's notation) for b-f processes in simple atoms. */
  double *recomb_simple_upweight;       /* multiplicative factor to account for ratio of total to "cooling" energy for b-f processes in simple atoms. */

/* Beginning of macro information */
  double kpkt_emiss;            /*This is the specific emissivity due to the conversion k-packet -> r-packet in the cell
                                   in the frequency range that is required for the final spectral synthesis. (SS) */

  double kpkt_abs;              /* k-packet equivalent of matom_abs. (SS) */

  int *kbf_use;                 /* List of the indices of the photoionization processes to be used for kappa_bf. (SS) */
  int kbf_nuse;                 /* Total number of photoionization processes to be used for kappa_bf. (SS) */

/* End of macro information */


  double t_r, t_r_old;          /*radiation temperature of cell */
  double t_e, t_e_old;          /*electron temperature of cell */
  double dt_e, dt_e_old;        /*How much t_e changed in the previous iteration */
  double heat_tot, heat_tot_old;        /* heating from all sources */
  double abs_tot;
  double heat_lines, heat_ff;
  double heat_comp;             /* 1108 NSH The compton heating for the cell */
  double heat_ind_comp;         /* 1205 NSH The induced compton heatingfor the cell */
  double heat_lines_macro, heat_photo_macro;    /* bb and bf heating due to macro atoms. Subset of heat_lines 
                                                   and heat_photo. SS June 04. */
  double heat_photo, heat_z;    /*photoionization heating total and of metals */
  double heat_auger;            /* photoionization heating due to inner shell ionizations */
  double abs_photo, abs_auger;  /* this is the energy absorbed from the photon due to these processes - different from 
                                   the heating rate because of the binding energy */
  double w;                     /*The dilution factor of the wind */

  int ntot;                     /*Total number of photon passages */

  /*  counters of the number of photon passages by origin */

  int ntot_star;
  int ntot_bl;
  int ntot_disk;
  int ntot_wind;
  int ntot_agn;


  int nscat_es;                 /* The number of electrons scatters in the cell */
  int nscat_res;                /* The number of resonant line scatters in the cell */

  double mean_ds;               /* NSH 6/9/12 Added to allow a check that a thin shell is really optically thin */
  int n_ds;                     /* NSH 6/9/12 Added to allow the mean ds to be computed */
  int nrad;                     /* Total number of photons created within the cell */
  int nioniz;                   /* Total number of photon passages by photons capable of ionizing H */
  double *ioniz, *recomb;       /* Number of ionizations and recombinations for each ion.
                                   The sense is ionization from ion[n], and recombinations 
                                   to each ion[n].  */
  double *inner_ioniz, *inner_recomb;
  int *scatters;                /* The number of scatters in this cell for each ion. */
  double *xscatters;            /* Diagnostic measure of energy scattered out of beam on extract. */
  double *heat_ion;             /* The amount of energy being transferred to the electron pool
                                   by this ion via photoionization. */
  double *heat_inner_ion;       /* The amount of energy being transferred to the electron pool
                                   by this ion via photoionization. */
  double *cool_rr_ion;          /* The amount of energy being released from the electron pool
                                   by this ion via recombination. */
  double *lum_rr_ion;           /* The recombination luminosity
                                   by this ion via recombination. */

  double *cool_dr_ion;
  double j, ave_freq;           /*Respectively mean intensity, intensity_averaged frequency, 
                                   luminosity and absorbed luminosity of shell */
  double xj[NXBANDS], xave_freq[NXBANDS];       /* 1108 NSH frequency limited versions of j and ave_freq */
  double fmin[NXBANDS];         /* the minimum frequency photon seen in a band - this is incremented during photon flight */
  double fmax[NXBANDS];         /* the maximum frequency photon seen in a band - this is incremented during photon flight */
  double fmin_mod[NXBANDS];     /* the minimum freqneucy that the model should be applied for */
  double fmax_mod[NXBANDS];     /* the maximum frequency that the model should be applied for */

  /* banded, directional fluxes - last element is used for the sum of magnitude of (flux)*/
  double F_vis[4];
  double F_UV[4];
  double F_Xray[4];

  /* The term direct here means from photons which have not been scattered. These are photons which have been
     created by the central object, or the disk, or in the simple case the wind, but which have not undergone
     any kind of interaction which would change their direction
   */
  double j_direct, j_scatt;     /* 1309 NSH mean intensity due to direct photons and scattered photons */
  double ip_direct, ip_scatt;   /* 1309 NSH mean intensity due to direct photons and scattered photons */
  double xsd_freq[NXBANDS];     /* 1208 NSH the standard deviation of the frequency in the band */
  int nxtot[NXBANDS];           /* 1108 NSH the total number of photon passages in frequency bands */
  double max_freq;              /* 1208 NSH The maximum frequency photon seen in this cell */
  double cool_tot;              /*The total cooling in a cell */
  /* The total luminosity of all processes in the cell, basically the emissivity of the cell times it volume. Not the same 
     as what escapes the cell, since photons can interact within the cell and lose weight or even be destroyed */
  double lum_lines, lum_ff, cool_adiabatic;
  double lum_rr, lum_rr_metals; /* the radiative recombination luminosity - not the same as the cooling rate */
  double cool_comp;             /* The compton luminosity of the cell */
  double cool_di;               /* The direct ionization luminosity */
  double cool_dr;               /* The dielectronic recombination luminosity of the cell */
  double cool_rr, cool_rr_metals;       /*fb luminosity & fb of metals metals */
  double lum_tot, lum_tot_old;  /* The specific radiative luminosity in frequencies defined by freqmin
                                   and freqmax.  This will depend on the last call to total_emission */

  double cool_tot_ioniz;
  double lum_lines_ioniz, lum_ff_ioniz, cool_adiabatic_ioniz;
  double lum_rr_ioniz;
  double cool_comp_ioniz;       /* The compton luminosity of the cell */
  double cool_di_ioniz;         /* The direct ionization luminosity */
  double cool_dr_ioniz;         /* The dielectronic recombination luminosity of the cell */
  double cool_rr_ioniz, cool_rr_metals_ioniz;   /*fb luminosity & fb of metals metals */
  double lum_tot_ioniz;         /* The specfic radiative luminosity in frequencies defined by freqmin
                                   and freqmax.  This will depend on the last call to total_emission */

  double heat_shock;            /*1805 ksl - An extra heating term added to allow for shock heating of the plasma (Implementef for FU Ori Project */

  /* JM 1807 -- these routines are for the BF_SIMPLE_EMISSIVITY_APPROACH
     they allow one to inspect the net flow of energy into and from the simple ion 
     ionization pool */
  double bf_simple_ionpool_in, bf_simple_ionpool_out;

  double comp_nujnu;            /* 1701 NSH The integral of alpha(nu)nuj(nu) used to 
                                   compute compton cooling-  only needs computing once per cycle 
                                 */

  double dmo_dt[3];             /*Radiative force of wind */
  double rad_force_es[4];       /*Radiative force of wind - 4th element is sum of magnitudes*/
  double rad_force_ff[4];       /*Radiative force of wind - 4th element is sum of magnitudes*/
  double rad_force_bf[4];       /*Radiative force of wind - 4th element is sum of magnitudes*/



  double gain;                  /* The gain being used in iterations of the structure */
  double converge_t_r, converge_t_e, converge_hc;       /* Three measures of whether the program believes the grid is converged.
                                                           The first two are the fractional changes in t_r, t_e between this and the last cycle. The third
                                                           number is the fraction between heating and cooling divided by the sum of the 2       */
  int trcheck, techeck, hccheck;        /* NSH the individual convergence checks used to calculate converge_whole.  Each of these values
                                           is 0 if the fractional change or in the case of the last check error is less than a value, currently
                                           set to 0.05.  ksl 111126   
                                           NSH 130725 - this number is now also used to say if the cell is over temperature - it is set to 2 in this case   */
  int converge_whole, converging;       /* converge_whole is the sum of the individual convergence checks.  It is 0 if all of the convergence checks indicated
                                           convergence. converging is an indicator of whether the program thought the cell is on the way to convergence 0
                                           implies converging */

#define CELL_CONVERGING 0       /* Indicator for a cell which is considered converging - temperature is oscillating and decreasing */
#define CELL_NOT_CONVERGING 1   /* Indicator for a cell which is considered not converging (temperature is shooting off in one direction) */
#define CONVERGENCE_CHECK_PASS 0        /* Indicator for that the cell has passed a convergence check */
#define CONVERGENCE_CHECK_FAIL 1        /* Indicator for that the cell has failed a convergence check */
#define CONVERGENCE_CHECK_OVER_TEMP 2   /* Indicator for a cell that its electron temperature is more than TMAX */

  /* 1108 Increase sim estimators to cover all of the bands */
  /* 1208 Add parameters for an exponential representation, and a switch to say which we prefer. */
  enum spec_mod_type_enum
  {
    SPEC_MOD_PL = 1,
    SPEC_MOD_EXP = 2,
    SPEC_MOD_FAIL = -1
  } spec_mod_type[NXBANDS];     /* A switch to say which type of representation we are using for this band in this cell. 
                                   Negative means we have no useful representation, 0 means power law, 1 means exponential */

  double pl_alpha[NXBANDS];     /*Computed spectral index for a power law spectrum representing this cell */
  double pl_log_w[NXBANDS];     /*This is the log version of the power law weight. It is in an attempt to allow very large 
                                   values of alpha to work with the PL spectral model to avoide NAN problems. 
                                   The pl_w version can be deleted once testing is complete */


  double exp_temp[NXBANDS];     /* The effective temperature of an exponential representation of the radiation field in a cell */
  double exp_w[NXBANDS];        /* The prefactor of an exponential representation of the radiation field in a cell */
  double ip;                    /* Ionization parameter calculated as number of photons over the lyman limit entering a cell, divided by the number density of hydrogen for the cell */
  double xi;                    /* Ionization parameter as defined by Tartar et al 1969 and described in Hazy. Its the ionizing flux over the number of hydrogen atoms */
} plasma_dummy, *PlasmaPtr;

PlasmaPtr plasmamain;

/* A storage area for photons.  The idea is that it is sometimes time-consuming to create the
cumulative distribution function for a process, but trivial to create more than one photon 
of a particular type once one has the cdf,  This appears to be case for f fb photons.  But 
the same procedure could be used for line photons */

#define NSTORE 10
typedef struct photon_store
{
  int n;                        /* This is the photon number that was last used */
  double t, f1, f2, freq[NSTORE];

} photon_store_dummy, *PhotStorePtr;

PhotStorePtr photstoremain;

/* A second photon store: this is very similar to photon_store above but for use in generating macro atom bf photons from cfds*/
typedef struct matom_photon_store
{
  int n;                        /* This is the photon number that was last used */
  double t, nconf, freq[NSTORE];

} matom_photon_store_dummy, *MatomPhotStorePtr;

MatomPhotStorePtr matomphotstoremain;
#define MATOM_BF_PDF 1000       //number of points to use in a macro atom bf PDF

typedef struct macro
{
  double *jbar;
  /* This will store the Sobolev mean intensity in transitions which is needed 
     for Macro Atom jumping probabilities. The indexing is by configuration (the 
     NLTE_LEVELS) and then by the upward bound-bound jumps from that level 
     (the NBBJUMPS) (SS) */

  double *jbar_old;

  double *gamma;
  /* This is similar to the jbar but for bound-free transitions. It records the 
     appropriate photoionisation rate co-efficient. (SS) */

  double *gamma_old;

  double *gamma_e;
  /* This is Leon's gamma_e: very similar to gamma but energy weighted. Needed
     for division of photoionisation energy into excitation and k-packets. (SS) */

  double *gamma_e_old;

  double *alpha_st;
  /* Same as gamma but for stimulated recombination rather than photoionisation. (SS) */

  double *alpha_st_old;

  double *alpha_st_e;
  /* Same as gamma_e but for stimulated recombination rather than photoionisation. (SS) */

  double *alpha_st_e_old;

  double *recomb_sp;
  /* Spontaneous recombination. (SS) */

  double *recomb_sp_e;
  /* "e" version of the spontaneous recombination coefficient. (SS) */

  double *matom_emiss;
  /* This is the specific emissivity due to the de-activation of macro atoms in the cell
     in the frequency range that is required for the final spectral synthesis. (SS) */

  double *matom_abs;
  /* This is the energy absorbed by the macro atom levels - recorded during the ionization 
     cycles and used to get matom_emiss (SS) */

  /* This portion of the macro structure  is not written out by windsave */
  int kpkt_rates_known;

  double *cooling_bf;
  double *cooling_bf_col;
  double *cooling_bb;

  /* set of cooling rate stores, which are calculated for each macro atom each cycle,
     and used to select destruction rates for kpkts */
  double cooling_normalisation;
  double cooling_bbtot, cooling_bftot, cooling_bf_coltot;
  double cooling_ff, cooling_ff_lofreq;
  double cooling_adiabatic;     // this is just cool_adiabatic / vol / ne


} macro_dummy, *MacroPtr;

MacroPtr macromain;

int xxxpdfwind;                 // When 1, line luminosity calculates pdf

int size_Jbar_est, size_gamma_est, size_alpha_est;

#define TMAX_FACTOR			1.5     /*Factor by which t_e can exceed
                                                   t_r in order for absorbed to 
                                                   match emitted flux */

#define TMAX    5e8             /*NSH 130725 - this is the maximum temperature permitted - this was 
                                   introduced following problems with adiabatically heated cells increasing 
                                   forever. The value was suggested by DP as a sensible compton teperature for the PK05/P05 Zeus models. */
#define TMIN   100              /* A minimum temperature below which various emission processes are set to 0 */


//These constants are used in the various routines which compute ionization state
#define SAHA 4.82907e15         /* 2* (2.*PI*MELEC*k)**1.5 / h**3  (Calculated in constants) */
#define MAXITERATIONS	200     //The number of loops to do to try to converge in ne
#define FRACTIONAL_ERROR 0.03   //The change in n_e which causes a break out of the loop for ne
#define THETAMAX	 1e4    //Used in initial calculation of n_e
#define MIN_TEMP	100.    //  ??? this is another minimum temperature, which is used in saha.c and variable_temperature.c )

// these definitions are for various ionization modes
#define IONMODE_ML93_FIXTE 0    // Lucy Mazzali using existing t_e (no HC balance)
#define IONMODE_LTE_TR 1        // LTE using t_r
#define IONMODE_LTE_TE 4        // LTE using t_e
#define IONMODE_FIXED 2         // Hardwired concentrations
#define IONMODE_ML93 3          // Lucy Mazzali
#define IONMODE_MATRIX_BB 8     // matrix solver BB model
#define IONMODE_MATRIX_SPECTRALMODEL 9  // matrix solver spectral model based on power laws
#define IONMODE_MATRIX_ESTIMATORS 10    // matrix solver spectral model based on power laws

// and the corresponding modes in nebular_concentrations
#define NEBULARMODE_TR 0        // LTE using t_r
#define NEBULARMODE_TE 1        // LTE using t_e
#define NEBULARMODE_ML93 2      // ML93 using correction
#define NEBULARMODE_NLTE_SIM 3  // Non_LTE with SS modification (Probably could be removed)
#define NEBULARMODE_LTE_GROUND 4        // A test mode which forces all levels to the GS (Probably could be removed)
#define NEBULARMODE_PAIRWISE_ML93 6     // pairwise ML93 (diffuse BB)
#define NEBULARMODE_PAIRWISE_SPECTRALMODEL 7    // pairwise spectral models (power law or expoentials)
#define NEBULARMODE_MATRIX_BB 8 // matrix solver BB model
#define NEBULARMODE_MATRIX_SPECTRALMODEL 9      // matrix solver spectral model
#define NEBULARMODE_MATRIX_ESTIMATORS 10        // matrix solver spectral model



typedef struct photon
{
  double x[3];                                  /* The position of packet */
  double lmn[3];                                /* Direction cosines of the packet */
  double freq, freq_orig, freq_orig_loc;        /* current, original frequency redshifted and unredshifted) of this packet */
  double w, w_orig;                             /* current and original weight of this packet */
  double tau;                                   /* optical depth of the photon since its creation or last interaction */
  enum istat_enum
  {
    P_INWIND = 0,               //in wind,
    P_SCAT = 1,                 //in process of scattering,
    P_ESCAPE = 2,               //Escaped to reach the universe,
    P_HIT_STAR = 3,             //absorbed by photosphere of star,
    P_TOO_MANY_SCATTERS = 4,    //in wind after MAXSCAT scatters
    P_ERROR = 5,                //Trying to scatter a photon in a location where it should not scatter
    P_ABSORB = 6,               //Photoabsorbed within wind
    P_HIT_DISK = 7,             //Banged into disk
    P_SEC = 8,                  //Photon hit secondary
    P_ADIABATIC = 9,            //records that a photon created a kpkt which was destroyed by adiabatic cooling
    P_ERROR_MATOM = 10,         //Some kind of error in processing of a photon which excited a macroattom
    P_LOFREQ_FF = 11,           //records a photon that had too low a frequency
    P_REPOSITION_ERROR = 12     //A photon passed through the disk due to dfudge pushing it through incorrectly
  } istat;                      /*status of photon. */

  enum frame
  {
    F_LOCAL = 0,
    F_OBSERVER = 1
  } frame;

  int nscat;                    /*Number of scatters for this photon */
  int nres;                     /*For line scattering, indicates the actual transition; 
                                   for continuum scattering, meaning 
                                   depends on matom vs non-matom. See headers of emission.c 
                                   or matom.c for details. */
  int nnscat;                   /* Used for the thermal trapping model of
                                   anisotropic scattering to carry the number of
                                   scattering to "extract" when needed for wind
                                   generated photons SS05. */
  int nrscat;                   /* number of resonance scatterings */
  int grid;                     /*grid position of the photon in the wind, if
                                   the photon is in the wind.  If the photon is not
                                   in the wind, then -1 implies inside the wind cone and  
                                   -2 implies outside the wind */

  enum origin_enum
  { PTYPE_STAR = 0,
    PTYPE_BL = 1,
    PTYPE_DISK = 2,
    PTYPE_WIND = 3,
    PTYPE_AGN = 4,
    PTYPE_STAR_MATOM = 10,
    PTYPE_BL_MATOM = 11,
    PTYPE_DISK_MATOM = 12,
    PTYPE_WIND_MATOM = 13,
    PTYPE_AGN_MATOM = 14
  } origin, origin_orig;        /* Where this photon originated.  If the photon has
                                   scattered its "origin" may be changed to "wind". */
  /* note that we add 10 to origin when processed by a macro-atom
     which means we need these values in the enum list.  In making spectra in spectrum_create
     10 is subtracted from the types.  If ever this logic is changed one must the be careful
     that it is fixed in create_spectra as well.  

     Comment - ksl - 180712 - The logic for all of this is obscure to me, since we keep track of the
     photons origin separately.  At some point one might want to revisit the necessity for this
   */
  int np;                       /* The photon number, which used ease tracking a photon for diagnostic
                                   purposes */
  double path;                  /* The total path length of a photon (used for reverberation calcuations) */
  double ds;                    /* the distance a photon has moved since its creation or last interaction */
}
p_dummy, *PhotPtr;

PhotPtr photmain;               /* A pointer to all of the photons that have been created in a subcycle. Added to ease 
                                   breaking the main routine of python into separate rooutines for inputs and running the
                                   program */

    /* minimum value for tau for p_escape_from_tau function- below this we 
       set to p_escape_ to 1 */
#define TAU_MIN 1e-6


    /* 68b - ksl - This is a structure in which the history of a single photon bundle can be recorded
     * See phot_util   phot_hist().  It needs to be used carefully.  if phot_hist_on is true
     * then photon histories will be attempted.
     */

#define MAX_PHOT_HIST	1000
int n_phot_hist, phot_hist_on, phot_history_spectrum;
struct photon xphot_hist[MAX_PHOT_HIST];

struct basis
{
  double a[3][3];

}
basis_cartesian;


    /* The next section defines the spectrum arrays.  The spectrum structure contains
       the selection criteria for each spectrum as well as the array in which to store the
       spectrum.  The first MSPEC spectra are reserved for the total generated spectrum,
       total emitted spectrum, the total spectrum of photons which are scattered and 
       the total spectrum of photons which are absorbed.  The remainder of the spectra pertain 
       to the spectrum as seen from various directions. Note that the space for the spectra 
       are allocated using a calloc statement in spectrum_init.  1409-ksl-A new spectrum
       was added.  This will be the first of the spectra.  It is simply the generated spectrum
       before passing through the wind.  It has the orginal weights as generated.  */



#define SPECTYPE_RAW        0   // As written this produces L_nu and so to get a Luminosity one needs to integrate
#define SPECTYPE_FLAMBDA    1
#define SPECTYPE_FNU        2
#define D_SOURCE 100.0          // distance to the source in parsecs for genearating spectra

int nspectra;                   /* After create_spectrum, the number of elements allocated for s, or
                                   alternatively the number of spectra one has to work with.  Note that
                                   general s[0],s[1] and s[2] are the escaping, scattered and absorbed photons,
                                   while elements higher than this will contain spectra as seen by different observers */

#define MSPEC               7   /* The number of standard spectra - i.e. not user defined angles */
#define SPEC_CREATED        0   /* The spectrum of the original weight before transmission through the wind */
#define SPEC_EMITTED        1   /* The emitted spectrum - i.e. photons with their weights changed by transmission through the wind */
#define SPEC_CENSRC         2   /* The emitted spectrum from photons emitted from the central source (if there is one) */
#define SPEC_DISK           3   /* The emitted spectrum from photons emitted from the disk (if there is one) */
#define SPEC_WIND           4   /* The emitted spectrum from photons emitted from the wind itself */
#define SPEC_HITSURF        5   /* The spectrum for photons which hit the a surface and were absorbed - should be zero for when reflection
                                 * is turned on */
#define SPEC_SCATTERED      6   /* The spectrum of photons which were scattered at least once in the wind - the weight used is the final
                                 * weight after transmission through the wind */

int nscat[MAXSCAT + 1], nres[MAXSCAT + 1], nstat[NSTAT];

typedef struct spectrum
{
  char name[40];
  float freqmin, freqmax, dfreq;
  float lfreqmin, lfreqmax, ldfreq;     /* NSH 1302 - values for logarithmic spectra */
  double lmn[3];
  double mmax, mmin;            /* Used only in live or die situations, mmax=cos(angle-DANG_LIVE_OR_DIE)
                                   and mmim=cos(angle+DANG_LIVE_OR_DIE).   In actually defining this
                                   one has to worry about signs and exactly whether you are above
                                   or below the plane */
  double renorm;                /* Factor used in Live or die option where it is 4*PI/domega;
                                   1.0 otherwise */
  int nphot[NSTAT];             /* Used to help define what happened to photons when trying to construct this
                                   particular spectrum */
  int nscat;                    /* nscat>999 -> select everything
                                   0 < nscat < MAXScat select only those photons which have
                                   scattered nscat times, number of scattering photons in spectrum, if nscat is negative
                                   then all photons with more than |nscat| are included.  */
  int top_bot;                  /* 0 ->select photons regardless of location 
                                   >0     -> select only photons whose "last" position is above the disk
                                   <0    -> select only photons whose last position is below the disk */
  double x[3], r;               /* The position and radius of a special region from which to extract spectra. 
                                   x is taken to be the center of the region and r is taken to be the radius of
                                   the region.   */
  double f[NWAVE];              /* The spectrum in linear (wavelength or frequency) units */
  double lf[NWAVE];             /* The specturm in log (wavelength or frequency)  units  */
  double lfreq[NWAVE];          /* We need to hold what freqeuncy intervals our logarithmic spectrum has been taken over */

  double f_wind[NWAVE];         /* The spectrum of photons created in the wind or scattered in the wind. Created for 
                                   reflection studies but possibly useful for other reasons as well. */
  double lf_wind[NWAVE];        /* The logarithmic version of this */
}
spectrum_dummy, *SpecPtr;


SpecPtr xxspec;

/* Parameters used only by py_wind
 * py_wind_projecti	0 -> simply print the various parameters without 
 * 			atempting toproject onto a yz plane
 * 			1 -> project onto a yz plane.
 */

int py_wind_min, py_wind_max, py_wind_delta, py_wind_project;
double *aaa;                    // A pointer to an array used by py_wind

/* This is the structure for storing cumulative distribution functions. The CDFs are
generated from a function which is usually only proportional to the probability density
function or from an array.  It is sometimes useful, e.g. in calculating the reweighting function to
have access to the proper normalization.  


*/


#define NCDF 30000              //The default size for these arrays.  This needs to be greater than
                                //the size of any model that is read in, hence larger than NWAVE in models.h
#define FUNC_CDF  200           //The size for CDFs made from functional form CDFs
#define ARRAY_PDF 1000          //The size for PDFs to be turned into CDFs from arrays


typedef struct Cdf
{
  double x[NCDF];               /* Positions for which the CDF is calculated */
  double y[NCDF];               /* The value of the CDF at x */
  double d[NCDF];               /* 57i -- the rate of change of the CDF at x */
  double limit1, limit2;        /* Limits (running from 0 to 1) that define a portion
                                   of the CDF to sample */
  double x1, x2;                /* limits if they exist on what is returned */
  double norm;                  //The scaling factor which would renormalize the CDF
  int ncdf;                     /* Size of this CDF */
}
 *CdfPtr, cdf_dummy;

struct Cdf cdf_ff;
struct Cdf cdf_fb;
struct Cdf cdf_vcos;
struct Cdf cdf_bb;
struct Cdf cdf_brem;




/* Variable used to allow something to be printed out the first few times
   an event occurs */
int itest, jtest;

char hubeny_list[132];          //Location of listing of files representing hubeny atmospheres






/* These variables are stored or used by the routines for anisotropic scattering */
/* Allow for the transfer of tau info to scattering routine */



/* N.B. cdf_randwind and phot_randwind are used in the routine anisowind for 
as part of effort to incorporate anisotropic scattering in to python.  
Added for python_43.2 */


/* Provide generally for having arrays which descibe the 3 xyz axes. 
these are initialized in main, and used in anisowind  */


double x_axis[3];
double y_axis[3];
double z_axis[3];

/* These are structures associated with frequency limits/s used for photon
generation and for calculating heating and cooling */

#define NBANDS 20
struct xbands
{
  double f1[NBANDS], f2[NBANDS];
  double alpha[NBANDS];
  double pl_const[NBANDS];
  double min_fraction[NBANDS];
  double nat_fraction[NBANDS];  // The fraction of the accepted luminosity in this band
  double used_fraction[NBANDS];
  double flux[NBANDS];          //The "luminosity" within a band
  double weight[NBANDS];
  int nphot[NBANDS];
  int nbands;                   // Actual number of bands in use
}
xband;


/* The next section contains the freebound structures that can be used for both the
 * specific emissivity of a free-bound transition, and for the recombination coefficient
 * assuming the array has been initialized, which can take a few minutes
*/

#define NTEMPS	60              // The number of temperatures which are stored in each fbstruct
                                /* NSH this was increased from 30 to 60 to take account of 3 extra OOM 
                                   intemperature we wanted to have in fb */
#define NFB	20              // The maximum number of frequency intervals for which the fb emission is calculated

struct fbstruc
{
  double f1, f2;
  double cool[NIONS][NTEMPS];   //cooling rate due to radiative recombination
  double lum[NIONS][NTEMPS];    //emissivity due to radiative recombinaion
  double cool_inner[NIONS][NTEMPS];     //cooling rate due to recombinations to inner shells
}
freebound[NFB];

double xnrecomb[NIONS][NTEMPS]; // There is only one set of recombination coefficients
double xninnerrecomb[NIONS][NTEMPS];    // There is only one set of recombination coefficients

double fb_t[NTEMPS];
int nfb;                        // Actual number of freqency intervals calculated




#include "version.h"            /*54f -- Added so that version can be read directly */
#include "templates.h"
#include "recipes.h"

/* kap_bf stores opacities for a single cell and as calculated by the routine kappa_bf. 
 * It was made an external array to avoid having to pass it between various calling routines
 * but this means that one has to be careful that data is not stale.  It is required for 
 * macro-atoms where bf is a scattering process, but not for the simple case.
 */

double kap_bf[NLEVELS];



// 12jun nsh - some commands to enable photon logging in given cells. There is also a pointer in the geo

FILE *pstatptr;                 //NSH 120601 - pointer to a diagnostic file that will contain photon data for given cells
int cell_phot_stats;            //1=do  it, 0=dont do it
#define  NCSTAT 10              //The maximum number of cells we are going to log
int ncstat;                     // the actual number we are going to log
int ncell_stats[NCSTAT];        //the numbers of the cells we are going to log


/* Added variables which count number of times two situations occur (See #91) */
int nerr_no_Jmodel;
int nerr_Jmodel_wrong_freq;



// advanced mode variables
struct advanced_modes
{
  /* these are all 0=off, 1=yes */
  int iadvanced;                // this is controlled by the -d flag, global mode control.
  int extra_diagnostics;        // when set various extra files will be written out depending what one wants to check
  int save_cell_stats;          // want to save photons statistics by cell
  int keep_ioncycle_windsaves;  // want to save wind file each ionization cycle
  int make_tables;              // create tables showing various parameters for each cycle
  int track_resonant_scatters;  // want to track resonant scatters
  int save_photons;             // want to track photons (in photon2d)
  int save_extract_photons;     // we want to save details on extracted photons
  int adjust_grid;              // the user wants to adjust the grid scale
  int diag_on_off;              // extra diagnostics
  int use_debug;                // print out debug statements
  int print_dvds_info;          // print out information on the velocity gradients
  int keep_photoabs;            // keep photoabsorption in final spectrum
  int quit_after_inputs;        // quit after inputs read in, testing mode
  int fixed_temp;               // do not alter temperature from that set in the parameter file
  int zeus_connect;             // We are connecting to zeus, do not seek new temp and output a heating and cooling file
  int rand_seed_usetime;        // default random number seed is fixed, not based on time
  int photon_speedup;
}
modes;


FILE *optr;                     //pointer to a diagnostic file that will contain dvds information



/* Structure containing all of the file and directory names created */
struct filenames
{
  char root[LINELENGTH];        // main rootname
  char windsave[LINELENGTH];    // wind save filename
  char old_windsave[LINELENGTH];        // old windsave name
  char specsave[LINELENGTH];    // spec save filename
  char diag[LINELENGTH];        // diag file
  char diagfolder[LINELENGTH];  // diag folder
  char input[LINELENGTH];       // input name if creating new pf file
  char new_pf[LINELENGTH];      // name of generated pf file
  char wspec[LINELENGTH];       // .spec_tot file (spectrum from last ionization cycle) 
  char lwspec[LINELENGTH];      // .log_spec_tot file (spectra from last ionization cycle in log wavelength scale)  
  char wspec_wind[LINELENGTH];  // .spec_tot_wind (spectra of wind photons in last ionization cycle on a linear scale limited to wind photons
  char lwspec_wind[LINELENGTH]; // .log_spec_tot_wind (same as above but in log units) 
  char spec[LINELENGTH];        // .spec file (extracted spectra on linear scale)
  char lspec[LINELENGTH];       // .spec file (extracted spectra on a log scale)
  char spec_wind[LINELENGTH];   // .spec file (extracted spectra limited to wind photons on a linear scale)
  char lspec_wind[LINELENGTH];  // .spec file (extracted spectra limited to wind photons on a log scale)
  char disk[LINELENGTH];        // disk diag file name
  char tprofile[LINELENGTH];    // non standard tprofile fname
  char phot[LINELENGTH];        // photfile e.g. python.phot
  char windrad[LINELENGTH];     // wind rad file
}
files;



#define NMAX_OPTIONS 20

/* these two variables are used by xdefine_phot() in photon_gen.c 
   to set the mode for get_matom_f()in matom.c and tell it 
   whether it has already calculated the matom emissivities or not. */
#define CALCULATE_MATOM_EMISSIVITIES 0
#define USE_STORED_MATOM_EMISSIVITIES 1


/* modes for kpkt calculations */
#define KPKT_MODE_CONTINUUM  0  /* only account for k->r processes */
#define KPKT_MODE_ALL        1  /* account for all cooling processes */

/* this variable controls whether to use the 
   Altered mode for bound-free in "simple-macro mode" */
#define BF_SIMPLE_EMISSIVITY_APPROACH 1


/* Variable introducted to cut off macroatom / estimator integrals when exponential function reaches extreme values. Effectivevly a max limit imposed on x = hnu/kT terms */
#define ALPHA_MATOM_NUMAX_LIMIT 30      /* maximum value for h nu / k T to be considered in integrals */
#define ALPHA_FF 100.           // maximum h nu / kT to create the free free CDF


/* non-radiative heat flow mode */
#define KPKT_NET_HEAT_MODE 0


/* DIAGNOSTIC for understanding problems with imported models
 */


#define BOUND_NONE   0
#define BOUND_INNER_CONE  1
#define BOUND_OUTER_CONE  2
#define BOUND_RMAX 3
#define BOUND_RMIN 4
#define BOUND_ZMIN 5
#define BOUND_ZMAX 6
#define BOUND_INNER_RHO 7
#define BOUND_OUTER_RHO 8

int xxxbound;


/* Structures associated with rdchoice.  This 
 * shtructure is required only in cases where one 
 * wants to use rdchoice multiple times with different
 * options.  But it can, or course be used for 
 * any such call.  It allows one to assoicated
 * a input word with an output value, using the routine
 * get_choices.  There needs to be one structure
 * for each input variable.  At present, this is only
 * used for the selection of spec_types
 */

typedef struct rdpar_choices
{
  char choices[10][LINELENGTH];
  int vals[10];
  int n;
} dummy_choices, *ChoicePtr;

struct rdpar_choices zz_spec;
