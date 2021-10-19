#include "atomic.h" 

 
int nelements;                  /* The actual number of ions read from the data file */
int nions;                      /*The actual number of ions read from the datafile */
int nlevels;                    /*These are the actual number of levels which were read in */
int nlte_levels;                /* Actual number of levels to treat explicityly */
int nlevels_macro;              /* Actual number of macro atom levels. (SS, June 04) */
int nlines;                     /* Actual number of lines that were read in */
int nlines_macro;               /* Actual number of Macro Atom lines that were read in.  New version of get_atomic
                                   data assumes that macro lines are read in before non-macro lines */
int n_inner_tot;                /*The actual number of inner shell ionization cross sections in total */

int nauger;                     /*Actual number of innershell edges for which autoionization is to be computed */

ElemPtr ele;

double rho2nh;                  /* Conversion constant from rho to nh the number of H atoms per unit volume */

IonPtr ion;

ConfigPtr config;

LinePtr line, lin_ptr[NLINES];  /* line[] is the actual structure array that contains all the data, *lin_ptr
                                   is an array which contains a frequency ordered set of ptrs to line */
struct lines fast_line;

int nline_min, nline_max, nline_delt;   /* Used to select a range of lines in a frequency band from the lin_ptr array 
                                           in situations where the frequency range of interest is limited, including for defining which
                                           lines come into play for resonant scattering along a line of sight, and in
                                           calculating band_limit luminosities.  The limits are established by the
                                           routine limit_lines.
                                         */

int n_coll_stren;

Coll_stren coll_stren[NLINES];  

int nxphot;                     /*The actual number of ions for which there are VFKY photoionization x-sections */
double phot_freq_min;           /*The lowest frequency for which photoionization can occur */
double inner_freq_min;          /*The lowest frequency for which inner shell ionization can take place */

int ntop_phot;                  /* The actual number of TopBase photoionzation x-sections */
int nphot_total;                /* total number of photoionzation x-sections = nxphot + ntop_phot */

Topbase_phot phot_top[NLEVELS];
TopPhotPtr phot_top_ptr[NLEVELS];       /* Pointers to phot_top in threshold frequency order - this */

Topbase_phot inner_cross[N_INNER * NIONS];
TopPhotPtr inner_cross_ptr[N_INNER * NIONS];

Inner_elec_yield inner_elec_yield[N_INNER * NIONS];

Inner_fluor_yield inner_fluor_yield[N_INNER * NIONS];

int ndrecomb;                   //This is the actual number of DR parameters

Drecomb drecomb[NIONS];         //set up the actual structure

double dr_coeffs[NIONS];        //this will be an array to temprarily store the volumetric dielectronic recombination rate coefficients for the current cell under interest. We may want to make this 2D and store the coefficients for a range of temperatures to interpolate.

