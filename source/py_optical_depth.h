/** ************************************************************************* */
/**
 * @file     py_optical_depth.h
 * @author   Edward Parkinson
 * @date     February 2021
 *
 * @brief    Globals used during optical depth diagnostics.
 *
 * @details
 *
 * Do not include this in any file other than py_optical_depth_sub.c, otherwise
 * there will be double definition errors. If you want to include it in
 * multiple, place MAXDIFF and N_PI_EDGES into py_optical_depth_sub.c to avoid
 * the errors. The reason they are in this header file, and why it exists, is
 * to keep things a bit cleaner in py_optical_depth_sub.c
 *
 * ************************************************************************** */

#define N_FREQ_BINS 10000

// Inclination angles structure

typedef struct SightLines_s
{
  char name[50];
  double lmn[3];
} SightLines_t;

// PI edges structure

typedef struct Edges_s
{
  char name[50];
  double freq;
} Edges_t;

// Control the column density extracted

enum COLUMN_DENSITY
{
  COLUMN_MODE_RHO,
  COLUMN_MODE_ION,
};

int COLUMN_MODE;
int COLUMN_MODE_ION_NUMBER;

// Control which domain to send photons from

int N_DOMAIN;

// Control how the optical depth integration is done

enum {
  RUN_MODE_OUTWARD,
  RUN_MODE_PHOTOSPHERE,
} MODE;

double TAU_DEPTH;

/*
 * Functions from other files
 */

void control_program(void);
SightLines_t *initialize_inclination_angles (int *n_angles);
void print_optical_depths (SightLines_t *inclinations, int n_inclinations, Edges_t edges[], int n_edges, double *optical_depth_values,
                           double *column_density_values);
void write_optical_depth_spectrum (SightLines_t * inclinations, int n_inclinations, double *tau_spectrum, double freq_min, double d_freq);
int create_photon (PhotPtr p_out, double freq, double *lmn);
