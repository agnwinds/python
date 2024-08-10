
/***********************************************************/
/** @file  windsave2fits.c 
 * @author ksl
 * @date   October, 2021
 *
 * @brief  Routines to inspect variables in a wind structure 
 *
 *###Notes###
 * This is intended just as a diagnostic routine 
 * so that one can print out whatever variables in
 * a windstrucutre one wants in order to diagnose
 * a problem.  It was written so that we could inspect
 * some of the macro atom variables in paralell mode
 * in diagnosing issue #898 and #910, but anyon 
 * should change it so that other problems might 
 * be addressed.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "fitsio.h"


char inroot[LINELENGTH], outroot[LINELENGTH], model_file[LINELENGTH], folder[LINELENGTH];
int model_flag, ksl_flag, cmf2obs_flag, obs2cmf_flag;

double line_matom_lum_single (double lum[], PlasmaPtr xplasma, int uplvl);
int line_matom_lum (int uplvl);
int create_matom_level_map ();

// Define a structure to hold spectral data
typedef struct
{
  int num_wavelengths;
  int num_spectra;
  float **data;                 // 2D array to hold spectra data
} Spectra;

/* Next fucntions are utilities that it should be possitle to call repeatedly*/

// Function to prepare 2D image data from a 2D array
float *
prepare_image_data (float **data, int width, int height)
{
  // Allocate memory for a contiguous 1D array to store the 2D image data
  float *image_data = calloc (width * height, sizeof (float));
  if (!image_data)
  {
    fprintf (stderr, "Memory allocation error\n");
    return NULL;
  }

  // Flatten the 2D data into a 1D array for FITS writing
  for (int i = 0; i < height; i++)
  {
    for (int j = 0; j < width; j++)
    {
      image_data[i * width + j] = data[i][j];
    }
  }

  return image_data;
}

// Function to prepare binary table data
void
prepare_binary_table (int num_rows, int num_columns, float **table_data)
{
  // Allocate memory for each row of the table
  for (int i = 0; i < num_rows; i++)
  {
    table_data[i] = malloc (num_columns * sizeof (float));
    if (!table_data[i])
    {
      fprintf (stderr, "Memory allocation error\n");
      return;
    }

    // Fill each row with some sample data
    for (int j = 0; j < num_columns; j++)
    {
      table_data[i][j] = (float) (i * num_columns + j); // Example data
    }
  }
}



// Function to write a 2D image to a FITS extension with a given name
int
write_image_extension (fitsfile *fptr, float *image_data, int width, int height, const char *extname)
{
  int status = 0;               // CFITSIO status value MUST be initialized to zero
  long naxes[2];                // Array to store dimension sizes

  // Set the dimensions of the 2D image: [width, height]
  naxes[0] = width;             // Number of columns
  naxes[1] = height;            // Number of rows

  // Create a new image extension of type FLOAT_IMG
  if (fits_create_img (fptr, FLOAT_IMG, 2, naxes, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Set the extension name using the EXTNAME keyword
  if (fits_update_key (fptr, TSTRING, "EXTNAME", (void *) extname, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the image data to the FITS extension
  if (fits_write_img (fptr, TFLOAT, 1, width * height, image_data, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  return 0;
}

// Function to write a binary table to a FITS extension with a given name
int
write_binary_table_extension (fitsfile *fptr, int num_rows, int num_columns, float **table_data, const char *extname)
{
  int status = 0;
  char *ttype[] = { "COL1" };   // Name for each column
  char *tform[] = { "E" };      // Format for each column (E = float)
  char *tunit[] = { "unit" };   // Unit for each column

  // Create a new binary table extension
  if (fits_create_tbl (fptr, BINARY_TBL, num_rows, num_columns, ttype, tform, tunit, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Set the extension name using the EXTNAME keyword
  if (fits_update_key (fptr, TSTRING, "EXTNAME", (void *) extname, NULL, &status))
  {
    fits_report_error (stderr, status);
    return status;
  }

  // Write the table data
  for (int i = 0; i < num_rows; i++)
  {
    if (fits_write_col (fptr, TFLOAT, 1, i + 1, 1, num_columns, table_data[i], &status))
    {
      fits_report_error (stderr, status);
      return status;
    }
  }

  return 0;
}






int
make_spec (inroot)
     char *inroot;
{

  fitsfile *fptr;               // Pointer to the FITS file
  int status = 0;               // CFITSIO status value MUST be initialized to zero

  // Create a new FITS file
  if (fits_create_file (&fptr, "!multi_extension.fits", &status))
  {
    fits_report_error (stderr, status); // Print any error message
    return status;
  }

// Create an empty primary HDU
  if (fits_create_img (fptr, FLOAT_IMG, 0, NULL, &status))
  {
    fits_report_error (stderr, status);
  }

  printf ("Hello World %s \n", inroot);
  printf ("Plasma  %d \n", NPLASMA);
  printf ("NBINS in spec %d \n", NBINS_IN_CELL_SPEC);


  Spectra spectra;
  spectra.num_spectra = NPLASMA;
  spectra.num_wavelengths = NBINS_IN_CELL_SPEC;
  spectra.data = calloc (spectra.num_spectra, sizeof (float *));

  for (int i = 0; i < spectra.num_spectra; i++)
  {
    spectra.data[i] = calloc (spectra.num_wavelengths, sizeof (float));

    // Fill with sample data
    for (int j = 0; j < spectra.num_wavelengths; j++)
    {
      // spectra.data[i][j] = (float) (i + j);
      spectra.data[i][j] = (float) plasmamain[i].cell_spec_flux[j];
    }
  }

  // Prepare image data for the first extension
  float *image_data1 = prepare_image_data (spectra.data, spectra.num_wavelengths, spectra.num_spectra);
  if (!image_data1)
  {
    // Free allocated memory before exiting
    for (int i = 0; i < spectra.num_spectra; i++)
    {
      free (spectra.data[i]);
    }
    free (spectra.data);
    return 1;
  }


  status = write_image_extension (fptr, image_data1, spectra.num_wavelengths, spectra.num_spectra, "Jnu");

  Spectra freq;
  freq.num_spectra = 1;
  freq.num_wavelengths = NBINS_IN_CELL_SPEC;
  freq.data = calloc (freq.num_spectra, sizeof (float *));
  freq.data[0] = calloc (freq.num_wavelengths, sizeof (float));

  for (int i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    freq.data[0][i] = (float) pow (10., geo.cell_log_freq_min + i * geo.cell_delta_lfreq);
  }


  float *image_data2 = prepare_image_data (freq.data, NBINS_IN_CELL_SPEC, 1);
/*
  double freq[2000][1];
  int i;

  for (i = 0; i < NBINS_IN_CELL_SPEC; i++)
  {
    freq[i][0] = (float) pow (10., geo.cell_log_freq_min + i * geo.cell_delta_lfreq);
  }

  float *image_data2 = prepare_image_data (freq, NBINS_IN_CELL_SPEC, 1);
*/
  status = write_image_extension (fptr, image_data2, NBINS_IN_CELL_SPEC, 1, "nu");

  if (fits_close_file (fptr, &status))
  {
    fits_report_error (stderr, status);
  }

  return 0;
}

/**********************************************************/
/**
 * @brief      parses the command line options
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * The general purpose of each of the command line options
 * should be fairly obvious from reading the code.
 *
 *
 * Although this routine uses the standard Log and Error commands
 * the diag files have not been created yet and so this information
 * is really simply written to the terminal.
 *
 **********************************************************/

int
xparse_command_line (argc, argv)
     int argc;
     char *argv[];
{
  int j = 0;
  int i;
  char dummy[LINELENGTH];
  int mkdir ();
  char *fgets_rc;


  sprintf (outroot, "%s", "new");

  model_flag = ksl_flag = obs2cmf_flag = cmf2obs_flag = 0;

  if (argc == 1)
  {
    printf ("Parameter file name (e.g. my_model.pf, or just my_model):");
    fgets_rc = fgets (dummy, LINELENGTH, stdin);
    if (!fgets_rc)
    {
      printf ("Input rootname is NULL or invalid\n");
      exit (1);
    }
    get_root (inroot, dummy);
  }
  else
  {

    for (i = 1; i < argc; i++)
    {
      if (strcmp (argv[i], "-out_root") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("python: Expected out_root after -out_root switch\n");
          exit (0);
        }

        get_root (outroot, dummy);
        i++;
        j = i;

      }
      if (strcmp (argv[i], "-model_file") == 0)
      {
        if (sscanf (argv[i + 1], "%s", dummy) != 1)
        {
          printf ("python: Expected a model file containing density, velocity and temperature after -model_file switch\n");
          exit (0);
        }
        get_root (model_file, dummy);
        i++;
        j = i;
        printf ("got a model file %s\n", model_file);
        model_flag = 1;
      }
      else if (strcmp (argv[i], "-ksl") == 0)
      {
        printf ("Carrying out a simple hard wired ion modification\n");
        ksl_flag = 1;
      }
      else if (strcmp (argv[i], "--dry-run") == 0)
      {
        modes.quit_after_inputs = 1;
        j = i;
      }
      else if (strcmp (argv[i], "-cmf") == 0)
      {
        obs2cmf_flag = 1;
      }
      else if (strcmp (argv[i], "-obs") == 0)
      {
        cmf2obs_flag = 1;
      }
      else if (strncmp (argv[i], "-", 1) == 0)
      {
        printf ("python: Unknown switch %s\n", argv[i]);
        exit (0);
      }
    }

    /* The last command line variable is always the windsave file */

    if (j + 1 == argc)
    {
      printf ("All of the command line has been consumed without specifying a parameter file name, so exiting\n");
      exit (1);
    }
    strcpy (dummy, argv[argc - 1]);
    get_root (inroot, dummy);

  }

  return (0);
}



/**********************************************************/
/**
 * @brief      the main routine which carries out the effort
 *
 * @param [in]  int  argc   the number of command line arguments
 * @param [in]  char *  argv[]   The command line arguments
 *
 *
 * ###Notes###
 *
 * This routine oversees the effort.  The basic steps are
 *
 * - parse the command line to get the names of files
 * - read the old windsave file
 * - read the densities from in this case H
 * - modify the densities
 * - write out the new windsave file
 *
 *
 **********************************************************/



int
main (argc, argv)
     int argc;
     char *argv[];
{

  char infile[LINELENGTH], outfile[LINELENGTH];
  FILE *fopen ();
  int mkdir ();


  xparse_command_line (argc, argv);

  sprintf (infile, "%.150s.wind_save", inroot);
  sprintf (outfile, "%.150s.txt", inroot);

  printf ("Reading %s and writing to %s\n", infile, outfile);

  zdom = calloc (MAX_DOM, sizeof (domain_dummy));
  if (zdom == NULL)
  {
    printf ("Unable to allocate memory for domain\n");
    return EXIT_FAILURE;
  }

  wind_read (infile);

  make_spec (inroot);

  return (0);
}
