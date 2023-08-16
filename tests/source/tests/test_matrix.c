/* ****************************************************************************************************************** */
/**
 *  @file test_matrix.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define LENGTH 256
#define FRAC_TOLERANCE 0.01

int solve_matrix (double *a_matrix, double *b_matrix, int size, double *x_matrix, int nplasma);
int invert_matrix (double *matrix, double *inverse_matrix, int matrix_size);

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

static int
get_invert_matrix_test_data (const char *matrix_path, const char *inverse_path, double **matrix, double **inverse, int *size)
{
  FILE *fp_matrix = fopen (matrix_path, "r");
  if (!fp_matrix)
  {
    fprintf (stderr, "unable to open matrix file %s\n", matrix_path);
    return EXIT_FAILURE;
  }

  int size_from_matrix = 0;
  fscanf (fp_matrix, "%d", &size_from_matrix);
  *size = size_from_matrix;

  *matrix = malloc (size_from_matrix * size_from_matrix * sizeof (double));
  *inverse = malloc (size_from_matrix * sizeof (double));

  int i;
  for (i = 0; i < size_from_matrix * size_from_matrix; ++i)
  {
    fscanf (fp_matrix, "%le", &(*matrix)[i]);
  }
  fclose (fp_matrix);

  FILE *fp_inverse = fopen (inverse_path, "r");
  if (!fp_inverse)
  {
    fprintf (stderr, "unable to open matrix file %s\n", inverse_path);
    return EXIT_FAILURE;
  }

  fscanf (fp_inverse, "%*d");
  for (i = 0; i < size_from_matrix; ++i)
  {
    fscanf (fp_inverse, "%le", &(*inverse)[i]);
  }

  fclose (fp_inverse);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

static int
get_solve_matrix_test_data (const char *a_path, const char *b_path, const char *x_path, double **a_matrix, double **b_vector,
                            double **x_vector, int *size)
{
  FILE *fp_a = fopen (a_path, "r");
  if (!fp_a)
  {
    fprintf (stderr, "unable to open matrix file %s\n", a_path);
    return EXIT_FAILURE;
  }

  int size_from_a = 0;
  fscanf (fp_a, "%d", &size_from_a);
  *size = size_from_a;

  *a_matrix = malloc (size_from_a * size_from_a * sizeof (double));
  *b_vector = malloc (size_from_a * sizeof (double));
  *x_vector = malloc (size_from_a * sizeof (double));

  int i;
  for (i = 0; i < size_from_a * size_from_a; ++i)
  {
    fscanf (fp_a, "%le", &(*a_matrix)[i]);
  }
  fclose (fp_a);

  FILE *fp_b = fopen (b_path, "r");
  if (!fp_b)
  {
    fprintf (stderr, "unable to open matrix file %s\n", b_path);
    return EXIT_FAILURE;
  }

  FILE *fp_x = fopen (x_path, "r");
  if (!fp_x)
  {
    fprintf (stderr, "unable to open matrix file %s\n", x_path);
    return EXIT_FAILURE;
  }

  fscanf (fp_b, "%*d");
  fscanf (fp_x, "%*d");
  for (i = 0; i < size_from_a; ++i)
  {
    fscanf (fp_b, "%le", &(*b_vector)[i]);
    fscanf (fp_x, "%le", &(*x_vector)[i]);
  }

  fclose (fp_b);
  fclose (fp_x);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

static int
internal_test_invert (int test_num)
{
  const char *python_path = getenv ((const char *) "PYTHON");
  if (!python_path)
  {
    fprintf (stderr, "$PYTHON env variable not set\n");
    return EXIT_FAILURE;
  }

  double *matrix;
  double *inverse;
  char matrix_filepath[LENGTH];
  char inverse_filepath[LENGTH];

  sprintf (matrix_filepath, "%s/tests/test_data/matrices/invert_test%d/matrix.txt", python_path, test_num);
  sprintf (inverse_filepath, "%s/tests/test_data/matrices/invert_test%d/inverse.txt", python_path, test_num);

  int size;
  const int get_err = get_invert_matrix_test_data (matrix_filepath, inverse_filepath, &matrix, &inverse, &size);

  if (get_err)
    return EXIT_FAILURE;

  double *test_inverse = malloc (size * sizeof (double));
  const int matrix_err = invert_matrix (matrix, test_inverse, size);

  int i;
  for (i = 0; i < size; ++i)
  {
    const double frac = fabs ((inverse[i] - test_inverse[i]) / test_inverse[i]);
    if (frac > FRAC_TOLERANCE)
    {
      fprintf (stderr, "internal_test_invert failure");
      return EXIT_FAILURE;
    }
  }

  free (matrix);
  free (inverse);
  free (test_inverse);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

static int
internal_test_solve (const int test_num)
{
  const char *python_path = getenv ((const char *) "PYTHON");
  if (!python_path)
  {
    fprintf (stderr, "$PYTHON env variable not set\n");
    return EXIT_FAILURE;
  }

  double *matrix_a;
  double *vector_b;
  double *vector_x;
  char matrix_a_filepath[LENGTH];
  char vector_b_filepath[LENGTH];
  char vector_x_filepath[LENGTH];

  sprintf (matrix_a_filepath, "%s/tests/test_data/matrices/solve_test%d/A.txt", python_path, test_num);
  sprintf (vector_b_filepath, "%s/tests/test_data/matrices/solve_test%d/b.txt", python_path, test_num);
  sprintf (vector_x_filepath, "%s/tests/test_data/matrices/solve_test%d/x.txt", python_path, test_num);

  int size;
  const int get_err =
    get_solve_matrix_test_data (matrix_a_filepath, vector_b_filepath, vector_x_filepath, &matrix_a, &vector_b, &vector_x, &size);

  if (get_err)
    return EXIT_FAILURE;

  double *test_vector_x = malloc (size * sizeof (double));
  const int matrix_err = solve_matrix (matrix_a, vector_b, size, test_vector_x, -1);

  int i;
  for (i = 0; i < size; ++i)
  {
    const double frac = fabs ((vector_x[i] - test_vector_x[i]) / test_vector_x[i]);
    if (frac > FRAC_TOLERANCE)
    {
      fprintf (stderr, "internal_test_solve failure\n");
      return EXIT_FAILURE;
    }
  }

  free (matrix_a);
  free (vector_b);
  free (vector_x);
  free (test_vector_x);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

int
test_solve_matrix (void)
{
  int error;

  error = internal_test_solve (1);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

int
test_invert_matrix (void)
{
  int error;

  error = internal_test_invert (1);

  return EXIT_SUCCESS;
}
