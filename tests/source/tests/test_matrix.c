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
#include <CUnit/CUnit.h>

#include "../source/matrix.h"

#define LENGTH 256
#define EPSILON 1.0e-6

/** *******************************************************************************************************************
 *
 * @brief Load in test data for the matrix inversion tests
 *
 * @param [in] matrix_path filepath to the input matrix for the tests
 * @param [in] inverse_path filepath to the input matrix for verification
 * @param [out] matrix pointer to pointer for matrix input
 * @param [out] inverse point to pointer for inverse used for verification
 * @param [out] size the size of the matrices
 *
 * @return an error status
 *
 * @details
 *
 * The matrix and inverse variables are pointers to pointers, as the memory is allocated for those variables in this
 * function.
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
  *inverse = malloc (size_from_matrix * size_from_matrix * sizeof (double));

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
  for (i = 0; i < size_from_matrix * size_from_matrix; ++i)
  {
    fscanf (fp_inverse, "%le", &(*inverse)[i]);
  }

  fclose (fp_inverse);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Load in test data for the solve matrix tests
 *
 * @param [in] a_path the filepath to the A matrix for input
 * @param [in] b_path the filepath to the b vector for input
 * @param [in] x_path the filepath to the x vector for verification
 * @param [out] a_matrix pointer to pointer for A matrix for input
 * @param [out] b_vector pointer to pointer for b vector for input
 * @param [out] x_vector pointer to pointer for x vector for verification
 * @param [out] size the size of the vectors
 *
 * @return an error status
 *
 * @details
 *
 * Pointer to pointers are used for a_matrix and etc. because the memory for those pointers/arrays are allocated here.
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
 * @brief Test case for invert_matrix
 *
 * @param [in] test_name the test name
 *
 * @return an error status
 *
 * @details
 *
 * Calls invert matrix with given input and checks the absolute tolerance.
 *
 *  ***************************************************************************************************************** */

static int
call_invert_matrix (const char *test_name)
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

  sprintf (matrix_filepath, "%s/tests/test_data/matrix/%s/matrix.txt", python_path, test_name);
  sprintf (inverse_filepath, "%s/tests/test_data/matrix/%s/inverse.txt", python_path, test_name);

  int matrix_size;
  const int get_err = get_invert_matrix_test_data (matrix_filepath, inverse_filepath, &matrix, &inverse, &matrix_size);
  CU_ASSERT_EQUAL_FATAL (get_err, EXIT_SUCCESS);

  double *test_inverse = malloc (matrix_size * matrix_size * sizeof (double));
  const int matrix_err = invert_matrix (matrix, test_inverse, matrix_size);
  CU_ASSERT_EQUAL_FATAL (matrix_err, EXIT_SUCCESS);

  int i;
  for (i = 0; i < matrix_size * matrix_size; ++i)
  {
    CU_ASSERT_DOUBLE_EQUAL_FATAL (test_inverse[i], inverse[i], EPSILON);
  }

  free (matrix);
  free (inverse);
  free (test_inverse);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief Test case for solve_matrix
 *
 * @param [in] test_name the test name
 *
 * @return an error status
 *
 * @details
 *
 * Calls solve_matrix with given input and checks the absolute tolerance.
 *
 *  ***************************************************************************************************************** */

int
call_solve_matrix (const char *test_name)
{
  const char *python_path = getenv ((const char *) "PYTHON");
  if (!python_path)
  {
    fprintf (stderr, "$PYTHON variable not set\n");
    exit (EXIT_FAILURE);
  }

  double *matrix_a;
  double *vector_b;
  double *vector_x;
  char matrix_a_filepath[LENGTH];
  char vector_b_filepath[LENGTH];
  char vector_x_filepath[LENGTH];

  sprintf (matrix_a_filepath, "%s/tests/test_data/matrix/%s/A.txt", python_path, test_name);
  sprintf (vector_b_filepath, "%s/tests/test_data/matrix/%s/b.txt", python_path, test_name);
  sprintf (vector_x_filepath, "%s/tests/test_data/matrix/%s/x.txt", python_path, test_name);

  int vector_size;
  const int get_err =
    get_solve_matrix_test_data (matrix_a_filepath, vector_b_filepath, vector_x_filepath, &matrix_a, &vector_b, &vector_x, &vector_size);
  CU_ASSERT_EQUAL_FATAL (get_err, EXIT_SUCCESS);

  double *test_vector_x = malloc (vector_size * sizeof (double));
  const int matrix_err = solve_matrix (matrix_a, vector_b, vector_size, test_vector_x, -1);
  CU_ASSERT_EQUAL_FATAL (matrix_err, EXIT_SUCCESS);

  int i;
  for (i = 0; i < vector_size; ++i)
  {
    CU_ASSERT_DOUBLE_EQUAL_FATAL (test_vector_x[i], vector_x[i], EPSILON);
  }

  free (matrix_a);
  free (vector_b);
  free (vector_x);
  free (test_vector_x);

  return EXIT_SUCCESS;
}

/** *******************************************************************************************************************
 *
 * @brief
 *
 * @details
 *
 *  ***************************************************************************************************************** */

void
test_solve_matrix (void)
{
  call_solve_matrix ("small_matrix");
  call_solve_matrix ("matrix_ion");
}

/** *******************************************************************************************************************
 *
 * @brief
 *
 * @details
 *
 *  ***************************************************************************************************************** */

void
test_invert_matrix (void)
{
  call_invert_matrix ("inverse_small");
}

/** *******************************************************************************************************************
 *
 * @brief
 *
 * @details
 *
 *  ***************************************************************************************************************** */

void
create_matrix_test_suite (void)
{
  CU_pSuite suite;

  /* It's a pain to have two suites for GSL and cuSolver, so let's keep it conditionally
     compiled */
#ifdef CUDA_ON
  suite = CU_add_suite ("Matrix Functions Suite: cuSolver", cusolver_create, cusolver_destroy);
#else
  suite = CU_add_suite ("Matrix Functions Suite: GSL", NULL, NULL);
#endif

  if (suite == NULL)
  {
    fprintf (stderr, "Failed to create `Matrix Functions Suite` suite\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }

  /* Add CPU tests to suite */
  if ((CU_add_test (suite, "Solve Matrix", test_solve_matrix) == NULL) ||
      (CU_add_test (suite, "Invert Matrix", test_invert_matrix) == NULL))
  {
    fprintf (stderr, "Failed to add tests to `Matrix Functions Suite: CPU`\n");
    CU_cleanup_registry ();
    exit (CU_get_error ());
  }
}
