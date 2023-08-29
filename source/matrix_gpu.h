/** ******************************************************************************************************************
 *
 *  @file matrix.h
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#ifndef MATRIX_H
#define MATRIX_H

#ifdef CUDA_ON
/* initialise and destroy cuda/cusolver */
int cusolver_create (void);
int cusolver_destroy (void);

/* gpu matrix routines */
int gpu_solve_matrix (double *a_matrix, double *b_vector, int size, double *x_vector);
int gpu_invert_matrix (double *matrix, double *inverse, int num_rows);

/* error string */
const char *cusolver_get_error_string (int error);
#endif

#endif
