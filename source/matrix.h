/* ****************************************************************************************************************** */
/**
 *  @file matrix.h
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#ifdef CUDA_ON

int gpu_solve_linear_system (double *a_matrix, double *b_vector, int size, double *x_vector);
int gpu_invert_matrix (double *matrix, double *inverse, int num_rows);

#endif
