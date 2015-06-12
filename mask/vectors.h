#include <stdlib.h>
#include <stdio.h>

/* Some vector routines by Scott Ransom since he hates the non-zero */
/* offset general vector routines in Numerical Recipes.  ;)         */
/* The routines are in vectors.c                                    */

float *gen_fvect(long length);
/* Generate a floating point vector */

double *gen_dvect(long length);
/* Generate a double precision vector */

short *gen_svect(long length);
/* Generate an short integer vector */

int *gen_ivect(long length);
/* Generate an integer vector */

long *gen_lvect(long length);
/* Generate a long integer vector */

unsigned char *gen_bvect(long length);
/* Generate a 'byte' or unsigned character vector */

unsigned char **gen_bmatrix(long nrows, long ncols);
/* Generate a 'byte' or unsigned char matrix (2 dimensions) */

short **gen_smatrix(long nrows, long ncols);
/* Generate a short int matrix (2 dimensions) */

int **gen_imatrix(long nrows, long ncols);
/* Generate an integer matrix (2 dimensions) */

float **gen_fmatrix(long nrows, long ncols);
/* Generate a floating point matrix (2 dimensions) */

double **gen_dmatrix(long nrows, long ncols);
/* Generate a double precision matrix (2 dimensions) */

void vect_free(void *vect);
/* Free a generated vector */ 

/*  Note:  To free memory allocated by these routines simply use  */
/*         the free() function.                                   */
/*                                                                */
/*  Example:                                                      */
/*                                                                */
/*  x = gen_fvect(100);   // Generate a 100 point float vector.   */
/*  free(x);              // Free the vector.                     */
