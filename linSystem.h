#ifndef __linSystem_H__
#define __linSystem_H__

#define MAXIT 10
#define EPS 1.0e-4

#include <stdio.h>

typedef float real_t; /*! \type typedef for switching between float/double */

/*! Structure for storing a system of linear equations */
typedef struct {
  real_t **A;     /*!< linear equation system's coeficient matrix */
  real_t *b;      /*!< linear equation system's right-hand array  */
  unsigned int n; /*!< dimension of the system                    */
} linSystem_t;

int gaussSeidel(linSystem_t *S, real_t *x, double *avg_time, real_t *norm);
real_t normL2(linSystem_t *S, real_t *x);
#endif
