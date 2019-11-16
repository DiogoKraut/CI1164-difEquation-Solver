#ifndef __linSystem_H__
#define __linSystem_H__

#define MAXIT 10
#define EPS 1.0e-4

#include <stdio.h>

typedef double real_t;

/*! Structure for storing a system of linear equations */
typedef struct {
  real_t *md, *sd, *id, *rid, *rsd; /*!< linear equation system's diagonals        */
  real_t *b;                        /*!< linear equation system's right-hand array */
  unsigned int n;                   /*!< dimension of the system                   */
} linSystem_t;

int gaussSeidel(linSystem_t *S, int nx, real_t *x, double *avg_time, real_t *norm);
real_t normL2(linSystem_t *S, real_t *x, int nx);
#endif
