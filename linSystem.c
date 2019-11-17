  /*! \file linSystem.c
      \brief Implementation of linear equation system related functions

      Uses the Gauss-Seidel method for solving systems of linear equations.
  */

    #include <math.h>
    #include <stdlib.h>
    #include <stdio.h>
    #include <likwid.h>
    #include "linSystem.h"
    #include "utils.h"

  /*!
    \brief Implementation of the Gauss-Seidel method

    \param S structure that stores a linear system of equations
    \param x solution array (contains initial guess on function call)
    \param[out] avg_time average iteration time
    \param[out] norm array conatining the norm for each iteration
  */
  int gaussSeidel(linSystem_t *S, int nx, real_t *x, double *avg_time, real_t *norm) {
    real_t *res = malloc(sizeof(real_t) * S->n);
    int k, i;
    double t1;

    *avg_time = 0.0;

    for (k = 0; k < MAXIT; k++) {
        t1 = timestamp();
      	LIKWID_MARKER_START("GS");

        /* First iteration */
        res[0] = ((S->sd[0] * x[1]) + (S->rsd[0] * x[nx]));

        /* Next nx iteration */
        for(i = 1; i < nx; i++) {
        	res[i] = ((S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]) + (S->rsd[i] * x[i+nx]));
        }

        /* From nx to S->n - nx */
        for(i = nx; i < S->n - nx-1; i++) {
        	res[i] = ((S->rid[i] * x[i-nx]) + (S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]) + (S->rsd[i] * x[i+nx]));
        }

        /* From S->n - nx to S->n */
        for(i = S->n - nx-1; i < S->n-1; i++) {
        	res[i] = ((S->rid[i-nx]) + (S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]));
        }

        /* Last iteration */
        res[S->n-1] = ((S->rid[S->n-1] * x[S->n - nx-1]) + (S->id[S->n-1] * x[S->n-1]));

        for(i = 0; i < S->n; i++) {
            res[i] = S->b[i] - res[i];
            x[i] = res[i] / S->md[i];
        }
	LIKWID_MARKER_STOP("GS");
    LIKWID_MARKER_START("Norm");
        norm[k] = normL2(S, x, res, nx);
    LIKWID_MARKER_STOP("Norm");

        *avg_time += timestamp() - t1;
    }
    *avg_time = *avg_time / MAXIT;
    free(res);
    return 0;
  }

  /*!
      \brief Calculates the norm of the array r = S->b - S->A*x
      \param S structure that stores a linear system of equations
      \param x solution array
  */

real_t normL2(linSystem_t *S, real_t *x, real_t *res, int nx) {
    int i, j;
    real_t part_norm;
    part_norm = 0.0;

    /* First iteration */
    res[0] = ((S->sd[0] * x[1]) + (S->rsd[0] * x[nx]) + (S->md[0] * x[0]));

    /* Next nx iteration */
    for(i = 1; i < nx; i++) {
        res[i] = ((S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]) + (S->rsd[i] * x[i+nx]) + (S->md[i] * x[i]));
    }

    /* From nx to S->n - nx */
    for(i = nx; i < S->n - nx-1; i++) {
        res[i] = ((S->rid[i] * x[i-nx]) + (S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]) + (S->rsd[i] * x[i+nx]) + (S->md[i] * x[i]));
    }

    /* From S->n - nx to S->n */
    for(i = S->n - nx-1; i < S->n-1; i++) {
        res[i] = ((S->rid[i-nx]) + (S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]) + (S->md[i] * x[i]));
    }

    /* Last iteration */
    res[S->n-1] = ((S->rid[S->n-1] * x[S->n - nx-1]) + (S->id[S->n-1] * x[S->n-1]) + (S->md[S->n-1] * x[S->n-1]));

    for(i = 0; i < S->n; i++) {
        res[i] = S->b[i] - res[i];
        res[i] = res[i]*res[i];
    }
    for(i = 0; i < S->n; i++) {
        part_norm += res[i];
    }
    return (real_t)sqrt(part_norm);         // returns norm of residue array
}
