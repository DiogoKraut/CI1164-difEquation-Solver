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
  	LIKWID_MARKER_INIT;
    real_t *res = malloc(sizeof(real_t) * S->n);
    int k, i, n = S->n;
    double t1;

    *avg_time = 0.0;

    for (k = 0; k < MAXIT; k++) {
        t1 = timestamp();
      	LIKWID_MARKER_START("GS");
        /* First iteration */
        res[0]  = ((S->diag[0].sd * x[1]) + (S->diag[0].rsd * x[nx]));

        /* Next nx iteration */
        for(i = 1; i < nx; i++) {
        	res[i]  = ((S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]));
        }

        /* From nx to n - nx */
        for(i = nx; i < n - nx-1; i++) {
        	res[i]  = ((S->diag[i].rid * x[i-nx]) + (S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]));
        }

        /* From n - nx to n */
        for(i = n - nx-1; i < n-1; i++) {
        	res[i] = ((S->diag[i-nx].rid) + (S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]));
        }

        /* Last iteration */
        res[n-1] = ((S->diag[n-1].rid * x[n - nx-1]) + (S->diag[n-1].id * x[n-1]));


        for(i = 0; i < n; i++) {
            res[i] = S->b[i] - res[i];
            x[i] = res[i] / S->diag[i].md;
        }
    	LIKWID_MARKER_STOP("GS");
        norm[k] = normL2(S, x, nx);
        *avg_time += timestamp() - t1;
    }
    *avg_time = *avg_time / MAXIT;
  	LIKWID_MARKER_CLOSE;
    free(res);
    return 0;
  }

/*!
  \brief Calculates the norm of the array r = S->b - S->A*x
  \param S structure that stores a linear system of equations
  \param x solution array
*/

real_t normL2(linSystem_t *S, real_t *x, int nx) {
    int i, n = S->n;
    real_t part_norm, *res = malloc(sizeof(real_t) * S->n);
    // Calculate residue array
    LIKWID_MARKER_START("Norm");
    part_norm = 0.0;
    /* First iteration */
    res[0] = ((S->diag[0].sd * x[1]) + (S->diag[0].rsd * x[nx]) + (S->diag[0].md * x[0]));
    res[0] = S->b[0] - res[0];
    part_norm += res[0]*res[0];

    /* Next nx iteration */
    for(i = 1; i < nx; i++) {
        res[i] = ((S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]) + (S->diag[i].md * x[i]));
        res[i] = S->b[i] - res[i];
        part_norm += res[i]*res[i];
    }

    /* From nx to n - nx */
    for(i = nx; i < n - nx-1; i++) {
        res[i] = ((S->diag[i].rid * x[i-nx]) + (S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]) + (S->diag[i].md * x[i]));
        res[i] = S->b[i] - res[i];
        part_norm += res[i]*res[i];
    }

    /* From n - nx to n */
    for(i = n - nx-1; i < n-1; i++) {
        res[i]  = S->diag[i].rid * x[i-nx];
        res[i] += S->diag[i].id  * x[i-1];
        res[i] += S->diag[i].sd  * x[i+1];
        res[i] += S->diag[i].md  * x[i];
        res[i] = S->b[i] - res[i];
        part_norm += res[i]*res[i];
    }

    /* Last iteration */
    res[n-1] = ((S->diag[n-1].rid * x[n - nx-1]) + (S->diag[n-1].id * x[n-2]) + (S->diag[n-1].md * x[n-1]));
    res[n-1] = S->b[n-1] - res[i];
    part_norm += res[i]*res[i];

    LIKWID_MARKER_STOP("Norm");
    free(res);
    return (real_t)sqrt(part_norm);         // returns partial norm
}