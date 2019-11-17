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
    real_t res;
    int k, i, n = S->n;
    double t1;

    *avg_time = 0.0;

    for (k = 0; k < MAXIT; k++) {
        t1 = timestamp();
      	LIKWID_MARKER_START("GS");
        /* First iteration */
        res  = ((S->diag[0].sd * x[1]) + (S->diag[0].rsd * x[nx]));
        x[0] = (S->b[0] - res) / S->diag[0].md;

        /* Next nx iteration */
        for(i = 1; i < nx; i++) {
        	res  = ((S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]));
            x[i] = (S->b[i] - res) / S->diag[i].md;
        }

        /* From nx to n - nx */
        for(i = nx; i < n - nx-1; i++) {
        	res  = ((S->diag[i].rid * x[i-nx]) + (S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]));
            x[i] = (S->b[i] - res) / S->diag[i].md;
        }

        /* From n - nx to n */
        for(i = n - nx-1; i < n-1; i++) {
        	res = ((S->diag[i-nx].rid) + (S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]));
            x[i] = (S->b[i] - res) / S->diag[i].md;
        }

        /* Last iteration */
        res = ((S->diag[n-1].rid * x[n - nx-1]) + (S->diag[n-1].id * x[n-1]));
        x[n-1] = (S->b[n-1] - res) / S->diag[n-1].md;

    	LIKWID_MARKER_STOP("GS");
        norm[k] = normL2(S, x, nx);
        *avg_time += timestamp() - t1;
    }
    *avg_time = *avg_time / MAXIT;
  	LIKWID_MARKER_CLOSE;
    return 0;
  }

/*!
  \brief Calculates the norm of the array r = S->b - S->A*x
  \param S structure that stores a linear system of equations
  \param x solution array
*/

real_t normL2(linSystem_t *S, real_t *x, int nx) {
    int i, n = S->n;
    real_t part_norm, res;
    // Calculate residue array
    LIKWID_MARKER_START("Norm");
    part_norm = 0.0;
    /* First iteration */
    res = ((S->diag[0].sd * x[1]) + (S->diag[0].rsd * x[nx]) + (S->diag[0].md * x[0]));
    // printf(" 0 %.3lf", res);
    res = S->b[0] - res;
    part_norm += res*res;

    /* Next nx iteration */
    for(i = 1; i < nx; i++) {
        res = ((S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]) + (S->diag[i].md * x[i]));
        // printf(" %d %.3lf",i, res);
        res = S->b[i] - res;
        part_norm += res*res;
    }

    /* From nx to n - nx */
    for(i = nx; i < n - nx-1; i++) {
        res = ((S->diag[i].rid * x[i-nx]) + (S->diag[i].id * x[i-1]) + (S->diag[i].sd * x[i+1]) + (S->diag[i].rsd * x[i+nx]) + (S->diag[i].md * x[i]));
        // printf(" %d %.3lf",i, res);
        res = S->b[i] - res;
        part_norm += res*res;
    }

    /* From n - nx to n */
    for(i = n - nx-1; i < n-1; i++) {
        res  = S->diag[i].rid * x[i-nx];
        res += S->diag[i].id  * x[i-1];
        res += S->diag[i].sd  * x[i+1];
        res += S->diag[i].md  * x[i];
        // printf(" %d %.3lf",i, res);
        res = S->b[i] - res;
        part_norm += res*res;
    }

    /* Last iteration */
    res = ((S->diag[n-1].rid * x[n - nx-1]) + (S->diag[n-1].id * x[n-2]) + (S->diag[n-1].md * x[n-1]));
    // printf(" l %.3lf", res);
    res = S->b[n-1] - res;
    part_norm += res*res;

    LIKWID_MARKER_STOP("Norm");
    return (real_t)sqrt(part_norm);         // returns partial norm
}