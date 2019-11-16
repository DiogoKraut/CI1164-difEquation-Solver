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
    real_t res, part_norm;
    int k, i;
    double t1;

    *avg_time = 0.0;

    for (k = 0; k < MAXIT; k++) {
        t1 = timestamp();
  	LIKWID_MARKER_START("GS");
        part_norm = 0.0;
        /* First iteration */
        res  = ((S->sd[0] * x[1]) + (S->rsd[0] * x[nx]));
        x[0] = (S->b[0] - res) / S->md[0];
        res += S->md[0] * x[0];
        // printf(" 0 %.3lf", res);
        res  = S->b[0] - res;
        part_norm += res*res;

        /* Next nx iteration */
        for(i = 1; i < nx; i++) {
        	res  = ((S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]) + (S->rsd[i] * x[i+nx]));
            x[i] = (S->b[i] - res) / S->md[i];
            res += S->md[i] * x[i];
            // printf(" %d %.3lf",i, res);
            res  = S->b[i] - res;
            part_norm += res*res;
        }

        /* From nx to S->n - nx */
        for(i = nx; i < S->n - nx-1; i++) {
        	res  = ((S->rid[i] * x[i-nx]) + (S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]) + (S->rsd[i] * x[i+nx]));
            x[i] = (S->b[i] - res) / S->md[i];
            res += S->md[i] * x[i];
            // printf(" %d %.3lf",i, res);
            res  = S->b[i] - res;
            part_norm += res*res;
        }

        /* From S->n - nx to S->n */
        for(i = S->n - nx-1; i < S->n-1; i++) {
        	res = ((S->rid[i-nx]) + (S->id[i] * x[i-1]) + (S->sd[i] * x[i+1]));
            x[i] = (S->b[i] - res) / S->md[i];
            res += S->md[i] * x[i];
            // printf(" %d %.3lf",i, res);
            res  = S->b[i] - res;
            part_norm += res*res;
        }

        /* Last iteration */
        res = ((S->rid[S->n-1] * x[S->n - nx-1]) + (S->id[S->n-1] * x[S->n-1]));
        x[S->n-1] = (S->b[S->n-1] - res) / S->md[S->n-1];
        res += S->md[S->n-1] * x[S->n-1];
        // printf(" l %.3lf", res);
        res  = S->b[S->n-1] - res;
        part_norm += res*res;

        // for(i = 0; i < S->n; i++)
        //     printf("%.3lf ", x[i]);
        // printf("%lf ", part_norm);
        // printf("\n");

        norm[k] = sqrt(part_norm);
        *avg_time += timestamp() - t1;
	LIKWID_MARKER_STOP("GS");
    }
    *avg_time = *avg_time / MAXIT;
  	LIKWID_MARKER_CLOSE;
    return 0;
  }

  // /*!
  //     \brief Calculates the norm of the array r = S->b - S->A*x
  //     \param S structure that stores a linear system of equations
  //     \param x solution array
  // */

  // real_t normL2(linSystem_t *S, real_t *x) {
  //     int i, j;
  //     real_t norm, res;
  //     norm = 0.0;
  //     // Calculate residue array
  //     for(i = 0; i < S->n; i++) {
  //         res = 0.0;
  //         for(j = 0; j < S->n; j++) {     // ith line of matrix A times x
  //             res += (S->A[i][j] * x[j]);
  //         }

  //         res = S->b[i] - res;            // partial residue
  //         norm += res*res;                // sum of all partial residues

  //     }
  //     return (real_t)sqrt(norm);         // returns norm of residue array
  // }
