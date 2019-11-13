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
  int gaussSeidel(linSystem_t *S, real_t *x, double *avg_time, real_t *norm) {
  	LIKWID_MARKER_INIT;
    int k, i, j;
    double t1;
    real_t sum;

    *avg_time = 0.0;

    for (k = 0; k < MAXIT; k++) {
    for(int i = 0; i < S->n; i++) {
        printf("%lf ", x[i]);
    }printf("\n");
  	LIKWID_MARKER_START("GS");
        t1 = timestamp();
        for(i = 0; i < S->n; i++) {
            sum = 0.0;
            for(j = 0; j < S->n; j++) {
                if(j != i)
                    sum += S->A[i][j] * x[j];
            }
            x[i] = S->b[i] - sum;
            x[i] /= S->A[i][i];
        }

        *avg_time += timestamp() - t1;
        norm[k] = normL2(S, x);
	LIKWID_MARKER_STOP("GS");
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

  real_t normL2(linSystem_t *S, real_t *x) {
      int i, j;
      real_t norm, res;
      norm = 0.0;
      // Calculate residue array
      for(i = 0; i < S->n; i++) {
          res = 0.0;
          for(j = 0; j < S->n; j++) {     // ith line of matrix A times x
              res += (S->A[i][j] * x[j]);
          }

          res = S->b[i] - res;            // partial residue
          norm += res*res;                // sum of all partial residues

      }
      return (real_t)sqrt(norm);         // returns norm of residue array
  }
