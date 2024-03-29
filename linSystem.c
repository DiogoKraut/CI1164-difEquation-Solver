/*! \file linSystem.c
    \brief Implementation of linear equation system related functions

    Uses the Gauss-Seidel method for solving systems of linear equations.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "linSystem.h"
#include "utils.h"

/*!
  \brief Implementation of the Gauss-Seidel method

  \param S structure that stores a linear system of equations
  \param x solution array (contains initial guess on function call)
  \param error threshold for acceptable answer
  \param[out] avg_time average iteration time
  \param[out] norm array conatining the norm for each iteration
  \param[out] int_count counter of the number of iterations
  \return 0 if method converged, 1 otherwise
*/
int gaussSeidel(linSystem_t *S, real_t *x, double *avg_time, real_t *norm) {
    int k, i, j;
    double t1;
    real_t sum;

    *avg_time = 0.0;
    real_t *c = calloc(S->n, sizeof(real_t));


    for (k = 0; k < MAXIT; k++) {
        t1 = timestamp();
        for(i = 0; i < S->n; i++) {
            c[i] = S->b[i];
            sum = 0.0;
            for(j = 0; j < S->n; j++) {
                if(j < i)
                    sum += S->A[i][j] * c[j]; // use previously calculated value
                else if(j > i)
                    sum += S->A[i][j] * x[j]; // use current value
            }
            c[i] -= sum;
            c[i] /= S->A[i][i];
        }
        *avg_time += timestamp() - t1;
        for(i = 0; i < S->n; i++)
                x[i] = c[i];               // update solution array

        norm[k] = normL2(S, x);
    }
    *avg_time = *avg_time / MAXIT;
    free(c);
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
    return (real_t)sqrtf(norm);         // returns norm of residue array
}
