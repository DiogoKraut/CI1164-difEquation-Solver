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
int gaussSeidel(linSystem_t *S, real_t *x, real_t error, double *avg_time, real_t *norm, int *int_count) {
    int k, i, j;
    int end; // Flag indicates when error is within acceptable range
    double t1;
    real_t sum, diff, res;

    *avg_time = 0.0;
    real_t *c = calloc(S->n, sizeof(real_t));

    for (k = 0, end = 1; k < MAXIT && end == 1; k++) {
        t1 = timestamp();
        end = 0;
        res = 0.0;
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

            diff = fabsf(c[i] - x[i]);
            if(diff > res)
                res = diff;
            /* Check stop condition */
            if(fabsf(c[i] - x[i]) > error)
                end = 1;
        }
        if(end == 1)
            for(i = 0; i < S->n; i++)
                x[i] = c[i];               // update solution array
        *avg_time += timestamp() - t1;
        norm[k] = res;
    }
    *int_count = k;
    t1 = t1 / k;
    free(c);
    // printf("%d\n", k);
    // printf("Solution:\n");
    // for(i = 0; i < S->n; i++){
    //     printf("x%d: %4.2f  ", i, x[i]);
    // }
    if(k >= MAXIT)
        return 1;
    return 0;
}
