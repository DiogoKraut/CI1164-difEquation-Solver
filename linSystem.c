#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "linSystem.h"

int gaussSeidel(linSystem_t *S, real_t *x, real_t error) {
    int k, i, j, end;
    real_t sum;
    real_t *c = calloc(S->n, sizeof(real_t));

    for (k = 0, end = 1; k < MAXIT && end == 1; k++) {
        end = 0;
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

            /* Check stop condition */

            if(fabsf(c[i] - x[i]) > error)
                end = 1;
        }
        if(end == 1)
            for(i = 0; i < S->n; i++)
                x[i] = c[i];               // update solution array
    }
    free(c);
    if(k >= MAXIT)
        return 1;
    return 0;
}
