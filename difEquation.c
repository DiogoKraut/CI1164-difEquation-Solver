/*! \file difEquation.c */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "difEquation.h"
#include "linSystem.h"

/*!
    \brief Generates the system of equations related to a partial diferential
    equation.

    Using the finite diferences method, generate the system of equations realted
    to the partial diferential equation defined in the Project 1 especification.
    \param[out] S structure that stores a linear system of equations
    \param nx number of points in the X axis of the mesh
    \param ny number of points in the Y axis of the mesh
*/
int difEquation(linSystem_t *S, int nx, int ny) {
    int i, j, k;
    real_t hx, hy;
    real_t num_pts = nx*ny; // Number of points in the mesh
    S->A = malloc(num_pts*sizeof(real_t*));
    for(i = 0; i < num_pts; i++)
        S->A[i] = calloc(num_pts, sizeof(real_t));

    S->b = calloc(num_pts, sizeof(real_t));
    S->n = num_pts;
    /* Find step sizes hx and hy */
    hx = M_PI / (nx - 1);
    hy = M_PI / (ny - 1);

    /* Find the diagonal constants */
    real_t shx = powf(hx, 2.0);
    real_t chx = powf(hx, 3.0);
    real_t shy = powf(hy, 2.0);
    real_t chy = powf(hy, 3.0);

    // Main diagonal
    real_t md;
    md = 8*(hx*chy + chx*hy + 2*chx*chy*powf(M_PI, 2.0));

    // Superior diagonal
    real_t sd;
    sd = 2*(shx*chy - 2*hx*chy);

    // Inferior diagonal
    real_t id;
    id = 2*(-2*hx*chy - shx*chy);

    // Removed superior diagonal
    real_t rsd;
    rsd = 2*(chx*shy - 2*chx*hy);

    // Removed inferior diagonal
    real_t rid;
    rid = 2*(-2*chx*hy - chx*shy);

    /* Fill main diagonal */
    for(i = 0; i < num_pts; i++) {
        S->A[i][i] = md;
        printf("%d\n", i);
    }

    /* Fill superior diagonal */
    for(i = 0; i < num_pts - 1; i++) {
        S->A[i][i+1] = sd;
    }

    /* Fill inferior diagonal */
    for(i = 1; i < num_pts; i++) {
        S->A[i][i-1] = id;
    }

    /* Fill removed superior diagonal */
    for(i = 0; i < nx*ny; i++) {
        S->A[i][i+nx] = rsd;
    }

    /* Fill removed inferior diagonal */
    for(i = nx; i < num_pts; i++) {
        S->A[i][i-nx] = rid;
    }
    /* Iterate through the mesh */
    int k = 0;
    for(j = 0; j < ny ; j++) {
        for(i = 0; i < nx; i++) {
            real_t x = i*hx;
            real_t y = j*hy;
            // b[k] = 16*powf(M_PI, 2.0)*chx*chy*(sin(2*M_PI*x)*sinh(M_PI*y) + sin(2*M_PI*(M_PI-x))sinh(M_PI*(M_PI-y)));

            if(i = 0) {
                S->A[i][j] = 1.0;
                b[k] = 0;
            } else if()

        }
    }
    printf("md:%7.2f\n", md);
    printf("sd:%7.2f\n", sd);
    printf("id:%7.2f\n", id);
    printf("rsd:%7.2f\n", rsd);
    printf("rid:%7.2f\n", rid);

    for(i = 0; i < num_pts; i++) {
        for(j = 0; j < num_pts-1; j++) {
            printf("%7.2f ", S->A[i][j]);
        }
        printf("%7.2f\n", S->A[i][j]);
    }
    return 0;
}
