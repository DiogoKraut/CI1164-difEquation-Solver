/*! \file difEquation.c
    \brief Implementation of the PDE solver */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <likwid.h>
#include "difEquation.h"
#include "linSystem.h"

/*!
    \brief Generates the system of equations related to a partial diferential
    equation.

    Using the finite diferences method, generate the system of equations related
    to the partial diferential equation defined in the Project 1 especification.
    \param[out] S structure that stores a linear system of equations
    \param nx number of points in the X axis of the mesh
    \param ny number of points in the Y axis of the mesh
*/
int difEquation(linSystem_t *S, int nx, int ny) {
    int i, j;
    real_t hx, hy, x, y;

    /* Find step sizes hx and hy */
    hx = M_PI / (nx - 1);
    hy = M_PI / (ny - 1);

    /* For simplification */
    real_t shx = hx*hx;
    real_t shy = hy*hy;;

    /* Fill diagonals */

    real_t md;   // main
    md = 4*shx + 4*shy + 8*shx*shy*M_PI*M_PI;
    for(i = 0; i < S->n; i++)
        S->md[i] = md;

    real_t sd;   // superior
    sd = hx*shy - 2*shy;
    for(i = 0; i < S->n - 1; i++)
        S->sd[i] = sd;

    real_t id;   // inferior
    id = -hx*shy - 2*shy;
    for(i = 1; i < S->n; i++)
        S->id[i] = id;

    real_t rsd;  // removed superior
    rsd = -2*shx + shx*hy;
    for(i = 0; i < (nx*ny)-nx; i++)
        S->rsd[i] = rsd;

    real_t rid;  // removed inferior
    rid = -2*shx - shx*hy;
    for(i = nx; i < S->n; i++)
        S->rid[i] = rid;

    /* Iterate through the mesh */
    int k = 0;
    for(j = 0; j < ny ; j++) {
        for(i = 0; i < nx; i++) {
            x = i*hx;
            y = j*hy;
            /* Set b array */
            S->b[k] = 8*M_PI*M_PI*shx*shy*(sin(2*M_PI*x)*sinh(M_PI*y) + sin(2*M_PI*(M_PI-x))*sinh(M_PI*(M_PI-y)) );

            /* Border conditions */
            if(i == 0 && k >= 1)
                S->id[k] = 0.0;

            if(i == nx-1 && k < S->n-1)
                S->sd[k] = 0.0;

            if(j == 0 && k >= nx) {
                S->rid[k]= 0.0;
                S->b[k] -= sin(2*M_PI*(M_PI-x))*sinh(M_PI*M_PI);
            }

            if(j == ny-1 && k < S->n-nx) {
                S->rsd[k]= 0.0;
                S->b[k] -= sin(2*M_PI*x)*sinh(M_PI*M_PI);
            }

            k++;
        }
    }
    return 0;
}