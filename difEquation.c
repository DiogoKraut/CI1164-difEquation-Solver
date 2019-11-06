/*! \file difEquation.c
    \brief Implementation of the PDE solver */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
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
    real_t num_pts = nx*ny; // Number of points in the mesh

    /* Find step sizes hx and hy */
    hx = M_PI / (nx + 1);
    hy = M_PI / (ny + 1);

    /* For simplification */
    real_t shx = hx*hx;
    real_t shy = hy*hy;;

    /* Calculate diagonal constants */
    // Main diagonal
    real_t md;
    md = 4*shx + 4*shy + 8*shx*shy*M_PI*M_PI;

    // Superior diagonal
    real_t sd;
    sd = hx*shy - 2*shy;

    // Inferior diagonal
    real_t id;
    id = -hx*shy - 2*shy;

    // Removed superior diagonal
    real_t rsd;
    rsd = -2*shx + shx*hy;

    // Removed inferior diagonal
    real_t rid;
    rid = -2shx - shx*hy;

    /* Fill main diagonal */
    for(i = 0; i < num_pts; i++)
        S->A[i][i] = md;

    /* Fill superior diagonal */
    for(i = 0; i < num_pts - 1; i++)
        S->A[i][i+1] = sd;

    /* Fill inferior diagonal */
    for(i = 1; i < num_pts; i++)
        S->A[i][i-1] = id;

    /* Fill removed superior diagonal */
    for(i = 0; i < (nx*ny)-nx; i++)
        S->A[i][i+nx] = rsd;

    /* Fill removed inferior diagonal */
    for(i = nx; i < num_pts; i++)
        S->A[i][i-nx] = rid;

    /* Iterate through the mesh */
    int k = 0;
    for(j = 0; j < ny ; j++) {
        for(i = 0; i < nx; i++) {
            x = i*hx;
            y = j*hy;
            S->b[k] = 8*M_PI*M_PI*shx*shy*(sin(2*M_PI*x)*sinh(M_PI*y) + sin(2*M_PI*(M_PI-x))*sinh(M_PI*(M_PI-y)) );

            /* Border conditions */
            if(i == 0 && k >= 1)
                S->A[k][k-1]  = 0.0;

            if(i == nx-1 && k < num_pts-1)
                S->A[k][k+1]  = 0.0;

            if(j == 0 && k >= nx) {
                S->A[k][k-nx] = 0.0;
                S->b[k]      -=  sin(2*M_PI*(M_PI-x))*sinh(M_PI*M_PI);
            }

            if(j == ny-1 && k < num_pts-nx) {
                S->A[k][k+nx] = 0.0;
                S->b[k]      -= sin(2*M_PI*x)*sinh(M_PI*M_PI);
            }

            k++;
        }
    }
    return 0;
}
