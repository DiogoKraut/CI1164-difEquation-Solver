/*! \file difEquation.c */

#include <stdlib.h>
#include <math.h>
#include "difEquation.h"
#include "linSystem.h"

/*!
    \brief Generates the system of equations related to a partial diferential
    equation.

    Using the finite diferences method, generate the system of equations realted
    to the partial diferential equation defined in the Project 1 especification.
    \param[out] A structure that stores a linear system of equations
    \param nx number of points in the X axis of the mesh
    \param ny number of points in the Y axis of the mesh
*/
int difEquation(linSystem_t *A, int nx, int ny) {
    int i, j, k;
    real_t hx, hy;
    real_t M[nx+1][ny+1];
    /* Find step sizes hx and hy */
    hx = M_PI / (nx + 1);
    hy = M_PI / (ny + 1);

    /* Iterate through the mesh */
    for(i = 0; i < nx + 1; i++) {
        for(j = 0; j < ny + 1; j++) {

        }
    }
}
