/*! \file gnuprint.c
    \brief File used to generate gnuplot compatible data on the PDE
*/
#include <stdio.h>
#include <math.h>
#include "gnuprint.h"
#include "linSystem.h"

/*!
    \brief Function prints the results of the pde

    \param sol solution array
    \param n number of items in sol
    \param avg_time average iteration time of gaussSeidel method
    \param norm array containing the norm of each gaussSeidel iteration
    \param fp out put file
*/
int gnuplot(real_t *sol, int nx, int ny, double avg_time, real_t *norm, FILE *fp) {
    int i, j;
    real_t hx, hy, x, y;
    hx = M_PI / (nx - 1);
    hy = M_PI / (ny - 1);
    /* Print norm and avg time */
    fprintf(fp, "#################################\n");
    fprintf(fp, "#Tempo Método GS: %lf ms\n#\n", avg_time);
    fprintf(fp, "#################################\n");

    fprintf(fp, "#Norma L2 do Residuo\n");
    for(i = 0; i < MAXIT; i++)
        fprintf(fp, "#i=%d: %f\n", i+1, norm[i]);
    fprintf(fp, "#################################\n");

    /* Print data */
    // x = y = 0.0;
    // for(i = 0; i < ny; i++) {
    //     for(j = 0; j < nx; j ++) {
    //         fprintf(fp, "%f %f %f\n", x, y, sol[i*nx + j]);
    //         x += hx;
    //     }
    //     x = 0.0;
    //     y += hy;
    // }
    // return 0;

}
