/*! \file gnuprint.c
    \brief File used to generate gnuplot compatible data on the PDE
*/
#include <stdio.h>
#include "gnuprint.h"
#include "linSystem.h"

/*!
    \brief Function prints the results of the pde

    \param sol solution array
    \param n number of items in sol
    \param avg_time average iteration time of gaussSeidel method
    \param norm array containing the norm of each gaussSeidel iteration
    \param int_count gaussSeidel iteration count
    \param fp out put file
*/
int gnuplot(real_t *sol, unsigned int n, double avg_time, real_t *norm, int int_count, FILE *fp) {
    int i;
    /* Print norm and avg time */
    fprintf(fp, "#######################\n");
    fprintf(fp, "#Tempo Método GS: %lf\n#\n", avg_time);
    fprintf(fp, "#Norma L2 do Residuo\n");
    for(i = 0; i < int_count; i ++)
        fprintf(fp, "#i=%d: %f\n", i+1, norm[i]);
    fprintf(fp, "#######################\n");

    /* Print data */
    for(i = 0; i < n - 1; i++) {
        fprintf(fp, "%f,", sol[i]);
    }
    fprintf(fp, "%f\n", sol[n - 1]);
    return 0;

}
