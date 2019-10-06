#include <stdlib.h>
#include "gnuprint.h"
#include "linSystem.h"
#include "difEquation.h"
#include "parse.h"

int main(int argc, char *const argv[]) {
    OPT_ARGS_t opt;
    double avg_time;
    int int_count, i;

    initOPTS(&opt);

    parseMain(argc, argv, &opt);

    /* Memory allocation */
    linSystem_t *S = malloc(sizeof(linSystem_t));
    S->A = malloc(opt.NX*opt.NY*sizeof(real_t*));
    for(i = 0; i < opt.NX*opt.NY; i++)
        S->A[i] = calloc(opt.NX*opt.NY, sizeof(real_t));
    S->b = calloc(opt.NX*opt.NY, sizeof(real_t));
    S->n = opt.NX*opt.NY;

    real_t *x = calloc(S->n, sizeof(real_t));
    real_t *norm = calloc(MAXIT, sizeof(real_t));

    FILE *fp;
    if(*opt.OUTPUT_FILE == '\0')
        fp = stdout;
    else
        fp = fopen(opt.OUTPUT_FILE, "w");

    difEquation(S, opt.NX, opt.NY);

    if(gaussSeidel(S, x, 0.000001, &avg_time, norm, &int_count) == 1) {
        fprintf(stderr, "Method didnt converge\n");
        return 1;
    }

    gnuplot(x, S->n, avg_time, norm, int_count, fp);

    /* Free */
    for(i = 0; i < opt.NX*opt.NY; i++)
        free(S->A[i]);
    free(S->b);
    free(S->A);
    free(S);
    free(x);
    free(norm);



    return 0;
}
