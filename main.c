#include <stdlib.h>
#include "linSystem.h"
#include "difEquation.h"
#include "parse.h"

int main(int argc, char *const argv[]) {
    OPT_ARGS_t opt;
    initOPTS(&opt);

    parseMain(argc, argv, &opt);
    printf("NX:%d \nNY:%d\nMI:%d\n", opt.NX, opt.NY, opt.MAXITER);

    linSystem_t *S = malloc(sizeof(linSystem_t));
    difEquation(S, opt.NX, opt.NY);
    real_t *x = calloc(S->n, sizeof(real_t));
    gaussSeidel(S, x, 0.000001);
    return 0;
}
