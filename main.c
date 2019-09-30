#include <stdlib.h>
#include "linSystem.h"
#include "difEquation.h"
#include "parse.h"

int main(int argc, char *const argv[]) {
    OPT_ARGS_t opt;
    initOPTS(&opt);

    parseMain(argc, argv, &opt);
    printf("NX:%d \nNT:%d\nMI:%d\n", opt.NX, opt.NY, opt.MAXITER);
    return 0;
}
